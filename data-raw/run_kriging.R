# data-raw/run_ata_kriging.R
#
# Area-to-area kriging with observation-level measurement error:
#
#   y_i = average_{D_i} z(s) + eps_i
#   eps_i ~ N(0, sigma_i^2)
#
#   Cov(y_i, y_j) =
#     C_bar_spatial(D_i, D_j) + sigma_i^2 * 1(i = j)
#
# Variogram fitting uses:
#
#   E[0.5 * (y_i - y_j)^2]
#     = gamma_A(D_i, D_j) + 0.5 * (sigma_i^2 + sigma_j^2)
#
# where:
#
#   gamma_A(D_i, D_j)
#     = 0.5 * {K_ii + K_jj - 2 K_ij}
#
# and K_ab is the average latent covariance between supports D_a and D_b.
#
# Important:
#   - The latent covariance has no pointwise nugget.
#   - Measurement error is observation-level, not prediction-level.
#   - Prediction variances do NOT include sigma_i^2 unless predicting a future noisy observation.
#
# Memory-safe: N_WORKERS=4, K_CHUNK=20

library(spatintegrate)
library(goebel2026)
library(sf)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------

N_PER_TRIANGLE  <- 16L
N_WORKERS       <- 4L
N_BINS          <- 20L
N_PAIRS_PER_BIN <- 200L
MAX_DIST        <- 15000
CHUNK           <- 500L
K_CHUNK         <- 20L
T_CHUNK         <- 50L

timing <- list()

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

message("== Loading data ==")

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo

soundings_sf  <- d_shared$soundings_proj
fine_grid_sf  <- d_shared$fine_grid_buffered
y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true
X_obs_water   <- d_shared$X_obs_water

target_idx <- which(fine_grid_sf$n_intersects > 0)
target_sf  <- fine_grid_sf[target_idx, ]

n        <- length(y_alb)
n_target <- length(target_idx)
crs_proj <- sf::st_crs(soundings_sf)

message(sprintf("  Soundings: %d    Target cells: %d", n, n_target))

cent_mat <- sf::st_coordinates(
  sf::st_centroid(sf::st_geometry(soundings_sf))
)

# ------------------------------------------------------------------------------
# 2. Pre-compute QMC sample points for each sounding footprint
# ------------------------------------------------------------------------------

message("\n== Step 1: Pre-computing QMC sample points ==")
t0 <- proc.time()

qmc_base      <- spatintegrate::generate_qmc_unit_square(N_PER_TRIANGLE)
sounding_geom <- sf::st_geometry(soundings_sf)

get_poly_pts <- function(poly_sfg, crs) {
  tris <- spatintegrate::triangulate_sf(sf::st_sfc(poly_sfg, crs = crs))

  if (length(tris) == 0L) {
    return(matrix(numeric(0), ncol = 2L))
  }

  do.call(rbind, lapply(seq_along(tris), function(k) {
    spatintegrate::map_unit_square_to_triangle(
      qmc_base,
      spatintegrate::get_triangle_coords(tris[[k]])
    )
  }))
}

future::plan(future::multisession, workers = N_WORKERS)

sounding_pts <- future.apply::future_lapply(
  seq_len(n),
  function(i) get_poly_pts(sounding_geom[[i]], crs_proj),
  future.seed     = TRUE,
  future.packages = "spatintegrate",
  future.globals  = list(
    sounding_geom               = sounding_geom,
    crs_proj                    = crs_proj,
    qmc_base                    = qmc_base,
    get_poly_pts                = get_poly_pts,
    triangulate_sf              = spatintegrate::triangulate_sf,
    get_triangle_coords         = spatintegrate::get_triangle_coords,
    map_unit_square_to_triangle = spatintegrate::map_unit_square_to_triangle
  )
)

future::plan(future::sequential)

n_pts_vec <- vapply(sounding_pts, nrow, integer(1L))
n_pts     <- as.integer(median(n_pts_vec))

message(sprintf(
  "  Points per sounding: median=%d  range=[%d,%d]",
  n_pts, min(n_pts_vec), max(n_pts_vec)
))

S_all <- do.call(rbind, sounding_pts)
bb    <- rowSums(S_all^2)

timing$sample_pts_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$sample_pts_sec))

# ------------------------------------------------------------------------------
# 2b. Observation-error variance factor
#
# sigma2_i = SIGMA2_BASE * obs_var_factor_i
#
# Default below: constant per sounding.
#
# If you have known retrieval uncertainty, footprint quality scores, or an
# effective number of independent sub-measurements, replace obs_var_factor with
# that quantity, normalized to have median 1.
# ------------------------------------------------------------------------------

obs_var_factor <- rep(1, n)

# Optional geometry-based alternative:
#
# obs_area <- as.numeric(sf::st_area(soundings_sf))
# obs_var_factor <- median(obs_area, na.rm = TRUE) / obs_area
# obs_var_factor <- obs_var_factor / median(obs_var_factor, na.rm = TRUE)

stopifnot(length(obs_var_factor) == n)
stopifnot(all(is.finite(obs_var_factor)))
stopifnot(all(obs_var_factor > 0))

# ------------------------------------------------------------------------------
# 3. Empirical areal variogram
# ------------------------------------------------------------------------------

message("\n== Step 2: Empirical areal variogram ==")
t0 <- proc.time()

breaks    <- seq(0, MAX_DIST, length.out = N_BINS + 1L)
bin_mids  <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
emp_sum   <- numeric(N_BINS)
emp_count <- integer(N_BINS)

for (i_start in seq(1, n, by = CHUNK)) {
  i_end <- min(i_start + CHUNK - 1L, n)
  i_idx <- i_start:i_end

  dx <- outer(cent_mat[i_idx, 1], cent_mat[, 1], "-")
  dy <- outer(cent_mat[i_idx, 2], cent_mat[, 2], "-")
  h  <- sqrt(dx^2 + dy^2)

  gv <- outer(y_alb[i_idx], y_alb, function(a, b) 0.5 * (a - b)^2)

  for (k in seq_along(i_idx)) {
    i_global <- i_idx[k]
    h[k,  seq_len(i_global)] <- NA_real_
    gv[k, seq_len(i_global)] <- NA_real_
  }

  bin_idx <- findInterval(h, breaks, rightmost.closed = TRUE)
  valid   <- !is.na(h) & bin_idx >= 1L & bin_idx <= N_BINS

  for (b in seq_len(N_BINS)) {
    sel_b <- valid & bin_idx == b

    if (any(sel_b)) {
      emp_sum[b]   <- emp_sum[b] + sum(gv[sel_b], na.rm = TRUE)
      emp_count[b] <- emp_count[b] + sum(sel_b)
    }
  }
}

emp_gamma <- ifelse(emp_count > 0L, emp_sum / emp_count, NA_real_)

timing$emp_vgm_sec <- as.numeric((proc.time() - t0)["elapsed"])

message(sprintf("  Elapsed: %.1f sec", timing$emp_vgm_sec))
message("  Bin counts:      ", paste(emp_count, collapse = " "))
message("  Empirical gamma: ", paste(round(emp_gamma, 4), collapse = " "))

# ------------------------------------------------------------------------------
# 4. Select representative pairs, stratified by distance bin
#
# Each representative pair gets its own empirical value and its own geometry-exact
# ATA theoretical value.
# ------------------------------------------------------------------------------

message("\n== Step 3: Selecting representative pairs ==")

set.seed(2026L)

all_rep_i <- integer(0)
all_rep_j <- integer(0)

for (i_start in seq(1, n, by = CHUNK)) {
  i_end <- min(i_start + CHUNK - 1L, n)
  i_idx <- i_start:i_end

  dx <- outer(cent_mat[i_idx, 1], cent_mat[, 1], "-")
  dy <- outer(cent_mat[i_idx, 2], cent_mat[, 2], "-")
  h  <- sqrt(dx^2 + dy^2)

  for (k in seq_along(i_idx)) {
    h[k, seq_len(i_idx[k])] <- NA_real_
  }

  idx <- which(!is.na(h) & h < MAX_DIST, arr.ind = TRUE)

  if (nrow(idx) > 0L) {
    all_rep_i <- c(all_rep_i, i_idx[idx[, 1]])
    all_rep_j <- c(all_rep_j, idx[, 2])
  }
}

h_all <- sqrt(
  (cent_mat[all_rep_i, 1] - cent_mat[all_rep_j, 1])^2 +
    (cent_mat[all_rep_i, 2] - cent_mat[all_rep_j, 2])^2
)

bin_all <- findInterval(h_all, breaks, rightmost.closed = TRUE)

sel <- unlist(lapply(seq_len(N_BINS), function(b) {
  idx_b <- which(bin_all == b)

  if (length(idx_b) == 0L) {
    return(integer(0))
  }

  if (length(idx_b) <= N_PAIRS_PER_BIN) {
    idx_b
  } else {
    sample(idx_b, N_PAIRS_PER_BIN)
  }
}))

rep_i    <- all_rep_i[sel]
rep_j    <- all_rep_j[sel]
rep_h    <- h_all[sel]
rep_bin  <- bin_all[sel]
emp_pair <- 0.5 * (y_alb[rep_i] - y_alb[rep_j])^2

message(sprintf("  Total representative pairs: %d", length(rep_i)))
message(sprintf(
  "  Empirical pair values: mean=%.5f  range=[%.5f, %.5f]",
  mean(emp_pair), min(emp_pair), max(emp_pair)
))

# ------------------------------------------------------------------------------
# 5. Joint WLS fit of latent ATA variogram and observation error
#
# Model:
#
#   emp_pair_ij ~= gamma_A(D_i, D_j; psill, range)
#                  + 0.5 * (sigma_i^2 + sigma_j^2)
#
# where:
#
#   sigma_i^2 = SIGMA2_BASE * obs_var_factor_i
#
# This is the kriging-style nugget/error term, but placed at the observation
# support rather than as a pointwise latent-field nugget.
# ------------------------------------------------------------------------------

areal_covariance_pair <- function(pts_i, pts_j, psill, range) {
  if (nrow(pts_i) == 0L || nrow(pts_j) == 0L) {
    return(NA_real_)
  }

  dx <- outer(pts_i[, 1], pts_j[, 1], "-")
  dy <- outer(pts_i[, 2], pts_j[, 2], "-")
  h  <- sqrt(dx^2 + dy^2)

  mean(psill * exp(-h / range))
}

areal_variogram_pair <- function(pts_i, pts_j, psill, range) {
  K_ij <- areal_covariance_pair(pts_i, pts_j, psill, range)
  K_ii <- areal_covariance_pair(pts_i, pts_i, psill, range)
  K_jj <- areal_covariance_pair(pts_j, pts_j, psill, range)

  0.5 * (K_ii + K_jj - 2 * K_ij)
}

# Optional simple WLS weights. This keeps every representative pair equally
# weighted by default. You can replace this later with bin-count or robust weights.
wls_weights <- rep(1, length(emp_pair))

iter_count <- 0L

wls_objective <- function(log_theta) {
  iter_count <<- iter_count + 1L

  psill      <- exp(log_theta[1L])
  range      <- exp(log_theta[2L])
  sigma2base <- exp(log_theta[3L])

  sigma2_i <- sigma2base * obs_var_factor

  future::plan(future::multisession, workers = N_WORKERS)

  latent_gamma <- future.apply::future_mapply(
    function(i, j) {
      areal_variogram_pair(
        pts_i = sounding_pts[[i]],
        pts_j = sounding_pts[[j]],
        psill = psill,
        range = range
      )
    },
    rep_i, rep_j,
    future.seed     = TRUE,
    future.packages = character(0),
    future.globals  = list(
      sounding_pts          = sounding_pts,
      psill                 = psill,
      range                 = range,
      areal_covariance_pair = areal_covariance_pair,
      areal_variogram_pair  = areal_variogram_pair
    )
  )

  future::plan(future::sequential)

  noise_pair <- 0.5 * (sigma2_i[rep_i] + sigma2_i[rep_j])
  theo_vals  <- latent_gamma + noise_pair

  valid <- is.finite(theo_vals) & is.finite(emp_pair) & is.finite(wls_weights)

  if (!any(valid)) {
    return(.Machine$double.xmax)
  }

  resid <- emp_pair[valid] - theo_vals[valid]
  obj   <- sum(wls_weights[valid] * resid^2)

  message(sprintf(
    "  iter %3d: psill=%.5f  range=%.0fm  sigma2base=%.8f  obj=%.8f",
    iter_count, psill, range, sigma2base, obj
  ))

  obj
}

message("\n== Step 4: Joint ATA variogram + observation-error fit ==")
t0 <- proc.time()

psill_init      <- var(y_alb, na.rm = TRUE)
range_init      <- MAX_DIST * 0.3
sigma2base_init <- max(0.05 * var(y_alb, na.rm = TRUE), 1e-8)

opt <- optim(
  par     = c(log(psill_init), log(range_init), log(sigma2base_init)),
  fn      = wls_objective,
  method  = "Nelder-Mead",
  control = list(maxit = 250L, reltol = 1e-4)
)

COV_NUGGET  <- 0
COV_PSILL   <- exp(opt$par[1L])
COV_RANGE   <- exp(opt$par[2L])
SIGMA2_BASE <- exp(opt$par[3L])
sigma2_i    <- SIGMA2_BASE * obs_var_factor

timing$vgm_fit_sec <- as.numeric((proc.time() - t0)["elapsed"])

message(sprintf("  Elapsed: %.1f sec", timing$vgm_fit_sec))
message(sprintf(
  "  Fitted latent covariance: psill=%.5f  range=%.1fm",
  COV_PSILL, COV_RANGE
))
message(sprintf(
  "  Fitted observation error: sigma2_base=%.8f  median sigma2_i=%.8f",
  SIGMA2_BASE, median(sigma2_i)
))
message(sprintf(
  "  True injected noise: %.8f  ratio median/true=%.2f",
  d_albedo$sigma_eps^2,
  median(sigma2_i) / d_albedo$sigma_eps^2
))

saveRDS(
  list(
    COV_NUGGET     = COV_NUGGET,
    COV_PSILL      = COV_PSILL,
    COV_RANGE      = COV_RANGE,
    SIGMA2_BASE    = SIGMA2_BASE,
    sigma2_i       = sigma2_i,
    obs_var_factor = obs_var_factor,
    opt            = opt,
    emp_gamma      = emp_gamma,
    emp_count      = emp_count,
    bin_mids       = bin_mids
  ),
  file = "ata_vgm_params.rds"
)

message("  Params saved to ata_vgm_params.rds")

# ------------------------------------------------------------------------------
# 6. Diagnostic: compare empirical pairs to fitted ATA + obs-error model
# ------------------------------------------------------------------------------

message("\n== Step 5: Variogram fit diagnostics ==")
t0 <- proc.time()

future::plan(future::multisession, workers = N_WORKERS)

latent_gamma_fit <- future.apply::future_mapply(
  function(i, j) {
    areal_variogram_pair(
      pts_i = sounding_pts[[i]],
      pts_j = sounding_pts[[j]],
      psill = COV_PSILL,
      range = COV_RANGE
    )
  },
  rep_i, rep_j,
  future.seed     = TRUE,
  future.packages = character(0),
  future.globals  = list(
    sounding_pts          = sounding_pts,
    COV_PSILL            = COV_PSILL,
    COV_RANGE            = COV_RANGE,
    areal_covariance_pair = areal_covariance_pair,
    areal_variogram_pair  = areal_variogram_pair
  )
)

future::plan(future::sequential)

noise_pair_fit <- 0.5 * (sigma2_i[rep_i] + sigma2_i[rep_j])
theo_pair_fit  <- latent_gamma_fit + noise_pair_fit
fit_resid_pair <- emp_pair - theo_pair_fit

diag_by_bin <- lapply(seq_len(N_BINS), function(b) {
  idx_b <- which(rep_bin == b)

  if (length(idx_b) == 0L) {
    return(data.frame(
      bin = b,
      h_mid = bin_mids[b],
      n_pairs = 0L,
      emp = NA_real_,
      latent = NA_real_,
      noise = NA_real_,
      total = NA_real_,
      resid = NA_real_
    ))
  }

  data.frame(
    bin = b,
    h_mid = bin_mids[b],
    n_pairs = length(idx_b),
    emp = mean(emp_pair[idx_b], na.rm = TRUE),
    latent = mean(latent_gamma_fit[idx_b], na.rm = TRUE),
    noise = mean(noise_pair_fit[idx_b], na.rm = TRUE),
    total = mean(theo_pair_fit[idx_b], na.rm = TRUE),
    resid = mean(fit_resid_pair[idx_b], na.rm = TRUE)
  )
})

diag_by_bin <- do.call(rbind, diag_by_bin)

timing$vgm_diag_sec <- as.numeric((proc.time() - t0)["elapsed"])

message("  First few fitted variogram diagnostics:")
print(utils::head(diag_by_bin, 8L), row.names = FALSE)

message(sprintf("  Elapsed: %.1f sec", timing$vgm_diag_sec))

# ------------------------------------------------------------------------------
# 7. Build K_AA: observation-support latent covariance matrix
#
# K[i, j] = average covariance between footprint D_i and D_j.
# This matrix does NOT include observation error.
# ------------------------------------------------------------------------------

message("\n== Step 6: Building K matrix ==")
t0 <- proc.time()

cov_spatial <- function(h) {
  COV_PSILL * exp(-h / COV_RANGE)
}

n_chunks <- ceiling(n / K_CHUNK)

chunk_list <- lapply(seq_len(n_chunks), function(chunk_idx) {
  i_start <- (chunk_idx - 1L) * K_CHUNK + 1L
  i_end   <- min(chunk_idx * K_CHUNK, n)
  i_start:i_end
})

future::plan(future::multisession, workers = N_WORKERS)

K_chunks <- future.apply::future_lapply(
  chunk_list,
  function(i_idx) {
    np <- n_pts

    row_idx <- as.vector(outer(seq_len(np), (i_idx - 1L) * np, "+"))

    S_chunk <- S_all[row_idx, , drop = FALSE]
    n_chunk <- length(i_idx)

    aa <- rowSums(S_chunk^2)
    ab <- tcrossprod(S_chunk, S_all)

    H2 <- outer(aa, bb, "+") - 2 * ab
    H  <- sqrt(pmax(H2, 0))

    CV <- cov_spatial(H)

    col_grp  <- rep(seq_len(n), each = np)
    CV_j_avg <- t(rowsum(t(CV), col_grp) / np)

    row_grp <- rep(seq_len(n_chunk), each = np)
    K_chunk <- rowsum(CV_j_avg, row_grp) / np

    list(i_idx = i_idx, K_chunk = K_chunk)
  },
  future.seed    = TRUE,
  future.globals = list(
    S_all       = S_all,
    bb          = bb,
    n           = n,
    n_pts       = n_pts,
    K_CHUNK     = K_CHUNK,
    COV_PSILL   = COV_PSILL,
    COV_RANGE   = COV_RANGE,
    cov_spatial = cov_spatial
  )
)

future::plan(future::sequential)

K <- matrix(0, nrow = n, ncol = n)

for (res in K_chunks) {
  K[res$i_idx, ] <- res$K_chunk
}

K <- 0.5 * (K + t(K))

timing$K_build_sec <- as.numeric((proc.time() - t0)["elapsed"])

message(sprintf("  K build elapsed: %.1f sec", timing$K_build_sec))
message(sprintf(
  "  K diagonal summary: min=%.6f median=%.6f max=%.6f",
  min(diag(K)), median(diag(K)), max(diag(K))
))

saveRDS(K, file = "ata_K.rds")
message("  K saved to ata_K.rds")

# ------------------------------------------------------------------------------
# 8. Cholesky of observation covariance
#
# Observation covariance:
#
#   K_obs = K_AA + diag(sigma2_i)
#
# Measurement error appears here only. It is not part of the latent prediction
# target self-variance.
# ------------------------------------------------------------------------------

message("\n== Step 7: Cholesky of observation covariance ==")
t0 <- proc.time()

K_obs <- K
diag(K_obs) <- diag(K_obs) + sigma2_i + 1e-10 * mean(diag(K))

CK <- tryCatch(
  chol(K_obs),
  error = function(e) {
    message("  Cholesky failed. Adding jitter...")
    diag(K_obs) <<- diag(K_obs) + 1e-4 * mean(diag(K_obs))
    chol(K_obs)
  }
)

timing$chol_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Cholesky elapsed: %.1f sec", timing$chol_sec))

ones   <- rep(1, n)
Kinv_1 <- backsolve(CK, forwardsolve(t(CK), ones))
mu_ok  <- sum(Kinv_1 * y_alb) / sum(Kinv_1)

e      <- y_alb - mu_ok
Kinv_e <- backsolve(CK, forwardsolve(t(CK), e))

message(sprintf("  Ordinary kriging mean: %.4f", mu_ok))

# ------------------------------------------------------------------------------
# 9. Prediction covariance vectors
#
# Current version preserves your centroid target approximation:
#
#   k_mat[target, obs] = Cov(z(target centroid), average over D_i)
#
# This is latent covariance only; no observation noise is added to k_mat.
# ------------------------------------------------------------------------------

message("\n== Step 8: Prediction covariance vectors ==")
t0 <- proc.time()

target_cents <- sf::st_coordinates(
  sf::st_centroid(sf::st_geometry(target_sf))
)

bb2   <- rowSums(target_cents^2)
k_mat <- matrix(0, nrow = n_target, ncol = n)
idx   <- rep(seq_len(n), each = n_pts)

n_t_chunks <- ceiling(n_target / T_CHUNK)

for (chunk_idx in seq_len(n_t_chunks)) {
  j_start <- (chunk_idx - 1L) * T_CHUNK + 1L
  j_end   <- min(chunk_idx * T_CHUNK, n_target)
  j_idx   <- j_start:j_end

  tc_chunk  <- target_cents[j_idx, , drop = FALSE]
  bb2_chunk <- bb2[j_idx]

  ab <- tcrossprod(S_all, tc_chunk)

  H2 <- outer(bb, bb2_chunk, "+") - 2 * ab
  H  <- sqrt(pmax(H2, 0))

  CV <- cov_spatial(H)

  k_mat[j_idx, ] <- t(rowsum(CV, idx) / n_pts)

  if (chunk_idx %% 20L == 0L || chunk_idx == n_t_chunks) {
    message(sprintf("  target chunk %d/%d", chunk_idx, n_t_chunks))
  }
}

timing$pred_cov_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$pred_cov_sec))

# ------------------------------------------------------------------------------
# 10. Predictions and latent prediction variances
#
# No observation-error variance is added to C_self. That would only be appropriate
# if predicting a future noisy measurement, not the latent target.
# ------------------------------------------------------------------------------

message("\n== Step 9: Predictions ==")
t0 <- proc.time()

mu_t <- mu_ok + as.numeric(k_mat %*% Kinv_e)

# Because the current prediction target is a point/centroid latent value,
# the self-covariance is the latent point sill.
#
# If you later switch to target-cell averages, replace C_self with the
# target-cell areal self-covariance K_**.
C_self <- COV_PSILL

Kinv_kt <- backsolve(CK, forwardsolve(t(CK), t(k_mat)))

kKinvk <- colSums(Kinv_kt * t(k_mat))
kKinv1 <- as.numeric(k_mat %*% Kinv_1)

ok_corr <- (1 - kKinv1)^2 / sum(Kinv_1)

var_t <- pmax(C_self - kKinvk + ok_corr, 0)
se_t  <- sqrt(var_t)

timing$predict_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$predict_sec))

# ------------------------------------------------------------------------------
# 11. Evaluate
# ------------------------------------------------------------------------------

tr_t <- y_latent_true[target_idx]
ns_t <- as.integer(Matrix::colSums(d_shared$A_flat[, target_idx] > 0))

ci_lo <- mu_t - 1.96 * se_t
ci_hi <- mu_t + 1.96 * se_t

coverage <- function(idx) {
  if (length(idx) == 0L) {
    return(NA_real_)
  }

  mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
}

resid <- mu_t - tr_t

rmse <- sqrt(mean(resid^2, na.rm = TRUE))

r2 <- 1 - sum(resid^2, na.rm = TRUE) /
  sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

message(sprintf(
  "\n  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  mean_SE=%.4f",
  rmse, r2, coverage(which(ns_t >= 1L)), mean(se_t, na.rm = TRUE)
))

message(sprintf(
  "  median sigma2_i estimated: %.8f  true injected: %.8f  ratio: %.2f",
  median(sigma2_i),
  d_albedo$sigma_eps^2,
  median(sigma2_i) / d_albedo$sigma_eps^2
))

# ------------------------------------------------------------------------------
# 12. Save
# ------------------------------------------------------------------------------

results_kriging_ata_albedo <- list(
  run_name = "kriging_ata_heteroskedastic_obs_error_albedo_330m",

  tags = list(
    resolution     = 330L,
    method         = "ata_kriging_joint_wls_obs_error",
    covariates     = "none",
    support        = "areal_obs_centroid_pred",
    cov_model      = "exponential_latent_no_point_nugget",
    obs_error      = "heteroskedastic_diagonal",
    obs_error_fit  = "joint_wls_variogram"
  ),

  timestamp    = Sys.time(),
  resolution_m = 330L,

  posterior_mean = mu_t,
  posterior_se   = se_t,
  ci_lower       = ci_lo,
  ci_upper       = ci_hi,

  mu_ok = mu_ok,

  sigma2_i       = sigma2_i,
  sigma2_base    = SIGMA2_BASE,
  obs_var_factor = obs_var_factor,

  cov_params = list(
    nugget = 0,
    psill  = COV_PSILL,
    range  = COV_RANGE
  ),

  emp_gamma   = emp_gamma,
  emp_count   = emp_count,
  bin_mids    = bin_mids,
  diag_by_bin = diag_by_bin,

  rmse = rmse,
  r2   = r2,

  coverage_95_all   = coverage(seq_along(mu_t)),
  coverage_95_obs   = coverage(which(ns_t >= 1L)),
  coverage_95_dense = coverage(which(ns_t >= 20L)),

  n_soundings_per_pixel = ns_t,
  timing                = timing
)

usethis::use_data(results_kriging_ata_albedo, overwrite = TRUE)

# ------------------------------------------------------------------------------
# 13. Summary
# ------------------------------------------------------------------------------

message("\n=== ATA kriging timing summary ===")

total <- 0

for (nm in names(timing)) {
  message(sprintf("  %-25s  %6.1f sec", nm, timing[[nm]]))
  total <- total + timing[[nm]]
}

message(sprintf("  %-25s  %6.1f sec", "TOTAL", total))

message("\n=== Final parameter summary ===")
message(sprintf("  psill:              %.8f", COV_PSILL))
message(sprintf("  range:              %.2f", COV_RANGE))
message(sprintf("  sigma2_base:        %.8f", SIGMA2_BASE))
message(sprintf("  median sigma2_i:    %.8f", median(sigma2_i)))
message(sprintf("  true sigma2_e:      %.8f", d_albedo$sigma_eps^2))
message(sprintf("  RMSE:               %.6f", rmse))
message(sprintf("  R2:                 %.6f", r2))
message(sprintf("  coverage obs:       %.6f", coverage(which(ns_t >= 1L))))
