# data-raw/run_ata_kriging_matched.R
#
# ATA kriging with covariance parameters matched to the CV-tuned SAR model.
#
# Argument: the A matrix (area-fraction aggregation) is a good approximation
# to full QMC areal integration. We verify this by running ATA kriging with
# the same Matérn covariance structure implied by the SAR model. If the
# posteriors are nearly identical, it validates that A captures the
# change-of-support geometry well enough that expensive QMC integration
# adds nothing.
#
# Covariance model: Matérn nu=1
#   C(h) = SIGMA2 * (1 + KAPPA*h) * exp(-KAPPA*h)
#
# Parameter matching:
#   KAPPA from Table 5.1 (SAR -> Matérn SPDE, normalised [0,1]^2 units,
#   then converted to metres):
#     kappa2_unit = (2*d / h_unit^2) * (1 - rho) / rho
#     KAPPA       = sqrt(kappa2_unit) / L
#
#   SIGMA2 = C(0), the pointwise marginal variance in albedo units, matched
#   directly to the SAR marginal variance:
#     SIGMA2 = sigma2_SAR * mean([Q^{-1}]_{jj}) for interior cells
#   This avoids the SPDE sigma2_M parameterisation and the queen vs rook
#   adjacency correction factor entirely.
#
#   SIGMA2_E = sigma2e from CV fit (profiled analytically, not method-of-moments).
#
# Prerequisites:
#   results_10m_no_cov_rho_cv  (run_main_results_10m.R, section 2b)
#
# Outputs:
#   results_ata_kriging_matched  -- saved via usethis::use_data()
#   fig_ata_vs_sar_matched.pdf   -- scatter plots (means vs means, SEs vs SEs)

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

N_PER_TRIANGLE <- 128L
N_WORKERS      <- 4L
K_CHUNK        <- 20L
J_CHUNK        <- 20L   # inner chunk over j-soundings to limit memory
T_CHUNK        <- 50L

timing <- list()

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

message("== 1. Loading data ==")

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

soundings_sf  <- d_shared$soundings_proj
fine_grid_sf  <- d_shared$fine_grid_buffered
W_queen       <- d_shared$W_queen
y_latent_true <- d_albedo$y_latent_true
sigma_eps     <- d_albedo$sigma_eps

keep_idx     <- d_sif$keep_idx
mean_alb_10m <- goebel2026::soundings_augmented$mean_albedo[keep_idx]
set.seed(2026L)
y_alb_10m <- mean_alb_10m + rnorm(length(mean_alb_10m), sd = sigma_eps)

target_idx <- which(fine_grid_sf$n_intersects > 0)
target_sf  <- fine_grid_sf[target_idx, ]

n        <- length(y_alb_10m)
n_target <- length(target_idx)
crs_proj <- sf::st_crs(soundings_sf)

message(sprintf("  Soundings: %d    Target cells: %d", n, n_target))

# ------------------------------------------------------------------------------
# 2. Parameter conversion
# ------------------------------------------------------------------------------

message("== 2. Parameter conversion ==")

e_cv <- new.env()
data("results_10m_no_cov_rho_cv", package = "goebel2026", envir = e_cv)
cv <- e_cv$results_10m_no_cov_rho_cv

rho_cv     <- cv$rho_opt
phi_cv     <- cv$phi
sigma2e_cv <- cv$sigma2e
sigma2_SAR <- phi_cv * sigma2e_cv

message(sprintf("  SAR: rho=%.4f  phi=%.4f  sigma2e=%.4g  sigma2_SAR=%.4g",
                rho_cv, phi_cv, sigma2e_cv, sigma2_SAR))

# --- KAPPA from Table 5.1 (normalised [0,1]^2 units, then to metres) ----------
d_dim  <- 2L
nu     <- 1L
alpha  <- nu + d_dim / 2   # = 2 for nu=1, d=2
h_grid <- 330              # grid spacing in metres

# Grid domain size (used to normalise h to [0,1]^2 units)
bbox   <- sf::st_bbox(fine_grid_sf)
L      <- as.numeric(max(bbox["xmax"] - bbox["xmin"],
                         bbox["ymax"] - bbox["ymin"]))  # longer dimension in metres
h_unit <- h_grid / L   # h in [0,1]^2 units

kappa2_unit <- (2 * d_dim / h_unit^2) * (1 - rho_cv) / rho_cv
KAPPA       <- sqrt(kappa2_unit) / L   # convert to metres^{-1}
KAPPA2      <- KAPPA^2

message(sprintf("  L=%.0fm  h_unit=%.6f  kappa2_unit=%.4g  KAPPA=%.6g m^-1",
                L, h_unit, kappa2_unit, KAPPA))
message(sprintf("  Effective range: %.1f m", sqrt(8 * nu) / KAPPA))

# --- SIGMA2: pointwise marginal variance matched to SAR, in albedo units ------
#
# Instead of the SPDE sigma2_M parameterisation (which lives in normalised
# [0,1]^2 units and requires an eigensum to recover physical variance), we
# use the standard Matérn parameterisation directly:
#
#   C(h) = SIGMA2 * (1 + KAPPA*h) * exp(-KAPPA*h)
#
# where SIGMA2 = C(0) is the pointwise marginal variance in albedo units.
# We match SIGMA2 to the SAR marginal variance:
#   SAR marginal var = sigma2_SAR * mean([Q^{-1}]_{jj}) for interior cells
#
# This avoids the queen vs rook adjacency correction factor issue entirely --
# we just directly equate the pointwise variances.

# --- SIGMA2: sounding-level variance matching --------------------------------
#
# We match SIGMA2 so that the ATA kriging observation covariance K matches
# A Q^{-1} A' from the SAR model at the sounding level.
#
# At the pointwise level, sigma2_SAR * [Q^{-1}]_jj matches C(0) = SIGMA2.
# But the observation covariance K_ii (areal integral over sounding footprint)
# differs from [A Q^{-1} A']_ii by a scale factor, because:
#   - K is built with a continuous Matérn integral
#   - A Q^{-1} A' uses discrete area fractions
#
# We find SIGMA2 such that K_ii = [A Q^{-1} A']_ii by:
#   1. Computing [A Q^{-1} A']_ii for a sample of soundings
#   2. Computing K_ii at SIGMA2=1 for the same soundings
#   3. Setting SIGMA2 = [A Q^{-1} A']_ii / K_ii(SIGMA2=1)
#
# This ensures the two observation covariance matrices match in scale,
# which is the correct level for comparing the two models.

message("  Computing sounding-level variance match...")

# Step 1: A Q^{-1} A' for first 50 soundings
Q_sar    <- Matrix::forceSymmetric(
  Matrix::crossprod(Matrix::Diagonal(nrow(W_queen)) - rho_cv * W_queen)
)
idx_snd  <- seq_len(min(50L, n))
A_sub    <- d_shared$A_flat[idx_snd, , drop = FALSE]
QinvAt   <- Matrix::solve(Q_sar, t(A_sub))
AQinvAt  <- sigma2_SAR * as.matrix(A_sub %*% QinvAt)
mean_AQA <- mean(diag(AQinvAt))

# Step 2: K_ii at SIGMA2=1 for same soundings using QMC points
# (computed after sounding_pts is built -- deferred, use analytical approx)
# K_ii(SIGMA2=1) ≈ avg_{s,t in D_i} (1 + KAPPA*||s-t||) * exp(-KAPPA*||s-t||)
# For now use the ratio from the target cell Qinv approach as a fast estimate
deg      <- as.integer(Matrix::rowSums(W_queen > 0))
interior <- which(deg == max(deg))[seq_len(min(20L, sum(deg == max(deg))))]
Qinv_diag_int <- vapply(interior, function(j) {
  ej <- rep(0, nrow(Q_sar)); ej[j] <- 1
  as.numeric(Matrix::solve(Q_sar, ej)[j])
}, numeric(1L))

# SIGMA2 at unit scale (will be corrected after K is built)
# Initial estimate using pointwise matching -- corrected post-K-build
SIGMA2_pointwise <- sigma2_SAR * mean(Qinv_diag_int)

message(sprintf("  mean(AQA'_ii) = %.6g", mean_AQA))
message(sprintf("  SIGMA2 (pointwise) = %.6g", SIGMA2_pointwise))
message(sprintf("  Sounding scale factor (AQA'/K will be computed post-K-build)"))

# We'll set SIGMA2 correctly after K is built at SIGMA2=1
# For now use pointwise -- will rescale
SIGMA2   <- SIGMA2_pointwise
SIGMA2_E <- sigma2e_cv
C_SELF   <- SIGMA2

message(sprintf("  SIGMA2_E = %.6g", SIGMA2_E))

# ------------------------------------------------------------------------------
# 3. Pre-compute QMC sample points
# ------------------------------------------------------------------------------

message("\n== 3. QMC sample points ==")
t0 <- proc.time()

qmc_base      <- spatintegrate::generate_qmc_unit_square(N_PER_TRIANGLE)
sounding_geom <- sf::st_geometry(soundings_sf)

get_poly_pts <- function(poly_sfg, crs) {
  tris <- spatintegrate::triangulate_sf(sf::st_sfc(poly_sfg, crs = crs))
  if (length(tris) == 0L) return(matrix(numeric(0), ncol = 2L))
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
message(sprintf("  Points per sounding: median=%d  range=[%d,%d]  uniform=%s",
                n_pts, min(n_pts_vec), max(n_pts_vec),
                if (length(unique(n_pts_vec)) == 1L) "yes" else "no"))

S_all <- do.call(rbind, sounding_pts)
bb    <- rowSums(S_all^2)

timing$qmc_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$qmc_sec))

grp_all <- rep(seq_len(n), each = n_pts)

# ------------------------------------------------------------------------------
# 4. Build observation covariance matrix K
#
# K[i,j] = avg_{s in D_i, t in D_j} C(||s - t||)
# Chunked over sounding rows i to keep memory bounded.
# KAPPA and SIGMA2 passed explicitly as globals -- no closure issues.
# ------------------------------------------------------------------------------

message("\n== 4. Building K matrix ==")
t0 <- proc.time()

n_chunks   <- ceiling(n / K_CHUNK)
chunk_list <- lapply(seq_len(n_chunks), function(k) {
  i_start <- (k - 1L) * K_CHUNK + 1L
  i_end   <- min(k * K_CHUNK, n)
  i_start:i_end
})

future::plan(future::multisession, workers = N_WORKERS)
K_chunks <- future.apply::future_lapply(
  chunk_list,
  function(i_idx) {
    n_chunk <- length(i_idx)
    np      <- n_pts

    row_idx_i <- as.vector(outer(seq_len(np), (i_idx - 1L) * np, "+"))
    S_i       <- S_all[row_idx_i, , drop = FALSE]
    aa_i      <- rowSums(S_i^2)

    # Inner j-chunking to limit memory: process J_CHUNK soundings at a time
    K_chunk <- matrix(0, nrow = n_chunk, ncol = n)
    j_starts <- seq(1L, n, by = J_CHUNK)

    for (j_start in j_starts) {
      j_end   <- min(j_start + J_CHUNK - 1L, n)
      j_idx   <- j_start:j_end
      n_j     <- length(j_idx)

      row_idx_j <- as.vector(outer(seq_len(np), (j_idx - 1L) * np, "+"))
      S_j       <- S_all[row_idx_j, , drop = FALSE]
      bb_j      <- rowSums(S_j^2)

      # Distances: S_i (n_chunk*np rows) x S_j (n_j*np rows)
      ab  <- tcrossprod(S_i, S_j)          # (n_chunk*np) x (n_j*np)
      H2  <- outer(aa_i, bb_j, "+") - 2 * ab
      H   <- sqrt(pmax(H2, 0))

      kh  <- KAPPA * H
      CV  <- SIGMA2 * (1 + kh) * exp(-kh)  # (n_chunk*np) x (n_j*np)

      # Average over j-sounding cols
      grp_j    <- rep(seq_len(n_j), each = np)
      CV_j_avg <- t(rowsum(t(CV), grp_j) / np)   # (n_chunk*np) x n_j

      # Average over i-sounding rows
      grp_i    <- rep(seq_len(n_chunk), each = np)
      K_chunk[, j_idx] <- rowsum(CV_j_avg, grp_i) / np   # n_chunk x n_j
    }

    list(i_idx = i_idx, K_chunk = K_chunk)
  },
  future.seed    = TRUE,
  future.globals = list(
    S_all    = S_all,
    n        = n,
    n_pts    = n_pts,
    J_CHUNK  = J_CHUNK,
    KAPPA    = KAPPA,
    SIGMA2   = SIGMA2
  )
)
future::plan(future::sequential)

K <- matrix(0, nrow = n, ncol = n)
for (res in K_chunks) K[res$i_idx, ] <- res$K_chunk
K <- 0.5 * (K + t(K))

# --- Rescale SIGMA2 to match A Q^{-1} A' at sounding level ------------------
# K was built with SIGMA2_pointwise. Now rescale so that
# mean(diag(K_rescaled)) = mean(diag(A Q^{-1} A')) for the sample soundings.
# K scales linearly with SIGMA2, so:
#   SIGMA2_final = SIGMA2_pointwise * mean_AQA / mean(diag(K[idx_snd, idx_snd]))
mean_K_diag <- mean(diag(K[idx_snd, idx_snd]))
scale_factor <- mean_AQA / mean_K_diag
SIGMA2       <- SIGMA2_pointwise * scale_factor
C_SELF       <- SIGMA2
K            <- K * scale_factor   # rescale K in place

message(sprintf("  K diagonal (pre-scale): mean=%.4g", mean_K_diag))
message(sprintf("  A Q^{-1} A' diagonal:   mean=%.4g", mean_AQA))
message(sprintf("  Scale factor: %.4f", scale_factor))
message(sprintf("  SIGMA2 (sounding-matched) = %.6g", SIGMA2))
message(sprintf("  SNR = SIGMA2/SIGMA2_E = %.1f", SIGMA2 / SIGMA2_E))

timing$K_build_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  K build elapsed: %.1f sec", timing$K_build_sec))
message(sprintf("  K diagonal (post-scale): mean=%.4g  range=[%.4g, %.4g]",
                mean(diag(K)), min(diag(K)), max(diag(K))))

# ------------------------------------------------------------------------------
# 5. Cholesky of K_obs = K + SIGMA2_E * I
# ------------------------------------------------------------------------------

message("\n== 5. Cholesky ==")
t0 <- proc.time()

K_obs       <- K
diag(K_obs) <- diag(K_obs) + SIGMA2_E + 1e-10 * mean(diag(K))

CK <- tryCatch(
  chol(K_obs),
  error = function(e) {
    message("  Cholesky failed, adding jitter...")
    diag(K_obs) <<- diag(K_obs) + 1e-4 * mean(diag(K_obs))
    chol(K_obs)
  }
)

timing$chol_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Cholesky elapsed: %.1f sec", timing$chol_sec))

ones   <- rep(1, n)
Kinv_1 <- backsolve(CK, forwardsolve(t(CK), ones))
mu_ok  <- sum(Kinv_1 * y_alb_10m) / sum(Kinv_1)
resid  <- y_alb_10m - mu_ok
Kinv_r <- backsolve(CK, forwardsolve(t(CK), resid))
message(sprintf("  Ordinary kriging mean: %.4f", mu_ok))

# ------------------------------------------------------------------------------
# 6. Prediction cross-covariance k(s*, D_i)
#
# k[j,i] = avg_{t in D_i} C(||s*_j - t||)
# Target side: centroid of target cell (point support).
# Observation side: QMC integration over sounding footprints.
# ------------------------------------------------------------------------------

message("\n== 6. Prediction cross-covariance ==")
t0 <- proc.time()

target_cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(target_sf)))
bb2          <- rowSums(target_cents^2)
k_mat        <- matrix(0, nrow = n_target, ncol = n)

n_t_chunks <- ceiling(n_target / T_CHUNK)

for (chunk_idx in seq_len(n_t_chunks)) {
  j_start   <- (chunk_idx - 1L) * T_CHUNK + 1L
  j_end     <- min(chunk_idx * T_CHUNK, n_target)
  j_idx     <- j_start:j_end
  tc        <- target_cents[j_idx, , drop = FALSE]
  bb2_chunk <- bb2[j_idx]

  ab  <- tcrossprod(S_all, tc)
  H2  <- outer(bb, bb2_chunk, "+") - 2 * ab
  H   <- sqrt(pmax(H2, 0))

  kh  <- KAPPA * H
  CV  <- SIGMA2 * (1 + kh) * exp(-kh)   # (n*np) x T_CHUNK

  k_mat[j_idx, ] <- t(rowsum(CV, grp_all) / n_pts)

  if (chunk_idx %% 20L == 0L || chunk_idx == n_t_chunks)
    message(sprintf("  target chunk %d/%d", chunk_idx, n_t_chunks))
}

timing$pred_cov_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$pred_cov_sec))

# ------------------------------------------------------------------------------
# 7. Predictions and kriging variances
# ------------------------------------------------------------------------------

message("\n== 7. Predictions ==")
t0 <- proc.time()

mu_t <- mu_ok + as.numeric(k_mat %*% Kinv_r)

Kinv_kt <- backsolve(CK, forwardsolve(t(CK), t(k_mat)))
kKinvk  <- colSums(Kinv_kt * t(k_mat))
kKinv1  <- as.numeric(k_mat %*% Kinv_1)
ok_corr <- (1 - kKinv1)^2 / sum(Kinv_1)

var_t <- pmax(C_SELF - kKinvk + ok_corr, 0)
se_t  <- sqrt(var_t)

timing$predict_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing$predict_sec))

# ------------------------------------------------------------------------------
# 8. Evaluate
# ------------------------------------------------------------------------------

tr_t  <- y_latent_true[target_idx]
ns_t  <- as.integer(Matrix::colSums(d_shared$A_flat[, target_idx] > 0))
ci_lo <- mu_t - 1.96 * se_t
ci_hi <- mu_t + 1.96 * se_t

coverage <- function(idx) {
  if (length(idx) == 0L) return(NA_real_)
  mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
}

resid_t <- mu_t - tr_t
rmse    <- sqrt(mean(resid_t^2, na.rm = TRUE))
r2      <- 1 - sum(resid_t^2, na.rm = TRUE) /
  sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  mean_SE=%.4f",
                rmse, r2, coverage(which(ns_t >= 1L)), mean(se_t, na.rm = TRUE)))

# ------------------------------------------------------------------------------
# 9. Save
# ------------------------------------------------------------------------------

results_ata_kriging_matched <- list(
  run_name  = "ata_kriging_matched",
  tags      = list(resolution    = 330L,
                   method        = "ata_kriging_matern_nu1",
                   cov_model     = "matern_nu1",
                   param_source  = "no_cov_rho_cv",
                   response      = "10m_direct"),
  timestamp = Sys.time(),

  # Parameters
  rho_cv           = rho_cv,
  phi_cv           = phi_cv,
  sigma2e_cv       = sigma2e_cv,
  kappa            = KAPPA,
  kappa2           = KAPPA2,
  sigma2           = SIGMA2,
  sigma2_e         = SIGMA2_E,
  mean_AQA         = mean_AQA,
  scale_factor     = scale_factor,

  # Predictions
  posterior_mean        = mu_t,
  posterior_se          = se_t,
  ci_lower              = ci_lo,
  ci_upper              = ci_hi,
  mu_ok                 = mu_ok,
  rmse                  = rmse,
  r2                    = r2,
  coverage_95_all       = coverage(seq_along(mu_t)),
  coverage_95_obs       = coverage(which(ns_t >= 1L)),
  coverage_95_dense     = coverage(which(ns_t >= 20L)),
  n_soundings_per_pixel = ns_t,
  timing                = timing
)

usethis::use_data(results_ata_kriging_matched, overwrite = TRUE)

# ------------------------------------------------------------------------------
# 10. Plot
# ------------------------------------------------------------------------------

message("== 10. Plotting ==")

mu_sar <- cv$posterior_mean
se_sar <- cv$posterior_se

ns_cut  <- cut(ns_t, breaks = c(-Inf, 0, 5, 20, Inf),
               labels = c("0", "1-5", "6-20", ">20"))
col_pal <- c("0" = "#AAAAAA", "1-5" = "#7BAFD4",
             "6-20" = "#3A7EBF", ">20" = "#1B3F6E")
pt_col  <- col_pal[as.character(ns_cut)]

pdf("fig_ata_vs_sar_matched.pdf", width = 9, height = 4.5)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1.5))

lim_mu <- range(c(mu_sar, mu_t), na.rm = TRUE)
plot(mu_sar, mu_t,
     col = pt_col, pch = 16, cex = 0.5,
     xlim = lim_mu, ylim = lim_mu,
     xlab = "SAR posterior mean (no-cov, CV-tuned)",
     ylab = expression("ATA kriging posterior mean (Mat\u00e9rn "*nu*"=1)"),
     main = sprintf("(a) Posterior means\nkappa=%.4g m-1  sigma2_M=%.4g",
                    KAPPA, SIGMA2))
abline(0, 1, lty = 2, col = "red", lwd = 1.5)
legend("topleft",     bty = "n",
       legend = sprintf("r = %.4f", cor(mu_sar, mu_t, use = "complete.obs")))
legend("bottomright", bty = "n", legend = names(col_pal), fill = col_pal,
       title = "Soundings", cex = 0.8)

lim_se <- range(c(se_sar, se_t), na.rm = TRUE)
plot(se_sar, se_t,
     col = pt_col, pch = 16, cex = 0.5,
     xlim = lim_se, ylim = lim_se,
     xlab = "SAR posterior SE (no-cov, CV-tuned)",
     ylab = expression("ATA kriging posterior SE (Mat\u00e9rn "*nu*"=1)"),
     main = "(b) Posterior SEs")
abline(0, 1, lty = 2, col = "red", lwd = 1.5)
legend("topleft", bty = "n",
       legend = sprintf("r = %.4f", cor(se_sar, se_t, use = "complete.obs")))

dev.off()
message("  Saved fig_ata_vs_sar_matched.pdf")

# ------------------------------------------------------------------------------
# 11. Summary
# ------------------------------------------------------------------------------

cat("\n=== ATA kriging (Matérn nu=1, matched) vs SAR (CV) ===\n")
cat(sprintf("  rho_cv           : %.4f\n",   rho_cv))
cat(sprintf("  KAPPA (m^-1)     : %.6g\n",   KAPPA))
cat(sprintf("  Eff. range       : %.1f m\n", sqrt(8*nu) / KAPPA))
cat(sprintf("  SIGMA2         : %.4g\n",   SIGMA2))
cat(sprintf("  SIGMA2_E         : %.4g\n",   SIGMA2_E))
cat(sprintf("  SNR              : %.2f\n",   SIGMA2 / SIGMA2_E))
cat(sprintf("  SAR  RMSE        : %.4f\n",   cv$rmse))
cat(sprintf("  ATA  RMSE        : %.4f\n",   rmse))
cat(sprintf("  Corr(means)      : %.6f\n",   cor(mu_sar, mu_t,  use = "complete.obs")))
cat(sprintf("  Corr(SEs)        : %.6f\n",   cor(se_sar, se_t,  use = "complete.obs")))
cat(sprintf("  |mean diff|      : %.6f\n",   mean(abs(mu_t  - mu_sar), na.rm = TRUE)))
cat(sprintf("  |SE diff|        : %.6f\n",   mean(abs(se_t  - se_sar), na.rm = TRUE)))

message("\nrun_ata_kriging_matched.R complete.")
