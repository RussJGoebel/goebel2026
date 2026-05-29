# check_centroid_variogram.R
#
# Fits a point-support variogram using sounding centroids, then runs
# ATA kriging using those parameters with the full areal covariance K.
#
# The centroid variogram uses ALL n*(n-1)/2 pairs -- no areal integration,
# just pairwise distances between centroids. Fast to compute.
#
# Key question: does the centroid variogram give different (psill, range)
# than the areal variogram? Specifically, does it give psill closer to
# the true field variance (~0.001)?
#
# If yes: use those parameters for ATA kriging and check coverage.
# If no: centroid approximation gives same answer as areal variogram.

library(Matrix)
library(goebel2026)
library(sf)
library(spatintegrate)
library(future)
library(future.apply)

# ------------------------------------------------------------------------------
# 0. Load data
# ------------------------------------------------------------------------------

message("== Loading data ==")
d_shared      <- goebel2026::setup_shared
d_albedo      <- goebel2026::setup_albedo
y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true
fine_grid_sf  <- d_shared$fine_grid_buffered
soundings_sf  <- d_shared$soundings_proj
target_idx    <- which(fine_grid_sf$n_intersects > 0)
target_sf     <- fine_grid_sf[target_idx, ]

n          <- length(y_alb)
n_target   <- length(target_idx)
SIGMA2_E   <- d_albedo$sigma_eps^2
TRUE_VAR   <- var(y_latent_true[target_idx], na.rm=TRUE)

message(sprintf("  n=%d  sigma2_e=%.8f  true_field_var=%.6f",
                n, SIGMA2_E, TRUE_VAR))

cent_mat <- sf::st_coordinates(
  sf::st_centroid(sf::st_geometry(soundings_sf))
)

# ------------------------------------------------------------------------------
# 1. Centroid-based empirical variogram -- ALL pairs
#
# Uses centroid distance as proxy for spatial separation.
# No areal integration -- pure point-support variogram.
# Sill constrained to var(y_alb) ≈ 0.00034 at long lags.
# ------------------------------------------------------------------------------

message("\n== Step 1: Centroid-based empirical variogram (all pairs) ==")
t0 <- proc.time()

N_BINS   <- 20L
MAX_DIST <- 15000
CHUNK    <- 500L

breaks   <- seq(0, MAX_DIST, length.out = N_BINS + 1L)
bin_mids <- 0.5 * (breaks[-1] + breaks[-length(breaks)])
emp_sum  <- numeric(N_BINS)
emp_cnt  <- integer(N_BINS)

for (i_start in seq(1, n, by = CHUNK)) {
  i_end <- min(i_start + CHUNK - 1L, n)
  i_idx <- i_start:i_end
  dx <- outer(cent_mat[i_idx,1], cent_mat[,1], "-")
  dy <- outer(cent_mat[i_idx,2], cent_mat[,2], "-")
  h  <- sqrt(dx^2 + dy^2)
  gv <- outer(y_alb[i_idx], y_alb, function(a,b) 0.5*(a-b)^2)
  for (k in seq_along(i_idx)) {
    ig <- i_idx[k]
    h[k,  seq_len(ig)] <- NA_real_
    gv[k, seq_len(ig)] <- NA_real_
  }
  bin_idx <- findInterval(h, breaks, rightmost.closed=TRUE)
  valid   <- !is.na(h) & bin_idx >= 1L & bin_idx <= N_BINS
  for (b in seq_len(N_BINS)) {
    sel <- valid & bin_idx == b
    if (any(sel)) {
      emp_sum[b] <- emp_sum[b] + sum(gv[sel], na.rm=TRUE)
      emp_cnt[b] <- emp_cnt[b] + sum(sel)
    }
  }
}

emp_gamma_cent <- ifelse(emp_cnt > 0L, emp_sum / emp_cnt, NA_real_)

timing_vgm <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing_vgm))
message(sprintf("  Bin counts range: [%d, %d]", min(emp_cnt), max(emp_cnt)))
message(sprintf("  Short-lag gamma: %.6f  Long-lag gamma: %.6f",
                emp_gamma_cent[1], emp_gamma_cent[N_BINS]))
message(sprintf("  Sill: %.6f  (obs var=%.6f  true field var=%.6f)",
                max(emp_gamma_cent, na.rm=TRUE), var(y_alb), TRUE_VAR))

# ------------------------------------------------------------------------------
# 2. Fit point-support variogram to centroid variogram
#
# Simple exponential WLS -- no areal integration needed.
# Fit nugget + psill + range (3 parameters) since centroid variogram
# can see the nugget at short lags.
# ------------------------------------------------------------------------------

message("\n== Step 2: Fitting point-support variogram (exponential, nugget free) ==")
t0 <- proc.time()

valid_bins <- !is.na(emp_gamma_cent) & emp_cnt > 0L

# With nugget
wls_3param <- function(log_theta) {
  nugget <- exp(log_theta[1L])
  psill  <- exp(log_theta[2L])
  range  <- exp(log_theta[3L])
  theo   <- nugget + psill * (1 - exp(-bin_mids / range))
  sum(emp_cnt[valid_bins] * (emp_gamma_cent[valid_bins] - theo[valid_bins])^2)
}

# Without nugget (sigma2_e fixed)
wls_2param <- function(log_theta) {
  psill <- exp(log_theta[1L])
  range <- exp(log_theta[2L])
  theo  <- SIGMA2_E + psill * (1 - exp(-bin_mids / range))
  sum(emp_cnt[valid_bins] * (emp_gamma_cent[valid_bins] - theo[valid_bins])^2)
}

# Fit with nugget free
opt_3 <- optim(
  par     = c(log(0.0001), log(0.0003), log(3000)),
  fn      = wls_3param,
  method  = "Nelder-Mead",
  control = list(maxit=500L, reltol=1e-6)
)
CENT_NUGGET <- exp(opt_3$par[1L])
CENT_PSILL  <- exp(opt_3$par[2L])
CENT_RANGE  <- exp(opt_3$par[3L])

# Fit with sigma2_e fixed
opt_2 <- optim(
  par     = c(log(0.0003), log(3000)),
  fn      = wls_2param,
  method  = "Nelder-Mead",
  control = list(maxit=500L, reltol=1e-6)
)
CENT_PSILL_FIXED <- exp(opt_2$par[1L])
CENT_RANGE_FIXED <- exp(opt_2$par[2L])

timing_fit <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Elapsed: %.1f sec", timing_fit))

cat(sprintf("\n=== Centroid variogram fit results ===\n"))
cat(sprintf("  Free nugget:   nugget=%.6f  psill=%.6f  range=%.0fm\n",
            CENT_NUGGET, CENT_PSILL, CENT_RANGE))
cat(sprintf("  Fixed sigma2e: nugget=%.8f  psill=%.6f  range=%.0fm\n",
            SIGMA2_E, CENT_PSILL_FIXED, CENT_RANGE_FIXED))
cat(sprintf("\n  Compare:\n"))
cat(sprintf("    Areal WLS (exp):      psill=%.6f  range=%.0fm\n", 0.00038, 1644))
cat(sprintf("    Areal WLS (Matern1):  psill=%.6f  range=%.0fm\n", 0.00034, 768))
cat(sprintf("    MLE Matern1 peak:     psill=%.6f\n", 0.001))
cat(sprintf("    True field variance:  psill=%.6f\n", TRUE_VAR))

# ------------------------------------------------------------------------------
# 3. Load saved K and rescale to centroid-fitted parameters
#
# K_m1 was built with Matern nu=1, psill=M1_PSILL=0.00034, range=768m.
# For centroid parameters we need to rebuild K unless range is similar.
#
# Strategy: if CENT_RANGE_FIXED is close to 768m, rescale K_m1.
# Otherwise note that a new K build would be needed.
# ------------------------------------------------------------------------------

message("\n== Step 3: Coverage check with centroid parameters ==")

if (file.exists("ata_K_matern1.rds")) {
  message("  Loading saved Matern nu=1 K...")
  K_m1     <- readRDS("ata_K_matern1.rds")
  M1_PSILL <- 0.00034
  M1_RANGE <- 768

  # Check if range is close enough to rescale
  range_ratio <- CENT_RANGE_FIXED / M1_RANGE
  message(sprintf("  Centroid range=%.0fm  vs  Matern K range=%.0fm  (ratio=%.2f)",
                  CENT_RANGE_FIXED, M1_RANGE, range_ratio))

  if (abs(range_ratio - 1) < 0.3) {
    message("  Range close enough -- rescaling K by psill ratio")
    K_cent <- K_m1 * (CENT_PSILL_FIXED / M1_PSILL)
  } else {
    message("  Range too different -- using K_m1 as approximation")
    message("  (For exact result would need to rebuild K with new range)")
    K_cent <- K_m1 * (CENT_PSILL_FIXED / M1_PSILL)
  }

  # Add sigma2_e and Cholesky
  K_obs <- K_cent
  diag(K_obs) <- diag(K_obs) + SIGMA2_E + 1e-10 * mean(diag(K_obs))

  message("  Computing Cholesky...")
  t0 <- proc.time()
  CK <- tryCatch(chol(K_obs),
                 error = function(e) {
                   diag(K_obs) <<- diag(K_obs) + 1e-4*mean(diag(K_obs))
                   chol(K_obs)
                 })
  message(sprintf("  Cholesky: %.1f sec", as.numeric((proc.time()-t0)["elapsed"])))

  CKt    <- as(CK, "dtrMatrix")
  ones   <- rep(1, n)
  Kinv_1 <- as.numeric(solve(CKt, solve(t(CKt), ones)))
  mu_ok  <- sum(Kinv_1 * y_alb) / sum(Kinv_1)
  e      <- y_alb - mu_ok
  Kinv_e <- as.numeric(solve(CKt, solve(t(CKt), e)))

  # Prediction covariance vectors
  message("  Computing prediction covariance...")
  target_cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(target_sf)))

  # Load QMC points if available
  if (!exists("S_all") || !exists("bb") || !exists("n_pts")) {
    message("  Rebuilding QMC points...")
    qmc_base      <- spatintegrate::generate_qmc_unit_square(16L)
    sounding_geom <- sf::st_geometry(soundings_sf)
    crs_proj      <- sf::st_crs(soundings_sf)
    get_poly_pts  <- function(poly_sfg, crs) {
      tris <- spatintegrate::triangulate_sf(sf::st_sfc(poly_sfg, crs=crs))
      if (length(tris)==0L) return(matrix(numeric(0), ncol=2L))
      do.call(rbind, lapply(seq_along(tris), function(k)
        spatintegrate::map_unit_square_to_triangle(
          qmc_base, spatintegrate::get_triangle_coords(tris[[k]]))))
    }
    sounding_pts <- vector("list", n)
    for (i in seq_len(n)) {
      if (i %% 1000L == 0L) message(sprintf("  QMC: %d/%d", i, n))
      sounding_pts[[i]] <- get_poly_pts(sounding_geom[[i]], crs_proj)
    }
    n_pts <- as.integer(median(vapply(sounding_pts, nrow, integer(1L))))
    S_all <- do.call(rbind, sounding_pts)
    bb    <- rowSums(S_all^2)
  }

  # Use Matern nu=1 covariance for prediction vectors
  cov_m1_cent <- function(h) {
    x <- h / CENT_RANGE_FIXED
    out <- CENT_PSILL_FIXED * x * besselK(x + 1e-10, nu=1)
    out[h==0] <- CENT_PSILL_FIXED
    out
  }

  T_CHUNK  <- 50L
  bb2      <- rowSums(target_cents^2)
  k_mat    <- matrix(0, nrow=n_target, ncol=n)
  idx_grp  <- rep(seq_len(n), each=n_pts)

  for (chunk_idx in seq_len(ceiling(n_target/T_CHUNK))) {
    j_start <- (chunk_idx-1L)*T_CHUNK+1L
    j_end   <- min(chunk_idx*T_CHUNK, n_target)
    j_idx   <- j_start:j_end
    tc      <- target_cents[j_idx,,drop=FALSE]
    ab      <- tcrossprod(S_all, tc)
    H2      <- outer(bb, rowSums(tc^2), "+") - 2*ab
    H       <- sqrt(pmax(H2, 0))
    CV      <- cov_m1_cent(H)
    k_mat[j_idx,] <- t(rowsum(CV, idx_grp) / n_pts)
    if (chunk_idx %% 20L==0L || chunk_idx==ceiling(n_target/T_CHUNK))
      message(sprintf("  pred cov chunk %d/%d", chunk_idx, ceiling(n_target/T_CHUNK)))
  }

  # Predictions
  mu_t    <- mu_ok + as.numeric(k_mat %*% Kinv_e)
  C_self  <- CENT_PSILL_FIXED

  k_t     <- t(k_mat)
  Kinv_kt <- as.matrix(solve(CKt, solve(t(CKt), k_t)))
  kKinvk  <- colSums(Kinv_kt * k_t)
  kKinv1  <- as.numeric(k_mat %*% Kinv_1)
  ok_corr <- (1 - kKinv1)^2 / sum(Kinv_1)
  var_t   <- pmax(C_self - kKinvk + ok_corr, 0)
  se_t    <- sqrt(var_t)

  # Evaluate
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- as.integer(Matrix::colSums(d_shared$A_flat[,target_idx] > 0))
  ci_lo <- mu_t - 1.96*se_t
  ci_hi <- mu_t + 1.96*se_t

  coverage <- function(idx) {
    if (length(idx)==0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm=TRUE)
  }

  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm=TRUE))
  r2    <- 1 - sum(resid^2,na.rm=TRUE)/sum((tr_t-mean(tr_t,na.rm=TRUE))^2,na.rm=TRUE)

  cat(sprintf("\n=== Centroid variogram -> ATA kriging results ===\n"))
  cat(sprintf("  Parameters: nugget=0  psill=%.6f  range=%.0fm  sigma2_e=%.8f\n",
              CENT_PSILL_FIXED, CENT_RANGE_FIXED, SIGMA2_E))
  cat(sprintf("  RMSE=%.4f  R2=%.4f\n", rmse, r2))
  cat(sprintf("  coverage_obs=%.3f  mean_SE=%.4f\n",
              coverage(which(ns_t>=1L)), mean(se_t,na.rm=TRUE)))
  cat(sprintf("  C_self=%.6f  mean(kKinvk)=%.6f  mean(var_t)=%.6f\n",
              C_self, mean(kKinvk), mean(var_t)))

  cat(sprintf("\n=== Comparison ===\n"))
  cat(sprintf("  %-40s  %6s  %6s  %8s\n", "Method", "RMSE", "R2", "Coverage"))
  cat(sprintf("  %-40s  %6.4f  %6.4f  %8.3f\n",
              "Exp WLS (fitted)", 0.0194, 0.6249, 0.728))
  cat(sprintf("  %-40s  %6.4f  %6.4f  %8.3f\n",
              "Matern1 WLS (fitted)", 0.0208, 0.5670, 0.739))
  cat(sprintf("  %-40s  %6.4f  %6.4f  %8.3f\n",
              "Oracle (true params)", 0.0196, 0.6153, 0.938))
  cat(sprintf("  %-40s  %6.4f  %6.4f  %8.3f\n",
              "Centroid variogram -> ATA", rmse, r2, coverage(which(ns_t>=1L))))
  cat(sprintf("  %-40s  %6.4f  %6.4f  %8.3f\n",
              "Your SAR model", 0.0208, 0.569, 0.970))

} else {
  message("  ata_K_matern1.rds not found -- cannot run coverage check")
  message("  Run check_matern1_kriging.R first to generate the K matrix")
}
