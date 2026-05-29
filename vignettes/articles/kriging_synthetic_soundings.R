# kriging_profile_subset.R
#
# Profile likelihood + LOO-CV over psill for ATA kriging.
# Tests two Matern nu=1 range values:
#   - RANGE_VARIO   : variogram-fitted (768m)
#   - RANGE_MATCHED : SAR-matched for rho=0.99, h=330m (Table 5.1)
#
# Checkpointed: each expensive step saves an RDS and is skipped on rerun.
# Checkpoint files are keyed by N_SUB so subset and full runs don't collide.
#
# To rerun a step, delete its checkpoint file:
#   qmc      -> ckpt_qmc_{N_SUB}.rds
#   K vario  -> ckpt_K0_vario_{N_SUB}.rds
#   K matched-> ckpt_K0_matched_{N_SUB}.rds
#   profile  -> ckpt_profile_{N_SUB}.rds
#   predict  -> ckpt_predict_{N_SUB}.rds

library(spatintegrate)
library(goebel2026)
library(sf)
library(Matrix)
library(future)
library(future.apply)

# ==============================================================================
# 0. Configuration
# ==============================================================================

N_SUB          <- 6850L   # set to smaller number for test run
N_PER_TRIANGLE <- 16L
N_WORKERS      <- 4L
K_CHUNK        <- 20L
T_CHUNK        <- 50L
SEED           <- 42L

H       <- 330
D       <- 2
RHO_SAR <- 0.99

PSILL_VARIO <- 0.00034
RANGE_VARIO <- 768

# Checkpoint file paths (keyed by N_SUB so subset/full don't collide)
ckpt <- function(name) sprintf("ckpt_%s_%d.rds", name, N_SUB)

# ==============================================================================
# 1. SAR <-> Matern nu=1 range matching
# ==============================================================================

kappa2_matched <- (2 * D / H^2) * ((1 - RHO_SAR) / RHO_SAR)
RANGE_MATCHED  <- sqrt(2) / sqrt(kappa2_matched)
kappa_vario    <- sqrt(2) / RANGE_VARIO
rho_vario      <- (2 * D / H^2) / (kappa_vario^2 + 2 * D / H^2)

cat("== SAR <-> Matern nu=1 range matching ==\n")
cat(sprintf("  rho=%.2f  h=%dm  d=%d\n", RHO_SAR, H, D))
cat(sprintf("  SAR-matched range : %.1fm\n", RANGE_MATCHED))
cat(sprintf("  Variogram range   : %.1fm  (equivalent rho~%.4f)\n\n",
            RANGE_VARIO, rho_vario))

# ==============================================================================
# 2. Load data
# ==============================================================================

cat("== Loading data ==\n")
d_shared      <- goebel2026::setup_shared
d_albedo      <- goebel2026::setup_albedo

soundings_sf  <- d_shared$soundings_proj
fine_grid_sf  <- d_shared$fine_grid_buffered
y_full        <- d_albedo$y
sigma2_e      <- d_albedo$sigma_eps^2
y_latent_true <- d_albedo$y_latent_true
target_idx    <- which(fine_grid_sf$n_intersects > 0)
true_var      <- var(y_latent_true[target_idx], na.rm = TRUE)
crs_proj      <- sf::st_crs(soundings_sf)

n_full <- length(y_full)
n_sub  <- min(N_SUB, n_full)

if (n_sub >= n_full) {
  sub_idx <- seq_len(n_full)
} else {
  set.seed(SEED)
  sub_idx <- sort(sample(seq_len(n_full), n_sub, replace = FALSE))
}
y_sub             <- y_full[sub_idx]
sounding_geom_sub <- sf::st_geometry(soundings_sf)[sub_idx]

cat(sprintf("  n_full=%d  n_sub=%d  sigma2_e=%.8f\n", n_full, n_sub, sigma2_e))
cat(sprintf("  True field variance : %.6f\n", true_var))
cat(sprintf("  Variogram psill     : %.6f\n\n", PSILL_VARIO))

# ==============================================================================
# Helper functions (defined before use)
# ==============================================================================

get_poly_pts <- function(poly_sfg, crs, qmc_base) {
  tris <- spatintegrate::triangulate_sf(sf::st_sfc(poly_sfg, crs = crs))
  if (length(tris) == 0L) return(matrix(numeric(0), ncol = 2L))
  do.call(rbind, lapply(seq_along(tris), function(k) {
    spatintegrate::map_unit_square_to_triangle(
      qmc_base, spatintegrate::get_triangle_coords(tris[[k]])
    )
  }))
}

build_K0 <- function(S_all, bb, n, n_pts, range_val, n_workers, k_chunk) {
  n_chunks   <- ceiling(n / k_chunk)
  chunk_list <- lapply(seq_len(n_chunks), function(ci) {
    i_start <- (ci - 1L) * k_chunk + 1L
    i_end   <- min(ci * k_chunk, n)
    i_start:i_end
  })
  future::plan(future::multisession, workers = n_workers)
  on.exit(future::plan(future::sequential), add = TRUE)
  K_chunks <- future.apply::future_lapply(
    chunk_list,
    function(i_idx) {
      n_chunk <- length(i_idx)
      np      <- n_pts
      row_idx <- as.vector(outer(seq_len(np), (i_idx - 1L) * np, "+"))
      S_chunk <- S_all[row_idx, , drop = FALSE]
      aa      <- rowSums(S_chunk^2)
      H2      <- outer(aa, bb, "+") - 2 * tcrossprod(S_chunk, S_all)
      H       <- sqrt(pmax(H2, 0))
      u       <- H / range_val
      CV      <- (1 + u) * exp(-u)
      col_grp  <- rep(seq_len(n), each = np)
      CV_j_avg <- t(rowsum(t(CV), col_grp) / np)
      row_grp  <- rep(seq_len(n_chunk), each = np)
      K_chunk  <- rowsum(CV_j_avg, row_grp) / np
      list(i_idx = i_idx, K_chunk = K_chunk)
    },
    future.seed    = TRUE,
    future.globals = list(
      S_all = S_all, bb = bb, n = n, n_pts = n_pts, range_val = range_val
    )
  )
  K <- matrix(0, nrow = n, ncol = n)
  for (res in K_chunks) K[res$i_idx, ] <- res$K_chunk
  0.5 * (K + t(K))
}

eval_profile_grid <- function(K0, y, psill_grid, sigma2_e, range_label) {
  # Returns lls, loo_mse, loo_cvg for each psill.
  # LOO-CV via influence matrix (Sundararajan & Keerthi 2001):
  #   loo_resid_i = (Sigma^{-1} y)_i / (Sigma^{-1})_ii
  #   loo_var_i   = 1 / (Sigma^{-1})_ii
  # Kinv_ii computed as colSums(CKinv^2) where CKinv = solve(chol(Sigma)).
  n       <- length(y)
  lls     <- numeric(length(psill_grid))
  loo_mse <- numeric(length(psill_grid))
  loo_cvg <- numeric(length(psill_grid))
  cat(sprintf("\n  -- %s --\n", range_label))
  for (k in seq_along(psill_grid)) {
    psill <- psill_grid[k]
    t0    <- proc.time()
    Sigma <- psill * K0
    diag(Sigma) <- diag(Sigma) + sigma2_e + 1e-10 * mean(diag(Sigma))
    CK <- tryCatch(chol(Sigma), error = function(e) {
      diag(Sigma) <<- diag(Sigma) + 1e-4 * mean(diag(Sigma)); chol(Sigma)
    })
    log_det    <- 2 * sum(log(diag(CK)))
    Kinv_y     <- backsolve(CK, forwardsolve(t(CK), y))
    quad       <- as.numeric(crossprod(y, Kinv_y))
    lls[k]     <- -0.5 * (n * log(2 * pi) + log_det + quad)
    # CK is upper triangular (R's chol returns U where Sigma = U'U)
    # U^{-1} = backsolve(CK, I)
    # Sigma^{-1}_ii = ||row_i of U^{-1}||^2 = rowSums(U^{-1}^2)
    CKinv      <- backsolve(CK, diag(n))
    Kinv_ii    <- rowSums(CKinv^2)          # colSums was wrong
    loo_resid  <- Kinv_y / Kinv_ii
    loo_var    <- 1 / Kinv_ii
    loo_mse[k] <- mean(loo_resid^2)
    loo_cvg[k] <- mean(abs(loo_resid) <= 1.96 * sqrt(pmax(loo_var, 0)))
    elapsed <- as.numeric((proc.time() - t0)["elapsed"])
    cat(sprintf("    psill=%.5f  ll=%10.2f  loo_rmse=%.5f  loo_cvg=%.3f  %.1fs\n",
                psill, lls[k], sqrt(loo_mse[k]), loo_cvg[k], elapsed))
  }
  list(lls = lls, loo_mse = loo_mse, loo_cvg = loo_cvg)
}

profile_with_extension <- function(K0, y, psill_grid, sigma2_e,
                                   range_label, max_extensions = 6L) {
  res     <- eval_profile_grid(K0, y, psill_grid, sigma2_e, range_label)
  lls     <- res$lls
  loo_mse <- res$loo_mse
  loo_cvg <- res$loo_cvg
  n_ext   <- 0L
  while (which.max(lls) == length(lls) && n_ext < max_extensions) {
    n_ext      <- n_ext + 1L
    psill_hi   <- max(psill_grid)
    new_psills <- sort(c(psill_hi * 1.5, psill_hi * 2.0, psill_hi * 3.0))
    cat(sprintf("\n  [ext %d] MLE at boundary (%.5f), extending to %.5f\n",
                n_ext, psill_hi, max(new_psills)))
    new_res    <- eval_profile_grid(K0, y, new_psills, sigma2_e,
                                    sprintf("%s [ext %d]", range_label, n_ext))
    psill_grid <- c(psill_grid, new_psills)
    lls        <- c(lls,     new_res$lls)
    loo_mse    <- c(loo_mse, new_res$loo_mse)
    loo_cvg    <- c(loo_cvg, new_res$loo_cvg)
    ord        <- order(psill_grid)
    psill_grid <- psill_grid[ord]
    lls        <- lls[ord]
    loo_mse    <- loo_mse[ord]
    loo_cvg    <- loo_cvg[ord]
  }
  if (which.max(lls) == length(lls))
    cat(sprintf("  WARNING: MLE still at boundary after %d extensions.\n", n_ext))
  best_ll <- which.max(lls)
  best_cv <- which.min(loo_mse)
  cat(sprintf("\n  MLE psill: %.5f  (ll=%.2f)\n",
              psill_grid[best_ll], lls[best_ll]))
  cat(sprintf("  CV  psill: %.5f  (loo_rmse=%.5f  loo_cvg=%.3f)\n",
              psill_grid[best_cv], sqrt(loo_mse[best_cv]), loo_cvg[best_cv]))
  list(psill_grid = psill_grid, lls = lls, loo_mse = loo_mse, loo_cvg = loo_cvg,
       best_ll_idx = best_ll, best_cv_idx = best_cv, n_extensions = n_ext)
}

predict_kriging <- function(K0, psill, sigma2_e, y, S_all, bb,
                            idx_pts, n_sub, n_pts, n_target,
                            target_cents, bb2, tr_t, ns_t,
                            range_val, t_chunk, label) {
  cat(sprintf("\n  -- %s --\n", label))
  Sigma <- psill * K0
  diag(Sigma) <- diag(Sigma) + sigma2_e + 1e-10 * mean(diag(Sigma))
  CK     <- chol(Sigma)
  CK_tri <- as(as(as(CK, "dMatrix"), "triangularMatrix"), "unpackedMatrix")
  ones   <- rep(1, n_sub)
  Kinv_1 <- as.numeric(solve(CK_tri, solve(t(CK_tri), ones)))
  mu_ok  <- sum(Kinv_1 * y) / sum(Kinv_1)
  Kinv_e <- as.numeric(solve(CK_tri, solve(t(CK_tri), y - mu_ok)))
  cat(sprintf("    Kriging global mean: %.4f\n", mu_ok))
  cat(sprintf("    Building k* (%d targets x %d soundings)...\n", n_target, n_sub))
  t0    <- proc.time()
  k_mat <- matrix(0, nrow = n_target, ncol = n_sub)
  for (chunk_idx in seq_len(ceiling(n_target / t_chunk))) {
    j_start  <- (chunk_idx - 1L) * t_chunk + 1L
    j_end    <- min(chunk_idx * t_chunk, n_target)
    j_idx    <- j_start:j_end
    tc_chunk <- target_cents[j_idx, , drop = FALSE]
    H2 <- outer(bb, rowSums(tc_chunk^2), "+") - 2 * tcrossprod(S_all, tc_chunk)
    H  <- sqrt(pmax(H2, 0))
    u  <- H / range_val
    CV <- psill * (1 + u) * exp(-u)
    k_mat[j_idx, ] <- t(rowsum(CV, idx_pts) / n_pts)
  }
  cat(sprintf("    k* done. %.1fs\n", as.numeric((proc.time() - t0)["elapsed"])))
  mu_t    <- mu_ok + as.numeric(k_mat %*% Kinv_e)
  k_t     <- t(k_mat)
  Kinv_kt <- as.matrix(solve(CK_tri, solve(t(CK_tri), k_t)))
  kKinvk  <- colSums(Kinv_kt * k_t)
  kKinv1  <- as.numeric(k_mat %*% Kinv_1)
  ok_corr <- (1 - kKinv1)^2 / sum(Kinv_1)
  var_t   <- pmax(psill - kKinvk + ok_corr, 0)
  se_t    <- sqrt(var_t)
  ci_lo   <- mu_t - 1.96 * se_t
  ci_hi   <- mu_t + 1.96 * se_t
  cov_fn  <- function(ii) {
    if (length(ii) == 0L) return(NA_real_)
    mean(tr_t[ii] >= ci_lo[ii] & tr_t[ii] <= ci_hi[ii], na.rm = TRUE)
  }
  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)
  cat(sprintf("    RMSE=%.5f  R2=%.4f  Cvg(all)=%.3f  Cvg(obs)=%.3f  Cvg(>=20)=%.3f\n",
              rmse, r2, cov_fn(seq_along(mu_t)),
              cov_fn(which(ns_t >= 1L)), cov_fn(which(ns_t >= 20L))))
  list(mu_t=mu_t, se_t=se_t, ci_lo=ci_lo, ci_hi=ci_hi,
       rmse=rmse, r2=r2,
       cov_all=cov_fn(seq_along(mu_t)), cov_obs=cov_fn(which(ns_t >= 1L)),
       cov_dense=cov_fn(which(ns_t >= 20L)))
}

# ==============================================================================
# 3. QMC quadrature points  [checkpoint: ckpt_qmc_{N_SUB}.rds]
# ==============================================================================

if (file.exists(ckpt("qmc"))) {
  cat(sprintf("== QMC: loading from %s ==\n", ckpt("qmc")))
  qmc_ckpt <- readRDS(ckpt("qmc"))
  S_all    <- qmc_ckpt$S_all
  bb       <- qmc_ckpt$bb
  n_pts    <- qmc_ckpt$n_pts
  cat(sprintf("  S_all: %d x %d  n_pts=%d\n", nrow(S_all), ncol(S_all), n_pts))
} else {
  cat("== Building QMC quadrature points ==\n")
  t0       <- proc.time()
  qmc_base <- spatintegrate::generate_qmc_unit_square(N_PER_TRIANGLE)
  sounding_pts <- vector("list", n_sub)
  for (i in seq_len(n_sub)) {
    if (i %% 500L == 0L) cat(sprintf("  QMC: %d / %d\n", i, n_sub))
    sounding_pts[[i]] <- get_poly_pts(sounding_geom_sub[[i]], crs_proj, qmc_base)
  }
  n_pts_vec <- vapply(sounding_pts, nrow, integer(1L))
  n_pts     <- as.integer(median(n_pts_vec))
  S_all     <- do.call(rbind, lapply(sounding_pts, function(pts) {
    if (nrow(pts) >= n_pts) pts[seq_len(n_pts), , drop = FALSE]
    else pts[sample(nrow(pts), n_pts, replace = TRUE), , drop = FALSE]
  }))
  bb <- rowSums(S_all^2)
  cat(sprintf("  Done. %.1fs  n_pts=%d  S_all: %d x %d\n",
              as.numeric((proc.time() - t0)["elapsed"]), n_pts,
              nrow(S_all), ncol(S_all)))
  saveRDS(list(S_all=S_all, bb=bb, n_pts=n_pts, sub_idx=sub_idx), ckpt("qmc"))
  cat(sprintf("  Saved to %s\n\n", ckpt("qmc")))
}

# ==============================================================================
# 4. K0 formation  [checkpoints: ckpt_K0_vario_{N_SUB}.rds, ckpt_K0_matched_{N_SUB}.rds]
# ==============================================================================

if (file.exists(ckpt("K0_vario"))) {
  cat(sprintf("== K0 vario: loading from %s ==\n", ckpt("K0_vario")))
  K0_vario <- readRDS(ckpt("K0_vario"))
  t_vario  <- NA
} else {
  cat("== Forming K0: variogram range (768m) ==\n")
  t0       <- proc.time()
  K0_vario <- build_K0(S_all, bb, n_sub, n_pts,
                       range_val=RANGE_VARIO, n_workers=N_WORKERS, k_chunk=K_CHUNK)
  t_vario  <- as.numeric((proc.time() - t0)["elapsed"])
  cat(sprintf("  Done. %.1fs  diag: min=%.4f  mean=%.4f  max=%.4f\n",
              t_vario, min(diag(K0_vario)), mean(diag(K0_vario)), max(diag(K0_vario))))
  saveRDS(K0_vario, ckpt("K0_vario"))
  cat(sprintf("  Saved to %s\n", ckpt("K0_vario")))
}

if (file.exists(ckpt("K0_matched"))) {
  cat(sprintf("== K0 matched: loading from %s ==\n", ckpt("K0_matched")))
  K0_matched <- readRDS(ckpt("K0_matched"))
  t_matched  <- NA
} else {
  cat(sprintf("== Forming K0: SAR-matched range (%.0fm) ==\n", RANGE_MATCHED))
  t0         <- proc.time()
  K0_matched <- build_K0(S_all, bb, n_sub, n_pts,
                         range_val=RANGE_MATCHED, n_workers=N_WORKERS, k_chunk=K_CHUNK)
  t_matched  <- as.numeric((proc.time() - t0)["elapsed"])
  cat(sprintf("  Done. %.1fs  diag: min=%.4f  mean=%.4f  max=%.4f\n",
              t_matched, min(diag(K0_matched)), mean(diag(K0_matched)), max(diag(K0_matched))))
  saveRDS(K0_matched, ckpt("K0_matched"))
  cat(sprintf("  Saved to %s\n", ckpt("K0_matched")))
}

# ==============================================================================
# 5. Profile likelihood + LOO-CV  [checkpoint: ckpt_profile_{N_SUB}.rds]
# ==============================================================================

psill_grid_init <- sort(unique(c(
  0.00010, 0.00020, PSILL_VARIO, 0.00050,
  0.00100, true_var, 0.00200, 0.00500
)))

if (file.exists(ckpt("profile"))) {
  cat(sprintf("== Profile: loading from %s ==\n", ckpt("profile")))
  prof_ckpt   <- readRDS(ckpt("profile"))
  res_vario   <- prof_ckpt$res_vario
  res_matched <- prof_ckpt$res_matched
} else {
  cat(sprintf("\n== Profile likelihood + LOO-CV ==\n"))
  cat(sprintf("  Initial psill grid: %s\n",
              paste(sprintf("%.5f", psill_grid_init), collapse=", ")))
  cat(sprintf("  sigma2_e = %.8f (true injected, fixed)\n\n", sigma2_e))
  res_vario <- profile_with_extension(
    K0_vario, y_sub, psill_grid_init, sigma2_e,
    sprintf("Range=%.0fm (rho~%.4f)", RANGE_VARIO, rho_vario)
  )
  res_matched <- profile_with_extension(
    K0_matched, y_sub, psill_grid_init, sigma2_e,
    sprintf("Range=%.0fm (SAR-matched rho=%.2f)", RANGE_MATCHED, RHO_SAR)
  )
  saveRDS(list(res_vario=res_vario, res_matched=res_matched), ckpt("profile"))
  cat(sprintf("\n  Saved to %s\n", ckpt("profile")))
}

psill_grid_v <- res_vario$psill_grid
psill_grid_m <- res_matched$psill_grid
ll_vario     <- res_vario$lls
ll_matched   <- res_matched$lls
best_vi      <- res_vario$best_ll_idx
best_mi      <- res_matched$best_ll_idx
best_vi_cv   <- res_vario$best_cv_idx
best_mi_cv   <- res_matched$best_cv_idx

# ==============================================================================
# 6. Profile results table
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat(sprintf("  Profile likelihood + LOO-CV  [n=%d  seed=%d]\n", n_sub, SEED))
cat(sprintf("  Covariance: Matern nu=1  sigma2_e=%.8f (fixed)\n", sigma2_e))
cat("================================================================\n\n")

for (rng in c("vario", "matched")) {
  if (rng == "vario") {
    pg <- psill_grid_v; ll <- ll_vario
    lm <- res_vario$loo_mse; lc <- res_vario$loo_cvg
    bi <- best_vi; bc <- best_vi_cv
    rl <- sprintf("Range=%.0fm (rho~%.4f)", RANGE_VARIO, rho_vario)
  } else {
    pg <- psill_grid_m; ll <- ll_matched
    lm <- res_matched$loo_mse; lc <- res_matched$loo_cvg
    bi <- best_mi; bc <- best_mi_cv
    rl <- sprintf("Range=%.0fm (SAR rho=%.2f)", RANGE_MATCHED, RHO_SAR)
  }
  cat(sprintf("  -- %s --\n", rl))
  cat(sprintf("  %-12s  %-12s  %-12s  %-10s\n",
              "psill", "log-lik", "loo_rmse", "loo_cvg"))
  cat(strrep("-", 52), "\n")
  for (k in seq_along(pg)) {
    fl <- if (k == bi) "L" else " "
    fc <- if (k == bc) "C" else " "
    cat(sprintf("  %-12.5f  %-12.2f  %-12.5f  %-10.3f  %s%s\n",
                pg[k], ll[k], sqrt(lm[k]), lc[k], fl, fc))
  }
  cat("  (L=MLE  C=CV)\n\n")
}

cat(sprintf("  %-28s  %-12s  %-12s\n", "",
            sprintf("%.0fm", RANGE_VARIO), sprintf("%.0fm", RANGE_MATCHED)))
cat(strrep("-", 56), "\n")
cat(sprintf("  %-28s  %-12.5f  %-12.5f\n", "MLE psill",
            psill_grid_v[best_vi], psill_grid_m[best_mi]))
cat(sprintf("  %-28s  %-12.5f  %-12.5f\n", "CV  psill",
            psill_grid_v[best_vi_cv], psill_grid_m[best_mi_cv]))
cat(sprintf("  %-28s  %-12.5f  %-12.5f\n", "True field variance",
            true_var, true_var))
cat(sprintf("  %-28s  %-12.5f  %-12.5f\n", "Variogram psill",
            PSILL_VARIO, PSILL_VARIO))
cat("================================================================\n")

# ==============================================================================
# 7. Predictions  [checkpoint: ckpt_predict_{N_SUB}.rds]
# ==============================================================================

target_sf    <- fine_grid_sf[target_idx, ]
n_target     <- length(target_idx)
tr_t         <- y_latent_true[target_idx]
ns_t         <- as.integer(Matrix::colSums(d_shared$A_flat[sub_idx, target_idx] > 0))
target_cents <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(target_sf)))
bb2          <- rowSums(target_cents^2)
idx_pts      <- rep(seq_len(n_sub), each = n_pts)

configs <- list(
  list(K0=K0_vario,   psill=psill_grid_v[best_vi],    rv=RANGE_VARIO,   label="768m  MLE"),
  list(K0=K0_vario,   psill=psill_grid_v[best_vi_cv], rv=RANGE_VARIO,   label="768m  CV"),
  list(K0=K0_matched, psill=psill_grid_m[best_mi],    rv=RANGE_MATCHED,
       label=sprintf("%.0fm MLE", RANGE_MATCHED)),
  list(K0=K0_matched, psill=psill_grid_m[best_mi_cv], rv=RANGE_MATCHED,
       label=sprintf("%.0fm CV",  RANGE_MATCHED))
)

if (file.exists(ckpt("predict"))) {
  cat(sprintf("\n== Predictions: loading from %s ==\n", ckpt("predict")))
  pred_results <- readRDS(ckpt("predict"))
} else {
  cat("\n== Predictions at MLE and CV psill ==\n")
  pred_results <- lapply(configs, function(cfg) {
    predict_kriging(
      K0=cfg$K0, psill=cfg$psill, sigma2_e=sigma2_e,
      y=y_sub, S_all=S_all, bb=bb,
      idx_pts=idx_pts, n_sub=n_sub, n_pts=n_pts,
      n_target=n_target, target_cents=target_cents, bb2=bb2,
      tr_t=tr_t, ns_t=ns_t, range_val=cfg$rv,
      t_chunk=T_CHUNK,
      label=sprintf("%s  psill=%.5f", cfg$label, cfg$psill)
    )
  })
  names(pred_results) <- sapply(configs, `[[`, "label")
  saveRDS(pred_results, ckpt("predict"))
  cat(sprintf("  Saved to %s\n", ckpt("predict")))
}

# ==============================================================================
# 8. Prediction results table
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat(sprintf("  Prediction metrics  [n=%d]\n", n_sub))
cat(sprintf("  True field variance: %.5f\n", true_var))
cat("================================================================\n")
cat(sprintf("  %-20s  %-9s  %-8s  %-7s  %-10s  %-10s\n",
            "Config", "psill", "RMSE", "R2", "Cvg(obs)", "Cvg(>=20)"))
cat(strrep("-", 72), "\n")
for (i in seq_along(configs)) {
  nm <- configs[[i]]$label
  p  <- pred_results[[nm]]
  cat(sprintf("  %-20s  %-9.5f  %-8.4f  %-7.4f  %-10.3f  %-10.3f\n",
              nm, configs[[i]]$psill, p$rmse, p$r2, p$cov_obs, p$cov_dense))
}
cat("================================================================\n")

# ==============================================================================
# 9. Save final results
# ==============================================================================

final <- list(
  n_sub           = n_sub,
  sub_idx         = sub_idx,
  seed            = SEED,
  psill_grid_v    = psill_grid_v,
  psill_grid_m    = psill_grid_m,
  ll_vario        = ll_vario,
  ll_matched      = ll_matched,
  loo_mse_vario   = res_vario$loo_mse,
  loo_mse_matched = res_matched$loo_mse,
  loo_cvg_vario   = res_vario$loo_cvg,
  loo_cvg_matched = res_matched$loo_cvg,
  best_vario_mle  = psill_grid_v[best_vi],
  best_vario_cv   = psill_grid_v[best_vi_cv],
  best_matched_mle = psill_grid_m[best_mi],
  best_matched_cv  = psill_grid_m[best_mi_cv],
  true_var        = true_var,
  sigma2_e        = sigma2_e,
  range_vario     = RANGE_VARIO,
  range_matched   = RANGE_MATCHED,
  rho_vario       = rho_vario,
  n_pts           = n_pts,
  pred_results    = pred_results
)

saveRDS(final, "kriging_profile_results.rds")
cat("\nFinal results saved to kriging_profile_results.rds\n")
cat(sprintf("Checkpoints: %s\n",
            paste(sapply(c("qmc","K0_vario","K0_matched","profile","predict"),
                         ckpt), collapse=", ")))
