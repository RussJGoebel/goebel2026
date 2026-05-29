# kriging_profile_subset.R
#
# Self-contained profile likelihood comparison for ATA kriging.
# Tests two Matern nu=1 range values:
#
#   - 768m  : variogram-fitted range (from existing work)
#   - 2333m : SAR-matched range for rho=0.99, h=330m, d=2
#
# And a grid of psill values at each range, with sigma2_e fixed
# at the true injected value (as in profile_likelihood_psill.R).
#
# Uses n=500 random subsample. K formation is the only real cost;
# Cholesky on 500x500 is negligible.
#
# Covariance: Matern nu=1 (Whittle kernel):
#   C(h) = psill * (1 + h/range) * exp(-h/range)
# This is the correct Matern nu=1, NOT the exponential.
# The SAR <-> Matern matching table assumes this family.
#
# Expected runtime: 2 * K_formation_time + negligible Cholesky
# K formation scales as (n_sub/n_full)^2, so ~3s per K at n=500
# if the full n=6850 K takes 20 min with same parallelism.

library(Matrix)
library(future)
library(future.apply)
library(goebel2026)

# ==============================================================================
# 0. Configuration -- edit these
# ==============================================================================

N_SUB     <- 500     # subsample size; set to Inf or n_full to use all
N_WORKERS <- parallel::detectCores() - 1L
K_CHUNK   <- 50L     # soundings per parallel chunk
N_PTS     <- 50L     # quadrature points per sounding (match your main script)
SEED      <- 42L

# SAR/Matern matching parameters
H         <- 330     # latent pixel size in meters
D         <- 2       # spatial dimension
RHO_SAR   <- 0.99   # SAR rho we are matching to

# Psill grid anchor from existing work
M1_PSILL_VARIO <- 0.00034

# ==============================================================================
# 1. SAR <-> Matern range matching
# ==============================================================================

# From Table 5.1:
#   kappa^2 = (2d / h^2) * (1 - rho) / rho
#   range   = sqrt(2 * nu) / kappa     [nu=1]
kappa2_matched <- (2 * D / H^2) * ((1 - RHO_SAR) / RHO_SAR)
kappa_matched  <- sqrt(kappa2_matched)
RANGE_MATCHED  <- sqrt(2 * 1) / kappa_matched

# What rho does the variogram range correspond to?
kappa_vario  <- sqrt(2 * 1) / 768
kappa2_vario <- kappa_vario^2
rho_vario    <- (2 * D / H^2) / (kappa2_vario + 2 * D / H^2)

RANGE_VARIO <- 768

cat("== SAR <-> Matern range matching ==\n")
cat(sprintf("  rho=%.2f  h=%dm  d=%d\n", RHO_SAR, H, D))
cat(sprintf("  SAR-matched range : %.1fm\n", RANGE_MATCHED))
cat(sprintf("  Variogram range   : %.1fm  (equivalent rho ~ %.4f)\n\n",
            RANGE_VARIO, rho_vario))

# ==============================================================================
# 2. Load data
# ==============================================================================

cat("== Loading data ==\n")
d_shared      <- goebel2026::setup_shared
d_albedo      <- goebel2026::setup_albedo

y_full        <- d_albedo$y
sigma2_e      <- d_albedo$sigma_eps^2
y_latent_true <- d_albedo$y_latent_true
fine_grid_sf  <- d_shared$fine_grid_buffered
target_idx    <- which(fine_grid_sf$n_intersects > 0)
true_var      <- var(y_latent_true[target_idx], na.rm = TRUE)

# Quadrature points: S_all is (n * n_pts) x 2, stacked by sounding.
# Sounding i occupies rows ((i-1)*n_pts + 1) : (i*n_pts).
# Adjust the name below if yours differs (e.g. d_shared$quad_pts).
S_all_full <- d_shared$S_all
n_full     <- length(y_full)
n_sub      <- min(N_SUB, n_full)

set.seed(SEED)
sub_idx <- sort(sample(seq_len(n_full), n_sub, replace = FALSE))
y_sub   <- y_full[sub_idx]

# Subset S_all to keep only rows for the sampled soundings
row_idx_sub <- as.vector(
  outer(seq_len(N_PTS), (sub_idx - 1L) * N_PTS, "+")
)
S_all_sub <- S_all_full[row_idx_sub, , drop = FALSE]
bb_sub    <- rowSums(S_all_sub^2)

cat(sprintf("  n_full=%d  n_sub=%d  n_pts=%d  sigma2_e=%.8f\n",
            n_full, n_sub, N_PTS, sigma2_e))
cat(sprintf("  True field variance : %.6f\n", true_var))
cat(sprintf("  Variogram psill     : %.6f\n\n", M1_PSILL_VARIO))

# ==============================================================================
# 3. Matern nu=1 covariance (unit psill)
# ==============================================================================

# Matern nu=1 (Whittle): C(h) = (1 + h/r) * exp(-h/r)
matern1_unit <- function(h, range) {
  u <- h / range
  (1 + u) * exp(-u)
}

# ==============================================================================
# 4. ATA K formation (unit psill, fixed range)
# ==============================================================================

build_K0 <- function(S_all, bb, n, n_pts, range, n_workers, k_chunk) {
  # Returns n x n unit-psill ATA Matern nu=1 covariance matrix.
  # Parallelizes over chunks of rows (soundings).

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

      # Quadrature rows for soundings in this chunk
      row_idx <- as.vector(outer(seq_len(n_pts), (i_idx - 1L) * n_pts, "+"))
      S_chunk <- S_all[row_idx, , drop = FALSE]
      aa      <- rowSums(S_chunk^2)

      # Pairwise distances between chunk quadrature points and all points
      H2  <- outer(aa, bb, "+") - 2 * tcrossprod(S_chunk, S_all)
      H   <- sqrt(pmax(H2, 0))
      CV  <- matern1_unit(H, range)   # (n_chunk*n_pts) x (n*n_pts)

      # Average over j quadrature points (columns)
      col_grp  <- rep(seq_len(n), each = n_pts)
      CV_j_avg <- t(rowsum(t(CV), col_grp) / n_pts)   # (n_chunk*n_pts) x n

      # Average over i quadrature points (rows)
      row_grp <- rep(seq_len(n_chunk), each = n_pts)
      K_chunk <- rowsum(CV_j_avg, row_grp) / n_pts    # n_chunk x n

      list(i_idx = i_idx, K_chunk = K_chunk)
    },
    future.seed    = TRUE,
    future.globals = list(
      S_all        = S_all,
      bb           = bb,
      n            = n,
      n_pts        = n_pts,
      matern1_unit = matern1_unit,
      range        = range
    )
  )

  K <- matrix(0, nrow = n, ncol = n)
  for (res in K_chunks) K[res$i_idx, ] <- res$K_chunk
  0.5 * (K + t(K))   # symmetrize for numerical safety
}

# ==============================================================================
# 5. Profile likelihood over psill grid (fixed sigma2_e)
# ==============================================================================

eval_profile_grid <- function(K0, y, psill_grid, sigma2_e, range_label) {
  n   <- length(y)
  lls <- numeric(length(psill_grid))
  cat(sprintf("\n  -- %s --\n", range_label))

  for (k in seq_along(psill_grid)) {
    psill <- psill_grid[k]
    t0    <- proc.time()

    Sigma <- psill * K0
    diag(Sigma) <- diag(Sigma) + sigma2_e + 1e-10 * mean(diag(Sigma))

    CK <- tryCatch(
      chol(Sigma),
      error = function(e) {
        message(sprintf("    chol failed at psill=%.5f, adding jitter", psill))
        diag(Sigma) <<- diag(Sigma) + 1e-4 * mean(diag(Sigma))
        chol(Sigma)
      }
    )

    log_det <- 2 * sum(log(diag(CK)))
    Kinv_y  <- backsolve(CK, forwardsolve(t(CK), y))
    quad    <- as.numeric(crossprod(y, Kinv_y))
    ll      <- -0.5 * (n * log(2 * pi) + log_det + quad)

    elapsed <- as.numeric((proc.time() - t0)["elapsed"])
    cat(sprintf("    psill=%.5f  ll=%10.2f  chol=%.3fs\n", psill, ll, elapsed))
    lls[k] <- ll
  }
  lls
}

# ==============================================================================
# 6. Psill grid
# ==============================================================================

psill_grid <- sort(unique(c(
  0.00010,
  0.00020,
  M1_PSILL_VARIO,
  0.00050,
  0.00100,
  true_var,
  0.00200,
  0.00500
)))

cat(sprintf("Psill grid: %s\n",
            paste(sprintf("%.5f", psill_grid), collapse = ", ")))

# ==============================================================================
# 7. Form K0 for each range
# ==============================================================================

cat("\n== Forming K0: variogram range (768m) ==\n")
t0       <- proc.time()
K0_vario <- build_K0(S_all_sub, bb_sub, n_sub, N_PTS,
                     range     = RANGE_VARIO,
                     n_workers = N_WORKERS,
                     k_chunk   = K_CHUNK)
t_vario  <- as.numeric((proc.time() - t0)["elapsed"])
cat(sprintf("  Done. %.1fs  |  diag: min=%.4f  max=%.4f\n",
            t_vario, min(diag(K0_vario)), max(diag(K0_vario))))

cat(sprintf("\n== Forming K0: SAR-matched range (%.0fm) ==\n", RANGE_MATCHED))
t0         <- proc.time()
K0_matched <- build_K0(S_all_sub, bb_sub, n_sub, N_PTS,
                       range     = RANGE_MATCHED,
                       n_workers = N_WORKERS,
                       k_chunk   = K_CHUNK)
t_matched  <- as.numeric((proc.time() - t0)["elapsed"])
cat(sprintf("  Done. %.1fs  |  diag: min=%.4f  max=%.4f\n",
            t_matched, min(diag(K0_matched)), max(diag(K0_matched))))

# ==============================================================================
# 8. Profile likelihood
# ==============================================================================

cat("\n== Profile likelihood ==\n")
cat(sprintf("  sigma2_e = %.8f (true injected, fixed)\n", sigma2_e))

ll_vario   <- eval_profile_grid(
  K0_vario, y_sub, psill_grid, sigma2_e,
  sprintf("Range=768m  (rho~%.4f)", rho_vario)
)
ll_matched <- eval_profile_grid(
  K0_matched, y_sub, psill_grid, sigma2_e,
  sprintf("Range=%.0fm (rho=%.2f, SAR-matched)", RANGE_MATCHED, RHO_SAR)
)

# ==============================================================================
# 9. Results table
# ==============================================================================

best_vi <- which.max(ll_vario)
best_mi <- which.max(ll_matched)

cat("\n")
cat("================================================================\n")
cat(sprintf("  Profile likelihood results  [n=%d subsample]\n", n_sub))
cat("================================================================\n\n")
cat(sprintf("  %-10s  %-16s  %-16s\n",
            "psill",
            "ll (768m)",
            sprintf("ll (%.0fm)", RANGE_MATCHED)))
cat(strrep("-", 48), "\n")
for (k in seq_along(psill_grid)) {
  m1 <- if (k == best_vi) "*" else " "
  m2 <- if (k == best_mi) "*" else " "
  cat(sprintf("  %-10.5f  %s%-14.2f  %s%-14.2f\n",
              psill_grid[k], m1, ll_vario[k], m2, ll_matched[k]))
}
cat("  (* = MLE for that range)\n\n")

cat(sprintf("  MLE psill  768m    : %.5f\n", psill_grid[best_vi]))
cat(sprintf("  MLE psill  %.0fm  : %.5f\n",  RANGE_MATCHED, psill_grid[best_mi]))
cat(sprintf("  True field variance: %.5f\n", true_var))
cat(sprintf("  Variogram psill    : %.5f\n", M1_PSILL_VARIO))
cat(sprintf("\n  Max ll  768m    : %.2f\n", max(ll_vario)))
cat(sprintf("  Max ll  %.0fm  : %.2f\n",   RANGE_MATCHED, max(ll_matched)))
cat(sprintf("  Diff at MLEs (matched - vario): %.2f\n",
            max(ll_matched) - max(ll_vario)))
cat(sprintf("\n  K formation: %.1fs (768m)  %.1fs (%.0fm)\n",
            t_vario, t_matched, RANGE_MATCHED))

cat("\n  Interpretation:\n")
cat("  [psill recovery]\n")
cat("    MLE psill >> variogram psill -> MLE recovers field variance\n")
cat("    that variogram WLS cannot. Supports likelihood-based fitting.\n")
cat("    MLE psill ~ variogram psill  -> Fundamental identification\n")
cat("    failure. Neither method can recover true field variance from\n")
cat("    the observation likelihood alone.\n")
cat("  [range sensitivity]\n")
cat("    Max ll similar at both ranges -> range doesn't matter in\n")
cat("    this observation regime, consistent with rho insensitivity\n")
cat("    in the SAR model. Neither range nor rho needs careful tuning.\n")
cat("    Max ll differs substantially  -> range matters and requires\n")
cat("    tuning, which means repeated K formation -- the core cost\n")
cat("    argument against kriging in your paper.\n")
cat("================================================================\n")

# ==============================================================================
# 10. Save
# ==============================================================================

results <- list(
  n_sub          = n_sub,
  sub_idx        = sub_idx,
  seed           = SEED,
  psill_grid     = psill_grid,
  ll_vario       = ll_vario,
  ll_matched     = ll_matched,
  best_vario     = psill_grid[best_vi],
  best_matched   = psill_grid[best_mi],
  true_var       = true_var,
  sigma2_e       = sigma2_e,
  range_vario    = RANGE_VARIO,
  range_matched  = RANGE_MATCHED,
  rho_vario      = rho_vario,
  t_K_vario      = t_vario,
  t_K_matched    = t_matched
)

saveRDS(results, "kriging_profile_subset_results.rds")
cat("\nSaved to kriging_profile_subset_results.rds\n")
cat("To scale up: set N_SUB <- Inf at top of script.\n")
