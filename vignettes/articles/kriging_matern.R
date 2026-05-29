# kriging_profile_subset.R
#
# Profile likelihood over psill for ATA kriging, testing two range values:
#
#   - RANGE_VARIO   : variogram-fitted range from existing pipeline
#   - RANGE_MATCHED : SAR-matched range for rho=0.99, h=330m, d=2
#                     Table 5.1: kappa^2 = (2d/h^2)(1-rho)/rho
#                                range   = sqrt(2*nu)/kappa  [nu=1]
#
# Covariance: Matern nu=1 (Whittle kernel)
#   C(h) = psill * (1 + h/range) * exp(-h/range)
#
# sigma2_e fixed at true injected value.
# n=500 subsample for fast test run.

library(spatintegrate)
library(goebel2026)
library(sf)
library(Matrix)
library(future)
library(future.apply)

# ==============================================================================
# 0. Configuration
# ==============================================================================

N_SUB          <- 500L
N_PER_TRIANGLE <- 16L
N_WORKERS      <- 4L
K_CHUNK        <- 20L
SEED           <- 42L

H       <- 330
D       <- 2
RHO_SAR <- 0.99

PSILL_VARIO <- 0.00034
RANGE_VARIO <- 768

# ==============================================================================
# 1. SAR <-> Matern nu=1 range matching
# ==============================================================================

kappa2_matched <- (2 * D / H^2) * ((1 - RHO_SAR) / RHO_SAR)
kappa_matched  <- sqrt(kappa2_matched)
RANGE_MATCHED  <- sqrt(2 * 1) / kappa_matched   # sqrt(2*nu)/kappa, nu=1

kappa_vario <- sqrt(2 * 1) / RANGE_VARIO
rho_vario   <- (2 * D / H^2) / (kappa_vario^2 + 2 * D / H^2)

cat("== SAR <-> Matern nu=1 range matching ==\n")
cat(sprintf("  rho=%.2f, h=%dm, d=%d\n", RHO_SAR, H, D))
cat(sprintf("  SAR-matched range : %.1fm\n", RANGE_MATCHED))
cat(sprintf("  Variogram range   : %.1fm  (equivalent rho ~ %.4f)\n\n",
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

set.seed(SEED)
sub_idx           <- sort(sample(seq_len(n_full), n_sub, replace = FALSE))
y_sub             <- y_full[sub_idx]
sounding_geom_sub <- sf::st_geometry(soundings_sf)[sub_idx]

cat(sprintf("  n_full=%d  n_sub=%d  sigma2_e=%.8f\n",
            n_full, n_sub, sigma2_e))
cat(sprintf("  True field variance : %.6f\n", true_var))
cat(sprintf("  Variogram psill     : %.6f\n\n", PSILL_VARIO))

# ==============================================================================
# 3. Build QMC quadrature points for subsample
# ==============================================================================

cat("== Building QMC quadrature points ==\n")
t0 <- proc.time()

qmc_base <- spatintegrate::generate_qmc_unit_square(N_PER_TRIANGLE)

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

sounding_pts <- vector("list", n_sub)
for (i in seq_len(n_sub)) {
  if (i %% 100L == 0L) cat(sprintf("  QMC: %d / %d\n", i, n_sub))
  sounding_pts[[i]] <- get_poly_pts(sounding_geom_sub[[i]], crs_proj)
}

n_pts_vec <- vapply(sounding_pts, nrow, integer(1L))
n_pts     <- as.integer(median(n_pts_vec))
cat(sprintf("  Points per sounding: median=%d  range=[%d,%d]\n",
            n_pts, min(n_pts_vec), max(n_pts_vec)))

# Stack into (n_sub * n_pts) x 2, padding/trimming to median n_pts
S_all <- do.call(rbind, lapply(sounding_pts, function(pts) {
  if (nrow(pts) >= n_pts) pts[seq_len(n_pts), , drop = FALSE]
  else pts[sample(nrow(pts), n_pts, replace = TRUE), , drop = FALSE]
}))
bb <- rowSums(S_all^2)

t_qmc <- as.numeric((proc.time() - t0)["elapsed"])
cat(sprintf("  Done. %.1fs  S_all: %d x %d\n\n", t_qmc, nrow(S_all), ncol(S_all)))

# ==============================================================================
# 4. ATA K formation
# ==============================================================================
#
# Matern nu=1: C(h) = (1 + h/r) * exp(-h/r)   [unit psill]
# K0[i,j] = (1/n_pts^2) * sum_{s in D_i, t in D_j} C(||s-t||)
#
# Key fix: all variables used inside future workers must be listed
# explicitly in future.globals -- closures over the enclosing environment
# do not work reliably across multisession workers.

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

      H2 <- outer(aa, bb, "+") - 2 * tcrossprod(S_chunk, S_all)
      H  <- sqrt(pmax(H2, 0))

      # Matern nu=1
      u  <- H / range_val
      CV <- (1 + u) * exp(-u)

      # Average over j quadrature points
      col_grp  <- rep(seq_len(n), each = np)
      CV_j_avg <- t(rowsum(t(CV), col_grp) / np)

      # Average over i quadrature points
      row_grp <- rep(seq_len(n_chunk), each = np)
      K_chunk <- rowsum(CV_j_avg, row_grp) / np

      list(i_idx = i_idx, K_chunk = K_chunk)
    },
    future.seed    = TRUE,
    future.globals = list(
      S_all     = S_all,
      bb        = bb,
      n         = n,
      n_pts     = n_pts,
      range_val = range_val   # explicitly passed, not closed over
    )
  )

  K <- matrix(0, nrow = n, ncol = n)
  for (res in K_chunks) K[res$i_idx, ] <- res$K_chunk
  0.5 * (K + t(K))
}

cat("== Forming K0: variogram range (768m) ==\n")
t0       <- proc.time()
K0_vario <- build_K0(S_all, bb, n_sub, n_pts,
                     range_val = RANGE_VARIO,
                     n_workers = N_WORKERS,
                     k_chunk   = K_CHUNK)
t_vario  <- as.numeric((proc.time() - t0)["elapsed"])
cat(sprintf("  Done. %.1fs  |  diag: min=%.4f  mean=%.4f  max=%.4f\n",
            t_vario,
            min(diag(K0_vario)), mean(diag(K0_vario)), max(diag(K0_vario))))

cat(sprintf("\n== Forming K0: SAR-matched range (%.0fm) ==\n", RANGE_MATCHED))
t0         <- proc.time()
K0_matched <- build_K0(S_all, bb, n_sub, n_pts,
                       range_val = RANGE_MATCHED,
                       n_workers = N_WORKERS,
                       k_chunk   = K_CHUNK)
t_matched  <- as.numeric((proc.time() - t0)["elapsed"])
cat(sprintf("  Done. %.1fs  |  diag: min=%.4f  mean=%.4f  max=%.4f\n",
            t_matched,
            min(diag(K0_matched)), mean(diag(K0_matched)), max(diag(K0_matched))))

# ==============================================================================
# 5. Profile likelihood (fixed sigma2_e, vary psill)
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

# Initial psill grid -- will be extended if MLE hits the upper boundary
psill_grid_init <- sort(unique(c(
  0.00010, 0.00020, PSILL_VARIO, 0.00050,
  0.00100, true_var, 0.00200, 0.00500
)))

cat(sprintf("\nInitial psill grid: %s\n",
            paste(sprintf("%.5f", psill_grid_init), collapse = ", ")))
cat(sprintf("sigma2_e fixed at %.8f (true injected)\n", sigma2_e))

# ------------------------------------------------------------------------------
# Profile likelihood with dynamic upper boundary extension.
# If the MLE lands on the top grid point, double the upper end and evaluate
# only the new points (reuse existing evaluations). Repeat until interior.
# MAX_EXTENSIONS caps the search to avoid runaway loops.
# ------------------------------------------------------------------------------

profile_with_extension <- function(K0, y, psill_grid, sigma2_e,
                                   range_label, max_extensions = 6L) {

  # Evaluate initial grid
  lls <- eval_profile_grid(K0, y, psill_grid, sigma2_e, range_label)

  n_ext <- 0L
  while (which.max(lls) == length(lls) && n_ext < max_extensions) {
    n_ext    <- n_ext + 1L
    psill_hi <- max(psill_grid)

    # Add 3 new points between current max and 2x current max
    new_psills <- sort(unique(c(
      psill_hi * 1.5,
      psill_hi * 2.0,
      psill_hi * 3.0
    )))

    cat(sprintf(
      "\n  [boundary hit, extension %d] MLE at upper edge (%.5f). Extending to %.5f\n",
      n_ext, psill_hi, max(new_psills)
    ))

    new_lls <- eval_profile_grid(K0, y, new_psills, sigma2_e,
                                 sprintf("%s [ext %d]", range_label, n_ext))

    psill_grid <- c(psill_grid, new_psills)
    lls        <- c(lls, new_lls)

    # Re-sort
    ord        <- order(psill_grid)
    psill_grid <- psill_grid[ord]
    lls        <- lls[ord]
  }

  if (which.max(lls) == length(lls)) {
    cat(sprintf(
      "\n  WARNING: MLE still at upper boundary (%.5f) after %d extensions.\n",
      max(psill_grid), n_ext
    ))
    cat("  Likelihood may be monotonically increasing -- psill not identified.\n")
  }

  list(psill_grid = psill_grid, lls = lls, n_extensions = n_ext)
}

cat("\n== Profile likelihood ==\n")

res_vario <- profile_with_extension(
  K0_vario, y_sub, psill_grid_init, sigma2_e,
  sprintf("Range=%.0fm  (rho_equiv~%.4f)", RANGE_VARIO, rho_vario)
)

res_matched <- profile_with_extension(
  K0_matched, y_sub, psill_grid_init, sigma2_e,
  sprintf("Range=%.0fm (SAR-matched rho=%.2f)", RANGE_MATCHED, RHO_SAR)
)

# Align to a common grid for the results table
ll_vario      <- res_vario$lls
ll_matched    <- res_matched$lls
psill_grid_v  <- res_vario$psill_grid
psill_grid_m  <- res_matched$psill_grid

# ==============================================================================
# 6. Results
# ==============================================================================

best_vi <- which.max(ll_vario)
best_mi <- which.max(ll_matched)

cat("\n")
cat("================================================================\n")
cat(sprintf("  Profile likelihood  [n=%d subsample, seed=%d]\n", n_sub, SEED))
cat(sprintf("  Covariance: Matern nu=1  C(h) = (1+h/r)*exp(-h/r)\n"))
cat(sprintf("  sigma2_e = %.8f (true injected, fixed)\n", sigma2_e))
cat("================================================================\n\n")

# Print each range separately since grids may differ after extension
cat(sprintf("  -- Range=%.0fm (rho~%.4f) --\n", RANGE_VARIO, rho_vario))
cat(sprintf("  %-12s  %s\n", "psill", "log-lik"))
cat(strrep("-", 30), "\n")
for (k in seq_along(psill_grid_v)) {
  m <- if (k == best_vi) " *" else "  "
  cat(sprintf("  %-12.5f  %.2f%s\n", psill_grid_v[k], ll_vario[k], m))
}

cat(sprintf("\n  -- Range=%.0fm (SAR-matched rho=%.2f) --\n",
            RANGE_MATCHED, RHO_SAR))
cat(sprintf("  %-12s  %s\n", "psill", "log-lik"))
cat(strrep("-", 30), "\n")
for (k in seq_along(psill_grid_m)) {
  m <- if (k == best_mi) " *" else "  "
  cat(sprintf("  %-12.5f  %.2f%s\n", psill_grid_m[k], ll_matched[k], m))
}
cat("  (* = MLE)\n\n")

cat(sprintf("  MLE psill  %.0fm   : %.5f\n", RANGE_VARIO,   psill_grid_v[best_vi]))
cat(sprintf("  MLE psill  %.0fm : %.5f\n",   RANGE_MATCHED, psill_grid_m[best_mi]))
cat(sprintf("  True field variance  : %.5f\n", true_var))
cat(sprintf("  Variogram psill      : %.5f\n", PSILL_VARIO))
cat(sprintf("\n  Max ll  %.0fm   : %.2f\n", RANGE_VARIO,   max(ll_vario)))
cat(sprintf("  Max ll  %.0fm : %.2f\n",   RANGE_MATCHED, max(ll_matched)))
cat(sprintf("  Delta ll at MLEs (matched - vario): %.2f\n",
            max(ll_matched) - max(ll_vario)))

cat("\n  Interpretation:\n")
cat("  [psill recovery]\n")
cat("    MLE psill >> variogram psill:\n")
cat("      Likelihood identifies field variance that variogram WLS\n")
cat("      cannot. MLE kriging gives better-calibrated intervals.\n")
cat("    MLE psill ~ variogram psill:\n")
cat("      Fundamental identification failure. Your CV tuning\n")
cat("      accesses prediction error information that the marginal\n")
cat("      likelihood alone cannot.\n")
cat("  [range sensitivity]\n")
cat("    |Delta ll| small (< ~2 per free parameter):\n")
cat("      Range doesn't matter in this observation regime --\n")
cat("      consistent with rho insensitivity in your SAR model.\n")
cat("    |Delta ll| large:\n")
cat("      Range matters, requiring repeated K formation at ~20min\n")
cat("      each. Core computational cost argument against kriging.\n")
cat("================================================================\n")

# ==============================================================================
# 7. Save
# ==============================================================================

results <- list(
  n_sub         = n_sub,
  sub_idx       = sub_idx,
  seed          = SEED,
  psill_grid_v  = psill_grid_v,
  psill_grid_m  = psill_grid_m,
  ll_vario      = ll_vario,
  ll_matched    = ll_matched,
  best_vario    = psill_grid_v[best_vi],
  best_matched  = psill_grid_m[best_mi],
  true_var      = true_var,
  sigma2_e      = sigma2_e,
  range_vario   = RANGE_VARIO,
  range_matched = RANGE_MATCHED,
  rho_vario     = rho_vario,
  t_qmc         = t_qmc,
  t_K_vario     = t_vario,
  t_K_matched   = t_matched,
  n_pts         = n_pts
)

saveRDS(results, "kriging_profile_subset_results.rds")
cat("\nSaved to kriging_profile_subset_results.rds\n")
cat("To scale to full data: set N_SUB <- n_full at top of script.\n")
