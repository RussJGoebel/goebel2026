# kriging_profile_subset.R
#
# Profiles the marginal log-likelihood over psill for ATA kriging,
# testing two Matern nu=1 range values:
#
#   - 768m  : variogram-fitted range (from existing script)
#   - 2333m : SAR-matched range for rho=0.99, h=330m, d=2
#             derived from kappa^2 = (2d/h^2)(1-rho)/rho
#                          range   = sqrt(2*nu) / kappa = sqrt(2) / kappa
#
# Uses a random subsample of n=500 soundings for a fast test run.
# Each Cholesky on a 500x500 matrix is negligible; K formation is the
# only real cost (~500^2 / 6850^2 * 20min ~ 3 seconds if parallelized).
#
# Outputs:
#   - Profile likelihood table for each range x psill combination
#   - MLE psill for each range
#   - Comparison of the two ranges at their respective MLEs
#   - Timing for K formation and Cholesky solves

library(Matrix)
library(goebel2026)

set.seed(42)

# ------------------------------------------------------------------------------
# 0. SAR <-> Matern range matching
# ------------------------------------------------------------------------------

H     <- 330    # latent pixel size in meters
D     <- 2      # spatial dimension
RHO   <- 0.99   # SAR rho

# From table: kappa^2 = (2d / h^2) * (1 - rho) / rho
kappa2_matched <- (2 * D / H^2) * ((1 - RHO) / RHO)
kappa_matched  <- sqrt(kappa2_matched)

# Matern range = sqrt(2 * nu) / kappa; nu = 1 here
range_matched  <- sqrt(2 * 1) / kappa_matched

# Verify: rho from variogram range 768m
kappa_vario  <- sqrt(2 * 1) / 768
kappa2_vario <- kappa_vario^2
rho_vario    <- (2 * D / H^2) / (kappa2_vario + 2 * D / H^2)

message("== SAR <-> Matern range matching ==")
message(sprintf("  rho=%.2f, h=%dm, d=%d", RHO, H, D))
message(sprintf("  Matched kappa^2 : %.8f", kappa2_matched))
message(sprintf("  Matched range   : %.1fm  (SAR-equivalent to rho=%.2f)", range_matched, RHO))
message(sprintf("  Variogram range : 768m  (equivalent rho ~ %.4f)", rho_vario))
message("")

RANGE_VARIO   <- 768
RANGE_MATCHED <- round(range_matched)

# ------------------------------------------------------------------------------
# 1. Load data and subsample
# ------------------------------------------------------------------------------

message("== Loading data ==")

d_shared      <- goebel2026::setup_shared
d_albedo      <- goebel2026::setup_albedo
y_full        <- d_albedo$y
sigma2_e_true <- d_albedo$sigma_eps^2

n_full <- length(y_full)
n_sub  <- 500

sub_idx <- sample(seq_len(n_full), n_sub, replace = FALSE)
y_sub   <- y_full[sub_idx]

message(sprintf("  Full n: %d  ->  Subsample n: %d", n_full, n_sub))
message(sprintf("  sigma2_e (true injected): %.8f", sigma2_e_true))
message(sprintf("  y_sub mean: %.4f  sd: %.4f", mean(y_sub), sd(y_sub)))
message("")

# ------------------------------------------------------------------------------
# 2. Form K for each range
#
# Assumes a function build_ata_K(footprint_list, range, nu) exists that
# returns a dense n x n ATA Matern covariance matrix with unit psill.
# Replace this block with your actual integration call.
# ------------------------------------------------------------------------------

build_K0 <- function(footprints, range_m, nu = 1) {
  # Wrapper around your ATA integration architecture.
  # Should return the n x n covariance matrix with psill = 1
  # (i.e., the normalized K0 such that K(psill) = psill * K0).
  #
  # Replace with your actual call, e.g.:
  #   ata_matern_K(footprints, range = range_m, nu = nu, psill = 1)
  stop("Replace build_K0() with your actual ATA integration call.")
}

footprints_sub <- d_shared$footprints[sub_idx]   # adjust to your geometry object

message("== Forming K0 for variogram range (768m) ==")
t0 <- proc.time()
K0_vario <- build_K0(footprints_sub, range_m = RANGE_VARIO)
t_vario  <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Done. Elapsed: %.1fs", t_vario))
message(sprintf("  K0 diag range: [%.4f, %.4f]  (should be ~1.0)",
                min(diag(K0_vario)), max(diag(K0_vario))))

message(sprintf("\n== Forming K0 for SAR-matched range (%dm) ==", RANGE_MATCHED))
t0 <- proc.time()
K0_matched <- build_K0(footprints_sub, range_m = RANGE_MATCHED)
t_matched  <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Done. Elapsed: %.1fs", t_matched))
message(sprintf("  K0 diag range: [%.4f, %.4f]  (should be ~1.0)",
                min(diag(K0_matched)), max(diag(K0_matched))))

# ------------------------------------------------------------------------------
# 3. Profile likelihood function
# ------------------------------------------------------------------------------

# Profile likelihood over psill with sigma2_e fixed.
# Sigma = psill * K0 + sigma2_e * I
# Profile out psill analytically? No -- psill and sigma2_e are not
# separable in the same way as sigma2_e alone. We fix sigma2_e at the
# true injected value (as in the existing script) and profile over psill.
#
# If you want to also profile over sigma2_e, set sigma2_e = NULL and
# this function will estimate it by MLE at each psill via:
#   sigma2_e_hat = (y' Sigma^{-1} y) / n   [only valid if psill -> 0]
# That's not right for joint profiling -- for now fix sigma2_e.

eval_profile <- function(K0, y, psill, sigma2_e, label = "") {
  n   <- length(y)
  t0  <- proc.time()

  Sigma <- psill * K0
  diag(Sigma) <- diag(Sigma) + sigma2_e

  # Jitter for numerical stability
  jitter <- 1e-10 * mean(diag(Sigma))
  diag(Sigma) <- diag(Sigma) + jitter

  CK <- tryCatch(
    chol(Sigma),
    error = function(e) {
      message(sprintf("    Cholesky failed at %s, adding larger jitter", label))
      diag(Sigma) <<- diag(Sigma) + 1e-4 * mean(diag(Sigma))
      chol(Sigma)
    }
  )

  log_det <- 2 * sum(log(diag(CK)))
  Kinv_y  <- backsolve(CK, forwardsolve(t(CK), y))
  quad    <- as.numeric(crossprod(y, Kinv_y))
  ll      <- -0.5 * (n * log(2 * pi) + log_det + quad)

  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  %-45s  ll=%10.2f  elapsed=%.2fs", label, ll, elapsed))
  list(ll = ll, elapsed = elapsed)
}

# ------------------------------------------------------------------------------
# 4. Psill grid
# ------------------------------------------------------------------------------

# Use the same grid as the existing script plus the true field variance
y_latent_true <- d_albedo$y_latent_true
fine_grid_sf  <- d_shared$fine_grid_buffered
target_idx    <- which(fine_grid_sf$n_intersects > 0)
true_field_var <- var(y_latent_true[target_idx], na.rm = TRUE)

psill_grid <- sort(unique(c(
  0.00010,
  0.00020,
  0.00034,   # variogram-fitted psill from existing script
  0.00050,
  0.00100,   # oracle / true field variance approx
  true_field_var,
  0.00200,
  0.00500
)))

message(sprintf("\n  True field variance: %.6f", true_field_var))
message(sprintf("  psill grid: %s\n",
                paste(sprintf("%.5f", psill_grid), collapse = ", ")))

# ------------------------------------------------------------------------------
# 5. Evaluate profile likelihood for both ranges
# ------------------------------------------------------------------------------

message("== Profile likelihood: variogram range (768m) ==\n")
ll_vario <- sapply(psill_grid, function(p) {
  res <- eval_profile(K0_vario, y_sub, psill = p,
                      sigma2_e = sigma2_e_true,
                      label = sprintf("range=768m  psill=%.5f", p))
  res$ll
})

message("\n== Profile likelihood: SAR-matched range (2333m) ==\n")
ll_matched <- sapply(psill_grid, function(p) {
  res <- eval_profile(K0_matched, y_sub, psill = p,
                      sigma2_e = sigma2_e_true,
                      label = sprintf("range=%dm  psill=%.5f", RANGE_MATCHED, p))
  res$ll
})

# ------------------------------------------------------------------------------
# 6. Results
# ------------------------------------------------------------------------------

best_vario   <- psill_grid[which.max(ll_vario)]
best_matched <- psill_grid[which.max(ll_matched)]

cat("\n")
cat("=================================================================\n")
cat("  Profile likelihood results -- n=500 subsample\n")
cat("=================================================================\n\n")

cat(sprintf("  %-10s  %-12s  %-12s\n", "psill", "ll (768m)", sprintf("ll (%dm)", RANGE_MATCHED)))
cat(strrep("-", 40), "\n")
for (k in seq_along(psill_grid)) {
  m1 <- if (k == which.max(ll_vario))   " <- MLE" else ""
  m2 <- if (k == which.max(ll_matched)) " <- MLE" else ""
  cat(sprintf("  %-10.5f  %-12.2f%s  %-12.2f%s\n",
              psill_grid[k], ll_vario[k], m1, ll_matched[k], m2))
}

cat(sprintf("\n  MLE psill (768m):    %.5f\n", best_vario))
cat(sprintf("  MLE psill (%dm): %.5f\n", RANGE_MATCHED, best_matched))
cat(sprintf("  True field variance: %.5f\n", true_field_var))
cat(sprintf("  Variogram psill:     %.5f\n", 0.00034))

cat(sprintf("\n  K formation time: %.1fs (768m)  %.1fs (%dm)\n",
            t_vario, t_matched, RANGE_MATCHED))

cat("\n  Interpretation guide:\n")
cat("  - If MLE psill >> variogram psill: MLE recovers field variance\n")
cat("    that variogram WLS cannot (supports your method's advantage).\n")
cat("  - If MLE psill ~ variogram psill: fundamental identification\n")
cat("    failure -- observations don't constrain psill well regardless\n")
cat("    of fitting method.\n")
cat("  - If ll(768m) ~ ll(2333m) at their respective MLEs: range\n")
cat("    doesn't matter much in this observation regime (consistent\n")
cat("    with rho insensitivity in your model).\n")
cat("  - If ll(768m) >> ll(2333m): range matters and variogram range\n")
cat("    is meaningfully better than SAR-matched range.\n")
cat("=================================================================\n")

# ------------------------------------------------------------------------------
# 7. Save results for full run
# ------------------------------------------------------------------------------

results <- list(
  n_sub         = n_sub,
  sub_idx       = sub_idx,
  psill_grid    = psill_grid,
  ll_vario      = ll_vario,
  ll_matched    = ll_matched,
  best_vario    = best_vario,
  best_matched  = best_matched,
  true_field_var = true_field_var,
  sigma2_e_true = sigma2_e_true,
  range_vario   = RANGE_VARIO,
  range_matched = RANGE_MATCHED,
  t_K_vario     = t_vario,
  t_K_matched   = t_matched
)

saveRDS(results, "kriging_profile_subset_results.rds")
message("\nResults saved to kriging_profile_subset_results.rds")
message("To expand to full n=6850, change n_sub <- n_full and re-run.")
