# profile_likelihood_psill.R
#
# Profiles the marginal log-likelihood over psill using the saved
# Matern nu=1 K matrix, without rebuilding K.
#
# Since K(psill) = psill * K0 where K0 is K with unit psill,
# we rescale the saved K_m1 by M1_PSILL to get K0, then
# evaluate the likelihood at a grid of psill values.
#
# Each evaluation requires one Cholesky (~5 min) so use a small grid.
# Key question: does MLE prefer psill ~ 0.00034 (variogram) or ~ 0.001 (oracle)?

library(Matrix)
library(goebel2026)

# ------------------------------------------------------------------------------
# 0. Load saved objects
# ------------------------------------------------------------------------------

message("== Loading saved objects ==")

d_shared      <- goebel2026::setup_shared
d_albedo      <- goebel2026::setup_albedo
y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true
fine_grid_sf  <- d_shared$fine_grid_buffered
target_idx    <- which(fine_grid_sf$n_intersects > 0)

n            <- length(y_alb)
SIGMA2_E     <- d_albedo$sigma_eps^2   # true injected -- used as known
M1_PSILL     <- 0.00034               # fitted Matern nu=1 psill
M1_RANGE     <- 768                   # fitted Matern nu=1 range

message(sprintf("  n=%d  sigma2_e=%.8f", n, SIGMA2_E))
message(sprintf("  Matern nu=1 fitted: psill=%.5f  range=%.0fm",
                M1_PSILL, M1_RANGE))
message(sprintf("  True field variance: %.6f",
                var(y_latent_true[target_idx], na.rm=TRUE)))

message("\nLoading K_m1 from ata_K_matern1.rds...")
K_m1 <- readRDS("ata_K_matern1.rds")
message(sprintf("  K_m1 dim: %d x %d", nrow(K_m1), ncol(K_m1)))
message(sprintf("  K_m1 diag range: [%.6f, %.6f]",
                min(diag(K_m1)), max(diag(K_m1))))

# ------------------------------------------------------------------------------
# 1. Get unit-psill K0
# ------------------------------------------------------------------------------

# K_m1 was built with psill = M1_PSILL
# K(psill) = psill * K0 where K0 = K_m1 / M1_PSILL
K0 <- K_m1 / M1_PSILL

message(sprintf("\n  K0 diag range: [%.6f, %.6f]  (should be ~1.0)",
                min(diag(K0)), max(diag(K0))))

# ------------------------------------------------------------------------------
# 2. Profile likelihood function
# ------------------------------------------------------------------------------

eval_loglik_psill <- function(psill, sigma2_e, label = NULL) {
  if (is.null(label)) label <- sprintf("psill=%.6f", psill)
  t0 <- proc.time()

  K_obs <- psill * K0
  diag(K_obs) <- diag(K_obs) + sigma2_e + 1e-10 * mean(diag(K_obs))

  CK <- tryCatch(
    chol(K_obs),
    error = function(e) {
      diag(K_obs) <<- diag(K_obs) + 1e-4 * mean(diag(K_obs))
      chol(K_obs)
    }
  )

  CKt     <- as(CK, "dtrMatrix")
  log_det <- 2 * sum(log(diag(CK)))
  Kinv_y  <- as.numeric(solve(CKt, solve(t(CKt), y_alb)))
  quad    <- as.numeric(crossprod(y_alb, Kinv_y))
  ll      <- -0.5 * (n * log(2 * pi) + log_det + quad)

  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  message(sprintf("  %-35s  ll=%10.2f  elapsed=%.1fs", label, ll, elapsed))
  ll
}

# ------------------------------------------------------------------------------
# 3. Evaluate at key psill values
# ------------------------------------------------------------------------------

message("\n== Profiling likelihood over psill (Matern nu=1, range=768m) ==")
message(sprintf("  sigma2_e fixed at: %.8f (true injected)", SIGMA2_E))
message("")

# Key values to test:
#   - variogram fitted psill
#   - true field variance
#   - several intermediate values
psill_grid <- c(
  0.00010,                                    # very small
  0.00020,                                    # intermediate
  M1_PSILL,                                   # variogram fitted (0.00034)
  0.00050,                                    # intermediate
  0.00100,                                    # oracle / true field variance
  0.00200,                                    # larger
  0.00500                                     # very large
)

labels <- c(
  "psill=0.00010 (very small)",
  "psill=0.00020 (intermediate)",
  sprintf("psill=%.5f (variogram fitted)", M1_PSILL),
  "psill=0.00050 (intermediate)",
  sprintf("psill=%.5f (true field var)", var(y_latent_true[target_idx], na.rm=TRUE)),
  "psill=0.00200 (larger)",
  "psill=0.00500 (very large)"
)

ll_vals <- mapply(eval_loglik_psill,
                  psill  = psill_grid,
                  label  = labels,
                  MoreArgs = list(sigma2_e = SIGMA2_E))

# ------------------------------------------------------------------------------
# 4. Results
# ------------------------------------------------------------------------------

best_idx  <- which.max(ll_vals)

cat(sprintf("\n=== Profile likelihood results ===\n"))
cat(sprintf("  %-35s  %10s  %10s\n", "Label", "psill", "log-lik"))
cat(strrep("-", 60), "\n")
for (k in seq_along(psill_grid)) {
  marker <- if (k == best_idx) " <-- MLE" else ""
  cat(sprintf("  %-35s  %10.5f  %10.2f%s\n",
              labels[k], psill_grid[k], ll_vals[k], marker))
}

cat(sprintf("\n  MLE psill:        %.5f\n", psill_grid[best_idx]))
cat(sprintf("  Variogram psill:  %.5f\n", M1_PSILL))
cat(sprintf("  True field var:   %.5f\n", var(y_latent_true[target_idx], na.rm=TRUE)))
cat(sprintf("  ll at MLE:        %.2f\n", max(ll_vals)))
cat(sprintf("  ll at variogram:  %.2f\n", ll_vals[which.min(abs(psill_grid - M1_PSILL))]))
cat(sprintf("  ll at oracle:     %.2f\n", ll_vals[which.min(abs(psill_grid - 0.001))]))
cat(sprintf("  ll diff (oracle - variogram): %.2f\n",
            ll_vals[which.min(abs(psill_grid - 0.001))] -
              ll_vals[which.min(abs(psill_grid - M1_PSILL))]))

cat(sprintf("\n  -> If MLE psill >> variogram psill: full likelihood\n"))
cat(sprintf("     identifies field variance that variogram WLS cannot.\n"))
cat(sprintf("  -> If MLE psill ≈ variogram psill: fundamental identification\n"))
cat(sprintf("     failure -- neither variogram nor MLE can recover field variance.\n"))
