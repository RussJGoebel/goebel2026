set.seed(42)

# use soundings_augmented for synthetic albedo response
truth_vals <- goebel2026::target_grid$mean_albedo
noise_sd   <- sd(truth_vals, na.rm = TRUE) / 20

y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

# precompute shared objects
W   <- goebel2026::make_W_matrix(goebel2026::target_grid)
A   <- spatintegrate::compute_overlap_fractions(
  spatintegrate::ensure_projected(goebel2026::soundings),
  spatintegrate::ensure_projected(goebel2026::target_grid)
)
A <- as(A, "dgCMatrix")

IminusW <- Matrix::Diagonal(nrow(W)) - 1.0 * W
Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
Q       <- Matrix::drop0(Q)
Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  Q       <- Matrix::drop0(Q)
  list(Q = Q)
}

### CV ###
future::plan(future::multisession())
tuned_cv <- fastblm::tune_cv(
  y          = y_alb,
  A          = A,
  Q_fun      = Q_fun,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 0.999,
  k          = 10L,
  solver     = "cholesky",
  parallel   = TRUE,
  verbose    = TRUE
)
future::plan(future::sequential())

### REML ###
tuned_reml <- fastblm::tune_reml(
  y             = y_alb,
  A             = A,
  Q_fun         = Q_fun,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 0.999,
  solver        = "cholesky",
  logdet_method = "cholesky",
  verbose       = TRUE
)

### fit both at their respective phi ###
fit_cv <- fastblm::fit_fastblm(y_alb, A, Q, tuned_cv$phi,   solver = "cholesky")
fit_reml <- fastblm::fit_fastblm(y_alb, A, Q, tuned_reml$phi, solver = "cholesky")

se_cv   <- fastblm::posterior_se(fit_cv,   A_new = diag(dim(A)[2]))
se_reml <- fastblm::posterior_se(fit_reml, A_new = diag(dim(A)[2]))

### compare against known truth ###
truth <- goebel2026::target_grid$mean_albedo
na_idx <- is.na(truth)

# RMSE of posterior mean
rmse_cv   <- sqrt(mean((fit_cv$posterior_mean[!na_idx]   - truth[!na_idx])^2))
rmse_reml <- sqrt(mean((fit_reml$posterior_mean[!na_idx] - truth[!na_idx])^2))

# coverage of 95% credible intervals
in_ci <- function(fit, se, truth) {
  lower <- fit$posterior_mean - 1.96 * se
  upper <- fit$posterior_mean + 1.96 * se
  mean(truth >= lower & truth <= upper, na.rm = TRUE)
}
cov_cv   <- in_ci(fit_cv,   se_cv,   truth)
cov_reml <- in_ci(fit_reml, se_reml, truth)

cat("\n--- CV ---\n")
cat(sprintf("phi:      %.4f\n", tuned_cv$phi))
cat(sprintf("sigma2e:  %.4f\n", tuned_cv$sigma2e))
cat(sprintf("RMSE:     %.4f\n", rmse_cv))
cat(sprintf("Coverage: %.4f\n", cov_cv))

cat("\n--- REML ---\n")
cat(sprintf("phi:      %.4f\n", tuned_reml$phi))
cat(sprintf("sigma2e:  %.4f\n", tuned_reml$sigma2e))
cat(sprintf("RMSE:     %.4f\n", rmse_reml))
cat(sprintf("Coverage: %.4f\n", cov_reml))
