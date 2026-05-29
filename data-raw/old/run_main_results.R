# data-raw/run_main_results.R
#
# Produces the main albedo ablation table (Table 2) and the canonical SIF run.
# Requires setup objects saved by data-raw/setup.R via usethis::use_data().
#
# Albedo ablations (all on semi-synthetic albedo, rho = 1 unless noted):
#   1. OLS baseline
#   2. No covariates, rho = 1
#   3. Water covariate, rho = 1  (no RSR)
#   4. Water covariate + RSR, rho = 1
#   5. Water covariate + RSR + R_inv, rho = 1
#   6. Water covariate + RSR, rho CV-tuned
#
# SIF:
#   7. Water covariate + RSR + R_inv, rho = 1  (canonical application run)
#
# Model structure:
#   Spatial-only:    fit on A (n x p),       posterior_mean = r (length p)
#   With covariates: fit on A_aug = [A|X_obs] (n x (p+q)),
#                    Q_aug = bdiag(Q, lambda*I_q),
#                    posterior_mean = [r; beta] (length p+q)
#   RSR:             constrain(fit_aug, C_aug) where C_aug = [t(X_obs)%*%A | 0]
#
# Outputs saved to data/ via usethis::use_data().

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

# Set to TRUE to force rerun all steps even if results already exist in data/
FORCE_RERUN <- FALSE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

soundings_proj     <- d_shared$soundings_proj
fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water    # n x 2: [1, water_obs]
X_latent_water     <- d_shared$X_latent_water # m x 2: [1, water_grid]

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true
sigma_eps     <- d_albedo$sigma_eps

y_sif     <- d_sif$y
R_inv_sif <- d_sif$R_inv

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)   # 2: intercept + water
lambda_beta           <- 0.01               # weak ridge prior on beta
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Shared model components
# ------------------------------------------------------------------------------

# Augmented design matrix for covariate runs: n x (p+q)
A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

# Q_fun factories returning function(theta) -> list(Q, log_det_Q)

make_Q_fun_fixed <- function(W) {
  S <- Matrix::Diagonal(nrow(W)) - W
  Q <- Matrix::forceSymmetric(Matrix::crossprod(S))
  Q <- Matrix::drop0(Q)
  function(theta) list(Q = Q, log_det_Q = NULL)
}

make_Q_fun_rho <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q   <- Matrix::drop0(Q)
    list(Q = Q, log_det_Q = NULL)
  }
}

# Augmented Q_fun: spatial Q_fun + lambda*I_q block for beta
make_Q_fun_aug <- function(Q_fun_spatial, q, lambda) {
  force(Q_fun_spatial); force(q); force(lambda)
  function(theta) {
    sp  <- Q_fun_spatial(theta)
    Q_a <- Matrix::forceSymmetric(
      Matrix::bdiag(sp$Q, lambda * Matrix::Diagonal(q))
    )
    list(Q = Q_a, log_det_Q = NULL)
  }
}

Q_fun_fixed     <- make_Q_fun_fixed(W_queen)
Q_fun_rho       <- make_Q_fun_rho(W_queen)
Q_fun_fixed_aug <- make_Q_fun_aug(Q_fun_fixed, q, lambda_beta)
Q_fun_rho_aug   <- make_Q_fun_aug(Q_fun_rho,   q, lambda_beta)

# RSR constraint function for augmented system.
# Per-fold: C_aug = [t(X_obs[train,]) %*% A_train | 0_{q x q}]
# Forces spatial block r orthogonal to covariates; beta block unconstrained.
make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    C_zeros   <- matrix(0, nrow = q, ncol = q)
    cbind(C_spatial, C_zeros)                   # q x (p+q), base matrix
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

# Full-data RSR constraint matrix (used after tuning)
C_spatial_full <- as.matrix(t(X_obs_water) %*% A)         # q x p
C_zeros_full   <- matrix(0, nrow = q, ncol = q)
C_aug_full     <- cbind(C_spatial_full, C_zeros_full)     # q x (p+q), base matrix

# ------------------------------------------------------------------------------
# 3. Fitting helpers
# ------------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Helper: TRUE if result object already loaded or saved in package data
.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name))
    return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error   = function(e) NULL,
    warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}


tune_and_fit <- function(y, A_fit, Q_fun_fit,
                         rsr      = FALSE,
                         R_inv    = NULL,
                         tune_rho = FALSE,
                         k        = 10L,
                         seed     = 2026L) {

  theta_init <- if (tune_rho) c(rho = 0.9) else numeric(0)
  lower      <- if (tune_rho) 0.01         else numeric(0)
  upper      <- if (tune_rho) 0.99         else numeric(0)
  constraint <- if (rsr) rsr_constraint else NULL

  tuned <- fastblm::tune_cv(
    y          = y,
    A          = A_fit,
    Q_fun      = Q_fun_fit,
    R_inv      = R_inv,
    theta_init = theta_init,
    lower      = lower,
    upper      = upper,
    k          = k,
    constraint = constraint,
    seed       = seed,
    parallel   = TRUE,
    verbose    = TRUE
  )

  fit <- fastblm::fit_fastblm(
    y      = y,
    A      = A_fit,
    Q      = tuned$Q,
    phi    = tuned$phi,
    R_inv  = R_inv,
    solver = "cholesky"
  )

  if (rsr) fit <- fastblm::constrain(fit, C_aug_full)

  list(fit = fit, tuned = tuned)
}

compute_outputs <- function(result, run_name, tags,
                            with_covariates = FALSE,
                            y_true          = y_latent_true) {
  fit   <- result$fit
  tuned <- result$tuned

  if (with_covariates) {
    r_hat    <- fit$posterior_mean[seq_len(p)]
    beta_hat <- fit$posterior_mean[p + seq_len(q)]
    mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
    A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
    se       <- fastblm::posterior_se(fit, A_new = A_pred)
    se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]
  } else {
    mu       <- fit$posterior_mean
    se       <- fastblm::posterior_se(fit)
    beta_hat <- NULL
    se_beta  <- NULL
  }

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  tr_t  <- y_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t
  ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }

  resid  <- mu_t - tr_t
  rmse   <- sqrt(mean(resid^2, na.rm = TRUE))
  ss_res <- sum(resid^2, na.rm = TRUE)
  ss_tot <- sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)
  r2     <- 1 - ss_res / ss_tot

  list(
    run_name              = run_name,
    tags                  = tags,
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = if ("rho" %in% names(tuned$theta)) tuned$theta[["rho"]] else 1,
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

# ------------------------------------------------------------------------------
# 4. Run ablations
# ------------------------------------------------------------------------------

# --- 1. OLS baseline ----------------------------------------------------------

if (!.should_skip("results_ols_baseline")) {
  message("\n== 1. OLS baseline ==")

  lm_fit   <- lm(y_alb ~ X_obs_water[, 2])
  beta_ols <- coef(lm_fit)
  mu_ols   <- beta_ols[1] + beta_ols[2] * X_latent_water[, 2]

  resid_ols <- mu_ols[target_idx] - y_latent_true[target_idx]
  rmse_ols  <- sqrt(mean(resid_ols^2, na.rm = TRUE))
  ss_res    <- sum(resid_ols^2, na.rm = TRUE)
  ss_tot    <- sum((y_latent_true[target_idx] -
                      mean(y_latent_true[target_idx], na.rm = TRUE))^2, na.rm = TRUE)
  r2_ols    <- 1 - ss_res / ss_tot

  results_ols_baseline <- list(
    run_name              = "ols_baseline",
    tags                  = list(tuning = "ols", covariates = "water",
                                 constraint = "none", W = "none", rho = NA),
    timestamp             = Sys.time(),
    posterior_mean        = mu_ols[target_idx],
    posterior_se          = rep(NA_real_, length(target_idx)),
    ci_lower              = rep(NA_real_, length(target_idx)),
    ci_upper              = rep(NA_real_, length(target_idx)),
    beta_hat              = beta_ols,
    se_beta               = summary(lm_fit)$coefficients[, 2],
    sigma2e               = summary(lm_fit)$sigma^2,
    phi                   = NA_real_,
    rho_opt               = NA_real_,
    cv_curve              = NULL,
    rmse                  = rmse_ols,
    r2                    = r2_ols,
    coverage_95_all       = NA_real_,
    coverage_95_obs       = NA_real_,
    coverage_95_dense     = NA_real_,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
  usethis::use_data(results_ols_baseline, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f", rmse_ols, r2_ols))
}

# --- 2. No covariates, rho = 1 ------------------------------------------------

if (!.should_skip("results_no_cov_rho1")) {
  message("\n== 2. No covariates, rho=1 ==")

  r2   <- tune_and_fit(y_alb, A_fit = A, Q_fun_fit = Q_fun_fixed)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r2$tuned$phi, r2$tuned$sigma2e))
  out2 <- compute_outputs(r2, "no_cov_rho1",
                          tags = list(tuning = "cv", covariates = "none",
                                      constraint = "none", W = "queen", rho = 1))
  results_no_cov_rho1 <- out2
  usethis::use_data(results_no_cov_rho1, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out2$rmse, out2$r2, out2$coverage_95_obs))
}

# --- 3. Water covariate, rho = 1, no RSR -------------------------------------

if (!.should_skip("results_water_rho1")) {
  message("\n== 3. Water covariate, rho=1, no RSR ==")

  r3   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r3$tuned$phi, r3$tuned$sigma2e))
  out3 <- compute_outputs(r3, "water_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "none", W = "queen", rho = 1),
                          with_covariates = TRUE)
  results_water_rho1 <- out3
  usethis::use_data(results_water_rho1, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out3$rmse, out3$r2, out3$coverage_95_obs))
}

# --- 4. Water + RSR, rho = 1 --------------------------------------------------

if (!.should_skip("results_water_rsr_rho1")) {
  message("\n== 4. Water + RSR, rho=1 ==")

  r4   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug, rsr = TRUE)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r4$tuned$phi, r4$tuned$sigma2e))
  out4 <- compute_outputs(r4, "water_rsr_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = 1),
                          with_covariates = TRUE)
  results_water_rsr_rho1 <- out4
  usethis::use_data(results_water_rsr_rho1, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out4$rmse, out4$r2, out4$coverage_95_obs))
}

# --- 5. Water + RSR + R_inv, rho = 1 -----------------------------------------
# Albedo noise is homoskedastic so R_inv = NULL (identity weights) here.
# Included to verify the R_inv code path before the SIF run.

if (!.should_skip("results_water_rsr_rinv_rho1")) {
  message("\n== 5. Water + RSR + R_inv, rho=1 (albedo) ==")

  r5   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug,
                       rsr = TRUE, R_inv = NULL)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r5$tuned$phi, r5$tuned$sigma2e))
  out5 <- compute_outputs(r5, "water_rsr_rinv_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = 1,
                                      R_inv = "identity_albedo"),
                          with_covariates = TRUE)
  results_water_rsr_rinv_rho1 <- out5
  usethis::use_data(results_water_rsr_rinv_rho1, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out5$rmse, out5$r2, out5$coverage_95_obs))
}

# --- 6. Water + RSR, rho CV-tuned --------------------------------------------

if (!.should_skip("results_water_rsr_rho_cv")) {
  message("\n== 6. Water + RSR, rho CV-tuned ==")

  r6   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                       rsr = TRUE, tune_rho = TRUE)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r6$tuned$phi, r6$tuned$sigma2e))
  out6 <- compute_outputs(r6, "water_rsr_rho_cv",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = "cv_tuned"),
                          with_covariates = TRUE)
  results_water_rsr_rho_cv <- out6
  usethis::use_data(results_water_rsr_rho_cv, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  rho_opt=%.4f  coverage_obs=%.3f",
                  out6$rmse, out6$r2, out6$rho_opt, out6$coverage_95_obs))
}

# --- 7. SIF canonical run: Water + RSR + R_inv, rho = 1 ----------------------

if (!.should_skip("results_sif_canonical")) {
  message("\n== 7. SIF canonical run ==")

  r7      <- tune_and_fit(y_sif, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug,
                          rsr = TRUE, R_inv = R_inv_sif)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", r7$tuned$phi, r7$tuned$sigma2e))
  fit7    <- r7$fit
  r_hat7  <- fit7$posterior_mean[seq_len(p)]
  beta7   <- fit7$posterior_mean[p + seq_len(q)]
  mu7     <- r_hat7 + as.numeric(X_latent_water %*% beta7)
  A_pred7 <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se7     <- fastblm::posterior_se(fit7, A_new = A_pred7)
  se_beta7 <- fastblm::posterior_se(fit7)[p + seq_len(q)]

  results_sif_canonical <- list(
    run_name              = "sif_canonical",
    tags                  = list(tuning = "cv", covariates = "water",
                                 constraint = "RSR", W = "queen", rho = 1,
                                 R_inv = "sif_uncertainty"),
    timestamp             = Sys.time(),
    posterior_mean        = mu7[target_idx],
    posterior_se          = se7[target_idx],
    ci_lower              = (mu7 - 1.96 * se7)[target_idx],
    ci_upper              = (mu7 + 1.96 * se7)[target_idx],
    beta_hat              = beta7,
    se_beta               = se_beta7,
    sigma2e               = fit7$sigma2e,
    phi                   = r7$tuned$phi,
    rho_opt               = 1,
    cv_curve              = r7$tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
  usethis::use_data(results_sif_canonical, overwrite = TRUE)
  message(sprintf("  phi=%.4f  sigma2e=%.4g  beta=[%.4f, %.4f]",
                  r7$tuned$phi, fit7$sigma2e, beta7[1], beta7[2]))
}

# --- DIAGNOSTIC: Water, no RSR, rho CV-tuned --------------------------------
# Sanity check: rho should tune near 1 without RSR constraint.
# Compare with run 6 (RSR + rho CV-tuned) to understand RSR's effect on rho.

if (!.should_skip("results_diagnostic_rho_norsr")) {
  message("\n== DIAGNOSTIC: Water, no RSR, rho CV-tuned ==")
  r_diag <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                         rsr = FALSE, tune_rho = TRUE)
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  r_diag$tuned$theta[["rho"]], r_diag$tuned$phi,
                  r_diag$tuned$sigma2e))
  out_diag <- compute_outputs(r_diag, "diagnostic_rho_norsr",
                              tags = list(tuning = "cv", covariates = "water",
                                          constraint = "none", W = "queen", rho = "cv_tuned"),
                              with_covariates = TRUE)
  results_diagnostic_rho_norsr <- out_diag
  usethis::use_data(results_diagnostic_rho_norsr, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out_diag$rho_opt, out_diag$rmse, out_diag$r2,
                  out_diag$coverage_95_obs))
}


# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nAll main results saved to data/.")
