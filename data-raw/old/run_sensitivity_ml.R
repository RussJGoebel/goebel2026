# data-raw/run_sensitivity_ml.R
#
# Supplement: Marginal likelihood tuning sensitivity.
#
# Compares ML-tuned phi (and rho) against the CV baseline (results_water_rho1).
# Both use: water covariate, no RSR, augmented system [A | X_obs_water].
# ML handles covariates via a flat (zero-precision) prior on the beta block;
# no RSR is applied since ML selects phi under the unconstrained likelihood.
#
# Outputs saved to data/ via usethis::use_data():
#   results_water_ml       -- ML-tuned phi and rho, water covariate, no RSR
#   results_water_ml_rho1  -- ML-tuned phi, rho fixed at 0.99, water covariate, no RSR
#   results_nocov_ml       -- ML-tuned phi and rho, no covariates
#   results_water_ml_zerolog -- ML-tuned phi and rho, water covariate, logdet(Q)=0

library(fastblm)
library(goebel2026)
library(Matrix)
library(usethis)

FORCE_RERUN <- FALSE

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

# Q_fun: rho-tunable SAR prior, spatial block only.
# tune_ml handles augmentation with X_fixed internally.
make_Q_fun_rho <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q   <- Matrix::drop0(Q)
    list(Q = Q, log_det_Q = NULL)
  }
}

Q_fun_rho <- make_Q_fun_rho(W_queen)

# Q_fun: rho-tunable SAR prior with logdet(Q) forced to zero.
# Used to test whether the logdet(Q) singularity near rho=1 is responsible
# for ML selecting a smaller phi than CV. Setting log_det_Q = 0 drops the
# logdet(Q) term from the likelihood entirely, matching the improper-prior
# convention used for intrinsic SAR models.
make_Q_fun_rho_zerolog <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q   <- Matrix::drop0(Q)
    list(Q = Q, log_det_Q = 0)
  }
}

Q_fun_rho_zerolog <- make_Q_fun_rho_zerolog(W_queen)

# Q_fun: near-intrinsic prior (rho = 0.99).
# rho = 1 is singular so logdet(Q) is undefined; 0.99 is positive definite,
# logdet is finite and constant w.r.t. phi so it does not distort phi selection.
make_Q_fun_fixed <- function(W, rho = 0.99) {
  S <- Matrix::Diagonal(nrow(W)) - rho * W
  Q <- Matrix::forceSymmetric(Matrix::crossprod(S))
  Q <- Matrix::drop0(Q)
  function(theta) list(Q = Q, log_det_Q = NULL)
}

Q_fun_fixed <- make_Q_fun_fixed(W_queen, rho = 0.99)
lambda_beta <- 0.01

# ------------------------------------------------------------------------------
# Shared: final fit and output computation given a tuned result
# ------------------------------------------------------------------------------

fit_and_compute <- function(tuned, rho_val) {
  A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")
  Q_aug <- Matrix::forceSymmetric(
    Matrix::bdiag(tuned$Q, lambda_beta * Matrix::Diagonal(q))
  )

  fit <- fastblm::fit_fastblm(
    y      = y_alb,
    A      = A_aug,
    Q      = Q_aug,
    phi    = tuned$phi,
    solver = "cholesky"
  )

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)

  A_pred  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se      <- fastblm::posterior_se(fit, A_new = A_pred)
  se_beta <- fastblm::posterior_se(fit)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f", beta_hat[1], beta_hat[2]))

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  tr_t  <- y_latent_true[target_idx]
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
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

# ------------------------------------------------------------------------------
# 3. Skip helper
# ------------------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# 4. Run
# ------------------------------------------------------------------------------

# --- ML with rho tuned --------------------------------------------------------

if (!.should_skip("results_water_ml")) {
  message("\n== ML tuning: water covariate, rho ML-tuned, no RSR ==")

  tuned_ml <- fastblm::tune_ml(
    y          = y_alb,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = X_obs_water,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )

  rho_hat <- tuned_ml$theta[["rho"]]
  message(sprintf("  ML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned_ml$phi, tuned_ml$sigma2e))

  out <- fit_and_compute(tuned_ml, rho_hat)

  results_water_ml <- c(
    list(run_name  = "water_ml",
         tags      = list(tuning     = "ml",
                          covariates = "water",
                          constraint = "none",
                          W          = "queen",
                          rho        = rho_hat),
         timestamp  = Sys.time(),
         ml_history = tuned_ml$history),
    out
  )

  usethis::use_data(results_water_ml, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat, tuned_ml$phi, out$rmse, out$r2, out$coverage_95_obs))
}

# --- ML with rho fixed at 1 ---------------------------------------------------

if (!.should_skip("results_water_ml_rho1")) {
  message("\n== ML tuning: water covariate, rho=0.99 fixed, no RSR ==")

  tuned_ml_rho1 <- fastblm::tune_ml(
    y       = y_alb,
    A       = A,
    Q_fun   = Q_fun_fixed,
    X_fixed = X_obs_water,
    verbose = TRUE
  )

  message(sprintf("  ML optimum: rho=0.99 (fixed)  phi=%.4f  sigma2e=%.4g",
                  tuned_ml_rho1$phi, tuned_ml_rho1$sigma2e))

  out1 <- fit_and_compute(tuned_ml_rho1, rho_val = 0.99)

  results_water_ml_rho1 <- c(
    list(run_name  = "water_ml_rho1",
         tags      = list(tuning     = "ml",
                          covariates = "water",
                          constraint = "none",
                          W          = "queen",
                          rho        = 0.99),
         timestamp  = Sys.time(),
         ml_history = tuned_ml_rho1$history),
    out1
  )

  usethis::use_data(results_water_ml_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  tuned_ml_rho1$phi, out1$rmse, out1$r2, out1$coverage_95_obs))
}

# --- ML no covariates, rho tuned ----------------------------------------------

# Separate fit helper: no X_fixed, so A is spatial-only, posterior_mean is
# length p with no beta block to extract.
fit_and_compute_nocov <- function(tuned, rho_val) {
  fit <- fastblm::fit_fastblm(
    y      = y_alb,
    A      = A,
    Q      = tuned$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )

  mu      <- fit$posterior_mean          # length p
  se      <- fastblm::posterior_se(fit)  # length p

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  tr_t  <- y_latent_true[target_idx]
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
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = NULL,
    se_beta               = NULL,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_val,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

if (!.should_skip("results_nocov_ml")) {
  message("\n== ML tuning: no covariates, rho ML-tuned ==")

  tuned_nocov <- fastblm::tune_ml(
    y          = y_alb,
    A          = A,
    Q_fun      = Q_fun_rho,
    X_fixed    = NULL,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )

  rho_hat_nocov <- tuned_nocov$theta[["rho"]]
  message(sprintf("  ML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_nocov, tuned_nocov$phi, tuned_nocov$sigma2e))

  out_nocov <- fit_and_compute_nocov(tuned_nocov, rho_hat_nocov)

  results_nocov_ml <- c(
    list(run_name  = "nocov_ml",
         tags      = list(tuning     = "ml",
                          covariates = "none",
                          constraint = "none",
                          W          = "queen",
                          rho        = rho_hat_nocov),
         timestamp  = Sys.time(),
         ml_history = tuned_nocov$history),
    out_nocov
  )

  usethis::use_data(results_nocov_ml, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat_nocov, tuned_nocov$phi,
                  out_nocov$rmse, out_nocov$r2, out_nocov$coverage_95_obs))
}

# --- Water covariate, rho ML-tuned, logdet(Q) = 0 ----------------------------

if (!.should_skip("results_water_ml_zerolog")) {
  message("\n== ML tuning: water covariate, rho ML-tuned, logdet(Q)=0 ==")
  message("  (tests whether logdet(Q) singularity near rho=1 distorts phi)")

  tuned_zl <- fastblm::tune_ml(
    y          = y_alb,
    A          = A,
    Q_fun      = Q_fun_rho_zerolog,
    X_fixed    = X_obs_water,
    theta_init = c(rho = 0.9),
    lower      = 0.01,
    upper      = 0.999,
    verbose    = TRUE
  )

  rho_hat_zl <- tuned_zl$theta[["rho"]]
  message(sprintf("  ML optimum: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat_zl, tuned_zl$phi, tuned_zl$sigma2e))

  out_zl <- fit_and_compute(tuned_zl, rho_hat_zl)

  results_water_ml_zerolog <- c(
    list(run_name  = "water_ml_zerolog",
         tags      = list(tuning     = "ml_zerolog",
                          covariates = "water",
                          constraint = "none",
                          W          = "queen",
                          rho        = rho_hat_zl),
         timestamp  = Sys.time(),
         ml_history = tuned_zl$history),
    out_zl
  )

  usethis::use_data(results_water_ml_zerolog, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat_zl, tuned_zl$phi, out_zl$rmse, out_zl$r2,
                  out_zl$coverage_95_obs))
}

# ------------------------------------------------------------------------------

message("\nML sensitivity results saved to data/.")
