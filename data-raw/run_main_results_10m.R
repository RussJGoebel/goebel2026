# data-raw/run_main_results_10m.R
#
# Replicates the canonical albedo ablation (run_01_main.R) but uses a response
# variable derived by direct 10m->coarse aggregation, i.e. soundings_augmented$mean_albedo
# with noise added, rather than the 330m-intermediate version in setup_albedo$y.
#
# The only difference from run_01_main.R is the response y_alb_10m:
#   Standard:  y_alb    = A * y_latent_330m + noise   (330m intermediate step)
#   This file: y_alb_10m = mean_albedo_per_sounding + noise  (direct 10m->coarse)
#
# Runs the same five albedo models as run_01_main.R, plus two additional:
#   1. OLS baseline
#   2. No covariates, rho=1
#   2b. No covariates, rho CV-tuned
#   3. Water, rho=1, no RSR
#   3b. Water, rho CV-tuned, no RSR  <-- used by run_matern_comparison.R
#   4. Water + RSR, rho=1
#   5. Water + RSR, rho CV-tuned  <-- primary comparison with results_water_rsr_rho_cv
#
# Outputs saved to data/:
#   results_10m_ols_baseline
#   results_10m_no_cov_rho1
#   results_10m_no_cov_rho_cv
#   results_10m_water_rho1
#   results_10m_water_rho_cv       <-- new: free-rho, no RSR, for Matern comparison
#   results_10m_water_rsr_rho1
#   results_10m_water_rsr_rho_cv

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- FALSE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water
fine_grid_buffered <- d_shared$fine_grid_buffered

y_latent_true <- d_albedo$y_latent_true
sigma_eps     <- d_albedo$sigma_eps

# Build 10m response: direct 10m->coarse aggregation + same noise level
# soundings_augmented$mean_albedo is exactextractr mean of 10m albedo raster
# over each sounding footprint -- no 330m intermediate step.
keep_idx     <- d_sif$keep_idx   # outlier already removed in setup
mean_alb_10m <- goebel2026::soundings_augmented$mean_albedo[keep_idx]

set.seed(2026L)
y_alb_10m <- mean_alb_10m + rnorm(length(mean_alb_10m), sd = sigma_eps)

message(sprintf("y_alb_10m: n=%d  mean=%.4f  sd=%.4f",
                length(y_alb_10m), mean(y_alb_10m), sd(y_alb_10m)))

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Shared model components (identical to run_01_main.R)
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

make_Q_fun_fixed <- function(W) {
  S <- Matrix::Diagonal(nrow(W)) - W
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
}

make_Q_fun_rho <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
    list(Q = Q, log_det_Q = NULL)
  }
}

make_Q_fun_aug <- function(Q_fun_spatial, q, lambda) {
  force(Q_fun_spatial); force(q); force(lambda)
  function(theta) {
    sp  <- Q_fun_spatial(theta)
    Q_a <- Matrix::drop0(Matrix::forceSymmetric(
      Matrix::bdiag(sp$Q, lambda * Matrix::Diagonal(q))
    ))
    list(Q = Q_a, log_det_Q = NULL)
  }
}

Q_fun_fixed     <- make_Q_fun_fixed(W_queen)
Q_fun_rho       <- make_Q_fun_rho(W_queen)
Q_fun_fixed_aug <- make_Q_fun_aug(Q_fun_fixed, q, lambda_beta)
Q_fun_rho_aug   <- make_Q_fun_aug(Q_fun_rho,   q, lambda_beta)

make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    cbind(C_spatial, matrix(0, nrow = q, ncol = q))
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

C_aug_full <- cbind(
  as.matrix(t(X_obs_water) %*% A),
  matrix(0, nrow = q, ncol = q)
)

# ------------------------------------------------------------------------------
# 3. Helpers (identical to run_01_main.R)
# ------------------------------------------------------------------------------

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name)); return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error   = function(e) NULL,
    warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name)); return(TRUE)
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
  lower      <- if (tune_rho) c(rho = 0.5) else numeric(0)
  upper      <- if (tune_rho) c(rho = 0.999) else numeric(0)
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
    se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
    se_beta  <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]
  } else {
    mu       <- fit$posterior_mean
    se       <- fastblm::posterior_se(fit, n_probes = 200L)
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
# 4. Albedo ablations using 10m response
# ------------------------------------------------------------------------------

# --- 1. OLS baseline ----------------------------------------------------------
if (!.should_skip("results_10m_ols_baseline")) {
  message("\n== 1. OLS baseline (10m) ==")

  lm_fit   <- lm(y_alb_10m ~ X_obs_water[, 2])
  beta_ols  <- coef(lm_fit)
  mu_ols    <- beta_ols[1] + beta_ols[2] * X_latent_water[, 2]

  resid_ols <- mu_ols[target_idx] - y_latent_true[target_idx]
  rmse_ols  <- sqrt(mean(resid_ols^2, na.rm = TRUE))
  ss_res    <- sum(resid_ols^2, na.rm = TRUE)
  ss_tot    <- sum((y_latent_true[target_idx] -
                      mean(y_latent_true[target_idx], na.rm = TRUE))^2, na.rm = TRUE)
  r2_ols    <- 1 - ss_res / ss_tot

  results_10m_ols_baseline <- list(
    run_name              = "10m_ols_baseline",
    tags                  = list(tuning = "ols", covariates = "water",
                                 constraint = "none", W = "none", rho = NA,
                                 response = "10m_direct"),
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
  usethis::use_data(results_10m_ols_baseline, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f", rmse_ols, r2_ols))
}

# --- 2. No covariates, rho=1 --------------------------------------------------
if (!.should_skip("results_10m_no_cov_rho1")) {
  message("\n== 2. No covariates, rho=1 (10m) ==")

  r2   <- tune_and_fit(y_alb_10m, A_fit = A, Q_fun_fit = Q_fun_fixed)
  out2 <- compute_outputs(r2, "10m_no_cov_rho1",
                          tags = list(tuning = "cv", covariates = "none",
                                      constraint = "none", W = "queen", rho = 1,
                                      response = "10m_direct"))
  results_10m_no_cov_rho1 <- out2
  usethis::use_data(results_10m_no_cov_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out2$phi, out2$rmse, out2$r2, out2$coverage_95_obs))
}

# --- 2b. No covariates, rho CV-tuned ------------------------------------------
if (!.should_skip("results_10m_no_cov_rho_cv")) {
  message("\n== 2b. No covariates, rho CV-tuned (10m) ==")

  r2b <- tune_and_fit(
    y         = y_alb_10m,
    A_fit     = A,
    Q_fun_fit = Q_fun_rho,
    tune_rho  = TRUE
  )
  out2b <- compute_outputs(r2b, "10m_no_cov_rho_cv",
                           tags = list(tuning     = "cv",
                                       covariates = "none",
                                       constraint = "none",
                                       W          = "queen",
                                       rho        = "cv_tuned",
                                       response   = "10m_direct"))
  results_10m_no_cov_rho_cv <- out2b
  usethis::use_data(results_10m_no_cov_rho_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out2b$phi, out2b$rho_opt, out2b$rmse, out2b$r2,
                  out2b$coverage_95_obs))
}

# --- 3. Water covariate, rho=1, no RSR ----------------------------------------
if (!.should_skip("results_10m_water_rho1")) {
  message("\n== 3. Water covariate, rho=1, no RSR (10m) ==")

  r3   <- tune_and_fit(y_alb_10m, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug)
  out3 <- compute_outputs(r3, "10m_water_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "none", W = "queen", rho = 1,
                                      response = "10m_direct"),
                          with_covariates = TRUE)
  results_10m_water_rho1 <- out3
  usethis::use_data(results_10m_water_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out3$phi, out3$rmse, out3$r2,
                  out3$beta_hat[2], out3$coverage_95_obs))
}

# --- 3b. Water covariate, rho CV-tuned, no RSR --------------------------------
#
# Used by run_matern_comparison.R to extract rho_opt for the SAR->Matern range
# conversion. No RSR so the free-rho SAR is directly comparable to the Matern
# fit (RSR is a SAR-specific constraint with no natural Matern equivalent).
# ------------------------------------------------------------------------------
if (!.should_skip("results_10m_water_rho_cv")) {
  message("\n== 3b. Water covariate, rho CV-tuned, no RSR (10m) ==")

  r3b  <- tune_and_fit(y_alb_10m, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                       rsr = FALSE, tune_rho = TRUE)
  out3b <- compute_outputs(r3b, "10m_water_rho_cv",
                           tags = list(tuning     = "cv",
                                       covariates = "water",
                                       constraint = "none",
                                       W          = "queen",
                                       rho        = "cv_tuned",
                                       response   = "10m_direct"),
                           with_covariates = TRUE)
  results_10m_water_rho_cv <- out3b
  usethis::use_data(results_10m_water_rho_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out3b$phi, out3b$rho_opt, out3b$rmse, out3b$r2,
                  out3b$beta_hat[2], out3b$coverage_95_obs))
}

# --- 4. Water + RSR, rho=1 ----------------------------------------------------
if (!.should_skip("results_10m_water_rsr_rho1")) {
  message("\n== 4. Water + RSR, rho=1 (10m) ==")

  r4   <- tune_and_fit(y_alb_10m, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug,
                       rsr = TRUE)
  out4 <- compute_outputs(r4, "10m_water_rsr_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = 1,
                                      response = "10m_direct"),
                          with_covariates = TRUE)
  results_10m_water_rsr_rho1 <- out4
  usethis::use_data(results_10m_water_rsr_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out4$phi, out4$rmse, out4$r2,
                  out4$beta_hat[2], out4$coverage_95_obs))
}

# --- 5. Water + RSR, rho CV-tuned ---------------------------------------------
if (!.should_skip("results_10m_water_rsr_rho_cv")) {
  message("\n== 5. Water + RSR, rho CV-tuned (10m) ==")

  r5   <- tune_and_fit(y_alb_10m, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                       rsr = TRUE, tune_rho = TRUE)
  out5 <- compute_outputs(r5, "10m_water_rsr_rho_cv",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen",
                                      rho = "cv_tuned",
                                      response = "10m_direct"),
                          with_covariates = TRUE)
  results_10m_water_rsr_rho_cv <- out5
  usethis::use_data(results_10m_water_rsr_rho_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out5$phi, out5$rho_opt, out5$rmse, out5$r2, out5$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_main_results_10m.R complete. Objects saved to data/:")
message("  results_10m_ols_baseline")
message("  results_10m_no_cov_rho1")
message("  results_10m_no_cov_rho_cv")
message("  results_10m_water_rho1")
message("  results_10m_water_rho_cv")
message("  results_10m_water_rsr_rho1")
message("  results_10m_water_rsr_rho_cv")
