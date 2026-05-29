# data-raw/run_01_main.R
#
# Main ablation results for the semi-synthetic albedo validation (Table 2)
# and the canonical SIF downscaling result (main paper Figure 4).
#
# Supplement sections covered:
#   - Section 1 (baseline): results_ols_baseline, results_no_cov_rho1
#   - Section 1 (RSR table): results_water_rho1, results_water_rsr_rho1
#   - Section 2 (neighbour baseline): results_water_rsr_rho1
#   - Section 4 (tuning table): results_water_rsr_rho1, results_water_rsr_rho_cv
#   - Main paper: results_sif_canonical
#
# Albedo runs (semi-synthetic, ground truth available):
#   1. OLS baseline           -- coarse-scale regression, no spatial structure
#   2. No covariates, rho=1   -- spatial only, intrinsic SAR prior
#   3. Water, rho=1, no RSR   -- water covariate, no confounding correction
#   4. Water + RSR, rho=1     -- RSR constraint, fixed rho (comparison table)
#   5. Water + RSR, rho CV    -- RSR constraint, rho CV-tuned (best albedo)
#
# SIF run (real data, no ground truth):
#   6. Water + RSR + R_inv, rho CV  -- canonical result, rho CV-tuned
#
# All tuning uses 10-fold CV. Posterior SE via Hutchinson (n_probes=200).
# Outputs saved to data/ via usethis::use_data().

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

A                  <- d_shared$A_flat          # n x p sparse
W_queen            <- d_shared$W_queen         # p x p queen adjacency
X_obs_water        <- d_shared$X_obs_water     # n x 2: [1, water_obs]
X_latent_water     <- d_shared$X_latent_water  # p x 2: [1, water_grid]
fine_grid_buffered <- d_shared$fine_grid_buffered

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

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

# Augmented design matrix: [A | X_obs_water], n x (p+q)
A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

# Q_fun factories: return function(theta) -> list(Q, log_det_Q)

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

# Augmented Q: spatial block + lambda*I_q ridge on beta
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

# Per-fold RSR constraint for augmented system.
# C_aug = [t(X_obs[train,]) %*% A_train | 0_{q x q}]
# Forces spatial component r orthogonal to water covariate; beta unconstrained.
make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    cbind(C_spatial, matrix(0, nrow = q, ncol = q))
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

# Full-data RSR constraint for final fit after tuning
C_aug_full <- cbind(
  as.matrix(t(X_obs_water) %*% A),   # q x p
  matrix(0, nrow = q, ncol = q)      # q x q zero block
)

# ------------------------------------------------------------------------------
# 3. Helpers
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

# Tune phi (and optionally rho) via 10-fold CV, then refit on full data.
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

# Extract posterior summaries and compute albedo metrics (RMSE, R2, coverage).
# For covariate runs: posterior_mean = r + X_latent %*% beta on full grid.
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
# 4. Albedo ablations
# ------------------------------------------------------------------------------

# --- 1. OLS baseline ----------------------------------------------------------
# Coarse-scale OLS: fit on upscaled observations, predict directly to grid.
# No spatial structure. Serves as lower bound in Table 2.

if (!.should_skip("results_ols_baseline")) {
  message("\n== 1. OLS baseline ==")

  lm_fit   <- lm(y_alb ~ X_obs_water[, 2])
  beta_ols  <- coef(lm_fit)
  mu_ols    <- beta_ols[1] + beta_ols[2] * X_latent_water[, 2]

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

# --- 2. No covariates, rho=1 --------------------------------------------------
# Spatial-only model with intrinsic SAR prior. No water covariate.

if (!.should_skip("results_no_cov_rho1")) {
  message("\n== 2. No covariates, rho=1 ==")

  r2   <- tune_and_fit(y_alb, A_fit = A, Q_fun_fit = Q_fun_fixed)
  out2 <- compute_outputs(r2, "no_cov_rho1",
                          tags = list(tuning = "cv", covariates = "none",
                                      constraint = "none", W = "queen", rho = 1))
  results_no_cov_rho1 <- out2
  usethis::use_data(results_no_cov_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out2$phi, out2$rmse, out2$r2, out2$coverage_95_obs))
}

# --- 3. Water covariate, rho=1, no RSR ----------------------------------------
# Adds water covariate but no spatial confounding correction.
# Expected: beta_water positive or near zero (wrong sign -- absorbed by spatial).
# Used in supplement Section 1 naive fit table.

if (!.should_skip("results_water_rho1")) {
  message("\n== 3. Water covariate, rho=1, no RSR ==")

  r3   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug)
  out3 <- compute_outputs(r3, "water_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "none", W = "queen", rho = 1),
                          with_covariates = TRUE)
  results_water_rho1 <- out3
  usethis::use_data(results_water_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out3$phi, out3$rmse, out3$r2,
                  out3$beta_hat[2], out3$coverage_95_obs))
}

# --- 4. Water + RSR, rho=1 ----------------------------------------------------
# RSR constraint enforced, rho held fixed at 1 (intrinsic prior).
# Used in supplement Section 1 RSR table and Section 2 neighbour baseline.
# Fixing rho=1 here keeps the RSR comparison clean (only constraint varies).

if (!.should_skip("results_water_rsr_rho1")) {
  message("\n== 4. Water + RSR, rho=1 ==")

  r4   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_fixed_aug, rsr = TRUE)
  out4 <- compute_outputs(r4, "water_rsr_rho1",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = 1),
                          with_covariates = TRUE)
  results_water_rsr_rho1 <- out4
  usethis::use_data(results_water_rsr_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out4$phi, out4$rmse, out4$r2,
                  out4$beta_hat[2], out4$coverage_95_obs))
}

# --- 5. Water + RSR, rho CV-tuned ---------------------------------------------
# Best albedo model: RSR constraint with rho selected by CV.
# Used in supplement Section 4 tuning comparison table.

if (!.should_skip("results_water_rsr_rho_cv")) {
  message("\n== 5. Water + RSR, rho CV-tuned ==")

  r5   <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                       rsr = TRUE, tune_rho = TRUE)
  out5 <- compute_outputs(r5, "water_rsr_rho_cv",
                          tags = list(tuning = "cv", covariates = "water",
                                      constraint = "RSR", W = "queen", rho = "cv_tuned"),
                          with_covariates = TRUE)
  results_water_rsr_rho_cv <- out5
  usethis::use_data(results_water_rsr_rho_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out5$phi, out5$rho_opt, out5$rmse, out5$r2, out5$coverage_95_obs))
}

# ------------------------------------------------------------------------------
# 5. SIF canonical run
# ------------------------------------------------------------------------------

# --- 5b. Water, no RSR, rho CV-tuned ------------------------------------------
# Fair no-RSR comparison: same rho freedom as results_water_rsr_rho_cv but
# without the RSR constraint. Shows beta_water goes wrong direction without RSR
# even when rho is CV-tuned. Needed for kriging comparison plots.

if (!.should_skip("results_water_rho_cv")) {
  message("\n== 5b. Water, no RSR, rho CV-tuned ==")

  r5b  <- tune_and_fit(y_alb, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                       rsr = FALSE, tune_rho = TRUE)
  out5b <- compute_outputs(r5b, "water_rho_cv",
                           tags = list(tuning = "cv", covariates = "water",
                                       constraint = "none", W = "queen",
                                       rho = "cv_tuned"),
                           with_covariates = TRUE)
  results_water_rho_cv <- out5b
  usethis::use_data(results_water_rho_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  RMSE=%.4f  R2=%.4f  beta_water=%.4f  coverage_obs=%.3f",
                  out5b$phi, out5b$rho_opt, out5b$rmse, out5b$r2,
                  out5b$beta_hat[2], out5b$coverage_95_obs))
}

# --- 6. SIF: Water + RSR + R_inv, rho CV-tuned --------------------------------
# Main application result. Uses per-sounding noise weights (R_inv) from
# SIF_Uncertainty_757nm. RSR constraint corrects for water/spatial confounding.
# rho selected by 10-fold CV (not fixed at 1) for best predictive performance.
# No ground truth: no RMSE/R2/coverage reported.

if (!.should_skip("results_sif_canonical")) {
  message("\n== 6. SIF canonical run: Water + RSR + R_inv, rho CV-tuned ==")

  r6 <- tune_and_fit(y_sif, A_fit = A_aug, Q_fun_fit = Q_fun_rho_aug,
                     rsr = TRUE, R_inv = R_inv_sif, tune_rho = TRUE)

  rho_hat <- r6$tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, r6$tuned$phi, r6$tuned$sigma2e))

  fit6     <- r6$fit
  r_hat6   <- fit6$posterior_mean[seq_len(p)]
  beta6    <- fit6$posterior_mean[p + seq_len(q)]
  mu6      <- r_hat6 + as.numeric(X_latent_water %*% beta6)
  A_pred6  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se6      <- fastblm::posterior_se(fit6, A_new = A_pred6, n_probes = 200L)
  se_beta6 <- fastblm::posterior_se(fit6, n_probes = 200L)[p + seq_len(q)]

  results_sif_canonical <- list(
    run_name              = "sif_canonical",
    tags                  = list(tuning     = "cv",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 R_inv      = "SIF_Uncertainty_757nm"),
    timestamp             = Sys.time(),
    posterior_mean        = mu6[target_idx],
    posterior_se          = se6[target_idx],
    ci_lower              = (mu6 - 1.96 * se6)[target_idx],
    ci_upper              = (mu6 + 1.96 * se6)[target_idx],
    beta_hat              = beta6,
    se_beta               = se_beta6,
    sigma2e               = fit6$sigma2e,
    phi                   = r6$tuned$phi,
    rho_opt               = rho_hat,
    cv_curve              = r6$tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
  usethis::use_data(results_sif_canonical, overwrite = TRUE)
  message(sprintf("  phi=%.4f  rho_opt=%.4f  sigma2e=%.4g  beta=[%.4f, %.4f]",
                  r6$tuned$phi, rho_hat, fit6$sigma2e, beta6[1], beta6[2]))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_01_main.R complete. Objects saved to data/:")
message("  results_ols_baseline")
message("  results_no_cov_rho1")
message("  results_water_rho1")
message("  results_water_rsr_rho1")
message("  results_water_rsr_rho_cv")
message("  results_water_rho_cv")
message("  results_sif_canonical")
