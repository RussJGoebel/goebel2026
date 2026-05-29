# data-raw/run_sensitivity_softrsr.R
#
# Supplement: Soft RSR penalty sensitivity.
# Instead of hard RSR (exact projection via constrain()), penalizes spatial
# variation in the covariate direction via a rank-q addition to Q:
#
#   Q_soft(v) = Q_sp(rho)(v) + alpha * t(C_n)(C_n v)
#
# where C_n = t(X_obs_train) %*% A_train / ||t(X_obs_train) %*% A_train||_F
# is computed from the TRAINING fold only, so the penalty direction is
# correctly specified for each fold's training problem.
#
# alpha = 0: unconstrained (= run 3, water no RSR)
# alpha -> Inf: approaches hard RSR (= run 4)
# Both rho and log(alpha) are CV-tuned jointly; phi is profiled out.
# No preconditioner -- PCG converges in ~22 iterations with relative tolerance.
#
# Q_fun_soft takes two arguments (theta, A_train), making it fold-aware.
# tune_cv detects this via length(formals(Q_fun)) >= 2 and calls it per fold.
#
# Baseline comparisons:
#   results_water_rho1      -- no RSR (alpha = 0)
#   results_water_rsr_rho1  -- hard RSR
#
# Outputs saved to data/ via usethis::use_data():
#   results_softrsr_albedo  -- soft RSR, alpha CV-tuned

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multicore, workers = parallel::detectCores() - 1L)

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
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Model objects
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_alb), 10L)

# ------------------------------------------------------------------------------
# 3. Q_fun -- fold-aware: takes (theta, A_train)
#
# C is computed from the training rows of A_aug only, so the soft RSR penalty
# direction is correctly specified for each fold. The Frobenius normalisation
# is also per-fold so alpha is on a consistent scale across folds.
# Q_rho depends only on theta, not on the fold.
# ------------------------------------------------------------------------------

Q_fun_soft <- function(theta, A_train) {
  rho        <- plogis(theta[["logit_rho"]])   # logit_rho = log(rho/(1-rho)), rho in (0,1)
  alpha      <- exp(theta[["log_alpha"]])
  alpha_safe <- max(alpha, 1e-10)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q_rho <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))

  # Fold-specific C from training data only
  A_sp_train <- A_train[, seq_len(p), drop = FALSE]
  X_tr       <- A_train[, p + seq_len(q), drop = FALSE]
  C_raw      <- as.matrix(Matrix::crossprod(X_tr, A_sp_train))  # q x p
  C_scale    <- norm(C_raw, "F")
  C_n        <- if (C_scale > 0) C_raw / C_scale else C_raw

  # Capture locals for closure
  .Q_rho      <- Q_rho
  .alpha_safe <- alpha_safe
  .C          <- C_n
  .p          <- p
  .q          <- q
  .lb         <- lambda_beta

  apply_Q_aug <- function(v) {
    v_sp   <- v[seq_len(.p)]
    v_beta <- v[.p + seq_len(.q)]
    Cv     <- as.numeric(.C %*% v_sp)
    CtCv   <- as.numeric(t(.C) %*% Cv)
    c(as.numeric(.Q_rho %*% v_sp) + .alpha_safe * CtCv,
      .lb * v_beta)
  }

  list(Q = apply_Q_aug, log_det_Q = NULL, Q_rho = Q_rho, precond = NULL)
}

# ------------------------------------------------------------------------------
# 4. Sanity check timing
# ------------------------------------------------------------------------------

message("\nTiming sanity check (sequential, k=10, phi=82, fold-aware Q_fun)...")

t_check <- system.time({
  cv_check <- fastblm:::.eval_cv(
    y_alb, A_aug,
    Q_fun_soft, c(logit_rho = qlogis(0.9), log_alpha = 0), NULL,
    82, fold_assignments,
    fastblm:::.make_score_fn("mse"),
    "pcg", 1e-6, NULL, NULL,
    fold_C_list = NULL,
    precond_fun = NULL,
    parallel    = FALSE
  )
})
message(sprintf("  %.2fs total  %.2fs/fold  cv_mse=%.4e",
                t_check["elapsed"], t_check["elapsed"] / 10, cv_check))

# ------------------------------------------------------------------------------
# 5. Run
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

if (!.should_skip("results_softrsr_albedo")) {
  message("\n== Soft RSR: rho + alpha CV-tuned, fold-aware Q_fun, no preconditioner ==")

  tuned <- fastblm::tune_cv(
    y           = y_alb,
    A           = A_aug,
    Q_fun       = Q_fun_soft,
    theta_init  = c(logit_rho = qlogis(0.75), log_alpha = 1),
    lower       = c(logit_rho = qlogis(0.5),  log_alpha = log(1e-6)),
    upper       = c(logit_rho = qlogis(0.999), log_alpha = log(1e4)),
    k           = 10L,
    solver      = "pcg",
    constraint      = NULL,
    precond_fun     = NULL,
    optim_method    = "coordinate",
    lbfgsb_control  = list(ndeps = c(0.1, 0.5), factr = 1e3),
    parallel        = TRUE,
    verbose         = TRUE,
    folds           = fold_assignments
  )

  alpha_hat <- exp(tuned$theta[["log_alpha"]])
  rho_hat   <- plogis(tuned$theta[["logit_rho"]])
  message(sprintf("  tuning done: rho=%.4f  alpha=%.6f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, alpha_hat, tuned$phi, tuned$sigma2e))

  # Final fit using full A_aug as A_train
  prior_final <- Q_fun_soft(c(logit_rho = qlogis(rho_hat), log_alpha = log(alpha_hat)), A_aug)
  C_scale_final <- norm(as.matrix(Matrix::crossprod(
    A_aug[, p + seq_len(q), drop = FALSE],
    A_aug[, seq_len(p),     drop = FALSE]
  )), "F")

  # Full-data preconditioner at tuned phi
  AtA_full      <- Matrix::forceSymmetric(Matrix::crossprod(A_aug))
  Q_approx_full <- Matrix::forceSymmetric(
    Matrix::bdiag(prior_final$Q_rho * tuned$phi,
                  lambda_beta * tuned$phi * Matrix::Diagonal(q))
  ) + AtA_full / 3e-6
  .chol_full    <- Matrix::Cholesky(Q_approx_full, LDL = FALSE, perm = TRUE)
  precond_final <- function(v) as.numeric(Matrix::solve(.chol_full, v))

  fit <- fastblm::fit_fastblm(
    y           = y_alb,
    A           = A_aug,
    Q           = prior_final$Q,
    phi         = tuned$phi,
    solver      = "pcg",
    pcg_precond = precond_final
  )

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)

  A_pred  <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se      <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

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

  results_softrsr_albedo <- list(
    run_name              = "softrsr_albedo",
    tags                  = list(tuning     = "cv",
                                 covariates = "water",
                                 constraint = "soft_RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 alpha      = alpha_hat,
                                 C_scale    = C_scale_final),
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    alpha_opt             = alpha_hat,
    C_scale               = C_scale_final,
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_softrsr_albedo, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  alpha_opt=%.6f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat, alpha_hat, rmse, r2,
                  results_softrsr_albedo$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSoft RSR sensitivity results saved to data/.")
