# data-raw/run_sensitivity_softrsr.R
#
# Supplement: Soft RSR penalty sensitivity.
# Instead of hard RSR (exact projection via constrain()), penalizes spatial
# variation in the covariate direction via a rank-q addition to Q:
#
#   Q_soft(v) = Q_sp(rho)(v) + alpha * t(C_n)(C_n v)
#
# where C_n = t(X_obs) %*% A / ||t(X_obs) %*% A||_F  (Frobenius-normalised).
# alpha = 0: unconstrained (= run 3, water no RSR)
# alpha -> Inf: approaches hard RSR (= run 4)
# Both rho and log(alpha) are CV-tuned jointly; phi is profiled out.
#
# Uses plain global scope -- no worker_env -- so future's automatic global
# detection serializes everything correctly to parallel workers.
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

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

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
# 2. Model objects -- plain global scope
# ------------------------------------------------------------------------------

A_aug    <- as(cbind(A, X_obs_water), "dgCMatrix")
C_raw    <- as.matrix(t(X_obs_water) %*% A)
C_scale  <- norm(C_raw, "F")
C_rsr_n  <- C_raw / C_scale
CtC_diag <- as.numeric(Matrix::colSums(C_rsr_n^2))

# Precompute fold Cholesky factors at rho=0.9, phi=82
set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_alb), 10L)
fold_sizes       <- tabulate(fold_assignments)
n_total_obs      <- length(y_alb)

message("Precomputing fold Cholesky factors...")
t_pre <- system.time({
  S_pre <- Matrix::Diagonal(nrow(W_queen)) - 0.9 * W_queen
  Q_pre <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_pre)))

  fold_chol_list <- lapply(seq_len(10L), function(f) {
    train_idx <- which(fold_assignments != f)
    A_f       <- A_aug[train_idx, ]
    AtA_f     <- Matrix::forceSymmetric(Matrix::crossprod(A_f))
    Q_approx  <- Matrix::forceSymmetric(
      Matrix::bdiag(Q_pre * 82, lambda_beta * 82 * Matrix::Diagonal(q))
    ) + AtA_f / 3e-6
    Matrix::Cholesky(Q_approx, LDL = FALSE, perm = TRUE)
  })
})
message(sprintf("  done in %.1fs", t_pre["elapsed"]))

# ------------------------------------------------------------------------------
# 3. Q_fun and precond_fun -- plain functions, globals captured lexically
# ------------------------------------------------------------------------------

Q_fun_soft <- function(theta) {
  rho        <- min(theta[["rho"]], 0.98)
  alpha      <- exp(theta[["log_alpha"]])
  alpha_safe <- max(alpha, 1e-10)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q_rho <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))

  # Capture everything needed by apply_Q_aug in local variables
  # so the returned closure is fully self-contained
  .Q_rho       <- Q_rho
  .alpha_safe  <- alpha_safe
  .C           <- C_rsr_n
  .p           <- p
  .q           <- q
  .lb          <- lambda_beta

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

precond_fun <- function(phi, prior, A_train, y_train) {
  n_train <- nrow(A_train)
  fold_id <- which(n_total_obs - fold_sizes == n_train)[1L]
  if (is.na(fold_id)) {
    message("  WARNING: fold not matched, building from scratch")
    AtA_f    <- Matrix::forceSymmetric(Matrix::crossprod(A_train))
    Q_approx <- Matrix::forceSymmetric(
      Matrix::bdiag(prior$Q_rho * phi, lambda_beta * phi * Matrix::Diagonal(q))
    ) + AtA_f / 3e-6
    .chol <- Matrix::Cholesky(Q_approx, LDL = FALSE, perm = TRUE)
  } else {
    .chol <- fold_chol_list[[fold_id]]
  }
  force(.chol)
  function(v) as.numeric(Matrix::solve(.chol, v))
}

# ------------------------------------------------------------------------------
# 4. Sanity check timing
# ------------------------------------------------------------------------------

message("\nTiming sanity check (sequential, k=10, phi=82)...")
prior_check <- Q_fun_soft(c(rho = 0.9, log_alpha = 0))
t_check <- system.time({
  cv_check <- fastblm:::.eval_cv(
    y_alb, A_aug, prior_check, 82, fold_assignments,
    fastblm:::.make_score_fn("mse"),
    "pcg", 1e-6, NULL, NULL,
    fold_C_list = NULL,
    precond_fun = precond_fun,
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
  message("\n== Soft RSR: rho + alpha CV-tuned, PCG solver ==")

  tuned <- fastblm::tune_cv(
    y           = y_alb,
    A           = A_aug,
    Q_fun       = Q_fun_soft,
    theta_init  = c(rho = 0.9, log_alpha = 0),
    lower       = c(rho = 0.5,  log_alpha = log(1e-6)),
    upper       = c(rho = 0.98, log_alpha = log(1e4)),
    k           = 10L,
    solver      = "pcg",
    constraint  = NULL,
    precond_fun = precond_fun,
    parallel    = TRUE,
    verbose     = TRUE,
    seed        = 2026L
  )

  alpha_hat <- exp(tuned$theta[["log_alpha"]])
  rho_hat   <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  alpha=%.6f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, alpha_hat, tuned$phi, tuned$sigma2e))

  # Final fit with full-data preconditioner
  prior_final   <- Q_fun_soft(c(rho = rho_hat, log_alpha = log(alpha_hat)))
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
                                 C_scale    = C_scale),
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
    C_scale               = C_scale,
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
