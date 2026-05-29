# data-raw/run_sensitivity_softrsr.R
#
# Supplement: Soft RSR penalty sensitivity -- grid search version.
#
# Evaluates cv_mse on a grid of (rho, log_alpha) values with phi profiled
# out at each point. Produces a clean surface for plotting and reporting.
#
# Grid:
#   rho        in {0.70, 0.90, 0.95, 0.99}
#   log_alpha  in {-4, -3, -2, -1, 0, 1, 2, 3, 4}
#
# Total: 36 grid points, phi profiled at each, ~1-2 hours parallel.
#
# Outputs saved to data/ via usethis::use_data():
#   results_softrsr_albedo  -- full results including grid, best fit

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
score_fn         <- fastblm:::.make_score_fn("mse")

# ------------------------------------------------------------------------------
# 3. Q_fun -- fold-aware: takes (theta, A_train)
# ------------------------------------------------------------------------------

Q_fun_soft <- function(theta, A_train) {
  rho        <- plogis(theta[["logit_rho"]])
  alpha      <- exp(theta[["log_alpha"]])
  alpha_safe <- max(alpha, 1e-10)

  S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q_rho <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))

  A_sp_train <- A_train[, seq_len(p), drop = FALSE]
  X_tr       <- A_train[, p + seq_len(q), drop = FALSE]
  C_raw      <- as.matrix(Matrix::crossprod(X_tr, A_sp_train))
  C_scale    <- norm(C_raw, "F")
  C_n        <- if (C_scale > 0) C_raw / C_scale else C_raw

  .Q_rho <- Q_rho; .alpha_safe <- alpha_safe
  .C <- C_n; .p <- p; .q <- q; .lb <- lambda_beta

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
# 4. Grid definition
# ------------------------------------------------------------------------------

rho_grid       <- c(0.70, 0.90, 0.95, 0.99)
log_alpha_grid <- seq(-4, 4, by = 1)
grid           <- expand.grid(rho = rho_grid, log_alpha = log_alpha_grid)
grid$logit_rho <- qlogis(grid$rho)

message(sprintf("\nGrid: %d rho values x %d log_alpha values = %d points",
                length(rho_grid), length(log_alpha_grid), nrow(grid)))

# ------------------------------------------------------------------------------
# 5. Evaluate grid
# ------------------------------------------------------------------------------

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name))
    return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error   = function(e) NULL, warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}

if (!.should_skip("results_softrsr_albedo")) {
  message("\n== Soft RSR grid search ==")
  t_grid <- system.time({
    grid_results <- lapply(seq_len(nrow(grid)), function(i) {
      theta <- c(logit_rho = grid$logit_rho[i], log_alpha = grid$log_alpha[i])

      phi_hat <- fastblm:::.profile_phi_cv(
        y_alb, A_aug, Q_fun_soft, theta, NULL, fold_assignments,
        score_fn, log(0.01), log(1000),
        "pcg", 1e-6, NULL, NULL,
        fold_C_list = NULL, precond_fun = NULL, parallel = TRUE
      )

      cv_mse <- fastblm:::.eval_cv(
        y_alb, A_aug, Q_fun_soft, theta, NULL, phi_hat, fold_assignments,
        score_fn, "pcg", 1e-6, NULL, NULL,
        fold_C_list = NULL, precond_fun = NULL, parallel = TRUE
      )

      message(sprintf("  rho=%.2f  log_alpha=%+.0f  phi=%.2f  cv_mse=%.6e",
                      grid$rho[i], grid$log_alpha[i], phi_hat, cv_mse))

      list(rho       = grid$rho[i],
           log_alpha = grid$log_alpha[i],
           phi       = phi_hat,
           cv_mse    = cv_mse)
    })
  })

  message(sprintf("\nGrid search done in %.1f min", t_grid["elapsed"] / 60))

  grid_df <- do.call(rbind, lapply(grid_results, as.data.frame))
  best_i  <- which.min(grid_df$cv_mse)
  rho_hat   <- grid_df$rho[best_i]
  alpha_hat <- exp(grid_df$log_alpha[best_i])
  phi_opt   <- grid_df$phi[best_i]

  message(sprintf("  Best: rho=%.2f  log_alpha=%.0f  alpha=%.4f  phi=%.2f  cv_mse=%.6e",
                  rho_hat, grid_df$log_alpha[best_i], alpha_hat, phi_opt,
                  grid_df$cv_mse[best_i]))

  # ------------------------------------------------------------------------------
  # 6. Final fit at best grid point
  # ------------------------------------------------------------------------------

  prior_final <- Q_fun_soft(
    c(logit_rho = qlogis(rho_hat), log_alpha = log(alpha_hat)), A_aug
  )
  C_scale_final <- norm(as.matrix(Matrix::crossprod(
    A_aug[, p + seq_len(q), drop = FALSE],
    A_aug[, seq_len(p),     drop = FALSE]
  )), "F")

  AtA_full      <- Matrix::forceSymmetric(Matrix::crossprod(A_aug))
  Q_approx_full <- Matrix::forceSymmetric(
    Matrix::bdiag(prior_final$Q_rho * phi_opt,
                  lambda_beta * phi_opt * Matrix::Diagonal(q))
  ) + AtA_full / 3e-6
  .chol_full    <- Matrix::Cholesky(Q_approx_full, LDL = FALSE, perm = TRUE)
  precond_final <- function(v) as.numeric(Matrix::solve(.chol_full, v))

  fit <- fastblm::fit_fastblm(
    y           = y_alb,
    A           = A_aug,
    Q           = prior_final$Q,
    phi         = phi_opt,
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
    tags                  = list(tuning     = "grid",
                                 covariates = "water",
                                 constraint = "soft_RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 alpha      = alpha_hat,
                                 C_scale    = C_scale_final),
    timestamp             = Sys.time(),
    grid                  = grid_df,
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = phi_opt,
    rho_opt               = rho_hat,
    alpha_opt             = alpha_hat,
    C_scale               = C_scale_final,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_softrsr_albedo, overwrite = TRUE)
  message(sprintf("  rho_opt=%.2f  alpha_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  rho_hat, alpha_hat, rmse, r2,
                  results_softrsr_albedo$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nSoft RSR sensitivity results saved to data/.")
