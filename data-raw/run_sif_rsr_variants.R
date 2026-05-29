# data-raw/run_02_sif_rsr_variants.R
#
# Supplement Section 1: Water Covariate and Spatial Confounding.
#
# Produces three results that together tell the RSR story:
#
#   1. results_sif_naive       -- Water, no RSR, rho CV-tuned (SIF)
#                                 Shows beta_water goes wrong direction without RSR.
#
#   2. results_sif_rsr_rho1    -- Water + hard RSR, rho=1 (SIF)
#                                 Fixed rho for clean comparison in the beta table.
#
#   3. results_softrsr_albedo  -- Soft RSR grid search (albedo, ground truth known)
#                                 2D grid over (rho, log_alpha): traces path from
#                                 unconstrained (alpha=0) to hard RSR (alpha->inf).
#                                 Supplement plot shows CV MSE surface; preferred
#                                 alpha is intermediate, validating hard RSR choice.
#
# Note: results_sif_canonical (hard RSR, rho CV-tuned) is produced in run_01_main.R
# and serves as the fourth column of the Section 1 beta table.
#
# Model structure reminder:
#   A_aug = [A | X_obs_water],  n x (p+q)
#   Q_aug = bdiag(Q_spatial, lambda*I_q)
#   Hard RSR: constrain(fit, C_aug_full) where C = [t(X_obs) %*% A | 0]
#   Soft RSR: Q_soft = Q_spatial + alpha * t(C_n) %*% C_n  (C_n = C/||C||_F)
#
# Outputs saved to data/ via usethis::use_data().

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- FALSE

# Runs 1-2 (SIF): multisession works fine (Cholesky solver, no nested parallelism)
# Run 3 (soft RSR grid): uses PCG with parallel=TRUE inside each grid point.
# multicore (fork-based) is required for this to work on Linux -- multisession
# does not support nested parallel calls and will produce Inf cv_mse silently.
# We switch plans around the soft RSR section below.
future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water     # n x 2: [1, water_obs]
X_latent_water     <- d_shared$X_latent_water  # p x 2: [1, water_grid]
fine_grid_buffered <- d_shared$fine_grid_buffered

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

y_sif     <- d_sif$y
R_inv_sif <- d_sif$R_inv

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Shared model components
# ------------------------------------------------------------------------------

A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

# Q_fun factories
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

# Hard RSR constraint helpers
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

# Fit SIF with covariates, extract posterior mean/SE on target pixels.
# No ground truth so no RMSE/R2/coverage.
fit_sif_covariate <- function(tuned, rho_hat) {
  fit      <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    R_inv  = R_inv_sif,
    solver = "cholesky"
  )
  list(fit = fit, tuned = tuned, rho_hat = rho_hat)
}

extract_sif_outputs <- function(fit_obj, run_name, tags, apply_rsr = FALSE) {
  fit   <- fit_obj$fit
  tuned <- fit_obj$tuned

  if (apply_rsr) fit <- fastblm::constrain(fit, C_aug_full)

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

  list(
    run_name              = run_name,
    tags                  = tags,
    timestamp             = Sys.time(),
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    ci_lower              = (mu - 1.96 * se)[target_idx],
    ci_upper              = (mu + 1.96 * se)[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = fit_obj$rho_hat,
    cv_curve              = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

# ------------------------------------------------------------------------------
# 4. SIF naive: water covariate, no RSR, rho CV-tuned
# ------------------------------------------------------------------------------
# Expected: beta_water near zero or positive -- spatial field absorbs the water
# signal when RSR is absent. Contrasts with results_sif_canonical to motivate RSR.

if (!.should_skip("results_sif_naive")) {
  message("\n== 1. SIF naive: water, no RSR, rho CV-tuned ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    R_inv      = R_inv_sif,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    k          = 10L,
    constraint = NULL,
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  fit_obj <- fit_sif_covariate(tuned, rho_hat)
  results_sif_naive <- extract_sif_outputs(
    fit_obj, "sif_naive",
    tags = list(tuning     = "cv",
                response   = "SIF_757nm",
                covariates = "water",
                constraint = "none",
                W          = "queen",
                rho        = rho_hat,
                R_inv      = "SIF_Uncertainty_757nm")
  )
  usethis::use_data(results_sif_naive, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  beta=[%.4f, %.4f]",
                  rho_hat, tuned$phi,
                  results_sif_naive$beta_hat[1],
                  results_sif_naive$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 5. SIF hard RSR, rho=1
# ------------------------------------------------------------------------------
# Fixed rho=1 (intrinsic prior) for the beta comparison table in Section 1.
# Keeping rho fixed isolates the effect of the RSR constraint vs. results_sif_naive.

if (!.should_skip("results_sif_rsr_rho1")) {
  message("\n== 2. SIF hard RSR, rho=1 ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    R_inv      = R_inv_sif,
    theta_init = numeric(0),
    k          = 10L,
    constraint = rsr_constraint,
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))

  fit_obj <- fit_sif_covariate(tuned, rho_hat = 1)
  results_sif_rsr_rho1 <- extract_sif_outputs(
    fit_obj, "sif_rsr_rho1",
    tags = list(tuning     = "cv",
                response   = "SIF_757nm",
                covariates = "water",
                constraint = "hard_RSR",
                W          = "queen",
                rho        = 1,
                R_inv      = "SIF_Uncertainty_757nm"),
    apply_rsr = TRUE
  )
  usethis::use_data(results_sif_rsr_rho1, overwrite = TRUE)
  message(sprintf("  phi=%.4f  beta=[%.4f, %.4f]",
                  tuned$phi,
                  results_sif_rsr_rho1$beta_hat[1],
                  results_sif_rsr_rho1$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 6. Soft RSR grid search (albedo)
# ------------------------------------------------------------------------------
# Grid search over (rho, log_alpha) using albedo where ground truth is known.
# Q_soft replaces the hard constraint with a penalty: Q_sp + alpha * t(C_n) %*% C_n
# where C_n is the per-fold RSR matrix normalized by its Frobenius norm.
# Normalization is fold-aware: C_n computed from A_train at each fold so the
# constraint only sees training data.
# At each grid point: tune phi by CV (PCG solver), record cv_mse and phi.
# Supplement plot: cv_mse surface showing preferred alpha is intermediate
# between unconstrained (alpha->0) and hard RSR (alpha->inf).
#
# Grid: rho in {0.70, 0.90, 0.95, 0.99} x log_alpha in seq(-4, 4, by=1) = 36 pts
# Expected runtime: ~1-2 hours parallel.

if (!.should_skip("results_softrsr_albedo")) {
  message("\n== 4. Soft RSR grid search (albedo) ==")

  # Grid
  rho_grid       <- c(0.70, 0.90, 0.95, 0.99)
  log_alpha_grid <- seq(-4, 4, by = 1)
  grid_params    <- expand.grid(rho = rho_grid, log_alpha = log_alpha_grid)
  grid_params$logit_rho <- qlogis(grid_params$rho)
  message(sprintf("  Grid: %d rho x %d log_alpha = %d points",
                  length(rho_grid), length(log_alpha_grid), nrow(grid_params)))

  # Fold-aware soft RSR Q_fun: takes (theta, A_train).
  # C_n is recomputed from A_train at each fold call so constraint only
  # sees training data. rho is parameterized on logit scale for stability.
  Q_fun_soft <- function(theta, A_train = NULL) {
    rho        <- plogis(theta[["logit_rho"]])
    alpha      <- exp(theta[["log_alpha"]])
    alpha_safe <- max(alpha, 1e-10)

    S_rho <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
    Q_rho <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))

    # Compute fold-specific normalized constraint matrix
    if (!is.null(A_train)) {
      A_sp_tr <- A_train[, seq_len(p), drop = FALSE]
      X_tr    <- A_train[, p + seq_len(q), drop = FALSE]
      C_raw   <- as.matrix(Matrix::crossprod(X_tr, A_sp_tr))
    } else {
      C_raw <- as.matrix(t(X_obs_water) %*% A)
    }
    C_scale <- norm(C_raw, "F")
    C_n     <- if (C_scale > 0) C_raw / C_scale else C_raw

    .Q_rho <- Q_rho; .alpha_safe <- alpha_safe
    .C_n <- C_n; .p <- p; .q <- q; .lb <- lambda_beta

    apply_Q_aug <- function(v) {
      v_sp   <- v[seq_len(.p)]
      v_beta <- v[.p + seq_len(.q)]
      Cv     <- as.numeric(.C_n %*% v_sp)
      CtCv   <- as.numeric(t(.C_n) %*% Cv)
      c(as.numeric(.Q_rho %*% v_sp) + .alpha_safe * CtCv,
        .lb * v_beta)
    }

    list(Q = apply_Q_aug, log_det_Q = NULL, Q_rho = Q_rho)
  }

  # Preconditioner: sparse approximation ignoring soft penalty term
  make_precond <- function(phi, prior, A_train, y_train) {
    AtA_f    <- Matrix::forceSymmetric(Matrix::crossprod(A_train))
    Q_approx <- Matrix::forceSymmetric(
      Matrix::bdiag(prior$Q_rho * phi, lambda_beta * phi * Matrix::Diagonal(q))
    ) + AtA_f / 3e-6
    ch <- Matrix::Cholesky(Q_approx, LDL = FALSE, perm = TRUE)
    force(ch)
    function(v) as.numeric(Matrix::solve(ch, v))
  }

  message("  Running grid search (sequential, one line per point)...")
  t_grid <- system.time({
    grid_results <- lapply(seq_len(nrow(grid_params)), function(i) {
      rho       <- grid_params$rho[i]
      log_alpha <- grid_params$log_alpha[i]
      logit_rho <- grid_params$logit_rho[i]

      Q_fun_fixed_point <- function(theta_empty, A_train = NULL) {
        Q_fun_soft(c(logit_rho = logit_rho, log_alpha = log_alpha), A_train)
      }

      tuned <- tryCatch(
        fastblm::tune_cv(
          y           = y_alb,
          A           = A_aug,
          Q_fun       = Q_fun_fixed_point,
          R_inv       = NULL,
          theta_init  = numeric(0),
          k           = 10L,
          solver      = "pcg",
          pcg_maxit   = 500L,        # cap so bad points fail fast
          precond_fun = make_precond,
          seed        = 2026L,
          parallel    = FALSE,
          verbose     = FALSE
        ),
        error = function(e) list(phi = NA_real_, value = Inf)
      )

      message(sprintf("  rho=%.2f  log_alpha=%+.0f  phi=%.2f  cv_mse=%.6e",
                      rho, log_alpha, tuned$phi, tuned$value))

      data.frame(rho       = rho,
                 log_alpha = log_alpha,
                 phi       = tuned$phi,
                 cv_mse    = tuned$value)
    })
  })

  message(sprintf("  Grid search done in %.1f min", t_grid["elapsed"] / 60))

  grid_df  <- do.call(rbind, grid_results)
  best_idx <- which.min(grid_df$cv_mse)
  best     <- grid_df[best_idx, ]
  message(sprintf("  Best: rho=%.2f  log_alpha=%.0f  alpha=%.4f  phi=%.2f  cv_mse=%.6e",
                  best$rho, best$log_alpha, exp(best$log_alpha),
                  best$phi, best$cv_mse))

  # Final fit at best grid point using PCG with preconditioner
  theta_best  <- c(logit_rho = qlogis(best$rho), log_alpha = best$log_alpha)
  prior_final <- Q_fun_soft(theta_best, A_train = NULL)

  # Build preconditioner on full data
  AtA_full      <- Matrix::forceSymmetric(Matrix::crossprod(A_aug))
  Q_approx_full <- Matrix::forceSymmetric(
    Matrix::bdiag(prior_final$Q_rho * best$phi,
                  lambda_beta * best$phi * Matrix::Diagonal(q))
  ) + AtA_full / 3e-6
  ch_full       <- Matrix::Cholesky(Q_approx_full, LDL = FALSE, perm = TRUE)
  precond_final <- function(v) as.numeric(Matrix::solve(ch_full, v))

  fit_best <- fastblm::fit_fastblm(
    y           = y_alb,
    A           = A_aug,
    Q           = prior_final$Q,
    phi         = best$phi,
    R_inv       = NULL,
    solver      = "pcg",
    pcg_precond = precond_final
  )

  r_hat_b    <- fit_best$posterior_mean[seq_len(p)]
  beta_hat_b <- fit_best$posterior_mean[p + seq_len(q)]
  mu_b       <- r_hat_b + as.numeric(X_latent_water %*% beta_hat_b)
  A_pred_b   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se_b       <- fastblm::posterior_se(fit_best, A_new = A_pred_b, n_probes = 200L)

  mu_t  <- mu_b[target_idx]
  se_t  <- se_b[target_idx]
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t
  ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }

  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

  # C_scale at best point (full data, for record-keeping)
  C_raw_full   <- as.matrix(t(X_obs_water) %*% A)
  C_scale_full <- norm(C_raw_full, "F")

  results_softrsr_albedo <- list(
    run_name              = "softrsr_albedo",
    tags                  = list(tuning     = "cv_grid",
                                 response   = "albedo",
                                 covariates = "water",
                                 constraint = "soft_RSR",
                                 W          = "queen",
                                 rho        = best$rho,
                                 alpha      = exp(best$log_alpha),
                                 C_scale    = C_scale_full),
    timestamp             = Sys.time(),
    grid                  = grid_df,       # cols: rho, log_alpha, phi, cv_mse
    rho_opt               = best$rho,
    log_alpha_opt         = best$log_alpha,
    alpha_opt             = exp(best$log_alpha),
    phi                   = best$phi,
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat_b,
    se_beta               = fastblm::posterior_se(fit_best, n_probes = 200L)[p + seq_len(q)],
    sigma2e               = fit_best$sigma2e,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
  usethis::use_data(results_softrsr_albedo, overwrite = TRUE)
  message(sprintf("  rho_opt=%.2f  alpha_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  best$rho, exp(best$log_alpha), rmse, r2,
                  results_softrsr_albedo$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_02_sif_rsr_variants.R complete. Objects saved to data/:")
message("  results_sif_naive")
message("  results_sif_rsr_rho1")
message("  results_softrsr_albedo")
