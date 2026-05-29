# data-raw/run_04_tuning.R
#
# Supplement Section 4: Hyperparameter Tuning.
#
# Produces results comparing random vs blocked CV, and CV vs ML tuning.
# All albedo runs have ground truth; SIF runs do not.
#
# Albedo:
#   1. results_cv_blocked       -- blocked CV, water+RSR, rho=1
#   2. results_water_ml         -- ML, water, rho free
#
# SIF:
#   3. results_sif_cv_blocked   -- blocked CV, water+RSR+R_inv, rho=1
#   4. results_sif_water_ml     -- ML, water+R_inv, rho free
#   5. results_sif_rho095_ml    -- ML, water+R_inv, rho=0.95 fixed
#   6. results_sif_rho095_cv    -- CV, water+RSR+R_inv, rho=0.95 fixed
#
# Baselines from run_01 (for comparison in supplement table):
#   results_water_rsr_rho1      -- albedo random CV, rho=1
#   results_sif_canonical       -- SIF random CV, rho free
#
# Key findings expected:
#   - Blocked CV selects much smaller phi than random CV (OCO-2 overlap effect)
#   - ML drifts toward rho~1 boundary (logdet(Q) distortion)
#   - At rho=0.95, ML and CV agree well
#
# Outputs saved to data/ via usethis::use_data().

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

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true
blocked_folds_alb <- d_albedo$blocked_folds   # pre-computed spatial blocks

y_sif     <- d_sif$y
R_inv_sif <- d_sif$R_inv
blocked_folds_sif <- d_sif$blocked_folds

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(A)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(A > 0))

# ------------------------------------------------------------------------------
# 2. Shared model components
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

Q_fun_fixed_aug <- make_Q_fun_aug(make_Q_fun_fixed(W_queen), q, lambda_beta)
Q_fun_rho_aug   <- make_Q_fun_aug(make_Q_fun_rho(W_queen),   q, lambda_beta)

# Fixed rho=0.95 Q_fun
make_Q_fun_fixed_rho <- function(W, rho_fixed) {
  S <- Matrix::Diagonal(nrow(W)) - rho_fixed * W
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
}

Q_fun_rho095_aug <- make_Q_fun_aug(
  make_Q_fun_fixed_rho(W_queen, 0.95), q, lambda_beta
)

# RSR constraint
rsr_constraint <- function(train_idx, A_aug_train) {
  A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
  C_spatial <- as.matrix(t(X_obs_water[train_idx, , drop = FALSE]) %*% A_train)
  cbind(C_spatial, matrix(0, nrow = q, ncol = q))
}

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
    error = function(e) NULL, warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}

# Albedo outputs: includes RMSE/R2/coverage
compute_albedo_outputs <- function(fit, tuned, run_name, tags, rho_opt = 1) {
  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

  mu_t  <- mu[target_idx]; se_t <- se[target_idx]
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t; ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }

  resid <- mu_t - tr_t
  rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
  r2    <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

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
    rho_opt               = rho_opt,
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

# SIF outputs: no ground truth
compute_sif_outputs <- function(fit, tuned, run_name, tags, rho_opt) {
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
    rho_opt               = rho_opt,
    cv_curve              = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

# ------------------------------------------------------------------------------
# 4. Albedo: blocked CV, water + RSR, rho=1
# ------------------------------------------------------------------------------
# Uses pre-computed k-means spatial blocks from setup_albedo$blocked_folds.
# Expected: phi much smaller than random CV (~7 vs ~82) due to OCO-2 overlap
# structure -- blocked folds create sparser train/test splits.

if (!.should_skip("results_cv_blocked")) {
  message("\n== 1. Albedo blocked CV, water + RSR, rho=1 ==")
  message(sprintf("  Blocked folds: k=%d", max(blocked_folds_alb)))

  tuned <- fastblm::tune_cv(
    y          = y_alb,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    theta_init = numeric(0),
    k          = max(blocked_folds_alb),
    folds      = blocked_folds_alb,
    constraint = rsr_constraint,
    parallel   = TRUE,
    verbose    = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_alb, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_cv_blocked <- compute_albedo_outputs(
    fit, tuned, "cv_blocked",
    tags = list(tuning = "cv_blocked", covariates = "water",
                constraint = "RSR", W = "queen", rho = 1)
  )
  usethis::use_data(results_cv_blocked, overwrite = TRUE)
  message(sprintf("  phi=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  tuned$phi, results_cv_blocked$rmse,
                  results_cv_blocked$r2, results_cv_blocked$coverage_95_obs))
}

# ------------------------------------------------------------------------------
# 5. Albedo: ML, water, rho free
# ------------------------------------------------------------------------------
# ML maximizes marginal likelihood over rho and phi jointly.
# Expected: ML wants more smoothing than CV (larger phi), rho near boundary.

if (!.should_skip("results_water_ml")) {
  message("\n== 2. Albedo ML, water, rho free ==")

  tuned <- fastblm::tune_ml(
    y          = y_alb,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    X_fixed    = NULL,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_alb, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_water_ml <- compute_albedo_outputs(
    fit, tuned, "water_ml",
    tags = list(tuning = "ml", covariates = "water",
                constraint = "RSR", W = "queen", rho = rho_hat),
    rho_opt = rho_hat
  )
  usethis::use_data(results_water_ml, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  RMSE=%.4f  R2=%.4f",
                  rho_hat, tuned$phi,
                  results_water_ml$rmse, results_water_ml$r2))
}

# ------------------------------------------------------------------------------
# 6. SIF: blocked CV, water + RSR + R_inv, rho=1
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_cv_blocked")) {
  message("\n== 3. SIF blocked CV, water + RSR + R_inv, rho=1 ==")
  message(sprintf("  Blocked folds: k=%d", max(blocked_folds_sif)))

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_fixed_aug,
    R_inv      = R_inv_sif,
    theta_init = numeric(0),
    k          = max(blocked_folds_sif),
    folds      = blocked_folds_sif,
    constraint = rsr_constraint,
    parallel   = TRUE,
    verbose    = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, R_inv = R_inv_sif, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_sif_cv_blocked <- compute_sif_outputs(
    fit, tuned, "sif_cv_blocked",
    tags = list(tuning = "cv_blocked", response = "SIF_757nm",
                covariates = "water", constraint = "RSR",
                W = "queen", rho = 1, R_inv = "SIF_Uncertainty_757nm"),
    rho_opt = 1
  )
  usethis::use_data(results_sif_cv_blocked, overwrite = TRUE)
  message(sprintf("  phi=%.6f  beta=[%.4f, %.4f]",
                  tuned$phi,
                  results_sif_cv_blocked$beta_hat[1],
                  results_sif_cv_blocked$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 7. SIF: ML, water + R_inv, rho free
# ------------------------------------------------------------------------------
# Expected: ML drifts toward rho~0.997, phi very small -- boundary degeneracy.

if (!.should_skip("results_sif_water_ml")) {
  message("\n== 4. SIF ML, water + R_inv, rho free ==")

  tuned <- fastblm::tune_ml(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    X_fixed    = NULL,
    R_inv      = R_inv_sif,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, R_inv = R_inv_sif, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_sif_water_ml <- compute_sif_outputs(
    fit, tuned, "sif_water_ml",
    tags = list(tuning = "ml", response = "SIF_757nm",
                covariates = "water", constraint = "RSR",
                W = "queen", rho = rho_hat,
                R_inv = "SIF_Uncertainty_757nm"),
    rho_opt = rho_hat
  )
  usethis::use_data(results_sif_water_ml, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  beta=[%.4f, %.4f]",
                  rho_hat, tuned$phi,
                  results_sif_water_ml$beta_hat[1],
                  results_sif_water_ml$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 8. SIF: ML, rho=0.95 fixed
# ------------------------------------------------------------------------------
# At rho=0.95 ML and CV should agree -- used to validate CV selection.

if (!.should_skip("results_sif_rho095_ml")) {
  message("\n== 5. SIF ML, rho=0.95 fixed ==")

  tuned <- fastblm::tune_ml(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho095_aug,
    X_fixed    = NULL,
    R_inv      = R_inv_sif,
    theta_init = numeric(0),
    verbose    = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, R_inv = R_inv_sif, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_sif_rho095_ml <- compute_sif_outputs(
    fit, tuned, "sif_rho095_ml",
    tags = list(tuning = "ml", response = "SIF_757nm",
                covariates = "water", constraint = "RSR",
                W = "queen", rho = 0.95,
                R_inv = "SIF_Uncertainty_757nm"),
    rho_opt = 0.95
  )
  usethis::use_data(results_sif_rho095_ml, overwrite = TRUE)
  message(sprintf("  phi=%.4f  beta=[%.4f, %.4f]",
                  tuned$phi,
                  results_sif_rho095_ml$beta_hat[1],
                  results_sif_rho095_ml$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 9. SIF: CV, rho=0.95 fixed
# ------------------------------------------------------------------------------
# Compare to ML at same rho -- should give similar phi if ML is behaving well.

if (!.should_skip("results_sif_rho095_cv")) {
  message("\n== 6. SIF CV, rho=0.95 fixed, water + RSR + R_inv ==")

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho095_aug,
    R_inv      = R_inv_sif,
    theta_init = numeric(0),
    k          = 10L,
    constraint = rsr_constraint,
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g", tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_sif, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, R_inv = R_inv_sif, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  results_sif_rho095_cv <- compute_sif_outputs(
    fit, tuned, "sif_rho095_cv",
    tags = list(tuning = "cv", response = "SIF_757nm",
                covariates = "water", constraint = "RSR",
                W = "queen", rho = 0.95,
                R_inv = "SIF_Uncertainty_757nm"),
    rho_opt = 0.95
  )
  usethis::use_data(results_sif_rho095_cv, overwrite = TRUE)
  message(sprintf("  phi=%.4f  beta=[%.4f, %.4f]",
                  tuned$phi,
                  results_sif_rho095_cv$beta_hat[1],
                  results_sif_rho095_cv$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------------------------

message("\n=== Section 4 tuning summary ===")

message("\nAlbedo (random CV baseline: results_water_rsr_rho1):")
for (nm in c("results_water_rsr_rho1", "results_cv_blocked", "results_water_ml")) {
  r <- tryCatch({
    e <- new.env(); data(list = nm, package = "goebel2026", envir = e); e[[nm]]
  }, error = function(e) NULL)
  if (!is.null(r))
    message(sprintf("  %-25s  rho=%.3f  phi=%7.3f  RMSE=%.4f  R2=%.4f",
                    nm, r$rho_opt, r$phi, r$rmse, r$r2))
}

message("\nSIF (random CV baseline: results_sif_canonical):")
for (nm in c("results_sif_canonical", "results_sif_cv_blocked",
             "results_sif_water_ml", "results_sif_rho095_ml",
             "results_sif_rho095_cv")) {
  r <- tryCatch({
    e <- new.env(); data(list = nm, package = "goebel2026", envir = e); e[[nm]]
  }, error = function(e) NULL)
  if (!is.null(r))
    message(sprintf("  %-30s  rho=%.3f  phi=%.4f  beta_water=%.4f",
                    nm, r$rho_opt, r$phi, r$beta_hat[2]))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_04_tuning.R complete. Objects saved to data/:")
message("  results_cv_blocked")
message("  results_water_ml")
message("  results_sif_cv_blocked")
message("  results_sif_water_ml")
message("  results_sif_rho095_ml")
message("  results_sif_rho095_cv")
