# data-raw/run_05_gA_cv.R
#
# Forward operator sensitivity: g-weighted A matrices with CV rho selection.
#
# The canonical SIF run (run_01) uses uniform A (aij = intersection area /
# footprint area) with rho CV-tuned to ~0.95. Here we ask: if A more
# accurately reflects the sensor response function g (a Gaussian-like
# weighting that down-weights footprint edges), does the CV-selected rho
# change? The hypothesis is that smoother A absorbs spatial structure that
# the SAR prior was compensating for, pushing rho toward smaller values.
#
# Setup: setup_g_A contains A_g (list: tau_0.2, tau_0.333, tau_0.5) where
# tau controls the spread of g, and A_uniform (the standard uniform A).
# Larger tau = more peaked g = sharper footprint weighting.
# Smaller tau = flatter g = smoother, more uniform-like weighting.
#
# All runs: SIF, water covariate + RSR + R_inv, rho CV-tuned.
# Compare rho_opt, phi, beta_water across tau values and vs canonical.
#
# Outputs saved to data/ via usethis::use_data():
#   results_sif_gA_tau02    -- tau=0.2  (smoothest g)
#   results_sif_gA_tau033   -- tau=0.333
#   results_sif_gA_tau05    -- tau=0.5  (sharpest g)

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
d_sif    <- goebel2026::setup_sif
d_gA     <- goebel2026::setup_g_A

fine_grid_buffered <- d_shared$fine_grid_buffered
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_sif     <- d_sif$y
R_inv_sif <- d_sif$R_inv

# g-weighted A matrices
A_list <- list(
  tau_0.2   = d_gA$A_g$tau_0.2,
  tau_0.333 = d_gA$A_g$tau_0.333,
  tau_0.5   = d_gA$A_g$tau_0.5
)

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(d_gA$A_g$tau_0.2)   # same p for all A
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(d_shared$A_flat > 0))

# ------------------------------------------------------------------------------
# 2. Shared model components
# ------------------------------------------------------------------------------

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

Q_fun_rho_aug <- make_Q_fun_aug(make_Q_fun_rho(W_queen), q, lambda_beta)

# Per-fold RSR constraint -- note: uses A_aug from the specific g-weighted A
make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    cbind(C_spatial, matrix(0, nrow = q, ncol = q))
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

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
    error   = function(e) NULL, warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already saved in data/)", obj_name))
    return(TRUE)
  }
  FALSE
}

fit_sif_gA <- function(A, obj_name, tau_val) {
  A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")

  # Full-data RSR constraint for final fit
  C_aug_full <- cbind(
    as.matrix(t(X_obs_water) %*% A),
    matrix(0, nrow = q, ncol = q)
  )

  tuned <- fastblm::tune_cv(
    y          = y_sif,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    R_inv      = R_inv_sif,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    k          = 10L,
    constraint = rsr_constraint,
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y      = y_sif,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    R_inv  = R_inv_sif,
    solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

  message(sprintf("  beta: intercept=%.4f  water=%.4f",
                  beta_hat[1], beta_hat[2]))

  list(
    run_name              = obj_name,
    tags                  = list(tuning     = "cv",
                                 response   = "SIF_757nm",
                                 covariates = "water",
                                 constraint = "RSR",
                                 W          = "queen",
                                 rho        = rho_hat,
                                 R_inv      = "SIF_Uncertainty_757nm",
                                 A          = "g_weighted",
                                 tau        = tau_val),
    timestamp             = Sys.time(),
    posterior_mean        = mu[target_idx],
    posterior_se          = se[target_idx],
    ci_lower              = (mu - 1.96 * se)[target_idx],
    ci_upper              = (mu + 1.96 * se)[target_idx],
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    cv_curve              = tuned$history,
    n_soundings_per_pixel = n_soundings_per_pixel[target_idx]
  )
}

# ------------------------------------------------------------------------------
# 4. Runs
# ------------------------------------------------------------------------------

if (!.should_skip("results_sif_gA_tau02")) {
  message("\n== 1. g-weighted A, tau=0.2 (smoothest) ==")
  results_sif_gA_tau02 <- fit_sif_gA(A_list$tau_0.2, "results_sif_gA_tau02", 0.2)
  usethis::use_data(results_sif_gA_tau02, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  beta_water=%.4f",
                  results_sif_gA_tau02$rho_opt,
                  results_sif_gA_tau02$phi,
                  results_sif_gA_tau02$beta_hat[2]))
}

if (!.should_skip("results_sif_gA_tau033")) {
  message("\n== 2. g-weighted A, tau=0.333 ==")
  results_sif_gA_tau033 <- fit_sif_gA(A_list$tau_0.333, "results_sif_gA_tau033", 0.333)
  usethis::use_data(results_sif_gA_tau033, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  beta_water=%.4f",
                  results_sif_gA_tau033$rho_opt,
                  results_sif_gA_tau033$phi,
                  results_sif_gA_tau033$beta_hat[2]))
}

if (!.should_skip("results_sif_gA_tau05")) {
  message("\n== 3. g-weighted A, tau=0.5 (sharpest) ==")
  results_sif_gA_tau05 <- fit_sif_gA(A_list$tau_0.5, "results_sif_gA_tau05", 0.5)
  usethis::use_data(results_sif_gA_tau05, overwrite = TRUE)
  message(sprintf("  rho_opt=%.4f  phi=%.4f  beta_water=%.4f",
                  results_sif_gA_tau05$rho_opt,
                  results_sif_gA_tau05$phi,
                  results_sif_gA_tau05$beta_hat[2]))
}

# ------------------------------------------------------------------------------
# 5. Summary comparison
# ------------------------------------------------------------------------------

message("\n=== Summary: CV-selected rho by A specification ===")
message(sprintf("  Uniform A (canonical):  rho=%.4f  (from results_sif_canonical)",
                tryCatch({
                  e <- new.env()
                  data("results_sif_canonical", package = "goebel2026", envir = e)
                  e$results_sif_canonical$rho_opt
                }, error = function(e) NA)))

for (nm in c("results_sif_gA_tau02", "results_sif_gA_tau033", "results_sif_gA_tau05")) {
  r <- tryCatch(get(nm), error = function(e) NULL)
  if (!is.null(r))
    message(sprintf("  g-weighted tau=%.3f:    rho=%.4f  phi=%.4f  beta_water=%.4f",
                    r$tags$tau, r$rho_opt, r$phi, r$beta_hat[2]))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_05_gA_cv.R complete. Objects saved to data/:")
message("  results_sif_gA_tau02")
message("  results_sif_gA_tau033")
message("  results_sif_gA_tau05")
