# data-raw/run_sensitivity_neighbor.R
#
# Supplement S4: Neighbor structure sensitivity analysis.
# Compares queen (baseline), rook, and landcover-aware W adjacency.
# All runs use: water covariate + RSR + rho=1 (canonical model).
# Baseline result (queen) is results_water_rsr_rho1 from run_main_results.R.
#
# Outputs saved to data/ via usethis::use_data():
#   results_neighbor_rook     -- rook adjacency
#   results_neighbor_lc       -- landcover-aware W

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

fine_grid_buffered <- d_shared$fine_grid_buffered
A                  <- d_shared$A_flat
W_queen            <- d_shared$W_queen
W_rook             <- d_shared$W_rook
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
# 2. Landcover-aware W builders
#
# Two variants, both with CV-tuned alpha:
#
#   Variant A -- epsilon floor (min spillover):
#     w_ij = (1[same_type] + eps) / Z_i
#     Hard zeros on cross-boundary edges, but eps > 0 prevents isolation.
#     alpha here is log(1/eps), tuned so eps = exp(-alpha).
#
#   Variant B -- soft alpha scaling:
#     w_ij = exp(-alpha * |water_i - water_j|) / Z_i
#     Continuous downweighting proportional to water-fraction difference.
#     alpha = 0 recovers uniform queen; alpha -> Inf approaches hard zero.
#
# Both are row-normalised. alpha is tuned via CV as part of Q_fun(theta).
# ------------------------------------------------------------------------------

water_p  <- fine_grid_buffered$proportion_water
nb_queen <- spdep::poly2nb(fine_grid_buffered, queen = TRUE, snap = 330 * 0.01)
m        <- nrow(fine_grid_buffered)
is_water <- water_p > 0.5

# Variant A: epsilon floor -- eps = exp(-alpha), alpha >= 0
make_W_eps <- function(nb, is_water, alpha) {
  eps <- exp(-alpha)   # alpha=0 -> eps=1 (uniform); alpha large -> eps~0
  entries <- lapply(seq_along(nb), function(i) {
    js  <- nb[[i]]; js <- js[js > 0L]
    if (length(js) == 0L) return(list(j = integer(0), x = numeric(0)))
    raw <- ifelse(is_water[js] == is_water[i], 1, eps)
    list(j = js, x = raw / sum(raw))
  })
  i_idx <- rep(seq_along(nb), lengths(lapply(entries, `[[`, "j")))
  j_idx <- unlist(lapply(entries, `[[`, "j"))
  x_val <- unlist(lapply(entries, `[[`, "x"))
  Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                       dims = c(length(nb), length(nb)))
}

# Variant B: soft alpha scaling on |water_i - water_j|
make_W_alpha <- function(nb, water, alpha) {
  entries <- lapply(seq_along(nb), function(i) {
    js  <- nb[[i]]; js <- js[js > 0L]
    if (length(js) == 0L) return(list(j = integer(0), x = numeric(0)))
    raw <- exp(-alpha * abs(water[js] - water[i]))
    # add tiny floor to prevent complete isolation at very high alpha
    raw <- raw + 1e-6
    list(j = js, x = raw / sum(raw))
  })
  i_idx <- rep(seq_along(nb), lengths(lapply(entries, `[[`, "j")))
  j_idx <- unlist(lapply(entries, `[[`, "j"))
  x_val <- unlist(lapply(entries, `[[`, "x"))
  Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                       dims = c(length(nb), length(nb)))
}

# ------------------------------------------------------------------------------
# 3. Shared model components
# ------------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

make_Q_fun_fixed <- function(W) {
  S <- Matrix::Diagonal(nrow(W)) - W
  Q <- Matrix::forceSymmetric(Matrix::crossprod(S))
  Q <- Matrix::drop0(Q)
  function(theta) list(Q = Q, log_det_Q = NULL)
}

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

make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    C_zeros   <- matrix(0, nrow = q, ncol = q)
    cbind(C_spatial, C_zeros)
  }
}

A_aug          <- as(cbind(A, X_obs_water), "dgCMatrix")
rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)
C_spatial_full <- as.matrix(t(X_obs_water) %*% A)
C_zeros_full   <- matrix(0, nrow = q, ncol = q)
C_aug_full     <- cbind(C_spatial_full, C_zeros_full)

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

# Q_fun with alpha-tunable W variant
# make_W_fn: function(alpha) -> sparse W matrix
make_Q_fun_lc <- function(make_W_fn, q, lambda) {
  force(make_W_fn); force(q); force(lambda)
  function(theta) {
    alpha <- theta[["alpha"]]
    W     <- make_W_fn(alpha)
    S     <- Matrix::Diagonal(nrow(W)) - W
    Q_sp  <- Matrix::forceSymmetric(Matrix::crossprod(S))
    Q_sp  <- Matrix::drop0(Q_sp)
    Q_a   <- Matrix::forceSymmetric(
      Matrix::bdiag(Q_sp, lambda * Matrix::Diagonal(q))
    )
    list(Q = Q_a, log_det_Q = NULL)
  }
}

tune_and_fit_lc <- function(y, make_W_fn, k = 10L, seed = 2026L) {
  Q_fun_lc <- make_Q_fun_lc(make_W_fn, q, lambda_beta)

  tuned <- fastblm::tune_cv(
    y          = y,
    A          = A_aug,
    Q_fun      = Q_fun_lc,
    theta_init = c(alpha = 1),   # start at alpha=1, tune upward
    lower      = 0,
    upper      = 20,
    k          = k,
    constraint = rsr_constraint,
    seed       = seed,
    parallel   = TRUE,
    verbose    = TRUE
  )

  # Rebuild W at optimal alpha for the full-data constraint
  alpha_opt <- tuned$theta[["alpha"]]
  W_opt     <- make_W_fn(alpha_opt)

  fit <- fastblm::fit_fastblm(
    y      = y,
    A      = A_aug,
    Q      = tuned$Q,
    phi    = tuned$phi,
    solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  list(fit = fit, tuned = tuned, alpha_opt = alpha_opt, W_opt = W_opt)
}

compute_outputs <- function(result, run_name, tags) {
  fit   <- result$fit
  tuned <- result$tuned

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred)
  se_beta  <- fastblm::posterior_se(fit)[p + seq_len(q)]

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
    rho_opt               = 1,
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
# 4. Runs
# ------------------------------------------------------------------------------

# --- Rook adjacency -----------------------------------------------------------

if (!.should_skip("results_neighbor_rook")) {
  message("\n== S4: Rook adjacency ==")
  Q_fun_rook_sp  <- make_Q_fun_fixed(W_rook)
  Q_fun_rook_aug <- make_Q_fun_aug(Q_fun_rook_sp, q, lambda_beta)
  tuned_rook <- fastblm::tune_cv(
    y = y_alb, A = A_aug, Q_fun = Q_fun_rook_aug,
    theta_init = numeric(0), k = 10L, seed = 2026L,
    constraint = rsr_constraint, parallel = TRUE, verbose = TRUE
  )
  fit_rook <- fastblm::fit_fastblm(
    y = y_alb, A = A_aug, Q = tuned_rook$Q,
    phi = tuned_rook$phi, solver = "cholesky"
  )
  fit_rook <- fastblm::constrain(fit_rook, C_aug_full)
  message(sprintf("  tuning done: phi=%.4f  sigma2e=%.4g",
                  tuned_rook$phi, tuned_rook$sigma2e))
  out_rook <- compute_outputs(list(fit = fit_rook, tuned = tuned_rook),
                              "neighbor_rook",
                              tags = list(tuning = "cv", covariates = "water",
                                          constraint = "RSR", W = "rook", rho = 1))
  results_neighbor_rook <- out_rook
  usethis::use_data(results_neighbor_rook, overwrite = TRUE)
  message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  out_rook$rmse, out_rook$r2, out_rook$coverage_95_obs))
}

# --- Landcover-aware W, variant A: epsilon floor (alpha = log(1/eps)) --------

if (!.should_skip("results_neighbor_lc_eps")) {
  message("\n== S4: Landcover-aware W, epsilon floor (alpha CV-tuned) ==")
  make_W_eps_fn <- function(alpha) make_W_eps(nb_queen, is_water, alpha)
  r_lc_eps <- tune_and_fit_lc(y_alb, make_W_eps_fn)
  message(sprintf("  tuning done: alpha=%.4f  phi=%.4f  sigma2e=%.4g",
                  r_lc_eps$alpha_opt, r_lc_eps$tuned$phi,
                  r_lc_eps$tuned$sigma2e))
  out_lc_eps <- compute_outputs(r_lc_eps, "neighbor_lc_eps",
                                tags = list(tuning = "cv", covariates = "water",
                                            constraint = "RSR", W = "lc_eps",
                                            rho = 1, alpha = r_lc_eps$alpha_opt))
  results_neighbor_lc_eps <- out_lc_eps
  usethis::use_data(results_neighbor_lc_eps, overwrite = TRUE)
  message(sprintf("  alpha_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  r_lc_eps$alpha_opt, out_lc_eps$rmse,
                  out_lc_eps$r2, out_lc_eps$coverage_95_obs))
}

# --- Landcover-aware W, variant B: soft alpha scaling ------------------------

if (!.should_skip("results_neighbor_lc_alpha")) {
  message("\n== S4: Landcover-aware W, soft alpha scaling (alpha CV-tuned) ==")
  make_W_alpha_fn <- function(alpha) make_W_alpha(nb_queen, water_p, alpha)
  r_lc_alpha <- tune_and_fit_lc(y_alb, make_W_alpha_fn)
  message(sprintf("  tuning done: alpha=%.4f  phi=%.4f  sigma2e=%.4g",
                  r_lc_alpha$alpha_opt, r_lc_alpha$tuned$phi,
                  r_lc_alpha$tuned$sigma2e))
  out_lc_alpha <- compute_outputs(r_lc_alpha, "neighbor_lc_alpha",
                                  tags = list(tuning = "cv", covariates = "water",
                                              constraint = "RSR", W = "lc_alpha",
                                              rho = 1, alpha = r_lc_alpha$alpha_opt))
  results_neighbor_lc_alpha <- out_lc_alpha
  usethis::use_data(results_neighbor_lc_alpha, overwrite = TRUE)
  message(sprintf("  alpha_opt=%.4f  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f",
                  r_lc_alpha$alpha_opt, out_lc_alpha$rmse,
                  out_lc_alpha$r2, out_lc_alpha$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nS4 neighbor sensitivity results saved to data/.")
