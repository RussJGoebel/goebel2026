# data-raw/run_08_likelihood_surface.R
#
# Computes ML likelihood and CV MSE profiles over a grid of phi values.
# Purpose: validate that random CV selects phi close to ML, while
# spatially blocked CV selects very different phi.
#
# Uses tune_cv and tune_ml for optima, and grid evaluation for curves.
# No RSR, no covariates. No R_inv for albedo, yes R_inv for SIF.
# rho=0.95 fixed throughout (avoids logdet singularity at rho=1).
#
# Outputs:
#   likelihood_surface_albedo
#   likelihood_surface_sif

library(fastblm)
library(goebel2026)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- TRUE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_sif    <- goebel2026::setup_sif

A       <- d_shared$A_flat
W_queen <- d_shared$W_queen
y_alb   <- d_albedo$y
y_sif   <- d_sif$y
R_inv   <- d_sif$R_inv

rho_val <- 0.95
S       <- Matrix::Diagonal(nrow(W_queen)) - rho_val * W_queen
Q_sp    <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
Q_fun   <- function(theta) list(Q = Q_sp)

log_phi_seq <- seq(log(1e-10), log(5000000), length.out = 150L)
phi_seq     <- exp(log_phi_seq)

set.seed(2026L)
folds_alb_random  <- sample(rep_len(1:10, length(y_alb)))
folds_sif_random  <- sample(rep_len(1:10, length(y_sif)))
folds_alb_blocked <- d_albedo$blocked_folds
folds_sif_blocked <- d_sif$blocked_folds

# ------------------------------------------------------------------------------
# 2. Helper: optima via tune_cv / tune_ml
# ------------------------------------------------------------------------------

get_cv_opt <- function(y, A, Q_fun, folds, R_inv = NULL, label = "") {
  message(sprintf("  tune_cv %s...", label))
  fastblm::tune_cv(
    y             = y,
    A             = A,
    Q_fun         = Q_fun,
    R_inv         = R_inv,
    theta_init    = numeric(0),
    k             = 10L,
    folds         = folds,
    solver        = "cholesky",
    log_phi_lower = log(0.0001),
    log_phi_upper = log(100000),
    parallel      = TRUE,
    verbose       = TRUE
  )
}

get_ml_opt <- function(y, A, Q_fun, R_inv = NULL, label = "") {
  message(sprintf("  tune_ml %s...", label))
  fastblm::tune_ml(
    y          = y,
    A          = A,
    Q_fun      = Q_fun,
    R_inv      = R_inv,
    theta_init = numeric(0),
    verbose    = TRUE
  )
}

# ------------------------------------------------------------------------------
# 3. Helper: CV MSE curve over phi grid
# ------------------------------------------------------------------------------

eval_cv_grid <- function(y, A, Q_sp, phi_seq, folds, R_inv = NULL) {
  n <- length(y)
  k <- max(folds)
  w <- if (is.null(R_inv)) rep(1.0, n) else Matrix::diag(R_inv)

  # Precompute fold indices
  fold_idx <- lapply(seq_len(k), function(fold)
    list(train = which(folds != fold), test = which(folds == fold)))

  # Capture everything needed explicitly
  y_     <- y
  A_     <- A
  Q_sp_  <- Q_sp
  w_     <- w

  unlist(future.apply::future_lapply(phi_seq, function(phi) {
    scores <- vapply(fold_idx, function(fi) {
      y_tr  <- y_[fi$train]
      A_tr  <- A_[fi$train, , drop = FALSE]
      A_te  <- A_[fi$test,  , drop = FALSE]
      y_te  <- y_[fi$test]
      Ri_tr <- Matrix::Diagonal(x = w_[fi$train])
      fit   <- tryCatch(
        fastblm::fit_fastblm(y_tr, A_tr, Q_sp_, phi = phi,
                             R_inv = Ri_tr, solver = "cholesky"),
        error = function(e) NULL)
      if (is.null(fit)) return(NA_real_)
      mean((y_te - as.numeric(A_te %*% fit$posterior_mean))^2,
           na.rm = TRUE)
    }, numeric(1L))
    mean(scores, na.rm = TRUE)
  },
  future.seed     = TRUE,
  future.packages = c("Matrix", "fastblm"),
  future.globals  = list(y_ = y_, A_ = A_, Q_sp_ = Q_sp_,
                         w_ = w_, fold_idx = fold_idx)))
}

# ------------------------------------------------------------------------------
# 4. Helper: ML likelihood curve over phi grid
# ------------------------------------------------------------------------------

eval_ml_grid <- function(y, A, Q_sp, phi_seq, R_inv = NULL) {
  n <- length(y)
  p <- ncol(A)

  CQ <- tryCatch(Matrix::Cholesky(Q_sp, LDL = FALSE, perm = TRUE),
                 error = function(e) NULL)
  logdet_Q <- if (!is.null(CQ))
    as.numeric(Matrix::determinant(CQ, logarithm = TRUE,
                                   sqrt = TRUE)$modulus) * 2
  else 0

  w    <- if (is.null(R_inv)) rep(1.0, n) else Matrix::diag(R_inv)
  wy   <- w * y
  yRy  <- sum(w * y^2)
  AtRy <- as.numeric(Matrix::crossprod(A, wy))
  AtRA <- Matrix::crossprod(A, A * w)

  vapply(phi_seq, function(phi) {
    K      <- Matrix::forceSymmetric(AtRA + (1/phi) * Q_sp)
    chol_K <- tryCatch(Matrix::Cholesky(K, LDL = FALSE, perm = TRUE),
                       error = function(e) NULL)
    if (is.null(chol_K)) return(NA_real_)
    mu      <- as.numeric(Matrix::solve(chol_K, AtRy))
    yHy     <- yRy - as.numeric(crossprod(AtRy, mu))
    if (yHy <= 0) return(NA_real_)
    sigma2e  <- yHy / n
    logdet_K <- as.numeric(
      Matrix::determinant(chol_K, logarithm = TRUE,
                          sqrt = TRUE)$modulus) * 2
    -n/2 * log(sigma2e) - 1/2 * logdet_K -
      p/2 * log(phi)    + 1/2 * logdet_Q
  }, numeric(1L))
}

# ==============================================================================
# 5. Albedo
# ==============================================================================

message("\n== Albedo ==")

tuned_alb_ml      <- get_ml_opt(y_alb, A, Q_fun, R_inv = NULL, "ML")
tuned_alb_cv_rand <- get_cv_opt(y_alb, A, Q_fun,
                                folds = folds_alb_random,
                                R_inv = NULL, "CV random")
tuned_alb_cv_blk  <- get_cv_opt(y_alb, A, Q_fun,
                                folds = folds_alb_blocked,
                                R_inv = NULL, "CV blocked")

message("  Computing albedo curves...")
ml_alb      <- eval_ml_grid(y_alb, A, Q_sp, phi_seq, R_inv = NULL)
cv_alb_rand <- eval_cv_grid(y_alb, A, Q_sp, phi_seq,
                            folds_alb_random,  R_inv = NULL)
cv_alb_blk  <- eval_cv_grid(y_alb, A, Q_sp, phi_seq,
                            folds_alb_blocked, R_inv = NULL)

message(sprintf("  Albedo: ML=%.2f  CV-rand=%.2f  CV-blk=%.4f",
                tuned_alb_ml$phi, tuned_alb_cv_rand$phi,
                tuned_alb_cv_blk$phi))

likelihood_surface_albedo <- list(
  log_phi_seq    = log_phi_seq,
  phi_seq        = phi_seq,
  rho            = rho_val,
  ml_ll          = ml_alb,
  cv_mse_random  = cv_alb_rand,
  cv_mse_blocked = cv_alb_blk,
  ml_phi_opt     = tuned_alb_ml$phi,
  cv_phi_random  = tuned_alb_cv_rand$phi,
  cv_phi_blocked = tuned_alb_cv_blk$phi,
  folds_random   = folds_alb_random,
  folds_blocked  = folds_alb_blocked,
  timestamp      = Sys.time()
)
usethis::use_data(likelihood_surface_albedo, overwrite = TRUE)

# ==============================================================================
# 6. SIF
# ==============================================================================

message("\n== SIF ==")

tuned_sif_ml      <- get_ml_opt(y_sif, A, Q_fun, R_inv = R_inv, "ML")
tuned_sif_cv_rand <- get_cv_opt(y_sif, A, Q_fun,
                                folds = folds_sif_random,
                                R_inv = R_inv, "CV random")
tuned_sif_cv_blk  <- get_cv_opt(y_sif, A, Q_fun,
                                folds = folds_sif_blocked,
                                R_inv = R_inv, "CV blocked")

message("  Computing SIF curves...")
ml_sif      <- eval_ml_grid(y_sif, A, Q_sp, phi_seq, R_inv = R_inv)
cv_sif_rand <- eval_cv_grid(y_sif, A, Q_sp, phi_seq,
                            folds_sif_random,  R_inv = R_inv)
cv_sif_blk  <- eval_cv_grid(y_sif, A, Q_sp, phi_seq,
                            folds_sif_blocked, R_inv = R_inv)

message(sprintf("  SIF: ML=%.4f  CV-rand=%.4f  CV-blk=%.4f",
                tuned_sif_ml$phi, tuned_sif_cv_rand$phi,
                tuned_sif_cv_blk$phi))

likelihood_surface_sif <- list(
  log_phi_seq    = log_phi_seq,
  phi_seq        = phi_seq,
  rho            = rho_val,
  ml_ll          = ml_sif,
  cv_mse_random  = cv_sif_rand,
  cv_mse_blocked = cv_sif_blk,
  ml_phi_opt     = tuned_sif_ml$phi,
  cv_phi_random  = tuned_sif_cv_rand$phi,
  cv_phi_blocked = tuned_sif_cv_blk$phi,
  folds_random   = folds_sif_random,
  folds_blocked  = folds_sif_blocked,
  timestamp      = Sys.time()
)
usethis::use_data(likelihood_surface_sif, overwrite = TRUE)

# ==============================================================================
# 7. Summary
# ==============================================================================

future::plan(future::sequential)

cat("\n=== Likelihood surface summary (rho=0.95, no covariates) ===\n")
cat(sprintf("\n%-10s  %-10s  %-12s  %-12s  %-8s  %-8s\n",
            "Response", "ML phi", "CV-rand phi", "CV-blk phi",
            "ML/CV-r", "ML/CV-b"))
cat(strrep("-", 72), "\n")
cat(sprintf("%-10s  %-10.2f  %-12.2f  %-12.4f  %-8.2f  %-8.2f\n",
            "Albedo",
            tuned_alb_ml$phi, tuned_alb_cv_rand$phi, tuned_alb_cv_blk$phi,
            tuned_alb_ml$phi / tuned_alb_cv_rand$phi,
            tuned_alb_ml$phi / tuned_alb_cv_blk$phi))
cat(sprintf("%-10s  %-10.4f  %-12.4f  %-12.4f  %-8.2f  %-8.2f\n",
            "SIF",
            tuned_sif_ml$phi, tuned_sif_cv_rand$phi, tuned_sif_cv_blk$phi,
            tuned_sif_ml$phi / tuned_sif_cv_rand$phi,
            tuned_sif_ml$phi / tuned_sif_cv_blk$phi))

message("\nrun_08 complete.")
message("  likelihood_surface_albedo saved to data/")
message("  likelihood_surface_sif     saved to data/")
