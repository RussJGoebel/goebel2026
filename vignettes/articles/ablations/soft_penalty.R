# =============================================================================
# Experiment: soft RSR penalty (alpha tuned), no hard constraint -- albedo
#
# The spatial field is penalized for correlating with covariates via a
# rank-q addition to Q:
#   Q_soft(v) = Q(v) + alpha * t(C) %*% (C %*% v)
# where C = t(X_fixed) %*% A  (the RSR constraint matrix).
# alpha=0 is unconstrained; alpha -> Inf recovers hard RSR.
# Applied via PCG (matrix-free) so Q_soft is never formed explicitly.
# =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)

set.seed(42)
t_start <- proc.time()

# -----------------------------------------------------------------------------
# Data
# -----------------------------------------------------------------------------
soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)
p <- ncol(A)

water_sounding <- goebel2026::soundings_augmented$proportion_water
X_fixed        <- cbind(intercept = 1, water = water_sounding)
q              <- ncol(X_fixed)

water_grid_p   <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0

albedo_grid    <- goebel2026::target_grid$mean_albedo
noise_sd       <- 0.05 * sd(albedo_grid, na.rm = TRUE)
y_albedo       <- as.numeric(A %*% albedo_grid) + rnorm(nrow(A), sd = noise_sd)

cat(sprintf("Albedo noise SD: %.6f\n", noise_sd))

# -----------------------------------------------------------------------------
# Standard W and base Q
# -----------------------------------------------------------------------------
lambda_beta <- 0.01
W_std <- goebel2026::make_W_matrix(goebel2026::target_grid)

IminusW  <- Matrix::Diagonal(p) - W_std
Q_sp     <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))

Q_sp_aug <- Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))

# RSR constraint matrix C = t(X_fixed) %*% A,  q x p
# Normalised by Frobenius norm so alpha is on a natural scale
C_rsr    <- as.matrix(t(X_fixed) %*% A)
C_scale  <- norm(C_rsr, "F")
C_rsr_n  <- C_rsr / C_scale          # normalised

# -----------------------------------------------------------------------------
# Augmented design matrix
# -----------------------------------------------------------------------------

A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

# -----------------------------------------------------------------------------
# Q_fun: base Q_sp + alpha * t(C) C, applied matrix-free via closure
# -----------------------------------------------------------------------------
# Diagonal of t(C_n) %*% C_n -- precomputed once for preconditioner
CtC_diag <- as.numeric(Matrix::colSums(C_rsr_n^2))   # p-vector

Q_fun <- function(theta) {
  alpha <- theta[["alpha"]]

  # apply_Q(v) = Q_sp v + alpha * t(C_n) (C_n v)
  # Augmented Q acts on vectors of length p + q:
  #   spatial block gets the soft penalty, beta block gets lambda_beta ridge
  apply_Q_aug <- function(v) {
    v_sp   <- v[seq_len(p)]
    v_beta <- v[p + seq_len(q)]

    Cv     <- as.numeric(C_rsr_n %*% v_sp)           # q-vector
    CtCv   <- as.numeric(t(C_rsr_n) %*% Cv)          # p-vector

    Q_sp_v <- as.numeric(Q_sp %*% v_sp) + alpha * CtCv
    Q_b_v  <- lambda_beta * v_beta

    c(Q_sp_v, Q_b_v)
  }

  # Diagonal preconditioner: diag(Q_sp) + alpha * diag(t(C_n) C_n)
  # Inverted once per theta eval, reused across all PCG solves.
  d_sp    <- Matrix::diag(Q_sp) + alpha * CtC_diag
  d_aug   <- c(d_sp, rep(lambda_beta, q))
  precond <- function(v) v / d_aug

  list(Q = apply_Q_aug, log_det_Q = NULL, precond = precond)
}

# -----------------------------------------------------------------------------
# CV
# -----------------------------------------------------------------------------
cat("Tuning soft RSR alpha (log scale) and phi, PCG solver ...\n")
t_cv_start <- proc.time()
# Preconditioner: base Cholesky (no penalty) at current phi and training fold
precond_fun <- function(phi, prior, A_train, y_train) {
  fit_base <- fastblm::fit_fastblm(y_train, A_train, Q_sp_aug,
                                   phi = phi, solver = "cholesky")
  function(v) as.numeric(Matrix::solve(fit_base$chol_factor, v))
}

tuned <- tune_cv(
  y             = y_albedo,
  A             = A_aug,
  Q_fun         = Q_fun,
  theta_init    = c(alpha = 1.0),
  lower         = c(alpha = 0.0),
  upper         = c(alpha = 1e4),
  k             = 10L,
  solver        = "pcg",
  constraint    = NULL,
  precond_fun   = precond_fun,
  verbose       = TRUE
)
t_cv <- proc.time() - t_cv_start
cat(sprintf("alpha = %.4f   phi = %.4f   CV-MSE = %.6f\n",
            tuned$theta[["alpha"]], tuned$phi, tuned$value))
cat(sprintf("CV time: %.1fs\n", t_cv["elapsed"]))

# -----------------------------------------------------------------------------
# Final fit
# -----------------------------------------------------------------------------
t_fit_start <- proc.time()
alpha_hat <- tuned$theta[["alpha"]]
phi_hat   <- tuned$phi

prior_final <- Q_fun(c(alpha = alpha_hat))
fit <- fastblm::fit_fastblm(
  y           = y_albedo,
  A           = A_aug,
  Q           = prior_final$Q,
  phi         = phi_hat,
  solver      = "pcg",
  pcg_precond = prior_final$precond
)
t_fit <- proc.time() - t_fit_start

beta_hat <- fit$posterior_mean[p + seq_len(q)]
cat(sprintf("beta: intercept = %.4f   water = %.4f\n",
            beta_hat[1], beta_hat[2]))
cat(sprintf("Fit time: %.1fs\n", t_fit["elapsed"]))

# -----------------------------------------------------------------------------
# Posterior mean on grid (PCG path: no posterior_se without Hutchinson)
# -----------------------------------------------------------------------------
X_grid    <- cbind(intercept = 1, water = water_grid_p)
A_pred    <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred_grid <- as.numeric(A_pred %*% fit$posterior_mean)

# -----------------------------------------------------------------------------
# Summary statistics on observed cells
# -----------------------------------------------------------------------------
observed  <- which(Matrix::colSums(A) > 0)
resid     <- pred_grid[observed] - albedo_grid[observed]
r2        <- 1 - sum(resid^2) / sum((albedo_grid[observed] - mean(albedo_grid[observed]))^2)
rmse      <- sqrt(mean(resid^2))

cat(sprintf("\n--- Summary (observed cells, n=%d) ---\n", length(observed)))
cat(sprintf("RMSE : %.6f\n", rmse))
cat(sprintf("R2   : %.4f\n", r2))

stats <- list(rmse = rmse, r2 = r2, n_observed = length(observed))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
t_total <- proc.time() - t_start
timings <- list(cv = t_cv, fit = t_fit, total = t_total)
cat(sprintf("Total time: %.1fs\n", t_total["elapsed"]))

results_softrsr_norsr_albedo <- list(
  tuned      = tuned,
  fit        = fit,
  alpha_hat  = alpha_hat,
  phi_hat    = phi_hat,
  beta_hat   = beta_hat,
  pred_grid  = pred_grid,
  stats      = stats,
  noise_sd   = noise_sd,
  y_albedo   = y_albedo,
  timings    = timings,
  p          = p,
  q          = q
)

usethis::use_data(results_softrsr_norsr_albedo, overwrite = TRUE)
cat("Saved via use_data\n")
