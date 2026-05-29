# =============================================================================
# Experiment: soft RSR penalty (alpha tuned) -- real SIF data
#
# Q_soft(v) = Q(v) + alpha * t(C_n) %*% (C_n %*% v)
# where C_n = t(X_fixed) %*% A / ||t(X_fixed) %*% A||_F
# alpha=0 unconstrained; alpha -> Inf recovers hard RSR.
# =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(ggplot2)

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

y_sif          <- goebel2026::soundings_augmented$SIF_757nm
water_sounding <- goebel2026::soundings_augmented$proportion_water
X_fixed        <- cbind(intercept = 1, water = water_sounding)
q              <- ncol(X_fixed)

water_grid_p   <- goebel2026::target_grid$proportion_water
water_grid_p[is.na(water_grid_p)] <- 0

# -----------------------------------------------------------------------------
# Standard W and base Q
# -----------------------------------------------------------------------------
W_std    <- goebel2026::make_W_matrix(goebel2026::target_grid)
IminusW  <- Matrix::Diagonal(p) - W_std
Q_sp     <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))

lambda_beta <- 0.01
Q_sp_aug    <- Matrix::bdiag(Q_sp, lambda_beta * Matrix::Diagonal(q))
A_aug       <- as(cbind(A, X_fixed), "dgCMatrix")

# RSR constraint matrix, normalised
C_rsr   <- as.matrix(t(X_fixed) %*% A)
C_rsr_n <- C_rsr / norm(C_rsr, "F")
CtC_diag <- as.numeric(Matrix::colSums(C_rsr_n^2))

# -----------------------------------------------------------------------------
# Q_fun
# -----------------------------------------------------------------------------
Q_fun <- function(theta) {
  alpha <- theta[["alpha"]]

  apply_Q_aug <- function(v) {
    v_sp   <- v[seq_len(p)]
    v_beta <- v[p + seq_len(q)]
    Cv     <- as.numeric(C_rsr_n %*% v_sp)
    CtCv   <- as.numeric(t(C_rsr_n) %*% Cv)
    c(as.numeric(Q_sp %*% v_sp) + alpha * CtCv,
      lambda_beta * v_beta)
  }

  d_aug   <- c(Matrix::diag(Q_sp) + alpha * CtC_diag, rep(lambda_beta, q))
  precond <- function(v) v / d_aug

  list(Q = apply_Q_aug, log_det_Q = NULL, precond = precond)
}

# -----------------------------------------------------------------------------
# CV
# -----------------------------------------------------------------------------
local({
  Q_sp_aug_captured <- Q_sp_aug
  precond_fun <<- function(phi, prior, A_train, y_train) {
    fit_base <- fastblm::fit_fastblm(y_train, A_train, Q_sp_aug_captured,
                                     phi = phi, solver = "cholesky")
    function(v) as.numeric(Matrix::solve(fit_base$chol_factor, v))
  }
})

cat("Tuning soft RSR alpha and phi ...\n")
t_cv_start <- proc.time()
tuned <- tune_cv(
  y           = y_sif,
  A           = A_aug,
  Q_fun       = Q_fun,
  theta_init  = c(alpha = 5.0),
  lower       = c(alpha = 0.0),
  upper       = c(alpha = 1e4),
  k           = 10L,
  solver      = "pcg",
  constraint  = NULL,
  precond_fun = precond_fun,
  verbose     = TRUE
)
t_cv <- proc.time() - t_cv_start
cat(sprintf("alpha = %.4f   phi = %.4f   CV-MSE = %.6f\n",
            tuned$theta[["alpha"]], tuned$phi, tuned$value))
cat(sprintf("CV time: %.1fs\n", t_cv["elapsed"]))

# -----------------------------------------------------------------------------
# Final fit
# -----------------------------------------------------------------------------
t_fit_start <- proc.time()
alpha_hat   <- tuned$theta[["alpha"]]
phi_hat     <- tuned$phi
prior_final <- Q_fun(c(alpha = alpha_hat))

fit <- fastblm::fit_fastblm(
  y           = y_sif,
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
# Posterior mean on grid
# -----------------------------------------------------------------------------
X_grid    <- cbind(intercept = 1, water = water_grid_p)
A_pred    <- as(cbind(Matrix::Diagonal(p), X_grid), "dgCMatrix")
pred_grid <- as.numeric(A_pred %*% fit$posterior_mean)

observed  <- which(Matrix::colSums(A) > 0)

# -----------------------------------------------------------------------------
# Map
# -----------------------------------------------------------------------------
plot_sf          <- target_proj
plot_sf$pred     <- ifelse(seq_len(p) %in% observed, pred_grid, NA_real_)
plot_sf$water    <- water_grid_p

ggplot(plot_sf) +
  geom_sf(aes(fill = pred), colour = NA) +
  scale_fill_gradient2(low  = "#2c7bb6", mid = "#ffffbf", high = "#d7191c",
                       midpoint = mean(plot_sf$pred, na.rm = TRUE),
                       na.value = "grey90", name = "SIF") +
  labs(title = sprintf("Soft RSR SIF  alpha=%.2f  phi=%.2f  beta_water=%.3f",
                       alpha_hat, phi_hat, beta_hat[2])) +
  theme_void(base_size = 11) +
  theme(plot.title     = element_text(face = "bold"),
        legend.key.height = unit(0.5, "cm"))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
t_total <- proc.time() - t_start
timings <- list(cv = t_cv, fit = t_fit, total = t_total)
cat(sprintf("Total time: %.1fs\n", t_total["elapsed"]))

results_softrsr_sif <- list(
  tuned     = tuned,
  fit       = fit,
  alpha_hat = alpha_hat,
  phi_hat   = phi_hat,
  beta_hat  = beta_hat,
  pred_grid = pred_grid,
  timings   = timings,
  p         = p,
  q         = q
)

usethis::use_data(results_softrsr_sif, overwrite = TRUE)
cat("Saved via use_data\n")
