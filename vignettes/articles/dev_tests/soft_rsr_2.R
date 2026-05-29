# test_pcg_tol_fix.R
#
# Benchmarks the effect of the relative PCG tolerance fix on:
#   (1) PCG iteration counts at a fixed (rho, phi) across multiple folds
#   (2) Wall time for a full .eval_cv call (k=10, sequential)
#
# Run this BEFORE and AFTER updating solvers.R to confirm the fix.
# The script prints a summary table and a timing comparison.
#
# Expected result after fix:
#   - PCG iterations collapse from ~maxit to O(10-50)
#   - .eval_cv wall time drops from 5+ min to ~2 min

library(fastblm)
library(goebel2026)
library(Matrix)

# ------------------------------------------------------------------------------
# 1. Load setup objects (same as run_sensitivity_softrsr.R)
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo

A           <- d_shared$A_flat
W_queen     <- d_shared$W_queen
X_obs_water <- d_shared$X_obs_water
y_alb       <- d_albedo$y

p <- ncol(A)
q <- ncol(X_obs_water)
lambda_beta <- 0.01

A_aug   <- as(cbind(A, X_obs_water), "dgCMatrix")
C_raw   <- as.matrix(t(X_obs_water) %*% A)
C_scale <- norm(C_raw, "F")
C_rsr_n <- C_raw / C_scale

# ------------------------------------------------------------------------------
# 2. Build Q_fun and precond_fun (identical to sensitivity script)
# ------------------------------------------------------------------------------

Q_fun_soft <- function(theta) {
  rho        <- min(theta[["rho"]], 0.98)
  alpha      <- exp(theta[["log_alpha"]])
  alpha_safe <- max(alpha, 1e-10)
  S_rho      <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
  Q_rho      <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_rho)))
  .Q_rho      <- Q_rho
  .alpha_safe <- alpha_safe
  .C          <- C_rsr_n
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

precond_fun <- function(phi, prior, A_train, y_train) {
  AtA_f    <- Matrix::forceSymmetric(Matrix::crossprod(A_train))
  Q_approx <- Matrix::forceSymmetric(
    Matrix::bdiag(prior$Q_rho * phi, lambda_beta * phi * Matrix::Diagonal(q))
  ) + AtA_f / 3e-6
  .chol <- Matrix::Cholesky(Q_approx, LDL = FALSE, perm = TRUE)
  force(.chol)
  function(v) as.numeric(Matrix::solve(.chol, v))
}

# ------------------------------------------------------------------------------
# 3. Fixed evaluation point
# ------------------------------------------------------------------------------

RHO   <- 0.9
PHI   <- 82
ALPHA <- 1.0   # log_alpha = 0

theta_test <- c(rho = RHO, log_alpha = log(ALPHA))
prior_test <- Q_fun_soft(theta_test)

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_alb), 10L)
score_fn         <- fastblm:::.make_score_fn("mse")

# ------------------------------------------------------------------------------
# 4. Per-fold PCG iteration count diagnostic
#
#    We instrument a single PCG call per fold to count iterations.
#    apply_K is built manually to mirror what fit_fastblm does internally.
# ------------------------------------------------------------------------------

message("\n=== PCG iteration count per fold at rho=", RHO,
        "  phi=", PHI, "  alpha=", ALPHA, " ===\n")

make_apply_K <- fastblm:::make_apply_K

iter_results <- lapply(seq_len(10L), function(fold) {
  train_idx <- which(fold_assignments != fold)
  A_tr      <- A_aug[train_idx, , drop = FALSE]
  y_tr      <- y_alb[train_idx]

  apply_A    <- function(v) as.numeric(A_tr %*% v)
  apply_At   <- function(v) as.numeric(Matrix::crossprod(A_tr, v))
  apply_Rinv <- function(v) v
  apply_K    <- make_apply_K(apply_A, apply_At, prior_test$Q, apply_Rinv, PHI)

  Rinvy   <- y_tr
  AtRinvy <- apply_At(Rinvy)
  b_norm  <- sqrt(sum(AtRinvy^2))

  precond <- precond_fun(PHI, prior_test, A_tr, y_tr)

  result <- fastblm:::pcg(apply_K, AtRinvy,
                          tol    = 1e-6,
                          maxit  = 2L * (p + q),
                          precond = precond)

  data.frame(
    fold       = fold,
    n_train    = length(train_idx),
    b_norm     = round(b_norm, 2),
    converged  = result$converged,
    iterations = result$iter,
    maxit      = 2L * (p + q)
  )
})

iter_df <- do.call(rbind, iter_results)
print(iter_df, row.names = FALSE)

message(sprintf(
  "\nMean iterations: %.1f  |  Max: %d  |  Converged: %d/%d",
  mean(iter_df$iterations),
  max(iter_df$iterations),
  sum(iter_df$converged),
  nrow(iter_df)
))

message(sprintf(
  "Mean b_norm: %.1f  -- absolute tol 1e-6 would require residual < %.2e,",
  mean(iter_df$b_norm),
  1e-6
))
message(sprintf(
  "  relative tol 1e-6 requires residual < %.2e  (%.0fx easier)",
  mean(iter_df$b_norm) * 1e-6,
  mean(iter_df$b_norm)
))

# ------------------------------------------------------------------------------
# 5. Wall-time comparison: .eval_cv sequential, k=10
# ------------------------------------------------------------------------------

message("\n=== .eval_cv wall-time (sequential, k=10) ===\n")

t_eval <- system.time({
  cv_val <- fastblm:::.eval_cv(
    y_alb, A_aug, prior_test, PHI, fold_assignments,
    score_fn,
    "pcg", 1e-6, NULL, NULL,
    fold_C_list       = NULL,
    fold_precond_list = NULL,
    precond_fun       = precond_fun,
    parallel          = FALSE
  )
})

message(sprintf("  Total:    %.2fs", t_eval["elapsed"]))
message(sprintf("  Per fold: %.2fs", t_eval["elapsed"] / 10))
message(sprintf("  cv_mse:   %.4e", cv_val))
message(sprintf(
  "\n  Target: ~2 min total (~12s/fold). Previous (broken): 5+ min (~30s/fold)."
))

# ------------------------------------------------------------------------------
# 6. Spot-check: does relative tol give same answer as a tighter absolute tol?
# ------------------------------------------------------------------------------

message("\n=== Solution accuracy spot-check (fold 1) ===\n")

fold       <- 1L
train_idx  <- which(fold_assignments != fold)
A_tr       <- A_aug[train_idx, , drop = FALSE]
y_tr       <- y_alb[train_idx]
apply_A    <- function(v) as.numeric(A_tr %*% v)
apply_At   <- function(v) as.numeric(Matrix::crossprod(A_tr, v))
apply_Rinv <- function(v) v
apply_K    <- make_apply_K(apply_A, apply_At, prior_test$Q, apply_Rinv, PHI)
AtRinvy    <- apply_At(y_tr)
precond    <- precond_fun(PHI, prior_test, A_tr, y_tr)

sol_rel <- fastblm:::pcg(apply_K, AtRinvy, tol = 1e-6,
                         maxit = 2L * (p + q), precond = precond)
sol_tight <- fastblm:::pcg(apply_K, AtRinvy, tol = 1e-10,
                           maxit = 2L * (p + q), precond = precond)

resid_rel   <- sqrt(sum((apply_K(sol_rel$x)   - AtRinvy)^2))
resid_tight <- sqrt(sum((apply_K(sol_tight$x) - AtRinvy)^2))
sol_diff    <- sqrt(sum((sol_rel$x - sol_tight$x)^2)) /
  sqrt(sum(sol_tight$x^2))

message(sprintf("  Relative tol (1e-6):  iters=%d  |Kx-b|=%.2e",
                sol_rel$iter, resid_rel))
message(sprintf("  Tight tol   (1e-10):  iters=%d  |Kx-b|=%.2e",
                sol_tight$iter, resid_tight))
message(sprintf("  Relative solution difference: %.2e  (should be << 1e-3)",
                sol_diff))

message("\nDone.")
