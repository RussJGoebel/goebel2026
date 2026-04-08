## =============================================================================
## test_pcg_tol.R
##
## Tests PCG convergence at different tolerances, with and without
## preconditioner, to find the best tradeoff for REML tuning.
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)

set.seed(42)

# =============================================================================
# 1. Setup (same as profile_pcg.R)
# =============================================================================

noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)
A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj,
                                           spatintegrate::ensure_projected(goebel2026::target_grid)),
  "dgCMatrix"
)
W <- goebel2026::make_W_matrix(goebel2026::target_grid)
Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  list(Q = Matrix::drop0(Q))
}

rho <- 0.95
phi <- 2.0
p   <- ncol(A)
n   <- nrow(A)

Q          <- Q_fun(c(rho = rho))$Q
apply_Q    <- function(v) as.numeric(Q %*% v)
apply_AtA  <- function(v) as.numeric(Matrix::crossprod(A, A %*% v))
apply_K    <- function(v) apply_AtA(v) + (1/phi) * apply_Q(v)
AtRinvy    <- as.numeric(Matrix::crossprod(A, y_alb))

# Build preconditioner
eps_prec  <- 1e-4 * mean(Matrix::diag(Q))
Q_reg     <- Matrix::forceSymmetric(Q + eps_prec * Matrix::Diagonal(p))
chol_prec <- Matrix::Cholesky(Q_reg, LDL = FALSE, perm = TRUE)
precond_Q <- function(v) as.numeric(Matrix::solve(chol_prec, v))

# Get reference solution at tight tolerance
cat("Computing reference solution (tol=1e-8)...\n")
x_ref <- fastblm:::pcg(apply_K, AtRinvy, tol = 1e-8, maxit = 4L * p)$x
cat(sprintf("Reference: ||x_ref|| = %.6f\n", sqrt(sum(x_ref^2))))

# =============================================================================
# 2. Test tolerances without preconditioner
# =============================================================================

cat("\n=== No preconditioner ===\n")
cat(sprintf("%-10s  %-8s  %-8s  %-10s  %-10s\n",
            "tol", "iters", "time(s)", "||err||/||x||", "converged"))

for (tol in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
  t0  <- proc.time()[["elapsed"]]
  res <- fastblm:::pcg(apply_K, AtRinvy, tol = tol, maxit = 4L * p)
  t1  <- proc.time()[["elapsed"]]
  rel_err <- sqrt(sum((res$x - x_ref)^2)) / sqrt(sum(x_ref^2))
  cat(sprintf("%-10.0e  %-8d  %-8.2f  %-10.2e  %-10s\n",
              tol, res$iter, t1 - t0, rel_err,
              ifelse(res$converged, "yes", "NO")))
}

# =============================================================================
# 3. Test tolerances with preconditioner
# =============================================================================

cat("\n=== With Q preconditioner ===\n")
cat(sprintf("%-10s  %-8s  %-8s  %-10s  %-10s\n",
            "tol", "iters", "time(s)", "||err||/||x||", "converged"))

for (tol in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
  t0  <- proc.time()[["elapsed"]]
  res <- fastblm:::pcg(apply_K, AtRinvy, precond = precond_Q,
                       tol = tol, maxit = 4L * p)
  t1  <- proc.time()[["elapsed"]]
  rel_err <- sqrt(sum((res$x - x_ref)^2)) / sqrt(sum(x_ref^2))
  cat(sprintf("%-10.0e  %-8d  %-8.2f  %-10.2e  %-10s\n",
              tol, res$iter, t1 - t0, rel_err,
              ifelse(res$converged, "yes", "NO")))
}

# =============================================================================
# 4. How much does solution accuracy affect the ll estimate?
# =============================================================================

cat("\n=== Effect on REML log-likelihood estimate ===\n")

n_eff  <- n
yRinvy <- as.numeric(crossprod(y_alb, y_alb))

.ll_from_x <- function(x, phi, p, n_eff, AtRinvy, yRinvy) {
  yHinvy <- yRinvy - as.numeric(crossprod(AtRinvy, x))
  if (yHinvy <= 0) return(NA)
  sigma2e <- yHinvy / n_eff
  # omit logdet terms -- just checking sigma2e sensitivity
  -n_eff/2 * log(sigma2e)
}

ll_ref <- .ll_from_x(x_ref, phi, p, n_eff, AtRinvy, yRinvy)
cat(sprintf("Reference ll component: %.4f\n\n", ll_ref))

cat(sprintf("%-10s  %-12s  %-10s\n", "tol", "ll component", "ll error"))
for (tol in c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6)) {
  res <- fastblm:::pcg(apply_K, AtRinvy, tol = tol, maxit = 4L * p)
  ll  <- .ll_from_x(res$x, phi, p, n_eff, AtRinvy, yRinvy)
  cat(sprintf("%-10.0e  %-12.4f  %-10.4f\n", tol, ll, ll - ll_ref))
}
