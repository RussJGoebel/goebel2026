## =============================================================================
## profile_pcg.R
##
## Profiles a single PCG solve in isolation to identify the bottleneck.
## Self-contained -- just run this script directly.
## =============================================================================

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)

set.seed(42)

# =============================================================================
# 1. Data setup
# =============================================================================

noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$synthetic_albedo_A_upscaled +
  rnorm(length(goebel2026::soundings_augmented$synthetic_albedo_A_upscaled),
        0, noise_sd)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

A <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)

W <- goebel2026::make_W_matrix(goebel2026::target_grid)

Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  list(Q = Matrix::drop0(Q))
}

# =============================================================================
# 2. Build operators at rho=0.95, phi=2  (near the slow region)
# =============================================================================

rho <- 0.95
phi <- 2.0
p   <- ncol(A)
n   <- nrow(A)

cat(sprintf("n=%d  p=%d  nnz(A)=%d\n", n, p, Matrix::nnzero(A)))

Q           <- Q_fun(c(rho = rho))$Q
apply_Q     <- function(v) as.numeric(Q %*% v)
apply_AtA   <- function(v) as.numeric(Matrix::crossprod(A, A %*% v))
apply_K     <- function(v) apply_AtA(v) + (1/phi) * apply_Q(v)
AtRinvy     <- as.numeric(Matrix::crossprod(A, y_alb))

# =============================================================================
# 3. Build Q preconditioner
# =============================================================================

eps_prec  <- 1e-4 * mean(Matrix::diag(Q))
Q_reg     <- Matrix::forceSymmetric(Q + eps_prec * Matrix::Diagonal(p))
chol_prec <- Matrix::Cholesky(Q_reg, LDL = FALSE, perm = TRUE)
precond_Q <- function(v) as.numeric(Matrix::solve(chol_prec, v))

# =============================================================================
# 4. Warmup run (avoid profiling JIT / first-call overhead)
# =============================================================================

cat("Warming up PCG...\n")
res_warm <- fastblm:::pcg(apply_K, AtRinvy, precond = precond_Q,
                          tol = 1e-6, maxit = 4L * p)
cat(sprintf("Warmup: %d iterations, converged=%s\n",
            res_warm$iter, res_warm$converged))

# =============================================================================
# 5. Timed run without profiler -- baseline
# =============================================================================

cat("\nTiming single PCG solve...\n")
t0 <- proc.time()
res <- fastblm:::pcg(apply_K, AtRinvy, precond = precond_Q,
                     tol = 1e-6, maxit = 4L * p)
t1 <- proc.time()

elapsed <- (t1 - t0)[["elapsed"]]
cat(sprintf("PCG: %d iterations, converged=%s, time=%.2fs\n",
            res$iter, res$converged, elapsed))
cat(sprintf("Time per iteration: %.2fms\n", 1000 * elapsed / res$iter))

# =============================================================================
# 6. Profile run
# =============================================================================

cat("\nProfiling PCG solve...\n")
prof_file <- tempfile(fileext = ".out")

Rprof(prof_file, interval = 0.005)
res_prof <- fastblm:::pcg(apply_K, AtRinvy, precond = precond_Q,
                          tol = 1e-6, maxit = 4L * p)
Rprof(NULL)

cat(sprintf("Profile written to: %s\n", prof_file))

# =============================================================================
# 7. Summary
# =============================================================================

prof <- summaryRprof(prof_file)

cat("\n=== Top 20 by SELF time (where time is actually spent) ===\n")
print(head(prof$by.self, 20))

cat("\n=== Top 20 by TOTAL time (call stack context) ===\n")
print(head(prof$by.total, 20))

# =============================================================================
# 8. Also time the individual components separately
# =============================================================================

cat("\n=== Component timings (100 calls each) ===\n")

v <- rnorm(p)

t_atA <- system.time(for (i in seq_len(100)) apply_AtA(v))[["elapsed"]]
cat(sprintf("apply_AtA (A'Av):       %.2fms per call\n", 1000 * t_atA / 100))

t_Q <- system.time(for (i in seq_len(100)) apply_Q(v))[["elapsed"]]
cat(sprintf("apply_Q (Qv):           %.2fms per call\n", 1000 * t_Q / 100))

t_K <- system.time(for (i in seq_len(100)) apply_K(v))[["elapsed"]]
cat(sprintf("apply_K (Kv):           %.2fms per call\n", 1000 * t_K / 100))

t_precond <- system.time(for (i in seq_len(100)) precond_Q(v))[["elapsed"]]
cat(sprintf("precond_Q (Q^{-1}v):    %.2fms per call\n", 1000 * t_precond / 100))

t_add <- system.time(for (i in seq_len(100)) { x <- v + 0.5 * v })[["elapsed"]]
cat(sprintf("vector add (x+alpha*d): %.2fms per call\n", 1000 * t_add / 100))

t_dot <- system.time(for (i in seq_len(100)) crossprod(v, v))[["elapsed"]]
cat(sprintf("crossprod (r'z):        %.2fms per call\n", 1000 * t_dot / 100))

cat(sprintf("\nPredicted time per PCG iter: %.2fms\n",
            1000 * (t_K + t_precond + 4 * t_add + 2 * t_dot) / 100))
cat(sprintf("Actual time per PCG iter:    %.2fms\n",
            1000 * elapsed / res$iter))

# flamegraph if profvis available
if (requireNamespace("profvis", quietly = TRUE)) {
  cat("\nOpening profvis flame graph...\n")
  print(profvis::profvis(prof_input = prof_file))
} else {
  cat("\nTip: install.packages('profvis') for a flame graph\n")
}
