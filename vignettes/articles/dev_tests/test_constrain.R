## =============================================================================
## test_constrain.R
##
## Simulation to verify that constrain() works correctly.
##
## We construct a small problem where:
##   1. We know the true beta and r
##   2. We can verify the constraint is satisfied
##   3. We can verify the predicted field is preserved after RSR
##   4. We can compare to a brute-force direct constrained solve
##
## Key test: the predicted field A_pred %*% mu_c should be close to
## A_pred %*% mu (unconstrained) -- RSR should be nearly a no-op
## on predictions when beta is correctly identified.
## =============================================================================

library(Matrix)
library(fastblm)

set.seed(42)

# =============================================================================
# 1. Small synthetic problem
# =============================================================================
# n = 50 soundings, p = 100 grid cells, q = 2 covariates (intercept + water)

n <- 50
p <- 100
q <- 2

# Random overlap matrix A (sparse, rows sum to ~1)
A_raw <- matrix(0, n, p)
for (i in seq_len(n)) {
  # Each sounding overlaps 3-5 adjacent cells
  start <- sample(seq_len(p - 4), 1)
  cells <- start:(start + sample(2:4, 1))
  weights <- runif(length(cells))
  A_raw[i, cells] <- weights / sum(weights)
}
A_sm <- as(A_raw, "dgCMatrix")

# Grid-level covariates
water_grid  <- runif(p, 0, 1)
X_grid_sim  <- cbind(1, water_grid)           # p x 2

# Sounding-level covariates (aggregated from grid)
water_snd   <- as.numeric(A_sm %*% water_grid)
X_fixed_sim <- cbind(1, water_snd)            # n x 2

# True parameters
beta_true <- c(0.2, -0.1)                     # intercept=0.2, water=-0.1
r_true    <- rnorm(p, 0, 0.05)               # spatial field, mean~0

# True field at grid level
y_grid_true <- as.numeric(X_grid_sim %*% beta_true) + r_true

# Observations
sigma2e_true <- 1e-4
y_obs <- as.numeric(A_sm %*% y_grid_true) + rnorm(n, 0, sqrt(sigma2e_true))

cat(sprintf("Setup: n=%d  p=%d  q=%d\n", n, p, q))
cat(sprintf("True beta: intercept=%.3f  water=%.3f\n",
            beta_true[1], beta_true[2]))
cat(sprintf("mean(r_true): %.6f  sd(r_true): %.6f\n",
            mean(r_true), sd(r_true)))

# =============================================================================
# 2. Augmented fit
# =============================================================================

A_aug_sim <- as(cbind(A_sm, X_fixed_sim), "dgCMatrix")

# SAR-like Q (simple second-order difference for 1D grid)
diags <- list(rep(2, p), rep(-1, p-1))
Q_sim <- bandSparse(p, k=c(0,1), diagonals=diags, symmetric=TRUE)
Q_sim <- as(forceSymmetric(Q_sim), "dgCMatrix")
# Fix boundary
Q_sim[1,1] <- 1; Q_sim[p,p] <- 1

phi_sim     <- 100
lambda_beta <- 0.01

Q_aug_sim <- forceSymmetric(bdiag(
  Q_sim,
  lambda_beta * Diagonal(q)
))

fit_sim <- fastblm::fit_fastblm(
  y      = y_obs,
  A      = A_aug_sim,
  Q      = Q_aug_sim,
  phi    = phi_sim,
  solver = "cholesky"
)

r_hat_sim    <- fit_sim$posterior_mean[seq_len(p)]
beta_hat_sim <- fit_sim$posterior_mean[p + seq_len(q)]

cat(sprintf("\nUnconstrained fit:\n"))
cat(sprintf("  beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_sim[1], beta_hat_sim[2]))
cat(sprintf("  mean(r_hat): %.6f\n", mean(r_hat_sim)))

# Predicted field before RSR
mu_pred_before_sim <- r_hat_sim + as.numeric(X_grid_sim %*% beta_hat_sim)
cat(sprintf("  cor(pred, truth): %.6f\n",
            cor(mu_pred_before_sim, y_grid_true)))

# =============================================================================
# 3. RSR constraint
# =============================================================================

C_spatial_sim <- as.matrix(t(X_fixed_sim) %*% A_sm)   # 2 x p
C_aug_sim     <- cbind(C_spatial_sim, matrix(0, q, q)) # 2 x (p+q)

cat(sprintf("\nBefore RSR: ||C r||_inf = %.4e\n",
            max(abs(C_spatial_sim %*% r_hat_sim))))

fit_sim_rsr <- fastblm::constrain(fit_sim, C_aug_sim)

r_hat_rsr_sim    <- fit_sim_rsr$posterior_mean[seq_len(p)]
beta_hat_rsr_sim <- fit_sim_rsr$posterior_mean[p + seq_len(q)]

cat(sprintf("\nRSR constrained fit:\n"))
cat(sprintf("  beta_hat: intercept=%.4f  water=%.4f\n",
            beta_hat_rsr_sim[1], beta_hat_rsr_sim[2]))
cat(sprintf("  mean(r_hat_rsr): %.6f\n", mean(r_hat_rsr_sim)))
cat(sprintf("  ||C r_rsr||_inf = %.4e\n",
            max(abs(C_spatial_sim %*% r_hat_rsr_sim))))

mu_pred_after_sim <- r_hat_rsr_sim + as.numeric(X_grid_sim %*% beta_hat_rsr_sim)
cat(sprintf("  cor(pred_rsr, truth): %.6f\n",
            cor(mu_pred_after_sim, y_grid_true)))
cat(sprintf("  cor(pred_before, pred_after): %.6f\n",
            cor(mu_pred_before_sim, mu_pred_after_sim)))
cat(sprintf("  ||pred_after - pred_before|| / ||pred_before||: %.4f\n",
            sqrt(sum((mu_pred_after_sim - mu_pred_before_sim)^2)) /
              sqrt(sum(mu_pred_before_sim^2))))

# =============================================================================
# 4. Brute-force constrained solve (ground truth)
# =============================================================================
# Solve the constrained system directly via null space of C_aug:
# Find basis B for null(C_aug), reparameterise gamma = B alpha,
# solve unconstrained problem in alpha space.

cat(sprintf("\n--- Brute force constrained solve ---\n"))

# Null space of C_aug (2 x (p+q)) via SVD
C_svd  <- svd(C_aug_sim, nu = 0, nv = ncol(C_aug_sim))
rank_C <- sum(C_svd$d > 1e-10 * max(C_svd$d))
B_null <- C_svd$v[, (rank_C+1):ncol(C_aug_sim), drop=FALSE]  # (p+q) x (p+q-2)

cat(sprintf("rank(C_aug): %d  null space dim: %d\n",
            rank_C, ncol(B_null)))

# Reparameterise: gamma = B_null %*% alpha
# Likelihood: y ~ N(A_aug B_null alpha, sigma2e I)
# Prior: alpha ~ N(0, sigma2e phi (B_null' Q_aug B_null)^{-1})
A_reduced  <- as.matrix(A_aug_sim %*% B_null)           # n x (p+q-2)
Q_reduced  <- t(B_null) %*% Q_aug_sim %*% B_null        # (p+q-2) x (p+q-2)

K_reduced  <- crossprod(A_reduced) + (1/phi_sim) * Q_reduced
mu_alpha   <- solve(K_reduced, t(A_reduced) %*% y_obs)
mu_bf      <- as.numeric(B_null %*% mu_alpha)           # back to original space

r_bf    <- mu_bf[seq_len(p)]
beta_bf <- mu_bf[p + seq_len(q)]

cat(sprintf("  beta_bf: intercept=%.4f  water=%.4f\n",
            beta_bf[1], beta_bf[2]))
cat(sprintf("  mean(r_bf): %.6f\n", mean(r_bf)))
cat(sprintf("  ||C r_bf||_inf = %.4e\n",
            max(abs(C_spatial_sim %*% r_bf))))

mu_pred_bf <- r_bf + as.numeric(X_grid_sim %*% beta_bf)
cat(sprintf("  cor(pred_bf, truth): %.6f\n",
            cor(mu_pred_bf, y_grid_true)))

# =============================================================================
# 5. Compare constrain() vs brute force
# =============================================================================

cat(sprintf("\n--- Comparison: constrain() vs brute force ---\n"))
cat(sprintf("  ||mu_rsr - mu_bf|| / ||mu_bf||: %.4e\n",
            sqrt(sum((fit_sim_rsr$posterior_mean - mu_bf)^2)) /
              sqrt(sum(mu_bf^2))))
cat(sprintf("  ||r_rsr - r_bf||:    %.4e\n",
            sqrt(sum((r_hat_rsr_sim - r_bf)^2))))
cat(sprintf("  |beta_rsr - beta_bf|: %.4e  %.4e\n",
            abs(beta_hat_rsr_sim[1] - beta_bf[1]),
            abs(beta_hat_rsr_sim[2] - beta_bf[2])))
cat(sprintf("  cor(pred_rsr, pred_bf): %.8f\n",
            cor(mu_pred_after_sim, mu_pred_bf)))

# =============================================================================
# 6. Repeat with intercept EXACTLY in col(A) -- your actual scenario
# =============================================================================

cat(sprintf("\n=== Test 2: intercept exactly in col(A) (your scenario) ===\n"))

# Make rows of A sum to exactly 1 (intercept exactly in col(A))
A_raw2 <- A_raw
A_raw2 <- A_raw2 / rowSums(A_raw2)   # normalise rows to sum to 1
A_sm2  <- as(A_raw2, "dgCMatrix")
cat(sprintf("rowSums range: [%.4f, %.4f]\n",
            min(rowSums(A_raw2)), max(rowSums(A_raw2))))

X_fixed_sim2 <- cbind(1, as.numeric(A_sm2 %*% water_grid))
A_aug_sim2   <- as(cbind(A_sm2, X_fixed_sim2), "dgCMatrix")

fit_sim2 <- fastblm::fit_fastblm(
  y=y_obs, A=A_aug_sim2, Q=Q_aug_sim, phi=phi_sim, solver="cholesky")

r2    <- fit_sim2$posterior_mean[seq_len(p)]
beta2 <- fit_sim2$posterior_mean[p + seq_len(q)]
cat(sprintf("Unconstrained: beta intercept=%.4f  water=%.4f  mean(r)=%.4f\n",
            beta2[1], beta2[2], mean(r2)))

C_spatial_sim2 <- as.matrix(t(X_fixed_sim2) %*% A_sm2)
C_aug_sim2     <- cbind(C_spatial_sim2, matrix(0, q, q))
cat(sprintf("||C r|| before RSR: %.4e\n",
            max(abs(C_spatial_sim2 %*% r2))))

fit_sim2_rsr <- fastblm::constrain(fit_sim2, C_aug_sim2)
r2_rsr    <- fit_sim2_rsr$posterior_mean[seq_len(p)]
beta2_rsr <- fit_sim2_rsr$posterior_mean[p + seq_len(q)]

mu_pred2_before <- r2    + as.numeric(X_grid_sim %*% beta2)
mu_pred2_after  <- r2_rsr + as.numeric(X_grid_sim %*% beta2_rsr)

cat(sprintf("RSR: beta intercept=%.4f  water=%.4f  mean(r)=%.4f\n",
            beta2_rsr[1], beta2_rsr[2], mean(r2_rsr)))
cat(sprintf("||C r_rsr||_inf: %.4e\n",
            max(abs(C_spatial_sim2 %*% r2_rsr))))
cat(sprintf("cor(pred_before, truth): %.6f\n",
            cor(mu_pred2_before, y_grid_true)))
cat(sprintf("cor(pred_after,  truth): %.6f\n",
            cor(mu_pred2_after,  y_grid_true)))
cat(sprintf("cor(pred_before, pred_after): %.6f\n",
            cor(mu_pred2_before, mu_pred2_after)))

# Brute force for this case too
C_svd2  <- svd(C_aug_sim2, nu = 0, nv = ncol(C_aug_sim2))
rank_C2  <- sum(C_svd2$d > 1e-10 * max(C_svd2$d))
B_null2  <- C_svd2$v[, (rank_C2+1):ncol(C_aug_sim2), drop=FALSE]
A_red2   <- as.matrix(A_aug_sim2 %*% B_null2)
Q_red2   <- t(B_null2) %*% Q_aug_sim %*% B_null2
K_red2   <- crossprod(A_red2) + (1/phi_sim) * Q_red2
mu_al2   <- solve(K_red2, t(A_red2) %*% y_obs)
mu_bf2   <- as.numeric(B_null2 %*% mu_al2)

r_bf2    <- mu_bf2[seq_len(p)]
beta_bf2 <- mu_bf2[p + seq_len(q)]
mu_pred_bf2 <- r_bf2 + as.numeric(X_grid_sim %*% beta_bf2)

cat(sprintf("\nBrute force RSR:\n"))
cat(sprintf("  beta: intercept=%.4f  water=%.4f\n",
            beta_bf2[1], beta_bf2[2]))
cat(sprintf("  cor(pred_bf, truth): %.6f\n",
            cor(mu_pred_bf2, y_grid_true)))
cat(sprintf("  ||constrain - brute_force|| / ||brute_force||: %.4e\n",
            sqrt(sum((fit_sim2_rsr$posterior_mean - mu_bf2)^2)) /
              sqrt(sum(mu_bf2^2))))
