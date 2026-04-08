## Fit downscaling model on OCO-2 SIF data
##
## Steps:
## 1. Compute aggregation matrix A
## 2. Precompute W (queen adjacency) for efficiency inside Q_fun
## 3. Use tune_cv to jointly select rho and phi via 10-fold CV,
##    with folds evaluated in parallel
## 4. Fit final model at optimum
## 5. Extract posterior mean and SE
###############################################################################

library(goebel2026)
library(spatintegrate)
library(fastblm)

### 1) Load data ###############################################################

target_grid   <- goebel2026::target_grid
soundings     <- goebel2026::soundings
soundings_aug <- goebel2026::soundings_augmented

### 2) Compute A matrix ########################################################

message("Computing A matrix...")
A <- spatintegrate::compute_overlap_fractions(
  spatintegrate::ensure_projected(soundings),
  spatintegrate::ensure_projected(target_grid)
)

### 3) Set up response #########################################################

y <- soundings$SIF_757nm

### 4) Precompute W (so poly2nb is only called once) ##########################

message("Computing W matrix...")
W <- goebel2026::make_W_matrix(target_grid)

Q_fun <- function(theta) {
  rho     <- theta[["rho"]]
  IminusW <- Matrix::Diagonal(nrow(W)) - rho * W
  Q       <- Matrix::forceSymmetric(Matrix::crossprod(IminusW))
  Q       <- Matrix::drop0(Q)
  list(Q = Q)
}

### 5) Tune rho and phi via 10-fold CV with parallel fold evaluation ###########

message("Starting CV tuning...")
future::plan(future::multisession())

tuned <- fastblm::tune_cv(
  y          = y,
  A          = A,
  Q_fun      = Q_fun,
  theta_init = c(rho = 0.9),
  lower      = 0.01,
  upper      = 1.0,
  k          = 10L,
  solver     = "cholesky",
  parallel   = TRUE,
  verbose    = TRUE
)

future::plan(future::sequential())

message(sprintf("Optimal rho: %.4f", tuned$theta[["rho"]]))
message(sprintf("Optimal phi: %.4f", tuned$phi))
message(sprintf("sigma2e:     %.4f", tuned$sigma2e))

### 6) Fit final model at optimum #############################################

message("Fitting final model...")
fit <- fastblm::fit_fastblm(
  y      = y,
  A      = A,
  Q      = tuned$Q,
  phi    = tuned$phi,
  solver = "cholesky"
)

### 7) Extract posterior mean and SE ##########################################

posterior_mean <- fit$posterior_mean
posterior_se   <- fastblm::posterior_se(fit, A_new = diag(dim(A)[2]))

# attach to target_grid
target_grid$SIF_posterior_mean <- posterior_mean
target_grid$SIF_posterior_se   <- posterior_se

### 8) Quick diagnostics #######################################################

cat("\nPosterior mean summary:\n")
print(summary(posterior_mean))

cat("\nPosterior SE summary:\n")
print(summary(posterior_se))

