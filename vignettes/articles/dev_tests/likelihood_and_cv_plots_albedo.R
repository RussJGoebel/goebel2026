# plot_likelihood_surface_albedo.R
#
# Plots 1D ML likelihood and CV MSE curves as a function of log-phi,
# at fixed rho values. No covariates, no RSR, no R_inv.
# Uses semi-synthetic albedo data where ground truth is known.

library(fastblm)
library(goebel2026)
library(Matrix)
library(ggplot2)
library(patchwork)
library(future)
library(future.apply)

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo

A       <- d_shared$A_flat
W_queen <- d_shared$W_queen
y_alb   <- d_albedo$y

p <- ncol(A)

set.seed(2026L)
fold_assignments <- fastblm:::.make_folds(length(y_alb), 10L)
score_fn         <- fastblm:::.make_score_fn("mse")

# ------------------------------------------------------------------------------
# 2. Grid
# ------------------------------------------------------------------------------

rho_values  <- c(0.7, 0.9, 0.95, 0.99)
log_phi_seq <- seq(log(0.005), log(10000), length.out = 35)
phi_seq     <- exp(log_phi_seq)
rho_cv      <- 0.99

# ------------------------------------------------------------------------------
# 3. ML likelihood grid
# ------------------------------------------------------------------------------

message("Evaluating ML likelihood grid (albedo)...")

if (file.exists("likelihood_surface_albedo_ml_grid.rds")) {
  message("  loading from cache...")
  ml_grid <- readRDS("likelihood_surface_albedo_ml_grid.rds")
} else {
  results <- lapply(rho_values, function(rho) {
    message(sprintf("\n  ML rho=%.2f", rho))

    S    <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
    Q_sp <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))

    CQ <- tryCatch(Matrix::Cholesky(Q_sp, LDL = FALSE, perm = TRUE),
                   error = function(e) NULL)
    logdet_Q <- if (!is.null(CQ))
      as.numeric(Matrix::determinant(CQ, logarithm = TRUE, sqrt = TRUE)$modulus) * 2
    else 0

    AtA     <- Matrix::crossprod(A)
    Aty     <- as.numeric(Matrix::crossprod(A, y_alb))
    yty     <- as.numeric(crossprod(y_alb))
    n       <- length(y_alb)

    ll_vals <- vapply(phi_seq, function(phi) {
      K      <- Matrix::forceSymmetric(AtA + (1/phi) * Q_sp)
      chol_K <- tryCatch(Matrix::Cholesky(K, LDL = FALSE, perm = TRUE),
                         error = function(e) NULL)
      if (is.null(chol_K)) return(NA_real_)
      mu     <- as.numeric(Matrix::solve(chol_K, Aty))
      yHy    <- yty - as.numeric(crossprod(Aty, mu))
      if (yHy <= 0) return(NA_real_)
      sigma2e  <- yHy / n
      logdet_K <- as.numeric(
        Matrix::determinant(chol_K, logarithm = TRUE, sqrt = TRUE)$modulus) * 2
      -n/2 * log(sigma2e) - 1/2 * logdet_K - p/2 * log(phi) + 1/2 * logdet_Q
    }, numeric(1L))

    message(sprintf("    done. ll range [%.1f, %.1f]",
                    min(ll_vals, na.rm=TRUE), max(ll_vals, na.rm=TRUE)))
    data.frame(rho = rho, log_phi = log_phi_seq, phi = phi_seq, ll = ll_vals)
  })

  ml_grid <- do.call(rbind, results)
  saveRDS(ml_grid, "likelihood_surface_albedo_ml_grid.rds")
  message("  ml_grid saved.")
}

# ------------------------------------------------------------------------------
# 4. CV grid at rho=0.99
# ------------------------------------------------------------------------------

message(sprintf("\nEvaluating CV MSE at rho=%.2f (albedo)...", rho_cv))

if (file.exists("likelihood_surface_albedo_cv_grid.rds")) {
  message("  loading from cache...")
  cv_grid <- readRDS("likelihood_surface_albedo_cv_grid.rds")
} else {
  S_cv    <- Matrix::Diagonal(nrow(W_queen)) - rho_cv * W_queen
  Q_sp_cv <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S_cv)))

  cv_vals <- vapply(phi_seq, function(phi) {
    prior <- list(Q = Q_sp_cv, log_det_Q = NULL)
    tryCatch(
      fastblm:::.eval_cv(
        y_alb, A, function(theta) prior, numeric(0), prior,
        phi, fold_assignments, score_fn,
        "cholesky", 1e-6, NULL, NULL,
        R_inv = NULL, fold_C_list = NULL, precond_fun = NULL,
        parallel = TRUE
      ),
      error = function(e) NA_real_
    )
  }, numeric(1L))

  message(sprintf("  CV done. mse range [%.4f, %.4f]",
                  min(cv_vals, na.rm=TRUE), max(cv_vals, na.rm=TRUE)))
  cv_grid <- data.frame(rho = rho_cv, log_phi = log_phi_seq,
                        phi = phi_seq, cv_mse = cv_vals)
  saveRDS(cv_grid, "likelihood_surface_albedo_cv_grid.rds")
  message("  cv_grid saved.")
}

# ------------------------------------------------------------------------------
# 5. CV optima at each rho (fast -- just phi profiling per rho)
# ------------------------------------------------------------------------------

message("\nFinding CV optima at each rho...")

if (file.exists("likelihood_surface_albedo_cv_opts.rds")) {
  message("  loading from cache...")
  cv_opts <- readRDS("likelihood_surface_albedo_cv_opts.rds")
} else {
  cv_opts <- do.call(rbind, lapply(rho_values, function(rho) {
    message(sprintf("  CV rho=%.2f", rho))
    S    <- Matrix::Diagonal(nrow(W_queen)) - rho * W_queen
    Q_sp <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
    Q_fun_fixed <- function(theta) list(Q = Q_sp, log_det_Q = NULL)

    tuned <- fastblm::tune_cv(
      y             = y_alb,
      A             = A,
      Q_fun         = Q_fun_fixed,
      theta_init    = numeric(0),
      k             = 10L,
      solver        = "cholesky",
      folds         = fold_assignments,
      log_phi_lower = log(0.005),
      log_phi_upper = log(10000),
      parallel      = TRUE,
      verbose       = FALSE
    )
    message(sprintf("    phi=%.4f  tau=%.4f", tuned$phi, 1/tuned$phi))
    data.frame(rho = rho, phi = tuned$phi, log_phi = log(tuned$phi),
               tau = 1/tuned$phi, cv_mse = tuned$value)
  }))
  saveRDS(cv_opts, "likelihood_surface_albedo_cv_opts.rds")
  message("  cv_opts saved.")
}

# ------------------------------------------------------------------------------
# 6. Optima
# ------------------------------------------------------------------------------

ml_opt <- do.call(rbind, lapply(rho_values, function(rho) {
  sub <- ml_grid[ml_grid$rho == rho & is.finite(ml_grid$ll), ]
  sub[which.max(sub$ll), c("rho", "phi", "log_phi", "ll")]
}))

cv_opt_099 <- cv_opts[cv_opts$rho == rho_cv, ]
ml_opt_099 <- ml_opt[ml_opt$rho == rho_cv, ]

# ------------------------------------------------------------------------------
# 7. Normalise ML per rho
# ------------------------------------------------------------------------------

ml_grid <- do.call(rbind, lapply(split(ml_grid, ml_grid$rho), function(d) {
  d$ll_norm <- d$ll - max(d$ll, na.rm = TRUE)
  d
}))

# ------------------------------------------------------------------------------
# 8. Combined 2-panel plot
# ------------------------------------------------------------------------------

ml_grid$rho_label <- factor(sprintf("rho=%.2f", ml_grid$rho))
ml_opt$rho_label  <- factor(sprintf("rho=%.2f", ml_opt$rho))
cv_opts$rho_label <- factor(sprintf("rho=%.2f", cv_opts$rho))

p_ml <- ggplot(ml_grid[is.finite(ml_grid$ll_norm), ],
               aes(x = log_phi, y = ll_norm, colour = rho_label)) +
  geom_line(linewidth = 0.8) +
  geom_point(data = ml_opt, aes(x = log_phi, y = 0),
             shape = 8, size = 3, show.legend = FALSE) +
  geom_vline(data = cv_opts,
             aes(xintercept = log_phi, colour = rho_label),
             linetype = "dashed", linewidth = 0.6, show.legend = FALSE) +
  scale_colour_brewer(palette = "RdYlBu", name = "rho") +
  labs(x = expression(log~phi), y = "ML log-likelihood (relative to max)",
       title = "ML likelihood profile (albedo)",
       subtitle = "Stars = ML optimum  |  Dashed = CV optimum (same colour = same rho)") +
  theme_minimal(base_size = 10)

p_cv <- ggplot(cv_grid[is.finite(cv_grid$cv_mse), ],
               aes(x = log_phi, y = cv_mse)) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  geom_vline(xintercept = cv_opt_099$log_phi, linetype = "dashed",
             colour = "steelblue", linewidth = 0.8) +
  geom_vline(xintercept = ml_opt_099$log_phi, linetype = "dashed",
             colour = "red", linewidth = 0.8) +
  annotate("text", x = cv_opt_099$log_phi, y = Inf,
           label = sprintf("CV\nphi=%.3f", cv_opt_099$phi),
           hjust = -0.1, vjust = 1.3, size = 3, colour = "steelblue") +
  annotate("text", x = ml_opt_099$log_phi, y = Inf,
           label = sprintf("ML\nphi=%.3f", ml_opt_099$phi),
           hjust = -0.1, vjust = 1.3, size = 3, colour = "red") +
  labs(x = expression(log~phi), y = "CV MSE",
       title = sprintf("CV MSE profile at rho=%.2f (albedo)", rho_cv),
       subtitle = "Blue dashed = CV optimum  |  Red dashed = ML optimum") +
  theme_minimal(base_size = 10)

print(
  p_ml / p_cv +
    plot_annotation(title = "Likelihood and CV surfaces: albedo, no covariates, no RSR")
)

# ------------------------------------------------------------------------------
# 9. 4-panel plot
# ------------------------------------------------------------------------------

panels <- lapply(rho_values, function(rho) {
  ml_sub <- ml_grid[ml_grid$rho == rho & is.finite(ml_grid$ll_norm), ]
  ml_o   <- ml_opt[ml_opt$rho == rho, ]
  cv_o   <- cv_opts[cv_opts$rho == rho, ]

  ml_sub$ll_scaled <- (ml_sub$ll_norm - min(ml_sub$ll_norm, na.rm=TRUE)) /
    diff(range(ml_sub$ll_norm, na.rm=TRUE))

  has_cv_curve <- (rho == rho_cv)

  p <- ggplot() +
    geom_line(data = ml_sub, aes(x = log_phi, y = ll_scaled),
              colour = "tomato", linewidth = 0.9)  +
    geom_point(data = ml_o, aes(x = log_phi, y = 1),
               colour = "tomato", shape = 8, size = 3)

  if (has_cv_curve) {
    cv_sub <- cv_grid[is.finite(cv_grid$cv_mse), ]
    cv_sub$cv_scaled <- 1 - (cv_sub$cv_mse - min(cv_sub$cv_mse, na.rm=TRUE)) /
      diff(range(cv_sub$cv_mse, na.rm=TRUE))
    p <- p + geom_line(data = cv_sub, aes(x = log_phi, y = cv_scaled),
                       colour = "steelblue", linewidth = 0.9)
  }

  p +
    geom_vline(xintercept = ml_o$log_phi, colour = "tomato",
               linetype = "dashed", linewidth = 0.6) +
    geom_vline(xintercept = cv_o$log_phi, colour = "steelblue",
               linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = ml_o$log_phi, y = 0.05,
             label = sprintf("ML\nphi=%.3f", ml_o$phi),
             hjust = 1.1, size = 2.8, colour = "tomato") +
    annotate("text", x = cv_o$log_phi, y = 0.15,
             label = sprintf("CV\nphi=%.3f", cv_o$phi),
             hjust = -0.1, size = 2.8, colour = "steelblue") +
    scale_y_continuous(limits = c(0, 1), labels = NULL) +
    labs(x = expression(log~phi), y = NULL,
         title = sprintf("rho = %.2f", rho),
         subtitle = if (has_cv_curve)
           "Red=ML  Blue=CV (scaled)" else "Red=ML  Blue dashed=CV opt") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank())
})

print(
  (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) +
    plot_annotation(
      title   = "ML likelihood vs CV MSE: albedo, no covariates, no RSR",
      caption = "Both scaled to [0,1]. CV curve at rho=0.99 only; CV optimum marked for all rho."
    )
)

# ------------------------------------------------------------------------------
# 10. Summary
# ------------------------------------------------------------------------------

cat("\n=== ML optima (albedo) ===\n")
ml_opt$tau <- round(1/ml_opt$phi, 2)
print(ml_opt[, c("rho", "phi", "tau")])

cat("\n=== CV optima (albedo) ===\n")
cv_opts$tau <- round(1/cv_opts$phi, 4)
print(cv_opts[, c("rho", "phi", "tau")])

future::plan(future::sequential)
