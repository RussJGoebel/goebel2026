# albedo_comparison.R
#
# Loads all semi-synthetic albedo downscaling results from the goebel2026
# package and prints a kable summary table of hyperparameters and predictive
# performance metrics (RMSE, R^2, approximate 95% coverage).
#
# Ground truth: goebel2026::target_grid$mean_albedo
# Target pixels: goebel2026::target_grid$n_intersects > 0

library(goebel2026)
library(knitr)

# ------------------------------------------------------------------------------
# 1. Result object names
# ------------------------------------------------------------------------------

result_names <- c(
  "results_ols_baseline",
  "results_no_cov_rho1",
  "results_water_rho1",
  "results_water_rsr_rho1",
  "results_water_rsr_rinv_rho1",
  "results_water_rsr_rho_cv",
  "results_water_zero",
  "results_softrsr_albedo",
  "results_cv_blocked",
  "results_neighbor_rook",
  "results_neighbor_lc",
  "results_neighbor_lc_eps",
  "results_neighbor_lc_alpha",
  "results_diagnostic_rho_nocov",
  "results_diagnostic_rho_norsr"
)

# ------------------------------------------------------------------------------
# 2. Load all results into a named list
# ------------------------------------------------------------------------------

results_list <- lapply(result_names, function(nm) {
  env <- new.env(parent = emptyenv())
  data(list = nm, package = "goebel2026", envir = env)
  env[[nm]]
})
names(results_list) <- result_names

# ------------------------------------------------------------------------------
# 3. Ground truth and target pixel index
# ------------------------------------------------------------------------------

target_idx <- which(goebel2026::target_grid$n_intersects > 0)
truth      <- goebel2026::target_grid$mean_albedo[target_idx]

# ------------------------------------------------------------------------------
# 4. Helper
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x)) x else y

summarise_result <- function(res) {
  tags <- res$tags %||% list()

  mu <- res$posterior_mean  # already subsetted to target pixels
  se <- res$posterior_se    # same; NULL for OLS baseline

  rmse <- sqrt(mean((truth - mu)^2))
  r2   <- summary(lm(truth ~ mu))$r.squared

  # Coverage only meaningful when posterior SE is available
  coverage <- if (!is.null(se)) {
    mean(abs(truth - mu) <= 1.96 * se)
  } else {
    NA_real_
  }

  data.frame(
    run_name   = res$run_name                %||% NA_character_,
    response   = tags$response               %||% NA_character_,
    covariates = tags$covariates             %||% NA_character_,
    constraint = tags$constraint             %||% NA_character_,
    W          = tags$W                      %||% NA_character_,
    tuning     = tags$tuning                 %||% NA_character_,
    rho_opt    = round(res$rho_opt %||% NA_real_, 4),
    phi        = round(res$phi     %||% NA_real_, 4),
    RMSE       = round(rmse,       5),
    R2         = round(r2,         4),
    Coverage95 = round(coverage,   4),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 5. Build table
# ------------------------------------------------------------------------------

summary_df <- do.call(rbind, lapply(results_list, summarise_result))
rownames(summary_df) <- NULL

# ------------------------------------------------------------------------------
# 6. Print kable
# ------------------------------------------------------------------------------

knitr::kable(
  summary_df,
  format    = "simple",
  digits    = 4,
  caption   = paste(
    "Semi-synthetic albedo downscaling: model comparison.",
    "Performance computed over pixels with at least one overlapping sounding.",
    "Coverage is NA for OLS baseline (no posterior SE)."
  ),
  col.names = c(
    "Run", "Response", "Covariates", "Constraint", "W", "Tuning",
    "rho_opt", "phi", "RMSE", "R2", "Coverage (95%)"
  )
)
