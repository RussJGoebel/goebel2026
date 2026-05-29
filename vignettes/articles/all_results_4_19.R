# print_results_table.R
#
# Prints a summary of all saved results across all run_ scripts.
# Missing objects print as "—". Run interactively to check progress.

library(goebel2026)

summarise_result <- function(obj_name) {
  r <- tryCatch(
    { e <- new.env(); data(list = obj_name, package = "goebel2026", envir = e); e[[obj_name]] },
    error = function(e) NULL, warning = function(e) NULL
  )
  if (is.null(r)) {
    return(data.frame(name = obj_name, rho = NA, phi = NA,
                      beta_water = NA, rmse = NA, r2 = NA,
                      coverage_obs = NA, stringsAsFactors = FALSE))
  }
  data.frame(
    name         = obj_name,
    rho          = if (!is.null(r$rho_opt)) round(r$rho_opt, 3) else NA,
    phi          = if (!is.null(r$phi))     round(r$phi, 3)     else NA,
    beta_water   = if (!is.null(r$beta_hat) && length(r$beta_hat) >= 2)
      round(r$beta_hat[2], 4) else NA,
    rmse         = if (!is.null(r$rmse))            round(r$rmse, 4)            else NA,
    r2           = if (!is.null(r$r2))              round(r$r2, 4)              else NA,
    coverage_obs = if (!is.null(r$coverage_95_obs)) round(r$coverage_95_obs, 3) else NA,
    stringsAsFactors = FALSE
  )
}

pr <- function(label, objects) {
  cat(sprintf("\n=== %s ===\n", label))
  df <- do.call(rbind, lapply(objects, summarise_result))
  print(df, row.names = FALSE, na.print = "—")
}

# run_01: Main ablation + SIF canonical
pr("run_01: ALBEDO ABLATION", c(
  "results_ols_baseline",
  "results_no_cov_rho1",
  "results_water_rho1",
  "results_water_rsr_rho1",
  "results_water_rsr_rho_cv"
))

pr("run_01: SIF CANONICAL", c(
  "results_sif_canonical"
))

# run_02: RSR variants
pr("run_02: SIF RSR VARIANTS", c(
  "results_sif_naive",
  "results_sif_rsr_rho1",
  "results_softrsr_albedo"
))

# run_03: Neighbour weights
pr("run_03: ALBEDO NEIGHBOUR WEIGHTS", c(
  "results_water_rsr_rho1",
  "results_neighbor_rook",
  "results_neighbor_lc_eps",
  "results_neighbor_lc_alpha"
))

pr("run_03: SIF NEIGHBOUR WEIGHTS", c(
  "results_sif_canonical",
  "results_sif_neighbor_rook",
  "results_sif_neighbor_lc_eps",
  "results_sif_neighbor_lc_alpha"
))

# run_04: Hyperparameter tuning
pr("run_04: ALBEDO TUNING", c(
  "results_water_rsr_rho1",
  "results_cv_blocked",
  "results_water_ml"
))

pr("run_04: SIF TUNING", c(
  "results_sif_canonical",
  "results_sif_cv_blocked",
  "results_sif_water_ml",
  "results_sif_rho095_ml",
  "results_sif_rho095_cv"
))

# run_05: g-weighted A (SIF)
pr("run_05: SIF g-WEIGHTED A", c(
  "results_sif_canonical",
  "results_sif_gA_tau05",
  "results_sif_gA_tau033",
  "results_sif_gA_tau02",
  "results_sif_gA_tau01",
  "results_sif_gA_tau001"
))

# run_05c: g-weighted A (albedo)
pr("run_05c: ALBEDO g-WEIGHTED A", c(
  "results_water_rsr_rho_cv",
  "results_albedo_gA_tau05",
  "results_albedo_gA_tau033",
  "results_albedo_gA_tau02",
  "results_albedo_gA_tau01",
  "results_albedo_gA_centroid"
))
