library(goebel2026)

load_result <- function(nm) {
  e <- new.env()
  tryCatch(
    { utils::data(list = nm, package = "goebel2026", envir = e); e[[nm]] },
    error = function(e) NULL
  )
}

# All known result object names
result_names <- c(
  # Albedo ablation (330m pipeline)
  "results_ols_baseline",
  "results_no_cov_rho1",
  "results_water_rho1",
  "results_water_rsr_rho1",
  "results_water_rsr_rho_cv",
  "results_water_rho_cv",
  # Albedo ablation (10m pipeline)
  "results_10m_ols_baseline",
  "results_10m_no_cov_rho1",
  "results_10m_water_rho1",
  "results_10m_water_rsr_rho1",
  "results_10m_water_rsr_rho_cv",
  # SIF
  "results_sif_canonical",
  "results_sif_naive",
  "results_sif_rsr_rho1",
  "results_sif_cv_blocked",
  # Tuning
  "results_cv_blocked",
  "results_water_ml",
  "results_sif_rho095_ml",
  "results_sif_rho095_cv",
  "results_sif_water_ml",
  # Neighbours
  "results_neighbor_rook",
  "results_neighbor_lc_eps",
  "results_neighbor_lc_alpha",
  "results_sif_neighbor_rook",
  "results_sif_neighbor_lc_eps",
  "results_sif_neighbor_lc_alpha",
  # Forward operator (albedo)
  "results_albedo_gA_tau05",
  "results_albedo_gA_tau033",
  "results_albedo_gA_tau02",
  "results_albedo_gA_tau01",
  "results_albedo_gA_centroid",
  # Forward operator (SIF)
  "results_sif_gA_tau05",
  "results_sif_gA_tau033",
  "results_sif_gA_tau02",
  "results_sif_gA_tau01",
  "results_sif_gA_tau001",
  # Kriging
  "results_kriging_albedo",
  "results_kriging_albedo_10m",
  "results_kriging_albedo_nocov_330m",
  "results_kriging_albedo_nocov_10m",
  "results_kriging_sif",
  # Contingency
  "results_contingency_uu",
  "results_contingency_ug",
  "results_contingency_gu",
  "results_contingency_gg",
  # Soft RSR
  "results_softrsr_albedo",
  # Rho sensitivity
  "results_water_rho_cv"
)

# Build summary table -- extract whatever fields exist
rows <- lapply(result_names, function(nm) {
  r <- load_result(nm)
  if (is.null(r)) {
    return(data.frame(name = nm, rmse = NA, r2 = NA,
                      coverage_obs = NA, phi = NA,
                      rho = NA, beta_water = NA,
                      mean_se = NA, status = "NOT FOUND",
                      stringsAsFactors = FALSE))
  }
  data.frame(
    name         = nm,
    rmse         = if (!is.null(r$rmse))           round(r$rmse, 5)    else NA,
    r2           = if (!is.null(r$r2))             round(r$r2, 4)      else NA,
    coverage_obs = if (!is.null(r$coverage_95_obs)) round(r$coverage_95_obs, 4) else NA,
    phi          = if (!is.null(r$phi))            round(r$phi, 4)     else NA,
    rho          = if (!is.null(r$rho_opt))        round(r$rho_opt, 4) else NA,
    beta_water   = if (!is.null(r$beta_hat) && length(r$beta_hat) >= 2)
      round(r$beta_hat[2], 4) else NA,
    mean_se      = if (!is.null(r$posterior_se))
      round(mean(r$posterior_se, na.rm = TRUE), 4) else NA,
    status       = "OK",
    stringsAsFactors = FALSE
  )
})

tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)

# Also print objects that aren't model results
message("\n--- Non-result objects ---")
other_names <- c(
  "quantile_fits_canonical",
  "likelihood_surface_albedo",
  "likelihood_surface_sif",
  "results_discretization_error"
)
for (nm in other_names) {
  r <- load_result(nm)
  if (is.null(r)) {
    message(sprintf("  %-40s  NOT FOUND", nm))
  } else {
    message(sprintf("  %-40s  OK  (class: %s)", nm, class(r)[1]))
  }
}
