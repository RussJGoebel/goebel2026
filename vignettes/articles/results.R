# In a plain .R script, just return the data frame
ablation_summary <- function() {
  result_names <- c(
    "results_ols_baseline",
    "results_no_cov_rho1",
    "results_water_rho1",
    "results_water_rsr_rho1",
    "results_water_rsr_rinv_rho1",
    "results_water_rsr_rho_cv",
    "results_diagnostic_rho_nocov",
    "results_diagnostic_rho_norsr"
  )

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  env <- new.env(parent = emptyenv())
  for (nm in result_names) {
    tryCatch(
      utils::data(list = nm, package = "goebel2026", envir = env),
      error   = function(e) NULL,
      warning = function(e) NULL
    )
  }

  rows <- lapply(result_names, function(nm) {
    obj <- get0(nm, envir = env)
    if (is.null(obj)) return(NULL)
    data.frame(
      run_name    = nm,
      covariates  = obj$tags$covariates %||% NA,
      constraint  = obj$tags$constraint %||% NA,
      rho         = if (is.numeric(obj$rho_opt)) round(obj$rho_opt, 3) else obj$rho_opt,
      phi         = round(obj$phi, 3),
      sigma2e     = signif(obj$sigma2e, 3),
      rmse        = round(obj$rmse, 4),
      r2          = round(obj$r2, 4),
      cov_all     = round(obj$coverage_95_all,   3),
      cov_obs     = round(obj$coverage_95_obs,   3),
      cov_dense   = round(obj$coverage_95_dense, 3),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), rows))
}

# Usage
tbl <- ablation_summary()
print(tbl, digits = 4)
