# sif_comparison.R
#
# Summary table of SIF run hyperparameter selections (phi, rho).
# No performance metrics -- SIF runs have no ground truth.

library(goebel2026)
library(knitr)

# ------------------------------------------------------------------------------
# 1. Result names
# ------------------------------------------------------------------------------

result_names <- c(
  # Main / canonical runs
  "results_sif_canonical",
  "results_sif_rsr",
  "results_sif_rsr_rho1",
  "results_sif_rsr_rho1_rinv",
  "results_sif_rsr_waterzero_rho1",
  "results_sif_water_zero",
  # ML runs
  "results_sif_water_ml",
  "results_sif_water_ml_rho1",
  "results_sif_nocov_ml"
)

# ------------------------------------------------------------------------------
# 2. Load
# ------------------------------------------------------------------------------

`%||%` <- function(x, y) if (!is.null(x)) x else y

results_list <- lapply(result_names, function(nm) {
  env <- new.env(parent = emptyenv())
  tryCatch(
    data(list = nm, package = "goebel2026", envir = env),
    warning = function(w) message(sprintf("WARNING loading %s: %s", nm, conditionMessage(w))),
    error   = function(e) message(sprintf("ERROR loading %s: %s",   nm, conditionMessage(e)))
  )
  obj <- env[[nm]]
  if (is.null(obj)) message(sprintf("  -> %s not found, skipping.", nm))
  obj
})
names(results_list) <- result_names
results_list <- Filter(Negate(is.null), results_list)
message(sprintf("Loaded %d of %d objects.", length(results_list), length(result_names)))

# ------------------------------------------------------------------------------
# 3. Summarise
# ------------------------------------------------------------------------------

fmt <- function(x, digits = 4) {
  if (is.null(x) || length(x) != 1L || !is.numeric(x)) NA_real_
  else round(x, digits)
}

summarise_sif <- function(res) {
  tags <- res$tags %||% list()
  data.frame(
    run_name   = res$run_name            %||% NA_character_,
    covariates = tags$covariates         %||% NA_character_,
    constraint = tags$constraint         %||% NA_character_,
    tuning     = tags$tuning             %||% NA_character_,
    rho_opt    = fmt(res$rho_opt,  4),
    phi        = fmt(res$phi,      2),
    sigma2e    = fmt(res$sigma2e,  6),
    stringsAsFactors = FALSE
  )
}

tab <- do.call(rbind, lapply(results_list, summarise_sif))
rownames(tab) <- NULL

# ------------------------------------------------------------------------------
# 4. Print
# ------------------------------------------------------------------------------

knitr::kable(
  tab,
  format    = "simple",
  col.names = c("Run", "Covariates", "Constraint", "Tuning",
                "rho_opt", "phi", "sigma2e"),
  caption   = "SIF downscaling runs: hyperparameter selections. No ground truth available."
)
