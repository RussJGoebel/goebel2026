# results.R
#
# Loads all available SIF and albedo results and prints summary tables.
# Run interactively to explore phi/rho selections across tuning methods.

library(goebel2026)
library(knitr)

`%||%` <- function(x, y) if (!is.null(x)) x else y

# ------------------------------------------------------------------------------
# 1. All result names
# ------------------------------------------------------------------------------

albedo_names <- c(
  # Main ablation
  "results_ols_baseline",
  "results_no_cov_rho1",
  "results_water_rho1",
  "results_water_rsr_rho1",
  "results_water_rsr_rinv_rho1",
  "results_water_rsr_rho_cv",
  "results_water_zero",
  # Supplement: blocking
  "results_cv_blocked",
  # Supplement: soft RSR
  "results_softrsr_albedo",
  # Supplement: neighbors
  "results_neighbor_rook",
  "results_neighbor_lc",
  "results_neighbor_lc_eps",
  "results_neighbor_lc_alpha",
  # Supplement: diagnostics
  "results_diagnostic_rho_nocov",
  "results_diagnostic_rho_norsr",
  # Supplement: ML
  "results_water_ml",
  "results_water_ml_rho1",
  "results_nocov_ml",
  "results_water_ml_zerolog"
)

sif_names <- c(
  # Main
  "results_sif_canonical",
  "results_sif_rsr",
  "results_sif_rsr_rho1",
  "results_sif_rsr_rho1_rinv",
  "results_sif_rsr_waterzero_rho1",
  "results_sif_water_zero",
  # Supplement: ML
  "results_sif_water_ml",
  "results_sif_water_ml_rho1",
  "results_sif_nocov_ml",
  "results_sif_water_ml_zerolog",
  # Supplement: ML + R_inv
  "results_sif_water_ml_rinv",
  "results_sif_water_ml_rho1_rinv",
  "results_sif_nocov_ml_rinv",
  # Supplement: rho=0.95 comparison
  "results_sif_rho095_ml",
  "results_sif_rho095_cv",
  "results_sif_rho095_cv_blocked",
  # Supplement: blocked CV
  "results_sif_cv_blocked",
  "results_sif_cv_blocked_rho_cv",
  # Supplement: CV water rho tuned
  "results_sif_water_rho_cv"
)

# ------------------------------------------------------------------------------
# 2. Load helper
# ------------------------------------------------------------------------------

load_result <- function(nm) {
  env <- new.env(parent = emptyenv())
  tryCatch(
    data(list = nm, package = "goebel2026", envir = env),
    warning = function(w) NULL,
    error   = function(e) NULL
  )
  env[[nm]]
}

# ------------------------------------------------------------------------------
# 3. Summarise helper
# ------------------------------------------------------------------------------

summarise_result <- function(res, dataset) {
  if (is.null(res)) return(NULL)
  tags <- res$tags %||% list()

  # tau = 1/phi (smoothing parameter in original parameterisation)
  phi <- res$phi %||% NA_real_
  tau <- if (!is.na(phi) && phi > 0) round(1/phi, 3) else NA_real_

  data.frame(
    dataset    = dataset,
    run_name   = res$run_name        %||% NA_character_,
    tuning     = tags$tuning         %||% NA_character_,
    covariates = tags$covariates     %||% NA_character_,
    constraint = tags$constraint     %||% NA_character_,
    R_inv      = !is.null(tags$R_inv),
    rho_opt    = round(res$rho_opt   %||% NA_real_, 4),
    phi        = round(phi,          4),
    tau        = tau,
    rmse       = round(res$rmse      %||% NA_real_, 5),
    r2         = round(res$r2        %||% NA_real_, 4),
    coverage   = round(res$coverage_95_obs %||% NA_real_, 4),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Load and build tables
# ------------------------------------------------------------------------------

albedo_list <- lapply(albedo_names, load_result)
names(albedo_list) <- albedo_names
albedo_list <- Filter(Negate(is.null), albedo_list)
message(sprintf("Loaded %d/%d albedo results", length(albedo_list), length(albedo_names)))

sif_list <- lapply(sif_names, load_result)
names(sif_list) <- sif_names
sif_list <- Filter(Negate(is.null), sif_list)
message(sprintf("Loaded %d/%d SIF results", length(sif_list), length(sif_names)))

albedo_df <- do.call(rbind, lapply(albedo_list, summarise_result, dataset = "albedo"))
sif_df    <- do.call(rbind, lapply(sif_list,    summarise_result, dataset = "sif"))
rownames(albedo_df) <- NULL
rownames(sif_df)    <- NULL

# ------------------------------------------------------------------------------
# 5. Print tables
# ------------------------------------------------------------------------------

col_names <- c("Dataset", "Run", "Tuning", "Covariates", "Constraint",
               "R_inv", "rho_opt", "phi", "tau", "RMSE", "R2", "Coverage")

cat("\n=== ALBEDO RESULTS ===\n\n")
knitr::kable(
  albedo_df,
  format    = "simple",
  col.names = col_names,
  caption   = "Albedo runs. tau = 1/phi."
)

cat("\n=== SIF RESULTS ===\n\n")
knitr::kable(
  sif_df[, setdiff(names(sif_df), c("rmse", "r2", "coverage"))],
  format    = "simple",
  col.names = col_names[!col_names %in% c("RMSE", "R2", "Coverage")],
  caption   = "SIF runs. No ground truth -- phi/rho only."
)

# ------------------------------------------------------------------------------
# 6. Quick ratio summaries
# ------------------------------------------------------------------------------

cat("\n=== PHI RATIOS (random CV as reference) ===\n\n")

if ("results_water_rho1" %in% names(albedo_list) &&
    "results_cv_blocked"  %in% names(albedo_list)) {
  r_random  <- albedo_list[["results_water_rho1"]]$phi
  r_blocked <- albedo_list[["results_cv_blocked"]]$phi
  r_ml      <- albedo_list[["results_water_ml"]]$phi %||% NA
  cat(sprintf("Albedo  random CV phi=%.4f  blocked CV phi=%.4f  ratio=%.1fx\n",
              r_random, r_blocked, r_random / r_blocked))
  if (!is.na(r_ml))
    cat(sprintf("Albedo  random CV phi=%.4f  ML phi=%.4f  ratio=%.1fx\n",
                r_random, r_ml, r_random / r_ml))
}

if ("results_sif_canonical" %in% names(sif_list) &&
    "results_sif_cv_blocked" %in% names(sif_list)) {
  s_random  <- sif_list[["results_sif_canonical"]]$phi
  s_blocked <- sif_list[["results_sif_cv_blocked"]]$phi
  s_ml      <- sif_list[["results_sif_water_ml"]]$phi %||% NA
  cat(sprintf("SIF     random CV phi=%.4f  blocked CV phi=%.4f  ratio=%.1fx\n",
              s_random, s_blocked, s_random / s_blocked))
  if (!is.na(s_ml))
    cat(sprintf("SIF     random CV phi=%.4f  ML phi=%.4f  ratio=%.1fx\n",
                s_random, s_ml, s_random / s_ml))
}
