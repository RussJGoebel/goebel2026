# print_all_results.R
#
# Prints all results from the goebel2026 package as formatted tables.
# Run interactively to check all results at once.
# Missing objects print as "—".

library(goebel2026)

# ------------------------------------------------------------------------------
# Helper
# ------------------------------------------------------------------------------

load_result <- function(nm) {
  tryCatch({
    e <- new.env()
    data(list = nm, package = "goebel2026", envir = e)
    e[[nm]]
  }, error = function(e) NULL, warning = function(e) NULL)
}

fmt <- function(x, digits = 4) {
  if (is.null(x) || all(is.na(x))) return("—")
  sprintf(paste0("%.", digits, "f"), x)
}

fmt2 <- function(x) fmt(x, 2)
fmt3 <- function(x) fmt(x, 3)
fmt4 <- function(x) fmt(x, 4)

hr  <- function(w = 85) cat(strrep("-", w), "\n")
hr2 <- function(w = 85) cat(strrep("=", w), "\n")

# ------------------------------------------------------------------------------
# Table printer: albedo results (with ground truth)
# ------------------------------------------------------------------------------

print_albedo_table <- function(label, objects) {
  hr2()
  cat(sprintf(" %s\n", label))
  hr2()
  cat(sprintf("  %-35s  %5s  %7s  %6s  %6s  %8s  %8s\n",
              "Model", "rho", "phi", "RMSE", "R2", "cov_obs", "beta_w"))
  hr()
  for (nm in objects) {
    r <- load_result(nm)
    if (is.null(r)) {
      cat(sprintf("  %-35s  %s\n", nm, "— NOT FOUND —"))
    } else {
      cat(sprintf("  %-35s  %5s  %7s  %6s  %6s  %8s  %8s\n",
                  nm,
                  fmt3(r$rho_opt),
                  fmt2(r$phi),
                  fmt4(r$rmse),
                  fmt4(r$r2),
                  fmt3(r$coverage_95_obs),
                  fmt4(if (!is.null(r$beta_hat) && length(r$beta_hat) >= 2)
                    r$beta_hat[2] else NA)))
    }
  }
  cat("\n")
}

# ------------------------------------------------------------------------------
# Table printer: SIF results (no ground truth)
# ------------------------------------------------------------------------------

print_sif_table <- function(label, objects) {
  hr2()
  cat(sprintf(" %s\n", label))
  hr2()
  cat(sprintf("  %-35s  %5s  %7s  %8s  %8s  %8s\n",
              "Model", "rho", "phi", "beta_int", "beta_w", "mean_SE"))
  hr()
  for (nm in objects) {
    r <- load_result(nm)
    if (is.null(r)) {
      cat(sprintf("  %-35s  %s\n", nm, "— NOT FOUND —"))
    } else {
      cat(sprintf("  %-35s  %5s  %7s  %8s  %8s  %8s\n",
                  nm,
                  fmt3(r$rho_opt),
                  fmt4(r$phi),
                  fmt4(if (!is.null(r$beta_hat) && length(r$beta_hat) >= 1)
                    r$beta_hat[1] else NA),
                  fmt4(if (!is.null(r$beta_hat) && length(r$beta_hat) >= 2)
                    r$beta_hat[2] else NA),
                  fmt4(mean(r$posterior_se, na.rm = TRUE))))
    }
  }
  cat("\n")
}

# ==============================================================================
# 1. Main paper table: albedo ablation (Table 2 in paper)
# ==============================================================================

print_albedo_table("TABLE 1: ALBEDO ABLATION (main paper Table 2)", c(
  "results_ols_baseline",
  "results_no_cov_rho1",
  "results_water_rho1",
  "results_water_rsr_rho1",
  "results_water_rsr_rho_cv"
))

# ==============================================================================
# 2. SIF canonical
# ==============================================================================

print_sif_table("TABLE 2: SIF CANONICAL", c(
  "results_sif_canonical"
))

# ==============================================================================
# 3. Supplement: rho sensitivity (albedo)
# ==============================================================================

print_albedo_table("TABLE 3 (SUPP): RHO SENSITIVITY -- ALBEDO", c(
  "results_water_rsr_rho1",
  "results_water_rsr_rho_cv",
  "results_water_rho1",
  "results_water_rho_cv"
))

# ==============================================================================
# 4. Supplement: RSR variants (SIF)
# ==============================================================================

print_sif_table("TABLE 4 (SUPP): RSR VARIANTS -- SIF", c(
  "results_sif_naive",
  "results_sif_rsr_rho1",
  "results_sif_canonical"
))

# ==============================================================================
# 5. Supplement: neighbour weights (albedo)
# ==============================================================================

print_albedo_table("TABLE 5 (SUPP): NEIGHBOUR WEIGHTS -- ALBEDO", c(
  "results_water_rsr_rho1",
  "results_neighbor_rook",
  "results_neighbor_lc_eps",
  "results_neighbor_lc_alpha"
))

# ==============================================================================
# 6. Supplement: neighbour weights (SIF)
# ==============================================================================

print_sif_table("TABLE 6 (SUPP): NEIGHBOUR WEIGHTS -- SIF", c(
  "results_sif_canonical",
  "results_sif_neighbor_rook",
  "results_sif_neighbor_lc_eps",
  "results_sif_neighbor_lc_alpha"
))

# ==============================================================================
# 7. Supplement: hyperparameter tuning (albedo)
# ==============================================================================

print_albedo_table("TABLE 7 (SUPP): HYPERPARAMETER TUNING -- ALBEDO", c(
  "results_water_rsr_rho1",       # random CV baseline
  "results_cv_blocked",           # blocked CV
  "results_water_ml"              # ML
))

# ==============================================================================
# 8. Supplement: hyperparameter tuning (SIF)
# ==============================================================================

print_sif_table("TABLE 8 (SUPP): HYPERPARAMETER TUNING -- SIF", c(
  "results_sif_canonical",        # random CV baseline
  "results_sif_cv_blocked",
  "results_sif_water_ml",
  "results_sif_rho095_ml",
  "results_sif_rho095_cv"
))

# ==============================================================================
# 9. Supplement: forward operator sensitivity -- g-weighted A (SIF)
# ==============================================================================

print_sif_table("TABLE 9 (SUPP): G-WEIGHTED A -- SIF", c(
  "results_sif_canonical",
  "results_sif_gA_tau05",
  "results_sif_gA_tau033",
  "results_sif_gA_tau02",
  "results_sif_gA_tau01",
  "results_sif_gA_tau001"
))

# ==============================================================================
# 10. Supplement: forward operator sensitivity -- g-weighted A (albedo)
# ==============================================================================

print_albedo_table("TABLE 10 (SUPP): G-WEIGHTED A -- ALBEDO", c(
  "results_water_rsr_rho_cv",
  "results_albedo_gA_tau05",
  "results_albedo_gA_tau033",
  "results_albedo_gA_tau02",
  "results_albedo_gA_tau01",
  "results_albedo_gA_centroid"
))

# ==============================================================================
# 11. Supplement: kriging comparison
# ==============================================================================

hr2()
cat(" TABLE 11 (SUPP): KRIGING COMPARISON\n")
hr2()
cat(sprintf("  %-30s  %6s  %6s  %8s  %8s\n",
            "Model", "RMSE", "R2", "cov_obs", "mean_SE"))
hr(75)

for (nm in c("results_water_rsr_rho_cv", "results_albedo_gA_centroid",
             "results_kriging_albedo")) {
  r <- load_result(nm)
  if (is.null(r)) {
    cat(sprintf("  %-30s  %s\n", nm, "— NOT FOUND —"))
  } else {
    cat(sprintf("  %-30s  %6s  %6s  %8s  %8s\n",
                nm,
                fmt4(r$rmse),
                fmt4(r$r2),
                fmt3(r$coverage_95_obs),
                fmt4(mean(r$posterior_se, na.rm = TRUE))))
  }
}

cat("\n")
cat("SIF comparison (no ground truth -- SE only):\n")
for (nm in c("results_sif_canonical", "results_kriging_sif")) {
  r <- load_result(nm)
  if (!is.null(r))
    cat(sprintf("  %-30s  mean_pred=%s  mean_SE=%s\n",
                nm,
                fmt4(mean(r$posterior_mean, na.rm = TRUE)),
                fmt4(mean(r$posterior_se,   na.rm = TRUE))))
}

# ==============================================================================
# 12. Soft RSR
# ==============================================================================

print_albedo_table("TABLE 12 (SUPP): SOFT RSR -- ALBEDO", c(
  "results_water_rsr_rho1",
  "results_softrsr_albedo"
))

# ==============================================================================
# Summary: what's missing
# ==============================================================================

hr2()
cat(" MISSING OBJECTS\n")
hr2()

all_objects <- c(
  "results_ols_baseline", "results_no_cov_rho1", "results_water_rho1",
  "results_water_rsr_rho1", "results_water_rsr_rho_cv", "results_water_rho_cv",
  "results_sif_naive", "results_sif_rsr_rho1", "results_sif_canonical",
  "results_softrsr_albedo",
  "results_neighbor_rook", "results_neighbor_lc_eps", "results_neighbor_lc_alpha",
  "results_sif_neighbor_rook", "results_sif_neighbor_lc_eps", "results_sif_neighbor_lc_alpha",
  "results_cv_blocked", "results_water_ml",
  "results_sif_cv_blocked", "results_sif_water_ml", "results_sif_rho095_ml", "results_sif_rho095_cv",
  "results_sif_gA_tau05", "results_sif_gA_tau033", "results_sif_gA_tau02",
  "results_sif_gA_tau01", "results_sif_gA_tau001",
  "results_albedo_gA_tau05", "results_albedo_gA_tau033", "results_albedo_gA_tau02",
  "results_albedo_gA_tau01", "results_albedo_gA_centroid",
  "results_kriging_albedo", "results_kriging_sif"
)

missing <- vapply(all_objects, function(nm) is.null(load_result(nm)), logical(1))
if (any(missing)) {
  cat("  The following objects were not found in goebel2026:\n")
  for (nm in names(missing)[missing]) cat(sprintf("    - %s\n", nm))
} else {
  cat("  All objects present!\n")
}
cat("\n")
