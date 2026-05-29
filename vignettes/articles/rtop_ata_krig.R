# data-raw/run_rtop_atak.R
#
# Area-to-area kriging (ATAK) experiment using rtop.
# Structured as a scaling diagnostic: variogram fitting first (obs-obs only),
# then a small prediction test before committing to the full grid.
#
# Mirrors the setup from run_150m_albedo.R (Method 3 block) but replaces
# point-support gstat kriging with rtop areal kriging.
#
# Outputs:
#   timing_rtop                  -- wall times at each stage
#   results_kriging_rtop_albedo  -- full results if prediction completes
#
# Requires: rtop (install.packages("rtop"))

library(rtop)
library(sf)
library(sp)
library(goebel2026)
library(usethis)

# ------------------------------------------------------------------------------
# 0. Config
# ------------------------------------------------------------------------------

# Number of target cells to use in the small prediction test.
# Set to NULL to skip the test and go straight to the full grid.
N_PRED_TEST <- 50L

# rtop integration parameters
RTOP_PARAMS <- list(
  gDist = TRUE,   # geodesic distances (appropriate for projected coords)
  cloud = FALSE
)

# Max neighbours used in kriging (caps memory/time per prediction cell)
RTOP_NMAX <- 200L

# Only attempt full prediction if test extrapolation is under this threshold
FEASIBILITY_THRESHOLD_MIN <- 30

timing <- list()

# ------------------------------------------------------------------------------
# 1. Load setup (same as run_150m_albedo.R)
# ------------------------------------------------------------------------------

message("== Loading setup ==")
d_shared <- goebel2026::setup_shared_150m
d_albedo <- goebel2026::setup_albedo_150m

soundings_sf  <- d_shared$soundings_proj        # sf POLYGON, n soundings
fine_grid_sf  <- d_shared$fine_grid_buffered    # sf POLYGON, m grid cells
y_alb         <- d_albedo$y                     # numeric n
y_latent_true <- d_albedo$y_latent_true         # numeric m
X_obs_water   <- d_shared$X_obs_water           # n x 2 (intercept, water)
target_idx    <- which(fine_grid_sf$n_intersects > 0)

n_soundings <- length(y_alb)
n_target    <- length(target_idx)
message(sprintf("  Soundings: %d    Target cells: %d", n_soundings, n_target))

# ------------------------------------------------------------------------------
# 2. Prepare sp objects for rtop
# ------------------------------------------------------------------------------

message("\n== Preparing sp objects ==")

water_obs    <- X_obs_water[, 2]
water_target <- fine_grid_sf$proportion_water[target_idx]

# Observations
obs_sf       <- soundings_sf
obs_sf$y     <- y_alb
obs_sf$water <- water_obs
obs_sp       <- as(obs_sf, "Spatial")

# Full target grid
target_sf_full       <- fine_grid_sf[target_idx, ]
target_sf_full$water <- water_target
target_sp_full       <- as(target_sf_full, "Spatial")

# Small prediction test subset
target_sp_test <- target_sp_full[seq_len(min(N_PRED_TEST, n_target)), ]

# ------------------------------------------------------------------------------
# 3. Create rtop object and compute empirical areal variogram
#
# O(n^2) areal integrations over sounding-sounding pairs.
# predictionLocations set to test subset for now; rebuilt for full grid later.
# ------------------------------------------------------------------------------

message("\n== Step 1: Empirical areal variogram ==")
t0 <- proc.time()

rtop_obj <- tryCatch({
  rtop::createRtopObject(
    observations        = obs_sp,
    predictionLocations = target_sp_test,
    formulaString       = y ~ water,
    params              = RTOP_PARAMS
  )
}, error = function(e) {
  message("  createRtopObject FAILED: ", conditionMessage(e))
  NULL
})

if (is.null(rtop_obj)) stop("createRtopObject failed -- see message above.")

rtop_obj <- tryCatch({
  rtop::rtopVariogram(rtop_obj)
}, error = function(e) {
  message("  rtopVariogram FAILED: ", conditionMessage(e))
  NULL
})

timing$variogram_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Variogram elapsed: %.1f sec", timing$variogram_sec))

if (is.null(rtop_obj)) stop("rtopVariogram failed -- see message above.")

# ------------------------------------------------------------------------------
# 4. Fit variogram model
# ------------------------------------------------------------------------------

message("\n== Step 2: Fit variogram model ==")
t0 <- proc.time()

rtop_obj <- tryCatch({
  rtop::rtopFitVariogram(rtop_obj)
}, error = function(e) {
  message("  rtopFitVariogram FAILED: ", conditionMessage(e))
  NULL
})

timing$vgm_fit_sec <- as.numeric((proc.time() - t0)["elapsed"])
message(sprintf("  Variogram fit elapsed: %.1f sec", timing$vgm_fit_sec))

if (is.null(rtop_obj)) stop("rtopFitVariogram failed -- see message above.")

vgm_fit <- rtop_obj$variogramModel
if (!is.null(vgm_fit)) {
  message(sprintf("  nugget=%.4f  psill=%.4f  range=%.1f",
                  vgm_fit$psill[1],
                  vgm_fit$psill[2],
                  vgm_fit$range[2]))
}

# ------------------------------------------------------------------------------
# 5. Small prediction test
#
# predictionLocations is already target_sp_test from object creation.
# Times per-cell cost and extrapolates to the full grid.
# ------------------------------------------------------------------------------

message(sprintf("\n== Step 3: Prediction test (%d cells) ==", N_PRED_TEST))
t0 <- proc.time()

rtop_test <- tryCatch({
  rtop::rtopKrige(rtop_obj)
}, error = function(e) {
  message("  rtopKrige test FAILED: ", conditionMessage(e))
  NULL
})

timing$pred_test_sec     <- as.numeric((proc.time() - t0)["elapsed"])
timing$pred_test_n       <- N_PRED_TEST
timing$pred_per_cell_sec <- timing$pred_test_sec / N_PRED_TEST
timing$pred_full_est_min <- timing$pred_per_cell_sec * n_target / 60

message(sprintf("  Test elapsed: %.1f sec  (%.2f sec/cell)",
                timing$pred_test_sec, timing$pred_per_cell_sec))
message(sprintf("  Estimated full-grid time: %.1f minutes",
                timing$pred_full_est_min))

# ------------------------------------------------------------------------------
# 6. Full prediction
#
# Rebuilds rtop object with full prediction locations but reuses the
# already-fitted variogram -- no need to rerun the expensive O(n^2) step.
# Only runs if test extrapolation is under FEASIBILITY_THRESHOLD_MIN.
# ------------------------------------------------------------------------------

run_full <- !is.null(rtop_test) &&
  timing$pred_full_est_min < FEASIBILITY_THRESHOLD_MIN

if (run_full) {
  message(sprintf("\n== Step 4: Full prediction (%d cells) ==", n_target))

  rtop_full <- tryCatch({
    rtop::createRtopObject(
      observations        = obs_sp,
      predictionLocations = target_sp_full,
      formulaString       = y ~ water,
      params              = RTOP_PARAMS
    )
  }, error = function(e) {
    message("  createRtopObject (full) FAILED: ", conditionMessage(e))
    NULL
  })

  if (!is.null(rtop_full)) {
    # Reuse fitted variogram from the test object
    rtop_full$variogram      <- rtop_obj$variogram
    rtop_full$variogramModel <- rtop_obj$variogramModel

    t0 <- proc.time()
    rtop_full <- tryCatch({
      rtop::rtopKrige(rtop_full)
    }, error = function(e) {
      message("  rtopKrige full FAILED: ", conditionMessage(e))
      NULL
    })
    timing$pred_full_sec <- as.numeric((proc.time() - t0)["elapsed"])
    message(sprintf("  Full prediction elapsed: %.1f sec",
                    timing$pred_full_sec))
  }

  if (!is.null(rtop_full) && !is.null(rtop_full$predicted)) {
    pred  <- rtop_full$predicted
    mu_t  <- pred$var1.pred
    se_t  <- sqrt(pmax(pred$var1.var, 0))
    tr_t  <- y_latent_true[target_idx]
    ns_t  <- as.integer(Matrix::colSums(d_shared$A_flat[, target_idx] > 0))
    ci_lo <- mu_t - 1.96 * se_t
    ci_hi <- mu_t + 1.96 * se_t

    coverage <- function(idx) {
      if (length(idx) == 0L) return(NA_real_)
      mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
    }
    resid <- mu_t - tr_t
    rmse  <- sqrt(mean(resid^2, na.rm = TRUE))
    r2    <- 1 - sum(resid^2, na.rm = TRUE) /
      sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

    message(sprintf("  RMSE=%.4f  R2=%.4f  coverage_obs=%.3f  mean_SE=%.4f",
                    rmse, r2,
                    coverage(which(ns_t >= 1L)),
                    mean(se_t, na.rm = TRUE)))

    results_kriging_rtop_albedo <- list(
      run_name              = "kriging_rtop_albedo_150m",
      tags                  = list(resolution = 150L,
                                   method     = "atak_rtop",
                                   covariates = "water",
                                   support    = "areal",
                                   vgm_model  = "rtop"),
      timestamp             = Sys.time(),
      resolution_m          = 150L,
      posterior_mean        = mu_t,
      posterior_se          = se_t,
      ci_lower              = ci_lo,
      ci_upper              = ci_hi,
      vgm_fit               = vgm_fit,
      rmse                  = rmse,
      r2                    = r2,
      coverage_95_all       = coverage(seq_along(mu_t)),
      coverage_95_obs       = coverage(which(ns_t >= 1L)),
      coverage_95_dense     = coverage(which(ns_t >= 20L)),
      n_soundings_per_pixel = ns_t,
      timing                = timing
    )
    usethis::use_data(results_kriging_rtop_albedo, overwrite = TRUE)
  }

} else {
  message(sprintf(
    "\nSkipping full prediction: estimated %.1f min exceeds threshold of %d min.",
    timing$pred_full_est_min, FEASIBILITY_THRESHOLD_MIN
  ))
}

# ------------------------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------------------------

message("\n=== rtop ATAK timing summary ===")
for (nm in names(timing)) {
  message(sprintf("  %-30s  %.2f", nm, timing[[nm]]))
}
