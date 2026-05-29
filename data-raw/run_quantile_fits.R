# data-raw/run_07_quantile_fits.R
#
# Computes binned quantile fits of posterior SIF vs NDVI for Fig 5.
# Uses results_sif_canonical (posterior mean + SE) and target_grid (NDVI).
#
# Approach: bin NDVI into equal-count bins, compute median SIF per bin
# across n_samples posterior draws. Produces uncertainty bands from
# sample-to-sample variation.
#
# Outputs saved via usethis::use_data():
#   quantile_fits_canonical  -- list with binned quantiles, bin edges, R2 samples
#                               including r2_mean_dense / r2_sd_dense for pixels
#                               with 20+ overlapping soundings

library(goebel2026)
library(usethis)

FORCE_RERUN <- TRUE

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  if (exists(obj_name, envir = .GlobalEnv)) {
    message(sprintf("  skipping %s (already in environment)", obj_name)); return(TRUE)
  }
  pkg_data <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    error = function(e) NULL, warning = function(e) NULL
  )
  if (!is.null(pkg_data)) {
    message(sprintf("  skipping %s (already in data/)", obj_name)); return(TRUE)
  }
  FALSE
}

if (.should_skip("quantile_fits_canonical")) {
  message("quantile_fits_canonical already saved -- done.")
  stop("Nothing to do.", call. = FALSE)
}

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

message("Loading data...")

sif    <- goebel2026::results_sif_canonical
ndvi   <- goebel2026::target_grid$mean_ndvi
target_idx <- which(goebel2026::target_grid$n_intersects > 0)

preds  <- sif$posterior_mean   # length = n target pixels
se     <- sif$posterior_se

ndvi_t <- ndvi[target_idx]     # NDVI at target pixels only

# Dense pixel mask: 20+ overlapping soundings
dense_mask <- goebel2026::target_grid$n_intersects[target_idx] > 20
message(sprintf("  n dense pixels (20+ soundings): %d / %d",
                sum(dense_mask), length(preds)))

# Also need sounding-level SIF and NDVI for left panel (raw observations)
y_sif         <- goebel2026::setup_sif$y
soundings_aug <- goebel2026::soundings_augmented

# filter

soundings_aug <- soundings_aug[-which.max(soundings_aug$SIF_757nm),]

# NDVI at each sounding -- use footprint-averaged mean_ndvi (same spatial support as SIF)
ndvi_soundings <- soundings_aug$mean_ndvi

message(sprintf("  n target pixels: %d  n soundings: %d",
                length(preds), length(y_sif)))

# ------------------------------------------------------------------------------
# 2. Binning setup
# ------------------------------------------------------------------------------

N_BINS    <- 25L      # number of equal-count NDVI bins
N_SAMPLES <- 50L      # posterior draws for uncertainty
QUANTILES <- c(0.025, 0.25, 0.5, 0.75, 0.975)  # inner + outer bands
set.seed(2026L)

# Equal-count bins on NDVI (target pixels)
ndvi_finite  <- ndvi_t[is.finite(ndvi_t)]
bin_breaks   <- quantile(ndvi_finite, probs = seq(0, 1, length.out = N_BINS + 1L),
                         na.rm = TRUE)
bin_breaks   <- unique(bin_breaks)   # remove duplicates at extremes
n_bins_actual <- length(bin_breaks) - 1L

bin_assign <- function(x) {
  findInterval(x, bin_breaks, rightmost.closed = TRUE)
}

bin_idx <- bin_assign(ndvi_t)
bin_mids <- 0.5 * (bin_breaks[-length(bin_breaks)] + bin_breaks[-1L])

# Same bins for sounding-level NDVI
bin_idx_snd <- bin_assign(ndvi_soundings)

# ------------------------------------------------------------------------------
# 3. Left panel: raw SIF soundings -- SIF bins, NDVI quantiles (horizontal bands)
# ------------------------------------------------------------------------------

message("Computing sounding-level binned quantiles (SIF bins, NDVI quantiles)...")

# Equal-count bins on SIF (soundings)
sif_finite   <- y_sif[is.finite(y_sif) & is.finite(ndvi_soundings)]
sif_breaks   <- quantile(sif_finite, probs = seq(0, 1, length.out = N_BINS + 1L),
                         na.rm = TRUE)
sif_breaks   <- unique(sif_breaks)
n_sif_bins   <- length(sif_breaks) - 1L
sif_bin_mids <- 0.5 * (sif_breaks[-length(sif_breaks)] + sif_breaks[-1L])

bin_idx_sif <- findInterval(y_sif, sif_breaks, rightmost.closed = TRUE)

snd_quantiles <- array(NA_real_,
                       dim = c(n_sif_bins, length(QUANTILES)),
                       dimnames = list(NULL, paste0("q", QUANTILES)))
snd_counts <- integer(n_sif_bins)

for (b in seq_len(n_sif_bins)) {
  idx_b <- which(bin_idx_sif == b & is.finite(ndvi_soundings))
  snd_counts[b] <- length(idx_b)
  if (length(idx_b) < 3L) next
  snd_quantiles[b, ] <- quantile(ndvi_soundings[idx_b], probs = QUANTILES, na.rm = TRUE)
}

# ------------------------------------------------------------------------------
# 4. Right panel: posterior SIF -- binned quantiles across samples
# ------------------------------------------------------------------------------

message(sprintf("Computing posterior binned quantiles (%d samples x %d bins)...",
                N_SAMPLES, n_bins_actual))

# Array: [sample, bin, quantile]
post_quantiles <- array(NA_real_,
                        dim = c(N_SAMPLES, n_bins_actual, length(QUANTILES)),
                        dimnames = list(NULL, NULL, paste0("q", QUANTILES)))

# R2 per sample (all pixels) and dense pixels (20+ soundings)
r2_samples       <- numeric(N_SAMPLES)
r2_samples_dense <- numeric(N_SAMPLES)

for (jj in seq_len(N_SAMPLES)) {
  if (jj %% 10 == 0) message(sprintf("  sample %d / %d", jj, N_SAMPLES))

  # Draw from posterior predictive
  pred_sample <- rnorm(length(preds), mean = preds, sd = se)

  # Binned quantiles
  for (b in seq_len(n_bins_actual)) {
    idx_b <- which(bin_idx == b)
    if (length(idx_b) < 3L) next
    post_quantiles[jj, b, ] <- quantile(pred_sample[idx_b],
                                        probs = QUANTILES, na.rm = TRUE)
  }

  # R2 of this sample vs NDVI (all pixels)
  finite_idx <- is.finite(pred_sample) & is.finite(ndvi_t)
  if (sum(finite_idx) > 2L) {
    lm_fit <- lm(pred_sample[finite_idx] ~ ndvi_t[finite_idx])
    ss_res <- sum(residuals(lm_fit)^2)
    ss_tot <- sum((pred_sample[finite_idx] -
                     mean(pred_sample[finite_idx]))^2)
    r2_samples[jj] <- 1 - ss_res / ss_tot
  }

  # R2 of this sample vs NDVI (dense pixels only: 20+ soundings)
  finite_dense <- is.finite(pred_sample) & is.finite(ndvi_t) & dense_mask
  if (sum(finite_dense) > 2L) {
    lm_dense <- lm(pred_sample[finite_dense] ~ ndvi_t[finite_dense])
    ss_res_d <- sum(residuals(lm_dense)^2)
    ss_tot_d <- sum((pred_sample[finite_dense] -
                       mean(pred_sample[finite_dense]))^2)
    r2_samples_dense[jj] <- 1 - ss_res_d / ss_tot_d
  }
}

# Summary: mean and SD across samples per bin per quantile
post_q_mean <- apply(post_quantiles, c(2, 3), mean, na.rm = TRUE)
post_q_sd   <- apply(post_quantiles, c(2, 3), sd,   na.rm = TRUE)

r2_mean <- mean(r2_samples, na.rm = TRUE)
r2_sd   <- sd(r2_samples,   na.rm = TRUE)

r2_mean_dense <- mean(r2_samples_dense, na.rm = TRUE)
r2_sd_dense   <- sd(r2_samples_dense,   na.rm = TRUE)

message(sprintf("  Posterior R2 (all)   = %.3f +/- %.3f", r2_mean,       r2_sd))
message(sprintf("  Posterior R2 (dense) = %.3f +/- %.3f", r2_mean_dense, r2_sd_dense))

# ------------------------------------------------------------------------------
# 5. R2 for raw soundings vs NDVI
# ------------------------------------------------------------------------------

finite_snd <- is.finite(y_sif) & is.finite(ndvi_soundings)
lm_snd <- lm(y_sif[finite_snd] ~ ndvi_soundings[finite_snd])
ss_res_snd <- sum(residuals(lm_snd)^2)
ss_tot_snd <- sum((y_sif[finite_snd] - mean(y_sif[finite_snd]))^2)
r2_soundings <- 1 - ss_res_snd / ss_tot_snd

message(sprintf("  Sounding R2 = %.3f", r2_soundings))

# ------------------------------------------------------------------------------
# 6. Save
# ------------------------------------------------------------------------------

quantile_fits_canonical <- list(
  # Binning
  bin_breaks    = bin_breaks,
  bin_mids      = bin_mids,
  n_bins        = n_bins_actual,
  quantiles     = QUANTILES,
  n_samples     = N_SAMPLES,

  # Left panel: raw soundings -- SIF bins, NDVI quantiles (horizontal bands)
  snd_quantiles  = snd_quantiles,   # [sif_bin x quantile] -- NDVI quantiles
  snd_counts     = snd_counts,
  sif_bin_mids   = sif_bin_mids,
  sif_breaks     = sif_breaks,
  r2_soundings   = r2_soundings,

  # Right panel: posterior
  post_q_mean   = post_q_mean,     # [bin x quantile] -- mean across samples
  post_q_sd     = post_q_sd,       # [bin x quantile] -- SD across samples
  post_quantiles = post_quantiles, # [sample x bin x quantile] -- full array
  r2_samples    = r2_samples,
  r2_mean       = r2_mean,
  r2_sd         = r2_sd,

  # Dense pixels (20+ overlapping soundings)
  r2_samples_dense = r2_samples_dense,
  r2_mean_dense    = r2_mean_dense,
  r2_sd_dense      = r2_sd_dense,

  # Water pixel flags (>50% water coverage shown as crosses in figure)
  is_water_target   = goebel2026::target_grid$proportion_water[target_idx] > 0.5,
  is_water_soundings = soundings_aug$proportion_water > 0.5,

  # Raw data for scatter
  ndvi_target   = ndvi_t,
  preds         = preds,
  se            = se,
  ndvi_soundings = ndvi_soundings,
  y_sif         = y_sif,

  timestamp     = Sys.time()
)

usethis::use_data(quantile_fits_canonical, overwrite = TRUE)

message("\nrun_07_quantile_fits.R complete.")
message("  quantile_fits_canonical saved to data/")
message(sprintf("  R2 soundings:        %.3f", r2_soundings))
message(sprintf("  R2 posterior (all):  %.3f +/- %.3f", r2_mean,       r2_sd))
message(sprintf("  R2 posterior (dense):%.3f +/- %.3f", r2_mean_dense, r2_sd_dense))
