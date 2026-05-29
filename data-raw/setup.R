# data-raw/setup.R
#
# Builds and saves all shared setup objects for the goebel2026 package.
# Run once before any run_*.R scripts.
#
# Sections:
#   1. Geometry         -- projected soundings and fine grid
#   2. Aggregation      -- A_flat (uniform g only; g-weighted A built in run_05*)
#   3. SAR priors       -- W_queen, W_rook, landcover-aware W; Q_fun factories
#   4. Covariates       -- X_latent_water, X_obs_water, y_albedo_latent
#   5. SIF data         -- y_sif, se_sif, R_inv_sif
#   6. Blocked folds    -- spatial k-means folds for supplement
#   7. Save             -- usethis::use_data() for all three setup objects
#
# Outputs (saved to data/):
#   setup_shared  -- geometry, A_flat, W, Q_fun, covariates
#   setup_albedo  -- semi-synthetic albedo response and truth
#   setup_sif     -- SIF response, uncertainty, R_inv
#
# Note: g-weighted A matrices (A_sg, A_sg_list) are NOT built here.
# They are constructed inline in run_05c_gA_albedo.R and run_05_gA_cv.R,
# which define build_g_A themselves and save results directly to the package.

library(goebel2026)
library(spatintegrate)
library(fastblm)
library(Matrix)
library(spdep)
library(sf)
library(future)
library(future.apply)
library(usethis)

future::plan(future::multisession, workers = parallel::detectCores() - 1L)

# ==============================================================================
# 1. Geometry
# ==============================================================================

message("== 1. Geometry ==")

soundings_proj     <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
fine_grid_buffered <- spatintegrate::ensure_projected(goebel2026::target_grid)
res_m              <- 330

# Remove outlier sounding before building any downstream objects (A_flat,
# X_obs_water, blocked_folds) so everything is consistent at n-1 soundings.
# The outlier has SIF = 4.91 (quality flag = 1, physically implausible).
sif_raw     <- sf::st_drop_geometry(goebel2026::soundings_augmented)$SIF_757nm
outlier_idx <- which.max(sif_raw)
keep_idx    <- setdiff(seq_along(sif_raw), outlier_idx)
message(sprintf("  Removing outlier sounding %d (SIF = %.3f) before building A",
                outlier_idx, sif_raw[outlier_idx]))
soundings_proj <- soundings_proj[keep_idx, ]

message(sprintf("  Soundings after filter: %d", nrow(soundings_proj)))
message(sprintf("  Grid cells: %d", nrow(fine_grid_buffered)))

# ==============================================================================
# 2. Aggregation matrices
# ==============================================================================

message("== 2. Aggregation matrices ==")

# A_flat: uniform g -- exact area intersection fractions (paper default)
message("  Building A_flat ...")
A_flat <- spatintegrate::compute_overlap_fractions(
  soundings  = soundings_proj,
  fine_grid  = fine_grid_buffered,
  parallel   = TRUE
)
message(sprintf("  A_flat: %d x %d, nnz = %d",
                nrow(A_flat), ncol(A_flat), Matrix::nnzero(A_flat)))

# ==============================================================================
# 3. SAR prior precision matrices
# ==============================================================================

message("== 3. SAR priors ==")

m <- nrow(fine_grid_buffered)

make_W_from_nb <- function(nb) {
  lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  W   <- spdep::listw2mat(lw)
  Matrix::Matrix(W, sparse = TRUE)
}

make_Q_fun <- function(W) {
  S <- Matrix::Diagonal(nrow(W)) - W
  Q <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
  function(theta) list(Q = Q, log_det_Q = NULL)
}

message("  Queen adjacency ...")
nb_queen <- spdep::poly2nb(fine_grid_buffered, queen = TRUE,  snap = res_m * 0.01)
W_queen  <- make_W_from_nb(nb_queen)
Q_fun_queen <- make_Q_fun(W_queen)

message("  Rook adjacency ...")
nb_rook  <- spdep::poly2nb(fine_grid_buffered, queen = FALSE, snap = res_m * 0.01)
W_rook   <- make_W_from_nb(nb_rook)
Q_fun_rook <- make_Q_fun(W_rook)

# ==============================================================================
# 4. Covariates
# ==============================================================================

message("== 4. Covariates ==")

X_latent_water  <- cbind(1, fine_grid_buffered$proportion_water)
X_obs_water     <- cbind(1, soundings_proj$proportion_water)
y_albedo_latent <- fine_grid_buffered$mean_albedo

set.seed(2026L)
sigma_eps <- 0.05 * sd(y_albedo_latent, na.rm = TRUE)
y_albedo  <- as.numeric(A_flat %*% y_albedo_latent) +
  rnorm(nrow(A_flat), sd = sigma_eps)

message(sprintf("  sigma_eps = %.4f", sigma_eps))

# Landcover-aware W (needs proportion_water, built here after X_latent_water)
message("  Landcover-aware W ...")
is_water <- fine_grid_buffered$proportion_water > 0.5
W_lc_entries <- lapply(seq_len(m), function(i) {
  js        <- nb_queen[[i]]
  js        <- js[js > 0L]
  same_type <- is_water[js] == is_water[i]
  js_keep   <- js[same_type]
  if (length(js_keep) == 0L) return(list(j = integer(0), x = numeric(0)))
  list(j = js_keep, x = rep(1 / length(js_keep), length(js_keep)))
})
i_lc <- rep(seq_len(m), lengths(lapply(W_lc_entries, `[[`, "j")))
j_lc <- unlist(lapply(W_lc_entries, `[[`, "j"))
x_lc <- unlist(lapply(W_lc_entries, `[[`, "x"))
W_landcover <- Matrix::sparseMatrix(i = i_lc, j = j_lc, x = x_lc,
                                    dims = c(m, m))
Q_fun_lc    <- make_Q_fun(W_landcover)

# ==============================================================================
# 5. SIF data
# ==============================================================================

message("== 5. SIF data ==")

sif_data <- sf::st_drop_geometry(goebel2026::soundings_augmented)

y_sif  <- sif_data$SIF_757nm
se_sif <- sif_data$SIF_Uncertainty_757nm

# outlier already removed from soundings_proj in section 1;
# keep_idx from section 1 applies here too
y_sif  <- y_sif[keep_idx]
se_sif <- se_sif[keep_idx]
message(sprintf("  n soundings after filter: %d", length(y_sif)))

# R_inv: diagonal with 1/se^2 weights, scaled so mean weight = 1
w_raw     <- 1 / se_sif^2
w_scaled  <- w_raw / mean(w_raw)
R_inv_sif <- Matrix::Diagonal(x = w_scaled)

message(sprintf("  SIF SE range: [%.4f, %.4f]", min(se_sif), max(se_sif)))
message(sprintf("  Weight range: [%.2f, %.2f]", min(w_scaled), max(w_scaled)))

# ==============================================================================
# 6. Spatially blocked CV folds
# ==============================================================================

message("== 6. Spatially blocked folds ==")

n_blocks    <- 10L
snd_coords  <- sf::st_coordinates(sf::st_centroid(soundings_proj))
set.seed(2026L)
km          <- kmeans(snd_coords, centers = n_blocks, nstart = 25L)
blocked_folds_sif    <- km$cluster
blocked_folds_albedo <- km$cluster

message(sprintf("  Block sizes: %s",
                paste(table(blocked_folds_sif), collapse = " ")))

# ==============================================================================
# 7. Assemble and save
# ==============================================================================

message("== 7. Saving ==")

setup_shared <- list(
  soundings_proj      = soundings_proj,
  fine_grid_buffered  = fine_grid_buffered,
  A_flat              = A_flat,
  W_queen             = W_queen,
  W_rook              = W_rook,
  W_landcover         = W_landcover,
  Q_fun_queen         = Q_fun_queen,
  Q_fun_rook          = Q_fun_rook,
  Q_fun_lc            = Q_fun_lc,
  X_latent_water      = X_latent_water,
  X_obs_water         = X_obs_water
)

setup_albedo <- list(
  y                   = y_albedo,
  y_latent_true       = y_albedo_latent,
  A                   = A_flat,
  X_obs               = X_obs_water,
  X_latent            = X_latent_water,
  blocked_folds       = blocked_folds_albedo,
  sigma_eps           = sigma_eps
)

setup_sif <- list(
  y                   = y_sif,
  se                  = se_sif,
  R_inv               = R_inv_sif,
  A                   = A_flat,
  keep_idx            = keep_idx,      # row indices into A_flat after outlier removal
  X_obs               = X_obs_water,
  X_latent            = X_latent_water,
  blocked_folds       = blocked_folds_sif
)

usethis::use_data(setup_shared, overwrite = TRUE)
usethis::use_data(setup_albedo, overwrite = TRUE)
usethis::use_data(setup_sif,    overwrite = TRUE)

future::plan(future::sequential)
message("Done. Objects saved to data/.")
