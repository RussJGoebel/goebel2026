# g_downscale_comparison.R
#
# Compares downscaling results using:
#   (1) uniform area-intersection A  (current paper approach)
#   (2) g-weighted A at tau = 1/3    (moderate Gaussian sensitivity)
#
# No covariates, intrinsic SAR prior (rho=1), phi tuned via ML.
# Albedo semi-synthetic data from goebel2026.

library(fastblm)
library(spatintegrate)
library(goebel2026)
library(Matrix)
library(sf)
library(ggplot2)
library(patchwork)
library(future)
library(future.apply)

set.seed(42)

# ------------------------------------------------------------------------------
# 1. Data
# ------------------------------------------------------------------------------

noise_sd <- sd(goebel2026::target_grid$mean_albedo, na.rm = TRUE) / 20
y_alb <- goebel2026::soundings_augmented$mean_albedo +
  rnorm(length(goebel2026::soundings_augmented$mean_albedo), 0, noise_sd)

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

truth <- goebel2026::target_grid$mean_albedo

# ------------------------------------------------------------------------------
# 2. Build A matrices
# ------------------------------------------------------------------------------

make_affine_g_fn <- function(sounding_sf, tau) {
  coords <- sf::st_coordinates(sounding_sf)
  xy     <- unique(coords[, c("X", "Y")])
  if (nrow(xy) != 4L) stop("Expected 4 unique vertices, got ", nrow(xy))
  cx <- mean(xy[, "X"]); cy <- mean(xy[, "Y"])
  xy_c   <- sweep(xy, 2, c(cx, cy), "-")
  xy_ord <- xy_c[order(atan2(xy_c[, "Y"], xy_c[, "X"])), , drop = FALSE]
  canon  <- matrix(c(1,-1, 1,1, -1,1, -1,-1), ncol = 2, byrow = TRUE)
  M_inv  <- solve(t(solve(crossprod(canon), crossprod(canon, xy_ord))))
  force(tau); force(cx); force(cy); force(M_inv)
  # Returns plain numeric vector (used as basis_fn via matrix wrapper)
  function(pts) {
    dx <- pts[, 1] - cx; dy <- pts[, 2] - cy
    uv <- cbind(dx, dy) %*% t(M_inv)
    exp(-0.5 * rowSums(uv^2) / tau^2)
  }
}

build_g_A <- function(soundings_sf, fine_grid_sf,
                      tau            = 1/3,
                      n_per_triangle = 500L) {

  n <- nrow(soundings_sf); m <- nrow(fine_grid_sf)
  message(sprintf("Building g-weighted A: n=%d, m=%d, tau=%.3f, workers=%d",
                  n, m, tau, future::nbrOfWorkers()))

  sounding_geom <- sf::st_geometry(soundings_sf)
  fine_geom     <- sf::st_geometry(fine_grid_sf)
  crs           <- sf::st_crs(soundings_sf)
  touches       <- sf::st_intersects(soundings_sf, fine_grid_sf, sparse = TRUE)

  worker <- function(i) {
    js <- as.integer(touches[[i]])
    if (length(js) == 0L)
      return(list(j = integer(0), x = numeric(0)))

    inter_geoms <- suppressWarnings(
      sf::st_intersection(sounding_geom[i], fine_geom[js])
    )
    areas <- as.numeric(sf::st_area(inter_geoms))
    keep  <- !sf::st_is_empty(inter_geoms) & areas > 0
    if (!any(keep))
      return(list(j = integer(0), x = numeric(0)))

    inter_sf   <- sf::st_sf(geometry = inter_geoms[keep], crs = crs)
    areas_keep <- areas[keep]
    js_keep    <- js[keep]

    g_fn    <- make_affine_g_fn(soundings_sf[i, ], tau = tau)
    g_basis <- function(coords) matrix(g_fn(coords), ncol = 1L)

    g_ij <- as.numeric(integrate_basis(
      basis_fn       = g_basis,
      polygons_sf    = inter_sf,
      n_per_triangle = n_per_triangle
    ))

    weights <- g_ij * areas_keep
    s       <- sum(weights, na.rm = TRUE)
    if (s > 0) weights <- weights / s

    list(j = js_keep, x = weights)
  }

  rows <- future.apply::future_lapply(
    seq_len(n),
    worker,
    future.seed     = TRUE,
    future.packages = c("sf", "spatintegrate"),
    future.globals  = list(
      soundings_sf            = soundings_sf,
      fine_geom               = fine_geom,
      crs                     = crs,
      touches                 = touches,
      tau                     = tau,
      n_per_triangle          = n_per_triangle,
      make_affine_g_fn        = make_affine_g_fn,
      integrate_basis         = spatintegrate:::integrate_basis,
      .integrate_one_polygon  = spatintegrate:::.integrate_one_polygon,
      .extract_polygon_pieces = spatintegrate:::.extract_polygon_pieces,
      .get_triangle_coords    = spatintegrate:::.get_triangle_coords,
      .tri_area_shoelace      = spatintegrate:::.tri_area_shoelace,
      .infer_k                = spatintegrate:::.infer_k,
      map_unit_square_to_triangle = spatintegrate:::map_unit_square_to_triangle,
      generate_qmc_unit_square    = spatintegrate:::generate_qmc_unit_square
    )
  )

  i_idx <- rep(seq_len(n), vapply(rows, function(r) length(r$j), integer(1)))
  j_idx <- unlist(lapply(rows, `[[`, "j"))
  x_val <- unlist(lapply(rows, `[[`, "x"))

  Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val, dims = c(n, m))
}

# Build A matrices
future::plan(future::multisession, workers = 16)

A_uniform <- as(
  spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
  "dgCMatrix"
)

message("Building g-weighted A (tau = 1/3)...")
A_g <- as(build_g_A(soundings_proj, target_proj, tau = 1/3), "dgCMatrix")

future::plan(future::sequential)

p <- ncol(A_uniform)

# ------------------------------------------------------------------------------
# 3. Prior: intrinsic SAR (rho = 1)
# ------------------------------------------------------------------------------

W       <- goebel2026::make_W_matrix(goebel2026::target_grid)
IminusW <- Matrix::Diagonal(nrow(W)) - W
Q_fixed <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))
Q_fun   <- function(theta) list(Q = Q_fixed, log_det_Q = 0)

# ------------------------------------------------------------------------------
# 4. Fit: tune phi via ML then fit
# ------------------------------------------------------------------------------

fit_model <- function(A, label) {
  message(sprintf("\n--- Fitting: %s ---", label))
  tuned <- tune_ml(
    y          = y_alb,
    A          = A,
    Q_fun      = Q_fun,
    X_fixed    = NULL,
    theta_init = numeric(0),
    verbose    = TRUE
  )
  fit <- fastblm::fit_fastblm(
    y      = y_alb,
    A      = A,
    Q      = Q_fixed,
    phi    = tuned$phi,
    solver = "cholesky"
  )
  list(fit = fit, tuned = tuned, label = label)
}

result_uniform <- fit_model(A_uniform, "Uniform A")
result_g       <- fit_model(A_g,       "g-weighted A (tau=1/3)")

# ------------------------------------------------------------------------------
# 5. Compare
# ------------------------------------------------------------------------------

r_uniform <- result_uniform$fit$posterior_mean
r_g       <- result_g$fit$posterior_mean

rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
r2   <- function(x, y) summary(lm(y ~ x))$r.squared

cat("\n=== Posterior mean comparison ===\n")
cat(sprintf("Uniform A:    RMSE = %.5f   R2 = %.4f   phi = %.4f\n",
            rmse(r_uniform, truth), r2(r_uniform, truth),
            result_uniform$tuned$phi))
cat(sprintf("g-weighted A: RMSE = %.5f   R2 = %.4f   phi = %.4f\n",
            rmse(r_g, truth), r2(r_g, truth),
            result_g$tuned$phi))
cat(sprintf("\nMax abs diff:  %.6f\n", max(abs(r_uniform - r_g), na.rm = TRUE)))
cat(sprintf("Mean abs diff: %.6f\n",  mean(abs(r_uniform - r_g), na.rm = TRUE)))

# ------------------------------------------------------------------------------
# 6. Maps
# ------------------------------------------------------------------------------

coords_df        <- as.data.frame(sf::st_coordinates(sf::st_centroid(target_proj)))
grid_df          <- sf::st_drop_geometry(target_proj)
grid_df$x        <- coords_df[, 1]
grid_df$y        <- coords_df[, 2]
grid_df$truth    <- truth
grid_df$uniform  <- r_uniform
grid_df$g        <- r_g
grid_df$diff     <- r_g - r_uniform

lims <- range(c(r_uniform, r_g, truth), na.rm = TRUE)

make_map <- function(col, title, fill_scale) {
  ggplot(grid_df, aes(x, y, fill = .data[[col]])) +
    geom_raster() + fill_scale + coord_equal() +
    labs(title = title) + theme_minimal()
}

viridis_scale <- scale_fill_viridis_c(limits = lims, name = "Albedo")
diff_scale    <- scale_fill_gradient2(name = "Diff", midpoint = 0,
                                      low = "blue", mid = "white", high = "red")

print(
  make_map("truth",   "Truth",                viridis_scale) +
    make_map("uniform", "Uniform A",            viridis_scale) +
    make_map("g",       "g-weighted (tau=1/3)", viridis_scale) +
    make_map("diff",    "g - uniform",          diff_scale) +
    plot_layout(ncol = 2)
)
