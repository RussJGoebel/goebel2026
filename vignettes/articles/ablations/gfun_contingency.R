# g_contingency.R
#
# 2x2 contingency experiment:
#   Rows    = data generating process (uniform A vs g-weighted A)
#   Columns = fitting model           (uniform A vs g-weighted A)
#
# For each of the 4 cells we report RMSE, R2, and a posterior mean map.
# The diagonal cells should outperform the off-diagonal ones if the
# choice of A matters.

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
# 1. Setup
# ------------------------------------------------------------------------------

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)
truth          <- goebel2026::target_grid$mean_albedo

W       <- goebel2026::make_W_matrix(goebel2026::target_grid)
IminusW <- Matrix::Diagonal(nrow(W)) - W
Q_fixed <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(IminusW)))
Q_fun   <- function(theta) list(Q = Q_fixed, log_det_Q = 0)

# ------------------------------------------------------------------------------
# 2. Build both A matrices
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
  function(pts) {
    dx <- pts[, 1] - cx; dy <- pts[, 2] - cy
    uv <- cbind(dx, dy) %*% t(M_inv)
    exp(-0.5 * rowSums(uv^2) / tau^2)
  }
}

build_g_A <- function(soundings_sf, fine_grid_sf,
                      tau = 1/3, n_per_triangle = 16L) {
  n <- nrow(soundings_sf); m <- nrow(fine_grid_sf)
  message(sprintf("Building g-weighted A: n=%d, m=%d, tau=%.3f, workers=%d",
                  n, m, tau, future::nbrOfWorkers()))
  sounding_geom <- sf::st_geometry(soundings_sf)
  fine_geom     <- sf::st_geometry(fine_grid_sf)
  crs           <- sf::st_crs(soundings_sf)
  touches       <- sf::st_intersects(soundings_sf, fine_grid_sf, sparse = TRUE)

  worker <- function(i) {
    js <- as.integer(touches[[i]])
    if (length(js) == 0L) return(list(j = integer(0), x = numeric(0)))
    inter_geoms <- suppressWarnings(
      sf::st_intersection(sounding_geom[i], fine_geom[js]))
    areas <- as.numeric(sf::st_area(inter_geoms))
    keep  <- !sf::st_is_empty(inter_geoms) & areas > 0
    if (!any(keep)) return(list(j = integer(0), x = numeric(0)))
    inter_sf   <- sf::st_sf(geometry = inter_geoms[keep], crs = crs)
    areas_keep <- areas[keep]; js_keep <- js[keep]
    g_fn    <- make_affine_g_fn(soundings_sf[i, ], tau = tau)
    g_basis <- function(coords) matrix(g_fn(coords), ncol = 1L)
    g_ij    <- as.numeric(integrate_basis(g_basis, inter_sf,
                                          n_per_triangle = n_per_triangle))
    weights <- g_ij * areas_keep
    s <- sum(weights, na.rm = TRUE)
    if (s > 0) weights <- weights / s
    list(j = js_keep, x = weights)
  }

  rows <- future.apply::future_lapply(
    seq_len(n), worker,
    future.seed = TRUE,
    future.packages = c("sf", "spatintegrate"),
    future.globals = list(
      soundings_sf = soundings_sf, fine_geom = fine_geom, crs = crs,
      touches = touches, tau = tau, n_per_triangle = n_per_triangle,
      make_affine_g_fn = make_affine_g_fn,
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

future::plan(future::multisession, workers = 16)

message("Building uniform A...")
A_unif <- as(spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
             "dgCMatrix")

message("Building g-weighted A (tau=1/3)...")
A_g <- as(build_g_A(soundings_proj, target_proj, tau = 1/3), "dgCMatrix")

future::plan(future::sequential)

# ------------------------------------------------------------------------------
# 3. Generate synthetic observations under each A
#    y = A %*% truth + noise,  noise ~ N(0, noise_sd^2)
# ------------------------------------------------------------------------------

noise_sd <- sd(truth, na.rm = TRUE) / 20

y_unif <- as.numeric(A_unif %*% truth) +
  rnorm(nrow(A_unif), 0, noise_sd)

y_g    <- as.numeric(A_g %*% truth) +
  rnorm(nrow(A_g), 0, noise_sd)

# ------------------------------------------------------------------------------
# 4. Fit function
# ------------------------------------------------------------------------------

fit_model <- function(y, A, label) {
  message(sprintf("  Fitting: %s", label))
  tuned <- tune_ml(y = y, A = A, Q_fun = Q_fun,
                   X_fixed = NULL, theta_init = numeric(0), verbose = FALSE)
  fit <- fastblm::fit_fastblm(y = y, A = A, Q = Q_fixed,
                              phi = tuned$phi, solver = "cholesky")
  list(fit = fit, tuned = tuned, label = label,
       rmse = sqrt(mean((fit$posterior_mean - truth)^2, na.rm = TRUE)),
       r2   = summary(lm(truth ~ fit$posterior_mean))$r.squared)
}

# ------------------------------------------------------------------------------
# 5. Run all 4 cells
# ------------------------------------------------------------------------------

message("\n--- Row 1: data generated with UNIFORM A ---")
r_unif_unif <- fit_model(y_unif, A_unif, "gen=uniform, fit=uniform")
r_unif_g    <- fit_model(y_unif, A_g,    "gen=uniform, fit=g-weighted")

message("\n--- Row 2: data generated with G-WEIGHTED A ---")
r_g_unif    <- fit_model(y_g, A_unif, "gen=g-weighted, fit=uniform")
r_g_g       <- fit_model(y_g, A_g,    "gen=g-weighted, fit=g-weighted")

# ------------------------------------------------------------------------------
# 6. Print summary table
# ------------------------------------------------------------------------------

# Intersecting cells only for metrics
covered <- Matrix::colSums(A_unif) > 0 & Matrix::colSums(A_g) > 0
truth_covered <- truth[covered]

recompute_metrics <- function(result) {
  pred <- result$fit$posterior_mean[covered]
  result$rmse <- sqrt(mean((pred - truth_covered)^2, na.rm = TRUE))
  result$r2   <- summary(lm(truth_covered ~ pred))$r.squared
  result
}

r_unif_unif <- recompute_metrics(r_unif_unif)
r_unif_g    <- recompute_metrics(r_unif_g)
r_g_unif    <- recompute_metrics(r_g_unif)
r_g_g       <- recompute_metrics(r_g_g)

# ------------------------------------------------------------------------------
# 6. Print summary table
# ------------------------------------------------------------------------------

cat("\n=== 2x2 Contingency Table ===\n")
cat(sprintf("%-35s  RMSE      R2      phi\n", ""))
for (r in list(r_unif_unif, r_unif_g, r_g_unif, r_g_g)) {
  cat(sprintf("%-35s  %.5f  %.4f  %.4f\n",
              r$label, r$rmse, r$r2, r$tuned$phi))
}

# ------------------------------------------------------------------------------
# 7. 2x2 map grid
# ------------------------------------------------------------------------------

coords_df <- as.data.frame(sf::st_coordinates(sf::st_centroid(target_proj)))

# Shared color scale across all four panels
all_vals <- c(
  r_unif_unif$fit$posterior_mean[covered],
  r_unif_g$fit$posterior_mean[covered],
  r_g_unif$fit$posterior_mean[covered],
  r_g_g$fit$posterior_mean[covered]
)
shared_lims <- range(all_vals, na.rm = TRUE)

make_map <- function(result, row_label, col_label) {
  df <- data.frame(
    x   = coords_df[covered, 1],
    y   = coords_df[covered, 2],
    val = result$fit$posterior_mean[covered]
  )
  subtitle <- sprintf("RMSE=%.4f  R²=%.3f", result$rmse, result$r2)
  ggplot(df, aes(x, y, fill = val)) +
    geom_raster() +
    scale_fill_viridis_c(name = "Albedo", limits = shared_lims) +
    coord_equal() +
    labs(title    = sprintf("Gen: %s | Fit: %s", row_label, col_label),
         subtitle = subtitle) +
    theme_minimal(base_size = 9) +
    theme(legend.position = "right")
}

p1 <- make_map(r_unif_unif, "Uniform",    "Uniform")
p2 <- make_map(r_unif_g,    "Uniform",    "g-weighted")
p3 <- make_map(r_g_unif,    "g-weighted", "Uniform")
p4 <- make_map(r_g_g,       "g-weighted", "g-weighted")

print(
  (p1 | p2) / (p3 | p4) +
    plot_annotation(
      title   = "2x2 Contingency: Data-generating A vs Fitting A",
      caption = "Diagonal = correctly specified; off-diagonal = misspecified"
    )
)
