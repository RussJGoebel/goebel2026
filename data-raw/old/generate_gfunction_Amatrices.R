# data-raw/build_g_A_matrices.R
#
# Builds Gaussian-weighted forward operator matrices A_g.
# Uses target_grid (not buffered grid) to match the comparison script.
#
# For each sounding i:
#   1. Find intersecting cells via st_intersects (precomputed once)
#   2. Compute intersection polygons D_i ∩ D_j
#   3. Integrate affine Gaussian g_i over each intersection polygon
#   4. Multiply by intersection area, row-normalise
#
# tau controls the Gaussian width in canonical [-1,1]^2 sounding space:
#   tau = 1/3  =>  sounding edge at ~1 sigma (recommended default)
#   tau -> inf =>  recovers uniform area-intersection A
#   tau -> 0   =>  concentrates at sounding centroid
#
# Outputs saved to data/ via usethis::use_data():
#   setup_g_A  -- list with A_g (named by tau), tau_values, A_uniform, timestamp

library(fastblm)
library(goebel2026)
library(spatintegrate)
library(sf)
library(Matrix)
library(future)
library(future.apply)
library(usethis)

FORCE_RERUN <- FALSE

future::plan(future::multisession, workers = parallel::detectCores() - 1L)
message(sprintf("Using %d workers", parallel::detectCores() - 1L))

# ------------------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------------------

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

message(sprintf("Soundings: %d  Target grid cells: %d",
                nrow(soundings_proj), nrow(target_proj)))

# ------------------------------------------------------------------------------
# 2. Helper functions
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
                      tau            = 1/3,
                      n_per_triangle = 500L,
                      verbose        = TRUE) {

  n <- nrow(soundings_sf)
  m <- nrow(fine_grid_sf)

  if (verbose) message(sprintf(
    "Building g-weighted A: n=%d, m=%d, tau=%.3f, n_per_triangle=%d, workers=%d",
    n, m, tau, n_per_triangle, future::nbrOfWorkers()
  ))

  sounding_geom <- sf::st_geometry(soundings_sf)
  fine_geom     <- sf::st_geometry(fine_grid_sf)
  crs           <- sf::st_crs(soundings_sf)

  # Precompute intersections once on main process
  touches <- sf::st_intersects(soundings_sf, fine_grid_sf, sparse = TRUE)

  # Pull internal spatintegrate functions for workers
  integrate_basis             <- spatintegrate:::integrate_basis
  .integrate_one_polygon      <- spatintegrate:::.integrate_one_polygon
  .extract_polygon_pieces     <- spatintegrate:::.extract_polygon_pieces
  .get_triangle_coords        <- spatintegrate:::.get_triangle_coords
  .tri_area_shoelace          <- spatintegrate:::.tri_area_shoelace
  .infer_k                    <- spatintegrate:::.infer_k
  map_unit_square_to_triangle <- spatintegrate:::map_unit_square_to_triangle
  generate_qmc_unit_square    <- spatintegrate:::generate_qmc_unit_square

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
      soundings_sf                = soundings_sf,
      sounding_geom               = sounding_geom,
      fine_geom                   = fine_geom,
      crs                         = crs,
      touches                     = touches,
      tau                         = tau,
      n_per_triangle              = n_per_triangle,
      make_affine_g_fn            = make_affine_g_fn,
      integrate_basis             = integrate_basis,
      .integrate_one_polygon      = .integrate_one_polygon,
      .extract_polygon_pieces     = .extract_polygon_pieces,
      .get_triangle_coords        = .get_triangle_coords,
      .tri_area_shoelace          = .tri_area_shoelace,
      .infer_k                    = .infer_k,
      map_unit_square_to_triangle = map_unit_square_to_triangle,
      generate_qmc_unit_square    = generate_qmc_unit_square
    )
  )

  i_idx <- rep(seq_len(n), vapply(rows, function(r) length(r$j), integer(1L)))
  j_idx <- unlist(lapply(rows, `[[`, "j"))
  x_val <- unlist(lapply(rows, `[[`, "x"))

  as(Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                          dims = c(n, m)), "dgCMatrix")
}

# ------------------------------------------------------------------------------
# 3. Skip helper
# ------------------------------------------------------------------------------

.should_skip <- function(obj_name) {
  if (FORCE_RERUN) return(FALSE)
  pkg <- tryCatch(
    utils::data(list = obj_name, package = "goebel2026", envir = new.env()),
    warning = function(w) NULL, error = function(e) NULL
  )
  if (!is.null(pkg)) {
    message(sprintf("  skipping %s (already saved)", obj_name))
    return(TRUE)
  }
  FALSE
}

# ------------------------------------------------------------------------------
# 4. Sanity check on small subset
# ------------------------------------------------------------------------------

message("\n--- Sanity check: 10 soundings ---")
t_test <- system.time({
  A_test <- build_g_A(soundings_proj[1:10, ], target_proj,
                      tau = 1/3, n_per_triangle = 500L)
})
rs <- Matrix::rowSums(A_test)
message(sprintf("  %.1fs  row sums [%.4f, %.4f]",
                t_test["elapsed"], min(rs), max(rs)))

# ------------------------------------------------------------------------------
# 5. Full build
# ------------------------------------------------------------------------------

if (!.should_skip("setup_g_A")) {
  tau_values <- c(1/5, 1/3, 1/2)

  # Uniform baseline -- built on same target_proj for fair comparison
  message("\n--- Building uniform A (baseline) ---")
  A_uniform <- as(
    spatintegrate::compute_overlap_fractions(soundings_proj, target_proj),
    "dgCMatrix"
  )
  rs <- Matrix::rowSums(A_uniform)
  message(sprintf("  A_uniform row sums [%.4f, %.4f]", min(rs), max(rs)))

  # g-weighted matrices
  A_g_list <- vector("list", length(tau_values))
  names(A_g_list) <- paste0("tau_", round(tau_values, 3))

  for (i in seq_along(tau_values)) {
    tau <- tau_values[i]
    nm  <- names(A_g_list)[i]
    message(sprintf("\n--- Building A_g: %s ---", nm))

    t_start <- proc.time()["elapsed"]
    A_g_list[[nm]] <- build_g_A(
      soundings_sf   = soundings_proj,
      fine_grid_sf   = target_proj,
      tau            = tau,
      n_per_triangle = 500L,
      verbose        = TRUE
    )
    elapsed <- proc.time()["elapsed"] - t_start

    rs   <- Matrix::rowSums(A_g_list[[nm]])
    diff <- norm(A_g_list[[nm]] - A_uniform, type = "F")
    message(sprintf("  done in %.1fs  row sums [%.4f, %.4f]  Frobenius diff = %.4f",
                    elapsed, min(rs), max(rs), diff))
  }

  setup_g_A <- list(
    A_g        = A_g_list,
    A_uniform  = A_uniform,
    tau_values = tau_values,
    timestamp  = Sys.time()
  )

  usethis::use_data(setup_g_A, overwrite = TRUE)
  message("\nsetup_g_A saved to data/.")
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nDone.")
