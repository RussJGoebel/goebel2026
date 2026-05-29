# data-raw/run_05c_gA_albedo.R
#
# Forward operator sensitivity on semi-synthetic albedo (ground truth known).
# Compares uniform A vs g-weighted A variants vs centroid limit.
# Having ground truth lets us compute RMSE/R2/coverage directly.
#
# All runs: albedo, water covariate + RSR + rho CV-tuned.
#
# Outputs saved to data/ via usethis::use_data():
#   results_albedo_gA_tau05    -- tau=0.5
#   results_albedo_gA_tau033   -- tau=0.333
#   results_albedo_gA_tau02    -- tau=0.2
#   results_albedo_gA_tau01    -- tau=0.1  (requires build, n_per_triangle=2000)
#   results_albedo_gA_centroid -- centroid limit (tau -> 0)

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

# ------------------------------------------------------------------------------
# 1. Load setup objects
# ------------------------------------------------------------------------------

d_shared <- goebel2026::setup_shared
d_albedo <- goebel2026::setup_albedo
d_gA     <- goebel2026::setup_g_A

fine_grid_buffered <- d_shared$fine_grid_buffered
W_queen            <- d_shared$W_queen
X_obs_water        <- d_shared$X_obs_water
X_latent_water     <- d_shared$X_latent_water

y_alb         <- d_albedo$y
y_latent_true <- d_albedo$y_latent_true

soundings_proj <- spatintegrate::ensure_projected(goebel2026::soundings_augmented)
target_proj    <- spatintegrate::ensure_projected(goebel2026::target_grid)

target_idx            <- which(fine_grid_buffered$n_intersects > 0)
p                     <- ncol(d_shared$A_flat)
q                     <- ncol(X_obs_water)
lambda_beta           <- 0.01
n_soundings_per_pixel <- as.integer(Matrix::colSums(d_shared$A_flat > 0))

# ------------------------------------------------------------------------------
# 2. g-A builder (same as run_05b)
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
                      tau = 1/3, n_per_triangle = 500L, verbose = TRUE) {
  n <- nrow(soundings_sf); m <- nrow(fine_grid_sf)
  if (verbose) message(sprintf(
    "Building g-weighted A: tau=%.4f, n_per_triangle=%d", tau, n_per_triangle))
  sounding_geom <- sf::st_geometry(soundings_sf)
  fine_geom     <- sf::st_geometry(fine_grid_sf)
  crs           <- sf::st_crs(soundings_sf)
  touches       <- sf::st_intersects(soundings_sf, fine_grid_sf, sparse = TRUE)
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
    if (length(js) == 0L) return(list(j = integer(0), x = numeric(0)))
    inter_geoms <- suppressWarnings(sf::st_intersection(sounding_geom[i], fine_geom[js]))
    areas <- as.numeric(sf::st_area(inter_geoms))
    keep  <- !sf::st_is_empty(inter_geoms) & areas > 0
    if (!any(keep)) return(list(j = integer(0), x = numeric(0)))
    inter_sf   <- sf::st_sf(geometry = inter_geoms[keep], crs = crs)
    areas_keep <- areas[keep]; js_keep <- js[keep]
    g_fn    <- make_affine_g_fn(soundings_sf[i, ], tau = tau)
    g_basis <- function(coords) matrix(g_fn(coords), ncol = 1L)
    g_ij <- as.numeric(integrate_basis(
      basis_fn = g_basis, polygons_sf = inter_sf, n_per_triangle = n_per_triangle))
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
      soundings_sf = soundings_sf, sounding_geom = sounding_geom,
      fine_geom = fine_geom, crs = crs, touches = touches,
      tau = tau, n_per_triangle = n_per_triangle,
      make_affine_g_fn = make_affine_g_fn,
      integrate_basis = integrate_basis,
      .integrate_one_polygon = .integrate_one_polygon,
      .extract_polygon_pieces = .extract_polygon_pieces,
      .get_triangle_coords = .get_triangle_coords,
      .tri_area_shoelace = .tri_area_shoelace,
      .infer_k = .infer_k,
      map_unit_square_to_triangle = map_unit_square_to_triangle,
      generate_qmc_unit_square = generate_qmc_unit_square
    )
  )
  i_idx <- rep(seq_len(n), vapply(rows, function(r) length(r$j), integer(1L)))
  j_idx <- unlist(lapply(rows, `[[`, "j"))
  x_val <- unlist(lapply(rows, `[[`, "x"))
  as(Matrix::sparseMatrix(i = i_idx, j = j_idx, x = x_val,
                          dims = c(n, m)), "dgCMatrix")
}

# ------------------------------------------------------------------------------
# 3. Shared model components
# ------------------------------------------------------------------------------

make_Q_fun_rho <- function(W) {
  function(theta) {
    rho <- theta[["rho"]]
    S   <- Matrix::Diagonal(nrow(W)) - rho * W
    Q   <- Matrix::drop0(Matrix::forceSymmetric(Matrix::crossprod(S)))
    list(Q = Q, log_det_Q = NULL)
  }
}

make_Q_fun_aug <- function(Q_fun_spatial, q, lambda) {
  force(Q_fun_spatial); force(q); force(lambda)
  function(theta) {
    sp  <- Q_fun_spatial(theta)
    Q_a <- Matrix::drop0(Matrix::forceSymmetric(
      Matrix::bdiag(sp$Q, lambda * Matrix::Diagonal(q))
    ))
    list(Q = Q_a, log_det_Q = NULL)
  }
}

Q_fun_rho_aug <- make_Q_fun_aug(make_Q_fun_rho(W_queen), q, lambda_beta)

make_rsr_constraint <- function(X_obs, p, q) {
  force(X_obs); force(p); force(q)
  function(train_idx, A_aug_train) {
    A_train   <- A_aug_train[, seq_len(p), drop = FALSE]
    C_spatial <- as.matrix(t(X_obs[train_idx, , drop = FALSE]) %*% A_train)
    cbind(C_spatial, matrix(0, nrow = q, ncol = q))
  }
}

rsr_constraint <- make_rsr_constraint(X_obs_water, p, q)

# ------------------------------------------------------------------------------
# 4. Helpers
# ------------------------------------------------------------------------------

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
    message(sprintf("  skipping %s (already saved in data/)", obj_name)); return(TRUE)
  }
  FALSE
}

fit_albedo_gA <- function(A, tau_val, obj_name, A_label = "g_weighted") {
  A_aug <- as(cbind(A, X_obs_water), "dgCMatrix")
  C_aug_full <- cbind(as.matrix(t(X_obs_water) %*% A),
                      matrix(0, nrow = q, ncol = q))

  tuned <- fastblm::tune_cv(
    y          = y_alb,
    A          = A_aug,
    Q_fun      = Q_fun_rho_aug,
    R_inv      = NULL,
    theta_init = c(rho = 0.9),
    lower      = c(rho = 0.5),
    upper      = c(rho = 0.999),
    k          = 10L,
    constraint = rsr_constraint,
    seed       = 2026L,
    parallel   = TRUE,
    verbose    = TRUE
  )

  rho_hat <- tuned$theta[["rho"]]
  message(sprintf("  tuning done: rho=%.4f  phi=%.4f  sigma2e=%.4g",
                  rho_hat, tuned$phi, tuned$sigma2e))

  fit <- fastblm::fit_fastblm(
    y = y_alb, A = A_aug, Q = tuned$Q,
    phi = tuned$phi, solver = "cholesky"
  )
  fit <- fastblm::constrain(fit, C_aug_full)

  r_hat    <- fit$posterior_mean[seq_len(p)]
  beta_hat <- fit$posterior_mean[p + seq_len(q)]
  mu       <- r_hat + as.numeric(X_latent_water %*% beta_hat)
  A_pred   <- as(cbind(Matrix::Diagonal(p), X_latent_water), "dgCMatrix")
  se       <- fastblm::posterior_se(fit, A_new = A_pred, n_probes = 200L)
  se_beta  <- fastblm::posterior_se(fit, n_probes = 200L)[p + seq_len(q)]

  mu_t  <- mu[target_idx]
  se_t  <- se[target_idx]
  tr_t  <- y_latent_true[target_idx]
  ns_t  <- n_soundings_per_pixel[target_idx]
  ci_lo <- mu_t - 1.96 * se_t
  ci_hi <- mu_t + 1.96 * se_t

  coverage <- function(idx) {
    if (length(idx) == 0L) return(NA_real_)
    mean(tr_t[idx] >= ci_lo[idx] & tr_t[idx] <= ci_hi[idx], na.rm = TRUE)
  }

  resid  <- mu_t - tr_t
  rmse   <- sqrt(mean(resid^2, na.rm = TRUE))
  r2     <- 1 - sum(resid^2, na.rm = TRUE) /
    sum((tr_t - mean(tr_t, na.rm = TRUE))^2, na.rm = TRUE)

  message(sprintf("  RMSE=%.4f  R2=%.4f  beta_water=%.4f",
                  rmse, r2, beta_hat[2]))

  list(
    run_name              = obj_name,
    tags                  = list(tuning = "cv", response = "albedo",
                                 covariates = "water", constraint = "RSR",
                                 W = "queen", rho = rho_hat,
                                 A = A_label, tau = tau_val),
    timestamp             = Sys.time(),
    posterior_mean        = mu_t,
    posterior_se          = se_t,
    ci_lower              = ci_lo,
    ci_upper              = ci_hi,
    beta_hat              = beta_hat,
    se_beta               = se_beta,
    sigma2e               = fit$sigma2e,
    phi                   = tuned$phi,
    rho_opt               = rho_hat,
    cv_curve              = tuned$history,
    rmse                  = rmse,
    r2                    = r2,
    coverage_95_all       = coverage(seq_along(mu_t)),
    coverage_95_obs       = coverage(which(ns_t >= 1L)),
    coverage_95_dense     = coverage(which(ns_t >= 20L)),
    n_soundings_per_pixel = ns_t
  )
}

# ------------------------------------------------------------------------------
# 5. Runs using pre-built A matrices from setup_g_A
# ------------------------------------------------------------------------------

if (!.should_skip("results_albedo_gA_tau05")) {
  message("\n== 1. Albedo g-weighted A, tau=0.5 ==")
  results_albedo_gA_tau05 <- fit_albedo_gA(d_gA$A_g$tau_0.5, 0.5, "results_albedo_gA_tau05")
  usethis::use_data(results_albedo_gA_tau05, overwrite = TRUE)
}

if (!.should_skip("results_albedo_gA_tau033")) {
  message("\n== 2. Albedo g-weighted A, tau=0.333 ==")
  results_albedo_gA_tau033 <- fit_albedo_gA(d_gA$A_g$tau_0.333, 0.333, "results_albedo_gA_tau033")
  usethis::use_data(results_albedo_gA_tau033, overwrite = TRUE)
}

if (!.should_skip("results_albedo_gA_tau02")) {
  message("\n== 3. Albedo g-weighted A, tau=0.2 ==")
  results_albedo_gA_tau02 <- fit_albedo_gA(d_gA$A_g$tau_0.2, 0.2, "results_albedo_gA_tau02")
  usethis::use_data(results_albedo_gA_tau02, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 6. tau=0.1 -- build inline with n_per_triangle=2000
# ------------------------------------------------------------------------------

if (!.should_skip("results_albedo_gA_tau01")) {
  message("\n== 4. Albedo g-weighted A, tau=0.1 (building inline) ==")
  t_build <- system.time({
    A_tau01 <- build_g_A(soundings_proj, target_proj,
                         tau = 0.1, n_per_triangle = 2000L)
  })
  message(sprintf("  built in %.1f min", t_build["elapsed"] / 60))
  results_albedo_gA_tau01 <- fit_albedo_gA(A_tau01, 0.1, "results_albedo_gA_tau01")
  usethis::use_data(results_albedo_gA_tau01, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 7. Centroid limit
# ------------------------------------------------------------------------------

if (!.should_skip("results_albedo_gA_centroid")) {
  message("\n== 5. Albedo centroid limit (tau -> 0) ==")
  sounding_centroids <- sf::st_centroid(soundings_proj)
  nearest_pixel      <- sf::st_nearest_feature(sounding_centroids, target_proj)
  A_centroid <- as(Matrix::sparseMatrix(
    i = seq_len(nrow(soundings_proj)), j = nearest_pixel, x = 1,
    dims = c(nrow(soundings_proj), nrow(target_proj))
  ), "dgCMatrix")
  message(sprintf("  row sums [%.4f, %.4f]",
                  min(Matrix::rowSums(A_centroid)),
                  max(Matrix::rowSums(A_centroid))))
  results_albedo_gA_centroid <- fit_albedo_gA(
    A_centroid, 0.0, "results_albedo_gA_centroid",
    A_label = "centroid_assignment"
  )
  usethis::use_data(results_albedo_gA_centroid, overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# 8. Summary table
# ------------------------------------------------------------------------------

message("\n=== Albedo summary: RMSE/R2 by A specification ===")

all_runs <- list(
  list(name = "Uniform (canonical)", obj = "results_water_rsr_rho_cv"),
  list(name = "g tau=0.500",         obj = "results_albedo_gA_tau05"),
  list(name = "g tau=0.333",         obj = "results_albedo_gA_tau033"),
  list(name = "g tau=0.200",         obj = "results_albedo_gA_tau02"),
  list(name = "g tau=0.100",         obj = "results_albedo_gA_tau01"),
  list(name = "Centroid (tau->0)",   obj = "results_albedo_gA_centroid")
)

for (run in all_runs) {
  r <- tryCatch({
    e <- new.env()
    data(list = run$obj, package = "goebel2026", envir = e); e[[run$obj]]
  }, error = function(e) NULL)
  if (!is.null(r))
    message(sprintf("  %-22s  rho=%.3f  phi=%6.2f  RMSE=%.4f  R2=%.4f  coverage=%.3f",
                    run$name, r$rho_opt, r$phi,
                    r$rmse, r$r2, r$coverage_95_obs))
}

# ------------------------------------------------------------------------------

future::plan(future::sequential)
message("\nrun_05c_gA_albedo.R complete. Objects saved to data/:")
message("  results_albedo_gA_tau05")
message("  results_albedo_gA_tau033")
message("  results_albedo_gA_tau02")
message("  results_albedo_gA_tau01")
message("  results_albedo_gA_centroid")
