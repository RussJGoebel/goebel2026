# run_spatial_residual_diagnostic.R
#
# Diagnostic: assess spatial dependence of standardized residuals from a
# full-data fastblm fit.
#
# The key quantity is the standardized prediction residual at each observed
# sounding:
#
#   r_i = (y_i - A_i mu) / sqrt(sigma2e * (A_i K^{-1} A_i' + 1))
#
# The denominator accounts for both posterior uncertainty (A_i K^{-1} A_i')
# and observation noise (+ 1), so soundings with few overlapping neighbors
# are not artificially flagged. Under a correctly specified model with no
# residual spatial correlation, r_i should be approximately i.i.d. N(0,1)
# with no spatial structure.
#
# USAGE
# -----
# Source your model setup, then call:
#
#   result <- spatial_residual_diagnostic(fit, y, A, sounding_coords)
#
# where sounding_coords is an n x 2 matrix or data.frame with columns x, y
# (projected coordinates, not lat/lon).
#
# The function returns a list with the standardized residuals and variogram,
# and saves a multi-panel PDF to disk.
#
# REQUIREMENTS
# ------------
#   fit             : fastblm_fit object from fit_fastblm(), solver = "cholesky"
#                     or "woodbury" (PCG path gives only stochastic variances).
#   y               : observed sounding values (length n)
#   A               : n x p observation matrix (sparse OK)
#   sounding_coords : n x 2 matrix/data.frame, columns named x and y,
#                     in a projected CRS (metres preferred).
#   gstat           : install.packages("gstat") -- for empirical variogram
#   sf              : install.packages("sf")    -- for spatial plotting
#
# OUTPUT FILES
# ------------
#   spatial_residual_diagnostic.pdf   -- four-panel diagnostic figure

# ─────────────────────────────────────────────────────────────────────────────
# 0.  Dependencies
# ─────────────────────────────────────────────────────────────────────────────
required_pkgs <- c("Matrix", "gstat", "sf", "ggplot2", "patchwork", "dplyr")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace,
                                       quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Please install missing packages: ",
       paste(missing_pkgs, collapse = ", "))
}

library(Matrix)
library(gstat)
library(sf)
library(ggplot2)
library(patchwork)
library(dplyr)

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Main diagnostic function
# ─────────────────────────────────────────────────────────────────────────────

#' Spatial residual diagnostic for a fastblm_fit
#'
#' @param fit          fastblm_fit object (cholesky or woodbury solver)
#' @param y            numeric observed response vector (length n)
#' @param A            n x p observation matrix
#' @param coords       n x 2 data.frame/matrix with columns x, y (projected)
#' @param n_lags       number of variogram distance lags (default 15)
#' @param max_dist     maximum variogram distance. NULL = half the bounding box
#'                     diagonal (a sensible default for most domains).
#' @param output_pdf   file path for output PDF
#' @param crs          EPSG code for coords (used only for the map; default 32619
#'                     = UTM zone 19N, covering Boston). Set to NULL to skip CRS.
#'
#' @return invisible list with fields:
#'   $std_resid   : numeric vector of standardized residuals
#'   $pred_sd     : posterior predictive SD at each sounding
#'   $variogram   : gstat variogram object
spatial_residual_diagnostic <- function(fit,
                                        y,
                                        A,
                                        coords,
                                        n_lags     = 15L,
                                        max_dist   = NULL,
                                        output_pdf = "spatial_residual_diagnostic.pdf",
                                        crs        = 32619L) {

  stopifnot(inherits(fit, "fastblm_fit"))
  stopifnot(fit$solver_type %in% c("cholesky", "woodbury"),
            "PCG solver does not store an exact Cholesky factor; use
            cholesky or woodbury for exact posterior variances.")

  n <- length(y)
  stopifnot(nrow(A) == n, nrow(coords) == n)

  coords <- as.data.frame(coords)
  if (!all(c("x", "y") %in% names(coords)))
    stop("coords must have columns named 'x' and 'y'.")

  # ── 1a. Fitted values ──────────────────────────────────────────────────────
  message("Computing fitted values ...")
  fitted_vals <- as.numeric(A %*% fit$posterior_mean)
  raw_resid   <- y - fitted_vals

  # ── 1b. Posterior predictive variance at each sounding ────────────────────
  # Var(y_i | y) = sigma2e * (A_i K^{-1} A_i'  +  1)
  #              = posterior variance of A_i gamma  +  sigma2e
  #
  # posterior_se() with A_new = A gives sqrt(sigma2e * A K^{-1} A') for each row
  message("Computing posterior predictive SDs (this may take a moment) ...")
  post_sd_gamma <- posterior_se(fit, A_new = A)   # length n
  pred_var      <- post_sd_gamma^2 + fit$sigma2e  # add observation noise
  pred_sd       <- sqrt(pred_var)

  # ── 1c. Standardized residuals ────────────────────────────────────────────
  std_resid <- raw_resid / pred_sd

  message(sprintf(
    "Standardized residuals: mean = %.3f, sd = %.3f (expect ~0 and ~1)",
    mean(std_resid), sd(std_resid)
  ))

  # ── 1d. Build sf object for spatial operations ────────────────────────────
  df <- data.frame(
    x          = coords$x,
    y          = coords$y,
    raw_resid  = raw_resid,
    std_resid  = std_resid,
    pred_sd    = pred_sd,
    n_soundings_approx = rowSums(A != 0)  # pixels each sounding touches
  )

  pts_sf <- sf::st_as_sf(df, coords = c("x", "y"),
                         crs = if (!is.null(crs)) crs else NA)

  # ── 1e. Empirical variogram of standardized residuals ────────────────────
  message("Computing empirical variogram ...")

  if (is.null(max_dist)) {
    bbox     <- sf::st_bbox(pts_sf)
    max_dist <- 0.5 * sqrt((bbox["xmax"] - bbox["xmin"])^2 +
                             (bbox["ymax"] - bbox["ymin"])^2)
  }

  sp_df <- as(pts_sf, "Spatial")   # gstat still prefers sp objects
  vgm_emp <- gstat::variogram(std_resid ~ 1,
                              data    = sp_df,
                              cutoff  = max_dist,
                              width   = max_dist / n_lags)

  # Fit a Matern/exponential model as a visual reference line
  # (we are NOT using this for inference, just as a guide to the eye)
  vgm_fit <- tryCatch(
    gstat::fit.variogram(vgm_emp,
                         gstat::vgm(psill = var(std_resid),
                                    model = "Exp",
                                    range = max_dist / 3,
                                    nugget = 0)),
    error = function(e) NULL
  )

  # ── 1f. Moran's I as a single summary statistic ───────────────────────────
  moran_result <- tryCatch({
    # inverse-distance weights, capped at 5 nearest neighbours
    coords_mat <- as.matrix(coords[, c("x", "y")])
    dmat       <- as.matrix(dist(coords_mat))
    diag(dmat) <- Inf
    k          <- 5L
    W          <- matrix(0, n, n)
    for (i in seq_len(n)) {
      nn_idx      <- order(dmat[i, ])[seq_len(k)]
      W[i, nn_idx] <- 1 / dmat[i, nn_idx]
    }
    W <- W / rowSums(W)   # row-standardise

    # Moran's I
    z    <- std_resid - mean(std_resid)
    I    <- (n / sum(W)) * as.numeric(z %*% W %*% z) / sum(z^2)
    E_I  <- -1 / (n - 1)
    list(I = I, E_I = E_I)
  }, error = function(e) {
    message("Moran's I skipped (", conditionMessage(e), ")")
    NULL
  })

  if (!is.null(moran_result)) {
    message(sprintf("Moran's I = %.4f  (expected under no autocorrelation: %.4f)",
                    moran_result$I, moran_result$E_I))
  }

  # ─────────────────────────────────────────────────────────────────────────
  # 2.  Plots
  # ─────────────────────────────────────────────────────────────────────────
  message("Generating plots ...")

  clamp <- function(x, lo = -3.5, hi = 3.5) pmax(pmin(x, hi), lo)

  # ── Panel A: spatial map of standardized residuals ────────────────────────
  p_map <- ggplot(df, aes(x = x, y = y, colour = clamp(std_resid))) +
    geom_point(size = 1.2, alpha = 0.8) +
    scale_colour_distiller(palette  = "RdBu",
                           limits   = c(-3.5, 3.5),
                           name     = "Std. residual\n(clamped ±3.5)") +
    coord_equal() +
    labs(title = "A: Standardized residuals (spatial)",
         x = "Easting", y = "Northing") +
    theme_bw(base_size = 10)

  # ── Panel B: residuals vs posterior predictive SD (information content) ───
  p_sd <- ggplot(df, aes(x = pred_sd, y = std_resid)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_hline(yintercept = c(-2, 0, 2), linetype = c("dashed","solid","dashed"),
               colour = "steelblue") +
    geom_smooth(method = "loess", se = TRUE, colour = "tomato", linewidth = 0.8) +
    labs(title = "B: Std. residual vs. predictive SD",
         subtitle = "High SD = poorly constrained soundings",
         x = "Posterior predictive SD", y = "Standardized residual") +
    theme_bw(base_size = 10)

  # ── Panel C: empirical variogram of standardized residuals ────────────────
  p_vgm <- ggplot(vgm_emp, aes(x = dist, y = gamma)) +
    geom_point(aes(size = np), colour = "grey30") +
    scale_size_continuous(name = "Pairs", range = c(1, 4)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "steelblue") +
    labs(title = "C: Empirical variogram of std. residuals",
         subtitle = "Flat at sill ≈ 1 → no residual spatial autocorrelation",
         x = "Distance", y = "Semivariance") +
    theme_bw(base_size = 10)

  if (!is.null(vgm_fit)) {
    vgm_line <- gstat::variogramLine(vgm_fit,
                                     maxdist = max_dist,
                                     n       = 200L)
    p_vgm <- p_vgm +
      geom_line(data = vgm_line, aes(x = dist, y = gamma),
                colour = "tomato", linewidth = 0.8, inherit.aes = FALSE)
  }

  # ── Panel D: Q-Q plot of standardized residuals ───────────────────────────
  p_qq <- ggplot(df, aes(sample = std_resid)) +
    stat_qq(alpha = 0.4, size = 0.8) +
    stat_qq_line(colour = "tomato", linewidth = 0.8) +
    labs(title = "D: Q-Q plot of std. residuals",
         subtitle = "Should follow N(0,1) if model is well-specified",
         x = "Theoretical quantiles", y = "Sample quantiles") +
    theme_bw(base_size = 10)

  # ── Assemble and save ─────────────────────────────────────────────────────
  combined <- (p_map | p_sd) / (p_vgm | p_qq) +
    patchwork::plot_annotation(
      title    = "Spatial residual diagnostic",
      subtitle = sprintf(
        "n = %d soundings  |  mean(r) = %.3f  |  sd(r) = %.3f%s",
        n, mean(std_resid), sd(std_resid),
        if (!is.null(moran_result))
          sprintf("  |  Moran's I = %.4f", moran_result$I)
        else ""
      )
    )

  ggsave(output_pdf, combined, width = 12, height = 9, device = "pdf")
  message("Saved: ", output_pdf)

  invisible(list(
    std_resid  = std_resid,
    raw_resid  = raw_resid,
    pred_sd    = pred_sd,
    variogram  = vgm_emp,
    vgm_fit    = vgm_fit,
    moran      = moran_result,
    data       = df
  ))
}


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Example usage (replace with your actual objects)
# ─────────────────────────────────────────────────────────────────────────────
if (FALSE) {

  # Source your model fitting infrastructure
  # source("fit.R"); source("posterior_se.R"); source("utils.R"); etc.

  # Fit the full model at the CV-selected hyperparameters
  fit <- fit_fastblm(y, A, Q, phi = tuned$phi, solver = "cholesky")

  # sounding_coords: n x 2 data.frame with columns x, y in metres
  # e.g. UTM easting/northing of sounding centroids
  result <- spatial_residual_diagnostic(
    fit        = fit,
    y          = y,
    A          = A,
    coords     = sounding_coords,   # your n x 2 coordinate data.frame
    n_lags     = 15L,
    max_dist   = NULL,              # auto: half bounding-box diagonal
    output_pdf = "spatial_residual_diagnostic.pdf",
    crs        = 32619L             # UTM zone 19N (Boston); adjust as needed
  )

  # Inspect
  hist(result$std_resid, breaks = 40, main = "Standardized residuals")
  print(result$variogram)
}
