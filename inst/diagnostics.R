## =============================================================================
## diagnostics.R
## Diagnostic functions for fastblm ablation experiments.
##
## Main entry point:
##   diag <- blm_diagnostics(fit, se, truth, grid_sf, label, ...)
##   plot(diag)                    # spatial maps panel
##   plot_cv_profile(tuned, ...)   # phi profile curve from tune_cv history
##   plot_fold_errors(soundings_sf, folds, fold_rmse, ...)  # spatial CV errors
##   constraint_check(fit, C)      # RSR sanity check
## =============================================================================

library(sf)
library(ggplot2)
library(patchwork)

# =============================================================================
# 1.  blm_diagnostics()  --  collect all scalars + spatial fields
# =============================================================================

#' Collect diagnostics for a fastblm fit
#'
#' @param fit       fastblm_fit object. May be an augmented fit over p+q
#'                  coefficients (spatial field + fixed effects). If so,
#'                  supply \code{X_grid} and \code{p_spatial} to get
#'                  diagnostics on the predicted field rather than r alone.
#' @param se        numeric vector of posterior SEs. For augmented fits this
#'                  should be the SE of the predicted field
#'                  \eqn{\hat{y}_i = r_i + x_{\text{grid},i}^\top \beta},
#'                  computed via \code{posterior_se(fit, A_new = A_pred)}
#'                  where \code{A_pred = [I_p | X_grid]}.
#' @param truth     numeric vector of true pixel values (length p, NA allowed)
#' @param grid_sf   sf object with p rows representing the target grid
#' @param label     character label for this model (used in plot titles)
#' @param A         n x p spatial overlap matrix. If supplied, coverage, CI
#'                  width, and outlier rate are restricted to pixels with at
#'                  least one overlapping sounding. RMSE and R^2 use all
#'                  non-NA pixels.
#' @param C_check   optional constraint matrix C (q x p) -- if supplied,
#'                  ||C mu||_inf is reported as a constraint residual. Should
#'                  use the spatial block of the posterior mean only.
#' @param X_grid    optional p x q matrix of fixed-effect covariates at the
#'                  grid level. If supplied together with \code{p_spatial},
#'                  the predicted field \eqn{r + X_{\rm grid} \beta} is used
#'                  for all diagnostics instead of the raw posterior mean.
#' @param p_spatial optional integer: number of spatial coefficients (p) when
#'                  fit is over an augmented p+q system. If NULL and
#'                  \code{X_grid} is NULL, assumes fit is spatial-only.
#'
#' @return list of class "blm_diagnostics"
blm_diagnostics <- function(fit, se, truth, grid_sf, label,
                            A = NULL, C_check = NULL,
                            X_grid = NULL, p_spatial = NULL) {

  stopifnot(inherits(fit, "fastblm_fit"))

  # -------------------------------------------------------------------------
  # Resolve the predicted field mu_pred (length p) and posterior SEs.
  #
  # Three cases:
  #   (a) Spatial-only fit: fit$posterior_mean is length p. mu_pred = mu.
  #   (b) Augmented fit, X_grid supplied: fit$posterior_mean is length p+q.
  #       mu_pred = r_hat + X_grid %*% beta_hat = A_pred %*% mu_aug.
  #       se should already be posterior_se(fit, A_new=A_pred) -- p-vector.
  #   (c) Augmented fit, X_grid not supplied: caller has pre-split r_hat
  #       and passed it as fit$posterior_mean (length p). Treat as (a).
  # -------------------------------------------------------------------------

  if (!is.null(X_grid) && !is.null(p_spatial)) {
    # Case (b): augmented fit, compute predicted field explicitly
    mu_aug  <- fit$posterior_mean            # length p+q
    q       <- length(mu_aug) - p_spatial
    r_hat   <- mu_aug[seq_len(p_spatial)]
    beta_hat <- mu_aug[p_spatial + seq_len(q)]
    mu      <- r_hat + as.numeric(X_grid %*% beta_hat)   # length p
    p       <- p_spatial
  } else {
    # Case (a) or (c): posterior_mean is already the p-vector we want
    mu <- fit$posterior_mean
    p  <- length(mu)
  }

  stopifnot(length(mu) == length(truth),
            length(se)  == length(truth))

  na_idx <- is.na(truth)

  # --- observed pixel mask --------------------------------------------------
  if (!is.null(A)) {
    obs_idx <- as.logical(Matrix::colSums(A > 0) > 0)
  } else {
    obs_idx <- rep(TRUE, p)
  }
  eval_idx <- !na_idx & obs_idx

  if (sum(eval_idx) == 0L)
    warning("No pixels pass both the NA and observation-overlap filters.")

  # --- scalar summaries -----------------------------------------------------
  resid <- mu - truth

  rmse    <- sqrt(mean(resid[!na_idx]^2))
  lm_fit  <- lm(truth[!na_idx] ~ mu[!na_idx])
  r2      <- summary(lm_fit)$r.squared
  lm_coef <- coef(lm_fit)

  lower    <- mu - 1.96 * se
  upper    <- mu + 1.96 * se
  coverage <- mean(truth[eval_idx] >= lower[eval_idx] &
                     truth[eval_idx] <= upper[eval_idx])
  ci_width     <- median((upper - lower)[eval_idx])
  outlier_rate <- mean(abs(resid[eval_idx]) > 2 * se[eval_idx])

  scalars <- list(
    label        = label,
    phi          = fit$phi,
    sigma2e      = fit$sigma2e,
    sigma2b      = fit$sigma2b,
    rmse         = rmse,
    r2           = r2,
    coverage     = coverage,
    ci_width     = ci_width,
    outlier_rate = outlier_rate,
    lm_intercept = lm_coef[[1]],
    lm_slope     = lm_coef[[2]]
  )

  # --- constraint residual --------------------------------------------------
  # C_check acts on the spatial field only
  r_for_constraint <- if (!is.null(X_grid) && !is.null(p_spatial))
    fit$posterior_mean[seq_len(p_spatial)]
  else
    mu

  constraint_linf <- if (!is.null(C_check)) {
    max(abs(as.numeric(C_check %*% r_for_constraint)))
  } else {
    NA_real_
  }
  scalars$constraint_linf <- constraint_linf
  scalars$n_obs_pixels    <- sum(obs_idx)
  scalars$n_eval_pixels   <- sum(eval_idx)

  # --- spatial fields -------------------------------------------------------
  in_ci_vec              <- truth >= lower & truth <= upper
  in_ci_vec[!obs_idx]    <- NA
  ci_width_vec           <- upper - lower
  ci_width_vec[!obs_idx] <- NA

  grid_sf$posterior_mean <- mu      # predicted field, not raw mu_aug
  grid_sf$posterior_se   <- se
  grid_sf$error          <- resid
  grid_sf$abs_error      <- abs(resid)
  grid_sf$in_ci          <- in_ci_vec
  grid_sf$ci_width       <- ci_width_vec
  grid_sf$truth          <- truth
  grid_sf$observed       <- obs_idx

  structure(
    list(
      label    = label,
      scalars  = scalars,
      grid     = grid_sf,
      fit      = fit,
      se       = se,
      truth    = truth,
      obs_idx  = obs_idx
    ),
    class = "blm_diagnostics"
  )
}

# =============================================================================
# 2.  print method
# =============================================================================

#' @export
print.blm_diagnostics <- function(x, ...) {
  s <- x$scalars
  cat(sprintf("\n=== %s ===\n", s$label))
  cat(sprintf("  phi         : %.4f\n", s$phi))
  cat(sprintf("  sigma2e     : %.4f\n", s$sigma2e))
  cat(sprintf("  sigma2b     : %.4f\n", s$sigma2b))
  cat(sprintf("  RMSE        : %.4f\n", s$rmse))
  cat(sprintf("  R^2         : %.4f   (lm slope=%.3f, intercept=%.3f)\n",
              s$r2, s$lm_slope, s$lm_intercept))
  cat(sprintf("  Coverage    : %.4f   (target 0.95, n=%d observed pixels)\n",
              s$coverage, s$n_eval_pixels))
  cat(sprintf("  Median CI width: %.4f\n", s$ci_width))
  cat(sprintf("  Outlier rate   : %.4f   (|err|>2*se, observed pixels)\n",
              s$outlier_rate))
  if (!is.na(s$constraint_linf))
    cat(sprintf("  ||C mu||_inf   : %.2e   (constraint residual)\n",
                s$constraint_linf))
  invisible(x)
}

# =============================================================================
# 3.  plot.blm_diagnostics()  --  four-panel spatial map
# =============================================================================

#' @export
plot.blm_diagnostics <- function(x, ...) {

  g  <- x$grid
  lb <- x$label

  th <- theme_minimal(base_size = 10) +
    theme(
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      panel.grid       = element_blank(),
      plot.title       = element_text(size = 9, face = "bold"),
      legend.key.width = unit(0.4, "cm"),
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7)
    )

  p_mean <- ggplot(g) +
    geom_sf(aes(fill = posterior_mean), colour = NA, size = 0) +
    scale_fill_viridis_c(option = "magma", name = "Post.\nmean",
                         na.value = "grey85") +
    labs(title = "Posterior mean") + th

  p_se <- ggplot(g) +
    geom_sf(aes(fill = posterior_se), colour = NA, size = 0) +
    scale_fill_viridis_c(option = "inferno", name = "Post. SE",
                         na.value = "grey85") +
    labs(title = "Posterior SE") + th

  err_lim <- max(abs(g$error), na.rm = TRUE)
  p_err <- ggplot(g) +
    geom_sf(aes(fill = error), colour = NA, size = 0) +
    scale_fill_distiller(palette = "RdBu", limits = c(-err_lim, err_lim),
                         name = "Error", na.value = "grey85") +
    labs(title = "Error (mean \u2212 truth)") + th

  p_cov <- ggplot(g) +
    geom_sf(aes(fill = in_ci), colour = NA, size = 0) +
    scale_fill_manual(values = c("TRUE" = "#2166ac", "FALSE" = "#d6604d"),
                      name = "In 95% CI", na.value = "grey85",
                      labels = c("TRUE" = "Yes", "FALSE" = "No")) +
    labs(title = sprintf("95%% CI coverage (%.1f%%)",
                         100 * x$scalars$coverage)) + th

  (p_mean | p_se) / (p_err | p_cov) +
    patchwork::plot_annotation(title = lb)
}

# =============================================================================
# 4.  scatter_truth_vs_mean()
# =============================================================================

scatter_truth_vs_mean <- function(diag, colour_var = NULL,
                                  colour_label = "Covariate") {
  mu    <- diag$grid$posterior_mean   # already the predicted field
  truth <- diag$truth
  ok    <- !is.na(truth)

  df <- data.frame(truth = truth[ok], mu = mu[ok])
  if (!is.null(colour_var)) df$col <- colour_var[ok]

  lim <- range(c(df$truth, df$mu), na.rm = TRUE)

  p <- ggplot(df, aes(x = truth, y = mu)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "black") +
    labs(
      x        = "Truth",
      y        = "Posterior mean",
      title    = diag$label,
      subtitle = sprintf("RMSE=%.4f  R\u00b2=%.3f  Coverage=%.3f",
                         diag$scalars$rmse, diag$scalars$r2,
                         diag$scalars$coverage)
    ) +
    coord_fixed(xlim = lim, ylim = lim) +
    theme_minimal(base_size = 11)

  if (!is.null(colour_var)) {
    p <- p +
      geom_point(aes(colour = col), size = 0.8, alpha = 0.6) +
      scale_colour_viridis_c(option = "viridis", name = colour_label)
  } else {
    p <- p + geom_point(size = 0.8, alpha = 0.4, colour = "#333333")
  }
  p
}

# =============================================================================
# 5.  plot_cv_profile()
# =============================================================================

plot_cv_profile <- function(tuned, n_phi = NULL) {

  h <- tuned$history
  if (is.null(h) || nrow(h) == 0) {
    message("No history found in tuned object.")
    return(invisible(NULL))
  }

  .flatten <- function(x) {
    if (is.list(x)) vapply(x, function(v) as.numeric(v)[[1L]], numeric(1L))
    else as.numeric(x)
  }

  phi_vals <- .flatten(h$phi)

  if (!is.null(h$ll)) {
    score_vals  <- .flatten(h$ll)
    score_label <- "log-likelihood"
    best_idx    <- which.max(score_vals)
  } else if (!is.null(h$cv_score)) {
    score_vals  <- .flatten(h$cv_score)
    score_label <- "CV score (MSE)"
    best_idx    <- which.min(score_vals)
  } else {
    message("History has neither 'll' nor 'cv_score' column.")
    return(invisible(NULL))
  }

  theta_name <- if (length(tuned$theta) == 1L) names(tuned$theta)[[1L]] else NULL
  has_theta  <- !is.null(theta_name) && !is.null(h$theta)

  if (has_theta) {
    theta_vals <- .flatten(h$theta)
    df <- data.frame(theta = theta_vals, log_phi = log(phi_vals),
                     score = score_vals)
    p <- ggplot(df, aes(x = theta, y = score)) +
      geom_point(aes(colour = log_phi), size = 2) +
      geom_point(data = df[best_idx, , drop = FALSE],
                 shape = 4, size = 4, colour = "red", stroke = 1.5) +
      scale_colour_viridis_c(name = "log(\u03d5)") +
      labs(x = theta_name, y = score_label,
           title    = sprintf("Optimisation trace (%s)", tuned$method),
           subtitle = sprintf("Optimum: %s=%.4f  phi=%.4f",
                              theta_name, tuned$theta[[1L]], tuned$phi)) +
      theme_minimal(base_size = 11)
  } else {
    df <- data.frame(log_phi = log(phi_vals), score = score_vals)
    p <- ggplot(df, aes(x = log_phi, y = score)) +
      geom_point(size = 2, colour = "#3182bd") +
      geom_point(data = df[best_idx, , drop = FALSE],
                 shape = 4, size = 4, colour = "red", stroke = 1.5) +
      labs(x = "log(\u03d5)", y = score_label,
           title    = sprintf("Optimisation trace (%s)", tuned$method),
           subtitle = sprintf("Optimum: phi=%.4f", tuned$phi)) +
      theme_minimal(base_size = 11)
  }
  p
}

# =============================================================================
# 6.  plot_fold_errors()
# =============================================================================

plot_fold_errors <- function(soundings_sf, folds, y, A, fit) {

  k     <- max(folds)
  mu    <- fit$posterior_mean
  y_hat <- as.numeric(A %*% mu[seq_len(ncol(A))])   # spatial block only
  resid <- y - y_hat

  fold_rmse <- vapply(seq_len(k), function(f) {
    idx <- which(folds == f)
    sqrt(mean(resid[idx]^2))
  }, numeric(1))

  cents           <- sf::st_centroid(soundings_sf)
  cents$fold      <- folds
  cents$fold_rmse <- fold_rmse[folds]
  cents$residual  <- resid

  th_map <- theme_minimal(base_size = 10) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())

  p_folds <- ggplot(cents) +
    geom_sf(aes(colour = factor(fold)), size = 0.6, alpha = 0.7) +
    scale_colour_manual(values = scales::hue_pal()(k), name = "Fold",
                        guide = guide_legend(override.aes = list(size = 2))) +
    labs(title = "Fold assignments") + th_map

  p_rmse <- ggplot(cents) +
    geom_sf(aes(colour = fold_rmse), size = 0.6, alpha = 0.9) +
    scale_colour_viridis_c(option = "plasma", name = "Fold\nRMSE") +
    labs(title = sprintf("Per-fold RMSE  (mean=%.4f)", mean(fold_rmse))) +
    th_map

  rlim  <- max(abs(resid))
  p_res <- ggplot(cents) +
    geom_sf(aes(colour = residual), size = 0.6, alpha = 0.9) +
    scale_colour_distiller(palette = "RdBu", limits = c(-rlim, rlim),
                           name = "Residual") +
    labs(title = "Sounding-level residuals") + th_map

  (p_folds | p_rmse | p_res) +
    patchwork::plot_annotation(title = "CV fold diagnostics")
}

# =============================================================================
# 7.  compare_diagnostics()
# =============================================================================

compare_diagnostics <- function(diag_list) {
  rows <- lapply(diag_list, function(d) {
    s <- d$scalars
    data.frame(
      model           = s$label,
      phi             = s$phi,
      sigma2e         = s$sigma2e,
      sigma2b         = s$sigma2b,
      rmse            = s$rmse,
      r2              = s$r2,
      coverage        = s$coverage,
      ci_width        = s$ci_width,
      outlier_rate    = s$outlier_rate,
      constraint_linf = s$constraint_linf,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  df
}

# =============================================================================
# 8.  plot_compare_maps()
# =============================================================================

plot_compare_maps <- function(diag_list, field = "posterior_mean") {

  stopifnot(field %in% c("posterior_mean", "posterior_se",
                         "error", "abs_error", "in_ci"))

  if (field != "in_ci") {
    all_vals <- unlist(lapply(diag_list, function(d) {
      v <- sf::st_drop_geometry(d$grid)[[field]]
      v[!is.na(v)]
    }))
    vmin <- min(all_vals); vmax <- max(all_vals)
    if (field == "error") { lim <- max(abs(c(vmin,vmax))); vmin <- -lim; vmax <- lim }
  }

  th <- theme_minimal(base_size = 9) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(size = 8, face = "bold"))

  plots <- lapply(diag_list, function(d) {
    g <- d$grid
    if (field == "in_ci") {
      ggplot(g) +
        geom_sf(aes(fill = in_ci), colour = NA, size = 0) +
        scale_fill_manual(values = c("TRUE"="#2166ac","FALSE"="#d6604d"),
                          name = "In CI", na.value = "grey85") +
        labs(title = d$label) + th
    } else if (field == "error") {
      ggplot(g) +
        geom_sf(aes_string(fill = field), colour = NA, size = 0) +
        scale_fill_distiller(palette = "RdBu", limits = c(vmin,vmax),
                             name = field, na.value = "grey85") +
        labs(title = d$label) + th
    } else {
      ggplot(g) +
        geom_sf(aes_string(fill = field), colour = NA, size = 0) +
        scale_fill_viridis_c(option = "magma", limits = c(vmin,vmax),
                             name = field, na.value = "grey85") +
        labs(title = d$label) + th
    }
  })

  patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(title = sprintf("Comparison: %s", field))
}

# =============================================================================
# 9.  save_ablation()
# =============================================================================

#' Save a single model's diagnostics and register it in the manifest
#'
#' @param model      list with: fit, se, and optionally tuned, colour_var,
#'                   colour_label, C_check, X_grid, p_spatial.
#'                   For augmented fits (spatial + fixed effects), supply
#'                   \code{X_grid} (p x q grid-level covariates) and
#'                   \code{p_spatial} (number of spatial coefficients) so
#'                   that diagnostics use the predicted field
#'                   \eqn{r + X_{\rm grid} \beta} rather than the raw
#'                   augmented posterior mean.
#'                   \code{se} should be the SE of that predicted field,
#'                   computed via \code{posterior_se(fit, A_new = A_pred)}
#'                   where \code{A_pred = [I_p | X_grid]}.
#' @param run_name   unique character key for this run
#' @param tags       named list of free-form metadata
#' @param truth      numeric vector of true pixel values (length p)
#' @param grid_sf    sf object with p rows
#' @param A          n x p spatial overlap matrix
#' @param overwrite  if FALSE and run exists, skip
#' @param width,height,dpi  plot dimensions
#'
#' @return invisibly returns the blm_diagnostics object
#' @export
save_ablation <- function(model, run_name, tags = list(),
                          truth, grid_sf, A,
                          overwrite = TRUE,
                          width = 10, height = 7, dpi = 150) {

  stopifnot(is.character(run_name), nchar(run_name) > 0)
  stopifnot(is.list(model), !is.null(model$fit), !is.null(model$se))

  pkg_root <- rprojroot::find_package_root_file(path = ".")
  abl_dir  <- file.path(pkg_root, "inst", "data", "ablations")
  run_dir  <- file.path(abl_dir, run_name)
  rds_path <- file.path(run_dir, "diagnostics.rds")
  manifest_path <- file.path(abl_dir, "manifest.csv")

  dir.create(run_dir, showWarnings = FALSE, recursive = TRUE)

  if (!overwrite && file.exists(rds_path)) {
    message(sprintf("[%s] already exists, skipping (overwrite=FALSE)", run_name))
    return(invisible(readRDS(rds_path)))
  }

  .save_plot <- function(p, fname, w = width, h = height) {
    path <- file.path(run_dir, fname)
    ggplot2::ggsave(path, plot = p, width = w, height = h, dpi = dpi)
    message("  saved: ", path)
  }

  message(sprintf("\n[%s] building diagnostics...", run_name))

  d <- blm_diagnostics(
    fit       = model$fit,
    se        = model$se,
    truth     = truth,
    grid_sf   = grid_sf,
    label     = run_name,
    A         = A,
    C_check   = model$C_check   %||% NULL,
    X_grid    = model$X_grid    %||% NULL,
    p_spatial = model$p_spatial %||% NULL
  )
  print(d)

  saveRDS(d, rds_path)
  message("  saved: ", rds_path)

  .save_plot(plot(d), "maps.png")
  .save_plot(scatter_truth_vs_mean(
    d,
    colour_var   = model$colour_var   %||% NULL,
    colour_label = model$colour_label %||% "Covariate"
  ), "scatter.png", w = 6, h = 6)

  if (!is.null(model$tuned)) {
    p_prof <- plot_cv_profile(model$tuned)
    if (!is.null(p_prof)) .save_plot(p_prof, "profile.png", w = 7, h = 5)
  }

  obs_density <- as.numeric(Matrix::colSums(A > 0))
  p_sed <- ggplot2::ggplot(
    data.frame(obs_density = obs_density, posterior_se = model$se),
    ggplot2::aes(x = obs_density, y = posterior_se)
  ) +
    ggplot2::geom_hex(bins = 40) +
    ggplot2::scale_fill_viridis_c(option = "magma", name = "count") +
    ggplot2::geom_smooth(method = "loess", se = FALSE,
                         colour = "red", linewidth = 0.8) +
    ggplot2::labs(x = "Observation density", y = "Posterior SE",
                  title = sprintf("%s \u2014 SE vs density", run_name)) +
    ggplot2::theme_minimal(base_size = 11)
  .save_plot(p_sed, "se_density.png", w = 7, h = 5)

  # --- update manifest ------------------------------------------------------
  s       <- d$scalars
  new_row <- data.frame(
    run_name     = run_name,
    timestamp    = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    phi          = s$phi,
    sigma2e      = s$sigma2e,
    sigma2b      = s$sigma2b,
    rmse         = s$rmse,
    r2           = s$r2,
    coverage     = s$coverage,
    ci_width     = s$ci_width,
    outlier_rate = s$outlier_rate,
    n_obs_pixels = s$n_obs_pixels,
    stringsAsFactors = FALSE
  )
  for (tag_nm in names(tags)) new_row[[tag_nm]] <- as.character(tags[[tag_nm]])

  if (file.exists(manifest_path)) {
    manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
    for (col in names(new_row))
      if (!col %in% names(manifest)) manifest[[col]] <- NA_character_
    for (col in names(manifest))
      if (!col %in% names(new_row)) new_row[[col]] <- NA_character_
    manifest <- manifest[manifest$run_name != run_name, , drop = FALSE]
    manifest <- rbind(manifest, new_row[, names(manifest)])
  } else {
    manifest <- new_row
  }

  write.csv(manifest, manifest_path, row.names = FALSE)
  message("  manifest updated: ", manifest_path)
  invisible(d)
}

# =============================================================================
# 10. load_ablations()
# =============================================================================

#' @export
load_ablations <- function(filter = NULL, run_names = NULL) {

  pkg_root <- rprojroot::find_package_root_file(path = ".")
  abl_dir  <- file.path(pkg_root, "inst", "data", "ablations")
  manifest_path <- file.path(abl_dir, "manifest.csv")

  if (!file.exists(manifest_path))
    stop("No manifest.csv found in ", abl_dir, ". Run save_ablation() first.")

  manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

  if (!is.null(run_names)) {
    manifest <- manifest[manifest$run_name %in% run_names, , drop = FALSE]
  } else if (!is.null(filter)) {
    keep <- rep(TRUE, nrow(manifest))
    for (tag_nm in names(filter)) {
      if (!tag_nm %in% names(manifest)) {
        warning(sprintf("Tag '%s' not in manifest; ignoring.", tag_nm)); next
      }
      keep <- keep & (manifest[[tag_nm]] == as.character(filter[[tag_nm]]))
    }
    manifest <- manifest[keep, , drop = FALSE]
  }

  if (nrow(manifest) == 0L) {
    message("No runs matched the filter.")
    return(invisible(list()))
  }

  diag_list <- lapply(manifest$run_name, function(nm) {
    rds_path <- file.path(abl_dir, nm, "diagnostics.rds")
    if (!file.exists(rds_path)) {
      warning(sprintf("diagnostics.rds not found for '%s', skipping.", nm))
      return(NULL)
    }
    message("loading: ", nm)
    readRDS(rds_path)
  })
  names(diag_list) <- manifest$run_name
  diag_list <- Filter(Negate(is.null), diag_list)
  message(sprintf("Loaded %d run(s).", length(diag_list)))
  diag_list
}

# =============================================================================
# 11. compare_ablations()
# =============================================================================

#' @export
compare_ablations <- function(diag_list, tag = "comparison",
                              width = 10, height = 7, dpi = 150) {

  stopifnot(is.list(diag_list), length(diag_list) >= 1L)

  pkg_root <- rprojroot::find_package_root_file(path = ".")
  cmp_dir  <- file.path(pkg_root, "inst", "data", "ablations", "compare", tag)
  dir.create(cmp_dir, showWarnings = FALSE, recursive = TRUE)

  .save_plot <- function(p, fname, w = width, h = height) {
    path <- file.path(cmp_dir, fname)
    ggplot2::ggsave(path, plot = p, width = w, height = h, dpi = dpi)
    message("  saved: ", path)
  }

  comp     <- compare_diagnostics(diag_list)
  csv_path <- file.path(cmp_dir, "summary.csv")
  write.csv(comp, csv_path, row.names = FALSE)
  message("summary saved: ", csv_path)
  print(comp, digits = 4)

  if (length(diag_list) > 1L) {
    for (field in c("posterior_mean", "error", "posterior_se", "in_ci")) {
      p_cmp <- plot_compare_maps(diag_list, field = field)
      .save_plot(p_cmp, sprintf("compare_%s.png", field),
                 w = 4 * length(diag_list), h = 5)
    }
  }

  message("Comparison outputs in: ", normalizePath(cmp_dir, mustWork = FALSE))
  invisible(comp)
}

# tiny helper
`%||%` <- function(a, b) if (!is.null(a)) a else b
