# data-raw/figures/figures.R
#
# Reproduces all main paper figures using plot_sf_on_google_satellite_fill_palette
# from SpatialBasis (function definition pasted below for portability).
#
# Usage:
#   1. Run shared_setup.R chunks interactively first
#   2. Run each figure section below
#
# Figures:
#   fig01_soundings.pdf         -- OCO-2 geometries + intersection counts
#   fig02_albedo_results.pdf    -- albedo true / observations / posterior mean
#   fig03_albedo_scatter.pdf    -- true vs predicted albedo, colored by water
#   fig04_sif_results.pdf       -- 6-panel SIF main result
#   fig05_sif_ndvi.pdf          -- SIF vs NDVI by landcover

if (!exists("plot_sf")) source("data-raw/figures_setup.R")

library(ggmap)
library(ggspatial)
library(patchwork)
library(dplyr)
library(grid)
library(PaperResults)

# ------------------------------------------------------------------------------
# paste plot_sf_on_google_satellite_fill_palette here or load from SpatialBasis
# ------------------------------------------------------------------------------

if (!exists("plot_sf_on_google_satellite_fill_palette")) {
  # Paste the function definition here if SpatialBasis not installed
  # (see document shared by Russell)
  stop("plot_sf_on_google_satellite_fill_palette not found.
       Either library(SpatialBasis) or paste the function definition above.")
}

# ------------------------------------------------------------------------------
# Shared parameters matching original paper (from original Rmd)
# ------------------------------------------------------------------------------

P <- list(
  api_key             = Sys.getenv("GOOGLE_MAPS_KEY"),
  zoom                = 10,
  frac                = 0.05,
  scalebar_pad_in     = c(0.08, 0.08),
  scalebar_text_cex   = 1.2,
  default_alpha       = 1.0,
  outline_color       = NA,
  legend_barheight_mm = 70,
  legend_barwidth_mm  = 6,
  legend_title_size   = 14,
  legend_text_size    = 11,
  # Palettes from original paper
  pal_gray            = c("black", "gray", "lightgray", "white"),
  pal_sif             = c("black", "lightblue", "white"),
  pal_ndvi            = c("black", "darkgreen", "lightgreen"),
  pal_se              = c("black", "red", "orange", "yellow")
)

# Thin wrapper replicating original map_plot() helper
map_plot <- function(data, column, palette,
                     lim            = NULL,
                     no_legend      = FALSE,
                     legend_title   = column,
                     discrete_levels = NULL,
                     alpha          = P$default_alpha,
                     outline_color  = P$outline_color) {

  p <- plot_sf_on_google_satellite_fill_palette(
    sf_poly           = data,
    fill_col          = column,
    api_key           = P$api_key,
    zoom              = P$zoom,
    palette           = palette,
    frac              = P$frac,
    scalebar_pad_in   = P$scalebar_pad_in,
    scalebar_text_cex = P$scalebar_text_cex,
    alpha             = alpha,
    outline_color     = outline_color,
    limits            = lim,
    no_legend         = no_legend,
    legend_title      = legend_title,
    discrete_levels   = discrete_levels
  )

  if (!no_legend) {
    p <- p +
      guides(fill = guide_colorbar(
        barheight      = unit(P$legend_barheight_mm, "mm"),
        barwidth       = unit(P$legend_barwidth_mm,  "mm"),
        title.position = "top"
      )) +
      theme(
        legend.title = element_text(size = P$legend_title_size),
        legend.text  = element_text(size = P$legend_text_size)
      )
  }
  p
}

# ------------------------------------------------------------------------------
# Shared data prep: add result columns to target_grid
# ------------------------------------------------------------------------------

target_grid_plot <- goebel2026::target_grid
idx <- target_grid_plot$n_intersects > 0

# Albedo
A_flat       <- goebel2026::setup_shared$A_flat
obs_at_pixel <- as.numeric(Matrix::crossprod(A_flat[, target_idx], y_alb)) /
  pmax(Matrix::colSums(A_flat[, target_idx]), 1e-10)

target_grid_plot$true_alb  <- NA_real_
target_grid_plot$obs_alb   <- NA_real_
target_grid_plot$pred_alb  <- NA_real_
target_grid_plot$true_alb[idx]  <- y_latent_true[target_idx]
target_grid_plot$obs_alb[idx]   <- obs_at_pixel
target_grid_plot$pred_alb[idx]  <- alb_canonical$posterior_mean

# SIF
target_grid_plot$sif_mean <- NA_real_
target_grid_plot$sif_se   <- NA_real_
target_grid_plot$sif_mean[idx] <- sif_canonical$posterior_mean
target_grid_plot$sif_se[idx]   <- sif_canonical$posterior_se

# NDVI already in target_grid as mean_ndvi

# Landcover 3-class
# dominant_class codes: 1=Urban, 8=Bare -> Built; 6=Water, 7=Wetlands -> Water;
#                       2=Croplands, 3=Forest, 4=Shrub, 5=Grass -> Green
dc <- target_grid_plot$dominant_class
target_grid_plot$dom_class <- factor(
  dplyr::case_when(
    dc %in% c(1, 8)       ~ "Built",
    dc %in% c(6, 7)       ~ "Water",
    dc %in% c(2, 3, 4, 5) ~ "Green",
    TRUE                  ~ NA_character_
  ),
  levels = c("Built", "Water", "Green")
)

# Shared limits
alb_lims <- range(c(target_grid_plot$true_alb,
                    target_grid_plot$pred_alb), na.rm = TRUE)
# SIF limits: match old figure range 0-4 for consistent color scale
# Posterior mean is smoother so has smaller range, but use same scale for comparability
sif_lims <- c(0, 4)
se_lims  <- range(target_grid_plot$sif_se, na.rm = TRUE)

tg <- target_grid_plot[idx, ]   # target pixels only

# ==============================================================================
# Fig 1: OCO-2 sounding geometries + intersection counts
# ==============================================================================

message("\n== Fig 1: Sounding geometries ==")

soundings_sf <- goebel2026::setup_shared$soundings_proj
soundings_sf$fill_col <- 1L

p1a <- plot_sf_on_google_satellite_fill_palette(
  sf_poly       = soundings_sf,
  fill_col      = "fill_col",
  api_key       = P$api_key,
  zoom          = P$zoom,
  palette       = c("blue", "blue"),
  frac          = P$frac,
  alpha         = 0.4,
  outline_color = "white",
  no_legend     = TRUE,
  scalebar_pad_in   = c(0.08, 0.08),
  scalebar_text_cex = 1.2
)

p1b <- map_plot(
  data         = tg,
  column       = "n_intersects",
  palette      = "viridis",
  legend_title = "Intersections"
)

save_fig(p1a, "fig01_1_soundings_geometries.pdf",  width = 5, height = 5)
save_fig(p1b, "fig01_2_intersection_counts.pdf",   width = 5, height = 5)

# ==============================================================================
# Fig 2: Semi-synthetic albedo results
# ==============================================================================

message("\n== Fig 2: Albedo results ==")

soundings_aug <- goebel2026::soundings_augmented
soundings_aug <- soundings_aug[-which.max(soundings_aug$SIF_757nm),]
soundings_aug <- sf::st_transform(soundings_aug, sf::st_crs(tg))
soundings_aug$synthetic_albedo <- y_alb

p2a <- map_plot(tg, "true_alb", P$pal_gray, lim = alb_lims,
                legend_title = "Mean\nAlbedo", no_legend = TRUE)

p2b <- plot_sf_on_google_satellite_fill_palette(
  sf_poly       = soundings_aug,
  fill_col      = "synthetic_albedo",
  api_key       = P$api_key,
  zoom          = P$zoom,
  palette       = P$pal_gray,
  frac          = P$frac,
  alpha         = 1.0,
  outline_color = NULL,
  limits        = alb_lims,
  no_legend     = TRUE,
  scalebar_pad_in   = P$scalebar_pad_in,
  scalebar_text_cex = P$scalebar_text_cex
)

p2c <- map_plot(tg, "pred_alb", P$pal_gray, lim = alb_lims,
                legend_title = "Mean\nAlbedo")

save_fig(p2a, "fig02_1_albedo_true.pdf",        width = 5, height = 5)
save_fig(p2b, "fig02_2_albedo_observations.pdf", width = 5, height = 5)
save_fig(p2c, "fig02_3_albedo_posterior.pdf",    width = 5, height = 5)

# Shared colorbar for albedo panels
p2_cbar <- map_plot(tg, "true_alb", P$pal_gray, lim = alb_lims,
                    legend_title = "Mean\nAlbedo")
save_fig(cowplot::get_legend(p2_cbar), "fig02_colorbar.pdf", width = 2, height = 4)

# ==============================================================================
# Fig 3: True vs predicted albedo scatterplot
# ==============================================================================

message("\n== Fig 3: Albedo scatterplot ==")

scatter_df <- data.frame(
  true  = y_latent_true[target_idx],
  pred  = alb_canonical$posterior_mean,
  water = goebel2026::target_grid$proportion_water[target_idx]
)
lims3 <- range(c(scatter_df$true, scatter_df$pred), na.rm = TRUE)

fig03 <- ggplot(scatter_df, aes(x = true, y = pred, color = water)) +
  geom_point(alpha = 0.4, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  scale_color_viridis_c(name = "Proportion\nWater",
                        limits = c(0, 1),
                        breaks = c(0, 0.25, 0.50, 0.75, 1.00)) +
  scale_x_continuous(limits = lims3) +
  scale_y_continuous(limits = lims3) +
  labs(x = "Mean Albedo", y = "Posterior Mean Albedo Prediction",
       subtitle = sprintf("R\u00b2 = %.2f", alb_canonical$r2)) +
  coord_fixed() +
  theme_minimal(base_size = 11)

save_fig(fig03, "fig03_albedo_scatter.pdf", width = 5.5, height = 5)

# ==============================================================================
# Fig 4: SIF main result -- 6 panels saved individually
# ==============================================================================

message("\n== Fig 4: SIF main result ==")

# Top-left: pure satellite -- use map_plot with alpha=0 so extent matches other panels
p4_sat <- map_plot(tg, "sif_mean", P$pal_sif,
                   lim = sif_lims, no_legend = TRUE, alpha = 0)

p4_mean <- map_plot(tg, "sif_mean", P$pal_sif,
                    lim = sif_lims, legend_title = "SIF")

p4_ndvi <- map_plot(tg, "mean_ndvi", P$pal_ndvi,
                    legend_title = "NDVI")

# SIF observations -- plot sounding polygons (coarse footprints), not target grid
soundings_aug$SIF_757nm <- goebel2026::soundings_augmented$SIF_757nm[-which.max(goebel2026::soundings_augmented$SIF_757nm)]
p4_obs <- plot_sf_on_google_satellite_fill_palette(
  sf_poly = soundings_aug, fill_col = "SIF_757nm",
  api_key = P$api_key, zoom = P$zoom,
  palette = P$pal_sif, frac = P$frac,
  limits = sif_lims, alpha = 1.0,
  outline_color = NULL, no_legend = TRUE,
  scalebar_pad_in = P$scalebar_pad_in,
  scalebar_text_cex = P$scalebar_text_cex
)

p4_lc <- map_plot(
  tg, "dom_class",
  palette         = c('#c4281b', '#429ae4', '#397e48'),
  discrete_levels = c("Built", "Water", "Green"),
  no_legend       = TRUE,
  legend_title    = "Landcover"
)

p4_se <- map_plot(tg, "sif_se", P$pal_se,
                  lim = se_lims, legend_title = "Std. Error")

save_fig(p4_sat,  "fig04_1_boston_satellite.pdf",  width = 5, height = 5)
save_fig(p4_mean, "fig04_2_sif_posterior_mean.pdf", width = 5, height = 5)
save_fig(p4_ndvi, "fig04_3_ndvi.pdf",               width = 5, height = 5)
save_fig(p4_obs,  "fig04_4_sif_observations.pdf",   width = 5, height = 5)
save_fig(p4_lc,   "fig04_5_landcover.pdf",           width = 5, height = 5)
save_fig(p4_se,   "fig04_6_posterior_se.pdf",        width = 5, height = 5)

# ==============================================================================
# Fig 5: SIF vs NDVI -- two panels with density scatter + binned quantile step
# ==============================================================================

message("\n== Fig 5: SIF vs NDVI ==")

library(cowplot)

qf <- goebel2026::quantile_fits_canonical

# --- Scatter data (unclipped first, used for axis limit computation) ----------
snd_df <- data.frame(
  sif   = qf$y_sif,
  ndvi  = qf$ndvi_soundings,
  water = qf$is_water_soundings
) |> dplyr::filter(is.finite(sif) & is.finite(ndvi))

post_df <- data.frame(
  ndvi  = qf$ndvi_target,
  sif   = qf$preds,
  water = qf$is_water_target
) |> dplyr::filter(is.finite(ndvi) & is.finite(sif))

# --- Shared axis limits -------------------------------------------------------
ndvi_lims <- c(
  min(c(snd_df$ndvi, post_df$ndvi), na.rm = TRUE) - 0.02,
  max(c(snd_df$ndvi, post_df$ndvi), na.rm = TRUE) + 0.02
)

SIF_YMAX <- max(
  quantile(c(snd_df$sif, post_df$sif), 0.995, na.rm = TRUE),
  max(qf$snd_quantiles[, "q0.975"], na.rm = TRUE),
  max(qf$post_q_mean[,  "q0.975"], na.rm = TRUE)
)

sif_ymin <- min(
  min(qf$snd_quantiles[, "q0.025"], na.rm = TRUE),
  min(qf$post_q_mean[,  "q0.025"], na.rm = TRUE)
) - 0.05

sif_lims_plot <- c(sif_ymin, SIF_YMAX)

message(sprintf("  Axis limits: NDVI=[%.2f, %.2f]  SIF=[%.2f, %.2f]",
                ndvi_lims[1], ndvi_lims[2], sif_lims_plot[1], sif_lims_plot[2]))

# --- Clipped scatter data for geom_point and marginals -----------------------
snd_df_clip <- dplyr::filter(snd_df,
                             ndvi >= ndvi_lims[1],    ndvi <= ndvi_lims[2],
                             sif  >= sif_lims_plot[1], sif  <= sif_lims_plot[2])

post_df_clip <- dplyr::filter(post_df,
                              ndvi >= ndvi_lims[1],    ndvi <= ndvi_lims[2],
                              sif  >= sif_lims_plot[1], sif  <= sif_lims_plot[2])

# --- Step median lines -------------------------------------------------------
step_snd_med <- data.frame(
  x = qf$bin_mids,
  y = qf$snd_quantiles[, "q0.5"],
  n = qf$snd_counts
)
step_snd_med <- step_snd_med[!is.na(step_snd_med$y) & step_snd_med$n >= 10, ]

step_post_med <- data.frame(
  x = qf$bin_mids,
  y = qf$post_q_mean[, "q0.5"]
)
step_post_med <- step_post_med[!is.na(step_post_med$y), ]

# --- Helper: geom_rect step bands --------------------------------------------
make_step_rects <- function(bin_mids, quant_mat, lo_col, hi_col,
                            bin_breaks, fill, alpha, counts = NULL, min_n = 10) {
  df <- data.frame(
    xmin = bin_breaks[-length(bin_breaks)],
    xmax = bin_breaks[-1L],
    ymin = quant_mat[, lo_col],
    ymax = quant_mat[, hi_col]
  )
  if (!is.null(counts)) df$n <- counts
  df <- df[complete.cases(df), ]
  if (!is.null(counts)) df <- df[df$n >= min_n, ]
  geom_rect(data = df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = fill, alpha = alpha)
}

# Larger base_size for Fig 5 panels
theme_sif_ndvi <- theme_minimal(base_size = 14) +
  theme(panel.grid.minor = element_blank(),
        axis.title       = element_text(size = 12))

ylab_sif      <- expression("SIF"[757]~"[mW m"^{-2}~"nm"^{-1}~"sr"^{-1}*"]")
ylab_sif_post <- expression("Posterior SIF"[757]~"[mW m"^{-2}~"nm"^{-1}~"sr"^{-1}*"]")

# ------------------------------------------------------------------------------
# Left panel: x=SIF, y=NDVI (flipped), bins on SIF axis
# ------------------------------------------------------------------------------
p5a_base <- ggplot() +
  # Outer band (2.5-97.5) light gray -- SIF bins, NDVI quantiles -> horizontal rects
  # xmin/xmax = NDVI quantiles, ymin/ymax = SIF bin edges
  {
    df_outer <- data.frame(
      xmin = qf$snd_quantiles[, "q0.025"],
      xmax = qf$snd_quantiles[, "q0.975"],
      ymin = qf$sif_breaks[-length(qf$sif_breaks)],
      ymax = qf$sif_breaks[-1L],
      n    = qf$snd_counts
    )
    df_outer <- df_outer[is.finite(df_outer$xmin) & is.finite(df_outer$xmax) &
                           df_outer$n >= 10, ]
    geom_rect(data = df_outer,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "grey70", alpha = 0.35, inherit.aes = FALSE)
  } +
  # Inner band (25-75)
  {
    df_inner <- data.frame(
      xmin = qf$snd_quantiles[, "q0.25"],
      xmax = qf$snd_quantiles[, "q0.75"],
      ymin = qf$sif_breaks[-length(qf$sif_breaks)],
      ymax = qf$sif_breaks[-1L],
      n    = qf$snd_counts
    )
    df_inner <- df_inner[is.finite(df_inner$xmin) & is.finite(df_inner$xmax) &
                           df_inner$n >= 10, ]
    geom_rect(data = df_inner,
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
              fill = "grey40", alpha = 0.35, inherit.aes = FALSE)
  } +
  # Non-water points
  geom_point(data = dplyr::filter(snd_df_clip, !water),
             aes(x = ndvi, y = sif),
             size = 0.3, alpha = 0.2, color = "grey20") +
  # Water crosses
  geom_point(data = dplyr::filter(snd_df_clip, water),
             aes(x = ndvi, y = sif),
             shape = 4, size = 0.8, color = "grey50", alpha = 0.6) +
  # Median step -- x=NDVI median, y=SIF bin mid -> horizontal step
  {
    df_med <- data.frame(
      x = qf$snd_quantiles[, "q0.5"],
      y = qf$sif_bin_mids,
      n = qf$snd_counts
    )
    df_med <- df_med[is.finite(df_med$x) & df_med$n >= 10, ]
    geom_step(data = df_med, aes(x = x, y = y),
              color = "red", linewidth = 0.8, direction = "vh")
  } +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("R\u00b2 = %.3f", qf$r2_soundings),
           size = 5, color = "steelblue", fontface = "italic") +
  scale_x_continuous(limits = ndvi_lims, expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  coord_cartesian(ylim = sif_lims_plot) +
  labs(x = "NDVI", y = ylab_sif) +
  theme_sif_ndvi

# ------------------------------------------------------------------------------
# Right panel: x=NDVI, y=posterior SIF, bins on NDVI axis
# ------------------------------------------------------------------------------
p5b_base <- ggplot() +
  # Outer band
  make_step_rects(qf$bin_mids, qf$post_q_mean,
                  "q0.025", "q0.975", qf$bin_breaks,
                  fill = "grey70", alpha = 0.35) +
  # Inner band
  make_step_rects(qf$bin_mids, qf$post_q_mean,
                  "q0.25", "q0.75", qf$bin_breaks,
                  fill = "grey40", alpha = 0.35) +
  # Non-water points
  geom_point(data = dplyr::filter(post_df_clip, !water),
             aes(x = ndvi, y = sif),
             size = 0.3, alpha = 0.2, color = "grey20") +
  # Water crosses
  geom_point(data = dplyr::filter(post_df_clip, water),
             aes(x = ndvi, y = sif),
             shape = 4, size = 0.8, color = "grey50", alpha = 0.6) +
  # Median step
  geom_step(data = step_post_med, aes(x = x, y = y),
            color = "red", linewidth = 0.8) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = sprintf("R\u00b2 = %.3f \u00b1 %.3f", qf$r2_mean, qf$r2_sd),
           size = 5, color = "steelblue", fontface = "italic") +
  scale_x_continuous(limits = ndvi_lims, expand = c(0.01, 0)) +
  scale_y_continuous(expand = c(0.01, 0)) +
  coord_cartesian(ylim = sif_lims_plot) +
  labs(x = "NDVI", y = ylab_sif_post) +
  theme_sif_ndvi

# Marginal histograms
# Filter to axis limits for marginals to align correctly
# Build marginals manually with cowplot for correct axis alignment
library(cowplot)

make_xhist <- function(df) {
  ggplot(df, aes(x = ndvi)) +
    geom_histogram(fill = "grey70", color = "white", bins = 30) +
    scale_x_continuous(limits = ndvi_lims, expand = c(0.01, 0)) +
    theme_void()
}

make_yhist <- function(df) {
  ggplot(df, aes(x = sif)) +
    geom_histogram(fill = "grey70", color = "white", bins = 30) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_flip(xlim = sif_lims_plot) +
    theme_void()
}

p5a_final <- cowplot::insert_xaxis_grob(p5a_base, make_xhist(snd_df_clip),
                                        position = "top", height = unit(0.2, "null"))
p5a_final <- cowplot::insert_yaxis_grob(p5a_final, make_yhist(snd_df_clip),
                                        position = "right", width = unit(0.2, "null"))
p5a <- cowplot::ggdraw(p5a_final)

p5b_final <- cowplot::insert_xaxis_grob(p5b_base, make_xhist(post_df_clip),
                                        position = "top", height = unit(0.2, "null"))
p5b_final <- cowplot::insert_yaxis_grob(p5b_final, make_yhist(post_df_clip),
                                        position = "right", width = unit(0.2, "null"))
p5b <- cowplot::ggdraw(p5b_final)

save_fig(p5a, "fig05_1_sif_vs_ndvi_soundings.pdf",  width = 5.5, height = 5.5)
save_fig(p5b, "fig05_2_sif_vs_ndvi_posterior.pdf",  width = 5.5, height = 5.5)

# ==============================================================================

message("\nAll panels saved to figures/:")
message("  fig01_1_soundings_geometries.pdf")
message("  fig01_2_intersection_counts.pdf")
message("  fig02_1_albedo_true.pdf")
message("  fig02_2_albedo_observations.pdf")
message("  fig02_3_albedo_posterior.pdf")
message("  fig02_colorbar.pdf")
message("  fig03_albedo_scatter.pdf")
message("  fig04_1_boston_satellite.pdf")
message("  fig04_2_sif_posterior_mean.pdf")
message("  fig04_3_ndvi.pdf")
message("  fig04_4_sif_observations.pdf")
message("  fig04_5_landcover.pdf")
message("  fig04_6_posterior_se.pdf")
message("  fig05_1_sif_vs_ndvi_soundings.pdf")
message("  fig05_2_sif_vs_ndvi_posterior.pdf")
