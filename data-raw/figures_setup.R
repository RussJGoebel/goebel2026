# data-raw/figures/shared_setup.R
#
# Shared setup for all figure scripts.
# Source this at the top of each figure script.
# Loads all required results and builds common spatial objects.

library(goebel2026)
library(ggplot2)
library(sf)
library(Matrix)

# ------------------------------------------------------------------------------
# Common spatial objects
# ------------------------------------------------------------------------------

target_sf          <- goebel2026::target_grid
fine_grid_buffered <- goebel2026::setup_shared$fine_grid_buffered
soundings_proj     <- goebel2026::setup_shared$soundings_proj
soundings_proj <- soundings_proj[-which.max(soundings_proj$SIF_757nm),]
target_idx         <- which(target_sf$n_intersects > 0)

# Filter soundings_proj

soundings_proj <- soundings_proj[-which.max(soundings_proj$SIF_757nm),]

# Target grid for plotting (non-buffered, target pixels only)
plot_sf <- target_sf[target_idx, ]

# ------------------------------------------------------------------------------
# Common data
# ------------------------------------------------------------------------------

y_latent_true <- goebel2026::setup_albedo$y_latent_true
y_alb         <- goebel2026::setup_albedo$y
y_sif         <- goebel2026::setup_sif$y

# ------------------------------------------------------------------------------
# Load results
# ------------------------------------------------------------------------------

# Albedo
alb_canonical  <- goebel2026::results_10m_water_rsr_rho_cv
alb_no_rsr     <- goebel2026::results_10m_water_rho1
alb_rsr_rho1   <- goebel2026::results_10m_water_rsr_rho1
alb_no_cov     <- goebel2026::results_10m_no_cov_rho1
alb_ols        <- goebel2026::results_ols_baseline

# SIF
sif_canonical  <- goebel2026::results_sif_canonical

# ------------------------------------------------------------------------------
# Common theme
# ------------------------------------------------------------------------------

theme_map <- function() {
  theme_void() +
    theme(
      legend.title     = element_text(size = 8),
      legend.text      = element_text(size = 7),
      plot.title       = element_text(size = 9, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = 7),
      plot.margin      = margin(2, 2, 2, 2)
    )
}

# ------------------------------------------------------------------------------
# Output directory
# ------------------------------------------------------------------------------

fig_dir <- file.path(getwd(), "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

save_fig <- function(p, name, width = 8, height = 5) {
  path <- file.path(fig_dir, name)
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300)
  message(sprintf("  saved: %s", path))
}

message("shared_setup.R loaded.")
