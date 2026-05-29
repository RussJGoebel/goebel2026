---
title: "Supplementary Materials: Bayesian Downscaling of Satellite-Derived Solar-Induced Chlorophyll Fluorescence"
author: "Russell Goebel, Leeza Moldavchuk, Taylor S. Jones, Jonathan Dooley, Lucy R. Hutyra, Luis Carvalho"
output:
  pdf_document:
    toc: true
    toc_depth: 2
    number_sections: true
    latex_engine: pdflatex
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo    = FALSE,
  message = FALSE,
  warning = FALSE,
  fig.align = "center",
  out.width = "\\textwidth"
)

library(goebel2026)
library(ggplot2)
library(patchwork)
library(dplyr)
library(knitr)
library(scales)
library(sf)

# Shared formatting helpers
fmt2 <- function(x) sprintf("%.2f", x)
fmt3 <- function(x) sprintf("%.3f", x)
fmt4 <- function(x) sprintf("%.4f", x)
pct  <- function(x) sprintf("%.1f\\%%", 100 * x)

# Load all result objects fresh from package (avoids stale session cache)
.pkg_get <- function(nm) get(nm, envir = asNamespace("goebel2026"))

r <- list(
  # Albedo ablation
  ols             = .pkg_get("results_ols_baseline"),
  no_cov          = .pkg_get("results_no_cov_rho1"),
  water           = .pkg_get("results_water_rho1"),
  water_rsr       = .pkg_get("results_water_rsr_rho1"),
  water_rsr_cv    = .pkg_get("results_water_rsr_rho_cv"),
  water_rho_cv    = .pkg_get("results_water_rho_cv"),
  # SIF
  sif             = .pkg_get("results_sif_canonical"),
  sif_naive       = .pkg_get("results_sif_naive"),
  sif_rsr_rho1    = .pkg_get("results_sif_rsr_rho1"),
  # Tuning
  cv_blocked      = .pkg_get("results_cv_blocked"),
  water_ml        = .pkg_get("results_water_ml"),
  sif_cv_blocked  = .pkg_get("results_sif_cv_blocked"),
  sif_water_ml    = .pkg_get("results_sif_water_ml"),
  sif_rho095_ml   = .pkg_get("results_sif_rho095_ml"),
  sif_rho095_cv   = .pkg_get("results_sif_rho095_cv"),
  # Neighbours
  nb_rook         = .pkg_get("results_neighbor_rook"),
  nb_lc_eps       = .pkg_get("results_neighbor_lc_eps"),
  nb_lc_alpha     = .pkg_get("results_neighbor_lc_alpha"),
  sif_nb_rook     = .pkg_get("results_sif_neighbor_rook"),
  sif_nb_lc_eps   = .pkg_get("results_sif_neighbor_lc_eps"),
  sif_nb_lc_alpha = .pkg_get("results_sif_neighbor_lc_alpha"),
  # Forward operator
  alb_gA_tau05    = .pkg_get("results_albedo_gA_tau05"),
  alb_gA_tau033   = .pkg_get("results_albedo_gA_tau033"),
  alb_gA_tau02    = .pkg_get("results_albedo_gA_tau02"),
  alb_gA_tau01    = .pkg_get("results_albedo_gA_tau01"),
  alb_gA_cent     = .pkg_get("results_albedo_gA_centroid"),
  sif_gA_tau05    = .pkg_get("results_sif_gA_tau05"),
  sif_gA_tau033   = .pkg_get("results_sif_gA_tau033"),
  sif_gA_tau02    = .pkg_get("results_sif_gA_tau02"),
  sif_gA_tau01    = .pkg_get("results_sif_gA_tau01"),
  sif_gA_tau001   = .pkg_get("results_sif_gA_tau001"),
  # Kriging
  krig_alb        = .pkg_get("results_kriging_albedo"),
  krig_sif        = .pkg_get("results_kriging_sif"),
  # Soft RSR
  softrsr         = .pkg_get("results_softrsr_albedo")
)

# Likelihood surfaces
ls_alb <- .pkg_get("likelihood_surface_albedo")
ls_sif <- .pkg_get("likelihood_surface_sif")

# Safe rho extraction -- fixed-rho results don't store rho_opt
get_rho <- function(x, default = 1.0) {
  v <- x$rho_opt
  if (is.null(v) || length(v) == 0) default else v
}

# Safe field extraction -- returns NA if field is NULL/missing
gf <- function(x, field, idx = NULL) {
  v <- x[[field]]
  if (is.null(v) || length(v) == 0) return(NA_real_)
  if (!is.null(idx)) v <- v[idx]
  if (length(v) == 0) NA_real_ else v
}

# Safe field extraction
get_field <- function(x, ...) {
  fields <- c(...)
  for (f in fields) {
    v <- x[[f]]
    if (!is.null(v) && length(v) > 0) return(v)
  }
  NA_real_
}

```

\newpage

# Hyperparameter Tuning {#sec:tuning}

The smoothing parameter $\phi$ controls the signal-to-noise ratio of the
spatial random effect and must be selected from data. We compare three
approaches: (1) random $K$-fold cross-validation (random CV), which assigns
soundings to folds uniformly at random; (2) spatially blocked cross-validation
(blocked CV), which assigns soundings to contiguous spatial blocks; and (3)
marginal likelihood maximisation (ML). Throughout this section we fix
$\rho = 0.95$ and include no covariates or RSR constraint, so that differences
in selected $\phi$ reflect the tuning method alone.

**Spatial blocking criterion.** For blocked CV, soundings are assigned to
folds by $k$-means clustering on their centroid coordinates, with $k = 10$.
This produces spatially contiguous blocks whose sizes reflect the local
sounding density — denser regions (the centre of the target-mode acquisition)
yield larger blocks. Figure 2 shows the resulting fold assignments overlaid
on the sounding geometries, coloured by fold index. The concentration of
soundings in the central strip means that a single held-out block removes a
disproportionately large and informationally dense region from the training
set, which is the root cause of the pathological $\phi$ selection documented
below.

```{r fold_map, fig.cap="Spatial fold assignments for blocked (left) and random (right) 10-fold cross-validation. Each sounding is coloured by its fold index. Blocked folds form contiguous spatial regions; random folds are scattered uniformly. The central high-density strip falls predominantly within one or two blocked folds, meaning that holding out a single block removes the most informative soundings from the training set.", fig.height=4, fig.width=9}

soundings_sf  <- goebel2026::setup_shared$soundings_proj
folds_blocked <- goebel2026::setup_albedo$blocked_folds
set.seed(2026L)
folds_random  <- sample(rep_len(1:10, nrow(soundings_sf)))

soundings_sf$fold_blocked <- factor(folds_blocked)
soundings_sf$fold_random  <- factor(folds_random)

# Bounding box for consistent extent
bb <- sf::st_bbox(soundings_sf)

make_fold_plot <- function(sf_obj, fold_col, title) {
  ggplot(sf_obj) +
    geom_sf(aes(fill = .data[[fold_col]]), color = NA, alpha = 0.8) +
    scale_fill_manual(
      values = setNames(
        scales::hue_pal()(10),
        as.character(1:10)
      ),
      name = "Fold", guide = guide_legend(ncol = 2)
    ) +
    coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
             ylim = c(bb["ymin"], bb["ymax"]),
             expand = FALSE) +
    labs(title = title) +
    theme_void(base_size = 9) +
    theme(legend.position  = "right",
          plot.title       = element_text(size = 10, hjust = 0.5))
}

p_blk <- make_fold_plot(soundings_sf, "fold_blocked", "Blocked CV")
p_rnd <- make_fold_plot(soundings_sf, "fold_random",  "Random CV")

p_blk | p_rnd
```

## Likelihood and CV objective profiles

Figure 2 shows the ML log-likelihood and CV MSE
profiles as functions of $\log\phi$ for both albedo (left) and SIF (right).
All three objectives are scaled to $[0, 1]$ and oriented so that higher is
better.

For both responses, the ML and random CV profiles peak at similar values of
$\phi$: the ratio of ML to random CV optima is
`r fmt2(ls_alb$ml_phi_opt / ls_alb$cv_phi_random)`x for albedo and
`r fmt2(ls_sif$ml_phi_opt / ls_sif$cv_phi_random)`x for SIF, indicating
strong agreement between the two methods. Both profiles have broad, flat
peaks, suggesting that predictions are relatively insensitive to the exact
value of $\phi$ in this region.

Spatially blocked CV selects a substantially different $\phi$: its optimum is
`r fmt2(ls_alb$cv_phi_blocked / ls_alb$cv_phi_random)`x larger than the
random CV optimum for albedo and
`r fmt2(ls_sif$cv_phi_blocked / ls_sif$cv_phi_random)`x larger for SIF.
This over-smoothing occurs because removing a contiguous spatial block removes
the most densely observed region of the domain, leaving the held-out pixels
far from any training observations and forcing the model to rely on heavy
spatial smoothing to achieve reasonable predictions. Random CV, by contrast,
holds out scattered soundings that always have nearby training neighbours,
producing a CV objective that more faithfully reflects the model's predictive
accuracy at the observed density.

```{r likelihood_surface, fig.cap="ML log-likelihood (red) and CV MSE profiles (random CV: blue; blocked CV: green) as functions of $\\log\\phi$, for albedo (left) and SIF (right). All curves are scaled to $[0,1]$ and oriented so that higher is better. Dashed vertical lines indicate the optimum for each method. Random CV and ML optima are in close agreement; blocked CV selects a substantially larger $\\phi$ in both cases.", fig.height=4, fig.width=9}

make_surface_panel <- function(surf, title) {
  lp  <- surf$log_phi_seq

  ml  <- surf$ml_ll
  ml  <- (ml  - max(ml,  na.rm=TRUE)) /
          abs(diff(range(ml,  na.rm=TRUE)))

  cvr <- surf$cv_mse_random
  cvr <- -(cvr - min(cvr, na.rm=TRUE)) /
           abs(diff(range(cvr, na.rm=TRUE)))

  cvb <- surf$cv_mse_blocked
  cvb <- -(cvb - min(cvb, na.rm=TRUE)) /
           abs(diff(range(cvb, na.rm=TRUE)))

  df <- data.frame(
    log_phi   = rep(lp, 3),
    value     = c(ml, cvr, cvb),
    method    = rep(c("ML", "CV random", "CV blocked"),
                    each = length(lp))
  )
  df$method <- factor(df$method,
                      levels = c("ML", "CV random", "CV blocked"))

  opts <- data.frame(
    method    = factor(c("ML", "CV random", "CV blocked"),
                       levels = c("ML", "CV random", "CV blocked")),
    log_phi   = log(c(surf$ml_phi_opt,
                      surf$cv_phi_random,
                      surf$cv_phi_blocked)),
    phi_label = c(
      sprintf("ML: %.3g",    surf$ml_phi_opt),
      sprintf("Rand: %.3g",  surf$cv_phi_random),
      sprintf("Block: %.3g", surf$cv_phi_blocked)
    )
  )

  # Build subtitle with phi values instead of overlapping text annotations
  subtitle_text <- sprintf(
    "ML: %.3g  |  Rand: %.3g  |  Block: %.3g",
    surf$ml_phi_opt, surf$cv_phi_random, surf$cv_phi_blocked
  )

  ggplot(df[is.finite(df$value), ],
         aes(x = log_phi, y = value, color = method)) +
    geom_line(linewidth = 0.9) +
    geom_vline(data = opts,
               aes(xintercept = log_phi, color = method),
               linetype = "dashed", linewidth = 0.6) +
    scale_color_manual(
      values = c("ML"         = "tomato",
                 "CV random"  = "steelblue",
                 "CV blocked" = "darkgreen"),
      name = NULL) +
    scale_y_continuous(limits = c(-1, 0.05)) +
    labs(x        = expression(log~phi),
         y        = "Scaled objective (higher = better)",
         title    = title,
         subtitle = paste0("phi_opt -- ML: ", round(surf$ml_phi_opt, 3),
                           "  |  Rand: ",     round(surf$cv_phi_random, 3),
                           "  |  Block: ",    round(surf$cv_phi_blocked, 3))) +
    theme_minimal(base_size = 10) +
    theme(legend.position   = "bottom",
          panel.grid.minor  = element_blank(),
          plot.subtitle     = element_text(size = 8, color = "grey40"))
}

p_alb <- make_surface_panel(ls_alb, "Albedo")
p_sif <- make_surface_panel(ls_sif, "SIF")
p_alb | p_sif
```

## Comparison with covariate model

The results above use no covariates or RSR constraint. When the water
proportion covariate and RSR are included (as in the main analysis),
the tuning methods still disagree. Table 1a shows CV results with RSR:
blocked CV selects $\phi$ approximately
`r fmt2(gf(r$water_rsr,"phi") / gf(r$cv_blocked,"phi"))`x smaller than random CV
for albedo and
`r fmt2(gf(r$sif_rsr_rho1,"phi") / gf(r$sif_cv_blocked,"phi"))`x smaller for SIF, with coverage
collapsing from
`r pct(gf(r$water_rsr,"coverage_95_obs"))` to
`r pct(gf(r$cv_blocked,"coverage_95_obs"))` for albedo.

Table 1b shows all three methods without RSR or covariates, at fixed
$\rho=0.95$ — the same conditions as Figure 1. Random CV and ML agree
closely (ratio `r fmt2(ls_alb$ml_phi_opt / ls_alb$cv_phi_random)`x for
albedo, `r fmt2(ls_sif$ml_phi_opt / ls_sif$cv_phi_random)`x for SIF),
while blocked CV selects a substantially larger $\phi$ in both cases.
ML cannot be combined with RSR since RSR is a posterior constraint with
no natural analogue in the marginal likelihood objective; the
`r fmt2(gf(r$sif_rho095_cv,"phi") / gf(r$sif_rho095_ml,"phi"))`x ratio
between ML and random CV at fixed $\rho=0.95$ with RSR+covariates
(not shown) is consistent with the no-RSR comparison.

Crucially, $\hat\beta_\text{water}$ is stable across both tables and both
responses (range:
`r fmt4(min(sapply(list(r$water_rsr, r$cv_blocked, r$water_ml), gf, "beta_hat", idx=2), na.rm=TRUE))`
to
`r fmt4(max(sapply(list(r$water_rsr, r$cv_blocked, r$water_ml), gf, "beta_hat", idx=2), na.rm=TRUE))`
for albedo;
`r fmt4(min(sapply(list(r$sif_rsr_rho1, r$sif_cv_blocked, r$sif_rho095_ml), gf, "beta_hat", idx=2), na.rm=TRUE))`
to
`r fmt4(max(sapply(list(r$sif_rsr_rho1, r$sif_cv_blocked, r$sif_rho095_ml), gf, "beta_hat", idx=2), na.rm=TRUE))`
for SIF), confirming that the tuning method does not affect the
substantive inference. Because random CV directly optimises prediction
accuracy, produces near-nominal coverage in the semi-synthetic
validation, and is consistent with ML at comparable model specifications,
we prefer it for the main analysis.

```{r tuning_comparison}
# Table 1a: CV methods with RSR
tuning_cv_alb <- list(r$water_rsr, r$cv_blocked)
tuning_cv_sif <- list(r$sif_rsr_rho1, r$sif_cv_blocked)

tab_cv <- data.frame(
  Response = c("Albedo", "Albedo", "SIF", "SIF"),
  Method   = rep(c("Random CV", "Blocked CV"), 2),
  rho      = c(sapply(tuning_cv_alb, get_rho),
               sapply(tuning_cv_sif, get_rho)),
  phi      = c(sapply(tuning_cv_alb, gf, "phi"),
               sapply(tuning_cv_sif, gf, "phi")),
  RMSE     = c(sapply(tuning_cv_alb, gf, "rmse"),
               NA_real_, NA_real_),
  Coverage = c(sapply(tuning_cv_alb, gf, "coverage_95_obs"),
               NA_real_, NA_real_),
  mean_SE  = c(NA_real_, NA_real_,
               sapply(tuning_cv_sif, function(x)
                 mean(x$posterior_se, na.rm = TRUE)))
)

kable(tab_cv,
  format    = "pipe",
  digits    = 4,
  caption   = "Table 1a: CV tuning methods with RSR constraint. Albedo results include RMSE and coverage (ground truth known); SIF reports mean posterior SE.",
  col.names = c("Response", "Method", "rho", "phi",
                "RMSE", "Coverage (95%)", "Mean SE"))

# Table 1b: all three methods without RSR, rho=0.95 fixed
# phi optima from likelihood surface objects (run_08)
tab_nrsr <- data.frame(
  Response = c("Albedo", "Albedo", "Albedo",
               "SIF",    "SIF",    "SIF"),
  Method   = rep(c("Random CV", "Blocked CV", "ML"), 2),
  rho      = 0.95,
  phi      = c(ls_alb$cv_phi_random,
               ls_alb$cv_phi_blocked,
               ls_alb$ml_phi_opt,
               ls_sif$cv_phi_random,
               ls_sif$cv_phi_blocked,
               ls_sif$ml_phi_opt)
)

kable(tab_nrsr,
  format    = "pipe",
  digits    = 4,
  caption   = "Table 1b: All three tuning methods without RSR constraint or covariates, at fixed rho=0.95. These are the phi optima shown in Figure 1. Random CV and ML agree closely for both responses; blocked CV selects a substantially larger phi due to the dense OCO-2 overlap structure.",
  col.names = c("Response", "Method", "rho", "phi (opt)"))
```
```

\newpage

# Neighbour Weight Sensitivity {#sec:neighbours}

The SAR prior is defined through a neighbourhood matrix $W$. We consider
four specifications, all using queen adjacency as the base structure.
Let $\mathcal{N}(i)$ denote the queen-adjacency neighbours of cell $i$,
$n_i = |\mathcal{N}(i)|$, $w_i \in [0,1]$ the water proportion of cell $i$,
and $\ell_{ij} = \mathbb{1}[w_i > 0.5 \neq w_j > 0.5]$ a binary indicator
that the edge crosses a land--water boundary. The four specifications are:

1. **Queen**: $\tilde W_{ij} = 1/n_i$ for all $j \in \mathcal{N}(i)$
2. **Rook**: restricts to horizontal/vertical neighbours only, $\tilde W_{ij} = 1/n_i^{(\text{rook})}$
3. **LC-eps**: $\tilde W_{ij} \propto (1 - \ell_{ij}) + \varepsilon\,\ell_{ij}$,\quad $\varepsilon = e^{-\alpha}$,\quad $\alpha$ CV-tuned
4. **LC-alpha**: $\tilde W_{ij} \propto e^{-\alpha |w_i - w_j|} + \delta$,\quad $\delta = 10^{-6}$,\quad $\alpha$ CV-tuned

In both LC-aware specifications $\alpha = 0$ recovers the queen baseline.
The CV-tuned values are
$\hat\alpha_\text{eps} = `r fmt3(gf(r$nb_lc_eps,"alpha_opt"))`$ and
$\hat\alpha_\text{alpha} = `r fmt3(gf(r$nb_lc_alpha,"alpha_opt"))`$,
both close to zero.

Tables 2 and 3 show that all four specifications give nearly identical
predictive performance for albedo and SIF respectively. Figure 3 shows the
unnormalised edge weight as a function of water fraction difference for both
LC-aware specifications at their tuned $\alpha$ values — both curves are
nearly flat, visually explaining the insensitivity. This insensitivity is
expected once the water proportion covariate is included: the covariate
already captures most of the land--water contrast, leaving little for the
neighbourhood structure to explain.

```{r nb_visualization, fig.cap="Left: unnormalised edge weight as a function of water fraction difference $|w_i - w_j|$ for the two LC-aware specifications at their CV-tuned $\\alpha$ values, alongside the flat queen baseline. Both LC-aware curves are nearly flat at the tuned $\\alpha$, explaining why results are insensitive to neighbourhood specification. Right: dominant landcover in the target grid, showing the land--water boundary where LC-aware weights would differ from queen.", fig.height=3.5, fig.width=9}

# Retrieve CV-tuned alpha values
alpha_eps   <- gf(r$nb_lc_eps,   "alpha_opt")
alpha_alpha <- gf(r$nb_lc_alpha, "alpha_opt")

# Weight as function of |w_i - w_j| (water fraction difference)
# lc_eps:   w_ij_unnorm = 1[same_type] + exp(-alpha)
#           same_type = 1 if both land or both water
#           For a continuous water diff d: same_type ~ 1 if d < 0.5 threshold
#           But actual implementation uses binary is_water (>0.5)
#           So: weight = 1 + exp(-alpha) if same type, exp(-alpha) if different
# lc_alpha: w_ij_unnorm = exp(-alpha * |w_i - w_j|) + 1e-6

d_seq <- seq(0, 1, by = 0.01)

# lc_eps: binary -- d >= 0.5 threshold means different type (land vs water)
eps_val <- exp(-alpha_eps)
w_eps   <- ifelse(d_seq < 0.5, 1 + eps_val, eps_val)

# lc_alpha: continuous
w_alpha_val <- exp(-alpha_alpha * d_seq) + 1e-6

# Queen: flat at 1 (before normalisation, all edges equal)
w_queen <- rep(1, length(d_seq))

df_w <- data.frame(
  d       = rep(d_seq, 3),
  weight  = c(w_queen, w_eps, w_alpha_val),
  Spec    = rep(c(
    "Queen (baseline)",
    sprintf("LC-eps (alpha=%.2f)", alpha_eps),
    sprintf("LC-alpha (alpha=%.2f)", alpha_alpha)
  ), each = length(d_seq))
)
df_w$Spec <- factor(df_w$Spec, levels = unique(df_w$Spec))

p_weights <- ggplot(df_w, aes(x = d, y = weight, color = Spec,
                               linetype = Spec)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey50") +
  annotate("text", x = 0.52, y = max(df_w$weight) * 0.95,
           label = "Land/water\nthreshold", hjust = 0, size = 2.5,
           color = "grey50") +
  scale_color_manual(
    values = c("steelblue", "orange", "tomato"),
    name = NULL) +
  scale_linetype_manual(
    values = c("solid", "dashed", "dotdash"),
    name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression("|"*w[i] - w[j]*"|"~"(water fraction difference)"),
       y = "Unnormalised edge weight",
       title = "Edge weight vs water fraction difference") +
  theme_minimal(base_size = 9) +
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10, hjust = 0.5))

# --- Right panel: landcover map ---
tg <- goebel2026::target_grid
dc <- tg$dominant_class
tg$dom_class <- dplyr::case_when(
  dc %in% c(1, 8)       ~ "Built/Bare",
  dc %in% c(6, 7)       ~ "Water",
  dc %in% c(2, 3, 4, 5) ~ "Green",
  TRUE                  ~ NA_character_
)

p_lc <- ggplot(tg[!is.na(tg$dom_class), ]) +
  geom_sf(aes(fill = dom_class), color = NA) +
  scale_fill_manual(
    values = c("Built/Bare" = "#c4281b",
               "Water"      = "#429ae4",
               "Green"      = "#397e48"),
    name = "Landcover", na.value = "grey90") +
  labs(title = "Dominant landcover (target grid)") +
  theme_void(base_size = 9) +
  theme(legend.position = "right",
        plot.title      = element_text(size = 10, hjust = 0.5))

p_weights | p_lc
```

```{r nb_albedo}
nbs_alb <- list(r$water_rsr, r$nb_rook, r$nb_lc_eps, r$nb_lc_alpha)

alpha_eps_hat   <- fmt3(gf(r$nb_lc_eps,   "alpha_opt"))
alpha_alpha_hat <- fmt3(gf(r$nb_lc_alpha, "alpha_opt"))

tab_nb_alb <- data.frame(
  Neighbourhood = c(
    "Queen (baseline)",
    "Rook",
    paste0("LC-eps (alpha=", alpha_eps_hat, ")"),
    paste0("LC-alpha (alpha=", alpha_alpha_hat, ")")
  ),
  rho      = sapply(nbs_alb, get_rho),
  phi      = sapply(nbs_alb, gf, "phi"),
  RMSE     = sapply(nbs_alb, gf, "rmse"),
  R2       = sapply(nbs_alb, gf, "r2"),
  Coverage = sapply(nbs_alb, gf, "coverage_95_obs"),
  beta_w   = sapply(nbs_alb, gf, "beta_hat", idx = 2)
)

tab_nb_alb %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "Albedo predictive performance across neighbourhood weight specifications. All models include the water proportion covariate and RSR constraint at $\\rho=1$.",
    label    = "tab:nb_albedo",
    col.names = c("Neighbourhood", "$\\hat\\rho$", "$\\hat\\phi$",
                  "RMSE", "$R^2$", "Coverage (95\\%)",
                  "$\\hat\\beta_\\text{water}$"),

  )
```

```{r nb_sif}
nbs_sif <- list(r$sif_rsr_rho1, r$sif_nb_rook, r$sif_nb_lc_eps, r$sif_nb_lc_alpha)

alpha_sif_eps_hat   <- fmt3(gf(r$sif_nb_lc_eps,   "alpha_opt"))
alpha_sif_alpha_hat <- fmt3(gf(r$sif_nb_lc_alpha, "alpha_opt"))

tab_nb_sif <- data.frame(
  Neighbourhood = c(
    "Queen (baseline)",
    "Rook",
    paste0("LC-eps (alpha=", alpha_sif_eps_hat, ")"),
    paste0("LC-alpha (alpha=", alpha_sif_alpha_hat, ")")
  ),
  rho     = sapply(nbs_sif, get_rho),
  phi     = sapply(nbs_sif, gf, "phi"),
  beta_w  = sapply(nbs_sif, gf, "beta_hat", idx = 2),
  mean_SE = sapply(nbs_sif, function(x)
              mean(x$posterior_se, na.rm = TRUE))
)

tab_nb_sif %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "SIF posterior summaries across neighbourhood weight specifications. All models include the water proportion covariate, RSR constraint, and $R^{-1}$ weighting.",
    label    = "tab:nb_sif",
    col.names = c("Neighbourhood", "$\\hat\\rho$", "$\\hat\\phi$",
                  "$\\hat\\beta_\\text{water}$", "Mean SE"),

  )
```

\newpage

# Restricted Spatial Regression {#sec:rsr}

Spatial confounding arises when both the covariate and the spatial random
effect can explain the same variation in the response, leading to unstable
and uninterpretable covariate estimates. We address this through restricted
spatial regression (RSR), which orthogonalises the spatial residual against
the covariate space, ensuring that $\hat\beta$ is identified by
the covariate variation alone.

The table below illustrates the importance of RSR for SIF. Without
any spatial random effect (OLS) or with an unconstrained spatial effect
(naive), the estimated water coefficient is
$\hat\beta_\text{water} = `r fmt4(gf(r$sif_naive,"beta_hat",idx=2))`$, which is
positive — implying that water pixels have higher SIF than land pixels, a
physically implausible result. Applying RSR reverses the sign:
$\hat\beta_\text{water} = `r fmt4(gf(r$sif,"beta_hat",idx=2))`$, consistent with
water having lower photosynthetic activity than vegetated land. This sign
reversal occurs because without RSR the spatial random effect absorbs the
land--water contrast, leaving the covariate to explain residual variation
that is spuriously positive.

For albedo, where the ground truth is known, RSR also improves $R^2$ and
coverage (Table 2 in the main paper), confirming that the orthogonalisation
produces better-calibrated uncertainty estimates as well as more
interpretable coefficients.

```{r rsr_sif}
rsrs <- list(r$sif_naive, r$sif_rsr_rho1, r$sif)

tab_rsr <- data.frame(
  Model    = c("Naive (no RSR)", "RSR, $\\rho=1$",
               "RSR, $\\rho$ tuned (canonical)"),
  rho      = sapply(rsrs, get_rho),
  phi      = sapply(rsrs, gf, "phi"),
  beta_int = sapply(rsrs, gf, "beta_hat", idx = 1),
  beta_w   = sapply(rsrs, gf, "beta_hat", idx = 2),
  mean_SE  = sapply(rsrs, function(x) mean(x$posterior_se, na.rm = TRUE))
)

tab_rsr %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "SIF posterior summaries with and without RSR. The naive model (no RSR) estimates a positive water coefficient, which is physically implausible. RSR corrects the sign and produces stable estimates across $\\rho$ specifications.",
    label    = "tab:rsr_sif",
    col.names = c("Model", "$\\hat\\rho$", "$\\hat\\phi$",
                  "$\\hat\\beta_0$", "$\\hat\\beta_\\text{water}$",
                  "Mean SE"),

  )
```

\newpage

# Forward Operator Sensitivity {#sec:forward_operator}

The aggregation matrix $A$ discretises the continuous change-of-support
integral by replacing each intersection integral with a weighted average
over the intersecting latent cell. The default specification uses uniform
weights within each intersection ($g_i \equiv 1$). Here we assess sensitivity
to this choice by considering a family of super-Gaussian sensor response
functions parameterised by a width $\tau$.

**Sensor response function.** For sounding $i$ with centroid $c_i$ and
affine map $M_i^{-1}$ transforming footprint coordinates to a canonical
unit square, the response function is:
$$g_i(s) = \exp\!\left(-\frac{\|M_i^{-1}(s - c_i)\|^2}{2\tau^2}\right)$$
As $\tau \to \infty$ the response approaches uniform ($g_i \equiv 1$,
the default). As $\tau \to 0$ the response concentrates at the footprint
centre, approaching the centroid approximation in the limit.
The aggregation weight for cell $j$ under sounding $i$ is then
$a_{ij} \propto \int_{D_i^{(o)} \cap D_j^{(\ell)}} g_i(s)\,ds$,
normalised so that $\sum_j a_{ij} = 1$.

Figure 4 shows $g_i$ evaluated on a representative rectangular OCO-2
footprint for each $\tau$ value considered. The five tau values tested
span from nearly uniform (tau=0.5) to sharply peaked (tau=0.01).

```{r gA_viz, fig.cap="Super-Gaussian sensor response $g_i(s)$ evaluated on a representative OCO-2 footprint for each $\\tau$ value. Uniform ($\\tau = \\infty$, top left) assigns equal weight everywhere; smaller $\\tau$ concentrates weight at the footprint centre. The centroid approximation (bottom right) is the limiting case.", fig.height=4, fig.width=9}

# Evaluate g on a canonical unit rectangle representing a footprint
# Use normalised coordinates: footprint spans [-1,1] x [-1,1]
# M_inv = identity (canonical), so ||M^{-1}(s-c)||^2 = x^2+y^2

nx <- 60; ny <- 40
xs <- seq(-1, 1, length.out = nx)
ys <- seq(-1, 1, length.out = ny)
grid_uv <- expand.grid(u = xs, v = ys)

tau_vals_viz <- c(Inf, 0.5, 1/3, 0.2, 0.1, 0.01)
tau_labs     <- c("Uniform (tau=inf)", "tau=0.5", "tau=0.33",
                  "tau=0.2",           "tau=0.1", "tau=0.01 (centroid)")

df_g <- do.call(rbind, lapply(seq_along(tau_vals_viz), function(k) {
  tau <- tau_vals_viz[k]
  if (is.infinite(tau)) {
    g <- rep(1, nrow(grid_uv))
  } else {
    g <- exp(-(grid_uv$u^2 + grid_uv$v^2) / (2 * tau^2))
  }
  data.frame(u     = grid_uv$u,
             v     = grid_uv$v,
             g     = g / max(g),   # normalise to [0,1] for display
             label = tau_labs[k])
}))
df_g$label <- factor(df_g$label, levels = tau_labs)

ggplot(df_g, aes(x = u, y = v, fill = g)) +
  geom_tile() +
  facet_wrap(~ label, nrow = 2) +
  scale_fill_gradient(low = "white", high = "tomato",
                      name = expression(g[i](s)~"(normalised)")) +
  coord_equal() +
  labs(x = "Footprint u-coordinate",
       y = "Footprint v-coordinate") +
  theme_minimal(base_size = 9) +
  theme(legend.position  = "bottom",
        panel.grid       = element_blank(),
        strip.text       = element_text(size = 8))
```

## Albedo validation

Because the albedo data has a known generating field, we can directly assess
the impact of forward operator misspecification on RMSE and coverage.
Figure 2 and the table below show results for
albedo as $\tau$ decreases from 0.5 (broad response, close to uniform) to
the centroid limit. Coverage degrades monotonically from
`r pct(gf(r$water_rsr_cv,"coverage_95_obs"))` at the uniform baseline to
`r pct(gf(r$alb_gA_cent,"coverage_95_obs"))` at the centroid limit, while RMSE
increases by `r fmt4(gf(r$alb_gA_cent,"rmse") - gf(r$water_rsr_cv,"rmse"))` absolute
units. This confirms that forward operator misspecification primarily damages
uncertainty quantification rather than point prediction accuracy.

```{r gA_albedo, fig.cap="Coverage (95\\%) and RMSE as functions of the sensor response width $\\tau$ for the albedo validation experiment. The uniform baseline ($\\tau = \\infty$, leftmost point) achieves near-nominal coverage; coverage degrades monotonically as $\\tau \\to 0$ (centroid limit, rightmost point).", fig.height=3.5, fig.width=7}

tau_vals <- c(Inf, 0.5, 1/3, 0.2, 0.1, 0)
cov_vals <- c(gf(r$water_rsr_cv, "coverage_95_obs"),
              gf(r$alb_gA_tau05,  "coverage_95_obs"),
              gf(r$alb_gA_tau033, "coverage_95_obs"),
              gf(r$alb_gA_tau02,  "coverage_95_obs"),
              gf(r$alb_gA_tau01,  "coverage_95_obs"),
              gf(r$alb_gA_cent,   "coverage_95_obs"))
rmse_vals <- c(gf(r$water_rsr_cv, "rmse"),
               gf(r$alb_gA_tau05,  "rmse"),
               gf(r$alb_gA_tau033, "rmse"),
               gf(r$alb_gA_tau02,  "rmse"),
               gf(r$alb_gA_tau01,  "rmse"),
               gf(r$alb_gA_cent,   "rmse"))

tau_labels <- c("Uniform\n(baseline)", "0.5", "0.33", "0.2", "0.1", "Centroid")
df_gA <- data.frame(
  tau       = factor(tau_labels, levels = tau_labels),
  coverage  = cov_vals,
  rmse      = rmse_vals
)

p_cov <- ggplot(df_gA, aes(x = tau, y = coverage, group = 1)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "gray50") +
  geom_line(color = "steelblue", linewidth = 0.9) +
  geom_point(color = "steelblue", size = 2.5) +
  scale_y_continuous(limits = c(0.5, 1.0),
                     labels = scales::percent_format(accuracy=1)) +
  labs(x = expression(tau), y = "Coverage (95%)") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_rmse <- ggplot(df_gA, aes(x = tau, y = rmse, group = 1)) +
  geom_line(color = "tomato", linewidth = 0.9) +
  geom_point(color = "tomato", size = 2.5) +
  labs(x = expression(tau), y = "RMSE") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank())

p_cov | p_rmse
```

```{r gA_albedo_table}
gAs_alb <- list(r$water_rsr_cv, r$alb_gA_tau05, r$alb_gA_tau033,
                r$alb_gA_tau02, r$alb_gA_tau01, r$alb_gA_cent)

tab_gA_alb <- data.frame(
  tau      = c("Uniform (baseline)", "0.5", "0.33", "0.2", "0.1", "Centroid"),
  rho      = sapply(gAs_alb, get_rho),
  phi      = sapply(gAs_alb, gf, "phi"),
  RMSE     = rmse_vals,
  R2       = sapply(gAs_alb, gf, "r2"),
  Coverage = cov_vals
)

tab_gA_alb %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "Albedo predictive performance under forward operator misspecification. $\\tau$ controls the width of the super-Gaussian sensor response function; $\\tau \\to 0$ corresponds to the centroid approximation. Coverage degrades monotonically while RMSE increases modestly.",
    label    = "tab:gA_albedo",
    col.names = c("$\\tau$", "$\\hat\\rho$", "$\\hat\\phi$",
                  "RMSE", "$R^2$", "Coverage (95\\%)"),

  )
```

\newpage

## SIF sensitivity

For SIF, ground truth is unavailable, so we assess sensitivity through
changes in posterior mean and uncertainty. The table below shows
that as $\tau$ decreases, $\hat\rho$ increases toward 1 and $\hat\phi$
decreases, reflecting the model compensating for a more peaked forward
operator by increasing spatial smoothness. The water coefficient
$\hat\beta_\text{water}$ remains stable across all specifications
(range: `r fmt4(min(sapply(list(r$sif, r$sif_gA_tau05, r$sif_gA_tau033, r$sif_gA_tau02, r$sif_gA_tau01, r$sif_gA_tau001), gf, "beta_hat", idx=2)))`
to
`r fmt4(max(sapply(list(r$sif, r$sif_gA_tau05, r$sif_gA_tau033, r$sif_gA_tau02, r$sif_gA_tau01, r$sif_gA_tau001), gf, "beta_hat", idx=2)))`),
suggesting that the posterior mean is robust to the exact form of the sensor
response function even as uncertainty estimates vary.

```{r gA_sif}
gAs_sif <- list(r$sif, r$sif_gA_tau05, r$sif_gA_tau033,
                r$sif_gA_tau02, r$sif_gA_tau01, r$sif_gA_tau001)

tab_gA_sif <- data.frame(
  tau     = c("Uniform (baseline)", "0.5", "0.33", "0.2", "0.1", "0.01"),
  rho     = sapply(gAs_sif, get_rho),
  phi     = sapply(gAs_sif, gf, "phi"),
  beta_w  = sapply(gAs_sif, gf, "beta_hat", idx = 2),
  mean_SE = sapply(gAs_sif, function(x) mean(x$posterior_se, na.rm = TRUE))
)

tab_gA_sif %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "SIF posterior summaries under forward operator misspecification. The water coefficient is stable across all specifications; mean posterior SE decreases as $\\tau \\to 0$ reflecting increased smoothing.",
    label    = "tab:gA_sif",
    col.names = c("$\\tau$", "$\\hat\\rho$", "$\\hat\\phi$",
                  "$\\hat\\beta_\\text{water}$", "Mean SE"),

  )
```

## Kriging comparison

The table below compares the Bayesian downscaling approach to
universal kriging on the albedo validation experiment. Both methods are
applied to the same semi-synthetic observations; kriging uses an empirically
fitted exponential variogram. The Bayesian approach achieves substantially
better coverage (`r pct(gf(r$water_rsr_cv,"coverage_95_obs"))` vs.\
`r pct(gf(r$krig_alb,"coverage_95_obs"))`) with similar RMSE, demonstrating
that the change-of-support formulation produces well-calibrated uncertainty
estimates that kriging — which treats observations as point-indexed — does not.

```{r kriging}
krigs <- list(r$water_rsr_cv, r$alb_gA_cent, r$krig_alb)

tab_krig <- data.frame(
  Method   = c("Bayesian (uniform $A$, baseline)",
               "Bayesian (centroid $A$)",
               "Universal kriging"),
  RMSE     = sapply(krigs, gf, "rmse"),
  R2       = sapply(krigs, gf, "r2"),
  Coverage = sapply(krigs, gf, "coverage_95_obs")
)

tab_krig %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "Comparison of Bayesian downscaling and universal kriging on the albedo validation experiment. The centroid approximation row shows performance when the aggregation matrix uses footprint centroids rather than intersection areas. Kriging coverage is substantially below nominal; the Bayesian approach with correct change-of-support aggregation achieves near-nominal coverage.",
    label    = "tab:kriging",
    col.names = c("Method", "RMSE", "$R^2$", "Coverage (95\\%)"),

  )
```

\newpage

# Rho Sensitivity {#sec:rho}

The SAR spatial range parameter $\rho$ controls the decay of spatial
correlation: $\rho = 1$ corresponds to an intrinsic (improper) prior with
no decay, while smaller values introduce exponential decay in correlations
with distance. In the main analysis we tune $\rho$ jointly with $\phi$ via
cross-validation.

The table below shows albedo results at $\rho = 1$ (intrinsic)
and at the CV-tuned value $\hat\rho =
`r fmt3(get_rho(r$water_rsr_cv))`$. The differences in RMSE, $R^2$, and
coverage are small (less than 1 percentage point in coverage), confirming
that the intrinsic prior is a reasonable default. We nevertheless prefer
the tuned $\rho$ in the main analysis as it provides marginally better
coverage calibration. For SIF, the canonical model uses
$\hat\rho = `r fmt3(get_rho(r$sif))`$; fixing $\rho = 0.95$ gives
essentially identical posterior means and $\hat\beta_\text{water}$ values
(Table 8 in the main results summary).

```{r rho_sensitivity}
rhos <- list(r$water_rsr, r$water_rsr_cv, r$water, r$water_rho_cv)

tab_rho <- data.frame(
  Model    = c("Water + RSR, $\\rho=1$",
               "Water + RSR, $\\rho$ tuned",
               "Water (no RSR), $\\rho=1$",
               "Water (no RSR), $\\rho$ tuned"),
  rho      = sapply(rhos, get_rho),
  phi      = sapply(rhos, gf, "phi"),
  RMSE     = sapply(rhos, gf, "rmse"),
  R2       = sapply(rhos, gf, "r2"),
  Coverage = sapply(rhos, gf, "coverage_95_obs")
)

tab_rho %>%
  kable(
    format = "simple",

    digits   = 4,
    caption  = "Albedo predictive performance at fixed $\\rho=1$ versus CV-tuned $\\rho$, with and without RSR. Differences are small in all metrics, confirming that the intrinsic prior ($\\rho=1$) is a reasonable default.",
    label    = "tab:rho_sensitivity",
    col.names = c("Model", "$\\hat\\rho$", "$\\hat\\phi$",
                  "RMSE", "$R^2$", "Coverage (95\\%)"),

  )
```
