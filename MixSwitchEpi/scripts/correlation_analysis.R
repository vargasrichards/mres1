#!/usr/bin/env Rscript
# Correlation, PCA and MDS analysis across mix-switch space
# This script loads output data from a parameter scan and computes how correlated
# epidemic metrics are, finds descriptor dimensionality via PCA, and clusters via MDS.

library(dplyr)
library(ggplot2)
library(corrplot)
library(tidyr)

# Load package
devtools::load_all()

# We can either load existing scan data or run a quick scan. Let's run a quick fine grid.
message("Running a high-resolution scan for correlation & PCA analysis...")

base_pars <- MixSwitchEpi::make_parms(
  n_activity = 5,
  epsilon = 0,
  switch_rate = 0,
  initial_condition_rule = "uniform",
  infected_fraction = 1e-3,
  activity_scheme = "gamma",
  gamma_dist = list(mean = 4.27, variance = 10.63),
  mixing_model = "garnett_mixing",
  switching_model = "simple_switching",
  pathogen_string = "sars-cov-2",
  target_r0 = 2
)

res_vals <- sort(unique(c(Inf, seq(1, 50, length.out = 20))))
switch_rates <- 1 / res_vals

# Simulate across both space and different R0s for robustness
scan_data <- list()
for (tr0 in c(1.5, 2.0, 3.0)) {
  sub_pars <- base_pars
  sub_pars$target_r0 <- tr0

  scan_results <- MixSwitchEpi::scan_mixswitch_framework(
    base_pars = sub_pars,
    n_epsilon = 20,
    n_switch = length(switch_rates),
    epsilon_vals = seq(0, 1, length.out = 20),
    switch_vals = switch_rates,
    classed_detail = FALSE,
    switch_model = "simple_switching",
    mix_model = "garnett_mixing",
    initial_condition_rule = "uniform",
    frac_infected_initial = 1e-3,
    calibration_mode = "switch"
  )

  df_sub <- dplyr::bind_rows(scan_results) |>
    dplyr::filter(class == -1) |>
    dplyr::mutate(target_r0 = tr0)

  scan_data[[as.character(tr0)]] <- df_sub
}

df <- dplyr::bind_rows(scan_data)

message("\nAnalysing metric dimensionality via correlation, PCA, and MDS...")
analysis <- MixSwitchEpi::analyse_metric_dimensionality(
  scan_df = df,
  metrics = c("hit", "final_size", "end_time", "r0", "psI", "ptI")
)

message("Results of Spearman Rank Correlation:")
print(analysis$correlation_spearman)

message("\nVariance explained by Principal Components:")
print(analysis$variance_explained)

out_dir <- "output/metric_dimensionality"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Plot raw correlation matrix
pdf(file.path(out_dir, "metric_correlation_plot.pdf"), width = 8, height = 8)
corrplot::corrplot(
  analysis$correlation_spearman,
  method = "ellipse",
  type = "upper",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  title = "Metric Correlation across Mix-Switch Space & varying R0",
  mar = c(0, 0, 2, 0)
)
dev.off()

# 2. Plot PCA Biplot
p_pca <- MixSwitchEpi::plot_metric_pca(analysis, color_var = "epsilon")
ggplot2::ggsave(file.path(out_dir, "metric_pca_biplot.pdf"), p_pca, width = 9, height = 7)

# 3. Plot MDS Space
p_mds <- MixSwitchEpi::plot_metric_mds(analysis)
ggplot2::ggsave(file.path(out_dir, "metric_mds_clustering.pdf"), p_mds, width = 8, height = 6)

# 4. (Optional) Pairwise grid
p_pairs <- MixSwitchEpi::plot_metric_pairs(analysis, color_var = "epsilon")
if (inherits(p_pairs, "patchwork") || inherits(p_pairs, "ggplot")) {
  ggplot2::ggsave(file.path(out_dir, "metric_pairwise_grid.pdf"), p_pairs, width = 15, height = 12)
}

message(sprintf("\nAnalysis complete! Visualisations saved to '%s'", out_dir))
