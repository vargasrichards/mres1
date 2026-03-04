library(dplyr)

switch_rates <- seq(0, 1, length.out = 20)
epsilon_vals <- seq(0, 1, length.out = 20)

num_cores <- 10

out_dir <- "output/heatmaps_switch"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
message("Switch rates to be used:")
print(switch_rates)
base_pars <- MixSwitchEpi::default_setup()

results <- MixSwitchEpi::scan_2d_parameter_space(
  base_pars = base_pars,
  epsilon_vals = epsilon_vals,
  switch_vals = switch_rates,
  num_cores = num_cores,
  initial_condition_rule = "uniform"
)

if (is.null(results) || nrow(results) == 0) stop("scan returned no results")

results <- results %>%
  mutate(residence_time = ifelse(switch_rate == 0, Inf, 1 / switch_rate))

plot_data <- results %>% filter(class == -1) # here only examining the overall metrics

metrics_to_plot <- intersect(c("final_size", "hit", "duration", "ptEI"), names(plot_data))

time_infectious <- 1 / base_pars$gamma
switch_at_infectious <- base_pars$gamma

# ensure ordering for nicer plots
plot_data <- plot_data %>% arrange(desc(switch_rate))

plots_list <- list()
for (metric in metrics_to_plot) {
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = epsilon, y = switch_rate, fill = .data[[metric]])) +
    ggplot2::geom_tile(color = NA) +
    viridis::scale_fill_viridis(option = "mako", name = metric, na.value = "grey50") +
    ggplot2::labs(
      x = expression(Assortativity ~ (epsilon)), y = "Switch rate (per day)",
      title = paste0(metric, " vs epsilon and switch rate")
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::geom_hline(yintercept = switch_at_infectious, color = "red", linetype = "dashed", linewidth = 0.6) +
    ggplot2::annotate("text",
      x = max(epsilon_vals, na.rm = TRUE), y = switch_at_infectious,
      label = "1/γ", hjust = 1.05, vjust = -0.2, color = "red", size = 3
    )

  plots_list[[metric]] <- p

  fname_png <- file.path(out_dir, paste0(metric, "_switch_heatmap.png"))
  fname_svg <- file.path(out_dir, paste0(metric, "_switch_heatmap.svg"))
  ggplot2::ggsave(fname_png, p, width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(fname_svg, p, width = 8, height = 6)
  message("Saved ", fname_png)
}

# Create a multi-panel figure suitable for a paper (2x2)
if (requireNamespace("patchwork", quietly = TRUE)) {
  sel <- metrics_to_plot[seq_len(min(4, length(metrics_to_plot)))]
  multi <- patchwork::wrap_plots(plots_list[sel], ncol = 2) &
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  out_multi_png <- file.path(out_dir, "multipanel_switch_heatmaps.png")
  out_multi_svg <- file.path(out_dir, "multipanel_switch_heatmaps.svg")
  ggplot2::ggsave(out_multi_png, multi, width = 10, height = 8, dpi = 300)
  ggplot2::ggsave(out_multi_svg, multi, width = 10, height = 8)
  message("Saved multipanel figure to: ", out_multi_png)
} else {
  message("Package 'patchwork' not available; skipped multipanel figure.")
}

message("Done. Plots saved to: ", normalizePath(out_dir))
