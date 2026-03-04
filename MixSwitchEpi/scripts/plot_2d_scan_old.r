library(readr)
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)

df <- read_csv("output/comparison/hit_final_size_end_time_r0_psI_ptI_scan_processed_long.csv")

# split metric and diff
metric_df <- df %>% filter(panel == "metric")
diff_df <- df %>% filter(panel == "diff")

# compute readable y breaks (same approach used in the function)
y_breaks <- pretty(range(df$res_time_plot, na.rm = TRUE), n = 6)
if (any(is.infinite(df$residence_time))) y_breaks <- sort(unique(c(y_breaks, max(df$res_time_plot, na.rm = TRUE))))

p_metrics <- ggplot(metric_df, aes(x = epsilon, y = res_time_plot, fill = value)) +
  geom_tile() +
  facet_wrap(~variable, nrow = 1) +
  scale_y_continuous(breaks = y_breaks) +
  scale_fill_viridis(option = "mako", na.value = "grey80") +
  theme_minimal()

lim_diff <- max(abs(diff_df$value), na.rm = TRUE)
p_diffs <- ggplot(diff_df, aes(x = epsilon, y = res_time_plot, fill = value)) +
  geom_tile() +
  facet_wrap(~variable, nrow = 1) +
  scale_y_continuous(breaks = y_breaks) +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b",
    midpoint = 0, limits = if (lim_diff > 0) c(-lim_diff, lim_diff) else NULL,
    na.value = "grey80"
  ) +
  theme_minimal()

(p_metrics / p_diffs) + plot_layout(heights = c(1, 1))
