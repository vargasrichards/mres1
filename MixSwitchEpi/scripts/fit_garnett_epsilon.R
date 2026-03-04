# Load package
devtools::load_all()

# Load necessary libraries for plotting and optimization
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(viridis)

# 1. Data Loading and Preparation
# Using the same activity scores and matrix as in compare_curves.R
raw_scores <- c(0.3146078, 1.1479220, 2.3361036, 4.2211382, 9.3466695)
contact_file <- "data/Haslemere.csv"
empirical_m <- as.matrix(utils::read.csv(contact_file, header = FALSE))
n_activity <- length(raw_scores)
class_sizes <- rep(1e6 / n_activity, n_activity)

# 2. Fitting Logic
# Objective function: Sum of Squared Errors (SSE) focused on whole matrix
sse_objective <- function(epsilon, target_m, acts, sizes) {
  pred_m <- garnett_mixing(
    n_activity = length(acts),
    epsilon_assort = epsilon,
    act_vector = acts,
    class_size_vector = sizes
  )
  # Fit to the whole matrix
  sum((target_m - pred_m)^2)
}

# Find the best epsilon in [0, 1]
fit <- stats::optimize(
  sse_objective,
  interval = c(0, 1),
  target_m = empirical_m,
  acts = raw_scores,
  sizes = class_sizes
)

best_epsilon <- fit$minimum
min_sse <- fit$objective

# Report Results
cat(sprintf("Fitted Epsilon (Whole Matrix): %.4f\n", best_epsilon))
cat(sprintf("Min SSE: %.6f\n", min_sse))

# 3. Visualization Preparation
pred_m_best <- garnett_mixing(n_activity, best_epsilon, raw_scores, class_sizes)
residuals_m <- empirical_m - pred_m_best

# Helper to convert matrix to long dataframe for ggplot
matrix_to_long <- function(m, name) {
  as.data.frame(m) %>%
    # Use mutate with row_number to avoid 1:n() warning
    mutate(row = row_number()) %>%
    pivot_longer(-row, names_to = "col", values_to = "value") %>%
    mutate(
      col = as.integer(gsub("V", "", col)),
      type = name
    )
}

df_emp <- matrix_to_long(empirical_m, "Empirical")
df_pred <- matrix_to_long(pred_m_best, "Fitted")
df_res <- matrix_to_long(residuals_m, "Residuals")

# Plotting function
# Note: Swapped axes per user request: Activity class (i) on X, Contacted class (j) on Y
plot_matrix <- function(df, title, palette = "mako", limits = NULL) {
  ggplot(df, aes(x = row, y = col, fill = value)) +
    geom_tile() +
    scale_fill_viridis(option = palette, name = "Value", limits = limits) +
    scale_y_continuous(breaks = 1:n_activity) +
    scale_x_continuous(breaks = 1:n_activity) +
    labs(title = title, x = "Activity class i", y = "Contacted class, j") +
    theme_minimal(base_size = 12) +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(size = 14),
      legend.position = "right", axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    )
}

# Define color limits for comparison
max_val <- max(c(as.numeric(empirical_m), as.numeric(pred_m_best)))
p_emp <- plot_matrix(df_emp, "A: Empirical Matrix", limits = c(0, max_val))
p_pred <- plot_matrix(df_pred, bquote("B: Fitted model (" * epsilon == .(round(best_epsilon, 3)) * ")"), limits = c(0, max_val))

# Residuals use a diverging palette
max_res <- max(abs(as.numeric(residuals_m)))
p_res <- ggplot(df_res, aes(x = row, y = col, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white",
    midpoint = 0, name = "Residuals",
    limits = c(-max_res, max_res)
  ) +
  scale_y_continuous(breaks = 1:n_activity) +
  scale_x_continuous(breaks = 1:n_activity) +
  labs(
    title = "C: Residuals (empirical - fitted)",
    x = "Activity class i",
    y = "Contacted class, j"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    aspect.ratio = 1,
    plot.title = element_text(size = 14),
    legend.position = "right", axis.text = element_text(size = 14),
    axis.title = element_text(size = 14)
  )

# Combine panels
combined_plot <- p_emp + p_pred + p_res +
  plot_layout(nrow = 1) +
  plot_annotation(
    title = "Parametric mixing fit to Haslemere contact matrix",
    theme = theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14)
    )
  )

output_dir <- "output/epsilon_fit"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "epsilon_fit_plots.pdf"),
  combined_plot,
  width = 14, height = 5,
  device = grDevices::cairo_pdf
)
ggsave(file.path(output_dir, "epsilon_fit_plots.png"),
  combined_plot,
  width = 14, height = 5, dpi = 300
)

cat(sprintf("Plots saved to %s\n", output_dir))
