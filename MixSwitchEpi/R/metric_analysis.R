#' Dimensionality and correlation analysis of epidemic metrics
#'
#' Evaluates how different epidemic outcome metrics correlate across
#' the mix-switch parameter space (and potentially different R0s),
#' and uses PCA to find a minimal set of descriptors.
#'
#' @param scan_df A data frame containing parameter scan results.
#' @param metrics A character vector of the metric column names to analyze.
#' @returns A list containing correlation matrices, PCA results, and summary statistics.
#' @export
analyse_metric_dimensionality <- function(scan_df,
                                          metrics = c("hit", "final_size", "end_time", "r0", "psI", "ptI")) {
  if (!all(metrics %in% names(scan_df))) {
    stop("Not all specified metrics are present in the scan_df.")
  }

  # Subset to numeric complete cases for the selected metrics
  df_sub <- as.data.frame(scan_df)[, metrics]
  valid_idx <- complete.cases(df_sub)
  df_clean <- df_sub[valid_idx, ]

  orig_data <- as.data.frame(scan_df)[valid_idx, ]

  # 1. Spearman Correlation (handles non-linear monotonic relationships)
  cor_mat <- cor(df_clean, method = "spearman")

  # 2. Pearson Correlation
  cor_pearson <- cor(df_clean, method = "pearson")

  # 3. PCA
  # Center and scale to unit variance since metrics have vastly different scales
  pca_res <- prcomp(df_clean, center = TRUE, scale. = TRUE)

  var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  cum_var_exp <- cumsum(var_exp)

  # Minimal descriptors info based on loadings
  loadings <- pca_res$rotation

  # 4. MDS (Multidimensional Scaling) on the columns (metrics) to see metric similarities
  # Distance between variables = 1 - abs(correlation)
  dist_metrics <- stats::as.dist(1 - abs(cor_mat))
  mds_metrics <- stats::cmdscale(dist_metrics, k = 2)

  list(
    correlation_spearman = cor_mat,
    correlation_pearson = cor_pearson,
    pca = pca_res,
    variance_explained = var_exp,
    cumulative_variance = cum_var_exp,
    loadings = loadings,
    mds_metrics = mds_metrics,
    data_clean = df_clean,
    orig_data = orig_data
  )
}

#' Plot Pairwise Scatterplots of Metrics
#'
#' @param analysis_res Result from analyse_metric_dimensionality
#' @param color_var A column name from the original data to color by (e.g. 'epsilon')
#' @export
plot_metric_pairs <- function(analysis_res, color_var = "epsilon") {
  df <- cbind(analysis_res$data_clean, analysis_res$orig_data[, color_var, drop = FALSE])

  # We construct a pairwise grid using ggplot2 and gridExtra/patchwork logic
  # Since GGally might not be a dependency, we build a simple melted version or use basic plots
  # Alternatively, just build a facet grid of all pairwise combinations

  metrics <- names(analysis_res$data_clean)
  n <- length(metrics)

  pairs_data <- expand.grid(x_var = metrics, y_var = metrics, stringsAsFactors = FALSE)
  pairs_data <- pairs_data[pairs_data$x_var != pairs_data$y_var, ]

  # Function to generate individual plots
  plot_list <- list()
  for (i in seq_len(nrow(pairs_data))) {
    xv <- pairs_data$x_var[i]
    yv <- pairs_data$y_var[i]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xv]], y = .data[[yv]], color = .data[[color_var]])) +
      ggplot2::geom_point(alpha = 0.6, size = 1) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::labs(x = xv, y = yv) +
      ggplot2::theme(legend.position = "none") +
      viridis::scale_color_viridis(option = "magma")

    plot_list[[paste0(xv, "_", yv)]] <- p
  }

  if (requireNamespace("patchwork", quietly = TRUE)) {
    patchwork::wrap_plots(plot_list, ncol = n - 1)
  } else {
    warning("patchwork package is required for full pairwise grid plotting.")
    plot_list
  }
}

#' Plot PCA biplot of metrics
#'
#' @param analysis_res Result from analyse_metric_dimensionality
#' @param color_var Variable to color the spatial points (e.g. 'epsilon', 'switch_rate')
#' @export
plot_metric_pca <- function(analysis_res, color_var = "epsilon") {
  pca <- analysis_res$pca
  var_exp <- analysis_res$variance_explained

  # Extract scores (points)
  pca_scores <- as.data.frame(pca$x)
  if (color_var %in% names(analysis_res$orig_data)) {
    pca_scores[[color_var]] <- analysis_res$orig_data[[color_var]]
  }

  # Extract loadings (arrows)
  loadings <- as.data.frame(pca$rotation)
  loadings$metric <- rownames(loadings)

  # Scale loadings for plotting to match scores range
  max_score <- max(abs(pca_scores[, 1:2]))
  max_loading <- max(abs(loadings[, 1:2]))
  scale_factor <- (max_score / max_loading) * 0.8
  loadings$PC1 <- loadings$PC1 * scale_factor
  loadings$PC2 <- loadings$PC2 * scale_factor

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = pca_scores,
      ggplot2::aes(x = PC1, y = PC2, color = .data[[color_var]]),
      alpha = 0.5, size = 2
    ) +
    viridis::scale_color_viridis(option = "mako") +
    ggplot2::geom_segment(
      data = loadings,
      ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "firebrick", linewidth = 1
    ) +
    ggplot2::geom_text(
      data = loadings,
      ggplot2::aes(x = PC1 * 1.1, y = PC2 * 1.1, label = metric),
      color = "firebrick", size = 5, fontface = "bold"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = "PCA Biplot of Epidemic Metrics",
      x = sprintf("PC1 (%.1f%% Variance)", var_exp[1] * 100),
      y = sprintf("PC2 (%.1f%% Variance)", var_exp[2] * 100)
    )

  p
}

#' Plot MDS of metrics
#'
#' Plots the multidimensional scaling projection of the metrics themselves
#' to observe clustering and find the minimal set of descriptors.
#' @param analysis_res Result from analyse_metric_dimensionality
#' @export
plot_metric_mds <- function(analysis_res) {
  mds_df <- as.data.frame(analysis_res$mds_metrics)
  colnames(mds_df) <- c("Dim1", "Dim2")
  mds_df$metric <- rownames(mds_df)

  ggplot2::ggplot(mds_df, ggplot2::aes(x = Dim1, y = Dim2, label = metric)) +
    ggplot2::geom_point(color = "darkblue", size = 3) +
    ggrepel::geom_text_repel(size = 5, fontface = "bold", color = "darkred") +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = "MDS Projection of Epidemic Metrics",
      subtitle = "Distance = 1 - abs(Spearman Correlation)",
      x = "MDS Dimension 1",
      y = "MDS Dimension 2"
    )
}
