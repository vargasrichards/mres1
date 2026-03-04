# 2D comparison heatmaps: epsilon × expected residence time (1/switch_rate)
#
# Scan output must have: epsilon, switch_rate, and response metric columns.
# Y-axis is 1/switch_rate (expected residence time). switch_rate==0 → Inf.


#' Make homogeneous parameters
#' @param base_params base parameters with potential heterogeneity
#' @export
make_homogeneous_parms <- function(base_params) {
  mean_activity <- sum(base_params$activity_scores * base_params$class_sizes) /
    base_params$n_activity
  homogeneous_parms <- base_params
  homogeneous_parms$act_vector <- mean_activity
  homogeneous_parms
}

#' Add y_plot = 1/switch_rate for plotting.
#' Inf (switch_rate==0) → finite sentinel above max finite value.
#' @keywords internal
.add_y_plot <- function(scan_df) {
  if (!"switch_rate" %in% names(scan_df)) {
    stop("scan_df must have switch_rate column")
  }
  scan_df$res_time <- 1 / scan_df$switch_rate
  finite_ert <- scan_df$res_time[is.finite(scan_df$res_time)]
  if (length(finite_ert) == 0) stop("no finite expected residence times")
  max_finite <- max(finite_ert, na.rm = TRUE)
  inf_plot_val <- max_finite * 1.15
  scan_df$y_plot <- ifelse(is.finite(scan_df$res_time),
    scan_df$res_time, inf_plot_val
  )
  attr(scan_df, "inf_plot_val") <- inf_plot_val
  scan_df
}

#' Build y-axis scale with Inf labeled as ∞
#' @keywords internal
.y_scale <- function(scan_df, n_breaks = 6) {
  inf_plot_val <- attr(scan_df, "inf_plot_val")
  finite_y <- sort(unique(scan_df$y_plot[is.finite(scan_df$res_time)]))
  if (length(finite_y) == 0) {
    return(ggplot2::scale_y_continuous())
  }

  brks <- pretty(finite_y, n = n_breaks)
  brks <- brks[brks >= min(finite_y) & brks <= max(finite_y)]
  if (!is.null(inf_plot_val) && any(!is.finite(scan_df$res_time))) {
    brks <- sort(unique(c(brks, inf_plot_val)))
  }
  labs <- ifelse(brks == inf_plot_val, expression(infinity),
    as.character(round(brks, 2))
  )
  ggplot2::scale_y_continuous(breaks = brks, labels = labs)
}

#' Multi-outcome 2-D comparison heatmaps
#'
#' @param base_params list from \code{make_parms()}
#' @param file scan output
#' @param response_vars character vector of metric column names
#' @param ref_epsilon,ref_switch reference point
#' @param ref_r0 optional R0
#' @param n_activity the (integer) number of activity classes in the model
#' @param output_dir save location
#' @param save write files?
#' @param width,height size per panel
#' @export
multi_outcome_comparison <- function(base_params, file,
                                     response_vars,
                                     ref_switch = NULL,
                                     ref_r0, n_activity = 5,
                                     output_dir = "output/comparison",
                                     save = TRUE, width = 6, height = 4,
                                     show_diff = TRUE,
                                     slice_res_times = NULL) {
  scan_df <- file
  eps_vals <- sort(unique(scan_df$epsilon))
  hom_metrics <- homogeneous_metrics(r0 = ref_r0)

  for (rv in response_vars) {
    if (!rv %in% names(scan_df)) {
      warning(rv, " not in scan data")
      next
    }
    scan_df[[paste0(rv, "_metric")]] <- as.numeric(scan_df[[rv]])
    if (!is.null(hom_metrics) && rv %in% names(hom_metrics)) {
      ref_val <- hom_metrics[[rv]]
    } else {
      stop(sprintf("Could not find homogeneous reference value
        for metric '%s' from homogeneous_metrics()", rv))
    }
    scan_df[[paste0(rv, "_diff")]] <- scan_df[[paste0(rv, "_metric")]] - as.numeric(ref_val)
  }
  plot_df <- scan_df


  # infectious-period annotation (added once per full patchwork, not per tile)
  gamma_time <- if (!is.null(base_params$gamma)) 1 / base_params$gamma else NULL
  ann <- if (!is.null(gamma_time)) {
    # Math label for infectious period: 1/gamma == value (days)
    lab <- paste0("1/gamma == ", signif(gamma_time, 3))
    list(
      ggplot2::geom_hline(
        yintercept = gamma_time,
        linetype = "dashed",
        colour = "#008b8b",
        linewidth = 0.7
      ),
      ggplot2::annotate(
        "text",
        x = min(eps_vals), y = gamma_time,
        label = lab,
        hjust = 0, vjust = -0.7,
        parse = TRUE, size = 3,
        colour = "#004c4c"
      )
    )
  }

  # Function to generate a single tile plot
  # x_show, y_show: logical controlling axis label/tick visibility
  make_tile_plot <- function(data, x_col, y_col, fill_col, title,
                             x_lab = NULL, y_lab = NULL, limits = NULL,
                             palette_type = "viridis",
                             x_show = TRUE, y_show = TRUE) {
    p <- ggplot2::ggplot(data, ggplot2::aes(
      x = .data[[x_col]], y = .data[[y_col]],
      fill = .data[[fill_col]]
    )) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        x = if (x_show) x_lab else NULL,
        y = if (y_show) y_lab else NULL,
        title = title,
        fill = if (palette_type == "divergent") "diff" else title
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        legend.position = "right",
        # Control axis labels/ticks
        axis.text.x = if (x_show) ggplot2::element_text(size = 13, angle = 45, hjust = 1) else ggplot2::element_blank(),
        axis.text.y = if (y_show) ggplot2::element_text(size = 13) else ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(), title = ggplot2::element_text(size = 14),
        # Remove whitespace padding around the plot area
        plot.margin = ggplot2::margin(2, 2, 2, 2),
        panel.spacing = ggplot2::unit(0, "lines")
      )

    # sparsify ticks if showing
    if (x_show) {
      levs <- levels(data[[x_col]])
      # show fewer breaks if many
      if (length(levs) > 10) {
        idx <- unique(round(seq(1, length(levs), length.out = 10)))
        p <- p + ggplot2::scale_x_discrete(breaks = levs[idx], expand = c(0, 0))
      }
    }
    if (y_show) {
      levs <- levels(data[[y_col]])
      if (length(levs) > 10) {
        idx <- unique(round(seq(1, length(levs), length.out = 10)))
        p <- p + ggplot2::scale_y_discrete(breaks = levs[idx], expand = c(0, 0))
      }
    }

    if (palette_type == "viridis") {
      # Use .data[["fill_col"]] inside aes, so we don't map fill_col as a variable name
      p <- p + viridis::scale_fill_viridis(option = "mako", na.value = "grey80")
    } else {
      p <- p + ggplot2::scale_fill_gradient2(
        low = "#2166ac", mid = "white", high = "#b2182b",
        midpoint = 0, limits = limits, na.value = "grey80"
      )
    }

    p
  }

  actual_plots <- list()
  diff_plots <- list()

  # Ensure plot_df is a data.frame to avoid data.table subsetting issues
  plot_df_safe <- as.data.frame(plot_df)

  # Epsilon: keep sorted unique values and show them as x factor levels.
  eps_levels_plot <- sort(unique(plot_df_safe$epsilon))
  plot_df_safe$x_fac <- factor(plot_df_safe$epsilon,
    levels = eps_levels_plot,
    labels = as.character(signif(eps_levels_plot, 3))
  )

  # Residence time (res_time): finite values in ascending order, with
  # an explicit "Inf" level appended if present.
  finite_res <- sort(unique(plot_df_safe$res_time[is.finite(plot_df_safe$res_time)]))
  y_levels_char <- as.character(finite_res)
  if (any(!is.finite(plot_df_safe$res_time))) y_levels_char <- c(y_levels_char, "Inf")
  # Labels: round finite numeric values for display; keep "Inf" literal.
  y_labels <- ifelse(y_levels_char == "Inf", "Inf", as.character(round(as.numeric(y_levels_char), 2)))

  plot_df_safe$y_fac <- factor(ifelse(is.finite(plot_df_safe$res_time), as.character(plot_df_safe$res_time), "Inf"),
    levels = y_levels_char, labels = y_labels
  )

  plots_combined_list <- list()
  n_vars <- length(response_vars)

  # If fewer than 3 vars, use that many columns; else cap at 3
  ncol_grid <- min(3, n_vars)

  for (i in seq_along(response_vars)) {
    rv <- response_vars[i]
    mc <- paste0(rv, "_metric")
    dc <- paste0(rv, "_diff")
    if (!mc %in% names(plot_df_safe)) next

    # Determine grid position visibility
    # Left-most column shows Y axis; diff panel shows X axis.
    col_idx <- (i - 1) %% ncol_grid + 1
    show_y <- (col_idx == 1)

    # Metric plot: hide x-axis when diff panel is shown below; show otherwise.
    metric_x_show <- if (isTRUE(show_diff)) FALSE else TRUE
    actual_plots[[rv]] <- make_tile_plot(plot_df_safe, "x_fac", "y_fac", mc,
      title = rv,
      x_lab = if (metric_x_show) expression(Assortativity ~ (epsilon)) else NULL,
      y_lab = if (show_y) "Expected residence time (days)" else NULL,
      palette_type = "viridis",
      x_show = metric_x_show,
      y_show = show_y
    )

    # Diff plot (Bottom of pair)
    lim_i <- max(abs(plot_df_safe[[dc]]), na.rm = TRUE)
    if (!is.finite(lim_i) || lim_i == 0) lim_i <- 1

    diff_plots[[rv]] <- make_tile_plot(plot_df_safe, "x_fac", "y_fac", dc,
      title = NULL,
      x_lab = expression(Assortativity ~ (epsilon)),
      y_lab = if (show_y) "Expected residence time (days)" else NULL,
      palette_type = "divergent",
      limits = c(-lim_i, lim_i),
      x_show = TRUE,
      y_show = show_y
    )

    # Combine them into one vertical unit
    if (show_diff) {
      unit <- patchwork::wrap_plots(actual_plots[[rv]], diff_plots[[rv]], ncol = 1, heights = c(1, 1))
    } else {
      unit <- actual_plots[[rv]]
    }
    plots_combined_list[[i]] <- unit
  }

  pw_heatmaps <- patchwork::wrap_plots(plots_combined_list, ncol = ncol_grid)

  # ---- 1D assortativity slice panels (optional) ----
  if (!is.null(slice_res_times) && length(slice_res_times) > 0) {
    # Match requested residence times to the closest available in the data
    avail_res <- sort(unique(plot_df_safe$res_time))
    match_res <- vapply(slice_res_times, function(rt) {
      if (is.infinite(rt)) {
        if (any(is.infinite(avail_res))) {
          return(Inf)
        }
        return(max(avail_res, na.rm = TRUE))
      }
      avail_res[which.min(abs(avail_res - rt))]
    }, numeric(1))
    match_res <- unique(match_res)

    slice_df <- plot_df_safe[plot_df_safe$res_time %in% match_res |
      (any(is.infinite(match_res)) & is.infinite(plot_df_safe$res_time)), ]
    slice_df$res_label <- factor(
      ifelse(is.finite(slice_df$res_time),
        paste0("Res = ", round(slice_df$res_time, 1), " d"),
        "No switching"
      ),
      levels = ifelse(is.finite(match_res),
        paste0("Res = ", round(match_res, 1), " d"),
        "No switching"
      )
    )

    # Build one line-plot per response variable, aligned to the heatmap columns
    slice_plots <- list()
    for (i in seq_along(response_vars)) {
      rv <- response_vars[i]
      mc <- paste0(rv, "_metric")
      if (!mc %in% names(slice_df)) next

      col_idx <- (i - 1) %% ncol_grid + 1
      show_y_slice <- (col_idx == 1)

      sp <- ggplot2::ggplot(slice_df, ggplot2::aes(
        x = epsilon,
        y = .data[[mc]],
        colour = res_label
      )) +
        ggplot2::geom_line(linewidth = 0.7) +
        ggplot2::geom_point(size = 0.8) +
        ggplot2::labs(
          x = expression(Assortativity ~ (epsilon)),
          y = if (show_y_slice) rv else NULL,
          colour = "Residence time"
        ) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          legend.position = if (i == length(response_vars)) "right" else "none",
          axis.text.y = if (show_y_slice) ggplot2::element_text(size = 8) else ggplot2::element_blank(),
          axis.title.y = if (show_y_slice) ggplot2::element_text(size = 10) else ggplot2::element_blank()
        )
      slice_plots[[i]] <- sp
    }

    pw_slices <- patchwork::wrap_plots(slice_plots, ncol = ncol_grid)
    # Stack: heatmaps on top (taller), 1D slices at bottom
    pw_combined <- patchwork::wrap_plots(pw_heatmaps, pw_slices,
      ncol = 1, heights = c(3, 1)
    )
  } else {
    pw_combined <- pw_heatmaps
  }

  # Add the infectious-period annotation once at the level of the full
  # patchwork object, to avoid repeated scale messages and duplicated lines.
  if (!is.null(ann)) {
    pw_combined <- pw_combined + ann
  }

  faceted_combined <- pw_combined

  # long tidy dataframe for CSV export
  # Reconstruct long format efficiently
  id_cols_vec <- c("epsilon", "switch_rate", "res_time")
  long_rows <- lapply(response_vars, function(rv) {
    mc <- paste0(rv, "_metric")
    dc <- paste0(rv, "_diff")
    if (!mc %in% names(plot_df_safe)) {
      return(NULL)
    }

    df_base <- plot_df_safe[, id_cols_vec, drop = FALSE]

    rbind(
      data.frame(df_base,
        variable = rv, panel = "metric",
        value = plot_df_safe[[mc]], stringsAsFactors = FALSE
      ),
      data.frame(df_base,
        variable = rv, panel = "diff",
        value = plot_df_safe[[dc]], stringsAsFactors = FALSE
      )
    )
  })
  long_df <- do.call(rbind, long_rows)
  long_df$panel <- factor(long_df$panel, levels = c("metric", "diff"))

  long_df$epsilon_fac <- factor(long_df$epsilon,
    levels = eps_levels_plot,
    labels = as.character(signif(eps_levels_plot, 3))
  )
  long_df$exp_res_fac <- factor(
    ifelse(is.finite(long_df$res_time),
      as.character(long_df$res_time), "Inf"
    ),
    levels = y_levels_char, labels = y_labels
  )


  if (save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    tag <- paste(response_vars, collapse = "_")

    # wide CSV/RDS
    wide_cols <- c("epsilon", "switch_rate", "res_time")
    for (rv in response_vars) wide_cols <- c(wide_cols, paste0(rv, c("_metric", "_diff")))
    wide_cols <- intersect(wide_cols, names(scan_df))
    wide_df <- if (data.table::is.data.table(scan_df)) {
      as.data.frame(scan_df)[, wide_cols, drop = FALSE]
    } else {
      scan_df[, wide_cols, drop = FALSE]
    }

    # Grid dimensions for saved PDF
    ncol_grid <- min(3, length(response_vars))
    nrow_grid <- ceiling(length(response_vars) / ncol_grid)
    panel_side <- width
    w_total <- min(panel_side * ncol_grid, 40)
    h_total <- min(panel_side * 2 * nrow_grid, 60)

    ggplot2::ggsave(file.path(output_dir, paste0(tag, "_multi_comparison.pdf")),
      pw_combined,
      width = w_total, height = h_total, device = grDevices::cairo_pdf
    )

    utils::write.csv(wide_df, file.path(output_dir, paste0(tag, "_scan_processed.csv")),
      row.names = FALSE, na = ""
    )

    utils::write.csv(long_df, file.path(output_dir, paste0(tag, "_scan_processed_long.csv")),
      row.names = FALSE, na = ""
    )

    message("Saved to ", output_dir)
  }

  list(
    patchwork = pw_combined, faceted = faceted_combined,
    wide_df = wide_df, long_df = long_df
  )
}


#' Example comparison scan and plot
#'
#' Scans
#'
#' @param res grid resolution
#' @param target_r0 the target r0 to calibrate against
#' @export
example_comp <- function(res = 50, target_r0, pathogen_string) {
  base_pars <- MixSwitchEpi::make_parms(
    n_activity = 5,
    epsilon = 0,
    switch_rate = 0,
    initial_condition_rule = "uniform",
    infected_fraction = 1e-3,
    activity_scheme = "gamma",
    gamma_dist = list(
      mean = 4.27,
      variance = 10.63
    ), # haslemere unique contacts
    mixing_model = "garnett_mixing",
    switching_model = "simple_switching",
    pathogen_string = pathogen_string, # can bne = "longer_recovery",
    target_r0 = target_r0
  )

  res_vals <- sort(unique(c(Inf, 1:50)))
  switch_rates <- 1 / res_vals

  scan <- MixSwitchEpi::scan_mixswitch_framework(
    base_pars = base_pars,
    n_epsilon = res,
    n_switch = length(switch_rates),
    epsilon_vals = seq(0, 1, length.out = res),
    switch_vals = switch_rates,
    classed_detail = FALSE, switch_model = "simple_switching",
    mix_model = "garnett_mixing", initial_condition_rule = "uniform",
    frac_infected_initial = 1e-3
  )

  scan_bound <- dplyr::bind_rows(scan)

  scan_bound_overall <- scan_bound |>
    dplyr::filter(class == -1) |>
    dplyr::mutate(res_time = 1 / .data$switch_rate)

  MixSwitchEpi::multi_outcome_comparison(base_pars,
    file = scan_bound_overall,
    ref_r0 = target_r0,
    response_vars = c(
      "hit",
      "final_size",
      "end_time",
      "r0",
      "psI",
      "ptI"
    ),
    save = TRUE,
    output_dir = "output/comparison",
    show_diff = FALSE
  )

  # Pull out specific residence times and plot 2D scan
  # (x=assortativity, y=outcome)
  target_res_times <- c(Inf, 1, 5, 20, 50)

  for (rtime in target_res_times) {
    switch_rate <- 1 / rtime
    specific_p <- base_pars
    specific_p$swch <- simple_switching(switch_rate = switch_rate, n_activity = n_activity)
  }

  plot_2d_data <- scan_bound_overall |>
    dplyr::filter(res_time %in% target_res_times) |>
    dplyr::mutate(res_time_label = factor(res_time,
      levels = target_res_times, labels = paste0("Res: ", target_res_times)
    ))

  p2d <- ggplot2::ggplot(plot_2d_data, ggplot2::aes(
    x = epsilon,
    y = hit,
    color = res_time_label
  )) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~res_time_label) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "HIT  Assortativity for selected residence times",
      x = "Assortativity",
      y = "Herd Immunity Threshold (HIT)",
      color = "Residence time"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 13)
    )

  ggplot2::ggsave("output/comparison/hit_2d_scan.pdf",
    p2d,
    width = 12,
    height = 8
  )

  final_size_2d <- ggplot2::ggplot(plot_2d_data, ggplot2::aes(
    x = "epsilon",
    y = "final_size",
    color = res_time_label
  )) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~res_time_label) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Final size against assortativity",
      x = "Assortativity",
      y = "Herd Immunity Threshold (HIT)",
      color = "Residence time (days)"
    ) +
    ggplot2::theme(
      text = ggplot2::element_text(size = 13),
      axis.text = ggplot2::element_text(size = 13)
    )

  ggplot2::ggsave("output/comparison/fs_2d_scan.pdf",
    p2d,
    width = 12,
    height = 8
  )

  # now we can look at the case where we repeatedly recalibrate.

  invisible(scan_bound_overall)
}


recalibration_scan <- function() {}
