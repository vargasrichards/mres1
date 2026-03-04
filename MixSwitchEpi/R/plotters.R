#' Plot 2D heatmaps of epidemic metrics
#'
#' @param scan_results Results from scan_2d_parameter_space()
#' @param metrics_to_plot Vector of metric names to plot (default: all)
#' @param overall_only If TRUE, only plot overall metrics (class == -1)
#' @returns List of ggplot objects
#' @export
plot_2d_heatmaps <- function(scan_results,
                             metrics_to_plot = NULL,
                             overall_only = TRUE) {
  if (overall_only) {
    plot_data <- scan_results[scan_results$class == -1, ]
  } else {
    plot_data <- scan_results
  }
  if (data.table::is.data.table(plot_data)) {
    plot_df <- as.data.frame(plot_data)
  } else {
    plot_df <- plot_data
  }

  plot_df <- plot_df |>
    dplyr::mutate(residence_time = 1 / switch_rate)

  plots <- list()

  for (metric in metrics_to_plot) {
    if (!metric %in% names(plot_data)) {
      message(glue::glue("Metric '{metric}' not found in data"))
      stop()
    }

    metric_labels <- list(
      duration = "Epidemic Duration",
      final_size = "Final Size",
      ptEI = "Peak Time (E+I)",
      ptI = "Peak Time (I)",
      psEI = "Peak Size (E+I)",
      psI = "Peak Size (I)",
      hit = "Herd immunity threshold",
      hit_time = "Time to HIT",
      rt_at_hit = "Rt at HIT"
    )

    plot_title <- metric_labels[[metric]] %||% metric
    plot_data_metric <- plot_data

    p <- ggplot2::ggplot(
      plot_data_metric,
      ggplot2::aes(x = epsilon, y = residence_time, fill = .data[[metric]])
    ) +
      ggplot2::geom_tile() +
      viridis::scale_fill_viridis(
        option = "mako",
        name = plot_title,
        na.value = "grey50"
      ) +
      ggplot2::labs(
        title = plot_title,
        x = "Assortativity (ε)",
        y = "Expected Residence Time"
      ) +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.position = "right"
      )
    plots[[metric]] <- p
  }
  plots
}

#' Save all heatmap plots
#'
#' @param plot_list List of ggplot objects from plot_2d_heatmaps()
#' @param output_dir Directory to save plots
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @export
save_heatmaps <- function(plot_list,
                          output_dir = "heatmaps",
                          width = 8,
                          height = 6) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  for (metric_name in names(plot_list)) {
    pngfile <- file.path(output_dir, paste0(metric_name, "_heatmap.png"))
    ggplot2::ggsave(
      pngfile,
      plot = plot_list[[metric_name]],
      width = width,
      height = height,
      dpi = 300
    )
    pdffile <- file.path(output_dir, paste0(metric_name, "_heatmap.pdf"))
    ggplot2::ggsave(
      pdffile,
      plot = plot_list[[metric_name]],
      width = width,
      height = height,
      device = cairo_pdf
    )
    message(glue::glue("Saved {pngfile} and {pdffile}"))
  }
}


plot_epidemics <- function(nt) {
  nt_long <- nt |>
    tidyr::pivot_longer(
      cols = tidyr::matches("^[SEIR][0-9]+$"),
      names_to = c("compartment", "class"),
      names_pattern = "([A-Z])([0-9]+)",
      values_to = "value"
    )

  nt_long$class <- as.factor(nt_long$class)

  a <- ggplot2::ggplot(nt_long, ggplot2::aes(
    x = time,
    y = value,
    color = class
  )) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~compartment, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = MetBrewer::met.brewer("Isfahan1",
        n = length(unique(nt_long$class))
      )
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::labs(
      title = glue::glue("SEIR Epidemic (n_activity = {length(unique(nt_long$class))})"),
      x = "Time",
      y = "Value (raw units)",
      color = "Activity class"
    )

  print(a)

  a
}

#' Plot epidemic course with annotated metrics
#'
#' @description
#' Plots the epidemic trajectory (S, E, I, R over time) and annotates
#' key epidemic metrics including peak timing/size and final size/time.
#'
#' @param model_output Neatened model output (from neaten_state)
#' @param metrics Output from compute_all() containing epidemic metrics
#' @param compartments  compartments to plot (default: c("S", "E", "I", "R"))
#' @param show_total If TRUE, plots total population dynamics (s_tot, etc.)
#' @param annotate_peak If TRUE, annotate peak timing and size
#' @param annotate_final If TRUE, annotate final epidemic size and time
#' @param title Plot title
#' @param x_truncation: scale_x_continuous used to cutoff the x axis for better visibility
#' @export
plot_epidemic_annotated <- function(model_output,
                                    metrics,
                                    compartments = c("I", "E"),
                                    show_total = TRUE,
                                    annotate_peak = TRUE,
                                    annotate_final = TRUE,
                                    title = "Example epidemic with annotated metrics",
                                    x_truncation = NULL,
                                    axis_text_size,
                                    axis_title_size) {
  if (inherits(model_output, "data.table")) {
    model_output <- as.data.frame(model_output)
  }
  if (inherits(metrics, "data.table")) {
    metrics <- as.data.frame(metrics)
  }
  model_output <- scale_results(model_output)
  overall_metrics <- metrics[metrics$class == -1, ]
  print(overall_metrics)
  plot_data <- data.frame(time = model_output$time)

  if (show_total) {
    for (comp in compartments) {
      col_name <- paste0(tolower(comp), "_tot")
      if (col_name %in% names(model_output)) {
        plot_data[[comp]] <- model_output[[col_name]]
      }
    }
  }


  plot_long <- data.frame()
  for (comp in compartments) {
    if (comp %in% names(plot_data)) {
      plot_long <- rbind(plot_long, data.frame(
        time = plot_data$time,
        value = plot_data[[comp]],
        compartment = comp
      ))
    }
  }

  plot_long$compartment <- factor(
    plot_long$compartment,
    levels = c("S", "E", "I", "R")
  )

  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(
      x = time,
      y = value,
      color = compartment
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "Time / days",
      y = "Fraction in state",
      color = "Infection state"
    ) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    )

  if (annotate_peak && nrow(overall_metrics) > 0) {
    if ("ptI" %in% names(overall_metrics) && "psI" %in% names(overall_metrics)) {
      peak_time_I <- overall_metrics$ptI
      peak_size_I <- overall_metrics$psI


      p <- p +
        ggplot2::geom_point(
          data = data.frame(x = peak_time_I, y = peak_size_I),
          ggplot2::aes(x = x, y = y),
          color = "turquoise4", size = 2, inherit.aes = FALSE
        ) +
        ggplot2::annotate("text",
          x = peak_time_I, y = peak_size_I,
          label = "list(p[t](I), p[s](I))", parse = TRUE,
          hjust = 0.5, vjust = 0,
          size = 3.5, color = "turquoise4"
        )
    }
  }

  if (annotate_final && nrow(overall_metrics) > 0) {
    if ("end_time" %in% names(overall_metrics) && "final_size" %in% names(overall_metrics)) {
      end_time <- overall_metrics |>
        dplyr::pull(end_time)
      final_size <- overall_metrics |>
        dplyr::pull(final_size)

      p <- p +
        ggplot2::geom_vline(
          xintercept = end_time, linetype = "dotted",
          color = "gray40", alpha = 0.6
        ) +
        ggplot2::annotate("text",
          x = end_time, y = 0.95,
          label = "list(f[t], f[s])", parse = TRUE,
          hjust = 1.1, vjust = 1, size = 3.5, color = "blue"
        ) + ggplot2::annotate("text",
          x = end_time, y = 0.95,
          label = "END", hjust = -0.7, size = 3.5, color = "blue"
        )
    }
  }

  if (annotate_peak && "hit" %in% names(overall_metrics)) {
    hit <- overall_metrics$hit
    p <- p +
      ggplot2::geom_hline(
        yintercept = hit, linetype = "dashed",
        color = "orange", alpha = 0.5
      ) +
      ggplot2::annotate("text",
        x = 0, y = hit,
        label = sprintf("HIT"),
        hjust = 0, size = 3.5, color = "orange"
      )
  }

  p <- p + ggplot2::theme(
    legend.position = "right",
    plot.title = ggplot2::element_text(size = axis_title_size),
    axis.text = ggplot2::element_text(size = axis_text_size),
    axis.title = ggplot2::element_text(size = axis_title_size),
    legend.text = ggplot2::element_text(size = axis_text_size)
  )

  if (!is.null(x_truncation)) {
    p <- p + ggplot2::scale_x_continuous(limits = c(0, x_truncation))
  }
  ggplot2::ggsave("annotated_epi.svg", plot = p, width = 8, height = 6)
  p
}

#' Plot all compartments (SEIR) separately
#'
#' @description
#' Creates a faceted plot showing S, E, I, R trajectories separately
#' with annotations on the relevant panels.
#'
#' @export
plot_epidemic_faceted <- function(model_output,
                                  metrics,
                                  show_total = TRUE,
                                  title = "Epidemic Dynamics by Compartment") {
  overall_metrics <- metrics[metrics$class == -1, ]

  # Prepare long format data
  compartments <- c("S", "E", "I", "R")
  plot_long <- data.frame()

  for (comp in compartments) {
    col_name <- paste0(tolower(comp), "_tot")
    if (col_name %in% names(model_output)) {
      plot_long <- rbind(plot_long, data.frame(
        time = model_output$time,
        value = model_output[[col_name]],
        compartment = comp
      ))
    }
  }

  p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::facet_wrap(~compartment, scales = "free_y", ncol = 2) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "Time",
      y = "Value (raw units)"
    ) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )

  if (nrow(overall_metrics) > 0 && "ptI" %in% names(overall_metrics)) {
    peak_time_I <- overall_metrics$ptI
    peak_size_I <- overall_metrics$psI

    p <- p +
      ggplot2::geom_vline(
        data = data.frame(compartment = "I", x = peak_time_I),
        ggplot2::aes(xintercept = x),
        linetype = "dashed", color = "red", alpha = 0.6
      )
  }

  if (nrow(overall_metrics) > 0 && "end_time" %in% names(overall_metrics)) {
    end_time <- overall_metrics$end_time

    p <- p +
      ggplot2::geom_vline(
        data = data.frame(
          compartment = rep(c("S", "E", "I", "R"), each = 1),
          x = end_time
        ),
        ggplot2::aes(xintercept = x),
        linetype = "dotted", color = "blue", alpha = 0.6
      )
  }

  p
}


#' Plot Rt over time with R0 reference line
#'
#' @description
#' Creates a plot of Rt over time, showing the epidemic threshold (Rt=1)
#' and the initial R0 value.
#'
#' @param model_output Model output with Rt column (from compute_rt_timeseries)
#' @param highlight_threshold If TRUE, highlights region where Rt > 1
#' @export
plot_rt_timeseries <- function(model_output,
                               highlight_threshold = TRUE,
                               x_truncation = NULL,
                               axis_text_size,
                               axis_title_size,
                               metrics = NULL) {
  if (inherits(model_output, "data.table")) {
    model_output <- as.data.frame(model_output)
  }

  r0 <- unique(model_output$R0)[1]

  p <- ggplot2::ggplot(model_output, ggplot2::aes(x = time, y = Rt)) +
    ggplot2::geom_line(linewidth = 1, color = "steelblue") +
    ggplot2::geom_hline(
      yintercept = 1, linetype = "dashed",
      color = "red", linewidth = 0.8
    ) +
    ggplot2::geom_hline(
      yintercept = r0, linetype = "dotted",
      color = "darkgreen", linewidth = 0.8
    ) +
    ggplot2::annotate("text",
      x = 0, y = 1,
      label = "list(R[t] == 1)", parse = TRUE,
      hjust = 0, vjust = -0.5, size = 3.5, color = "red"
    ) +
    ggplot2::annotate("text",
      x = 0, y = r0,
      label = sprintf("R[0] == %.1f", r0), parse = TRUE,
      hjust = 0, vjust = -0.5, size = 3.5, color = "darkgreen"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      x = "Time/days",
      y = expression(R[t])
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title = ggplot2::element_text(size = axis_title_size)
    )

  # Add HIT time vertical line if HIT time is found
  hit_time <- NA_real_
  if (!is.null(metrics)) {
    hit_time <- metrics$hit_time[metrics$class == -1]
  } else {
    # Fallback: find first time Rt <= 1
    idx <- which(model_output$Rt <= 1)
    if (length(idx) > 0) {
      hit_time <- model_output$time[idx[1]]
    }
  }

  if (!is.na(hit_time)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = hit_time, linetype = "dashed",
        color = "purple", alpha = 0.7
      ) +
      ggplot2::annotate("text",
        x = hit_time, y = max(model_output$Rt, na.rm = TRUE) * 0.9,
        label = "Time at HIT",
        angle = 90, vjust = -0.5, size = 3.5, color = "purple"
      )
  }


  if (highlight_threshold) {
    p <- p +
      ggplot2::geom_ribbon(
        data = model_output[model_output$Rt > 1, ],
        ggplot2::aes(ymin = 1, ymax = Rt),
        fill = "red", alpha = 0.1
      )
  }

  if (!is.null(x_truncation)) {
    p <- p + ggplot2::scale_x_continuous(limits = c(0, x_truncation))
  }

  p
}

#' Combined plot of epidemic trajectory and Rt
#'
#' @description
#' Creates a two-panel plot showing the epidemic course (I, E) and Rt
#'
#' @param model_output Model output with Rt column (from compute_rt_timeseries)
#' @param computed_metrics Computed epidemic metrics - this allows for
#' @param title the title of the plot
#' @param rt_only whether to plot only the evolution of Rt with time
#' @export
plot_epidemic_with_rt <- function(model_output, computed_metrics,
                                  compartments = c("I", "E"),
                                  title = "Epidemic Dynamics and Rt",
                                  rt_only = FALSE,
                                  axis_text_size,
                                  axis_title_size) {
  epi_endtime <- computed_metrics |>
    dplyr::filter(class == -1) |>
    dplyr::select("end_time") |>
    dplyr::pull()

  endtime_xaxis <- epi_endtime * 1.2

  overall_peakEI <- computed_metrics |>
    dplyr::filter(class == -1) |>
    dplyr::select("ptEI") |>
    dplyr::pull()

  model_output <- as.data.frame(model_output)

  if ("Rt" %in% names(model_output) && any(model_output$Rt < 1, na.rm = TRUE)) {
    rt_unity <- model_output$time[which(model_output$Rt < 1)[1]]
  } else {
    rt_unity <- NA_real_
  }

  plot_long <- data.frame()
  for (comp in compartments) {
    col_name <- paste0(tolower(comp), "_tot")
    if (col_name %in% names(model_output)) {
      plot_long <- rbind(plot_long, data.frame(
        time = model_output$time,
        value = model_output[[col_name]],
        compartment = comp
      ))
    }
  }

  plot_long$compartment <- factor(
    plot_long$compartment,
    levels = c("S", "E", "I", "R")
  )


  r0 <- unique(model_output$R0)[1]

  p2 <- ggplot2::ggplot(model_output, ggplot2::aes(x = time, y = Rt)) +
    ggplot2::geom_line(linewidth = 1, color = "#05524b") +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    ggplot2::geom_vline(
      xintercept = rt_unity,
      linetype = "dashed",
      color = "purple"
    ) +
    ggplot2::annotate("text",
      x = rt_unity,
      y = max(model_output$Rt) * 0.9,
      label = "Rt = 1",
      angle = 90,
      vjust = -0.5,
      color = "purple",
      size = 4
    ) +
    ggplot2::scale_x_continuous(limits = c(0, epi_endtime * 1.1)) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title = ggplot2::element_text(size = axis_title_size),
      legend.text = ggplot2::element_text(size = axis_text_size),
    ) +
    ggplot2::labs(
      x = "Time /days",
      y = expression(R[t])
    )
  if (rt_only == TRUE) {
    p <- p2
  } else if (rt_only == FALSE) {
    p <- p1 / p2 + patchwork::plot_layout(heights = c(1, 1))
  }
  p
}


plot_epidemic_ridges_quantitative <- function(model_output,
                                              compartment = "I",
                                              spacing = 1) {
  if (!inherits(model_output, "data.table")) {
    model_output <- data.table::as.data.table(model_output)
  }

  info <- get_info(model_output)
  n_classes <- info$num_classes
  class_sizes <- info$class_size_all

  plot_data_list <- vector("list", n_classes)

  for (i in seq_len(n_classes)) {
    class_size <- class_sizes[i]

    if (compartment == "I") {
      infected_col <- paste0("I", i)
      prevalence <- model_output[[infected_col]] / class_size
    } else if (compartment == "E") {
      exposed_col <- paste0("E", i)
      prevalence <- model_output[[exposed_col]] / class_size
    } else if (compartment == "EI") {
      exposed_col <- paste0("E", i)
      infected_col <- paste0("I", i)
      prevalence <- (model_output[[exposed_col]] + model_output[[infected_col]]) / class_size
    }

    baseline <- (n_classes - i) * spacing

    plot_data_list[[i]] <- data.frame(
      time = model_output$time,
      prevalence = prevalence,
      baseline = baseline,
      class = paste0("Class ", i),
      class_num = i
    )
  }

  plot_data <- do.call(rbind, plot_data_list)

  y_breaks <- seq(0, (n_classes - 1) * spacing, by = spacing)
  y_labels <- paste0("Class ", n_classes:1)

  ggplot(plot_data, aes(x = time, y = baseline, group = class, fill = class)) +
    geom_ridgeline(
      aes(height = prevalence),
      min_height = 0,
      alpha = 0.7,
      color = "black",
      size = 0.5
    ) +
    scale_fill_viridis_d(name = "Activity class", option = "D") +
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      sec.axis = sec_axis(~ . / spacing,
        breaks = seq(0, n_classes - 1, by = 0.2),
        labels = function(x) sprintf("%.1f", x %% 1),
        name = "Prevalence scale"
      )
    ) +
    labs(
      title = "Epidemic intensity across activity classes",
      x = "Time",
      y = "Activity class"
    ) +
    theme_minimal() +
    theme(
      legend.position = "right",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y.right = element_text(size = 8, color = "grey50")
    )
}

#' Plot 1D parameter sweeps
#'
#'
#'
#' @param results_tibble the results tibble containing parameter scan results
#' @export
plot_sweep_results_1d <- function(results_tibble) {
  all_res_bound <- results_tibble |>
    dplyr::mutate(class = as.character(class))

  classes <- unique(all_res_bound$class)

  linetypes <- setNames(ifelse(classes == "-1", "solid", "dashed"), classes)
  linewidths <- setNames(ifelse(classes == "-1", 1.2, 0.7), classes)
  linetypes <- ifelse(classes == -1, "solid", "dashed")
  linewidths <- ifelse(classes == -1, 1.5, 0.8)

  p_finalsize <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = final_size,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_finaltime <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = duration,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_peaktimeEI <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = ptEI,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_peaksizeEI <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = psEI,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_peaktimeI <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = ptI,
      colour = as.factor(claåss),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_peaksizeI <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = psI,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()


  p_hit <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = hit,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_r0 <- ggplot2::ggplot(
    all_res_bound,
    ggplot2::aes(
      x = epsilon, y = hit,
      colour = as.factor(class),
      group = class
    )
  ) +
    ggplot2::geom_line(aes(
      linetype = as.factor(class),
      linewidth = as.factor(class)
    )) +
    ggplot2::scale_linetype_manual(values = linetypes) +
    ggplot2::scale_linewidth_manual(values = linewidths) +
    ggplot2::theme_bw()

  p_all <- patchwork::wrap_plots(p_finalsize,
    p_finaltime,
    p_peaktimeEI,
    p_peaksizeEI,
    p_peaktimeI,
    p_peaksizeI,
    p_hit,
    guides = "collect"
  )


  ggplot2::ggtitle(label = glue::glue("initial cond {}"))
  print(p_all)

  ggplot2::ggsave(p_all, filename = "summary.pdf")
}
