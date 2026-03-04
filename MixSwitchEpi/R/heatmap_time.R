# Class × time heatmaps of infection prevalence
#
# Simulate one or two parameter-space points and visualise
# (E + I) / class_size through time as a heatmap.


#' Simulate and plot a class × time infection heatmap
#'
#' For each parameter set supplied, runs the SEIRS model and
#' produces a heatmap with activity class on the x-axis,
#' time on the y-axis, and \code{(E + I) / class_size} as the
#' fill (inferno palette).
#'
#' When two parameter sets are provided the heatmaps are placed
#' side-by-side for direct comparison.
#'
#' @param parms_list a named list of 1 or 2 parameter lists,
#'   each suitable for \code{run_system_wparms()}.
#'   Names are used as panel titles.
#' @param time_range optional numeric vector \code{c(t_start, t_end)}
#'   to restrict the y-axis.  Defaults to the full simulation range.
#' @param save logical; if \code{TRUE} saves a PDF to \code{output_file}.
#' @param output_file path for the PDF.
#' @param width,height PDF dimensions in inches.
#' @return A patchwork object (invisibly).
#' @export
heatmap_through_time <- function(parms_list,
                                 time_range = NULL,
                                 save = TRUE,
                                 output_file = "output/class_time_heatmap.pdf",
                                 width = 8, height = 6) {
  if (!is.list(parms_list[[1]])) {
    stop("parms_list must be a list of parameter lists (length 1 or 2)")
  }
  if (length(parms_list) > 2) {
    warning("Only the first two parameter sets will be plotted")
    parms_list <- parms_list[1:2]
  }
  # Default names

  if (is.null(names(parms_list))) {
    names(parms_list) <- paste("Scenario", seq_along(parms_list))
  }

  # Simulate each scenario and build a long data-frame
  long_list <- lapply(names(parms_list), function(nm) {
    pars <- parms_list[[nm]]
    out <- run_system_wparms(pars, model = seirs_act())
    n_act <- pars$n_activity
    cs <- pars$class_sizes

    # Extract E_i, I_i columns and compute prevalence
    e_cols <- paste0("E", seq_len(n_act))
    i_cols <- paste0("I", seq_len(n_act))

    # out is a data.table
    times <- as.numeric(out$time)
    prev_mat <- matrix(0, nrow = length(times), ncol = n_act)
    for (k in seq_len(n_act)) {
      prev_mat[, k] <- (as.numeric(out[[e_cols[k]]]) +
        as.numeric(out[[i_cols[k]]])) / cs[k]
    }

    # Thin to at most 500 time-points per scenario
    keep_idx <- unique(round(seq(1, length(times),
      length.out = min(500, length(times))
    )))
    times <- times[keep_idx]
    prev_mat <- prev_mat[keep_idx, , drop = FALSE]

    # To long format
    df <- data.frame(
      scenario = nm,
      time = rep(times, n_act),
      class = rep(seq_len(n_act), each = length(times)),
      prevalence = as.numeric(prev_mat),
      stringsAsFactors = FALSE
    )
    df
  })

  long_df <- do.call(rbind, long_list)
  long_df$class <- factor(long_df$class)
  long_df$scenario <- factor(long_df$scenario, levels = names(parms_list))

  # Optional time subsetting

  if (!is.null(time_range)) {
    long_df <- long_df[long_df$time >= time_range[1] &
      long_df$time <= time_range[2], ]
  }

  p <- ggplot2::ggplot(long_df, ggplot2::aes(
    x = class,
    y = time,
    fill = prevalence
  )) +
    ggplot2::geom_tile() +
    viridis::scale_fill_viridis(
      option = "inferno",
      name = "(E + I) / N_class",
      na.value = "grey80"
    ) +
    ggplot2::scale_y_reverse(expand = c(0, 0)) +
    ggplot2::labs(
      x = "Activity class",
      y = "Time (days)"
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 13, face = "bold")
    )

  if (length(parms_list) > 1) {
    p <- p + ggplot2::facet_wrap(~scenario, nrow = 1)
  } else {
    p <- p + ggplot2::ggtitle(names(parms_list)[1])
  }

  if (save) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggplot2::ggsave(output_file, p, width = width, height = height)
    message("Saved heatmap to ", output_file)
  }

  invisible(p)
}
