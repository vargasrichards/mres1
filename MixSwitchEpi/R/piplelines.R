#' Run and record a 2D parameter scan
#'
#' A wrapper for parameter set construction, calibration,
#' simulation and writing out results.
#'
#' @param n_activity the number of activity classes in the model.
#' @param switch_model switching model
#' @param mix_model mixing model
#' @param resol the number of values in each direction to scan.
#' so the total number of parameter combinations run will be resol^2.
#' @param gamma_prmt the gamma distribution according to which the activity
#' scores are calculated.
#' @param initial initial condition rule. allowed values include "uniform",
#'
#' @param calib_eps the epsilon at which to calibrate to target_r0
#' @param calib_swch the switch rate at which to calibrate t0 target_r0
#' @param target_r0 the R0 to calibrate to.
#' @export
run_combination <- function(n_activity,
                            switch_model,
                            mix_model,
                            resol,
                            gamma_prmt,
                            initial,
                            calib_eps, calib_swch,
                            target_r0) {
  pars <- MixSwitchEpi::make_parms(
    n_activity = n_activity,
    epsilon = calib_eps,
    switch_rate = calib_swch,
    initial_condition_rule = initial,
    infected_fraction = 1e-03,
    activity_scheme = "gamma",
    gamma_dist = gamma_prmt,
    mixing_model = mix_model,
    switching_model = switch_model,
    pathogen_string = "sars-cov-2",
    target_r0 = target_r0
  )

  desc_model <- paste0(
    switch_model, "_", mix_model, "_",
    initial, "r0=",
    as.character(target_r0), ".csv"
  )

  epsilon_vals <- seq(0, 1, length.out = resol)
  res_vals <- c(Inf, seq(20, 1, length.out = resol))
  # scanresidence out to 5 times infectious period
  switch_vals <- 1 / res_vals
  switch_vals <- seq(0, 1, length.out = resol)
  results_list <-
    MixSwitchEpi::scan_mixswitch_framework(
      base_pars = pars,
      epsilon_vals = epsilon_vals,
      switch_vals = switch_vals,
      classed_detail = TRUE,
      switch_model = switch_model,
      mix_model = mix_model
    )
  results_dt <- data.table::rbindlist(results_list, fill = TRUE)

  if (!"residence_time" %in% names(results_dt)) {
    results_dt[, residence_time := ifelse(is.na(switch_rate) | switch_rate == 0,
      Inf, 1 / switch_rate
    )]
  }

  write.csv(results_dt, file = desc_model)

  outcome_vars <- c(
    "final_size",
    "hit",
    "psEI",
    "rt_at_hit",
    "r0",
    "spectral_gap"
  )

  plot_dt <- dplyr::filter(results_dt, class == -1)


  make_heatmap <- function(var) {
    ggplot2::ggplot(
      plot_dt,
      ggplot2::aes(
        x = epsilon,
        y = switch_rate,
        fill = .data[[var]]
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_x_continuous(name = expression(epsilon)) +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(title = var, fill = var) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size = 11),
        axis.text.y = ggplot2::element_text(size = 9)
      )
  }
  plots <- lapply(outcome_vars, make_heatmap)
  p <- patchwork::wrap_plots(plots, ncol = 2)

  fig_name <- paste0(
    "ht_switchrate_",
    switch_model, "_",
    mix_model, "_",
    initial, "_",
    "r0=", target_r0, "_",
    "gamma=", paste(gamma_prmt, collapse = "-"), "_",
    "nact=", n_activity,
    ".png"
  )

  ggplot2::ggsave(
    filename = fig_name,
    plot = p,
    width = 12,
    height = 8,
    dpi = 300
  )
}

#' Run a standard set of routines
#'
#' Wrapper for testing combinations of parameters / conditions.
#'
#' @param n_activity the number of activity classes in the model.
#' @param resol the number of values in each direction to scan.
#' @param gamma_prmt the gamma distribution according to which the activity
#' scores are calculated.
#' @export
metaroutine <- function(n_activity = 5,
                        resol = 10,
                        gamma_prmt = list(mean = 30, variance = 2500)) {
  run_combination(n_activity,
    switch_model = "simple_switching",
    mix_model = "garnett_mixing",
    resol,
    gamma_prmt,
    initial = "uniform",
    calib_eps = 0.5,
    calib_swch = 0,
    target_r0 = 2
  )

  run_combination(n_activity,
    switch_model = "simple_switching",
    mix_model = "garnett_mixing",
    resol,
    gamma_prmt,
    initial = "uniform",
    calib_eps = 0.5, calib_swch = 0,
    target_r0 = 1.5
  )

  run_combination(n_activity,
    switch_model = "simple_switching",
    mix_model = "garnett_mixing",
    resol,
    gamma_prmt,
    initial = "uniform",
    calib_eps = 0.5,
    calib_swch = 0,
    target_r0 = 1
  )

  run_combination(n_activity,
    switch_model = "simple_switching",
    mix_model = "exponential",
    resol,
    gamma_prmt,
    initial = "uniform",
    calib_eps = 0.5,
    calib_swch = 0,
    target_r0 = 2
  )
}

#' Plot combined results
#'
#' Plot the results from a saved CSV file.
#'
#' Makes a heatmap by residency time.
#'
#' @param results_file the results file from
#' which the results to be plotted are read.
#' @export
plot_combin <- function(results_file) {
  rs <- read.csv(file = results_file)

  rs_whole <- rs |> # the overall metrics
    dplyr::filter(class == -1)


  p <- ggplot2::ggplot(data = )


  plot_combin <- function(results_file, out_png = NULL, out_pdf = NULL) {
    rs <- read.csv(file = results_file)

    # ensure residence_time exists
    if (!"residence_time" %in% names(rs)) {
      if ("switch_rate" %in% names(rs)) {
        rs$residence_time <- ifelse(is.na(rs$switch_rate) | rs$switch_rate == 0, Inf, 1 / rs$switch_rate)
      } else if ("res_vals" %in% names(rs)) {
        rs$residence_time <- rs$res_vals
      } else {
        stop("results file does not contain a switch_rate or res_vals column to reconstruct residence times")
      }
    }

    rs_whole <- rs |> dplyr::filter(class == -1)

    # pick outcome columns: numeric columns excluding parameter columns
    param_cols <- c("epsilon", "switch_rate", "residence_time", "class", "res_vals")
    candidate <- setdiff(names(rs_whole), param_cols)
    numeric_cols <- candidate[sapply(rs_whole[candidate], is.numeric)]
    if (length(numeric_cols) == 0) stop("no numeric outcome columns found in results file")

    plot_dt <- rs_whole |>
      dplyr::mutate(
        res_fac = forcats::fct_rev(factor(as.character(.data$residence_time), levels = as.character(sort(unique(.data$residence_time)))))
      ) |>
      dplyr::select(epsilon, residence_time, res_fac, dplyr::all_of(numeric_cols)) |>
      tidyr::pivot_longer(cols = dplyr::all_of(numeric_cols), names_to = "outcome", values_to = "value")

    p <- ggplot2::ggplot(plot_dt, ggplot2::aes(x = epsilon, y = res_fac, fill = value)) +
      ggplot2::geom_tile() +
      ggplot2::facet_wrap(~outcome, scales = "free") +
      ggplot2::scale_y_discrete(labels = function(x) ifelse(x == "Inf", "∞", x), name = "Residence time") +
      ggplot2::scale_x_continuous(name = expression(epsilon)) +
      ggplot2::scale_fill_viridis_c(na.value = "grey90") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.text = ggplot2::element_text(), axis.text.y = ggplot2::element_text(size = 9))

    if (!is.null(out_png)) ggplot2::ggsave(filename = out_png, plot = p, width = 12, height = 8, dpi = 300)
    if (!is.null(out_pdf)) ggplot2::ggsave(filename = out_pdf, plot = p, width = 12, height = 8, device = cairo_pdf)

    p
  }
}
