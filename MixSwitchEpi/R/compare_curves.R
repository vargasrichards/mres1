#' Homogeneous SEIR model
#' @export
seir_model <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    ds <- -beta * S * I
    de <- beta * S * I - sigma * E
    di <- sigma * E - gamma * I
    dr <- gamma * I
    list(c(ds, de, di, dr))
  })
}

#' Find the final size for a homogeneous SEIR epidemic given R0
#' Uses the standard SIR/SEIR final-size equation (no demography):
#'   \eqn{s_\infty = \exp(-R_0 (1 - s_\infty))}
#' where \eqn{s_\infty} is the final susceptible fraction and
#' the final size is \eqn{z = 1 - s_\infty}.
#' @param r0 the basic reproduction number
#' @returns the final size (fraction of population infected by end of epidemic)
#' @export
solve_r0_fs <- function(r0) {
  # Solve s_inf - exp(-r0 * (1 - s_inf)) = 0 for s_inf in (0, 1)
  final_size_eq <- function(s_inf) {
    s_inf - exp(-r0 * (1 - s_inf))
  }
  root <- stats::uniroot(final_size_eq, lower = 1e-12, upper = 1 - 1e-12)
  final_size <- 1 - root$root
  final_size
}


#' Homogeneous SEIR trajectory (mass-action)
#'
#' Convenience function that simulates a homogeneous SEIR epidemic
#' with fixed gamma, sigma and beta chosen to match the requested target R0.
#'
#' @param target_r0 Desired basic reproduction number R0
#' @param initial_frac Initial exposed fraction (default 1e-3)
#' @returns Data.frame with columns time, S, E, I, R
#' @export
scalar_reference <- function(target_r0, initial_frac = 1e-3) {
  stopifnot(0 < initial_frac && initial_frac < 1)
  parameters <- basic_scalar_parms("sars-cov-2")
  parameters$beta <- target_r0 * parameters$gamma

  assertthat::assert_that( # extra safety here
    abs(parameters$beta / parameters$gamma - target_r0) < 1e-6,
    msg = sprintf(
      "Basic scalar parameters do not match target R0 = %.2f",
      target_r0
    )
  )
  state <- c(
    S = 1 - initial_frac,
    E = initial_frac,
    I = 0,
    R = 0
  )
  times <- seq(0,
    1500,
    by = 0.1
  )

  out <- deSolve::ode(
    y = state,
    times = times,
    func = seir_model,
    parms = parameters
  )
  as.data.frame(out)
}

#' Metrics for homogeneous SEIR epidemic
#'
#' Computes peak times/sizes, final size, HIT and returns a metrics list
#' for the homogeneous (mass-action) SEIR epidemic defined by target R0.
#'
#' @param r0 Target basic reproduction number
#' @returns A named list with ptI, psI, ptEI, psEI,
#' hit, final_size, end_time, r0
#' @export
homogeneous_metrics <- function(r0) {
  A <- scalar_reference(target_r0 = r0)

  peak_EI_row <- A |>
    dplyr::mutate(EI = E + I) |>
    dplyr::filter(EI == max(EI))

  ptEI <- peak_EI_row |> dplyr::pull(time)
  psEI <- peak_EI_row |> dplyr::pull(EI)

  peak_I_row <- A |> dplyr::filter(I == max(I))
  ptI <- peak_I_row |> dplyr::pull(time)
  psI <- peak_I_row |> dplyr::pull(I)

  end_row <- A |>
    dplyr::filter(I < 1e-5) |>
    dplyr::filter(I == max(I))
  final_size <- 1 - (end_row |> dplyr::pull(S))
  final_time <- end_row |> dplyr::pull(time)

  hit <- 1 - (1 / r0)

  # sanity check via Rt crossing
  A <- A |> dplyr::mutate(Rt = r0 * S)
  Rt_cross_row <- A |>
    dplyr::filter(Rt <= 1)
  if (nrow(Rt_cross_row) > 0) {
    Rt_cross_row <- Rt_cross_row |> dplyr::filter(Rt == max(Rt))
    HIT_check <- 1 - (Rt_cross_row |> dplyr::pull(S))
    if (!is.na(HIT_check) && abs(HIT_check - hit) > 1e-3) {
      warning(
        "HIT check mismatch: analytic vs Rt crossing differ by ",
        abs(HIT_check - hit)
      )
    }
  }

  metrics <- list(
    ptI = ptI,
    psI = psI,
    ptEI = ptEI,
    psEI = psEI,
    hit = hit,
    final_size = final_size,
    end_time = final_time,
    r0 = r0
  )
  metrics
}


#' Pivot SEIR epidemic df + prepare for plotting
#'
#' Pivots the epidemic dataframe into suitable plotting format
#' and then binds an epidemic label
#'
#' @param df the epidemic dataframe
#' @param model_name string descriptor for the model
#' @export
pivot_seir <- function(df, model_name) {
  pivoted_df <- df |>
    tidyr::pivot_longer(!time, names_to = "infec_state", values_to = "value") |>
    dplyr::mutate(infec_state = forcats::fct_relevel(
      infec_state,
      "S",
      "E",
      "I",
      "R"
    )) |>
    dplyr::mutate(model_name = model_name)
  pivoted_df
}

#' Plot homogeneous mass action model
#'
#'
#' @param n_activity the number of activity classes
#' involved.
#' @param r0 the target r0 for the homogeneous model.
#' Note that this uniquely defines the epidemic given we have
#' fixed params
#' @return a ggplot2 object and saves output to pdf file
#' @export
homogeneous_curve <- function(r0) {
  epi_df1 <- scalar_reference(target_r0 = r0)
  epi_df2 <- scalar_reference(target_r0 = 2 * r0)

  # should get the max time around here... otherwise we're focussing too much on irrelevant bits.


  # here we construct the heterogeneously active population

  epi_df3 <- pivot_seir(epi_df1, model_name = "Homogeneous pop")
  epi_df4 <- pivot_seir(epi_df2, model_name = "double_r0")

  epi_df_long <- rbind(epi_df3, epi_df4)

  p <- ggplot2::ggplot(
    data = epi_df_long,
    ggplot2::aes(x = time, y = value, group = model_name, colour = model_name)
  ) +
    ggplot2::geom_line() +
    ggplot2::facet_grid(rows = dplyr::vars(infec_state)) +
    ggplot2::theme_bw() +
    MetBrewer::scale_colour_met_d("Isfahan2") +
    ggplot2::labs(
      title = bquote("SEIR epidemics: " ~ R[0] == .(r0)),
      x = "Time /days"
    ) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, size = 12))

  ggplot2::ggsave("epidemic_curve.pdf")
  print(p)
  p
}

#' Empirical Haslemere scan
#'
#' Runs a scan over a range of residence times (switching rates)
#' using empirical Haslemere contact data.
#'
#' @param target_r0 basic reproduction number to calibrate to.
#' @param calibrate_each logical (default \code{FALSE}). When \code{TRUE} the
#'   model is re-calibrated at every finite residence time so that R0 equals
#'   \code{target_r0} (the reference R0, computed at infinite residence time /
#'   no switching). The recalibrated results are returned alongside the
#'   fixed-beta results, tagged with a \code{calibration} column
#'   (\code{"Fixed p"} / \code{"Recalibrated"}).
#' @param switch_model the switching model to use. See \code{simple_switching}
#' for details.
#' @param contact_file path to the empirical contact matrix CSV file.
#' @param activity_scores vector of activity scores for each class.
#' @export
haslemere_empirical <- function(target_r0,
                                calibrate_each = FALSE,
                                switch_model = "simple_switching",
                                contact_file = "data/Haslemere.csv",
                                activity_scores = c(
                                  0.3146078,
                                  1.1479220,
                                  2.3361036,
                                  4.2211382,
                                  9.3466695
                                )) {
  raw_scores <- activity_scores

  empirical_m <- as.matrix(utils::read.csv(contact_file, header = FALSE))

  n_activity <- length(raw_scores)

  scalars <- basic_scalar_parms("sars-cov-2")
  base_pars <- list()
  base_pars$beta <- scalars$beta
  base_pars$sigma <- scalars$sigma
  base_pars$gamma <- scalars$gamma
  base_pars$n_activity <- n_activity
  base_pars$N <- 1e6
  base_pars$class_sizes <- rep(
    base_pars$N / n_activity,
    n_activity
  )
  base_pars$activity_scores <- raw_scores
  base_pars$t_end <- 3000
  base_pars$switching_model <- switch_model
  base_pars$omega <- 0

  ic <- make_initial_conditions(
    n_activity = n_activity,
    initial_condition_rule = "uniform",
    fraction_exposed = 1e-3
  )

  base_pars$E0 <- ic$E0
  base_pars$S0 <- ic$S0
  if (switch_model == "simple_switching") {
    base_pars$swch <- simple_switching(
      switch_rate = 0,
      num_classes = n_activity
    )
  } else if (switch_model == "adjacent_switching") {
    base_pars$swch <- adjacent_switching(
      switch_rate = 0,
      num_classes = n_activity
    )
  } else {
    stop("Switch model not recognised")
  }

  base_pars$m <- empirical_m
  residence_times <- 1:50
  switch_rates <- c(1 / residence_times, 0)

  # Calibrate at Inf residence time (no switching) — the reference beta
  calibrated <- calibrate_parms(
    preliminary_parms = base_pars,
    target_r0 = target_r0
  )
  assertthat::assert_that(abs(compute_r0(calibrated) - target_r0) < 1e-6,
    msg = "Calibration failed to achieve target R0"
  )

  # ── helper: simulate one switch_rate and return a metrics data.table row.
  # Adaptively doubles t_end (up to max_doublings times) if Rt never crosses
  # below 1 within the initial window — this can happen at fast switching
  # where the epidemic is very slow.
  run_one <- function(sr,
                      pars_in,
                      switch_model = "simple_switching",
                      max_doublings = 4L) {
    if (switch_model == "simple_switching") {
      pars_in$swch <- MixSwitchEpi::simple_switching(
        switch_rate = sr,
        num_classes = n_activity
      )
    } else if (switch_model == "adjacent_switching") {
      pars_in$swch <- MixSwitchEpi::adjacent_switching(
        switch_rate = sr,
        num_classes = n_activity
      )
    } else {
      stop("Switch model not recognised")
    }
    epi_out <- run_system_wparms(parms = pars_in, model = seirs_act())
    m <- compute_all(epi_out, classed_detail = FALSE, params = pars_in)

    # If HIT is NA  — extend t_end and retry
    n_doublings <- 0L
    while (is.na(m$hit[m$class == -1]) && n_doublings < max_doublings) {
      pars_in$t_end <- pars_in$t_end * 2L
      n_doublings <- n_doublings + 1L
      epi_out <- run_system_wparms(parms = pars_in, model = seirs_act())
      m <- compute_all(epi_out, classed_detail = FALSE, params = pars_in)
    }

    m$switch_rate <- sr
    m$residence_time <- ifelse(sr == 0, Inf, 1 / sr)
    m
  }

  # ── Fixed-beta scan ───────────────────────────────────────────────────────
  results_fixed <- lapply(switch_rates, run_one, pars_in = calibrated, switch_model = switch_model)
  dt_fixed <- data.table::rbindlist(results_fixed, fill = TRUE)
  dt_fixed$calibration <- "Fixed p"

  # ── Per-residence-time recalibration scan (optional) ─────────────────────
  if (calibrate_each) {
    results_recal <- vector("list", length(switch_rates))
    for (k in seq_along(switch_rates)) {
      sr <- switch_rates[k]
      if (sr == 0) {
        # Inf residence time: this is identical to the reference calibration
        results_recal[[k]] <- run_one(sr, calibrated, switch_model = switch_model)
      } else {
        # Build a fresh prelim with this switch rate and re-calibrate
        prelim_k <- calibrated # copy all fields incl. empirical m
        if (switch_model == "simple_switching") {
          prelim_k$swch <- MixSwitchEpi::simple_switching(
            switch_rate = sr,
            num_classes = n_activity
          )
        } else if (switch_model == "adjacent_switching") {
          prelim_k$swch <- MixSwitchEpi::adjacent_switching(
            switch_rate = sr,
            num_classes = n_activity
          )
        } else {
          stop("Switch model not recognised")
        }
        recal_k <- calibrate_parms(
          preliminary_parms = prelim_k,
          target_r0 = target_r0
        )
        results_recal[[k]] <- run_one(sr, recal_k)
      }
    }
    dt_recal <- data.table::rbindlist(results_recal, fill = TRUE)
    dt_recal$calibration <- "Recalibrated"
    results_dt <- data.table::rbindlist(list(dt_fixed, dt_recal), fill = TRUE)
  } else {
    results_dt <- dt_fixed
  }

  # calculate the R0 for a homogeneous scenario using the same transmission rate (p/beta)
  # as the calibrated structured model (at residence_time = Inf)
  hom_r0 <- (calibrated$beta * mean(activity_scores)) / calibrated$gamma
  results_dt$hom_r0 <- hom_r0

  outdir <- glue::glue("output/haslemere{target_r0}")
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  csvpath <- file.path(outdir, "haslemere_scan_results.csv")
  data.table::fwrite(results_dt, file = csvpath)

  overall <- results_dt[results_dt$class == -1, ]

  # single-panel plots for chosen metrics: include empirical 'hit' (lowercase) if present
  metrics_to_plot <- c("final_size", "psEI", "hit", "ptEI", "r0")
  metrics_to_plot <- metrics_to_plot[metrics_to_plot %in% names(overall)]

  # compute homogeneous expectations using the same p as the structured case
  hom_metrics <- homogeneous_metrics(hom_r0)

  plotdir <- file.path(outdir, "single_panels")
  if (!dir.exists(plotdir)) dir.create(plotdir, recursive = TRUE)

  summary_df <- as.data.frame(overall)
  finite_rt <- summary_df$residence_time[is.finite(summary_df$residence_time)]
  max_finite_rt <- if (length(finite_rt) > 0) {
    max(finite_rt,
      na.rm = TRUE
    )
  } else {
    50
  }
  inf_plot_val <- max_finite_rt * 1.25
  summary_df$res_time_plot <- ifelse(is.infinite(summary_df$residence_time),
    inf_plot_val, summary_df$residence_time
  )
  gamma_time <- 1 / base_pars$gamma

  list(results = results_dt, csv = csvpath)
}


#' Haslemere grid scan
#'
#' Produces a gridded plot of the effects of switching on epi outcome.
#' Each element of \code{r0_list} produces one row-panel in the output.
#'
#' @param r0_list list or vector of Target R0 values to scan
#' @param calibrate_each logical (default \code{FALSE}). Passed through to
#'   \code{\link{haslemere_empirical}}. When \code{TRUE} the model is
#'   re-calibrated at each finite residence time to match \code{target_r0}.
#' @export
haslemere_grid <- function(r0_list, calibrate_each = FALSE, plot = TRUE, save_path = "Haslemere_outcome_1d.pdf") {
  grid_results <- list()
  for (set_r0 in r0_list) {
    result_list <- haslemere_empirical(set_r0, calibrate_each = calibrate_each)
    result_df <- result_list$results
    result_elab <- result_df |>
      dplyr::mutate(calib_at_r0 = set_r0)

    grid_results[[length(grid_results) + 1]] <- result_elab
  }

  final_df <- data.table::rbindlist(grid_results)

  if (plot) {
    p <- plot_haslemere_grid(final_df, save_path = save_path)
    print(p)
  }

  final_df
}

#' Plot the grid of (e.g.,) Haslemere results
#'
#' Plots a grid where rows define the Calibration R0 and columns define the
#' metric.  Manually constructs rows to ensure consistent Y-scales for each
#' metric across R0 values (columns) while allowing scales to differ between
#' metrics.  Includes homogeneous reference lines with labels.
#'
#' When \code{grid_results} contains a \code{calibration} column (produced by
#' \code{\link{haslemere_empirical}} with \code{calibrate_each = TRUE}), lines
#' for each calibration series are drawn in distinct colours with a legend, so
#' the \emph{Fixed p} and \emph{Recalibrated} sweeps can be compared
#' directly on the same facets.
#'
#' @param grid_results a dataframe of results from \code{\link{haslemere_grid}}
#'   or \code{\link{haslemere_empirical}}.
#' @param gamma recovery rate, used only for annotations (default 1/4).
#' @param save_path if strictly character, saves the combined plot to this file. Set to NULL to disable.
#' @export
plot_haslemere_grid <- function(grid_results, gamma = 1 / 4, save_path = "Haslemere_outcome_1d.pdf") {
  # Detect whether a calibration-series column is present
  has_calib_col <- "calibration" %in% names(grid_results)

  # Palette for the calibration series (only used when has_calib_col is TRUE)
  calib_palette <- c("Fixed p" = "#1b6ca8", "Recalibrated" = "#c0392b")

  # Metric codes to plot
  desired_metrics <- c("final_size", "hit", "ptEI", "psEI", "r0")

  # metric codes to plotmath-parseable labels
  lab_map <- c(
    "r0"         = "R[0]",
    "final_size" = "Final~size",
    "hit"        = "HIT",
    "ptEI"       = "Peak~time~(days)",
    "psEI"       = "Peak~size"
  )
  lab_map <- lab_map[names(lab_map) %in% names(grid_results)]

  # Split: finite-residence rows for the main lines, Inf rows for reference lines.
  # When calibrate_each is used, keep only the "Fixed p" Inf row as the
  # no-switching reference (both series converge there by design).
  grid_results_df <- as.data.frame(grid_results)
  is_inf <- !is.finite(grid_results_df$residence_time)

  if (has_calib_col) {
    inf_rows <- grid_results_df[is_inf & grid_results_df$calibration == "Fixed p", ]
  } else {
    inf_rows <- grid_results_df[is_inf, ]
  }
  fin_rows <- grid_results_df[!is_inf, ]

  # Use the actual residence_time for the x-axis (1–50)
  fin_rows$res_time_plot <- fin_rows$residence_time

  # X-axis settings (1 to 50 days)
  x_range <- range(fin_rows$res_time_plot, na.rm = TRUE)
  x_breaks_all <- sort(unique(fin_rows$res_time_plot))
  if (length(x_breaks_all) > 10) {
    display_breaks <- c(1, seq(10, max(x_range), by = 10))
    display_breaks <- sort(unique(display_breaks))
  } else {
    display_breaks <- x_breaks_all
  }
  display_labels <- as.character(display_breaks)

  # Columns to pivot: metric codes + grouping fields
  group_cols <- c("res_time_plot", "calib_at_r0")
  if (has_calib_col) group_cols <- c(group_cols, "calibration")
  if ("hom_r0" %in% names(grid_results_df)) group_cols <- c(group_cols, "hom_r0")

  # Pivot main (finite) data to long format
  long_df <- fin_rows |>
    dplyr::select(dplyr::any_of(c(group_cols, names(lab_map)))) |>
    tidyr::pivot_longer(
      cols = dplyr::any_of(names(lab_map)),
      names_to = "metric_code",
      values_to = "value"
    ) |>
    dplyr::mutate(metric_label = factor(lab_map[metric_code], levels = lab_map))

  if (has_calib_col) {
    long_df$calibration <- factor(long_df$calibration,
      levels = names(calib_palette)
    )
  }

  # Pivot Inf rows to long format for reference lines
  if (nrow(inf_rows) > 0) {
    inf_long <- inf_rows |>
      dplyr::select(dplyr::any_of(c("calib_at_r0", names(lab_map)))) |>
      tidyr::pivot_longer(
        cols = dplyr::any_of(names(lab_map)),
        names_to = "metric_code",
        values_to = "inf_value"
      ) |>
      dplyr::mutate(
        metric_label = factor(lab_map[metric_code], levels = lab_map),
        x_label_pos  = max(display_breaks, na.rm = TRUE)
      )
  } else {
    inf_long <- NULL
  }

  # Global Y limits per metric (for consistent scales across R0 rows)
  
  # First, compute homogeneous references for all R0s to ensure they are inside limits
  hom_refs <- list()
  for (cur_r0 in sort(unique(long_df$calib_at_r0))) {
    sub_df_r0 <- long_df |> dplyr::filter(calib_at_r0 == cur_r0)
    cur_hom_r0 <- if ("hom_r0" %in% names(sub_df_r0)) {
        sub_df_r0$hom_r0[1]
    } else {
        cur_r0
    }
    # Fixed-p reference (cur_hom_r0 <= cur_r0 in general)
    hom <- homogeneous_metrics(cur_hom_r0)
    hom_lower <- stats::setNames(hom, tolower(names(hom)))
    for (mcode in names(lab_map)) {
      mcode_lower <- tolower(mcode)
      if (mcode_lower %in% names(hom_lower)) {
        hom_refs[[length(hom_refs) + 1]] <- data.frame(
          metric_label = factor(lab_map[mcode], levels = lab_map),
          value = hom_lower[[mcode_lower]]
        )
      }
    }
    # Matched-R0 reference (cur_r0 exactly)
    hom_matched <- homogeneous_metrics(cur_r0)
    hom_matched_lower <- stats::setNames(hom_matched, tolower(names(hom_matched)))
    for (mcode in names(lab_map)) {
      mcode_lower <- tolower(mcode)
      if (mcode_lower %in% names(hom_matched_lower)) {
        hom_refs[[length(hom_refs) + 1]] <- data.frame(
          metric_label = factor(lab_map[mcode], levels = lab_map),
          value = hom_matched_lower[[mcode_lower]]
        )
      }
    }
  }
  hom_refs_df <- if (length(hom_refs) > 0) data.table::rbindlist(hom_refs) else data.frame(metric_label = factor(levels = lab_map), value = numeric(0))
  
  # Combine data sources to find global limits
  all_values_df <- dplyr::bind_rows(
    long_df |> dplyr::select(metric_label, value),
    hom_refs_df
  )
  if (!is.null(inf_long)) {
    all_values_df <- dplyr::bind_rows(
      all_values_df,
      inf_long |> dplyr::transmute(metric_label, value = inf_value)
    )
  }

  limits_df <- all_values_df |>
    dplyr::group_by(metric_label) |>
    dplyr::summarise(
      min_val = min(value, na.rm = TRUE),
      max_val = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::pivot_longer(cols = c(min_val, max_val), values_to = "value") |>
    dplyr::mutate(res_time_plot = max(x_range))
  
  # Expansion factor to ensure line labels aren't cutoff
  limits_df <- limits_df |>
    dplyr::group_by(metric_label) |>
    dplyr::mutate(
      range_val = max(value) - min(value),
      range_val = ifelse(range_val == 0, max(abs(value)) * 0.1 + 1e-3, range_val),
      value = ifelse(value == max(value), value + 0.1 * range_val, value - 0.1 * range_val)
    ) |>
    dplyr::ungroup()

  # Build one plot per R0 row
  r0_vals <- sort(unique(long_df$calib_at_r0), decreasing = TRUE)
  plot_list <- list()

  for (i in seq_along(r0_vals)) {
    cur_r0 <- r0_vals[i]
    sub_df <- long_df |> dplyr::filter(calib_at_r0 == cur_r0)

    cur_hom_r0 <- if ("hom_r0" %in% names(sub_df)) {
        sub_df$hom_r0[1]
    } else {
        cur_r0
    }

    # --- Homogeneous reference lines ---
    # Definition 1: Fixed-p — same beta as structured model, uniform pop at mean activity
    hom <- homogeneous_metrics(cur_hom_r0)
    hom_lower <- stats::setNames(hom, tolower(names(hom)))
    hom_vals <- list()
    if ("final_size" %in% names(lab_map)) hom_vals$final_size <- hom_lower[["final_size"]]
    if ("hit" %in% names(lab_map)) hom_vals$hit <- hom_lower[["hit"]]
    if ("ptEI" %in% names(lab_map)) hom_vals$ptEI <- hom_lower[["ptei"]]
    if ("psEI" %in% names(lab_map)) hom_vals$psEI <- hom_lower[["psei"]]
    if ("r0" %in% names(lab_map)) hom_vals$r0 <- hom_lower[["r0"]]

    hom_df <- data.frame(
      metric_code = names(hom_vals),
      hom_value = as.numeric(hom_vals),
      stringsAsFactors = FALSE
    )
    hom_df$metric_label <- factor(lab_map[hom_df$metric_code], levels = lab_map)
    hom_df$x_label_pos <- max(display_breaks, na.rm = TRUE)

    # Definition 2: Matched-R0 — mass-action SIR recalibrated to the same observed R0 (= cur_r0)
    hom_matched <- homogeneous_metrics(cur_r0)
    hom_matched_lower <- stats::setNames(hom_matched, tolower(names(hom_matched)))
    hom_matched_vals <- list()
    if ("final_size" %in% names(lab_map)) hom_matched_vals$final_size <- hom_matched_lower[["final_size"]]
    if ("hit" %in% names(lab_map)) hom_matched_vals$hit <- hom_matched_lower[["hit"]]
    if ("ptEI" %in% names(lab_map)) hom_matched_vals$ptEI <- hom_matched_lower[["ptei"]]
    if ("psEI" %in% names(lab_map)) hom_matched_vals$psEI <- hom_matched_lower[["psei"]]
    if ("r0" %in% names(lab_map)) hom_matched_vals$r0 <- hom_matched_lower[["r0"]]

    hom_matched_df <- data.frame(
      metric_code = names(hom_matched_vals),
      hom_matched_value = as.numeric(hom_matched_vals),
      stringsAsFactors = FALSE
    )
    hom_matched_df$metric_label <- factor(lab_map[hom_matched_df$metric_code], levels = lab_map)
    hom_matched_df$x_label_pos <- max(display_breaks, na.rm = TRUE)

    # --- Inf-residence (no switching) reference lines ---
    if (!is.null(inf_long)) {
      inf_cur <- inf_long |> dplyr::filter(calib_at_r0 == cur_r0)
    } else {
      inf_cur <- data.frame(
        metric_label = factor(character(0), levels = lab_map),
        inf_value = numeric(0),
        x_label_pos = numeric(0)
      )
    }

    # Base plot ---------------------------------------------------------------
    if (has_calib_col) {
      # Colour lines by calibration series; include a legend
      p <- ggplot2::ggplot(
        sub_df,
        ggplot2::aes(x = res_time_plot, y = value, colour = calibration)
      ) +
        ggplot2::geom_line(linewidth = 0.6) +
        ggplot2::scale_colour_manual(
          name = "Calibration",
          values = calib_palette,
          guide = if (i == 1) {
            ggplot2::guide_legend(title.position = "top")
          } else {
            "none"
          }
        )
    } else {
      p <- ggplot2::ggplot(sub_df, ggplot2::aes(x = res_time_plot, y = value)) +
        ggplot2::geom_line(linewidth = 0.5)
    }

    p <- p +
      ggplot2::geom_blank(
        data = limits_df,
        ggplot2::aes(x = res_time_plot, y = value),
        inherit.aes = FALSE
      ) +
      ggplot2::facet_wrap(
        ~metric_label,
        scales = "free_y", nrow = 1, strip.position = "top",
        labeller = ggplot2::label_parsed
      ) +
      ggplot2::scale_x_continuous(
        breaks = display_breaks,
        labels = display_labels,
        limits = x_range
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(y = bquote(R[0] == .(cur_r0))) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_text(angle = 90, vjust = 0.5, size = 14),
        strip.text   = ggplot2::element_text(size = 14)
      )

    # Homogeneous reference lines — Definition 1: Fixed-p (red dashed)
    if (nrow(hom_df) > 0) {
      p <- p +
        ggplot2::geom_hline(
          data = hom_df,
          ggplot2::aes(yintercept = hom_value, linetype = "Hom. (fixed p)"),
          colour = "#E41A1C", linewidth = 0.6
        )
    }

    # Homogeneous reference lines — Definition 2: Matched-R0 (orange longdash)
    if (nrow(hom_matched_df) > 0) {
      p <- p +
        ggplot2::geom_hline(
          data = hom_matched_df,
          ggplot2::aes(yintercept = hom_matched_value, linetype = "Hom. (matched R0)"),
          colour = "#FF7F00", linewidth = 0.6
        )
    }

    # Inf-residence (no switching) reference lines + labels (green, dotted)
    if (nrow(inf_cur) > 0) {
      p <- p +
        ggplot2::geom_hline(
          data = inf_cur,
          ggplot2::aes(yintercept = inf_value, linetype = "Heterogeneous pop. \n no switching"),
          colour = "#4DAF4A", linewidth = 0.6
        )
    }

    # Build the linetype + colour overrides in the same order as the breaks.
    lt_breaks <- character(0)
    lt_values <- character(0)
    lt_colors <- character(0)

    if (nrow(hom_df) > 0) {
      lt_breaks <- c(lt_breaks, "Hom. (fixed p)")
      lt_values <- c(lt_values, "dashed")
      lt_colors <- c(lt_colors, "#E41A1C")
    }
    if (nrow(hom_matched_df) > 0) {
      lt_breaks <- c(lt_breaks, "Hom. (matched R0)")
      lt_values <- c(lt_values, "longdash")
      lt_colors <- c(lt_colors, "#FF7F00")
    }
    if (nrow(inf_cur) > 0) {
      lt_breaks <- c(lt_breaks, "Heterogeneous pop. \n no switching")
      lt_values <- c(lt_values, "dotted")
      lt_colors <- c(lt_colors, "#4DAF4A")
    }
    names(lt_values) <- lt_breaks

    p <- p +
      ggplot2::scale_linetype_manual(
        name = "Reference models",
        breaks = lt_breaks,
        values = lt_values,
        guide = if (i == 1) {
          ggplot2::guide_legend(title.position = "top", override.aes = list(colour = lt_colors))
        } else {
          "none"
        }
      )

    if (i > 1) {
      p <- p + ggplot2::theme(
        strip.text       = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
      )
    }
    if (i < length(r0_vals)) {
      p <- p + ggplot2::theme(
        axis.text.x  = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
      )
    } else {
      p <- p + ggplot2::labs(
        x = "Expected residence time in activity class (days)"
      ) + ggplot2::theme(
        axis.title.x = ggplot2::element_text(size = 14),
        axis.text.x  = ggplot2::element_text(size = 13)
      )
    }
    p <- p + ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 12),
      legend.position = if (has_calib_col && i == 1) "right" else "none"
    )
    plot_list[[i]] <- p
  }
  pgrid <- patchwork::wrap_plots(plot_list, ncol = 1)

  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = save_path,
      plot     = pgrid,
      width    = 12,
      height   = 3 * length(r0_vals)
    )
  }
  pgrid
}


make_model_palette <- function(num_models) {
  colors <- RColorBrewer::brewer.pal(n = num_models, name = "Set1")
  colors
}

#'
#' Make parameter sets for model fitting
#'
make_parameter_sets <- function(N,
                                r0,
                                num_activity,
                                haslemere_activity_scores,
                                path_to_empirical_mix = "data/Haslemere.csv") {
  build_prelim <- function(epsilon, switch_rate, empirical_mix = TRUE) {
    prelim <- basic_scalar_parms("sars-cov-2")

    prelim$n_activity <- num_activity
    prelim$N <- N
    prelim$t_end <- 1800
    prelim$class_sizes <- rep(N / num_activity, num_activity)
    prelim$activity_scores <- haslemere_activity_scores
    prelim$epsilon <- epsilon

    prelim$swch <- MixSwitchEpi::simple_switching(
      switch_rate = switch_rate,
      num_classes = num_activity
    )
    if (!empirical_mix) {
      prelim$m <- MixSwitchEpi::garnett_mixing(
        n_activity = num_activity,
        epsilon_assort = epsilon,
        act_vector = haslemere_activity_scores,
        class_size_vector = prelim$class_sizes
      )
    } else if (empirical_mix) {
      prelim$m <- as.matrix(utils::read.csv(path_to_empirical_mix,
        header = FALSE
      ))
    }
    # initial conditions: uniform small fraction exposed
    ic <- make_initial_conditions(
      n_activity = num_activity,
      initial_condition_rule = "uniform",
      fraction_exposed = 1e-3
    )
    prelim$E0 <- ic$E0
    prelim$S0 <- ic$S0
    prelim$omega <- 0
    prelim
  }

  # Build and calibrate each model variant
  # Homogeneous: truly uniform activity (mean of empirical scores),  proportionate
  # mixing, no switching.  This is the only variant where the analytical HIT
  # (1 - 1/R0) and final-size equation should match exactly.
  prelim_hom <- build_prelim(epsilon = 0, switch_rate = 0)
  mean_act <- mean(haslemere_activity_scores)
  prelim_hom$activity_scores <- rep(mean_act, num_activity)
  # rebuild the mixing matrix with the now-uniform activity scores
  prelim_hom$m <- MixSwitchEpi::garnett_mixing(
    n_activity = num_activity,
    epsilon_assort = 0,
    act_vector = prelim_hom$activity_scores,
    class_size_vector = prelim_hom$class_sizes
  )

  prelim_prop <- build_prelim(epsilon = 0, switch_rate = 0, empirical_mix = FALSE)
  prelim_np_ns <- build_prelim(epsilon = 0.3, switch_rate = 0, empirical_mix = TRUE)
  prelim_switchy <- build_prelim(epsilon = 0.3, switch_rate = 0.3, empirical_mix = FALSE)

  # Calibrate homogeneous to target R0
  hom_cal <- calibrate_parms(preliminary_parms = prelim_hom, target_r0 = r0)
  beta_hom <- hom_cal$beta

  prop_cal <- prelim_prop
  prop_cal$beta <- beta_hom

  # Keep the more complex cases tied to the homogeneous beta (fixed-beta family
  # anchored at the homogeneous reference)
  npns_cal <- prelim_np_ns
  npns_cal$beta <- beta_hom

  switch_cal <- prelim_switchy
  switch_cal$beta <- beta_hom

  parsets <- list(
    list(
      name = "Homogeneous population",
      calibrated_parms = hom_cal
    ),
    list(
      name = "Heterogeneous pop (proportionate mix)",
      calibrated_parms = prop_cal
    ),
    list(
      name = "Heterogeneous pop (empirical mix)",
      calibrated_parms = npns_cal
    ),
    list(
      name = "With switching",
      calibrated_parms = switch_cal
    )
  )

  parsets
}

#' Metric comparison across models
#'
#' Computes the comparison metrics for a list of parameter sets.
#'
#' @param parsets A list of parameter sets.
#' @export
make_metric_comparison <- function(parsets) {
  # Create a tibble to store the comparison metric variables along with model names
  # this is then used for plotting a bar chart where groups are model metrics and colour
  # is model type (e.g., homogeneous, heterogeneous with switching, etc.)
  comparison_df <- tibble::tibble(
    model_type = character(),
    ptEI = numeric(),
    psEI = numeric(),
    ptI = numeric(),
    psI = numeric(),
    hit = numeric(),
    final_size = numeric(),
    final_time = numeric(),
    r0 = numeric()
  )
  for (pars in parsets) {
    compute_all(pars)
    comparison_df <- comparison_df |>
      dplyr::add_row(
        model_type = pars$name,
        ptEI = pars$metrics$ptEI,
        psEI = pars$metrics$psEI,
        ptI = pars$metrics$ptI,
        psI = pars$metrics$psI,
        hit = pars$metrics$hit,
        final_size = pars$metrics$final_size,
        final_time = pars$metrics$end_time,
        r0 = pars$calibrated_parms$r0
      )
  }
}

#' Bar comparison of outcome metrics across models
#'
#'
#' @param comparison_df the dataframe of comparison metrics generated by make_metric_comparison,
#' which should include columns for model_type, ptEI, psEI, ptI, psI, hit, final_size, final_time, and r0.
#' The function pivots this dataframe into a long format suitable for ggplot2,
#' where each row corresponds to a specific metric for a given model type.
#'
#' @param palette gives the colour palette to use for each model type.
#'  Should be a vector of colours of length equal to the number of model
#' types (rows) in comparison_df. Can be generated with make_model_palette.
#'
#' @export
bar_comparison <- function(comparison_df, palette) {
  comparison_long <- comparison_df |>
    tidyr::pivot_longer(
      cols = c(ptEI, psEI, ptI, psI, hit, final_size, final_time),
      names_to = "metric",
      values_to = "value"
    )
  p <- ggplot2::ggplot(
    comparison_long,
    ggplot2::aes(x = metric, y = value, fill = model_type)
  ) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::theme_bw() +
    ggplot2::labs(x = "Metric", y = "Value", fill = "Model Type") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  print(p)
}

#' Compare SEIR curves across models, using the palette to colour the curves by model type.
#' This is a more direct comparison of the epidemic trajectories than the bar_comparison of summary
#' metrics, and can reveal differences in timing and shape of the epidemics that may not be
#' captured by summary metrics alone.
#' Facets by S, e, i ,r and plots the occupancy for each of those states for each model .
#' @param parsets the list of parameter sets to compare, which should include the calibrated
#' parameters and model names for each model type.
#' @param palette a vector of colours to use for each model type, of the same length
#' as the number
#' of model types in parsets. Can be generated with make_model_palette.
#' @param show_homlines a boolean indicating whether to show the homogeneous epidemic lines
#' in the bar plots
#' @export
curve_comparison <- function(parsets,
                             palette,
                             show_homlines = FALSE) {
  curve_dfs <- list()
  for (pars in parsets) {
    epi_out <- run_system_wparms(
      parms = pars$calibrated_parms,
      model = seirs_act()
    )
    curve_df <- as.data.frame(epi_out) |>
      dplyr::mutate(model_type = pars$name)
    curve_dfs[[length(curve_dfs) + 1]] <- curve_df
  }
  combined_df <- do.call(rbind, curve_dfs)
  p <- ggplot2::ggplot(
    combined_df,
    ggplot2::aes(x = time, y = value, colour = model_type, group = metric)
  ) +
    ggplot2::facet_wrap(~infection_state) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Time", y = "Infectious Proportion",
      colour = "Model Type"
    ) +
    ggplot2::scale_colour_manual(values = palette)
  print(p)
}

run_comparison <- function(N, r0, num_activity, haslemere_activity_scores) {
  parsets <- make_parameter_sets(
    N = N,
    r0 = r0,
    num_activity = num_activity,
    haslemere_activity_scores = haslemere_activity_scores
  )
  comparison_df <- make_metric_comparison(parsets)
  write.csv(comparison_df, "model_comparison_metrics.csv", row.names = FALSE)
  comparison_df
}

#' Main comparison: faceted SEIR curves (left) + grouped metric bars (right)
#'
#' Builds a combined figure where the left panel shows epidemic trajectories
#' (S, E, I, R) faceted and coloured by model type, and the right panel shows
#' grouped bars for key metrics (HIT, R0, final size, peak I) with matching
#' colours.
#'
#' @param N Total population size
#' @param r0 Target basic reproduction number for calibration
#' @param num_activity Number of activity classes
#' @param haslemere_activity_scores Numeric vector of activity scores (length == num_activity)
#' @param path_to_empirical_mix Path to empirical mixing CSV
#' @param output_dir Directory to save outputs (PDF/PNG); set NULL to skip saving
#' @param save If TRUE, write combined figure to disk
#' @param show_hom_ref If TRUE, overlays the homogeneous expectations (black horizontal lines)
#' @export
compare_curves_main <- function(
  N = 1e6,
  r0 = 2,
  num_activity = 5,
  haslemere_activity_scores = c(
    0.936825,
    2.182864,
    3.492136,
    5.271117,
    9.467058
  ),
  path_to_empirical_mix = "data/Haslemere.csv",
  output_dir = "output/curve_compare",
  save = TRUE,
  show_hom_ref = FALSE
) {
  parsets <- make_parameter_sets(
    N = N,
    r0 = r0,
    num_activity = num_activity,
    haslemere_activity_scores = haslemere_activity_scores,
    path_to_empirical_mix = path_to_empirical_mix
  )

  model_names <- vapply(parsets, function(p) p$name, character(1))
  # Canonical left-to-right order for bars and legend
  model_levels <- c(
    "Homogeneous population",
    "Heterogeneous pop (proportionate mix)",
    "Heterogeneous pop (empirical mix)",
    "With switching"
  )
  palette <- make_model_palette(length(parsets))
  names(palette) <- model_levels

  curve_dfs <- list()
  for (pars in parsets) {
    epi_out <- run_system_wparms(
      parms = pars$calibrated_parms,
      model = seirs_act()
    )
    scaled <- scale_results(epi_out)
    df_tot <- as.data.frame(scaled)[, c(
      "time",
      "s_tot",
      "e_tot",
      "i_tot",
      "r_tot"
    )]
    names(df_tot) <- c("time", "S", "E", "I", "R")
    df_long <- tidyr::pivot_longer(df_tot,
      cols = c("S", "E", "I", "R"),
      names_to = "infec_state", values_to = "value"
    )
    df_long$infec_state <- factor(df_long$infec_state, levels = c("S", "E", "I", "R"))
    df_long$model_type <- pars$name
    curve_dfs[[length(curve_dfs) + 1]] <- df_long
  }
  curves_df <- do.call(rbind, curve_dfs)
  curves_df$model_type <- factor(curves_df$model_type, levels = model_levels)

  p_curves <- ggplot2::ggplot(
    curves_df,
    ggplot2::aes(x = time, y = value, colour = model_type)
  ) +
    ggplot2::geom_line(alpha = 0.9) +
    ggplot2::facet_wrap(~infec_state, ncol = 2, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = bquote("SEIR trajectories (homogeneous population calibrated to " ~ R[0] == .(r0) ~ ")"),
      x = "Time (days)", y = "Fraction"
    ) +
    ggplot2::scale_colour_manual(name = "Model type", values = palette) +
    ggplot2::scale_x_continuous(limits = c(0, 150)) +
    ggplot2::theme(
      strip.text       = ggplot2::element_text(size = 12),
      axis.text        = ggplot2::element_text(size = 11),
      legend.position  = "none"
    )

  metric_rows <- list()
  predicted_from_r0 <- list()

  for (pars in parsets) {
    epi_out <- run_system_wparms(
      parms = pars$calibrated_parms,
      model = seirs_act()
    )
    res <- compute_all(epi_out,
      classed_detail = FALSE,
      params = pars$calibrated_parms
    )
    overall <- as.data.frame(res)
    overall <- overall[overall$class == -1, , drop = FALSE]
    metric_rows[[length(metric_rows) + 1]] <- data.frame(
      model_type = pars$name,
      hit = overall$hit,
      r0 = overall$r0,
      final_size = overall$final_size,
      psI = overall$psI,
      stringsAsFactors = FALSE
    )
  }
  metrics_df <- do.call(rbind, metric_rows)
  metrics_df$model_type <- factor(metrics_df$model_type, levels = model_levels)

  # Compute per-model homogeneous predictions from their respective R0s
  metrics_df_aug <- metrics_df |>
    dplyr::mutate(
      hit_hom_pred = 1 - 1 / r0,
      fs_hom_pred  = vapply(r0, solve_r0_fs, numeric(1))
    )

  metrics_long <- tidyr::pivot_longer(
    metrics_df,
    cols = c(hit, r0, final_size, psI),
    names_to = "metric",
    values_to = "value"
  )

  level_map <- c(
    hit        = "HIT",
    r0         = "R[0]",
    final_size = "Final~size",
    psI        = "Peak~size"
  )

  metrics_long$metric <- factor(
    level_map[metrics_long$metric],
    levels = unname(level_map)
  )

  dodge_w <- 0.8

  p_bars <- ggplot2::ggplot(
    metrics_long,
    ggplot2::aes(x = model_type, y = value, fill = model_type)
  ) +
    ggplot2::geom_bar(
      stat = "identity",
      position = ggplot2::position_dodge(width = dodge_w)
    )

  if (show_hom_ref) {
    # Build reference lines for HIT and final_size facets only.
    ref_lines_long <- dplyr::bind_rows(
      dplyr::transmute(metrics_df_aug,
        model_type,
        metric = "hit",
        hom_pred = hit_hom_pred,
        legend_label = "predicted from R0 via homogeneous eq."
      ),
      dplyr::transmute(metrics_df_aug,
        model_type,
        metric = "final_size",
        hom_pred = fs_hom_pred,
        legend_label = "predicted from R0 via homogeneous eq."
      )
    )
    ref_lines_long$model_type <- factor(ref_lines_long$model_type, levels = model_levels)
    ref_lines_long$metric <- factor(
      level_map[ref_lines_long$metric],
      levels = unname(level_map)
    )

    p_bars <- p_bars +
      ggplot2::geom_errorbar(
        data = ref_lines_long,
        ggplot2::aes(
          x        = model_type,
          ymin     = hom_pred,
          ymax     = hom_pred,
          linetype = legend_label
        ),
        inherit.aes = FALSE,
        colour = "black",
        linewidth = 0.8,
        width = 0.7,
        position = ggplot2::position_dodge(width = dodge_w),
        show.legend = TRUE
      ) +
      ggplot2::scale_linetype_manual(
        name = NULL,
        values = c("predicted from R0 via homogeneous eq." = "solid"),
        guide = ggplot2::guide_legend(
          override.aes = list(colour = "black") # force key to black line
        )
      )
  }

  p_bars <- p_bars +
    # Fix the fill legend so keys show pure colour with no black outline
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(colour = NA)
      )
    ) +
    ggplot2::facet_wrap(
      ~metric,
      scales   = "free_y",
      labeller = ggplot2::label_parsed # <-- replaces the named-vector labeller
    ) +
    ggplot2::scale_fill_manual(name = "Model type", values = palette) +
    ggplot2::theme_bw() +
    ggplot2::labs(x = NULL, y = "Value", fill = "Model type") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 14),
      strip.text = ggplot2::element_text(size = 14),
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 12)
    )

  combined <- patchwork::wrap_plots(p_curves + p_bars) +
    patchwork::plot_annotation(tag_levels = "A")

  if (!is.null(output_dir) && save) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    grDevices::cairo_pdf(
      file.path(
        output_dir,
        "curve_metrics_comparison.pdf"
      ),
      width = 14, height = 6
    )
    print(combined)
    grDevices::dev.off()
    ggplot2::ggsave(
      file.path(
        output_dir,
        "curve_metrics_comparison.png"
      ), combined,
      width = 14, height = 6, dpi = 300
    )
  }
  list(
    curves      = p_curves,
    bars        = p_bars,
    patchwork   = combined,
    metrics_df  = metrics_df,
    metrics_aug = metrics_df_aug # include hit_hom_pred and fs_hom_pred
  )
}
