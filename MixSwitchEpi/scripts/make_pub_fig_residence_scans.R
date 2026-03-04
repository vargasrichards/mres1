#!/usr/bin/env Rscript
# Generate publication-quality 1D residence time scans
# Columns: Epsilon (Assortativity)
# Rows: Target R0
# X-axis: Residence Time
# Metrics: HIT, Final Size, Peak Size (psI), R0
# Color: Calibration Mode (None vs Assortativity)

library(dplyr)
library(ggplot2)
library(tidyr)
devtools::load_all()

# Constants
eps_to_plot <- c(0, 0.2, 0.4, 0.6, 0.8, 1.0)
r0_targets <- c(1.5, 2.0, 3.0)
res_vals <- sort(unique(c(Inf, seq(1, 50, length.out = 30)))) # Increased resolution
switch_rates <- 1 / res_vals

metrics_dict <- c(
  "hit" = "Herd Immunity Threshold (HIT)",
  "final_size" = "Final Size (Attack Rate)",
  "psI" = "Peak Proportion Infectious (psI)",
  "r0" = "Functional R0"
)

run_scans <- function(calib_mode) {
  df_list <- list()
  for (tr0 in r0_targets) {
    message("Running R0 = ", tr0, " with calibration = ", calib_mode)
    base_pars <- MixSwitchEpi::make_parms(
      n_activity = 5, epsilon = 0, switch_rate = 0, initial_condition_rule = "uniform",
      infected_fraction = 1e-3, activity_scheme = "gamma", gamma_dist = list(mean = 4.27, variance = 10.63),
      mixing_model = "garnett_mixing", switching_model = "simple_switching",
      pathogen_string = "sars-cov-2", target_r0 = tr0
    )
    base_pars$t_end <- 4000 # Extended horizon
    
    # Calculate analytical homogeneous R0 given the calibrated structure
    hom_r0_val <- (base_pars$beta * mean(base_pars$activity_scores)) / base_pars$gamma

    scan_res <- MixSwitchEpi::scan_mixswitch_framework(
      base_pars = base_pars,
      n_epsilon = length(eps_to_plot),
      n_switch = length(switch_rates),
      epsilon_vals = eps_to_plot,
      switch_vals = switch_rates,
      classed_detail = FALSE, switch_model = "simple_switching",
      mix_model = "garnett_mixing", initial_condition_rule = "uniform",
      frac_infected_initial = 1e-3, calibration_mode = calib_mode
    )

    scan_df <- dplyr::bind_rows(scan_res) |>
      dplyr::filter(class == -1) |>
      dplyr::mutate(
        target_r0 = tr0,
        hom_r0 = hom_r0_val,        # Fixed-p: same beta, uniform pop at mean activity
        matched_r0 = tr0,            # Matched-R0: mass-action SIR with same observed R0
        res_time = 1 / switch_rate,
        calibration = calib_mode
      )
    df_list[[as.character(tr0)]] <- scan_df
  }
  dplyr::bind_rows(df_list)
}

message("=== SCANNING: Calibration Mode 'none' ===")
df_none <- run_scans("none")

message("=== SCANNING: Calibration Mode 'both' ===")
df_both <- run_scans("both")

df_combined <- bind_rows(df_none, df_both) |>
  mutate(calibration = factor(calibration,
    levels = c("none", "both"),
    labels = c(
      "Fixed p",
      "Recalibrated"
    )
  ))

out_dir <- "output/publication/residence_scans"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plot_pub_combined <- function(df, metric_col, metric_name) {
  df_plot <- df
  # Handle Inf plotting
  finite_rt <- df_plot$res_time[is.finite(df_plot$res_time)]
  max_finite <- max(finite_rt, na.rm = TRUE)
  inf_break_val <- max_finite * 1.15
  df_plot$res_time_plot <- ifelse(is.infinite(df_plot$res_time), inf_break_val, df_plot$res_time)

  df_plot$eps_factor <- factor(df_plot$epsilon, levels = eps_to_plot, labels = paste0("Assortativity ε = ", eps_to_plot))
  df_plot$r0_factor <- factor(df_plot$target_r0, levels = r0_targets, labels = paste0("Target R0 = ", r0_targets))

  # Define color scheme
  calib_palette <- c(
    "Fixed p" = "#1b6ca8",
    "Recalibrated" = "#c0392b"
  )

  # Parse labels for R0 and Epsilon to allow plotmath
  df_plot$eps_factor <- factor(df_plot$epsilon,
    levels = eps_to_plot,
    labels = paste0("epsilon == ", eps_to_plot)
  )

  df_plot$r0_factor <- factor(df_plot$target_r0,
    levels = r0_targets,
    labels = paste0("R[0]^{target} == ", r0_targets)
  )

  # Prepare data: separate finite and infinite residence times
  df_finite <- df_plot %>% filter(is.finite(res_time))
  
  # For infinite residence time (no switching), Fixed p and Recalibrated are identical. Filter to just Fixed p.
  df_inf <- df_plot %>% filter(is.infinite(res_time) & calibration == "Fixed p")

  # Prepare dynamic reference data (one per row/target_r0)
  # Definition 1 — Fixed-p: same beta as structured model, uniform pop at mean activity
  #   hom_r0 = beta * mean(activity) / gamma   (< target_r0 because variance inflates structured R0)
  # Definition 2 — Matched-R0: standard mass-action SIR recalibrated to the same observed R0
  #   matched_r0 = target_r0
  hom_refs <- df_plot |>
    dplyr::select(target_r0, r0_factor, hom_r0, matched_r0) |>
    dplyr::distinct() |>
    dplyr::mutate(
      hom_value_fixedp   = NA_real_,
      hom_value_matchedr0 = NA_real_
    )

  compute_ref_value <- function(r0_vec, metric_col) {
    vapply(r0_vec, function(r0) {
      if (metric_col == "hit") {
        1 - 1 / r0
      } else if (metric_col == "r0") {
        r0
      } else if (metric_col == "final_size") {
        MixSwitchEpi:::solve_r0_fs(r0)
      } else {
        h_mets <- MixSwitchEpi:::homogeneous_metrics(r0)
        if (metric_col %in% names(h_mets)) as.numeric(h_mets[metric_col]) else NA_real_
      }
    }, numeric(1))
  }

  hom_refs$hom_value_fixedp    <- compute_ref_value(hom_refs$hom_r0,    metric_col)
  hom_refs$hom_value_matchedr0 <- compute_ref_value(hom_refs$matched_r0, metric_col)
  hom_refs <- hom_refs |> dplyr::filter(!is.na(hom_value_fixedp))

  # Also handle the metric name if it involves R0
  metric_lab <- if (metric_col == "r0") expression(Functional ~ R[0]) else metric_name

  p <- ggplot() +
    # The main dynamic curves
    geom_line(data = df_finite, aes(x = res_time, y = .data[[metric_col]], color = calibration), linewidth = 1.2, alpha = 0.85) +
    
    # Add the target R0 horizontal line as purely aesthetic reference
    {
      if (metric_col == "r0") geom_hline(data = df_finite, aes(yintercept = target_r0), linetype = "dotted", alpha = 0.5, color = "black")
    } +
    
    # Homogeneous reference lines
    # Red dashed  = Fixed-p (same beta, mean activity — inherently lower R0 than target)
    geom_hline(
      data = hom_refs,
      aes(yintercept = hom_value_fixedp, linetype = "Hom. (fixed p)"),
      colour = "#E41A1C", linewidth = 0.8, alpha = 0.8
    ) +
    # Orange dashed = Matched-R0 (mass-action SIR recalibrated to same observed R0)
    geom_hline(
      data = hom_refs,
      aes(yintercept = hom_value_matchedr0, linetype = "Hom. (matched R0)"),
      colour = "#FF7F00", linewidth = 0.8, alpha = 0.8
    ) +
    
    # No switching reference line
    geom_hline(
      data = df_inf,
      aes(yintercept = .data[[metric_col]], linetype = "Heterogeneous pop. \n no switching"),
      colour = "#4DAF4A", linewidth = 0.8, alpha = 0.8
    ) +
    
    facet_grid(r0_factor ~ eps_factor, scales = "free_y", labeller = label_parsed) +
    theme_minimal(base_size = 18) +
    scale_color_manual(values = calib_palette) +
    
    # Configure the linetype legend — breaks + override colours must be in the same order
    scale_linetype_manual(
      name = "Reference models",
      breaks = c("Hom. (fixed p)", "Hom. (matched R0)", "Heterogeneous pop. \n no switching"),
      values = c(
        "Hom. (fixed p)"                       = "dashed",
        "Hom. (matched R0)"               = "longdash",
        "Heterogeneous pop. \n no switching"   = "dotted"
      ),
      guide = guide_legend(override.aes = list(colour = c("#E41A1C", "#FF7F00", "#4DAF4A")))
    ) +
    
    labs(
      title = metric_lab,
      x = "Expected residence time (days)",
      y = metric_lab,
      color = "Calibration"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
    theme(
      legend.position = "bottom",
      strip.background = element_rect(fill = "grey90", color = "grey80"),
      strip.text = element_text(size = 14, color = "black"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(fill = NA, color = "grey90"),
      plot.title = element_text(size = 26),
      axis.title.x = element_text(margin = margin(t = 12), size = 18),
      axis.title.y = element_text(margin = margin(r = 12), size = 18),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      panel.spacing = unit(1.5, "lines")
    )

  # Added "No switching" annotation for the dashed lines if desired,
  # but color-coded lines are usually cleaner.
  # We can add a small text label or rely on the legend.

  # Output filename logic
  base_name <- if (metric_col == "hit") "HIT_parameteric_contact_patterns" else paste0("publication_", metric_col, "_comparison")
  file_name <- file.path(out_dir, paste0(base_name, ".pdf"))

  ggsave(file_name, p, width = 20, height = 12, device = cairo_pdf)
  message("Saved: ", file_name)
}

for (m_col in names(metrics_dict)) {
  m_name <- metrics_dict[[m_col]]
  plot_pub_combined(df_combined, m_col, m_name)
}

message("All publication figures successfully generated!")
