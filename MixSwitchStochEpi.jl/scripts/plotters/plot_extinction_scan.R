# code for plotting the results of the parallised extinction scan
#
#first lets read in the data

library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyselect)
library(patchwork)
library(forcats)
library(cowplot)

ext_scan <- read.csv("/Users/alexis/Documents/sph_imperial/MixSwitchStochEpi.jl/output/haslemere_stochastic_extinction/combined_extinction_scans.csv")

finite_levels <- sort(unique(ext_scan$residence_time[is.finite(
                     ext_scan$residence_time)]),
                     decreasing = FALSE)
all_levels <- c(as.character(finite_levels), "Inf")

ext_scan <- ext_scan |>
  dplyr::mutate(
    residence_time_f = ifelse(is.infinite(residence_time), 
        "Inf", as.character(residence_time)),
    residence_time_f = factor(residence_time_f, levels = all_levels)
  )

split_by_r0 <- ext_scan |>
  dplyr::group_by(r0) |>
  dplyr::group_split()

#now we have the separated criticality regimes.
#we wish to compute the diff between the high activity and 
#low activity seeds


#now we can construct the plot.
#diff with white as neutral point; red can be lower pr_extinction, blue higher...
#share colourmaps across panels to give a good sense of the scale involved.

# prepare a list to collect per-r0 results
diff_by_r0 <- vector("list", length(split_by_r0))
names(diff_by_r0) <- sapply(split_by_r0, function(df) unique(df$r0))

# compute diff = p_ext(high) - p_ext(low) for each (residence_time, ξ)
for (i in seq_along(split_by_r0)) {
  crit_regime <- split_by_r0[[i]]

  wide <- crit_regime |>
    pivot_wider(
      names_from = seed,
      values_from = p_ext,
      names_prefix = "p_"
    )

  wide <- wide |>
    mutate(
      diff = p_high - p_low
    ) |>
    ungroup()

  # store result (optionally attach r0 value)
  wide$r0 <- unique(crit_regime$r0)
  diff_by_r0[[i]] <- wide
}

diff_all_r0 <- bind_rows(diff_by_r0) |>
  dplyr::mutate(r0_calib = forcats::as_factor(r0))

# print(head(diff_by_r0[[1]]))

write.csv(diff_all_r0, "output/extinction_diff_high_minus_low_by_r0.csv",
row.names = FALSE)


# global range for Pr(extinction) across all seeds and r0
p_min <- min(ext_scan$p_ext, na.rm = TRUE)
p_max <- max(ext_scan$p_ext, na.rm = TRUE)

# global range for diff (high - low)
diff_min <- min(diff_all_r0$diff, na.rm = TRUE)
diff_max <- max(diff_all_r0$diff, na.rm = TRUE)
# make diff scale symmetric around zero for meaningful diverging palette
diff_maxabs <- max(abs(diff_min), abs(diff_max), na.rm = TRUE)

# use mako palette and enlarge colorbar text/keys so labels render clearly in PDF
p_ext_scale <- scale_fill_viridis_c(
  name = expression(Pr~"(extinction)"),
  limits = c(p_min, p_max),
  option = "mako",
  na.value = "grey90",
  guide = guide_colorbar(
    title.theme = element_text(size = 14, face = "bold"),
    label.theme = element_text(size = 12),
    barwidth = unit(0.5, "cm"),
    barheight = unit(3, "cm")
  )
)

# diverging scale symmetric about zero — use plotmath for Delta in the legend title
diff_scale <- scale_fill_gradient2(
  low = "#2166ac", mid = "#f7f7f7", high = "#b2182b", midpoint = 0,
  limits = c(-diff_maxabs, diff_maxabs), na.value = "grey90",
  name = expression(Delta~Pr~"(extinction)"),
  guide = guide_colorbar(
    title.theme = element_text(size = 14, face = "bold"),
    label.theme = element_text(size = 12),
    barwidth = unit(0.5, "cm"),
    barheight = unit(3, "cm")
  )
)

# helper: format residence-time factor labels with up to 1 decimal place (drop trailing .0)
res_time_label <- function(x) {
   sapply(x, function(xx) {
    num <- suppressWarnings(as.numeric(as.character(xx)))
    if (is.na(num) || is.infinite(num)) return(xx)
    if (abs(num - round(num)) < 1e-8) return(as.character(as.integer(round(num))))
    return(format(round(num, 1), nsmall = 1))
  })
}

pdiff_plot <- ggplot(
  diff_all_r0,
  aes(
    x = ε,
    y = residence_time_f,
    fill = diff
  )
) +
  geom_tile() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), labels = res_time_label) +
  diff_scale +
  facet_wrap(~ r0_calib, nrow = 1) +
  theme_bw() +
  labs(
    x = expression(paste("Assortativity, ", epsilon)),
    y = "Expected residence time (days)"
  ) + theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14)
    )

    pbasal_plots <- lapply(split_by_r0, function(df) {
        # use unified p_ext scale for consistency across panels
        ggplot(df, aes(x = ε, y = residence_time_f, fill = p_ext)) +
          geom_tile() +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0), labels = res_time_label) +
         p_ext_scale +
         facet_wrap(~ seed) +
         theme_bw() +
         labs(x = expression(paste("Assortativity, ", epsilon)),
              y = "Expected residence time (days)") +
              theme(axis.text.x = element_text(size = 14),
                     axis.text.y = element_text(size = 14),
                     axis.title.x = element_text(size = 14),
                     axis.title.y = element_text(size = 14),
                     strip.text = element_text(size = 14))
      })


    # list of R0 calibration values (one column per R0)
    unique_r0s <- unique(diff_all_r0$r0)

    # --- build rows (one plot per R0) rather than column composites to avoid horizontal gutters ---
    low_row_plots <- vector("list", length(unique_r0s))
    high_row_plots <- vector("list", length(unique_r0s))
    diff_row_plots <- vector("list", length(unique_r0s))

    for (ii in seq_along(unique_r0s)) {
      r0v <- unique_r0s[ii]
      df_base <- ext_scan |> filter(r0 == r0v)
      df_low <- df_base |> filter(seed == "low")
      df_high <- df_base |> filter(seed == "high")
      df_diff <- diff_all_r0 |> filter(r0 == r0v)

      show_y <- ii == 1

      # base themes (depend on whether this is leftmost column)
      base_theme <- theme_minimal() + theme(legend.position = "none")
      theme_topmid <- base_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

      if (show_y) {
        theme_topmid <- theme_topmid + theme(
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 15),
          axis.ticks.y = element_line(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)
        )
        theme_left_full <- base_theme + theme(
          axis.text.y = element_text(size = 13),
          axis.title.y = element_text(size = 15),
          axis.ticks.y = element_line(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)
        )
      } else {
        theme_topmid <- theme_topmid + theme(
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)
        )
        theme_left_full <- base_theme + theme(
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12)
        )
      }

      # low
      p_low_i <- ggplot(df_low, aes(x = ε, y = residence_time_f, fill = p_ext)) +
        geom_tile() +
        p_ext_scale +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_discrete(labels = res_time_label, expand = c(0,0)) +
        theme_topmid + theme(aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "pt"),
                             axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             panel.background = element_rect(fill = NA, colour = NA), plot.background = element_rect(fill = NA, colour = NA))
      if (show_y) p_low_i <- p_low_i + labs(title = bquote(R[0] == .(r0v)), x = NULL, y = "Residence time") else p_low_i <- p_low_i + labs(title = bquote(R[0] == .(r0v)), x = NULL, y = NULL)

      # high
      p_high_i <- ggplot(df_high, aes(x = ε, y = residence_time_f, fill = p_ext)) +
        geom_tile() +
        p_ext_scale +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_discrete(labels = res_time_label, expand = c(0,0)) +
        theme_topmid + theme(aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "pt"),
                             axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                             panel.background = element_rect(fill = NA, colour = NA), plot.background = element_rect(fill = NA, colour = NA))
      if (show_y) p_high_i <- p_high_i + labs(y = "Residence time")

      # diff (bottom row keeps x axis)
      p_diff_i <- ggplot(df_diff, aes(x = ε, y = residence_time_f, fill = diff)) +
        geom_tile() +
        diff_scale +
        scale_x_continuous(expand = c(0,0)) +
        scale_y_discrete(labels = res_time_label, expand = c(0,0)) +
        labs(x = expression(paste("Assortativity, ", epsilon)), y = NULL) +
        theme_left_full + theme(aspect.ratio = 1, plot.margin = unit(c(0,0,0,0), "pt"),
                                panel.background = element_rect(fill = NA, colour = NA), plot.background = element_rect(fill = NA, colour = NA))
      if (show_y) p_diff_i <- p_diff_i + labs(y = "Residence time")

      low_row_plots[[ii]] <- p_low_i
      high_row_plots[[ii]] <- p_high_i
      diff_row_plots[[ii]] <- p_diff_i
    }

    # assemble main grid as three rows of N columns (this removes inter-column gutters)
    main_grid <- (wrap_plots(low_row_plots, ncol = length(low_row_plots)) /
                  wrap_plots(high_row_plots, ncol = length(high_row_plots)) /
                  wrap_plots(diff_row_plots, ncol = length(diff_row_plots))) + plot_layout(heights = c(1,1,1))

    # left-side label column (Low/High/Diff) - nudge labels rightwards (closer to heatmaps)
    lab_low <- ggplot() + theme_void() + annotate("text", x = 0.99, y = 0.5, label = "Low", size = 5, hjust = 1)
    lab_high <- ggplot() + theme_void() + annotate("text", x = 0.99, y = 0.5, label = "High", size = 5, hjust = 1)
    lab_diff <- ggplot() + theme_void() + annotate("text", x = 0.99, y = 0.5, label = "Diff", size = 5, hjust = 1)
    label_col <- (lab_low / lab_high / lab_diff) + plot_layout(heights = c(1,1,1))

    # build legend grobs for colorbars (used at right side)
    dummy_p <- data.frame(x = 1, y = 1, z = p_min)
    dummy_diff <- data.frame(x = 1, y = 1, z = 0)
    leg_p <- ggplot(dummy_p, aes(x = x, y = y, fill = z)) +
      geom_tile() + p_ext_scale + theme_void() +
      theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))
    leg_diff <- ggplot(dummy_diff, aes(x = x, y = y, fill = z)) +
      geom_tile() + diff_scale + theme_void() +
      theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 12))
    legend_p_grob <- cowplot::get_legend(leg_p)
    legend_diff_grob <- cowplot::get_legend(leg_diff)
    legend_p_plot <- wrap_elements(full = legend_p_grob)
    legend_diff_plot <- wrap_elements(full = legend_diff_grob)

    # right-side legend column
    legend_col <- (plot_spacer() / legend_p_plot / legend_diff_plot) + plot_layout(heights = c(0.15, 1, 0.75))

    # final assembly: labels | main grid | legends; keep horizontal spacing at zero, increase vertical gutter slightly
    final_plot <- label_col | main_grid | legend_col +
                  plot_layout(widths = c(0.03, 1, 0.12)) &
                  theme(panel.spacing.x = unit(0, "lines"), panel.spacing.y = unit(0.25, "lines"), plot.margin = unit(c(0.05,0.05,0.05,0.01), "cm"))

    print(final_plot)
    ggsave("output/extinction_basal_and_diff_grid.png", final_plot, width = 4 * length(unique_r0s) + 2, height = 10, dpi = 300)
    ggsave("output/extinction_basal_and_diff_grid.pdf", final_plot, width = 4 * length(unique_r0s) + 2, height = 10, dpi = 300)
