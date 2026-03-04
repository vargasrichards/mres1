# plot_pext_switch_scan.R
#
# Reads output/pext_switch_scan.csv produced by scripts/run_pext_switch_scan.jl
# and creates a faceted plot:
#   x = Number of initially exposed (E0)
#   y = Pr(extinction)
#   colour / linetype = init_class (lowest vs highest activity)
#   columns = Expected residence time (1 / switch rate)
#   rows    = R0 (calibrated at ξ=0)
#
# Usage:
#   Rscript scripts/plotters/plot_pext_switch_scan.R

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(cowplot)

# ── read data ──────────────────────────────────────────────────────
basedir <- file.path(
  "/Users/alexis/Documents/sph_imperial/MixSwitchStochEpi.jl"
)
dat <- read.csv(file.path(basedir, "output", "pext_switch_scan_quicker.csv"),
                stringsAsFactors = FALSE)


hom_case <- NULL

write.csv(dat, file.path(basedir, "output", "unified_pext"), row.names = FALSE)
# ── compute expected residence time and label helpers ──────────────
dat <- dat |>
  dplyr::mutate(
    residence_time = ifelse(switch_rate == 0, Inf, 1 / switch_rate),
    init_class_lab = ifelse(
      init_class == min(init_class),
      "Lowest activity",
      "Highest activity"
    ),
    residence_lab = ifelse(
      is.infinite(residence_time),
      "Expec.~residence~time~infinity",
      paste0("Expec.~residence~time~", round(residence_time, 1), "~days")
    ),
    r0_lab = paste0("R[0]~'='~", R0)
  )

# order facets
unique_res <- dat |>
  distinct(switch_rate, residence_time, residence_lab) |>
  arrange(residence_time)

# reverse so that the first level (left-most) is the largest residence time
dat$residence_lab <- factor(dat$residence_lab, levels = rev(unique_res$residence_lab))
dat$r0_lab <- factor(dat$r0_lab,
  levels = paste0("R[0]~'='~", rev(sort(unique(dat$R0))))
)

# Extract homogeneous model rows (from the main CSV or appended hom CSV).
# We'll collapse any duplicate homogeneous runs (e.g., different init_class)
# into a single homogeneous curve per panel by taking the mean p_extinction
# across init_class for each (R0, switch_rate, E0).
hom_rows <- NULL
if ("model_type" %in% names(dat)) {
  hom_rows <- dat |> filter(tolower(as.character(model_type)) == "hom")
}
if (!is.null(hom_case) && nrow(hom_case) > 0) {
  # also include any hom_case rows from the separate file
  hom_rows <- bind_rows(hom_rows, hom_case)
}
if (!is.null(hom_rows) && nrow(hom_rows) > 0) {
  hom_rows <- hom_rows |> mutate(R0 = as.numeric(R0))
  # ensure residence_lab exists on hom_rows (compute from switch_rate)
  hom_rows <- hom_rows |>
    mutate(residence_time = ifelse(switch_rate == 0, Inf, 1 / switch_rate)) |>
    mutate(residence_lab = ifelse(is.infinite(residence_time),
                                  "Expec.~residence~time~infinity",
                                  paste0("Expec.~residence~time~", round(residence_time, 1), "~days")))
  # average across init_class if needed
  hom_case_expanded <- hom_rows |>
    group_by(R0, residence_lab, E0) |>
    summarise(p_extinction = mean(p_extinction, na.rm = TRUE), .groups = "drop") |>
    ungroup() |>
    mutate(
      residence_lab = factor(residence_lab, levels = levels(dat$residence_lab)),
      r0_lab = factor(paste0("R[0]~'='~", R0), levels = levels(dat$r0_lab))
    )
} else {
  hom_case_expanded <- NULL
}

# Ensure numeric types and consistent ordering so geom_line connects points in
# increasing E0 order within each panel. This avoids lines zig-zagging when the
# input CSV is not pre-sorted.
dat <- dat |> mutate(E0 = as.integer(E0),
                     R0 = as.numeric(R0),
                     init_class = as.integer(init_class)) |>
  arrange(R0, residence_lab, init_class, E0)
if (!is.null(hom_case_expanded)) {
  hom_case_expanded <- hom_case_expanded |> mutate(E0 = as.integer(E0)) |>
    arrange(R0, residence_lab, E0)
}

# ── plot ───────────────────────────────────────────────────────────
# Use only inhomogeneous rows for the main coloured curves; homogeneous
# curves will be added separately so they appear as a single legend entry.
dat_inhom <- dat
if ("model_type" %in% names(dat)) {
  dat_inhom <- dat |> filter(tolower(as.character(model_type)) != "hom")
}

p <- ggplot(dat_inhom, aes(
  x = E0,
  y = p_extinction,
  colour = init_class_lab,
  shape = init_class_lab,
  group = init_class_lab
)) +
  geom_hline(aes(yintercept = 0.05, linetype = "Threshold (5%)"),
             colour = "#000000", linewidth = 0.5) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2, alpha = 0.7) +
  # homogeneous benchmark: single black solid line on every panel
  { if (!is.null(hom_case_expanded)) geom_line(data = hom_case_expanded,
    aes(x = E0, y = p_extinction, linetype = "Homogeneous population", group = 1),
    inherit.aes = FALSE, colour = "#000000",
      linewidth = 0.5) else NULL } +
  facet_grid(
    rows = vars(r0_lab),
    cols = vars(residence_lab),
    labeller = label_parsed
  ) +
  scale_colour_manual(
    values = c("Lowest activity"  = "#e41a1c",
               "Highest activity" = "#377eb8"),
    # Swap the label text while keeping colours tied to the same levels
    labels = c("Highest activity (class 5)", "Lowest activity (class 1)")
  ) +
  scale_linetype_manual(
    values = c("Homogeneous population" = "solid",
               "Threshold (5%)" = "dashed"),
    breaks = c("Homogeneous population", "Threshold (5%)"),
    labels = c("Homogeneous population", "Pr(ext.) = 0.05")
  ) +
  scale_shape_manual(
    values = c("Lowest activity" = 4,   # cross
               "Highest activity" = 16), # filled circle
    # Swap the label text while keeping shapes tied to the same levels
    labels = c("Highest activity (class 5)", "Lowest activity (class 1)")
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Number initially Exposed",
    y = expression(Pr(extinction))
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "bottom",
    strip.text       = element_text(size = 16),
    axis.title       = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    # panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "grey70", fill = NA)
  ) +
  guides(
    colour = guide_legend(title = "Initially-exposed class", order = 1),
    shape = guide_legend(title = "Initially-exposed class", order = 1),
    linetype = guide_legend(title = NULL, order = 2)
  )

outdir <- file.path(basedir, "output")
ggsave(file.path(outdir, "pext_switch_scan_quicker.png"), p,
       width = 11, height = 9, dpi = 300)
ggsave(file.path(outdir, "pext_switch_scan_quicker.pdf"), p,
       width = 11, height = 9, device = cairo_pdf)

cat("Saved to", file.path(outdir, "pext_switch_scan_quicker.{png,pdf}"), "\n")
