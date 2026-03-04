#!/usr/bin/env Rscript
# 2D epsilon × switching scan using cruise1 activity scores.
# Garnett parametric mixing is built from the cruise1 class-mean contact rates.
# Beta is calibrated to target_r0 at (epsilon=0, switch_rate=0); the same beta
# is held fixed across the grid (calibration_mode = "none").

library(glue)
library(dplyr)
library(data.table, quietly = TRUE)
devtools::load_all()

# ── Configuration ─────────────────────────────────────────────────────────
n_activity   <- 5
target_r0    <- 2
mix_model    <- "garnett_mixing"
switch_model <- "simple_switching"

epsilon_vals   <- seq(0, 1, length.out = 21)
residence_vals <- sort(unique(c(Inf, seq(1, 50, length.out = 30))))
switch_vals    <- ifelse(is.finite(residence_vals), 1 / residence_vals, 0)
sim_t_end      <- 3000

out_dir <- "output/cruise1_epsilon_switch"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load cruise1 data ──────────────────────────────────────────────────────
cruise1_activity_scores <- utils::read.csv(
  "data/cruise1_results/cruise1_activity_scores.csv"
)$activity_score

cruise1_contact_matrix <- as.matrix(utils::read.csv(
  "data/cruise1_results/cruise1_mean_contact_matrix.csv",
  row.names = 1
))

stopifnot(length(cruise1_activity_scores) == n_activity)
stopifnot(all(dim(cruise1_contact_matrix) == n_activity))

# ── Base parameters ────────────────────────────────────────────────────────
N           <- 1e6
class_sizes <- rep(N / n_activity, n_activity)
scalars     <- MixSwitchEpi:::basic_scalar_parms("sars-cov-2")

base_pars <- list(
  beta            = scalars$beta,
  sigma           = scalars$sigma,
  gamma           = scalars$gamma,
  omega           = 0,
  n_activity      = n_activity,
  N               = N,
  t_end           = sim_t_end,
  epsilon         = 0,
  switch_rate     = 0,
  target_r0       = target_r0,
  class_sizes     = class_sizes,
  activity_scores = cruise1_activity_scores
)

base_pars$swch <- MixSwitchEpi::simple_switching(
  switch_rate = 0, num_classes = n_activity
)
base_pars$m <- MixSwitchEpi::garnett_mixing(
  n_activity        = n_activity,
  epsilon_assort    = 0,
  act_vector        = cruise1_activity_scores,
  class_size_vector = class_sizes
)

ic <- MixSwitchEpi::make_initial_conditions(
  n_activity             = n_activity,
  initial_condition_rule = "uniform",
  fraction_exposed       = 1e-3
)
base_pars$E0 <- ic$E0
base_pars$S0 <- ic$S0

base_pars <- MixSwitchEpi::calibrate_parms(
  preliminary_parms = base_pars,
  target_r0         = target_r0
)

# ── Scan ───────────────────────────────────────────────────────────────────
scan_out <- MixSwitchEpi::scan_mixswitch_framework(
  base_pars              = base_pars,
  epsilon_vals           = epsilon_vals,
  switch_vals            = switch_vals,
  classed_detail         = TRUE,
  switch_model           = switch_model,
  mix_model              = mix_model,
  initial_condition_rule = "uniform",
  frac_infected_initial  = 1e-3,
  calibration_mode       = "none",
  t_end                  = sim_t_end
)

results_dt <- data.table::rbindlist(scan_out, fill = TRUE)
results_dt[, residence_time := ifelse(switch_rate == 0, Inf, 1 / switch_rate)]

write.csv(results_dt,
  file.path(out_dir, "cruise1_epsilon_switch_scan_raw.csv"),
  row.names = FALSE
)

# ── Multi-outcome heatmaps ─────────────────────────────────────────────────
results_overall <- results_dt |>
  dplyr::filter(class == -1) |>
  dplyr::mutate(res_time = ifelse(switch_rate == 0, Inf, 1 / switch_rate)) |>
  as.data.frame()

response_vars <- intersect(
  c("hit", "final_size", "end_time", "r0", "psI", "ptI"),
  names(results_overall)
)

# HIT is undefined where the effective Rt (incorporating switching) never
# drops below 1; final_size is well-defined everywhere.
na_hit <- sum(is.na(results_overall$hit))
if (na_hit > 0) message(glue("{na_hit} grid points have HIT = NA (fast-switching regime)"))

MixSwitchEpi::multi_outcome_comparison(
  base_params    = base_pars,
  file           = results_overall,
  ref_r0         = target_r0,
  response_vars  = response_vars,
  save           = TRUE,
  output_dir     = out_dir,
  show_diff      = FALSE,
  width          = 5,
  height         = 4
)

message("Done. Outputs: ", normalizePath(out_dir))
