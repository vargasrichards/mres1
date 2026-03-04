# Scan residence times (Inf, 20, 19, ..., 1) across assortativity and produce multi-outcome comparison
library(glue)
library(data.table, quietly = TRUE)

n_activity <- 5
mix_model <- "garnett_mixing"
switch_model <- "simple_switching"
resol_epsilon <- 21
epsilon_vals <- seq(0, 1, length.out = resol_epsilon)
residence_vals <- c(Inf, 20:1) # Inf, 20d,19d,...,1d
switch_vals <- ifelse(is.finite(residence_vals), 1 / residence_vals, 0)

# make base parameters and run scan
base_pars <- MixSwitchEpi::make_parms(
  n_activity = n_activity,
  epsilon = 0.0,
  switch_rate = 0.0,
  initial_condition_rule = "uniform",
  infected_fraction = 1e-03,
  activity_scheme = "gamma",
  gamma_dist = list(mean = 5, variance = 10),
  mixing_model = mix_model,
  switching_model = switch_model,
  pathogen_string = "sars-cov-2",
  target_r0 = 2
)

out <- MixSwitchEpi::scan_mixswitch_framework(
  base_pars = base_pars,
  epsilon_vals = epsilon_vals,
  switch_vals = switch_vals,
  classed_detail = TRUE,
  switch_model = switch_model,
  mix_model = mix_model
)

results_dt <- data.table::rbindlist(out, fill = TRUE)
results_dt <- results_dt[, residence_time := ifelse(switch_rate == 0, Inf, 1 / switch_rate)]

# save
ofile <- glue::glue("scan_residence_series_{mix_model}_{switch_model}.csv")
write.csv(results_dt, file = ofile, row.names = FALSE)
message("Saved scan results to ", ofile)

results_dt <- results_dt |>
  dplyr::filter(class == -1)

# pick numeric response variables to plot: exclude parameter columns
exclude_cols <- c(
  "spectral_gap", "rt_at_hit", "hit_time",
  "epsilon", "switch_rate",
  "residence_time", "class",
  "start_time", "end_time", "duration", "X"
)
num_cols <- names(results_dt)[sapply(results_dt, is.numeric)]
response_vars <- setdiff(num_cols, exclude_cols)

# produce comparison for all numeric outcomes
p <- MixSwitchEpi:::multi_outcome_comparison(
  base_params = base_pars,
  file = results_dt,
  response_vars = response_vars,
  ref_epsilon = NULL,
  ref_switch = NULL,
  ref_r0 = NULL,
  output_dir = "output/comparison",
  save = TRUE,
  width = 5,
  height = 4
)

message("Multi-outcome comparison produced and saved.")
