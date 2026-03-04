devtools::load_all()
base_pars <- MixSwitchEpi::make_parms(
  n_activity = 5, epsilon = 0, switch_rate = 0, initial_condition_rule = "uniform",
  infected_fraction = 1e-3, activity_scheme = "gamma", gamma_dist = list(mean = 4.27, variance = 10.63),
  mixing_model = "garnett_mixing", switching_model = "simple_switching",
  pathogen_string = "sars-cov-2", target_r0 = 2
)
res_vals <- sort(unique(c(Inf, 1, 5, 50)))
switch_rates <- 1 / res_vals
scan <- MixSwitchEpi::scan_mixswitch_framework(
  base_pars = base_pars, n_epsilon = 3, n_switch = length(switch_rates),
  epsilon_vals = seq(0, 1, length.out = 3), switch_vals = switch_rates,
  classed_detail = FALSE, switch_model = "simple_switching", mix_model = "garnett_mixing",
  initial_condition_rule = "uniform", frac_infected_initial = 1e-3, calibration_mode = "switch"
)
scan_bound <- dplyr::bind_rows(scan)
df <- scan_bound |> dplyr::filter(class == -1)
mat <- tapply(df$r0, list(epsilon = df$epsilon, switch = df$switch_rate), mean)
print("R0 values (rows=epsilon, cols=switch_rate):")
print(mat)
