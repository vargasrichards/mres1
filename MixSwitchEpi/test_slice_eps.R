devtools::load_all()

base_pars <- MixSwitchEpi::make_parms(
  n_activity = 5, epsilon = 0, switch_rate = 0, initial_condition_rule = "uniform",
  infected_fraction = 1e-3, activity_scheme = "gamma", gamma_dist = list(mean = 4.27, variance = 10.63),
  mixing_model = "garnett_mixing", switching_model = "simple_switching",
  pathogen_string = "sars-cov-2", target_r0 = 2
)

cat("Running calibration_mode = 'none'...\n")
res_none <- MixSwitchEpi::example_comp(res = 10, target_r0 = 2, pathogen_string = "sars-cov-2", calibration_mode = "none")

cat("Running calibration_mode = 'assort'...\n")
res_assort <- MixSwitchEpi::example_comp(res = 10, target_r0 = 2, pathogen_string = "sars-cov-2", calibration_mode = "assort")

cat("Done!\n")
