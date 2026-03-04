# Minimal working example for heatmap_through_time()
# Save in project root and run with: Rscript scripts/mwe_heatmap_time.R

pkgload::load_all(quiet = TRUE)

message("Building base parameters (calibrates to target R0)")
base_pars <- MixSwitchEpi::make_parms(
  n_activity = 5,
  epsilon = 0,
  switch_rate = 0,
  initial_condition_rule = "uniform",
  infected_fraction = 1e-3,
  activity_scheme = "gamma",
  gamma_dist = list(mean = 4.27, variance = 10.63),
  mixing_model = "garnett_mixing",
  switching_model = "simple_switching",
  pathogen_string = "sars-cov-2",
  target_r0 = 2
)

# Short simulation horizon for a quick example
base_pars$t_end <- 150

# Scenario A: No switching (calibration point)
pars_A <- base_pars

# Scenario B: Fast switching (residence ~1 day)
pars_B <- base_pars
pars_B$swch <- simple_switching(switch_rate = 1, num_classes = pars_B$n_activity)

out_file <- "output/mwe_heatmap.pdf"
message("Running heatmap_through_time() and saving to: ", out_file)

p <- heatmap_through_time(
  parms_list = list("No switching" = pars_A, "Fast switching (res=1d)" = pars_B),
  time_range = c(0, base_pars$t_end),
  save = TRUE,
  output_file = out_file,
  width = 10,
  height = 4
)

message("Done. Preview with: open ", normalizePath(out_file))
