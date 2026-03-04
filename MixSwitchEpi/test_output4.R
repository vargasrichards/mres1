devtools::load_all()
base_pars <- MixSwitchEpi::make_parms(
  n_activity = 5, epsilon = 0, switch_rate = 0, initial_condition_rule = "uniform",
  infected_fraction = 1e-3, activity_scheme = "gamma", gamma_dist = list(mean = 4.27, variance = 10.63),
  mixing_model = "garnett_mixing", switching_model = "simple_switching",
  pathogen_string = "sars-cov-2", target_r0 = 2
)

base_pars$beta <- 1.0 # Set beta = 1 to just see the raw lambda

eps_vals <- c(0, 0.5, 1)
sw_vals <- c(0, 0.02, 0.2, 1)

mat <- matrix(0, nrow = 3, ncol = 4)
for (i in 1:3) {
  for (j in 1:4) {
    p <- base_pars
    p$epsilon <- eps_vals[i]
    p$switch_rate <- sw_vals[j]
    p$m <- MixSwitchEpi::garnett_mixing(5, p$epsilon, p$activity_scores, p$class_sizes)
    p$swch <- MixSwitchEpi::simple_switching(p$switch_rate, 5)
    mat[i, j] <- MixSwitchEpi::compute_r0(p)
  }
}
rownames(mat) <- eps_vals
colnames(mat) <- sw_vals
print("Raw Eigenvalues (beta=1):")
print(mat)
