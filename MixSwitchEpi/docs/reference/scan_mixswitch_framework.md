# Scan 2D parameter space (epsilon × switch_rate.

Basic flexible framework function which can be adapted significantly as
we require. Runs a set of workers in parallel. We precompute the mixing
and switching matrices to reduce computational time.

## Usage

``` r
scan_mixswitch_framework(
  base_pars,
  epsilon_vals,
  switch_vals,
  n_epsilon = 11,
  n_switch = 11,
  classed_detail = FALSE,
  switch_model = "simple_switching",
  mix_model,
  power_param = NULL,
  initial_condition_rule = NULL,
  frac_infected_initial = 1e-04
)
```

## Arguments

- base_pars:

  Calibrated base parameters from default_setup().

- epsilon_vals:

  Vector of epsilon values to scan (default: 0 to 1)

- switch_vals:

  Vector of switching rate values to scan (default: 0 to 1)

- n_epsilon:

  Number of epsilon points (if epsilon_vals not provided)

- n_switch:

  Number of switch rate points (if switch_vals not provided)

- classed_detail:

  Whether to compute class-specific metrics

- switch_model:

  the model of switching to be used. allowed values include:

- mix_model:

  the model for

- initial_condition_rule:

  the rule to generate the initial cond.

- frac_infected_initial:

  the initial fraction infected of the population This will be put in
  the E (exposed) compartment.

## Value

data.table with all epidemic metrics for each parameter combination
