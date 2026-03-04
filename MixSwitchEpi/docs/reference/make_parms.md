# Make a default parameterisaiton

This function makes the parameterisation for mix-switch simulation. N
(total population) is set to 1 million people by default. Switching rate
has default 0.1. Initial conditions are set according to the specified
rule Function assumes that the class sizes are of equal size.

## Usage

``` r
make_parms(
  n_activity,
  epsilon,
  switch_rate,
  initial_condition_rule,
  infected_class = NULL,
  infected_fraction,
  activity_scheme = "gamma",
  gamma_dist = NULL,
  mixing_model,
  switching_model,
  target_r0,
  power_param = NULL,
  pathogen_string = "sars-cov-2"
)
```

## Arguments

- n_activity:

  number of activity classes

- epsilon:

  assortativity parameter for mixing matrix

- switch_rate:

  rate of switching between activity classes

- initial_condition_rule:

  rule for initial conditions. allowed values are `uniformly_infected`
  and `specific_class_infected`

- infected_class:

  if `initial_condition_rule` is `specific_class_infected`, this
  parameter specifies which class is initially infected.

- activity_scheme:

  scheme for activity scores. Allowed values are `linear`, `null` and
  `gamma`. Default is `gamma`. This specifies how activity scores are
  generated for each activity class.

- gamma_dist:

  if `activity_scheme` is `gamma`, this parameter is a list containing
  the parameters for the gamma distribution. It can contain either
  `mean` and `variance`, or `shape` and `rate`. (But not a mixture! :p)

- power_param:

  optional power parameter for the `polynomial` mixing model.

- pathogen_string:

  string describing the pathogen to be modelled. Currently only
  `sars-cov-2` is implemented.

## Value

a list of parameters for the simulation.
