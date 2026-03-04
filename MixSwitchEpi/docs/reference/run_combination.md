# Run and record a 2D parameter scan

A wrapper for parameter set construction, calibration, simulation and
writing out results.

## Usage

``` r
run_combination(
  n_activity,
  switch_model,
  mix_model,
  resol,
  gamma_prmt,
  initial,
  calib_eps,
  calib_swch,
  target_r0
)
```

## Arguments

- n_activity:

  the number of activity classes in the model.

- switch_model:

  switching model

- mix_model:

  mixing model

- resol:

  the number of values in each direction to scan. so the total number of
  parameter combinations run will be resol^2.

- gamma_prmt:

  the gamma distribution according to which the activity scores are
  calculated.

- initial:

  initial condition rule. allowed values include "uniform",

- calib_eps:

  the epsilon at which to calibrate to target_r0

- calib_swch:

  the switch rate at which to calibrate t0 target_r0

- target_r0:

  the R0 to calibrate to.
