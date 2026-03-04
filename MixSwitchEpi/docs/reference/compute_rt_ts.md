# Compute R0 and Rt time series from model output

Uses an effective NGM approach, accounting for the depletion of
susceptibles

## Usage

``` r
compute_rt_ts(model_output, pars, add_to_output = TRUE)
```

## Arguments

- model_output:

  Model output from dust/odin

- pars:

  Model parameters

- add_to_output:

  If TRUE, adds R0 and Rt columns to model_output

## Value

If add_to_output is TRUE, returns model_output with R0 and Rt columns
added. Otherwise, returns a data.frame with time, R0, Rt columns.
