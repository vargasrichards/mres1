# Multi-outcome 2-D comparison heatmaps

Multi-outcome 2-D comparison heatmaps

## Usage

``` r
multi_outcome_comparison(
  base_params,
  file,
  response_vars,
  ref_switch = NULL,
  ref_r0,
  n_activity = 5,
  output_dir = "output/comparison",
  save = TRUE,
  width = 6,
  height = 4,
  epsilon_levels = NULL,
  res_time_levels = NULL
)
```

## Arguments

- base_params:

  list from
  [`make_parms()`](https://vargasrichards.github.io/MixSwitchEpi/reference/make_parms.md)

- file:

  scan output

- response_vars:

  character vector of metric column names

- ref_r0:

  optional R0

- n_activity:

  the (integer) number of activity classes in the model

- output_dir:

  save location

- save:

  write files?

- width, height:

  size per panel

- ref_epsilon, ref_switch:

  reference point
