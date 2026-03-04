# Vary both mixing and switching

Vary both mixing and switching

Vary both mixing (epsilon) and switching (rate) in parallel

## Usage

``` r
vary_mixing_switching(
  n_activity,
  switch_rate_range = c(0, 1),
  assort_range = c(0, 1),
  num_switch_points = 10,
  num_assort_points = 10,
  num_cores = 10
)

vary_mixing_switching(
  n_activity,
  switch_rate_range = c(0, 1),
  assort_range = c(0, 1),
  num_switch_points = 10,
  num_assort_points = 10,
  num_cores = 10
)
```

## Arguments

- n_activity:

  Number of activity classes.

- switch_rate_range:

  Numeric length-2, range for switching rate (0..1).

- assort_range:

  Numeric length-2, range for epsilon (0..1).

- num_switch_points:

  Grid resolution for switching.

- num_assort_points:

  Grid resolution for epsilon.

- num_cores:

  Parallel worker count

- seed:

  Optional RNG seed for reproducibility.

## Value

description

Tibble containing results for each (switch_rate, epsilon) pair.
