# Calculate R0 across a grid of epsilon and switching values

Calculate R0 across a grid of epsilon and switching values

## Usage

``` r
compute_r0_grid(
  epsilon_vals,
  switch_vals,
  base_pars,
  switch_rule,
  increment_scalar = NULL
)
```

## Arguments

- epsilon_vals:

  Vector of assortativity values

- switch_vals:

  Vector of switching rates

- base_pars:

  Base parameter list (modified for each combination)
