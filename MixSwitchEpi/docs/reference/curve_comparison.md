# Compare SEIR curves across models, using the palette to colour the curves by model type. This is a more direct comparison of the epidemic trajectories than the bar_comparison of summary metrics, and can reveal differences in timing and shape of the epidemics that may not be captured by summary metrics alone. Facets by S, e, i ,r and plots the occupancy for each of those states for each model .

Compare SEIR curves across models, using the palette to colour the
curves by model type. This is a more direct comparison of the epidemic
trajectories than the bar_comparison of summary metrics, and can reveal
differences in timing and shape of the epidemics that may not be
captured by summary metrics alone. Facets by S, e, i ,r and plots the
occupancy for each of those states for each model .

## Usage

``` r
curve_comparison(parsets, palette)
```

## Arguments

- parsets:

  the list of parameter sets to compare, which should include the
  calibrated parameters and model names for each model type.

- palette:

  a vector of colours to use for each model type, of the same length as
  the number of model types in parsets. Can be generated with
  make_model_palette.
