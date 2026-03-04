# Plot the grid of (e.g.,) Haslemere results

Plots a grid where rows define the Calibration R0 and columns define the
metric. Manually constructs rows to ensure consistent Y-scales for each
metric across R0 values (columns) while allowing scales to differ
between metrics. Includes homogeneous reference lines with labels.

## Usage

``` r
plot_haslemere_grid(grid_results, gamma = 1/4)
```

## Arguments

- grid_results:

  a dataframe of results from haslemere_grid
