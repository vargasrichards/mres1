# Plot 2D heatmaps of epidemic metrics

Plot 2D heatmaps of epidemic metrics

## Usage

``` r
plot_2d_heatmaps(scan_results, metrics_to_plot = NULL, overall_only = TRUE)
```

## Arguments

- scan_results:

  Results from scan_2d_parameter_space()

- metrics_to_plot:

  Vector of metric names to plot (default: all)

- overall_only:

  If TRUE, only plot overall metrics (class == -1)

## Value

List of ggplot objects
