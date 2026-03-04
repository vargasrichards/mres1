# Plot Rt over time with R0 reference line

Creates a plot of Rt over time, showing the epidemic threshold (Rt=1)
and the initial R0 value.

## Usage

``` r
plot_rt_timeseries(model_output, highlight_threshold = TRUE)
```

## Arguments

- model_output:

  Model output with Rt column (from compute_rt_timeseries)

- highlight_threshold:

  If TRUE, highlights region where Rt \> 1
