# Combined plot of epidemic trajectory and Rt

Creates a two-panel plot showing the epidemic course (I, E) and Rt

## Usage

``` r
plot_epidemic_with_rt(
  model_output,
  computed_metrics,
  compartments = c("I", "E"),
  title = "Epidemic Dynamics and Rt",
  rt_only = FALSE
)
```

## Arguments

- model_output:

  Model output with Rt column (from compute_rt_timeseries)

- computed_metrics:

  Computed epidemic metrics - this allows for

- title:

  the title of the plot

- rt_only:

  whether to plot only the evolution of Rt with time
