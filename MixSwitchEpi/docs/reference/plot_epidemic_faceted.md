# Plot all compartments (SEIR) separately

Creates a faceted plot showing S, E, I, R trajectories separately with
annotations on the relevant panels.

## Usage

``` r
plot_epidemic_faceted(
  model_output,
  metrics,
  show_total = TRUE,
  title = "Epidemic Dynamics by Compartment"
)
```
