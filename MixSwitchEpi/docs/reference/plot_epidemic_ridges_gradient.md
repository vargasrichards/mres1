# Ridgeline-style plot with time gradient and adjustable spacing

Ridgeline-style plot with time gradient and adjustable spacing

## Usage

``` r
plot_epidemic_ridges_gradient(model_output, compartment = "I", spacing = 1)
```

## Arguments

- model_output:

  Output from the model

- compartment:

  Which compartment to plot: "I", "E", or "EI"

- spacing:

  Vertical spacing between classes. Default 1 = no overlap. Use \< 1 for
  more overlap (e.g., 0.5 = 50% overlap), \> 1 for more separation
