# Plot epidemic course with annotated metrics

Plots the epidemic trajectory (S, E, I, R over time) and annotates key
epidemic metrics including peak timing/size and final size/time.

## Usage

``` r
plot_epidemic_annotated(
  model_output,
  metrics,
  compartments = c("I", "E"),
  show_total = TRUE,
  annotate_peak = TRUE,
  annotate_final = TRUE,
  title = "Example epidemic with annotated metrics"
)
```

## Arguments

- model_output:

  Neatened model output (from neaten_state)

- metrics:

  Output from compute_all() containing epidemic metrics

- compartments:

  compartments to plot (default: c("S", "E", "I", "R"))

- show_total:

  If TRUE, plots total population dynamics (s_tot, etc.)

- annotate_peak:

  If TRUE, annotate peak timing and size

- annotate_final:

  If TRUE, annotate final epidemic size and time

- title:

  Plot title
