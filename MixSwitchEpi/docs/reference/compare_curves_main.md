# Main comparison: faceted SEIR curves (left) + grouped metric bars (right)

Builds a combined figure where the left panel shows epidemic
trajectories (S, E, I, R) faceted and coloured by model type, and the
right panel shows grouped bars for key metrics (HIT, R0, final size,
peak I) with matching colours. Returns a list with the individual plots,
combined patchwork, and the underlying metrics dataframe.

## Usage

``` r
compare_curves_main(
  N = 1e+06,
  r0 = 2,
  num_activity = 5,
  haslemere_activity_scores = c(0.936825, 2.182864, 3.492136, 5.271117, 9.467058),
  calibrate_all = TRUE,
  output_dir = "output/curve_compare",
  save = TRUE
)
```

## Arguments

- N:

  Total population size

- r0:

  Target basic reproduction number for calibration

- num_activity:

  Number of activity classes

- haslemere_activity_scores:

  Numeric vector of activity scores (length == num_activity)

- calibrate_all:

  If TRUE, calibrates each model to target r0; otherwise copies
  calibrated beta

- output_dir:

  Directory to save outputs (PDF/PNG); set NULL to skip saving

- save:

  If TRUE, write combined figure to disk
