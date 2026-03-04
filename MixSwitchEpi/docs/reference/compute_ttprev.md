# Compute time to threshold prevalence

Compute time to threshold prevalence

## Usage

``` r
compute_ttprev(
  model_output,
  prevalence_threshold,
  num_classes,
  classed_detail = TRUE
)
```

## Arguments

- model_output:

  Output from dust/odin

- prevalence_threshold:

  Numeric threshold for (E + I) / N_j

- num_classes:

  Number of activity classes

- classed_detail:

  Whether to compute per class

## See also

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`characterise_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_peaks.md),
[`compute_all()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_all.md),
[`compute_hit()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_hit.md),
[`compute_rt()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_rt.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`find_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/find_peaks.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
