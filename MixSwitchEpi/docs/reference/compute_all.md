# Wrapper function to compute all metrics

Wrapper function to compute all metrics

## Usage

``` r
compute_all(
  model_output,
  classed_detail = TRUE,
  prevalence_threshold = NULL,
  params
)
```

## Arguments

- model_output:

  Model output (must contain Rt column or params must be provided)

- prevalence_threshold:

  Prevalence threshold for ttprev

- params:

  Model parameters used to generate model output

## See also

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`characterise_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_peaks.md),
[`compute_hit()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_hit.md),
[`compute_rt()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_rt.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`compute_ttprev()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_ttprev.md),
[`find_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/find_peaks.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
