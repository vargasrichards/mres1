# Compute empirical herd immunity threshold (HIT)

Computes the empirical herd immunity threshold (HIT). This is defined as
the proportion of the population that has been infected when the
effective reproduction number Rt crosses below 1 for the first time.

## Usage

``` r
compute_hit(model_output_rt)
```

## Arguments

- model_output_rt:

  Model output data.table with Rt column already calculated

## Value

A list with elements:

- HIT: the empirical herd immunity threshold

- HIT_time: the time when Rt crosses below 1

- Rt_at_HIT: the value of Rt at the HIT time (should be close to 1)

## See also

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`characterise_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_peaks.md),
[`compute_all()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_all.md),
[`compute_rt()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_rt.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`compute_ttprev()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_ttprev.md),
[`find_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/find_peaks.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
