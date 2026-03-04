# Characterise epidemic peaks for total and by class

This takes in scaled epidemic output and computes peak sizes.

## Usage

``` r
characterise_peaks(model_output, info, classed_detail = TRUE, pars)
```

## Arguments

- model_output:

  Model output (scaled)

- info:

  Info list from get_info()

- classed_detail:

  If TRUE, compute for every activity class

## See also

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`compute_all()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_all.md),
[`compute_hit()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_hit.md),
[`compute_rt()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_rt.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`compute_ttprev()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_ttprev.md),
[`find_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/find_peaks.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
