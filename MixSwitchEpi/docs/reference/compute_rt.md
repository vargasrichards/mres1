# Computes the Rt using the population state at a particular time t

Computes the Rt using the population state at a particular time t

## Usage

``` r
compute_rt(pars, popstate)
```

## Arguments

- pars:

  the model parameters

- popstate:

  the state of the population at the time for Rt. Must have a list
  element \$s_sizes giving the number of people who are susceptible to
  infection in each activity class. This will then be divided through by
  the vector of activity class sizes.

## See also

`compute_r0`

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`characterise_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_peaks.md),
[`compute_all()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_all.md),
[`compute_hit()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_hit.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`compute_ttprev()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_ttprev.md),
[`find_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/find_peaks.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
