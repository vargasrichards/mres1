# Find epidemic peaks (for I or E + I)

This function requires that the data are scaled.

## Usage

``` r
find_peaks(epi_data, s_label, e_label, i_label, pop_size, class, pars)
```

## Arguments

- epi_data:

  data.table with epidemic trajectories

- s_label:

  Column name for susceptibles

- e_label:

  Column name for exposed

- i_label:

  Column name for infected

- pop_size:

  Total (sub)population size

- class:

  Activity class identifier

## See also

Other epidemic_metrics:
[`characterise_end()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_end.md),
[`characterise_peaks()`](https://vargasrichards.github.io/MixSwitchEpi/reference/characterise_peaks.md),
[`compute_all()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_all.md),
[`compute_hit()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_hit.md),
[`compute_rt()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_rt.md),
[`compute_spectral_gap()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_spectral_gap.md),
[`compute_ttprev()`](https://vargasrichards.github.io/MixSwitchEpi/reference/compute_ttprev.md),
[`scale_results()`](https://vargasrichards.github.io/MixSwitchEpi/reference/scale_results.md)
