# Compute R0 from a static mixing matrix and scalar parameters

Convenience wrapper used in tests to compute R0 when the mixing matrix
is supplied directly (no switching dynamics).

## Usage

``` r
compute_r0_static(mix_mat, scalar_parameters)
```

## Arguments

- mix_mat:

  mixing matrix (n x n) where rows sum to 1

- scalar_parameters:

  list containing at least `beta` and `gamma`. May also contain
  `activity_scores` and `class_sizes`.
