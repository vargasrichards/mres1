# Find expected values (means) within gamma distribution bins

Find expected values (means) within gamma distribution bins

## Usage

``` r
find_gamma_expectations(
  cutoffs,
  shape = NULL,
  rate = NULL,
  scale = NULL,
  mean = NULL,
  variance = NULL,
  num_bins = NULL
)
```

## Arguments

- cutoffs:

  Vector of bin boundaries

- shape:

  Shape parameter of gamma distribution (optional if mean/variance
  provided)

- rate:

  Rate parameter (optional if mean/variance provided)

- scale:

  Scale parameter (alternative to rate)

- mean:

  Mean of gamma distribution (alternative to shape/rate)

- variance:

  Variance of gamma distribution (alternative to shape/rate)

- num_bins:

  Number of bins (not used, just for compatibility)

## Value

Vector of conditional expected values for each bin
