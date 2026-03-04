# Get activity scores for each class according to a gamma distribution with equal-sized activity classes.

Get activity scores for each class according to a gamma distribution
with equal-sized activity classes.

## Usage

``` r
activity_gamma(shape, scale, n_classes)
```

## Arguments

- shape:

  Shape parameter of the gamma distribution

- scale:

  Scale parameter of the gamma distribution

- n_classes:

  Number of activity classes

## Value

A numeric vector of activity scores

## See also

[`activity_linear()`](https://vargasrichards.github.io/MixSwitchEpi/reference/activity_linear.md),
[`activity_null()`](https://vargasrichards.github.io/MixSwitchEpi/reference/activity_null.md)

## Examples

``` r
activity_gamma(shape = 2, scale = 2, n_classes = 5)
#> [1] 1.462099 2.377668 3.356694 4.578563 6.470374
```
