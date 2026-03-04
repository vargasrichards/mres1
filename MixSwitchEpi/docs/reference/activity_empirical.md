# Retrieve fitted activity scores for a given setting

This function retrieves the fitted activity scores for a specific
setting based on empirical data from (Pung et al., 2024) Interface

## Usage

``` r
activity_empirical(setting_string, n_classes)
```

## Arguments

- setting_string:

  A string indicating the setting

- n_classes:

  The number of activity classes

## Value

A numeric vector of activity scores

## Examples

``` r
activity_empirical("haslemere", n_classes = 5)
#> [1] 0.3146078 1.1479220 2.3361036 4.2211382 9.3466695
```
