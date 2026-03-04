# Metrics for homogeneous SEIR epidemic

Computes peak times/sizes, final size, HIT and returns a metrics list
for the homogeneous (mass-action) SEIR epidemic defined by target R0.

## Usage

``` r
homogeneous_metrics(r0)
```

## Arguments

- r0:

  Target basic reproduction number

## Value

A named list with ptI, psI, ptEI, psEI, hit, final_size, end_time, r0
