# Homogeneous SEIR trajectory (mass-action)

Convenience function that simulates a homogeneous SEIR epidemic with
fixed gamma, sigma and beta chosen to match the requested target R0.

## Usage

``` r
scalar_reference(target_r0, initial_frac = 0.001)
```

## Arguments

- target_r0:

  Desired basic reproduction number R0

- initial_frac:

  Initial exposed fraction (default 1e-3)

## Value

Data.frame with columns time, S, E, I, R
