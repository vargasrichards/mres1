# Changes the β value until the target R0 is met

Uses a very simple approach by dividing the target R0 by the current R0
and scaling β by that factor.

## Usage

``` r
calibrate_parms(preliminary_parms, target_r0)
```

## Arguments

- preliminary_parms:

  list of preliminary parameters including mixing matrix, activity
  scores, class sizes, etc.

- target_r0:

  desired R0 value

## Value

updated parameters with β set to achieve target R0
