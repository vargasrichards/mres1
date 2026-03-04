# Run a standard set of routines

Wrapper for testing combinations of parameters / conditions.

## Usage

``` r
metaroutine(
  n_activity = 5,
  resol = 10,
  gamma_prmt = list(mean = 30, variance = 2500)
)
```

## Arguments

- n_activity:

  the number of activity classes in the model.

- resol:

  the number of values in each direction to scan.

- gamma_prmt:

  the gamma distribution according to which the activity scores are
  calculated.
