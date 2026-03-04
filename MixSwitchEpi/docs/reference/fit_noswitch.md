# Fit non-switching system

Produce a fitted time-invariant system which is the 'best' approximation
to a time-varying epidemic system.

## Usage

``` r
fit_noswitch(
  parms,
  model,
  eps_init = NULL,
  metric = c("i_tot", "e_tot"),
  normalize = TRUE,
  verbose = TRUE,
  nl_opts = list(algorithm = "NLOPT_LN_BOBYQA", maxeval = 200, xtol_rel = 1e-06)
)
```

## Arguments

- parms:

  a parameterisation for the modelled system which includes switching
  rate \> 0 (can also be zero for testing)

- model:

  odin model used for simulation

- eps_init:

  initial value of epsilon which is used for the reference system.

- metric:

  the way

## Value

sol solution of the optimisation problem

## Details

This function produces, for a given parameterisation of the activity
class switching matrix, the mixing/contact matrix producing the most
similar results for the Exposed compartment.
