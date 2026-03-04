# Function to plot the approximated system against the original system

Allows for immediate visual comparison of the trajectory of the possible
epidemic, both with and without movemement between activity classes.
This function plots the trajectories next to each other.

## Usage

``` r
compare_approx(ref_parms, approx_parms, model)
```

## Arguments

- shared_parms:

  the parameters which are shared by both the approxi- mate system and
  the exact system. Usually this is just the scalar parameters

- original_parms:

  the original parameters: usually nonzero switching matrix and the
  contact matrix.

- approximate_parms:

  the optimised parameters giving similar epidemic trajectories as the
  original parameters

## Details

It also computes the
