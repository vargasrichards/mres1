# Initial conditions where each compartment is inoculated simultaenously with an equal fraction of exposed individuals (or, rather, stuff...)

This is a simple initial condition, Everything else is in the
Susceptible compartment.

## Usage

``` r
make_initial_uniform(n_activity, fraction_exposed)
```

## Arguments

- n_activity:

  the number of activity classes in the model

- fraction_exposed:

  the fraction of each activity class which is exposed to the infection.

## Value

list, whose items are `E0` and `S0`. These arrays are of length
(n_activity).
