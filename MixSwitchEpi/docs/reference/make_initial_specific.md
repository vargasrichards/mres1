# Specific activity class infected first.

initial condition where we only expose one activity class
(user-specified) rather than exposing all activity classes initially.
This has nontrivial consequences for the ensuing dynamics,

## Usage

``` r
make_initial_specific(n_activity, initial_class, fraction_exposed)
```

## Arguments

- n_activity:

  The number of activity classes in the model.

- initial_class:

  the initial class which contains the exposeds

- fraction_exposed:

  The fraction of the activity class which is moved to the Exposed (E)
  compartment.
