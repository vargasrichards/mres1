# Generate a mixing/contact matrix with selectable model

Wrapper to produce a mixing matrix. Defaults to the Garnett model.

## Usage

``` r
generate_mixingmat(
  n_activity,
  epsilon_assort,
  act_vector,
  class_size_vector,
  mixing_model = c("garnett", "polynomial", "exponential"),
  power_param = NULL
)
```

## Arguments

- n_activity:

  integer number of activity classes

- epsilon_assort:

  assortativity parameter

- act_vector:

  vector of activity scores

- class_size_vector:

  vector of class sizes

- mixing_model:

  character; one of 'garnett', 'polynomial', 'exponential'
