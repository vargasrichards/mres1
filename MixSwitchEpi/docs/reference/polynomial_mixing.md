# Parametrically generate mixing matrix for the desired

Note that this has an extra parameter compared to the other mixing
models. Thus the call signature needs to be adjusted accordingly

## Usage

``` r
polynomial_mixing(
  n_activity,
  power_param,
  epsilon_assort,
  act_vector,
  class_size_vector
)
```

## Arguments

- n_activity:

  the number of activity classes in the model

- power_param:

  the power parameter. if = 1

- epsilon_assort:

  the assortativity.

- act_vector:

  the vector of activity class scores. Gives a measure of the activity
  of each class.

- class_size_vector:

  the vector of activity class sizes
