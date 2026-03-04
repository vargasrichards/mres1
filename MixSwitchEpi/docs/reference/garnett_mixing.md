# parameterically generate mixing matrix for the desired number ofactivity classes

Produces a symmetric matrix of dimension n_activity x n_activity. The
matrix needs to be multiplied by the state vectors S(t) and I(t) in the
transmission term. Since the each entry in matrix m: m i,j is per capita
susceptible class per capita infectious class, then the rows of m sum to
1, making it a stochastic matrix.

## Usage

``` r
garnett_mixing(n_activity, epsilon_assort, act_vector, class_size_vector)
```

## Arguments

- n_activity:

  integer: the number of activity classes to account for

- act_vector:

  a vector containing the scale of the activity classes

- assortativity:

  float: the degree to which classes mix with themselves preferentially.
  when assortativity = 1 the function will return the identity matrix
