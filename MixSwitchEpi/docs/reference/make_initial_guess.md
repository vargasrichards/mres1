# Make an initial guess of an 'effective' mixing matrix

Function takes the mixing and switching matrices to be reproduced and
applies an operator to them to make the 'effective' mixing matrix. That
is, the matrix that if

## Usage

``` r
make_initial_guess(mix_mat, switch_mat, parameterisation)
```

## Arguments

- mix_mat:

  The contacts between different activity classes (and within) those
  activity classes.

- switch_mat:

  Describes the transitions between the different activity classes.

- parameterisation:

  The parameterisation (other than the mix/switch) matrices. So,
  recovery rate, transmission rate, etc.
