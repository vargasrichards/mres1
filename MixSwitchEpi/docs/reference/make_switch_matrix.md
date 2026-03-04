# Make the switching matrix according to supplied model

Wrapper - mostly for concision - which sanitises the switch rate input
into the function and uses the supplied model to produce the CTMC
generator matrix supplied to odin

## Usage

``` r
make_switch_matrix(n_activity, switch_rate, switching_model)
```

## Arguments

- n_activity:

  the number of activity classes in the odin model.

- switch_rate:

  the rate of switching. this = 1/Expected residence time in a
  particular activity class

- switching_model:

  allowed values are `simple_switching` and `adjacent_ switching`.

## Examples

``` r
# Example call
swch <- make_switch_matrix(
  n_activity = 5,
  switch_rate = 0.1,
  switching_model = "simple_switching"
)
#> ℹ using
#> switching model simple_switching at rate 0.1
print(swch)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] -0.100  0.025  0.025  0.025  0.025
#> [2,]  0.025 -0.100  0.025  0.025  0.025
#> [3,]  0.025  0.025 -0.100  0.025  0.025
#> [4,]  0.025  0.025  0.025 -0.100  0.025
#> [5,]  0.025  0.025  0.025  0.025 -0.100
# Example (adjacent wrap-around switching)
swch <- make_switch_matrix(
  n_activity = 5,
  switch_rate = 0.1,
  switching_model = "adjacent_switching"
)
#> ℹ using
#> switching model adjacent_switching at rate 0.1
print(swch)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,] -0.05  0.05  0.00  0.00  0.00
#> [2,]  0.05 -0.10  0.05  0.00  0.00
#> [3,]  0.00  0.05 -0.10  0.05  0.00
#> [4,]  0.00  0.00  0.05 -0.10  0.05
#> [5,]  0.00  0.00  0.00  0.05 -0.05
```
