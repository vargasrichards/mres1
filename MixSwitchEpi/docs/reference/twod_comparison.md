# 2D parameter scan & comparison w/ homogeneous

Idea here is that we compare the effect of mixing/switching on the
epidemic dynamics relative to the

## Usage

``` r
twod_comparison(
  base_params,
  file,
  response_var,
  ref_epsilon = NULL,
  ref_switch = NULL,
  ref_r0 = NULL,
  output_dir = "output/comparison",
  save = TRUE,
  width = 8,
  height = 6
)
```

## Arguments

- base_params:

  the basic parameters

- file:

  the path to a file containing the 2d scan data

- response_var:

  the name of a response variable
