# Generates a single entry of the rho matrix

The entry is derived from the assortative mixing model used in (Garrett
et al., 1999)

## Usage

``` r
rho_entry(i, j, params)
```

## Arguments

- i:

  activity class i

- j:

  activity class j

- params:

  the model parameter. Must contain epsilon, activity_scores,
  class_sizes

## Value

rho_ij the (i,j) entry of the rho matrix

## See also

`rho_matrix`
