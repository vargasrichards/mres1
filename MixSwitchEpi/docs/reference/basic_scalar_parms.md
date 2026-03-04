# Get the scalar epi parameters for given pathogen

Convenience function which returns the scalar parameters (ie
transmission rate (beta), rate of progression to I (sigma) and the
recovery rate (gamma).

## Usage

``` r
basic_scalar_parms(pathogen_string = "sars-cov-2")
```

## Arguments

- pathogen_string:

  string describing the desired pathogen for which the parameters are
  required. Default value is "sars-cov-2" for which we use the
  parameterisation from (Britton et al., 2020) Science

## Value

a list of scalars whose elements are `$beta`, `$sigma`, `$gamma`
