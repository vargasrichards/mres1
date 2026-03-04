# MixSwitchEpi

``` r

library(MixSwitchEpi)
```

MixSwitchEpi provides a simple model for parameterically varying
patterns of interpersonal contact (‘Mix’), and changes include contact
patterns through time (‘switch’). It allows us to explore the effect of
these changes on epidemic dynamics (‘Epi’).

The package provides tools to define and simulate compartmental epidemic
models with varying contact patterns, using odin2 and dust2 for
efficient simulation.

## A simple simulation

We illustrate the central principles of the package using a simple
simulated epidemic in a population.

### Defining contact patterns

Subdivide the population into 5 activity classes. In the most general
case, the activity classes represent the number of risky exposures per
unit time. Clearly, this is

## :::

:::
