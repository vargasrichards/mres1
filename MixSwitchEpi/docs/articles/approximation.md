# MixSwitchEpi

``` r

library(MixSwitchEpi)
library(cli)
library(ggplot2)
message("loaded MixSwitchEpi")
#> loaded MixSwitchEpi
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

It will be useful for our analysis to retrieve some

``` r

library(epiparameter)

epiparameters <- epiparameter_db()
#> Returning 125 results that match the criteria (100 are parameterised). 
#> Use subset to filter by entry variables or single_epiparameter to return a single entry. 
#> To retrieve the citation for each use the 'get_citation' function

influenza_incubation <- epiparameter_db(
  disease = "influenza",
  epi_name = "incubation period",
  single_epiparameter = TRUE
)
#> Using Virlogeux V, Li M, Tsang T, Feng L, Fang V, Jiang H, Wu P, Zheng J, Lau
#> E, Cao Y, Qin Y, Liao Q, Yu H, Cowling B (2015). "Estimating the
#> Distribution of the Incubation Periods of Human Avian Influenza A(H7N9)
#> Virus Infections." _American Journal of Epidemiology_.
#> doi:10.1093/aje/kwv115 <https://doi.org/10.1093/aje/kwv115>.. 
#> To retrieve the citation use the 'get_citation' function

covid_params <- epiparameter_db(
  disease = "COVID-19",
  epi_name = "incubation period", single_epiparameter = T
)
#> Using Linton N, Kobayashi T, Yang Y, Hayashi K, Akhmetzhanov A, Jung S, Yuan
#> B, Kinoshita R, Nishiura H (2020). "Incubation Period and Other
#> Epidemiological Characteristics of 2019 Novel Coronavirus Infections
#> with Right Truncation: A Statistical Analysis of Publicly Available
#> Case Data." _Journal of Clinical Medicine_. doi:10.3390/jcm9020538
#> <https://doi.org/10.3390/jcm9020538>.. 
#> To retrieve the citation use the 'get_citation' function
```

First let’s definine a basic set of parameters with switching
`parms_wswitch`, using the convenience function `default_setup()`.

To illustrate the method, we have chosen a switching rate of 0.25. This
corresponds to there being a

These baseline parms will be the

### Defining contact patterns

Subdivide the population into 5 activity classes. In the most general
case, the activity classes represent the number of risky exposures per
unit time. Clearly, this is

### Performing the static approximations

##### Constrained approximation:

``` r

main_demo <- function(n_activity,
                      switch_rate,
                      initial_condition_rule,
                      pathogen_string,
                      activity_scheme) {
  parms <- make_parms(
    n_activity = n_activity,
    switch_rate = switch_rate,
    initial_condition_rule = initial_condition_rule,
    pathogen_string = pathogen_string, activity_scheme = activity_scheme
  )
  model <- seirs_act()

  reference <- MixSwitchEpi::run_system_wparms(parms = parms, model = model)

  fit_res <- MixSwitchEpi::fit_noswitch(parms, model)
  approx_parms <- fit_res$best_parms

  cat(sprintf(
    "Best epsilon: %0.4f (cost=%0.6g)\n",
    fit_res$best_eps, fit_res$best_cost
  ))

  approx_sim <- MixSwitchEpi::run_system_wparms(
    parms = approx_parms,
    model = model
  )

  Ntot <- parms$N
  print(reference)
  df_ref <- data.frame(time = reference$time, i = reference$i_tot / Ntot, type = "reference")
  df_approx <- data.frame(time = approx_sim$time, i = approx_sim$i_tot / Ntot, type = "approximation")
  df_plot <- rbind(df_ref, df_approx)

  p <- ggplot(df_plot, aes(x = time, y = i, color = type)) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Reference system and no-switch approximation",
      x = "Time", y = "Prevalence (I / N)"
    ) +
    ggplot2::scale_color_manual(values = c("reference" = "#1f78b4", "approximation" = "#e31a1c"))

  ggplot2::ggsave(filename = "noswitch_approx.svg", plot = p, width = 9, height = 4)
}
```

``` r

# main_demo(switch_rate = 0.25,
#           n_activity = 5,
#           initial_condition_rule = "uniform",
#           pathogen_string = "sars-cov-2" )
```

#### Unconstrained approximation

Here we relax the constraint that the contact pattern must be generated
by the Garnett model.

## :::

:::
