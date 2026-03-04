# stochastic modelsuses Gillespie tau-leap Simulation Algorithm to approximately simulate thestochastic dynamics.A. Vargas Richards ICL. Oct. 2025.

stochastic modelsuses Gillespie tau-leap Simulation Algorithm to
approximately simulate thestochastic dynamics.A. Vargas Richards ICL.
Oct. 2025.

## Usage

``` r
switching_helpers(switching_matrix)
```

## Details

rate_SE\[\] \<- lambdai \# might simplify this, not really needed output
\<- (rate_SE) p_SE\[\] \<- 1 - exp(-sigma \* dt) \# array describing the
number of individuals in each activity class transitioning from the

## y2:length(y) \<- Binomial(size - sum(y1:(i - 1)), probi / sum(probi:length(y)))

s_ij, \<- mi, j \* Ii / (Ii + Ei + Ri + Si) lambda\[\] \<- beta \*
sum(s_iji, )

## Susceptible to Exposed classes

p_EI\[\] \<- 1 - exp(-sigma \* dt) p_IR\[\] \<- 1 - exp(-gamma \* dt)
p_RS\[\] \<- 1 - exp(-omega \* dt) \# for extensibility if/when we add
waning immunity.

n_SE\[\] \<- Binomial(Si, p_SEi) n_EI\[\] \<- Binomial(Ei, p_EIi)
n_IR\[\] \<- Binomial(Ii, p_IRi) n_RS\[\] \<- Binomial(Ri, p_RSi)

## first we compute the number exiting each activity class

n_EE\[\] \<- Binomial(Ei, p_EEi) \# these are the transitions between
activity classes n_II\[\] \<- Binomial(Ii, p_IIi) n_RR\[\] \<-
Binomial(Ri, p_RRi) n_SS\[\] \<- Binomial(Si, p_SSi)

switch_s_tmp, \<- swchi, j \* Si \# switch_s\[\] \<- sum(switch_s_tmp,
i) switch_e_tmp, \<- swchi, j \* Ei switch_e\[\] \<- sum(switch_e_tmp,
i) switch_i_tmp, \<- swchi, j \* Ii switch_i\[\] \<- sum(switch_i_tmp,
i) switch_r_tmp, \<- swchi, j \* Ri switch_r\[\] \<- sum(switch_r_tmp,
i)

p_SS\[\] \<- 1 - exp(-switch_si \* dt) p_EE\[\] \<- 1 - exp(-switch_ei
\* dt) p_II\[\] \<- 1 - exp(-switch_ii \* dt) p_RR\[\] \<- 1 -
exp(-switch_ri \* dt)

## n_XX\[\] is the array which details the flux within the X compartment

## due to the switching between activity classes.

## swch = no_switching

swch = sswch ) A semi-stochastic model with deterministic activity-class
switching and stochastic transmission events
