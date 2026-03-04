# Score difference between model outputs

Produces a score of the difference in model result between the original
mixing/switching parameterisation (`true_matrices`) and the candidate
mixing matrix.

## Usage

``` r
score_difference(x, reference_pattern, shared_parms, model, similar_matrices)
```

## Arguments

- x:

  the vector used by the optimiser to represent the mixing matrix.

- reference_pattern:

  gives the desired pattern of e.g., E(t) through time, which we're
  going to try and match as close as poss. with no switching.

- shared_parms:

  the shared parameters (e.g., transmission rate, recovery rate).

- model:

  the odin/monty model

- similar_matrices:

  bool. Indicates whether to reward effective contact matrices which are
  similar to the original contact matrix. Experimental feature.
