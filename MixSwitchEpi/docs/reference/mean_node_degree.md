# Computes the mean degree of a node over the observed time.

Need to be careful that systematically gapped data (like there is in the
Pung et al., 2024 data) does not influence results

## Usage

``` r
mean_node_degree(reconstructed_data, node)
```

## Details

Hence, comparisons should only be made between individuals from same
dataset. Be aware that some time might have all nodal degrees = 0 even
if this is not actually the case.
