# Report the parameter scan results

Writes out the summary statistics for the parameter scan, along with the
time at which started and completed. Makes a directory for this.

## Usage

``` r
report_scan_results(results_tibble, parameter_tibble, create_dir = FALSE)
```

## Arguments

- results_tibble:

  a tibble/df containing the results of a parameter scan

- parameter_tibble:

  the tibble containing the parameters used for that scan.

## Details

perhaps should be updated to measure the amount of cpu time the scan
takes?
