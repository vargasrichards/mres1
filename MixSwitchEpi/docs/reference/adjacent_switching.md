# Switching matrix where you can only move between adjacent classes

Hence there is twice expected residence time in the extreme classes.

## Usage

``` r
adjacent_switching(switch_rate, num_classes)
```

## Arguments

- switch_rate:

  switcjing rate. in general = 1/expected residence time for the
  non-extreme classes (ie, 2,3,4) for a 5-class model.

- num_classes:

  number of activity classes

## See also

[`simple_switching`](https://vargasrichards.github.io/MixSwitchEpi/reference/simple_switching.md)
which uses uniform movements where you can move more than distance of 1
class at once.
