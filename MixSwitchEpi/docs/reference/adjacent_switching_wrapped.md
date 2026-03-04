# Switching matrix where you can only move between adjacent classes but wraps lowest and highest activities.

Note that we can have min activity -\> max activity in a single step.

## Usage

``` r
adjacent_switching_wrapped(switch_rate, num_classes)
```

## Arguments

- switch_rate:

  switching rate. in general = 1/expected residence time for the
  non-extreme classes (ie, 2,3,4) for a 5-class model.

- num_classes:

  number of activity classes
