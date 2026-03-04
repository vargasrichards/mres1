# Make the simplest possible switching matrix.

User specifies the number of activity classes to be used and the
switching rate = 1 - rate(remain) in class i, conditional on not
recovering.

## Usage

``` r
simple_switching(switch_rate, num_classes)
```

## Arguments

- switch_rate:

  rate of switching between activity classes

- num_classes:

  number of activity classes

## Value

A square matrix of dimension num_classes x num_classes
