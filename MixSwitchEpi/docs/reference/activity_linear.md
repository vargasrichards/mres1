# Generate the activity-class scaler for each activity class Linear activity scores

Simple function which linearly interpolates between the lowest-activity
and the highest-activity class. Hence, for each class i we have a
positive number A_i which gives a measure of activity.

## Usage

``` r
activity_linear(min_act, max_act, n_classes)
```

## Details

(More complex models of how activity changes with the activity class are
possible.)
