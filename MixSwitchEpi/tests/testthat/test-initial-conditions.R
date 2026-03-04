library(testthat)
library(MixSwitchEpi)

test_that("make_initial_conditions scales E0 and S0 to class sizes
    for uniform exposure", {
  n_activity <- 5
  fraction_exposed <- 1e-4
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)

  ic <- MixSwitchEpi::make_initial_conditions(
    n_activity = n_activity,
    initial_condition_rule = "uniform",
    fraction_exposed = fraction_exposed
  )

  expect_true(is.numeric(ic$E0))
  expect_true(is.numeric(ic$S0))
  # total exposed should equal fraction_exposed * N
  total_exposed <- sum(ic$E0 * class_sizes)
  expect_equal(total_exposed, fraction_exposed * N)
  # each class susceptible fraction should be 1 - fraction_exposed
  expect_equal(ic$S0, rep(1 - fraction_exposed, n_activity))
})

test_that("make_initial_conditions specific class exposure behaves correctly", {
  n_activity <- 5
  initial_class <- 3
  fraction_exposed <- 1e-3
  N <- 1e6
  class_sizes <- rep(N / n_activity, n_activity)

  ic <- MixSwitchEpi::make_initial_conditions(
    n_activity = n_activity,
    initial_condition_rule = "specific",
    fraction_exposed = fraction_exposed,
    infected_class = initial_class
  )

  # exposed should only be in the specified class
  total_exposed <- sum(ic$E0 * class_sizes)
  expected_exposed <- fraction_exposed * class_sizes[initial_class]
  expect_equal(total_exposed, expected_exposed)
  # the infected class susceptible fraction should be 1 - fraction_exposed
  expect_equal(ic$S0[initial_class], 1 - fraction_exposed)
  # all other classes should be fully susceptible
  expect_true(all(ic$S0[-initial_class] == 1))
})
