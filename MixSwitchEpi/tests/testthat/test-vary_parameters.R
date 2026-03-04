library(testthat)
library(MixSwitchEpi)

test_that("make_initial_uniform produces valid output", {
  ic <- MixSwitchEpi::make_initial_uniform(
    n_activity = 5,
    fraction_exposed = 0.01
  )
  expect_equal(length(ic$S0), 5)
  expect_equal(length(ic$E0), 5)
  expect_true(all(ic$S0 + ic$E0 == 1))
})

test_that("simple_switching has zero row sums", {
  W <- MixSwitchEpi::simple_switching(
    switch_rate = 0.3,
    num_classes = 4
  )
  expect_equal(rowSums(W), rep(0, 4), tolerance = 1e-12)
})
