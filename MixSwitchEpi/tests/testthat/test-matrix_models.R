library(testthat)
library(MixSwitchEpi)

test_that("generate_mixingmat produces stochastic rows and reduces to homogeneous case", {
  n <- 4
  act <- rep(1, n)
  class_sizes <- rep(1 / n, n)
  eps <- 0.0

  m <- MixSwitchEpi::garnett_mixing(
    n_activity = n,
    epsilon_assort = eps,
    act_vector = act,
    class_size_vector = class_sizes
  )

  expect_true(is.matrix(m))
  # rows sum to 1
  expect_equal(rowSums(m), rep(1, n), tolerance = 1e-12)
  # with epsilon = 0 and homogeneous activity/class sizes, matrix should be constant rows
  expect_true(all(abs(m - matrix(1 / n, nrow = n, ncol = n)) < 1e-12))
})

test_that("activity generators behave as expected", {
  lin <- MixSwitchEpi::activity_linear(1, 5, 5)
  expect_equal(length(lin), 5)
  expect_equal(lin[1], 1)
  expect_equal(lin[5], 5)

  nul <- MixSwitchEpi::activity_null(3)
  expect_equal(nul, rep(1, 3))
})

test_that("simple_switching returns matrix with zero row sums and correct diagonal", {
  n <- 5
  rate <- 0.2
  W <- MixSwitchEpi::simple_switching(switch_rate = rate, num_classes = n)
  expect_true(is.matrix(W))
  expect_equal(dim(W), c(n, n))
  # row sums should be zero (outflows sum to zero)
  expect_equal(rowSums(W), rep(0, n), tolerance = 1e-12)
  # diagonal entries equal -rate
  expect_equal(diag(W), rep(-rate, n))
})
