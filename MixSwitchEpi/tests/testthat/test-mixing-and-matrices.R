library(testthat)
library(MixSwitchEpi)

test_that("rho_matrix rows sum to 1 and t_matrix non-negative", {
  params <- list(
    n_activity = 4,
    activity_scores = rep(1, 4),
    class_sizes = rep(1, 4) / 4,
    epsilon = 0.3,
    beta = 0.2,
    gamma = 0.1,
    swch = matrix(0, nrow = 4, ncol = 4)
  )

  rho <- MixSwitchEpi::rho_matrix(params)
  expect_true(is.matrix(rho))
  # each row of rho should sum to 1 (mixing probabilities)
  row_sums <- rowSums(rho)
  expect_equal(row_sums, rep(1, params$n_activity), tolerance = 1e-12)

  tmat <- MixSwitchEpi::t_matrix(params)
  expect_true(is.matrix(tmat))
  # transmission rates should be non-negative
  expect_true(all(tmat >= 0))

  sigma <- MixSwitchEpi::sigma_matrix(params)
  expect_true(is.matrix(sigma))
  # sigma diagonal should equal swch diag minus gamma
  expect_equal(diag(sigma), diag(params$swch) - params$gamma)
})


test_that("rho_matrix reduces to uniform matrix when activity scores and class sizes are uniform and epsilon=0", {
  params <- list(
    n_activity = 4,
    activity_scores = rep(1, 4),
    class_sizes = rep(1 / 4, 4),
    epsilon = 0
  )

  rho <- MixSwitchEpi::rho_matrix(params)
  expect_true(is.matrix(rho))
  # all entries should be 1/n_activity in the random mixing case
  expect_equal(rho, matrix(1 / params$n_activity,
    nrow = params$n_activity,
    ncol = params$n_activity
  ), tolerance = 1e-12)
})
