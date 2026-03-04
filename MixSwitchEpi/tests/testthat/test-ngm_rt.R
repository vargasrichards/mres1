library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)

test_that("compute_r0 matches spectral radius of -T %*% inv(Sigma)", {
  pars <- list(
    n_activity = 3,
    class_sizes = c(100, 200, 300),
    activity_scores = c(1, 2, 3),
    beta = 0.15,
    gamma = 0.05,
    epsilon = 0.2,
    swch = matrix(0, 3, 3)
  )

  tmat <- MixSwitchEpi::t_matrix(pars)
  sig <- MixSwitchEpi::sigma_matrix(pars)
  k_manual <- -tmat %*% solve(sig)
  r0_manual <- max(abs(eigen(k_manual, only.values = TRUE)$values))

  r0_fn <- MixSwitchEpi::compute_r0(pars)

  expect_equal(r0_fn, r0_manual, tolerance = 1e-12)
})

test_that("compute_rt returns R0 at t=0 for fully susceptible state", {
  pars <- list(
    n_activity = 2,
    class_sizes = c(50, 50),
    activity_scores = c(1, 1),
    beta = 0.2,
    gamma = 0.1,
    epsilon = 0.3,
    swch = matrix(0, 2, 2)
  )

  r0 <- MixSwitchEpi::compute_r0(pars)

  # model_output with S columns equal to class sizes (counts)
  model_counts <- data.table(time = 0, S1 = pars$class_sizes[1], S2 = pars$class_sizes[2])
  rt_counts <- MixSwitchEpi::compute_rt_ts(model_counts,
    pars,
    add_to_output = FALSE
  )
  expect_equal(rt_counts$R0[1], r0)
  expect_equal(rt_counts$Rt[1], r0, tolerance = 1e-12)

  # model_output with fractions
  model_frac <- data.table(time = 0, S1 = 1, S2 = 1)
  rt_frac <- MixSwitchEpi::compute_rt_ts(model_frac, pars, add_to_output = FALSE)
  expect_equal(rt_frac$R0[1], r0)
  expect_equal(rt_frac$Rt[1], r0, tolerance = 1e-12)
})
