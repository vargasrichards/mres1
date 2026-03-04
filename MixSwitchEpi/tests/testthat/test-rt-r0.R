library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("Rt equals R0 at t=0 for counts and fractions", {
  pars <- list(
    n_activity = 3,
    class_sizes = c(100, 200, 300),
    activity_scores = c(1, 1, 1),
    beta = 0.2,
    gamma = 0.1,
    epsilon = 0.3,
    swch = matrix(0, nrow = 3, ncol = 3)
  )

  r0 <- MixSwitchEpi::compute_r0(pars)

  # model output with absolute counts (S columns are counts)
  model_counts <- data.table(
    time = c(0, 1),
    S1 = c(pars$class_sizes[1], pars$class_sizes[1] - 1),
    S2 = c(pars$class_sizes[2], pars$class_sizes[2] - 2),
    S3 = c(pars$class_sizes[3], pars$class_sizes[3] - 3)
  )

  res_counts <- MixSwitchEpi::compute_rt_ts(model_counts,
    pars,
    add_to_output = FALSE
  )
  expect_equal(res_counts$R0[1], r0)
  expect_equal(res_counts$Rt[1], r0, tolerance = 1e-10)

  # model output with fractions (S columns are already scaled)
  model_frac <- copy(model_counts)
  model_frac[, S1 := S1 / pars$class_sizes[1]]
  model_frac[, S2 := S2 / pars$class_sizes[2]]
  model_frac[, S3 := S3 / pars$class_sizes[3]]

  res_frac <- MixSwitchEpi::compute_rt_ts(model_frac, pars, add_to_output = FALSE)
  expect_equal(res_frac$R0[1], r0)
  expect_equal(res_frac$Rt[1], r0, tolerance = 1e-10)

  # also test compute_rt for single-popstate inputs (counts and fractions)
  popstate_counts <- list(s_sizes = pars$class_sizes)
  rt_counts_single <- MixSwitchEpi::compute_rt(pars, popstate_counts)
  expect_equal(rt_counts_single, r0, tolerance = 1e-10)

  popstate_frac <- list(s_sizes = rep(1, pars$n_activity))
  rt_frac_single <- MixSwitchEpi::compute_rt(pars, popstate_frac)
  expect_equal(rt_frac_single, r0, tolerance = 1e-10)
})
