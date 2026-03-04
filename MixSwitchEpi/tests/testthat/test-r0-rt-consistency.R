library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("compute_r0 equals Rt at t=0 for unscaled and scaled inputs", {
  pars <- list(
    n_activity = 3,
    class_sizes = c(50, 50, 50),
    activity_scores = c(1, 1, 1),
    beta = 0.15,
    gamma = 0.05,
    epsilon = 0.4,
    swch = matrix(0, nrow = 3, ncol = 3)
  )

  r0 <- MixSwitchEpi::compute_r0(pars)

  # single-popstate counts equal to class sizes
  popstate_counts <- list(s_sizes = pars$class_sizes)
  rt_counts <- MixSwitchEpi::compute_rt(pars, popstate_counts)
  expect_equal(rt_counts, r0, tolerance = 1e-12)

  # model_output with S columns as counts
  model_counts <- data.table(
    time = c(0, 1),
    S1 = c(pars$class_sizes[1], pars$class_sizes[1] - 1),
    S2 = c(pars$class_sizes[2], pars$class_sizes[2] - 1),
    S3 = c(pars$class_sizes[3], pars$class_sizes[3] - 1)
  )
  res_counts <- MixSwitchEpi::compute_rt_ts(model_counts, pars, add_to_output = FALSE)
  expect_equal(res_counts$R0[1], r0)
  expect_equal(res_counts$Rt[1], r0, tolerance = 1e-12)

  # model_output with S columns as fractions
  model_frac <- copy(model_counts)
  model_frac[, S1 := S1 / pars$class_sizes[1]]
  model_frac[, S2 := S2 / pars$class_sizes[2]]
  model_frac[, S3 := S3 / pars$class_sizes[3]]
  res_frac <- MixSwitchEpi::compute_rt_ts(model_frac, pars, add_to_output = FALSE)
  expect_equal(res_frac$Rt[1], r0, tolerance = 1e-12)
})
