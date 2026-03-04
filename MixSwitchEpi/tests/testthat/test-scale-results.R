library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("scale_results converts totals and class columns to fractions", {
  # construct a toy output with totals and per-class counts
  model_output <- data.table(
    time = 0:1,
    class1 = c(100, 100),
    class2 = c(200, 200),
    class3 = c(700, 700),
    pop_tot = rep(1000, 2),
    s_tot = c(900, 800),
    e_tot = c(50, 100),
    i_tot = c(50, 100),
    r_tot = c(0, 0),
    S1 = c(90, 80), S2 = c(180, 160), S3 = c(630, 560),
    E1 = 0, E2 = 0, E3 = 0,
    I1 = 0, I2 = 0, I3 = 0,
    R1 = 0, R2 = 0, R3 = 0
  )

  scaled <- MixSwitchEpi::scale_results(model_output)
  # totals scaled by pop_tot
  expect_equal(scaled$s_tot[1], 900 / 1000)
  expect_equal(scaled$e_tot[1], 50 / 1000)
  # per-class S columns scaled by class sizes
  expect_equal(scaled$S1[1], 90 / model_output$class1[1])
  expect_equal(scaled$S2[1], 180 / model_output$class2[1])
  expect_equal(scaled$S3[1], 630 / model_output$class3[1])
})
