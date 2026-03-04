library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("get_info extracts population and class sizes correctly", {
  # create a minimal model_output with class columns and pop_tot
  model_output <- data.table(
    time = 0:1,
    class1 = c(100, 100),
    class2 = c(200, 200),
    class3 = c(300, 300),
    pop_tot = rep(600, 2)
  )

  info <- MixSwitchEpi::get_info(model_output)
  expect_equal(info$pop_total, 600)
  expect_equal(info$num_classes, 3)
  expect_equal(info$class_size_all, c(100, 200, 300))
})

test_that("compute_hit handles Rt never < 1 and scaled/unscaled s_tot", {
  # Case 1: Rt never < 1 -> returns NAs with warning
  dt1 <- data.table(time = 0:400, Rt = rep(1.5, 5), s_tot = rep(500, 5), pop_tot = rep(1000, 5))
  res1 <- MixSwitchEpi::compute_hit(dt1)
  expect_true(is.list(res1))
  expect_true(is.na(res1$HIT))

  # Case 2: s_tot in counts but pop_tot present -> scales correctly
  dt2 <- data.table(time = 0:2, Rt = c(1.2, 0.9, 0.8), s_tot = c(900, 700, 600), pop_tot = rep(1000, 3))
  res2 <- MixSwitchEpi::compute_hit(dt2)
  # At first Rt < 1 is at time=1 with s_tot=700 -> HIT = 1 - 700/1000 = 0.3
  expect_equal(res2$HIT, 0.3, tolerance = 1e-12)
  expect_equal(res2$HIT_time, 1)

  # Case 3: s_tot already scaled
  dt3 <- data.table(time = 0:2, Rt = c(1.2, 0.9, 0.8), s_tot = c(0.9, 0.7, 0.6))
  res3 <- MixSwitchEpi::compute_hit(dt3)
  expect_equal(res3$HIT, 0.3, tolerance = 1e-12)
})
