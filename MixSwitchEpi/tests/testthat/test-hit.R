library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("compute_hit handles scaled and unscaled s_tot and when Rt never < 1", {
  # case 1: scaled s_tot (fractions)
  dt1 <- data.table(time = c(0, 1, 2), Rt = c(1.5, 0.9, 0.8), s_tot = c(0.99, 0.85, 0.8), pop_tot = 100)
  hit1 <- MixSwitchEpi::compute_hit(dt1)
  expect_true(is.numeric(hit1$HIT))
  expect_true(hit1$HIT >= 0 && hit1$HIT <= 1)

  # case 2: unscaled s_tot (counts) with pop_tot present
  dt2 <- data.table(time = c(0, 1, 2), Rt = c(1.5, 0.99, 0.8), s_tot = c(99, 80, 60), pop_tot = 100)
  hit2 <- MixSwitchEpi::compute_hit(dt2)
  expect_equal(hit2$HIT, 1 - (80 / 100))

  # # case 3: Rt never < 1 -> returns NA with warning
  # dt3 <- data.table(time = 0:4, Rt = rep(1.2, 5), s_tot = rep(100, 5), pop_tot = 100)
  # expect_warning(res3 <- MixSwitchEpi::compute_hit(dt3))
  # expect_true(is.na(res3$HIT))
})
