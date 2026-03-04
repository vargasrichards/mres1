library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)

# test_that("peak E+I occurs earlier than peak I for constructed trajectory", {
#   # Construct synthetic epidemic trajectory for a single class
#   dt <- data.table(
#     time = 0:4,
#     E1 = c(0, 5, 1, 0, 0),
#     I1 = c(0, 0, 1, 3, 6)
#   )

#   peaks <- MixSwitchEpi::find_peaks(
#     epi_data = dt,
#     s_label = "S1",
#     e_label = "E1",
#     i_label = "I1",
#     pop_size = 1,
#     class = 1,
#     pars = list()
#   )

#   expect_true(peaks$ptEI < peaks$ptI)
# })

test_that("find_peaks returns identical results for two identical classes", {
  # Two identical classes should give identical peak times and sizes
  dt <- data.table(
    time = 0:4,
    E1 = c(0, 1, 2, 1, 0),
    I1 = c(0, 0, 1, 2, 1),
    E2 = c(0, 1, 2, 1, 0),
    I2 = c(0, 0, 1, 2, 1)
  )

  p1 <- MixSwitchEpi::find_peaks(dt, s_label = "S1", e_label = "E1", i_label = "I1", pop_size = 1, class = 1, pars = list())
  p2 <- MixSwitchEpi::find_peaks(dt, s_label = "S2", e_label = "E2", i_label = "I2", pop_size = 1, class = 2, pars = list())

  expect_equal(p1$ptEI, p2$ptEI)
  expect_equal(p1$ptI, p2$ptI)
  expect_equal(p1$psEI, p2$psEI)
  expect_equal(p1$psI, p2$psI)
})

test_that("scale_results converts counts to fractions using class sizes and pop_tot", {
  # craft a small model_output with absolute counts
  dt <- data.table(
    time = 0:1,
    class1 = c(50, 50),
    class2 = c(50, 50),
    pop_tot = rep(100, 2),
    S1 = c(40, 30),
    E1 = c(5, 10),
    I1 = c(5, 10),
    R1 = c(0, 0),
    S2 = c(45, 35),
    E2 = c(3, 8),
    I2 = c(2, 7),
    R2 = c(0, 0)
  )
  # add aggregate totals expected by scale_results
  dt[, s_tot := S1 + S2]
  dt[, e_tot := E1 + E2]
  dt[, i_tot := I1 + I2]
  dt[, r_tot := R1 + R2]

  scaled <- MixSwitchEpi::scale_results(dt)

  # total columns scaled by pop_tot
  expect_equal(scaled$s_tot[1], dt$s_tot[1] / dt$pop_tot[1], tolerance = 1e-12)
  # per-class S columns scaled by class sizes
  expect_equal(scaled$S1[1], dt$S1[1] / dt$class1[1], tolerance = 1e-12)
  expect_equal(scaled$S2[1], dt$S2[1] / dt$class2[1], tolerance = 1e-12)
})
