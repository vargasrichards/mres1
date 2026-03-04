library(testthat)
library(data.table, quietly = TRUE)
library(MixSwitchEpi)
test_that("get_info and scale_results behave as expected", {
  # create a tiny model output in counts
  dt <- data.table(
    time = 0:1,
    class1 = c(50, 50), class2 = c(50, 50),
    S1 = c(50, 40), S2 = c(50, 45),
    E1 = 0, E2 = 0, I1 = 0, I2 = 0, R1 = 0, R2 = 0,
    s_tot = c(100, 85), e_tot = c(0, 0), i_tot = c(0, 0), r_tot = c(0, 0), pop_tot = 100
  )

  info <- MixSwitchEpi::get_info(dt)
  expect_equal(info$pop_total, 100)
  expect_equal(info$num_classes, 2)
  expect_equal(info$class_size_all, c(50, 50))

  scaled <- MixSwitchEpi::scale_results(dt)
  # s_tot should be scaled to fraction of pop
  expect_equal(scaled$s_tot[1], dt$s_tot[1] / info$pop_total)
  # class columns should be scaled by their class sizes
  expect_equal(scaled$S1[1], dt$S1[1] / info$class_size_all[1])
})

test_that("characterise_end and find_peaks return sensible structures", {
  dt <- data.table(
    time = 0:4,
    class1 = rep(50, 5), class2 = rep(50, 5),
    S1 = c(49, 48, 47, 46, 45), E1 = c(0, 1, 0, 0, 0), I1 = c(0, 2, 3, 1, 0), R1 = c(0, 0, 0, 3, 5),
    S2 = c(50, 49, 48, 47, 46), E2 = 0, I2 = c(0, 1, 2, 1, 0), R2 = c(0, 0, 0, 1, 4),
    s_tot = c(99, 97, 95, 93, 91), e_tot = c(0, 2, 0, 0, 0), i_tot = c(0, 3, 5, 2, 0), r_tot = c(0, 0, 0, 4, 9), pop_tot = 100
  )

  ends <- MixSwitchEpi::characterise_end(dt, classed_detail = TRUE)
  expect_true(all(c("class", "start_time", "end_time", "duration", "final_size") %in% names(ends)))

  peaks <- MixSwitchEpi::characterise_peaks(dt, MixSwitchEpi::get_info(dt), classed_detail = TRUE, pars = list())
  expect_true(nrow(peaks) >= 1)
})
