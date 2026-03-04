test_that("make_homogeneous_parms works as expected", {
  base <- list(
    activity_scores = c(1, 2, 3),
    class_sizes = c(10, 10, 10),
    n_activity = 3
  )
  res <- make_homogeneous_parms(base)
  expect_equal(res$act_vector, 20)
})

test_that(".add_y_plot adds res_time correctly", {
  df <- data.frame(switch_rate = c(1, 0.5, 0))
  res <- .add_y_plot(df)
  expect_equal(res$res_time, c(1, 2, Inf))
  expect_true(!is.null(attr(res, "inf_plot_val")))
})

test_that(".add_y_plot throws error on missing switch_rate", {
  df <- data.frame(a = c(1, 2))
  expect_error(.add_y_plot(df), "scan_df must have switch_rate column")
})
