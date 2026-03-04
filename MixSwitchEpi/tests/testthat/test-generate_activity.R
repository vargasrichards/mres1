test_that("activity_empirical handles valid input correctly", {
  expected_scores <- c(0.3146078, 1.1479220, 2.3361036, 4.2211382, 9.3466695)
  scores <- activity_empirical("haslemere", n_classes = 5)
  expect_equal(scores, expected_scores)
})

test_that("activity_empirical warns when n_classes != 5 for haslemere", {
  expect_warning(activity_empirical("haslemere", n_classes = 4))
})

test_that("activity_empirical handles not implemented cruise setting", {
  expect_error(activity_empirical("cruise", n_classes = 5))
})

test_that("activity_empirical handles unknown setting string", {
  expect_error(activity_empirical("unknown", n_classes = 5), "Unknown setting string")
})
