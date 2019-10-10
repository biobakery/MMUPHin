test_that("construct_design works", {
  expect_equal(construct_design(NULL), NULL)

  df_test <- data.frame(a = 1:3)
  expect_equal(ncol(construct_design(df_test)), 2)

  df_test <- data.frame(a = factor(1:3))
  expect_equal(ncol(construct_design(df_test)), 3)
})
