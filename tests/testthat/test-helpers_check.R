test_that("check_feature_abd works", {
  mat1 <- matrix(c(0.0, 1.0, 2.0, 3.0), nrow = 2, ncol = 2)
  mat2 <- matrix(c(0, 0, 1, 1), nrow = 2, ncol = 2)
  mat_wrong <- matrix(c(0, 0, 1, 1.1), nrow = 2, ncol = 2)

  expect_equal(check_feature_abd(mat1), "counts")
  expect_equal(check_feature_abd(mat2), "proportions")
  expect_error(check_feature_abd(mat_wrong))
})

test_that("check_samples works", {
  mat_test <- matrix(1, nrow = 2, ncol = 2)
  df_test <- data.frame(a = c(1, 1), row.names = NULL)
  expect_error(check_samples(mat_test, df_test))

  sample_names <- c("1", "2")
  colnames(mat_test) <- sample_names
  expect_equal(check_samples(mat_test, df_test), sample_names)

  mat_test <- cbind(mat_test, "3" = c(1, 1))
  expect_error(check_samples(mat_test, df_test))
})

test_that("check_metadata works", {
  variables_test <- paste0("var", 1:3)
  df_test <- data.frame(var1 = c(1, 2),
                        var2 = c(1, 2),
                        var3 = c(1, NA))

  expect_error(check_metadata(df_test, variables_test))
  expect_error(check_metadata(df_test, c(variables_test, "var4")))
  expect_equal(check_metadata(df_test, variables_test[1:2]), df_test[, 1:2])
  expect_equal(check_metadata(df_test, NULL), NULL)
})

test_that("check_batch works", {
  batch_test <- 1:2

  expect_warning(check_batch(batch_test))
  expect_equal(check_batch(batch_test), factor(1:2))
  expect_error(check_batch(batch_test, min_n_batch = 3))
})
