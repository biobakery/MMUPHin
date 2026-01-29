test_that("set_pseudo works", {
  expect_error(set_pseudo(NA))
  expect_error(set_pseudo(-1))
  expect_error(set_pseudo(0))
  expect_equal(set_pseudo(1), 0.5)
})

test_that("match_control works", {
  default_test <- list("a" = 1, "b" = 1)
  control_test1 <- list("b" = 2)
  control_test2 <- list("b" = 2, "c" = 2)
  control_modified <- list("a" = 1, "b" = 2)

  expect_equal(match_control(default_test), default_test)
  expect_equal(match_control(default_test, control_test1), control_modified)
  expect_equal(match_control(default_test, control_test2), control_modified)
})

test_that("normalize_features works", {
  features_test <- matrix(1:4, 2, 2)

  expect_equal(normalize_features(features_test), features_test)
  features_norm <- normalize_features(features_test,
                                      normalization = "TSS",
                                      pseudo_count = 0.5)
  expect_equal(all(features_norm > 0), TRUE)
  expect_equal(all(apply(features_norm, 2, sum) == 1), TRUE)
})

test_that("transform_features works", {
  features_test <- matrix(1:4, 2, 2)

  expect_error(transform_features(features_test))
  features_norm <- normalize_features(features_test,
                                      normalization = "TSS",
                                      pseudo_count = 0.5)
  expect_equal(all(transform_features(features_norm, "LOG") < 0), TRUE)
})

test_that("fill_dimnames works", {
  x_test <- matrix(1:6, 2, 3)

  expect_error(fill_dimanes(x_test))

  x_filled <- fill_dimnames(x_test,
                            row_prefix = "a",
                            col_prefix = "b")
  expect_equal(rownames(x_filled), paste0("a", 1:2))
  expect_equal(colnames(x_filled), paste0("b", 1:3))

  # data frame by default always have dimnames so fill_dimnames won't do
  # anything
  expect_equal(fill_dimnames(as.data.frame(x_test), "a", "b"),
               as.data.frame(x_test))
})

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