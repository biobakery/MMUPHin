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
