test_that("check_exposure works", {
  df_test <- data.frame(exposure = as.factor(c(1, 2, 1, 2)),
                        batch = as.factor(c(1, 1, 2, 2)))
  expect_equal(check_exposure(exposure = df_test$exposure,
                              batch = df_test$batch), 
               c("1" = TRUE, "2" = TRUE))
  
  df_test$exposure <- c(1, NA, 1, 2)
  expect_equal(check_exposure(exposure = df_test$exposure,
                              batch = df_test$batch), 
               c("1" = FALSE, "2" = TRUE))
  
  df_test$exposure <- as.factor(c(1, 2, 1, 3))
  expect_error(check_exposure(exposure = df_test$exposure,
                              batch = df_test$batch))
})

test_that("check_covariates works", {
  df_test <- data.frame(a = c(1, 1, 1, 2))
  batch <- as.factor(c(1, 1, 2, 2))
  mat_ind_test <- matrix(c(FALSE, TRUE), nrow = 2, ncol = 1)
  dimnames(mat_ind_test) <- list(c("1", "2"), "a")
  expect_equal(check_covariates(data_covariates = df_test,
                                batch = batch), 
               mat_ind_test)
  
  mat_ind_test <- matrix(nrow = 2, ncol = 0)
  dimnames(mat_ind_test) <- list(c("1", "2"))
  expect_equal(check_covariates(data_covariates = df_test[NULL],
                                batch = batch), 
               mat_ind_test)
})

test_that("check_covariates_random works", {
  df_test <- data.frame(a = c(1, 1, 1, 1, 1, 2),
                        b = c(1, 2, 3, 1, 2, 3))
  batch <- as.factor(c(1, 1, 1, 2, 2, 2))
  mat_ind_test <- matrix(c(FALSE, TRUE, FALSE, FALSE), 
                         nrow = 2, ncol = 2)
  dimnames(mat_ind_test) <- list(c("1", "2"), c("a", "b"))
  expect_equal(check_covariates_random(data_covariates = df_test,
                                       batch = batch), 
               mat_ind_test)
  
  df_test$a <- 1
  expect_error(check_covariates_random(data_covariates = df_test,
                                       batch = batch))
})
