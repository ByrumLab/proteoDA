test_that("perseus_impute works on a simple matrix and preserves all-NA rows", {
  set.seed(42)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4)
  X[c(1,3), 2] <- NA_real_
  X[5, ] <- NA_real_  # entire row missing
  
  X_imp <- perseus_impute(X, seed = 1)
  expect_true(is.matrix(X_imp))
  # All-NA row should remain all NA
  expect_true(all(is.na(X_imp[5, ])))
  # Other NAs should be imputed
  expect_false(any(is.na(X_imp[-5, ])))
  # Mask present by default
  mask <- attr(X_imp, "imputed_mask")
  expect_true(is.matrix(mask))
  expect_true(mask[1,2] && mask[3,2])
  expect_true(!any(mask[5, ]))  # no mask for rows left NA
})

test_that("perseus_impute returns reproducible results with seed", {
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  X[1, c(2,3)] <- NA_real_
  
  A <- perseus_impute(X, seed = 99)
  B <- perseus_impute(X, seed = 99)
  expect_equal(A, B)
})

test_that("perseus_impute saves before/after in matrix mode when requested", {
  set.seed(7)
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  X[2, 1] <- NA_real_
  
  X_imp <- perseus_impute(X, seed = 5, save_before_after = TRUE)
  expect_equal(attr(X_imp, "before_log2"), X)
  expect_true(is.matrix(attr(X_imp, "imputed_mask")))
})

test_that("perseus_impute writes per-contrast data and saves before/after", {
  # Minimal DAList w/ per-contrast data
  X <- matrix(rnorm(12), nrow = 3, ncol = 4,
              dimnames = list(paste0("p",1:3), paste0("s",1:4)))
  X[c(1,3),2] <- NA_real_
  
  DA <- list(
    data = X,  # global fallback source
    data_per_contrast = list(
      ContrastA = list(log2 = X),
      ContrastB = list(log2 = X)
    )
  )
  
  out <- perseus_impute(DA, seed = 11, save_before_after = TRUE)
  
  # Per-contrast matrices exist and are imputed
  expect_true(is.matrix(out$data_per_contrast$ContrastA))
  expect_true(is.matrix(out$data_per_contrast$ContrastB))
  expect_false(any(is.na(out$data_per_contrast$ContrastA[-3, ])))  # row 3 not entirely missing here
  
  # Before/after + mask saved in imputation_per_contrast
  expect_true(is.matrix(out$imputation_per_contrast$ContrastA$before_log2))
  expect_true(is.matrix(out$imputation_per_contrast$ContrastA$after_log2))
  expect_true(is.matrix(out$imputation_per_contrast$ContrastA$imputed_mask))
})

test_that("perseus_impute falls back to global data when no per-contrast slot", {
  X <- matrix(rnorm(12), nrow = 3, ncol = 4)
  X[1,2] <- NA_real_
  
  DA <- list(data = X)
  out <- perseus_impute(DA, seed = 3, save_before_after = TRUE)
  expect_true(is.matrix(out$data))
  expect_false(any(is.na(out$data[-which(rowSums(is.na(X)) == ncol(X)), ])))
  expect_true(is.matrix(out$imputation$before_log2))
  expect_true(is.matrix(out$imputation$after_log2))
  expect_true(is.matrix(out$imputation$imputed_mask))
})

test_that("plot_perseus_imputation returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  X <- matrix(rnorm(20), nrow = 5, ncol = 4,
              dimnames = list(NULL, paste0("S",1:4)))
  X[c(1,3),2] <- NA_real_
  
  X_imp <- perseus_impute(X, seed = 2)
  p <- plot_perseus_imputation(X, X_imp, bins = 10, facet_ncol = 2)
  expect_s3_class(p, "ggplot")
})
