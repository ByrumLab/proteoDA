test_that("missing_to_zero converts missing values as expected", {

  input <- readRDS(test_path("fixtures", "missing_to_zero_input.rds"))
  output_NA <- readRDS(test_path("fixtures", "missing_to_zero_output_NA.rds"))

  expect_equal(missing_to_zero(input), output_NA)
})


test_that("missing_to_zero converts custom missing values as expected", {

  input <- readRDS(test_path("fixtures", "missing_to_zero_input.rds"))
  output_neg9 <- readRDS(test_path("fixtures", "missing_to_zero_output_neg9.rds"))
  output_both <- readRDS(test_path("fixtures", "missing_to_zero_output_both.rds"))

  expect_equal(missing_to_zero(input, missing_values = c(-9)), output_neg9)
  expect_equal(missing_to_zero(input, missing_values = c(NA, -9)), output_both)

})

