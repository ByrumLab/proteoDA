
test_that("filter_samples removes sample data as expected", {
  input <- readRDS(test_path("fixtures", "filter_samples_input.rds"))
  output <- readRDS(test_path("fixtures", "filter_samples_output.rds"))

  expect_equal(filter_samples(input, keeper == T), output)
  expect_equal(filter_samples(input, sample_ID %notin% paste0("sample", 6:10)), output)
  expect_equal(filter_samples(input, numeric <= 5), output)
})
