
test_that("filter_samples removes sample data as expected", {
  input <- readRDS(test_path("fixtures", "filter_samples_input.rds"))
  output <- readRDS(test_path("fixtures", "filter_samples_output.rds"))

  suppressMessages(expect_equal(filter_samples(input, keeper == T), output))
  suppressMessages(expect_equal(filter_samples(input, sample_ID %notin% paste0("sample", 6:10)), output))
  suppressMessages(expect_equal(filter_samples(input, numeric <= 5), output))
})

test_that("filter_samples checks that input is DAList", {
  expect_error(filter_samples("x"), "must be a DAList object")
})

test_that("filter_samples checks for existence of metadata", {
  input <- readRDS(test_path("fixtures", "filter_samples_input.rds"))
  input["metadata"] <- list(NULL)

  expect_error(filter_samples(input, keeper == T), "does not contain metadata")
})


test_that("filter samples gives error when input and output + filtered samples don't match", {
  input <- readRDS(test_path("fixtures", "filter_samples_input.rds"))
  input$metadata$keeper[2] <- NA

  expect_error(filter_samples(input, keeper == T), "Issue when filtering samples")
})



