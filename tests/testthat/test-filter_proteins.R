
test_that("filter_proteins_contaminants removes contaminants as expected", {
  input <- readRDS(test_path("fixtures", "filter_protein_contam_input.rds"))
  output <- readRDS(test_path("fixtures", "filter_protein_contam_output.rds"))

  expect_equal(filter_proteins_contaminants(input), output)
})
