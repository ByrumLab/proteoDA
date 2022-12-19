
test_that("filter_proteins_contaminants removes contaminants as expected", {
  input <- readRDS(test_path("fixtures", "filter_protein_contam_input.rds"))
  output <- readRDS(test_path("fixtures", "filter_protein_contam_output.rds"))

  expect_equal(filter_proteins_contaminants(input), output)
})


test_that("filter_proteins_by_annotation removes proteins as expected", {

  input <- readRDS(test_path("fixtures", "filter_protein_annotation_input.rds"))
  output_ID <- readRDS(test_path("fixtures", "filter_protein_annotation_output_ID.rds"))
  output_name <- readRDS(test_path("fixtures", "filter_protein_annotation_output_name.rds"))
  output_MW <- readRDS(test_path("fixtures", "filter_protein_annotation_output_MW.rds"))

  expect_equal(filter_proteins_by_annotation(input, protein_ID != "protein10"), output_ID)
  expect_equal(filter_proteins_by_annotation(input, !stringr::str_detect(Protein.Name, "keratin")), output_name)
  expect_equal(filter_proteins_by_annotation(input, molecular_weight < 10), output_MW)
})

test_that("filter_proteins_by_annotation gives error when filter cols contain NA", {
  input <- readRDS(test_path("fixtures", "filter_protein_annotation_input.rds"))
  expect_error(filter_proteins_by_annotation(input, extra > 5), "Issue")

})

test_that("filter_proteins_by_annotation gives error when filtering on nonexisting cols", {
  input <- readRDS(test_path("fixtures", "filter_protein_annotation_input.rds"))
  expect_error(filter_proteins_by_annotation(input, no_column > 5), "not found")
})
