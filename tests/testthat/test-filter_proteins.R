
# filter_protein_contaminants ----------------------------------------------


test_that("filter_proteins_contaminants removes contaminants as expected", {
  input <- readRDS(test_path("fixtures", "filter_protein_contam_input.rds"))
  output <- readRDS(test_path("fixtures", "filter_protein_contam_output.rds"))

  expect_equal(filter_proteins_contaminants(input), output)
})


# filter_proteins_by_annotation --------------------------------------------

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



# filter_proteins_by_groups -----------------------------------------------

test_that("filter_proteins_by_groups gives error when group column is not present", {
  input <- readRDS(test_path("fixtures", "filter_proteins_by_group_input.rds"))

  expect_error(filter_proteins_by_group(input, grouping_column = "xxx", min_reps = 10, min_groups = 3), "not found")
})

test_that("filter_proteins_by_groups gives error when min arguments aren't specifed", {
  input <- readRDS(test_path("fixtures", "filter_proteins_by_group_input.rds"))

  expect_error(filter_proteins_by_group(input), "min_reps")
  expect_error(filter_proteins_by_group(input, min_reps = 2), "min_groups")
})

test_that("filter_proteins_by_groups gives error when group and individual thresholds are too high", {
  input <- readRDS(test_path("fixtures", "filter_proteins_by_group_input.rds"))

  expect_error(filter_proteins_by_group(input, min_reps = 10, min_groups = 3), "Groups below threshold")
  expect_error(filter_proteins_by_group(input, min_reps = 4, min_groups = 4), "Not enough groups")
})

test_that("filter_proteins_by_groups removes proteins as expected", {
  input <- readRDS(test_path("fixtures", "filter_proteins_by_group_input.rds"))
  output_13 <- readRDS(test_path("fixtures", "filter_proteins_by_group_output13.rds"))
  output_21 <- readRDS(test_path("fixtures", "filter_proteins_by_group_output21.rds"))
  output_23 <- readRDS(test_path("fixtures", "filter_proteins_by_group_output23.rds"))
  output_33 <- readRDS(test_path("fixtures", "filter_proteins_by_group_output33.rds"))
  output_43 <- readRDS(test_path("fixtures", "filter_proteins_by_group_output43.rds"))


  expect_equal(filter_proteins_by_group(input, min_reps = 1, min_groups = 3), output_13)
  expect_equal(filter_proteins_by_group(input, min_reps = 2, min_groups = 1), output_21)
  expect_equal(filter_proteins_by_group(input, min_reps = 2, min_groups = 3), output_23)
  expect_equal(filter_proteins_by_group(input, min_reps = 3, min_groups = 3), output_33)
  expect_equal(filter_proteins_by_group(input, min_reps = 4, min_groups = 3), output_43)
})
