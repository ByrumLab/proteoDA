

# external DAList function ------------------------------------------------
# Just test for the errors
test_that("DAList() gives error when there is no uniprot_id column", {
  input <- readRDS(test_path("fixtures", "s3-class_input.rds"))

  colnames(input$annotation) <- "unipro_xx"

  # Assemble the DAList from existing raw data.
  expect_error(DAList(data = input$data,
                      annotation = input$annotation,
                      metadata = input$metadata),
               "uniprot_id")
})


test_that("DAList() gives error when uniprot_id column is not unique", {
  input <- readRDS(test_path("fixtures", "s3-class_input.rds"))

  input$annotation$uniprot_id[1] <- input$annotation$uniprot_id[2]

  # Assemble the DAList from existing raw data.
  expect_error(DAList(data = input$data,
                      annotation = input$annotation,
                      metadata = input$metadata),
               "not unique")

})

# Internal validator ------------------------------------------------------


