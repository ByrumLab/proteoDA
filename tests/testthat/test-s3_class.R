

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

test_that("validate_DAList checks for proper number and order of slots", {

  expect_error(validate_DAList(list(data = "data")),
               "missing the following slots: annotation")

  input <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = input$data,
                  annotation = input$annotation,
                  metadata = input$metadata)
  extra <- c(input, xxx = "xxx")
  expect_error(validate_DAList(extra),
               "xxx")

  out_of_order <- input[rev(names(input))]
  expect_message(validate_DAList(out_of_order),
                 "Reordering.")
})

test_that("validate_DAList checks data slot", {
  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata)

  input$data <- as.matrix(input$data)
  # Works if a matrix
  expect_s3_class(validate_DAList(input), class = "DAList")
  # Fails if not a matrix of df
  input$data <- "fail"
  expect_error(validate_DAList(input), " data frame or matrix")

  # Fails if df is not numeric
  input$data <- dfs$data
  input$data$sampleA <- as.character(input$data$sampleA)
  expect_error(validate_DAList(input), "numeric")
})

test_that("validate_DAList checks the annotation slot", {

  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata)
  # Fails if uniprot_id not unique
  input$annotation$uniprot_id[1] <- input$annotation$uniprot_id[2]

  expect_error(validate_DAList(input), "not unique")

  # Fails if uniprot_id not present
  colnames(input$annotation) <- c("x", "y")
  expect_error(validate_DAList(input), "uniprot_id")

})





