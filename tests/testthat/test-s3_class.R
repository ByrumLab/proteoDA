

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

test_that("DAList() gives useful error when data and annotation have different number of rows", {

  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))


  expect_error(DAList(data = dfs$data[1:2,],
                      annotation = dfs$annotation,
                      metadata = dfs$metadata),
               "same number of rows")
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

test_that("validate_DAList checks that annotation and data match", {
  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))


  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata)
  input$data <- input$data[1:2,]

  # Same number of rows
  expect_error(validate_DAList(input),
               "same number of rows")

  # same row names
  input$data <- dfs$data
  rownames(input$data) <- c("xxx", "M9NH73", "H9AZR4")
  expect_error(validate_DAList(input),
               "must match")

  # Rownames equal the uniprot_id col
  rownames(input$annotation)[1] <- "xxx"
  expect_error(validate_DAList(input),
               "uniprot_id")
})

test_that("validate_DAList checks the metadata", {

  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata)
  # rows of metadata match cols of data
  input$metadata <- input$metadata[1:2,]
  expect_error(validate_DAList(input),
               "number of samples in the data")

  # rownames metadata match colnames data
  input$metadata <- dfs$metadata
  rownames(input$metadata)[1] <- "xxx"
  expect_error(validate_DAList(input),
               "column names of the data")
})

test_that("validate_DAList checks the design", {
  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group)
  input$design <- list(x = 1)

  # Checks length of allowable design
  # Too short
  expect_error(validate_DAList(input),
               "of at least 2")

  # Checks for not-allowed names
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group)
  input$design <- c(input$design, extra = "extra")
  expect_error(validate_DAList(input),
               "extra")
  # checks for 2 required elements
  input$design$extra <- NULL
  input$design$random_factor <- "temp"
  input$design$design_formula <- NULL
  expect_error(validate_DAList(input),
               "must contain at least a")

  # Checks that design_matrix has proper attributes
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group)
  attributes(input$design$design_matrix) <- NULL
  expect_error(validate_DAList(input),
               "attribute")

  # checks that design matrix has same number of rows as metadata
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group)
  # Turns out delete ing rows is a little tricky,
  # have to manage the metadata too
  tmp <- attributes(input$design$design_matrix)
  tmp$dim[1] <- 2
  tmp$dimnames[[1]] <- c("sampleA", "sampleB")
  input$design$design_matrix <- input$design$design_matrix[1:2,]
  attributes(input$design$design_matrix) <- tmp
  expect_error(validate_DAList(input),
               "samples in the metadata")

  # Rownames of design matrix equal colnames of data
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group)
  rownames(input$design$design_matrix)[1] <- "xxx"
  expect_error(validate_DAList(input),
               "row names")

  # If there's a random effect, must be present in metadata
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group + (1 | sample_id))
  input$design$random_factor <- "xxx"
  expect_error(validate_DAList(input),
               "random factor term")

  # If there are contrasts, must have both contrast slots
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ 0 + group) |>
    add_contrasts(contrasts_vector = c("test = treatment - control"))
  input$design$contrast_matrix <- NULL
  expect_error(validate_DAList(input),
               "must have both")

  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ 0 + group) |>
    add_contrasts(contrasts_vector = c("test = treatment - control"))
  input$design$contrast_vector <- NULL
  expect_error(validate_DAList(input),
               "must have both")

  # Contrast matrix rows equal design matrix cols
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ 0 + group) |>
    add_contrasts(contrasts_vector = c("test = treatment - control"))
  rownames(input$design$contrast_matrix)[1] <- "xxx"
  expect_error(validate_DAList(input),
               "do not match")
})

test_that("validate_DAList checks the eBayes_fit", {

  # check class of eBayes fit
  dfs <- readRDS(test_path("fixtures", "s3-class_input.rds"))
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group) |>
    normalize_data("log2") |>
    fit_limma_model()
  input$eBayes_fit <- unclass(input$eBayes_fit)
  expect_error(validate_DAList(input),
               "class MArrayLM")

  # Random factor in design not in model
  input <- DAList(data = dfs$data,
                  annotation = dfs$annotation,
                  metadata = dfs$metadata) |>
    add_design(~ group) |>
    normalize_data("log2") |>
    fit_limma_model()
  input$design$random_factor <- "sample_id"
  expect_error(validate_DAList(input),
               "Design includes")

  # Random factor in model results no design

  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  d <- suppressMessages(input |>
                          normalize_data("log2") |>
                          add_design(~ 0 + sex + (1 | treatment)) |>
                          fit_limma_model())
  d$design$random_factor <- NULL
  expect_error(validate_DAList(d),
               "Model fit includes")
})

test_that("validate_DAList checks the results", {

  # Mismatch between rownames of data and results
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  treatment_noint_int <- suppressMessages(
    input |>
      normalize_data("log2") |>
      add_design(~ 0 + treatment) |>
      add_contrasts(contrasts_vector = c("Treatment_vs_Control= treatment - control")) |>
      fit_limma_model() |>
      extract_DA_results(extract_intercept = T)
    )

  group_int_noint <- suppressMessages(
    input |>
      normalize_data("log2") |>
      add_design(~ 0 + group) |>
      fit_limma_model() |>
      extract_DA_results(extract_intercept = F)
  )

  sex_noint_int <- suppressMessages(
    input |>
      normalize_data("log2") |>
      add_design(~ sex) |>
      fit_limma_model() |>
      extract_DA_results(extract_intercept = T)
  )

  sex_noint_noint <- suppressMessages(
    input |>
      normalize_data("log2") |>
      add_design(~ sex) |>
      fit_limma_model() |>
      extract_DA_results(extract_intercept = F)
  )

  # incorrect rownames
  bad_names <- sex_noint_int
  rownames(bad_names$results$M)[1] <- "xxx"
  expect_error(validate_DAList(bad_names),
               "between statistical results and data")

  # mismatch_contrast_noint <-

})



