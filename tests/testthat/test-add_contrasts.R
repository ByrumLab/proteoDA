
# Tests for add_contrasts -------------------------------------------------
# Main function is tested here, with some subfunctions tested below

test_that("add_contrasts gives useful error when DAList doesn't have a statistical design", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  good_input <- c("Treatment1_vs_Control= Treatment1 - Control",
                  "Treatment2_vs_Control= Treatment2-Control",
                  "Treatment2 vs Treatment1=Treatment2-Treatment1")

  expect_error(add_contrasts(input, contrasts_vector = good_input), "not have a statistical design")
})


test_that("add_contrasts gives useful error when both contrasts_vector and contrasts_file are provided", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))
  good_input <- c("Treatment1_vs_Control= Treatment1 - Control",
                  "Treatment2_vs_Control= Treatment2-Control",
                  "Treatment2 vs Treatment1=Treatment2-Treatment1")

  expect_error(
    add_contrasts(input,
                  contrasts_vector = good_input,
                  contrasts_file = "file.csv"), "Must supply either")
})

test_that("add_contrasts gives useful error when neither contrasts_vector and contrasts_file are provided", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  expect_error(
    add_contrasts(input), "Must supply either")

})

test_that("add_contrasts gives useful error when contrasts_file is the wrong type", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  expect_error(
    add_contrasts(input,
                  contrasts_file = "file.xxx"), "Input file must end")
})

test_that("add_contrasts gives useful error when contrasts_file is not found", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  expect_error(
    add_contrasts(input,
                  contrasts_file = "file.csv"), "does not exist")

})

test_that("add_contrasts gives useful error when contrasts_file has more than one column", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  expect_error(
    add_contrasts(input,
                  contrasts_file = test_path("fixtures", "bad_contrast_mult_col.csv")),
    "more than one column")
})

test_that("add_contrasts gives useful error when groups in contrast are not present in the design", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  expect_error(
    add_contrasts(input,
                  contrasts_file = test_path("fixtures", "bad_contrast_wrong_group.csv")),
    "XXX")
})

test_that("add_contrasts notifies user when overwriting an existing contrast", {

  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))
  existing_contrast <- input |>
    add_contrasts(contrasts_vector = c("Treatment_vs_Control= treatment - control"))


  expect_message(
    add_contrasts(existing_contrast,
                  contrasts_vector = c("Treatment_vs_Control= treatment - control")),
    "DAList already contains contrasts. Overwriting."
  )
})


# warning when deleting existing model fit
test_that("add_contrasts gives warning when deleting existing model fit", {
  input <- suppressMessages(
    readRDS(test_path("fixtures", "add_design_input.rds")) |>
      add_design(~ 0 + treatment) |>
      fit_limma_model()
  )
  suppressMessages(
    expect_message(
      add_contrasts(input,
                    contrasts_vector = c("Treatment_vs_Control= treatment - control")
                    ),
      "model fit"
    )
  )

  suppressMessages(
    expect_null(
      add_contrasts(input,
                    contrasts_vector = c("Treatment_vs_Control= treatment - control")
      )$eBayes_fit
    )
  )
})

# warning when deleting existing results
test_that("add_contrasts gives warning when deleting existing results", {
  input <- suppressMessages(
    readRDS(test_path("fixtures", "add_design_input.rds")) |>
      add_design(~ 0 + treatment) |>
      fit_limma_model() |>
      extract_DA_results()
  )
  suppressMessages(
    expect_message(
      add_contrasts(input,
                    contrasts_vector = c("Treatment_vs_Control= treatment - control")
      ),
      "results"
    )
  )

  suppressMessages(
    expect_null(
      add_contrasts(input,
                    contrasts_vector = c("Treatment_vs_Control= treatment - control")
      )$results
    )
  )
})



test_that("add_contrasts outputs proper contrast when contrasts_vector is supplied", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  contrast <- c("Treatment_vs_Control= treatment - control")

  output <- add_contrasts(input,
                          contrasts_vector = contrast)

  expected_matrix <-
    matrix(data = c(-1,1),
           nrow = 2,
           ncol = 1,
           dimnames = list(
             Levels = c("control", "treatment"),
             Contrasts = c("Treatment_vs_Control")
           )
    )

  expect_equal(output$design$contrast_vector, contrast)
  expect_equal(output$design$contrast_matrix, expected_matrix)
})

test_that("add_contrasts outputs proper contrast when contrasts_file is supplied", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  contrast <- c("Treatment_vs_Control= treatment - control")

  output <- add_contrasts(input,
                          contrasts_file = test_path("fixtures", "good_contrast.csv"))

  expected_matrix <-
    matrix(data = c(-1,1),
           nrow = 2,
           ncol = 1,
           dimnames = list(
             Levels = c("control", "treatment"),
             Contrasts = c("Treatment_vs_Control")
           )
    )

  expect_equal(output$design$contrast_vector, contrast)
  expect_equal(output$design$contrast_matrix, expected_matrix)

})

test_that("add_contrasts outputs proper contrast when contrasts_file is supplied as txt", {
  input <- readRDS(test_path("fixtures", "add_contrasts_input.rds"))

  contrast <- c("Treatment_vs_Control= treatment - control")

  output <- add_contrasts(input,
                          contrasts_file = test_path("fixtures", "good_contrast.txt"))

  expected_matrix <-
    matrix(data = c(-1,1),
           nrow = 2,
           ncol = 1,
           dimnames = list(
             Levels = c("control", "treatment"),
             Contrasts = c("Treatment_vs_Control")
           )
    )

  expect_equal(output$design$contrast_vector, contrast)
  expect_equal(output$design$contrast_matrix, expected_matrix)

})




# Tests for validate_contrasts --------------------------------------------
test_that("validate contrasts checks for the right number of equals signs", {
  good_input <- c("Treatment1_vs_Control= Treatment1 - Control",
                  "Treatment2_vs_Control= Treatment2-Control",
                  "Treatment2 vs Treatment1=Treatment2-Treatment1")
  no_equals <- c("Treatment1_vs_Control Treatment1 - Control",
                 "Treatment2_vs_Control = Treatment2-Control",
                 "Treatment2 vs Treatment1=Treatment2-Treatment1")
  multiple_equals <- c("Treatment1_vs_Control == Treatment1 - Control",
                       "Treatment2_vs_Control == Treatment2-Control",
                       "Treatment2 vs Treatment1=Treatment2-Treatment1")

  expect_equal(validate_contrasts(good_input), good_input)
  expect_error(validate_contrasts(no_equals), "1")
  expect_error(validate_contrasts(multiple_equals), "1 and 2")

})

# Tests for extract_contrast_groups ---------------------------------------

test_that("extract contrast groups extracts group names as expected", {

  # Right now, allow group names to be numbers, but only allowed to add and subtract them
  input <- c("Treatment1_vs_Control= Treatment1 - Control",
             "Treatment2_vs_Control= Treatment2-Control",
             "Treatment2 vs Treatment1=Treatment2-Treatment1")
  output <- c("Treatment1", "Control", "Treatment2")

  input2 <- c("something = 999 - 100",
              "something = XXX - YYY",
              "something else = (XXX-YYY)/2",
              "another one = (a+b+c+d+e+f+g+h+i+j)/10")

  output2 <- c("999", "100", "XXX", "YYY", letters[1:10])

  expect_equal(unique(unlist(extract_contrast_groups(input))), output)
  expect_equal(unique(unlist(extract_contrast_groups(input2))), output2)

})


