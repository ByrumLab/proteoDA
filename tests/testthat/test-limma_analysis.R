
# tests for limma_analysys ------------------------------------------------

# test fit_limma_model ----------------------------------------------------

test_that("fit_limma_model gives error when no design is present", {
  input_1 <- readRDS(test_path("fixtures", "fit_limma_model_input.rds"))
  input_1$design <- NULL
  input_2 <- readRDS(test_path("fixtures", "fit_limma_model_input_nodesign.rds"))

  expect_error(fit_limma_model(input_1), "statistical design")
  expect_error(fit_limma_model(input_2), "statistical design")

})


test_that("fit_limma_model gives message when input data are not normalized", {
  input <- readRDS(test_path("fixtures", "fit_limma_model_input_nonorm.rds"))

  expect_message(fit_limma_model(input), "not normalized")
  input$tags$normalized <- F
  expect_message(fit_limma_model(input), "not normalized")

})


test_that("fit_limma_model gives message on intra-block correlation with random effects models", {
  input <- suppressMessages(readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
    add_design(~ 0 + treatment + (1 | group)))

  expect_message(try(fit_limma_model(input), silent = T), # need a try here to ignore the error that happends because the cor is negative
                 "correlation =")

})

test_that("fit_limma_model gives error when intra-block correlation is negative", {

  input <- suppressMessages(readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
    add_design(~ 0 + treatment + (1 | group)))

  suppressMessages(expect_error(fit_limma_model(input),
                 "negative"))
})

test_that("fit_limma_model gives message when overwriting existing statistical results", {
  # should work for both eBayes fit and results

  input <- readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
    fit_limma_model()

  expect_message(fit_limma_model(input), "Overwriting")
  input2 <- suppressMessages(input |>
    extract_DA_results())

  expect_message(fit_limma_model(input2), "Overwriting")

  input3 <- input2
  input3["eBayes_fit"] <- list(NULL)

  expect_message(fit_limma_model(input3), "Overwriting")

})

test_that("fit_limma_model gives consistent output", {
  # Test a few different model types
  # Just using snapshot testing to check for changes

  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  a <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + group)

  b <- input |>
    normalize_data("log2") |>
    add_design(~ group)

  c <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + treatment) |>
    add_contrasts(contrasts_vector = c("Treatment_vs_Control= treatment - control"))

  d <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + sex + (1 | treatment))

  expect_snapshot(fit_limma_model(a))
  expect_snapshot(fit_limma_model(b))
  expect_snapshot(fit_limma_model(c))
  expect_snapshot(fit_limma_model(d))

})


# test extract_DA_results -------------------------------------------------

test_that("extract_DA_results checks arguments", {
  input <- readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
    fit_limma_model()

  expect_error(extract_DA_results(input, pval_thresh = -1), "pval_thresh")
  expect_error(extract_DA_results(input, pval_thresh = 1.1), "pval_thresh")
  expect_error(extract_DA_results(input, lfc_thresh = -1), "lfc_thresh")
  expect_error(extract_DA_results(input, adj_method = "xxx"), "xxx")

})

test_that("extract_DA_results errors when no eBayes fit is present", {
  input <- readRDS(test_path("fixtures", "fit_limma_model_input.rds"))

  suppressMessages(expect_error(extract_DA_results(input), "does not have a model fit"))

})

test_that("extract_DA_results informs user when overwriting old results", {
  input <- suppressMessages(readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
    fit_limma_model() |>
    extract_DA_results())

  suppressMessages(expect_message(extract_DA_results(input), "Overwriting."))

})


test_that("extract_DA_results handles intercept terms as expected", {
  input_int <- suppressMessages(
    readRDS(test_path("fixtures", "fit_limma_model_input.rds")) |>
      add_design(~ group) |>
      fit_limma_model() |>
      extract_DA_results()
    )

  no_int <- suppressMessages(extract_DA_results(input_int))
  with_int <- suppressMessages(extract_DA_results(input_int, extract_intercept = T))

  expect_false("Intercept" %in% names(no_int$results))
  expect_true("Intercept" %in% names(with_int$results))

})

test_that("extract_DA_results gives consistent results", {

  # Test a few different model types
  # Just using snapshot testing to check for changes

  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  a <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + group) |>
    fit_limma_model()

  b <- input |>
    normalize_data("log2") |>
    add_design(~ group) |>
    fit_limma_model()

  c <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + treatment) |>
    add_contrasts(contrasts_vector = c("Treatment_vs_Control= treatment - control")) |>
    fit_limma_model()

  d <- input |>
    normalize_data("log2") |>
    add_design(~ 0 + sex + (1 | treatment)) |>
    fit_limma_model()

  suppressMessages(expect_snapshot(extract_DA_results(a)))
  suppressMessages(expect_snapshot(extract_DA_results(b)))
  suppressMessages(expect_snapshot(extract_DA_results(c)))
  suppressMessages(expect_snapshot(extract_DA_results(d)))

})



# check_DA_perc -----------------------------------------------------------

test_that("check_DA_perc gives warning when above threshold", {
  fake_DA_outcomes <- data.frame(termA2 = c(1,-1,0,0,0,0,0,0,0,0),
                                 termB5 = c(1,-1,1,-1,1,0,0,0,0,0),
                                 termC8 = c(1,1,1,1,1,1,1,-1, 0,0)) |>
    as.matrix()

  # Default threshold
  # Not ideal: running into weird problems testing the term names specifically,
  # some issue with matching quotes and spaces in the regex.
  # Enough to know that this should match multiple terms,
  # and the threshold below should match only 1
  expect_message(check_DA_perc(fake_DA_outcomes,
                               pval_thresh = 0.055,
                               lfc_thresh = 1,
                               adj_method = "none"),
                 regexp = "Problematic terms:")

  # Highest threshold
  expect_message(check_DA_perc(fake_DA_outcomes,
                               DA_warn_threshold = 0.7,
                               pval_thresh = 0.055,
                               lfc_thresh = 1,
                               adj_method = "none"),
                 "termC8")

})

test_that("check_DA_perc is silent when below threshold", {
  fake_DA_outcomes <- data.frame(termA2 = c(1,-1,0,0,0,0,0,0,0,0),
                                 termB5 = c(1,-1,1,-1,1,0,0,0,0,0),
                                 termC8 = c(1,1,1,1,1,1,1,-1, 0,0)) |>
    as.matrix()

  expect_silent(check_DA_perc(fake_DA_outcomes,
                               DA_warn_threshold = 0.9,
                               pval_thresh = 0.055,
                               lfc_thresh = 1,
                               adj_method = "none"))

})

test_that("check_DA_perc returns correct percentage", {
  fake_DA_outcomes <- data.frame(termA2 = c(1,-1,0,0,0,0,0,0,0,0),
                                 termB5 = c(1,-1,1,-1,1,0,0,0,0,0),
                                 termC8 = c(1,1,1,1,1,1,1,-1, 0,0)) |>
    as.matrix()
  expect_equal(check_DA_perc(fake_DA_outcomes,
                              DA_warn_threshold = 0.9,
                              pval_thresh = 0.055,
                              lfc_thresh = 1,
                              adj_method = "none"),
               c(termA2 = 0.2, termB5 = 0.5, termC8 = 0.8))

})


test_that("check_DA_perc does not check intercept columns", {
  fake_DA_outcomes_with_int <- data.frame(Intercept = c(1,-1,1,1,1,1,1,1,1,1),
                                          intercept = c(1,-1,1,1,1,1,1,1,1,1),
                                          termB5 = c(1,-1,1,-1,1,0,0,0,0,0),
                                          termC8 = c(1,1,1,1,1,1,1,-1, 0,0)) |>
    as.matrix()

  expect_silent(check_DA_perc(fake_DA_outcomes_with_int,
                             DA_warn_threshold = 0.9,
                             pval_thresh = 0.055,
                             lfc_thresh = 1,
                             adj_method = "none"))

  expect_equal(check_DA_perc(fake_DA_outcomes_with_int,
                             DA_warn_threshold = 0.9,
                             pval_thresh = 0.055,
                             lfc_thresh = 1,
                             adj_method = "none"),
               c(termB5 = 0.5, termC8 = 0.8))

})
