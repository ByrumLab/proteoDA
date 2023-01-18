
# tests for limma_analysys ------------------------------------------------
# Main functions tested here, subfunctiosn below


# test fit_limma_model ----------------------------------------------------




# test extract_DA_results -------------------------------------------------




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
