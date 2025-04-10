test_that("compute_moving_sd_and_zscores works with valid input", {
  # Simulate a simple raw object
  set.seed(123)
  n <- 200
  raw <- list(
    data = data.frame(
      matrix(rnorm(n * 5), nrow = n, ncol = 5),
      Log2Fold = rnorm(n)
    ),
    results = list(
      treat_vs_control = list(logFC = rnorm(n)),
      drug_vs_control = list(logFC = rnorm(n))
    ),
    tags = list()
  )
  
  binsize <- 50
  raw <- compute_movingSD_zscores(raw, binsize)
  
  # Check output structure
  expect_true("movingSDs" %in% names(raw$tags))
  expect_true("logFC_z_scores" %in% names(raw$tags))
  
  # Check each comparison has expected output
  for (comp in names(raw$results)) {
    expect_equal(length(raw$tags$movingSDs[[comp]]), n)
    expect_equal(length(raw$tags$logFC_z_scores[[comp]]), n)
    expect_false(any(is.na(raw$tags$movingSDs[[comp]])))
    expect_false(any(is.na(raw$tags$logFC_z_scores[[comp]])))
  }
})
