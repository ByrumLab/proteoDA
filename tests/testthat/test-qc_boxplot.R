test_that("qc_boxplot works with default arguments", {
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  p <- qc_boxplot(data)
  
  expect_s3_class(p, "gg")
  expect_equal(p$labels$title, " ")
})


test_that("qc_boxplot validates input", {
  # If a list is passed that looks like a DAList but lacks data slots,
  # the function should error with an informative message.
  expect_error(
    qc_boxplot(list(a = 1:10)),
    "no usable data found in data_per_contrast or data"
  )
  
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  expect_error(qc_boxplot(data, groups = rep(1, 5)))
})

test_that("qc_boxplot applies colorblind palette", {
  data   <- matrix(rnorm(100), nrow = 10, ncol = 10)
  groups <- rep(1:2, each = 5)
  pal    <- c("#FF0000", "#00FF00")
  
  p <- qc_boxplot(data, groups = groups, colorblind_palette = pal)
  
  scale <- p$scales$scales[[1]]
  # Ask the scale what colors it will use for the two group levels
  colors_used <- scale$palette(2)
  expect_equal(colors_used, pal)
})

test_that("qc_boxplot handles custom margins", {
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  custom_margins <- unit(c(1, 1, 1, 1), "cm")
  
  p <- qc_boxplot(data, plot_margin = custom_margins)
  
  expect_equal(p$theme$plot.margin, custom_margins)
})

test_that("qc_boxplot handles custom width and alpha", {
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  p <- qc_boxplot(data, boxplot_width = 0.3, boxplot_alpha = 0.5)
  
  layer <- p$layers[[1]]
  
  expect_equal(layer$aes_params$alpha, 0.5)
  
  built <- ggplot2::ggplot_build(p)
  expect_true("alpha" %in% names(built$data[[1]]))
  expect_true(all(built$data[[1]]$alpha == 0.5))
})

### new test for before norm plot
test_that("qc_boxplot_beforeNorm works with matrix input", {
  data <- matrix(rnorm(100, mean = 20, sd = 5), nrow = 10, ncol = 10)
  
  p <- qc_boxplot_beforeNorm(data)
  
  expect_s3_class(p, "gg")
})

test_that("qc_boxplot_beforeNorm applies log2 transform when data not log scaled", {
  data <- matrix(runif(100, min = 1e3, max = 1e6), nrow = 10, ncol = 10)
  
  p <- qc_boxplot_beforeNorm(data)
  
  # Build the ggplot object
  built <- ggplot2::ggplot_build(p)
  
  # Defensive: ensure we have at least one layer
  expect_true(length(built$data) >= 1, info = "ggplot_build produced no data layers")
  
  layer_df <- built$data[[1]]
  
  # Collect numeric columns from the layer's data.frame
  num_cols <- vapply(layer_df, is.numeric, logical(1))
  numeric_vals <- unlist(layer_df[, num_cols, drop = FALSE], use.names = FALSE)
  
  # Ensure we have some finite numeric values to test
  expect_true(length(numeric_vals) > 0 && any(is.finite(numeric_vals)),
              info = "No finite numeric values produced by the first plot layer")
  
  # Use the finite numeric values to check the scale; use a generous cutoff
  max_val <- max(numeric_vals, na.rm = TRUE)
  expect_true(is.finite(max_val), info = "Max of numeric layer values is not finite")
  expect_true(max_val < 50)
})

test_that("qc_boxplot_beforeNorm accepts DAList-like input", {
  dat <- matrix(rnorm(100, mean = 10), nrow = 10, ncol = 10)
  da  <- list(data = dat)
  
  p <- qc_boxplot_beforeNorm(da)
  
  expect_s3_class(p, "gg")
})

test_that("qc_boxplot_beforeNorm errors on invalid input", {
  expect_error(
    qc_boxplot_beforeNorm(list(a = 1:10)),
    "no usable 'data' matrix"
  )
})

