test_that("qc_boxplot works with default arguments", {
  data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  p <- qc_boxplot(data)
  
  expect_s3_class(p, "gg")
  expect_equal(p$labels$title, " ")
})

test_that("qc_boxplot validates input", {
  expect_error(qc_boxplot(list(a = 1:10)), "is.matrix")
  
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
  
  # width should be in geom_params
  layer <- p$layers[[1]]
  expect_equal(layer$geom_params$width, 0.3)
  
  # alpha is only visible after building the plot
  built <- ggplot2::ggplot_build(p)
  alpha_vals <- built$data[[1]]$alpha
  
  expect_true(is.numeric(alpha_vals))
  expect_true(all(alpha_vals == 0.5))
})

