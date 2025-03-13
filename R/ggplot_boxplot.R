#' Generate a QC Boxplot
#'
#' @param data A matrix or data frame of intensities.
#' @param groups A vector indicating class membership (numeric, integer, character, factor, or NULL).
#' @param sample_labels A vector of names to display on the plot (character or NULL).
#' @param title A character string for the plot title (numeric, integer, character, or NULL).
#' @param text.sizes A numeric vector specifying text sizes for title, x-axis, y-axis, and legend.
#' @param legend.position A character string specifying the legend position (default: "right").
#' @param boxplot_width A numeric value specifying the width of the boxplots (default: dynamic based on groups).
#' @param boxplot_alpha A numeric value specifying the transparency level of the boxplots (default: dynamic based on groups).
#' @param plot_margin A grid unit object specifying the plot margins (default: dynamic based on groups).
#' @param colorblind_palette A character vector specifying colorblind-friendly colors to use (default: predefined palette).
#'
#' @return A ggplot2 object displaying the QC boxplot.
#'
#' @examples
#' \dontrun{
#'# example of boxplot
#'
#'}
#' @export

qc_boxplot <- function (data,
                        groups = NULL,
                        sample_labels = NULL,
                        title = NULL,
                        text.sizes = c(12,10,10,10),
                        legend.position="right",
                        boxplot_width = NULL,
                        boxplot_alpha = NULL,
                        plot_margin = NULL,
                        colorblind_palette = NULL){

  # Define a default colorblind-friendly palette
  base_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  if (is.null(colorblind_palette)) {
    colorblind_palette <- base_palette
  }

  # Validate input arguments
  stopifnot(is.matrix(data) || is.data.frame(data))
  if (!is.null(groups)) {
    stopifnot(length(groups) == ncol(data))
  }
  if (!is.null(sample_labels)) {
    stopifnot(length(sample_labels) == ncol(data))
  }
  stopifnot(length(text.sizes) == 4)

  # Prepare arguments
  if (is.null(groups)) {
    groups <- rep("1", ncol(data))
  }
  groups <- factor(groups)

  if (is.null(sample_labels)) {
    sample_labels <- colnames(data)
  }
  if (is.null(title)) {
    title <- " "
  }

  # Sort groups by levels
  group_order <- order(groups)

  # Reorder data accordingly
  sample_labels <- sample_labels[group_order]
  groups <- groups[group_order]
  data <- data[, group_order, drop = FALSE]

  # Create metadata frame
  plot.meta <- data.frame(
    ind = factor(colnames(data)),
    labels = factor(sample_labels),
    group = groups
  )

  # Convert data to long format
  plot.data <- merge(plot.meta, utils::stack(as.data.frame(data)), sort = FALSE)

  # Extend the colorblind-friendly palette to match the number of groups
  unique_groups <- length(unique(groups))
  colorblind_palette <- rep(colorblind_palette, length.out = unique_groups)

  # Dynamically adjust width, alpha, and plot margins if not provided
  if (is.null(boxplot_width)) {
    boxplot_width <- max(0.1, min(0.5, 1 / unique_groups))
  }
  if (is.null(boxplot_alpha)) {
    boxplot_alpha <- ifelse(unique_groups > 6, 0.4, 0.2)
  }
  if (is.null(plot_margin)) {
    plot_margin <- unit(rep(ifelse(unique_groups > 6, 0.3, 0.2), 4), "cm")
  }

########
# run testthat
#
#   Default arguments: This test checks if the function works with default arguments and returns a ggplot object. It also verifies that the title is set correctly and that the axis title text size is as expected.
#
#   Input validation: The test ensures that an error is thrown when the input data is not a matrix or data frame, or when the length of the groups argument does not match the number of columns in data.
#
#   Colorblind palette: This test checks if the colorblind_palette argument is correctly applied. The color palette should be applied to the plot as specified.
#
#   Custom margins: This test verifies if the custom plot_margin is correctly set and applied to the plot.
#
#   Boxplot width and transparency: This test checks if the boxplot_width and boxplot_alpha arguments are correctly passed into the plot's aesthetic parameters.
#
# These tests ensure that your function behaves correctly across different inputs and edge cases.
#########
  #  test_that("qc_boxplot works with default arguments", {
  #    data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  #    p <- qc_boxplot(data)
  #
  #    # Check if the output is a ggplot object
  #    expect_s3_class(p, "gg")
  #
  #    # Check if the plot title is correctly set
  #    expect_equal(p$labels$title, " ")
  #
  #    # Check if default text size is applied to the axis title
  # # #  expect_equal(p$theme$axis.title.y$size, 10)
  #  })

  ## check number arg values match no cols data
  testthat::expect_equal(object = length(groups), expected = ncol(data))
  testthat::expect_equal(object = length(sample_labels), expected = ncol(data))
  testthat::expect_equal(object = length(text.sizes), expected = 4)

  test_that("qc_boxplot throws error for non-matrix data", {
    data <- list(a = 1:10, b = 11:20)
    expect_error(qc_boxplot(data), "is.matrix(data) || is.data.frame(data)")
  })

  # test_that("qc_boxplot throws error when groups length doesn't match data columns", {
  #   data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  #   groups <- rep(1, 5)
  #   expect_error(qc_boxplot(data, groups = groups), "length(groups) == ncol(data)")
  # })
  #
  # test_that("qc_boxplot handles colorblind_palette correctly", {
  #   data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  #   groups <- rep(1:2, each = 5)
  #   color_palette <- c("#FF0000", "#00FF00")
  #   p <- qc_boxplot(data, groups = groups, colorblind_palette = color_palette)
  #
  #   # Check if the colors used in the plot match the colorblind_palette
  #   colors_used <- p$scales$scales[[1]]$palette(2)  # Should return 2 colors from the palette
  #   expect_equal(colors_used, color_palette)
  # })
  #
  # test_that("qc_boxplot handles custom margins", {
  #   data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  #   p <- qc_boxplot(data, plot_margin = unit(c(1, 1, 1, 1), "cm"))
  #
  #   # Check if plot margins have been set
  #   expect_equal(p$theme$plot.margin, unit(c(1, 1, 1, 1), "cm"))
  # })
  #
  # test_that("qc_boxplot handles boxplot width and alpha correctly", {
  #   data <- matrix(rnorm(100), nrow = 10, ncol = 10)
  #   p <- qc_boxplot(data, boxplot_width = 0.3, boxplot_alpha = 0.5)
  #
  #   # Check if width and alpha are passed correctly
  #   expect_equal(p$layers[[1]]$mapping$width, 0.3)
  #   expect_equal(p$layers[[1]]$aes_params$alpha, 0.5)
  # })

########
# Generate the plot
  p1 <- ggplot(plot.data, aes(x = labels, y = values, fill = group)) +
    geom_boxplot(width = boxplot_width, color = "black", alpha = boxplot_alpha) +
    scale_fill_manual(values = colorblind_palette, name = NULL) +
    labs(y = "Density", x = "", fill = "") +
    theme_gray() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = text.sizes[3]),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = text.sizes[2]),
      axis.text.y = element_text(size = text.sizes[3]),
      plot.title = element_text(size = text.sizes[1]),
      plot.margin = plot_margin,
      legend.position = legend.position,
      legend.text = element_text(size = text.sizes[4]),
      legend.key.size = unit(0.5, "cm")
    ) +
    ggtitle(title)

  return(p1)
}
