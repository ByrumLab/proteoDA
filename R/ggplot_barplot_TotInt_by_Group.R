#' Generate Total Intensity Barplots by Group
#'
#' This function creates barplots showing total intensity values for samples, grouped by a specified column.
#' It includes an optional percentile threshold to flag potential outliers and uses a color-blind friendly palette.
#'
#' @param DAList A list containing `data` (matrix) and `metadata` (data frame).
#' @param label_column A character string specifying the column name in metadata for sample labels.
#' @param grouping_column A character string specifying the column name in metadata for grouping.
#' @param percentile A numeric value (0-1) for the percentile threshold to flag low-intensity samples.
#' @param colors A vector of colors for groups. Defaults to a color-blind friendly palette (viridis).
#' @param nrow Number of rows for facet wrapping.
#' @param ncol Number of columns for facet wrapping.
#' @param legend.position Position of the legend ('none', 'right', 'top', 'left', 'bottom').
#'
#' @return A list containing ggplot objects for total intensity, total number, and total missing values barplots.
#' @examples
#' # example code
#'
#' \dontrun {
#'
#'# example of barplot colored by group
#'barplot1 <- qc_totInt_by_group(DAList = results,
#'                             label_column = "sample",
#'                             grouping_column = "group",
#'                             percentile = 0,
#'                             colors = NULL,   # or c("#E69F00", "#56B4E9", "#009E73")
#'                             legend.position = "right")
#'
#' }
#'
#' @export
qc_totInt_by_group <- function(DAList,
                               label_column,
                               grouping_column,
                               percentile,
                               colors = NULL,
                               nrow = NULL,
                               ncol = NULL,
                               legend.position = "none") {

  # Validate legend position argument
  legend.position <- rlang::arg_match(legend.position,
                                      c("none", "right", "top", "left", "bottom"))

  # Validate columns
  #label_column <- check_label_column(DAList$metadata, label_column)
  #grouping_column <- check_grouping_column(DAList$metadata, grouping_column)

  # Check percentile input
  if (!is.numeric(percentile) || any(percentile < 0 | percentile > 1)) {
    cli::cli_abort("Percentile should be a numeric value between 0 and 1.")
  }

  # Prepare data for plotting
  # bar_data <- data.frame(
  #   label = DAList$metadata[[label_column]],
  #   tot.int = rowSums(DAList$data, na.rm = TRUE),
  #   tot.num = rowSums(!is.na(DAList$data)),
  #   tot.na = rowSums(is.na(DAList$data)),
  #   group = DAList$metadata[[grouping_column]],
  #   stringsAsFactors = FALSE
  # )

  bar_data <- data.frame(
    label = DAList$metadata[[label_column]],
    tot.int = colSums(DAList$data, na.rm = TRUE),
    tot.num = colSums(!is.na(DAList$data)),
    tot.na = colSums(is.na(DAList$data)),
    group = DAList$metadata[[grouping_column]],
    stringsAsFactors = FALSE
  )

  # Convert group column to factor and sort
  bar_data$group <- factor(bar_data$group, levels = unique(bar_data$group))
  num_groups <- length(unique(bar_data$group))

  # Calculate percentile threshold
  perc_int_thresh <- ifelse(percentile == 0, 0, quantile(bar_data$tot.int, probs = percentile))
  bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh

  # Assign colors if not provided
  if (is.null(colors)) {
    colors <- viridis::viridis(num_groups, option = "D", direction = -1)
  }

  # Check color length matches number of groups
  if (length(colors) != num_groups) {
    cli::cli_abort("Number of colors does not match the number of unique groups.")
  }

  names(colors) <- levels(bar_data$group)

  # Generate ggplot objects
  p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    geom_hline(yintercept = perc_int_thresh, linetype = "dashed", color = "black") +
    labs(y = "Total Intensity", fill = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = legend.position)

  list(plot = p, data = bar_data, colors = colors, threshold = perc_int_thresh)
}
