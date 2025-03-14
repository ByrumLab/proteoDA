#' Generate a Barplot of Total Intensities by Sample Group
#'
#' This function creates a barplot where samples are organized and colored by group.
#' It also identifies potential outliers based on a specified percentile of total intensity values.
#'
#' @param DAList A list containing `metadata` and `data`. The `metadata` should include sample labels and grouping information.
#' @param label_column A string specifying the column name in `metadata` containing sample labels.
#' @param grouping_column A string specifying the column name in `metadata` containing grouping information.
#' @param percentile A numeric value (0-1) specifying the percentile threshold for identifying low-intensity samples.
#' @param colors A vector of colors (one per group). If NULL, a color-blind friendly palette (viridis) is used.
#' @param legend.position A string specifying legend position. Options: "none", "right", "top", "left", "bottom".
#'
#' @return A list containing:
#'   - `p`: A ggplot object of the barplot.
#'   - `bar_data`: A data frame with total intensity values and grouping information.
#'   - `colors`: The color palette used.
#'   - `percentile`: The input percentile value.
#'   - `perc_int_thresh`: The computed intensity threshold.
#'
#' @import ggplot2
#' @import viridisLite
#' @importFrom stats quantile
#'
#' @examples
#' \dontrun {
#'
#'barplot1 <- qc_totInt(DAList = results,
#'                  label_column = "sample",
#'                  grouping_column = "group",
#'                  percentile = 0,
#'                  colors = NULL,
#'                  legend.position = "right")
#'
#' }
#'
#'
#' @export
qc_totInt <- function(DAList, label_column, grouping_column, percentile, colors = NULL, legend.position = "right") {

  legend.position <- rlang::arg_match(legend.position, c("none", "right", "top", "left", "bottom"))

 # label_column <- check_label_column(metadata = DAList$metadata, label_column = label_column)
 # grouping_column <- check_grouping_column(metadata = DAList$metadata, grouping_column = grouping_column)

  # Validate percentile parameter
  if (!is.numeric(percentile) || percentile < 0 || percentile > 1) {
    cli::cli_abort("percentile should be a numeric value between 0 and 1")
  }

  ## create data.frame with sample labels, totalIntensity, and group column
  bar_data <- data.frame(label    = DAList$metadata[, label_column],
                         tot.int = as.numeric(colSums(DAList$data, na.rm = TRUE)),
                         tot.num = as.numeric(colSums(!is.na(DAList$data),
                                                      na.rm = FALSE)),
                         group    = DAList$metadata[, grouping_column],
                         check.names = FALSE,
                         fix.empty.names = FALSE,
                         stringsAsFactors = FALSE,
                         check.rows = FALSE)

  # Convert group column to factor and sort
  bar_data$group <- factor(bar_data$group, levels = unique(bar_data$group))
  bar_data <- bar_data[order(bar_data$group), ]

  # Compute intensity threshold for the given percentile
  ### percentile A double [0-1] used to identify the specified percentile of the data.
  ## takes the total intensity values for the samples, and identifies the intensity
  ## threshold for the say 10th percentile of the data. samples with values below
  ## the percentile are flagged as potential outliers. so it identifies samples with
  ## the lowest total intensities in the data. these samples may not have sequenced
  ## well and will likely have a lot of missing values. the horizontal line
  ## on the plots corresponds to the intensity threshold for the specified percentile
  perc_int_thresh <- ifelse(percentile == 0, 0, quantile(bar_data$tot.int, probs = percentile))

  # Flag potential outliers
  bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh

  # Assign colors
  num_groups <- length(unique(bar_data$group))
  if (is.null(colors)) {
    colors <- viridisLite::viridis(num_groups)
  }
  if (length(colors) != num_groups) {
    cli::cli_abort("Number of colors does not match the number of unique groups")
  }
  names(colors) <- levels(bar_data$group)

  # Generate the plot
  p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = colors) +
    geom_hline(yintercept = perc_int_thresh, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(y = "Total Intensity", fill = grouping_column) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = legend.position)

  # Return results
  return(list(
    p = p,
    bar_data = bar_data,
    colors = colors,
    percentile = percentile,
    perc_int_thresh = perc_int_thresh
  ))
}
