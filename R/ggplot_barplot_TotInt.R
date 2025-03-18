#' Generate a Barplot of Total Intensities by Sample Group
#'
#' This function creates a barplot where samples are organized and colored by group.
#' It also identifies potential outliers based on a specified percentile of total intensity values.
#'
#' @param DAList A list containing `metadata` and `data`. The `metadata` should include sample labels and grouping information.
#' @param label_column A string specifying the column name in `metadata` containing sample labels.
#' @param grouping_column A string specifying the column name in `metadata` containing grouping information.
#' @param percentile A numeric value (0-1) specifying the percentile threshold for identifying low-intensity samples.
#' @param colors A vector of colors (one per group). If NULL, a predefined color palette is used.
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
#' @import Polychrome
#' @importFrom stats quantile
#' @export
qc_totInt <- function(DAList, label_column, grouping_column, percentile, colors = NULL, legend.position = "right") {
  
  legend.position <- rlang::arg_match(legend.position, c("none", "right", "top", "left", "bottom"))
  
  if (!is.numeric(percentile) || percentile < 0 || percentile > 1) {
    cli::cli_abort("percentile should be a numeric value between 0 and 1")
  }
  
  bar_data <- data.frame(label    = DAList$metadata[, label_column],
                         tot.int = as.numeric(colSums(DAList$data, na.rm = TRUE)),
                         tot.num = as.numeric(colSums(!is.na(DAList$data),
                                                      na.rm = FALSE)),
                         group    = DAList$metadata[, grouping_column],
                         check.names = FALSE,
                         fix.empty.names = FALSE,
                         stringsAsFactors = FALSE,
                         check.rows = FALSE)
  
  bar_data$group <- factor(bar_data$group, levels = unique(bar_data$group))
  bar_data <- bar_data[order(bar_data$group), ]
  
  perc_int_thresh <- ifelse(percentile == 0, 0, quantile(bar_data$tot.int, probs = percentile))
  bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh
  
  num_groups <- length(unique(bar_data$group))
  
  cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                  "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", 
                  "#009E73", "#F0E442")
  
  if (is.null(colors)) {
    if (num_groups <= length(cb_palette)) {
      colors <- cb_palette[1:num_groups]
    } else {
      colors <- Polychrome::createPalette(N = num_groups, seedcolors = cb_palette)
    }
  }
  
  if (length(colors) != num_groups) {
    cli::cli_abort("Number of colors does not match the number of unique groups")
  }
  
  names(colors) <- levels(bar_data$group)
  
  p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = colors) +
    geom_hline(yintercept = perc_int_thresh, linetype = "dashed", color = "black", linewidth = 0.5) +
    labs(y = "Total Intensity", fill = grouping_column) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = legend.position)
  
  return(list(
    p = p,
    bar_data = bar_data,
    colors = colors,
    percentile = percentile,
    perc_int_thresh = perc_int_thresh
  ))
}