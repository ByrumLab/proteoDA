#' Generate Total Intensity Barplots by Group
#'
#' This function creates barplots showing total intensity values for samples, grouped by a specified column.
#' It includes an optional percentile threshold to flag potential outliers and uses a color-blind friendly palette.
#'
#' @param DAList A list containing `data` (matrix) and `metadata` (data frame).
#' @param label_column A character string specifying the column name in metadata for sample labels.
#' @param grouping_column A character string specifying the column name in metadata for grouping.
#' @param percentile A numeric value (0-1) for the percentile threshold to flag low-intensity samples.
#' @param colors A vector of colors for groups. Defaults to a color-blind friendly palette (max 12, >12 uses Polychrome).
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
#'                             nrow = NULL,
#'                             ncol = NULL,
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
  
  legend.position <- rlang::arg_match(legend.position,
                                      c("none", "right", "top", "left", "bottom"))
  
  if (!is.numeric(percentile) || any(percentile < 0 | percentile > 1)) {
    cli::cli_abort("Percentile should be a numeric value between 0 and 1.")
  }
  
  # Prepare data
  bar_data <- data.frame(
    label = DAList$metadata[[label_column]],
    tot.int = colSums(DAList$data, na.rm = TRUE),
    tot.num = colSums(!is.na(DAList$data)),
    tot.na = colSums(is.na(DAList$data)),
    group = factor(DAList$metadata[[grouping_column]], 
                   levels = unique(DAList$metadata[[grouping_column]])), 
    stringsAsFactors = FALSE
  )
  
  num_groups <- length(unique(bar_data$group))
  
  # Calculate percentile threshold
  perc_int_thresh <- ifelse(percentile == 0, 0, quantile(bar_data$tot.int, probs = percentile))
  bar_data$tot.int.outlier <- bar_data$tot.int <= perc_int_thresh
  
  # Define color-blind friendly palette
  cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                  "#D55E00", "#CC79A7", "#999999", "#E69F00", "#56B4E9", 
                  "#009E73", "#F0E442")
  
  # Assign colors to unique groups
  if (is.null(colors)) {
    if (num_groups <= length(cb_palette)) {
      colors <- setNames(cb_palette[1:num_groups], levels(bar_data$group))
    } else {
      colors <- setNames(Polychrome::createPalette(N = num_groups, seedcolors = cb_palette), levels(bar_data$group))
    }
  }
  
  if (length(colors) != num_groups) {
    cli::cli_abort("Number of colors does not match the number of unique groups.")
  }
  
  # Create the plot
  p <- ggplot(bar_data, aes(x = label, y = tot.int, fill = group)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors) +
    geom_hline(yintercept = perc_int_thresh, linetype = "dashed", color = "black") +
    labs(
     # title = paste("Total Intensity by", grouping_column),  # Add heading dynamically
      y = "Total Intensity",
      fill = grouping_column
    ) +
    facet_wrap(~ group, scales = "free_x", shrink = F, nrow = nrow, ncol = ncol) + # split by groups
    ggtitle(paste("Total Intensity by", grouping_column)) + # add title above the facets
    geom_text(aes(label=tot.num), position=position_dodge(width=0.9), vjust=-0.5, size = 3) +  # add number of proteins
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = legend.position) +
    guides(fill=guide_legend(title=grouping_column))
  
  return(list(
    plot = p,
    data = bar_data,
    colors = colors,
    threshold = perc_int_thresh
  ))
}
