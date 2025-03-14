#' Generate a Violin Plot for Intensity Data
#'
#' @param data A matrix or data frame containing intensity values.
#' @param groups A vector (numeric, integer, character, or factor) indicating class membership. Defaults to NULL.
#' @param sample_labels A character vector specifying sample names to display on the plot. Defaults to column names of `data`.
#' @param title A character string specifying the plot title. Defaults to an empty string.
#' @param text.sizes A numeric vector of length 4 specifying font sizes for title, x-axis, y-axis, and legend text, respectively.
#' @param legend.position A character string specifying legend position (e.g., "right", "left", "none"). Defaults to "right".
#' @param color_palette A character vector specifying a colorblind-friendly palette. If NULL, a default base palette is used based on RColorBrewer
#' @return A ggplot object representing the violin plot.
#' @import ggplot2 testthat utils
#' 
#' @example
#' \dontrun {
#' # example of violin plot colored by groups 
#' violin1 <- qc_violin(data            = results$data,
#'                      groups          = results$metadata$group,
#'                      sample_labels   = results$metadata$sample,
#'                      title           = "",
#'                      text.sizes      = c(12, 10, 10, 10),
#'                      legend.position = "right",
#'                      color_palette   = c("#66c2A5","#8DA0CB", "#FC8D62")  
#'                      )
#' 
#' 
#' }
#' 
#' 
#' @export
qc_violin <- function(data,
                      groups = NULL,
                      sample_labels = NULL,
                      title = "",
                      text.sizes = c(12, 10, 10, 10),
                      legend.position = "right",
                      color_palette = NULL) {

  # Argument checks
  stopifnot(is.matrix(data) || is.data.frame(data))
  stopifnot(is.null(groups) || is.numeric(groups) || is.integer(groups) ||
              is.character(groups) || is.factor(groups))
  stopifnot(is.null(sample_labels) || is.character(sample_labels))
  stopifnot(is.character(title) || is.numeric(title) || is.integer(title))
  stopifnot(is.numeric(text.sizes) && length(text.sizes) == 4)

  # Prepare arguments
  if (is.null(groups)) groups <- rep("1", ncol(data))
  groups <- proteoDA:::make_factor(as.character(groups), prefix = NULL)

  if (is.null(sample_labels)) sample_labels <- colnames(data)

  # Ensure correct lengths
  stopifnot(length(groups) == ncol(data))
  stopifnot(length(sample_labels) == ncol(data))

  # Reorder based on group levels
  group_order <- order(groups)
  sample_labels <- sample_labels[group_order]
  groups <- groups[group_order]
  data <- data[, group_order]

  # Prepare metadata for plotting
  plot.meta <- data.frame(
    ind = proteoDA:::make_factor(colnames(data)),
    labels = proteoDA:::make_factor(sample_labels),
    group = groups
  )

  # Convert data to long format
  plot.data <- merge(plot.meta, utils::stack(as.data.frame(data)), sort = FALSE)

  # Define color palette
  unique_groups <- unique(groups)
  num_colors <- length(unique_groups)
  base_palette <- RColorBrewer::brewer.pal(min(num_colors, 8), "Set2")
  if (is.null(color_palette)) {
    color_palette <- base_palette
  }

  # Create the violin plot
  p1 <- ggplot(plot.data, aes(x = labels, y = values, fill = group)) +
    geom_violin(draw_quantiles = c(0.5), na.rm = TRUE, color = "black") +
    scale_fill_manual(values = color_palette[1:num_colors], name = NULL) +
    labs(y = "Density", x = "", fill = "") +
    ggtitle(title) +
    theme_gray() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1, size = text.sizes[3]),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = text.sizes[2]),
      axis.text.y = element_text(size = text.sizes[3]),
      plot.title = element_text(size = text.sizes[1]),
      plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
      legend.position = legend.position,
      legend.text = element_text(size = text.sizes[4]),
      legend.key.size = unit(0.5, 'cm')
    )

  return(p1)
}
