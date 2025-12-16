#' Generate a QC Boxplot
#'
#' @param data A matrix or data frame of intensities.
#' @param groups A vector indicating class membership (numeric, integer, character, factor, or NULL).
#' @param sample_labels A vector of names to display on the plot (character or NULL).
#' @param title A character string for the plot title.
#' @param text.sizes A numeric vector specifying text sizes for title, x-axis, y-axis, and legend.
#' @param legend.position A character string specifying the legend position (default: "right").
#' @param boxplot_width A numeric value specifying the width of the boxplots (default: dynamic based on groups).
#' @param boxplot_alpha A numeric value specifying the transparency level of the boxplots (default: dynamic).
#' @param plot_margin A grid unit object specifying the plot margins (default: dynamic based on groups).
#' @param colorblind_palette A character vector specifying colorblind-friendly colors to use.
#'
#' @return A ggplot2 object displaying the QC boxplot.
#' @export
qc_boxplot <- function(data,
                       groups = NULL,
                       sample_labels = NULL,
                       title = NULL,
                       text.sizes = c(12, 10, 10, 10),
                       legend.position = "right",
                       boxplot_width = NULL,
                       boxplot_alpha = NULL,
                       plot_margin = NULL,
                       colorblind_palette = NULL) {
  
  # Default colorblind-friendly palette
  base_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                    "#0072B2", "#D55E00", "#CC79A7", "#999999")
  if (is.null(colorblind_palette)) {
    colorblind_palette <- base_palette
  }
  
  # Validate input
  stopifnot(is.matrix(data) || is.data.frame(data))
  
  if (!is.null(groups)) {
    stopifnot(length(groups) == ncol(data))
  }
  if (!is.null(sample_labels)) {
    stopifnot(length(sample_labels) == ncol(data))
  }
  stopifnot(length(text.sizes) == 4)
  
  # Ensure column names exist
  if (is.null(colnames(data))) {
    colnames(data) <- paste0("Sample", seq_len(ncol(data)))
  }
  
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
  groups        <- groups[group_order]
  data          <- data[, group_order, drop = FALSE]
  
  # Create metadata frame
  plot.meta <- data.frame(
    ind    = factor(colnames(data)),
    labels = factor(sample_labels),
    group  = groups
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
  
  # Generate the plot
  p1 <- ggplot(plot.data, aes(x = labels, y = values, fill = group)) +
    geom_boxplot(width = boxplot_width, color = "black", alpha = boxplot_alpha) +
    scale_fill_manual(values = colorblind_palette, name = NULL) +
    labs(y = "Density", x = "", fill = "") +
    theme_gray() +
    theme(
      axis.text.x      = element_text(angle = 45, vjust = 0.9, hjust = 1, size = text.sizes[3]),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = text.sizes[2]),
      axis.text.y      = element_text(size = text.sizes[3]),
      plot.title       = element_text(size = text.sizes[1]),
      plot.margin      = plot_margin,
      legend.position  = legend.position,
      legend.text      = element_text(size = text.sizes[4]),
      legend.key.size  = unit(0.5, "cm")
    ) +
    ggtitle(title)
  
  return(p1)
}
