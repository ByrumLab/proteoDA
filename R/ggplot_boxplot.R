#' Generate a QC Boxplot
#'
#' This function draws a boxplot of sample intensities for quality control.
#' The `data` argument may be a numeric matrix or data.frame of intensities,
#' or a DAList object produced by the proteoDA pipeline. If a DAList is
#' supplied, the function will prefer the normalized data stored in
#' `DAList$data_per_contrast[[contrast]]` when that element exists. If
#' 'contrast' is NULL and data_per_contrast has exactly one element, that
#' element will be used. Otherwise DAList$data is used.
#'
#' @param data A matrix or data.frame of intensities, or a DAList object.
#' @param contrast Optional character string naming the contrast to use when
#'   a DAList is provided. If NULL and data_per_contrast has exactly one
#'   element, that element is used.
#' @param groups A vector indicating class membership (numeric, integer,
#'   character, factor, or NULL). Length must match number of samples.
#' @param sample_labels A character vector of names to display on the plot,
#'   or NULL to use column names from the data.
#' @param title A character string for the plot title.
#' @param text.sizes Numeric vector of length 4 specifying text sizes for
#'   title, y-axis, x-axis labels, and legend text, in that order.
#' @param legend.position Character string specifying legend position
#'   (default: \"right\").
#' @param boxplot_width Numeric width for the boxplots. If NULL a dynamic
#'   value based on number of groups is used.
#' @param boxplot_alpha Numeric transparency for the boxplots. If NULL a
#'   dynamic value based on number of groups is used.
#' @param plot_margin A grid::unit object specifying the plot margins. If
#'   NULL a dynamic margin is used.
#' @param colorblind_palette Character vector of color values to use for
#'   group fills. If NULL a default colorblind friendly palette is used.
#'
#' @return A ggplot2 object displaying the QC boxplot.
#' @export
qc_boxplot <- function(data,
                       contrast = NULL,
                       groups = NULL,
                       sample_labels = NULL,
                       title = NULL,
                       text.sizes = c(12, 10, 10, 10),
                       legend.position = "right",
                       boxplot_width = NULL,
                       boxplot_alpha = NULL,
                       plot_margin = NULL,
                       colorblind_palette = NULL) {
  
  # If a DAList object is provided, select the appropriate data matrix
  # Accept any plain list-like DAList, but not a data.frame
  if (is.list(data) && !is.data.frame(data)) {
    # Try to detect DAList structure
    if (!is.null(data$data_per_contrast) && is.list(data$data_per_contrast) &&
        length(data$data_per_contrast) > 0) {
      if (!is.null(contrast) && nzchar(contrast) &&
          !is.null(data$data_per_contrast[[contrast]])) {
        data_to_plot <- data$data_per_contrast[[contrast]]
      } else if (length(data$data_per_contrast) == 1) {
        data_to_plot <- data$data_per_contrast[[1]]
      } else if (!is.null(data$data)) {
        data_to_plot <- data$data
      } else {
        stop("DAList provided but no usable data found in data_per_contrast or data.")
      }
    } else if (!is.null(data$data)) {
      data_to_plot <- data$data
    } else {
      stop("DAList provided but no usable data found in data_per_contrast or data.")
    }
    data <- data_to_plot
  }
  
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
  
  # Sort groups by levels (keep within-group sample order stable)
  group_order <- order(groups)
  
  # Reorder data accordingly
  sample_labels <- sample_labels[group_order]
  groups        <- groups[group_order]
  data          <- data[, group_order, drop = FALSE]
  
  # Create metadata frame
  plot.meta <- data.frame(
    ind    = factor(colnames(data)),
    labels = factor(sample_labels),
    group  = groups,
    stringsAsFactors = FALSE
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
    plot_margin <- grid::unit(rep(ifelse(unique_groups > 6, 0.3, 0.2), 4), "cm")
  }
  
  # Generate the plot
  p1 <- ggplot2::ggplot(plot.data, ggplot2::aes(x = labels, y = values, fill = group)) +
    ggplot2::geom_boxplot(width = boxplot_width, color = "black", alpha = boxplot_alpha) +
    ggplot2::scale_fill_manual(values = colorblind_palette, name = NULL) +
    ggplot2::labs(y = "Density", x = "", fill = "") +
    ggplot2::theme_gray() +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, vjust = 0.9, hjust = 1, size = text.sizes[3]),
      axis.title.x     = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_text(size = text.sizes[2]),
      axis.text.y      = ggplot2::element_text(size = text.sizes[3]),
      plot.title       = ggplot2::element_text(size = text.sizes[1]),
      plot.margin      = plot_margin,
      legend.position  = legend.position,
      legend.text      = ggplot2::element_text(size = text.sizes[4]),
      legend.key.size  = grid::unit(0.5, "cm")
    ) +
    ggplot2::ggtitle(title)
  
  return(p1)
}

#' QC boxplot for data before normalization
#'
#' Generate a QC boxplot from raw intensity data prior to normalization.
#' The input may be a numeric matrix or data.frame, or a DAList-like object
#' containing a top-level element named 'data'. If the data do not appear to
#' be on the log2 scale, a log2 transform is applied before plotting.
#'
#' Log2 scale detection uses a conservative heuristic based on data range and
#' value distribution. When the scale is ambiguous, the function will apply
#' a log2 transform to ensure comparable distributions.
#'
#' @param data A numeric matrix or data.frame of intensities, or a DAList-like
#'   object (list) containing a 'data' matrix.
#' @param groups Optional vector indicating class membership. Length must match
#'   the number of samples (columns).
#' @param sample_labels Optional character vector of sample labels. Length must
#'   match the number of samples (columns).
#' @param pseudo Optional positive numeric value added before log2 transform to
#'   avoid log2(0). If NULL, a small data-driven constant is chosen.
#' @param ... Additional arguments passed to qc_boxplot.
#'
#' @return A ggplot2 object displaying the QC boxplot of log2-scaled raw data.
#' @export
qc_boxplot_beforeNorm <- function(data,
                                  groups = NULL,
                                  sample_labels = NULL,
                                  pseudo = NULL,
                                  ...) {
  
  # Resolve DAList-like objects
  if (is.list(data) && !is.data.frame(data)) {
    if (!is.null(data$data) &&
        (is.matrix(data$data) || is.data.frame(data$data))) {
      mat <- as.matrix(data$data)
    } else {
      stop("DAList provided but no usable 'data' matrix found.")
    }
  } else {
    if (!(is.matrix(data) || is.data.frame(data))) {
      stop("data must be a matrix, data.frame, or a DAList-like object containing 'data'.")
    }
    mat <- as.matrix(data)
  }
  
  # Basic validation
  if (nrow(mat) < 1 || ncol(mat) < 1) {
    stop("data must have at least one row and one column.")
  }
  
  # Detect log2 scale using a conservative heuristic
  vals <- as.numeric(mat)
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    stop("data contains no finite values.")
  }
  
  q <- stats::quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)
  
  is_log2 <- TRUE
  if (any(vals <= 0, na.rm = TRUE)) {
    is_log2 <- FALSE
  } else if ((q[2] - q[1]) > 50) {
    is_log2 <- FALSE
  } else if (q[2] > 50) {
    is_log2 <- FALSE
  }
  
  # Apply log2 transform if needed
  if (!is_log2) {
    if (is.null(pseudo)) {
      pos_vals <- vals[vals > 0]
      if (length(pos_vals) > 0) {
        pseudo <- min(pos_vals, na.rm = TRUE) / 10
        pseudo <- max(pseudo, 1e-6)
      } else {
        pseudo <- 1
      }
    } else {
      if (!is.numeric(pseudo) || length(pseudo) != 1 || pseudo <= 0) {
        stop("pseudo must be a single positive numeric value.")
      }
    }
    mat <- log2(mat + pseudo)
  }
  
  # Delegate plotting to qc_boxplot
  qc_boxplot(
    data = mat,
    groups = groups,
    sample_labels = sample_labels,
    ...
  )
}

