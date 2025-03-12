#' @title PCA Plot for Subgroups
#' @description Generates PCA plots for different subgroups based on provided metadata.
#' @param DAList A list containing data and metadata for differential analysis. Should include a 'data' matrix and 'metadata' data frame.
#' @param grouping_column The name of the column in the metadata that defines the grouping for the samples.
#' @param label_column The name of the column in the metadata for sample labels. Default: NULL
#' @param group_color A color to represent the groups in the plot. Default: NULL (will use default color).
#' @param top The number of the most variable proteins to include in the PCA analysis. Default: 500
#' @param standardize Logical indicating whether to standardize the input data to have mean 0 and standard deviation 1. Default: TRUE
#' @param pca_axes A numeric vector of length 2 specifying which PCA axes to plot. Default: c(1, 2)
#' @param text.sizes A numeric vector specifying text sizes for different plot elements: title, axis titles, axis labels, and legend. Default: c(14, 14, 12, 10)
#' @param point.size Size of points in the PCA plot. Default: 3
#' @param label.size Size of sample labels in the plot. Default: 4
#' @param max.labels Maximum number of labels to show in the plot. Default: 100
#' @param legend.position Position of the legend in the plot. Default: 'right'
#' @param show.plot Logical indicating whether to display the plot. Default: TRUE
#' @return A list containing PCA plots for each subgroup and the PCA analysis results.
#' @details This function allows users to visualize the PCA results for different subgroups in the dataset, providing insights into the variance and structure of the data.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage
#'  qc_pca_plot7_subgroups(DAList, "group", "sample_label")
#'  }
#' }
#' @seealso
#'  \code{\link[proteoDA]{validate_DAList}}, \code{\link[proteoDA]{make_factor}}, \code{\link[cli]{cli_abort}}
#' @rdname qc_pca_plot7_subgroups
#' @export
#' @importFrom proteoDA validate_DAList make_factor binfcolors
#' @importFrom cli cli_abort
qc_pca_plot7_subgroups <- function(DAList,
                                   grouping_column,
                                   label_column = NULL,
                                   group_color = NULL,
                                   top = 500,
                                   standardize = TRUE,
                                   pca_axes = c(1, 2),
                                   text.sizes = c(14, 14, 12, 10),
                                   point.size = 3,
                                   label.size = 4,
                                   max.labels = 100,
                                   legend.position = "right",
                                   show.plot = TRUE,
                                   max_pc = NULL) {

  DAList <- proteoDA:::validate_DAList(x = DAList)

  grouping_column <- check_grouping_column(metadata = DAList$metadata,
                                           grouping_column = grouping_column)
  groups <- DAList$metadata[, grouping_column]

  if (!is.factor(groups)) {
    groups <- proteoDA:::make_factor(x = groups)
  } else {
    groups <- droplevels(x = groups)
  }

  if (is.null(label_column)) {
    sample_labels <- colnames(DAList$data)
  } else {
    label_column <- check_label_column(metadata = DAList$metadata,
                                       label_column = label_column)
    sample_labels <- DAList$metadata[, label_column]
  }

  if (is.null(group_color)) {
    group_color <- proteoDA:::binfcolors[1]
  }

  if (length(group_color) != 1L) {
    cli::cli_abort("The group_color argument should be one color.")
  }

  out <- lapply(levels(groups), function(x) {
    plot_colors <- rep("gray85", length(unique(groups)))
    names(plot_colors) <- levels(groups)

    idx <- which(names(plot_colors) %in% x)
    plot_colors[idx] <- group_color

    sample_labels2 <- rep("", length(sample_labels))
    sample_labels2[which(groups %in% x)] <- sample_labels[which(groups %in% x)]

    qc_pca_plot7(data            = DAList$data,
                 groups          = groups,
                 sample_labels   = sample_labels2,
                 text.sizes      = text.sizes,
                 point.size      = point.size,
                 label.size      = label.size,
                 max.labels      = max.labels,
                 top             = top,
                 standardize     = standardize,
                 pca_axes        = pca_axes,
                 colors          = plot_colors,
                 title           = paste0("colored by ", grouping_column, ": ", x),
                 legend.position = legend.position,
                 show.plot       = show.plot,
                 max_pc          = max_pc)
  })
  names(out) <- levels(groups)

  out
}

#' @title PCA Plot
#' @description Generates PCA plots from normalized intensity data.
#' @param data A data frame of intensity data, likely normalized. Rows should be proteins and columns should be samples.
#' @param groups A character or factor vector, listing the group(s) the samples belong to. Default: NULL
#' @param sample_labels Optional, a vector of sample labels to use. If not supplied, defaults to using the column names in the data. Default: NULL
#' @param text.sizes A numeric vector specifying text sizes for different plot elements: title, axis titles, axis labels, and legend. Default: c(14, 14, 12, 10)
#' @param point.size Size of points in the PCA plot. Default: 3
#' @param label.size Size of sample labels in the plot. Default: 4
#' @param max.labels Maximum number of labels to show in the plot. Default: 100
#' @param top The number of most variable proteins to use for the PCA analysis. Default: 500
#' @param standardize Logical indicating whether to standardize the input data before PCA. Default: TRUE
#' @param pca_axes A numeric vector of length 2 specifying which PCA axes to plot. Default: c(1, 2)
#' @param colors A vector of colors corresponding to groups. Default: NULL
#' @param title Title for the plot. Default: NULL
#' @param legend.position Position of the legend in the plot. Default: 'right'
#' @param show.plot Logical indicating whether to display the plot. Default: TRUE
#' @return A list containing the PCA plot, the PCA analysis results, and additional information.
#' @details This function performs PCA on the provided data and generates a corresponding plot.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage
#'  qc_pca_plot7(data, groups)
#'  }
#' }
#' @seealso
#'  \code{\link[rlang]{arg_match}}, \code{\link[cli]{cli_abort}},
#'  \code{\link[proteoDA]{make_factor}}, \code{\link[proteoDA]{colorGroup}}, \code{\link[proteoDA]{rowVars}},
#'  \code{\link[stats]{prcomp}}, \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{aes}},
#'  \code{\link[ggrepel]{geom_label_repel}}, \code{\link[ggplotify]{as.ggplot}}
#' @rdname qc_pca_plot7
#' @export
#' @importFrom rlang arg_match
#' @importFrom cli cli_abort
#' @importFrom proteoDA make_factor colorGroup rowVars
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplotify as.ggplot
qc_pca_plot7 <- function(data,
                         groups = NULL,
                         sample_labels = NULL,
                         text.sizes = c(14, 14, 12, 10),
                         point.size = 3,
                         label.size = 4,
                         max.labels = 100,
                         top = 500,
                         standardize = TRUE,
                         pca_axes = c(1, 2),
                         colors = NULL,
                         title = NULL,
                         max_pc = NULL,
                         legend.position = "right",
                         show.plot = TRUE) {

  # Legend handling
  legend.position <- rlang::arg_match(arg = legend.position,
                                      values = c("none", "top", "right", "bottom", "left"),
                                      multiple = FALSE)
  show.legend <- ifelse(legend.position == "none", FALSE, TRUE)

  if (is.null(title)) {
    title <- ""
  }

  if (is.null(groups)) {
    groups <- rep(1, ncol(data))
  }

  if (length(groups) != ncol(data)) {
    cli::cli_abort("groups should be the same length as the number of columns in data")
  }

  if (!is.factor(groups)) {
    groups <- proteoDA:::make_factor(as.character(groups), prefix = NULL)
  }
  groups <- droplevels(x = groups)

  # Plot colors
  if (is.null(colors)) {
    colors <- proteoDA:::colorGroup(levels(groups))
  }

  if (!all(length(colors) == length(levels(groups)))) {
    cli::cli_abort("length of colors and levels of groups do not match.")
  }
  names(colors) <- levels(groups)

  if (is.null(max.labels)) {
    max.labels = ncol(data)
  }

  if (is.null(sample_labels)) {
    sample_labels = colnames(data)
  }
  sample_labels <- proteoDA:::make_factor(as.character(sample_labels))

  # Top variable proteins
  dat <- data[!apply(is.na(data), 1, any), ]
  dat <- dat[order(proteoDA:::rowVars(as.matrix(dat)), decreasing = TRUE), ]

  if (is.null(top)) { top <- nrow(dat) }
  top2 <- ifelse(nrow(dat) >= top, top, nrow(dat))
  dat <- dat[1:top2, ]

  # Standardize the data
  if (standardize) {
    dat <- t(scale(x = t(dat), center = TRUE, scale = TRUE))
  }

  # Perform PCA
  pca <- stats::prcomp(t(dat), scale = FALSE)
  pca$summary <- summary(pca)$importance

  # Prepare data for plotting
  pca_data <- as.data.frame(pca$x[, pca_axes])
  pca_data$ind <- rownames(pca_data)
  colnames(pca_data) <- c("x", "y", "ind")

  plot_data <- merge(pca_data, sample_group_info, by = "ind")

  ########################
  ## PCA PLOT
  ########################
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$x,
                                               y = .data$y,
                                               color = .data$groups,
                                               label = .data$labels)) +
    geom_point(show.legend = show.legend, size = point.size) +
    scale_color_manual(values = colors[groups], name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    theme_gray() +
    theme(axis.title.x = element_text(size = text.sizes[2]),
          axis.text.x = element_text(size = text.sizes[3], hjust = 0.5,
                                     margin = margin(t = 0.3, b = 0.15, unit = "cm")),
          axis.title.y = element_text(size = text.sizes[2]),
          axis.text.y = element_text(size = text.sizes[3], vjust = 0.5,
                                     margin = margin(r = 0.3, l = 0.15, unit = "cm"))) +
    theme(axis.ticks.y = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.y = unit(0.25, "cm"),
          axis.ticks.x = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.x = unit(0.25, "cm")) +
    theme(panel.grid.major = element_line(linetype = "dotted", colour = "white", linewidth = 0.3),
          panel.grid.minor = element_line(linetype = "dotted", colour = "white", linewidth = 0.3),
          panel.border = element_rect(linetype = "solid", fill = NA, colour = "gray90", linewidth = 0.3)) +
    theme(plot.title = element_text(size = text.sizes[1], face = "bold", hjust = 0)) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
    theme(legend.position = legend.position,
          legend.text = element_text(size = text.sizes[4], vjust = 0.5)) +
    theme(legend.key.size = unit(0.5, 'cm'),
          legend.key = element_rect(fill = NA)) +
    theme(legend.margin = margin(c(legend.margin), unit = "cm")) +
    xlab(paste0("PC", pca_axes[1], " (", round(pca$summary["Proportion of Variance", pca_axes[1]] * 100, 2), " %)")) +
    ylab(paste0("PC", pca_axes[2], " (", round(pca$summary["Proportion of Variance", pca_axes[2]] * 100, 2), " %)")) +
    ggtitle(title)

  if (nrow(plot_data) < max.labels) {
    p <- p + ggrepel::geom_label_repel(show.legend = F,
                                       label.size = NA,
                                       fill = NA,
                                       size = label.size,
                                       max.time = 5,
                                       max.iter = 1e9,
                                       max.overlaps = 10,
                                       force = 1,
                                       force_pull = 0.5,
                                       label.padding = 0.25,
                                       point.padding = 0.25)
  }
  p <- ggplotify::as.ggplot(p)

  if (show.plot) { plot(p) }

  ####################
  ## SCREE PLOT
  ####################
  scree <- qc_pca_scree_plot7(pca = pca,
                              max_pc = max_pc,
                              text.sizes = text.sizes[1:3])

  invisible(list(p = p, pca = pca, plot_data = plot_data,
                 sample_group_info = sample_group_info,
                 param = list(top = top, standardize = standardize),
                 scree = scree))
}

#' @title Scree Plot for PCA
#' @description Generates a scree plot showing the proportion of variance for each principal component.
#' @param pca A prcomp object resulting from PCA analysis.
#' @param max_pc The maximum number of PCs to display in the scree plot. Default: NULL (shows all PCs)
#' @param text.sizes A numeric vector specifying text sizes for different plot elements: title, axis titles, axis labels. Default: c(16, 14, 14)
#' @return A list containing the scree plot and the data used for plotting.
#' @details This function visualizes the variance explained by each principal component in the PCA analysis.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage
#'  qc_pca_scree_plot7(pca)
#'  }
#' }
#' @seealso
#'  \code{\link[cli]{cli_abort}}
#' @rdname qc_pca_scree_plot7
#' @export
#' @importFrom cli cli_abort
qc_pca_scree_plot7 <- function(pca,
                               max_pc = NULL,
                               text.sizes = c(16, 14, 14)) {

  if (!is(pca, "prcomp")) {
    cli::cli_abort("prcomp argument must be of class prcomp. see qc_pca_plot7 pca output.")
  }

  if (length(text.sizes) != 3L) {
    cli::cli_abort(c("the text.sizes argument must be a vector of length 3
                   with positive numeric values (title, axis title, axis text)."))
  }

  is_cex <- sapply(text.sizes, check_num)
  if (!all(is_cex)) {
    cli::cli_abort("the text.sizes argument must be a vector of length 3
                 with positive numeric values (title, axis title, axis text)")
  }

  if (is.null(max_pc)) {
    max_pc <- ncol(pca$summary)
  }

  if (!check_num(x = max_pc)) {
    cli::cli_abort(c("{.arg max_pc} argument is invalid.",
                     "the max number of PCs displayed must be a number between
                  {.val {1}} and {.val {nrow(data)}}."))
  }

  num_pc <- ifelse(max_pc > ncol(pca$summary), ncol(pca$summary), max_pc)

  data <- data.frame(PC = c(1:ncol(pca$summary)),
                     PropVar = c(pca$summary["Proportion of Variance", ]),
                     check.rows = F,
                     check.names = F,
                     fix.empty.names = F,
                     stringsAsFactors = F)
  data$PC <- proteoDA:::make_factor(as.character(data$PC))

  scree_data <- data[1:num_pc, ]

  # Scree plot
  p <- ggplot2::ggplot(data = scree_data, ggplot2::aes(x = .data$PC, y = .data$PropVar, group = 1)) +
    geom_line(linetype = "solid", color = "black", linewidth = 1.2) +
    geom_point(size = 3) +
    theme(axis.title.x = element_text(size = text.sizes[2]),
          axis.text.x = element_text(size = text.sizes[3], hjust = 0.5,
                                     margin = margin(t = 0.3, b = 0.15, unit = "cm")),
          axis.title.y = element_text(size = text.sizes[2]),
          axis.text.y = element_text(size = text.sizes[3], vjust = 0.5,
                                     margin = margin(r = 0.3, l = 0.15, unit = "cm"))) +
    theme(axis.ticks.y = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.y = unit(0.25, "cm"),
          axis.ticks.x = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.x = unit(0.25, "cm")) +
    theme(panel.grid.major = element_line(linetype = "dotted", colour = "white", linewidth = 0.4),
          panel.grid.minor = element_line(linetype = "dotted", colour = "white", linewidth = 0.4),
          panel.border = element_rect(linetype = "solid", fill = NA, colour = "gray90",
                                      linewidth = 0.5)) +
    theme(plot.title = element_text(size = text.sizes[1], face = "bold", hjust = 0)) +
    theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")) +
    xlab(paste0("Principle Component (PC)")) +
    ylab(paste0("Proportion of Variance")) +
    ggtitle(paste0("Scree Plot: Showing PC1 - PC", num_pc, " of ", nrow(data), " PCs"))

  out <- list(p = p, scree_data = data, max_pc = max_pc)

  out
}
