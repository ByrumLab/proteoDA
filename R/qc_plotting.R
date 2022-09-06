#' Violin plot for QC report
#'
#' Makes a violin plot of per-sample intensities. Samples are grouped by the "groups"
#' argument on the x-axis.
#'
#' @param data A data frame of intensity data, likely normalized. Rows should be
#'   proteins and columns should be samples. Generally, a single element in the
#'   normList slot of the object returned by \code{\link{process_data}}.
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param sample_labels Optional, a vector of sample labels to use. If not supplied,
#'   defaults to using the column names in the data.
#'
#' @return A ggplot object of the plot.
#'
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_manual theme_bw theme element_text element_blank ylab
#' @export
#'
#' @examples
#' # No examples yet
qc_violin_plot <- function(data,
                           groups = NULL,
                           sample_labels = colnames(data)) {

  # Get alphabetical order of groups
  group_order <- sort.int(groups, index.return = T)$ix


  # Rename column names to sample names provided
  colnames(data) <- sample_labels

  # Match sample names to groups, for coloring
  sample_group_info <- data.frame(ind = factor(x = colnames(data)[group_order], levels = colnames(data)[group_order]),
                                  group = as.factor(groups[group_order]))

  # make and return plot
  # Must merge with sample info first, to get order correct
  merge(sample_group_info,
        utils::stack(as.data.frame(data)), sort = F) %>%
    ggplot(aes(x = as.factor(.data$ind), y = .data$values, fill = .data$group)) +
    geom_violin(draw_quantiles = c(0.5),
                na.rm = T,
                col = "black") +
    scale_fill_manual(values = colorGroup2(sample_group_info$group), name = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank()) +
    ylab("Normalized intensity")
}



#' PCA plot for QC report
#'
#' Performs and then plots a principal component analysis of sample intensities.
#' Samples are colored according to the groups argument, with colors
#' determined by the \code{\link{colorGroup2}} function. By default, uses
#' only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_violin_plot
#' @param top_proteins The number of most variable proteins to use for the analysis.
#'   Default is 500.
#' @param standardize Should input data be standardized to a mean of 0 and std.dev of
#'   1? If input data are not yet standardized, should be TRUE. Default is TRUE.
#' @param pca_axes A numeric vector of length 2 which lists the PC axes to plot.
#'   Default is c(1,2), to plot the first two principal components.
#'
#' @return A ggplot object of the plot.
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual theme_bw theme xlab ylab element_text
#'
#' @examples
#' # No examples yet
qc_pca_plot <- function(data,
                        groups = NULL,
                        sample_labels = colnames(data),
                        top_proteins = 500,
                        standardize = TRUE,
                        pca_axes = c(1, 2)) {

  # Rename column names to sample names provided
  colnames(data) <- sample_labels

  # Match sample names to groups, for coloring
  sample_group_info <- data.frame(ind = colnames(data),
                                  group = groups)

  # Get top variable proteins
  pca_data <- data[!apply(is.na(data), 1, any), ]
  o <- order(rowVars(as.matrix(pca_data)), decreasing = TRUE)
  pca_data <- pca_data[o, ]
  top <- ifelse(nrow(pca_data) >= top_proteins, top_proteins, nrow(pca_data))
  pca_data <- pca_data[1:top, ]

  ## center/scale rows (proteins) mean=0;stdev=1
  if (standardize) {
    pca_data <- t(scale(x = t(pca_data), center = TRUE, scale = TRUE))
  }

  ## Run the PCA
  pca <- stats::prcomp(t(pca_data), scale = FALSE)
  pca$summary <- summary(pca)$importance

  # Prep the data needed for plotting
  plot_data <- as.data.frame(pca$x[, pca_axes])
  plot_data$ind <- rownames(plot_data)
  colnames(plot_data) <- c("x", "y", "ind")

  # Make the plot
  plot <- merge(plot_data,
        sample_group_info) %>%
    ggplot(aes(x = .data$x, y = .data$y, color = .data$group, label = .data$ind)) +
    geom_point() +
    scale_color_manual(values = colorGroup2(groups)[groups], name = NULL) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(paste0("PC", pca_axes[1], " (", round(pca$summary["Proportion of Variance", pca_axes[1]] * 100, 2), " %)")) +
    ylab(paste0("PC", pca_axes[2], " (", round(pca$summary["Proportion of Variance", pca_axes[2]] * 100, 2), " %)"))

  # When fewer than 50 points,
  # add labels
  if (nrow(plot_data) < 50) {
    plot <- plot +
      ggrepel::geom_label_repel(show.legend = F, label.size = NA, fill = NA)
  }

  plot
}


#' Dendrogram for QC report
#'
#' Performs and then plots a hierarchical clustering analysis of sample intensities
#' across samples. Sample labels are colored according to the "groups" argument,
#' with colors determined by the \code{\link{colorGroup2}} function. By default,
#' uses only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_pca_plot
#' @param dist_metric The metric used to define distance for clustering. Default is
#'   "euclidean". See \code{\link[stats:dist]{stats::hclust}} for options.
#' @param clust_method The agglomeration method to use for clustering. Default
#'   is "complete", See \code{\link[stats:hclust]{stats::hclust}} for options.
#'
#' @return A ggplot object of the plot.
#' @export
#'
#' @importFrom ggplot2 aes theme element_text margin
#' @importFrom ggtree %<+%
#'
#' @examples
#' # No examples yet

qc_dendro_plot <- function(data,
                          groups = NULL,
                          sample_labels = NULL,
                          top_proteins = 500,
                          standardize = TRUE,
                          dist_metric = "euclidean",
                          clust_method = "complete") {

  # Check arguments
  dist_metric <- rlang::arg_match(
    arg = dist_metric,
    values = c(
      "euclidean", "maximum", "manhattan",
      "canberra", "binary", "minkowski"
    ),
    multiple = FALSE
  )

  clust_method <- rlang::arg_match(
    arg = clust_method,
    values = c(
      "ward.D", "ward.D2", "single",
      "complete", "average", "mcquitty",
      "median", "centroid"
    ), multiple = FALSE
  )

  # Get top variable proteins
  o <- order(rowVars(as.matrix(data)), decreasing = TRUE)
  data <- data[o, ]
  top <- ifelse(nrow(data) >= top_proteins, top_proteins, nrow(data))
  data <- data[1:top, ]

  ## center/scale rows (proteins) mean=0;stdev=1
  if (standardize) {
    data <- t(scale(x = t(data), center = TRUE, scale = TRUE))
  }

  # Do clustering
  hc <- stats::hclust(stats::dist(t(data), method = dist_metric), method = clust_method)

  # Get sample info for colors
  sample_group_info <- data.frame(label = colnames(data),
                                  group = groups)

  # Make plot
  ggtree::ggtree(hc) %<+% sample_group_info +
    ggtree::layout_dendrogram() +
    ggtree::geom_tippoint(aes(color = .data$group)) +
    ggtree::geom_tiplab(aes(color = .data$group), angle=90, hjust=1, offset = -0.5, show.legend=FALSE) +
    scale_color_manual(values = colorGroup2(groups)[groups], name = NULL) +
    ggtree::theme_dendrogram(plot.margin=margin(6,6,80,6)) +
    theme(plot.title = element_text(hjust = 0.5))

}



#' Correlation heatmap for QC report
#'
#' Plots a ComplexHeatmap object showing the pairwise correlations of intensities
#' across samples. Samples are grouped by similarity on the x- and y-axis, with
#' labels colored by group and batch.
#'
#' @inheritParams qc_pca_plot
#' @param batch A character or factor vector, listing the batch(es) the samples
#'   belong to.
#'
#' @return Invisibly returns the correlation matrix created by
#'   \code{\link[stats:cor]{stats::cor}}.
#' @export
#'
#' @examples
#' # No examples yet
qc_corr_hm <- function(data,
                      groups, batch = NULL, sampleLabels = NULL) {

  groups <- make_factor(groups)
  if (is.null(batch)) {
    batch <- c(rep("1", ncol(data)))
  }
  batch <- make_factor(as.character(batch), prefix = NULL)
  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }

  # Calculate correlation matrix
  cor_mat <- stats::cor(data, use = "pairwise.complete.obs", method = "pearson")

  # Get column annotations
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(
    Sample = groups,
    col = list(Sample = colorGroup2(groups)),
    annotation_legend_param = list(
      Sample = list(
        title = "Groups",
        at = levels(groups),
        labels = paste(levels(groups))
      )
    )
  )
  # Get row annotations
  RowAnn <- ComplexHeatmap::rowAnnotation(
    Batch = batch, col = list(Batch = colorBatch(batch)),
    annotation_legend_param = list(
      Batch = list(
        title = "Batch",
        at = levels(batch),
        labels = paste("Batch", levels(batch))
      )
    )
  )

  fontsize <- ifelse(length(groups) < 100, 12, 10)
  # make the complex heatmap plot object
  hm_corr <- ComplexHeatmap::Heatmap(cor_mat,
    name = "Pearson correlation", border = TRUE,
    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat)) / 7)),
      colors = c("#F7FBFF", "#DEEBF7",
                 "#C6DBEF", "#9ECAE1",
                 "#6BAED6", "#4292C6",
                 "#2171B5", "#084594"), # originally from RColorBrewer::brewer.pal(8, "Blues")
      transparency = 0.6
    ),
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "horizontal",
      legend_width = grid::unit(5, "cm"),
      title_position = "topcenter"
    ),
    column_names_gp = grid::gpar(fontsize = fontsize),
    row_names_gp = grid::gpar(fontsize = fontsize),
    top_annotation = ColAnn,
    left_annotation = RowAnn,
    column_labels = sampleLabels,
    row_labels = sampleLabels,

  )

  # Get original graphics pars
  orig_par <- graphics::par(no.readonly = TRUE)
  # Schedule them to be reset when we exit the function
  on.exit(graphics::par(orig_par), add = TRUE)

  # Set params, draw plot
  graphics::par(mar = c(10, 10, 10, 10))
  ComplexHeatmap::draw(hm_corr, heatmap_legend_side = "top",
                       ht_gap = grid::unit(2, "cm"))

  invisible(cor_mat)
}

