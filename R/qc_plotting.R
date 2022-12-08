#' Violin plot for QC report
#'
#' Makes a violin plot of per-sample intensities. Samples are grouped by the "groups"
#' argument on the x-axis.
#'
#' @param data A data frame of intensity data, likely normalized. Rows should be
#'   proteins and columns should be samples. Generally, a single element in the
#'   normList slot of the object returned by .
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param sample_labels Optional, a vector of sample labels to use. If not supplied,
#'   defaults to using the column names in the data.
#'
#' @return A ggplot object of the plot.
#'
#' @importFrom ggplot2 ggplot aes geom_violin scale_fill_manual theme_bw theme element_text element_blank ylab
#' @keywords internal
#'
#' @examples
#' # No examples yet
qc_violin_plot <- function(data,
                           groups = NULL,
                           sample_labels = colnames(data)) {



  # Get ordering of samples, by order in which groups are supplied
  group_order <- order(match(groups, unique(groups)))


  # Rename column names to sample names provided
  colnames(data) <- sample_labels

  # Match sample names to groups, for coloring
  sample_group_info <- data.frame(ind = factor(x = colnames(data)[group_order],
                                               levels = colnames(data)[group_order]),
                                  group = factor(x = groups[group_order]))

  # make and return plot
  # Must merge with sample info first, to get order correct
  merge(sample_group_info,
        utils::stack(as.data.frame(data)), sort = F) |>
    ggplot(aes(x = as.factor(.data$ind), y = .data$values, fill = .data$group)) +
    geom_violin(draw_quantiles = c(0.5),
                na.rm = T,
                col = "black") +
    scale_fill_manual(values = colorGroup(groups), limits = unique(groups), name = NULL) +
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
#' determined by the \code{\link{colorGroup}} function. By default, uses
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
#' @keywords internal
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
        sample_group_info) |>
    ggplot(aes(x = .data$x, y = .data$y, color = .data$group, label = .data$ind)) +
    geom_point() +
    scale_color_manual(values = colorGroup(groups), limits = unique(groups), name = NULL) +
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
#' with colors determined by the \code{\link{colorGroup}} function. By default,
#' uses only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_pca_plot
#' @param dist_metric The metric used to define distance for clustering. Default is
#'   "euclidean". See \code{\link[stats:dist]{stats::hclust}} for options.
#' @param clust_method The agglomeration method to use for clustering. Default
#'   is "complete", See \code{\link[stats:hclust]{stats::hclust}} for options.
#'
#' @return A ggplot object of the plot.
#' @keywords internal
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

  if (is.null(sample_labels)) {
    sample_labels <- colnames(data)
  } else {
    colnames(data) <- sample_labels
  }

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
  sample_group_info <- data.frame(text = sample_labels,
                                  group = groups)

  # Make plot
  ggtree::ggtree(hc) %<+% sample_group_info +
    ggtree::layout_dendrogram() +
    ggtree::geom_tippoint(aes(color = .data$group)) +
    ggtree::geom_tiplab(aes(color = .data$group), angle=90, hjust=1, offset = -0.5, show.legend=FALSE) +
    scale_color_manual(values = colorGroup(groups), limits = unique(groups), name = NULL) +
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
#'
#' @return A \code{\link[grid:gTree]{grid::gTree}} object of the ComplexHeatmap,
#'   which can be plotted with \code{\link[grid:grid.draw]{grid::grid.draw}}.
#' @keywords internal
#'
#' @examples
#' # No examples yet
qc_corr_hm <- function(data,
                       groups,
                       sample_labels = NULL) {

  groups <- make_factor(groups)


  if (is.null(sample_labels)) {
    sample_labels <- colnames(data)
  }

  # Calculate correlation matrix
  cor_mat <- stats::cor(data, use = "pairwise.complete.obs", method = "pearson")

  # Get column annotations
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(
    Sample = groups,
    col = list(Sample = colorGroup(groups)),
    annotation_legend_param = list(
      Sample = list(
        title = "Groups",
        at = levels(groups),
        labels = paste(levels(groups))
      )
    )
  )

  # Adjust plot size based on number of samples
  # This also occurs in the write_qc_report function and qc_corr_hm function.
  # If you make changes here, make corresponding changes in those functions
  fontsize <- ifelse(length(groups) < 100, 12, 10)
  if (ncol(data) > 50) {
    height <- 15
    width <- 15
  } else {
    height <- 6
    width <- 6
  }


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
    left_annotation = NULL,
    column_labels = sample_labels,
    row_labels = sample_labels,
    width = grid::unit(width, "in"),
    height = grid::unit(height, "in")
  )

  # Output the Heatmap as a gTree object,
  # Similar to our other plotting fxns, which output
  # ggplot objects
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      hm_corr, heatmap_legend_side = "top",
      ht_gap = grid::unit(2, "cm")
    )
  )
}




#' Create a missing values heatmap
#'
#' Makes a ComplexHeatmap object showing a heatmap of missing values in the
#' input data.
#'
#' @inheritParams qc_pca_plot
#' @param column_sort How should the columns of the heatmap be sorted? Options
#'   are: "cluster"- sort by similarity in missing values, "group"- sort samples
#'   by the grouping variable.
#' @param group_var_name The name of the variable being used for group sorting.
#' @param show_all_proteins Should all proteins be shown in missing value heatmap,
#'  of only those with missing data? Default is F (only those with missing data).
#'
#' @return A \code{\link[grid:gTree]{grid::gTree}} object of the ComplexHeatmap,
#'   which can be plotted with \code{\link[grid:grid.draw]{grid::grid.draw}}.
#' @keywords internal
#'
#' @examples
#' # No examples yet.
#'
qc_missing_hm <- function(data,
                          groups,
                          sample_labels = NULL,
                          column_sort = c("cluster", "group"),
                          group_var_name = "",
                          show_all_proteins = F) {
  # Check args
  column_sort <- rlang::arg_match(column_sort)

  if (is.null(sample_labels)) {
    sample_labels <- colnames(data)
  }

  # Set up the column ordering for the different options
  if (column_sort == "cluster") {
    # If clustering just order things as they are
    cluster <- T
    order <- seq_along(groups)
    title <- "Sorted by similarity"
  } else if (column_sort == "group") {
    # If ordering by group, sort by group
    cluster <- F
    order <- order(groups)
    title <- paste0("Sorted by ", group_var_name)
  } else {
    cli::cli_abort("Invalid value for {.arg column_sort}: cannot be {.val {column_sort}}")
  }

  #Prepare data
  missing <- !is.na(data)
  if (!show_all_proteins) {
    complete = apply(missing, 1, all)
    completeNA = apply(!missing, 1, all)
    missing <- missing[!complete & !completeNA,]
  }

  # Adjust plot size based on number of samples
  # This also occurs in the write_qc_report function and qc_corr_hm function.
  # If you make changes here, make corresponding changes in those functions
  if (ncol(data) > 50) {
    height <- 15
    width <- 15
  } else {
    height <- 6
    width <- 6
  }

  # Then order the groups/batches/etc
  sorted_groups <- groups[order]
  sorted_labels <- sample_labels[order]
  ordered_data <- missing[,order]

  # Set up heatmap annotation
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(
    sample = sorted_groups,
    col = list(sample = colorGroup(groups)),
    annotation_legend_param = list(sample = list(title = "Group",
                                                 at = unique(sorted_groups),
                                                 labels = paste("", unique(sorted_groups)))),
    show_legend = T
  )
  # Plot heatmap
  hm_clust <- ComplexHeatmap::Heatmap(
    ordered_data + 0,
    col = c("white", "black"),
    column_names_side = "top",
    column_title = title,
    show_row_names = FALSE,
    show_column_names = TRUE,
    name = "Status",
    column_names_gp = grid::gpar(fontsize=7),
    heatmap_legend_param = list(at = c(0, 1),
                                labels = c("Missing", "Valid")),
    show_heatmap_legend = T,
    top_annotation = ColAnn,
    cluster_columns = cluster,
    column_labels = sorted_labels,
    width = grid::unit(width, "in"),
    height = grid::unit(height, "in")
  )


  # Output the Heatmap as a gTree object,
  # Similar to our other plotting fxns, which output
  # ggplot objects
  grid::grid.grabExpr(
    ComplexHeatmap::draw(
      hm_clust,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      ht_gap = grid::unit(2, "cm"),
      column_title = "Missing Values"
      )
    )
}
