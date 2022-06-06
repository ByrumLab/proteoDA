#' Boxplot for QC report
#'
#' Makes a boxplot of per-sample intensities. Samples are grouped by the "groups"
#' argument on the x-axis. Colors are determined by the \code{\link{colorGroup2}}
#' function.
#'
#' @param data A data frame of intensity data, likely normalized. Rows should be
#'   proteins and columns should be samples.
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param sampleLabels Optional, a vector of sample labels to use. If not supplied,
#'   defaults to using the column names in the data.
#' @param title Optional, a plot title.
#' @param legend Include a legend in the plot? Default is TRUE.
#'
#' @return Invisibly returns the list object created by
#'   \code{\link[graphics:boxplot]{graphics::boxplot}}.
#' @export
#'
#' @examples
#' # No examples yet
qc_boxplot <- function(data,
                       groups = NULL, sampleLabels = NULL,
                       title = NULL, legend = TRUE) {
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)

  # Reorder group affiliations
  # Since groups is an ordered factor (see make_factor()),
  # this puts them in the order of their levels, then sorts by name.
  # Does not necessarily sort alphabetically.
  group_order <- sort.int(groups, index.return = T)$ix

  # Then, reorder the input data and groups
  data <- data[,group_order]
  sampleLabels <- sampleLabels[group_order]
  groups <- groups[group_order]

  # Get original graphics pars
  orig_par <- graphics::par(no.readonly = TRUE)
  # Schedule them to be reset when we exit the function
  on.exit(graphics::par(orig_par), add = TRUE)

  ## plot margins
  b_mar <- left_margin(x = sampleLabels)
  r_mar <- right_margin(x = groups)
  if (!legend) r_mar <- 1


  # Setup graphics parameters:
  graphics::par(mar = c(b_mar, 6, 3, r_mar),
                oma = c(0.5, 0, 0.5, r_mar/2 + 3))

  # Make the plot
  box <- graphics::boxplot(data,
                           col = colorGroup2(groups)[groups], names = sampleLabels,
                           notch = FALSE, horizontal = FALSE, outline = FALSE,
                           las = 2, cex.axis = 1, cex.labs = 1, cex = 1
  )
  graphics::mtext(side = 2, text = "Normalized Intensity", font = 1, line = 3, cex = 1.2)
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.1)
  if (legend) {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
           bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
           col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1)
  }

  invisible(box)
}

#' Violin plot for QC report
#'
#' Makes a violin plot of per-sample intensities. Samples are grouped by the "groups"
#' argument on the x-axis. Colors are determined by the \code{\link{colorGroup2}}
#' function.
#'
#' @inheritParams qc_boxplot
#'
#' @return  Invisibly returns the list object created by
#'   \code{\link[vioplot:vioplot]{vioplot::vioplot}}.
#' @export
#'
#' @examples
#' # No examples yet
qc_violin_plot <- function(data,
                       groups = NULL, sampleLabels = NULL,
                       title = NULL, legend = TRUE) {
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)

  ## convert data data.frame to a named list
  plotData <- list()
  for (k in 1:ncol(data)) {
    plotData[[k]] <- data[, k]
  }
  names(plotData) <- colnames(data)

  # Reorder group affiliations
  # Since groups is an ordered factor (see make_factor()),
  # this puts them in the order of their levels, then sorts by name.
  # Does not necessarily sort alphabetically.
  group_order <- sort.int(groups, index.return = T)$ix

  # Then, reorder the input data and groups
  plotData <- plotData[group_order]
  sampleLabels <- sampleLabels[group_order]
  groups <- groups[group_order]

  # Get original graphics pars
  orig_par <- graphics::par(no.readonly = TRUE)
  # Schedule them to be reset when we exit the function
  on.exit(graphics::par(orig_par), add = TRUE)

  ## plot margins
  b_mar <- left_margin(x = sampleLabels)
  r_mar <- right_margin(x = groups)
  if (!legend) r_mar <- 1

  # Setup graphics parameters:
  graphics::par(mar = c(b_mar, 6, 3, r_mar),
                oma = c(0.5, 0, 0.5, r_mar/2 + 3))



  # Make the plot
  vio <- vioplot::vioplot(plotData,
                          col = colorGroup2(groups)[groups], names = sampleLabels,
                          main = "", ylab = "", xlab = "",
                          las = 2, font = 1, cex.axis = 1, cex.labs = 1, cex = 1
  )
  graphics::mtext(side = 2, text = "Normalized Intensity", line = 3, cex = 1.2)
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.1)
  if (legend) {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
           bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
           col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1)
  }

  invisible(vio)
}



#' PCA plot for QC report
#'
#' Performs and then plots a principal component analysis of sample intensities.
#' Samples are colored according to the groups argument, with colors
#' determined by the \code{\link{colorGroup2}} function. By default, uses
#' only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_boxplot
#' @param top The number of most variable proteins to use for the analysis.
#'   Default is 500.
#' @param stdize Should input data be standardized to a mean of 0 and std.dev of
#'   1? If input data are not yet standardized, should be TRUE. Default is TRUE.
#' @param dims A numeric vector of length 2 which lists the PC axes to plot.
#'   Default is c(1,2), to plot the first two principal components.
#' @param cex.dot Character expansion factor for plotting points. Default is 2.
#' @param xlim Optional. Custom x-axis limits of the plot. By default, the min is
#'   10% below the min PC score on the x-axis, and the max is above the max PC
#'   score on the x-axis, with some padding for sample labels.
#' @param ylim Optional. Custom y-axis limits of the plot. Default is 10% above
#'   and below the min/max PC score on the y-axis.
#'
#' @return Invisibly returns a list with two elements: \enumerate{
#'   \item "pca"- The list returned by
#'     \code{\link[stats:princomp]{stats::princomp}}.
#'   \item "pca.data"- The data matrix that went into the PCA.
#' }
#' @export
#'
#' @examples
#' # No examples yet
qc_pca_plot <- function(data,
                    groups = NULL, sampleLabels = NULL,
                    title = NULL, legend = TRUE,
                    top = 500, stdize = TRUE,
                    dims = c(1, 2), cex.dot = 2,
                    xlim = NULL, ylim = NULL) {

  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)


  # Get top variable proteins
  data <- data[!apply(is.na(data), 1, any), ]
  o <- order(rowVars(as.matrix(data)), decreasing = TRUE)
  data <- data[o, ]
  top <- ifelse(nrow(data) >= top, top, nrow(data))
  data <- data[1:top, ]

  ## center/scale rows (proteins) mean=0;stdev=1
  if (stdize) {
    data <- t(scale(x = t(data), center = TRUE, scale = TRUE))
  }

  ## PCA
  pca <- stats::prcomp(t(data), scale = FALSE)
  pca$summary <- summary(pca)$importance


  # Get original graphics pars
  orig_par <- graphics::par(no.readonly = TRUE)
  # Schedule them to be reset when we exit the function
  on.exit(graphics::par(orig_par), add = TRUE)

  # Set plot margins
  b_mar <- left_margin(x = sampleLabels)
  r_mar <- right_margin(x = groups)
  if (!legend) r_mar <- 1

  graphics::par(mar = c(b_mar, 6, 3, r_mar),
                oma = c(0.5, 0, 0.5, r_mar/2 + 3))

  # Set plot limits
  longest_label <- max(nchar(sampleLabels))
  label_cex <- ifelse(longest_label <= 10, 1, 0.8)

  if (is.null(xlim)) { # if not already given
    # 10% lower than min x
    xmin <- min(pca$x[, dims[1]])*1.1
    # Xmax is based on length of the labels when we are adding labels
    # (when there are less than 50 points)
    if (ncol(data) <= 50) {
      offset <- 0.6 * ((longest_label / 10) * label_cex)
    } else {
      offset <- 0.05
    }
    xmax <- max(pca$x[, dims[1]]) * (1 + offset)
    xlim <- c(xmin, xmax)
  }

  ## ylim
  if (is.null(ylim)) { # if not already given
    ymin <- min(pca$x[,dims[2]])*1.1
    ymax <- max(pca$x[,dims[2]])*1.1
    ylim <- c(ymin, ymax)
  }
  # make plot
  plot(NULL,
       xlim = xlim, ylim = ylim,
       xlab = paste0("PC ", dims[1], " (", round(pca$summary["Proportion of Variance", dims[1]] * 100, 2), " %)"),
       ylab = paste0("PC ", dims[2], " (", round(pca$summary["Proportion of Variance", dims[2]] * 100, 2), " %)"))
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.2)
  graphics::grid(col = "grey65") # add the grid UNDER the points
  graphics::points(
    x = pca$x[, dims[1]], y = pca$x[, dims[2]],
    pch = 21, bg = colorGroup2(groups)[groups], col = colorGroup2(groups)[groups],
    cex = cex.dot, lwd = 1, las = 1, cex.axis = 1.2, cex.lab = 1.3
  )

  # When points less than 50, add labels
  if (ncol(data) <= 50) {
    if (longest_label < 3) {
      adjust <- -0.6
    } else if (longest_label < 5) {
      adjust <- -0.4
    } else {
      adjust <- -0.2
    }
    graphics::text(
      x = pca$x[, dims[1]], y = pca$x[, dims[2]],
      labels = sampleLabels, cex = label_cex,
      adj = adjust,
      col = colorGroup2(groups)[as.character(groups)]
    )
  }

  if (legend) {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
           bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
           col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1)
  }

  data2 <- list(pca = pca, pca.data = data)

  invisible(data2)
}


#' Dendrogram for QC report
#'
#' Performs and then plots a hierarchical clustering analysis of sample intensities
#' across samples. Sample labels are colored according to the "groups" argument,
#' with colors determined by the \code{\link{colorGroup2}} function. By default,
#' uses only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_pca_plot
#' @param clust.metric The cluster metric used to define distance. Default is
#'   "euclidean". See
#'   \code{\link[ClassDiscovery:distanceMatrix]{ClassDiscovery::distanceMatrix}}
#'   for options.
#' @param clust.method The agglomeration method to use for clustering. Default
#'   is "complete", See \code{\link[stats:hclust]{stats::hclust}} for options.
#' @param cex.names Character expansion factor for sample labels. Default is 1.
#'
#' @return Invisibly returns a list with seven elements: \enumerate{
#'   \item "hc"- The hclust object returned by
#'     \code{\link[stats:hclust]{stats::hclust}}.
#'   \item "dist"- The distance matrix returned by
#'     \code{\link[ClassDiscovery:distanceMatrix]{ClassDiscovery::distanceMatrix}}.
#'   \item "corr.data"- The data matrix that went into the correlation analysis.
#'   \item "clust.method"- The clust.method argument used.
#'   \item "clust.metric"- The clust.method argument used.
#'   \item "stdize"- The stdize argument used.
#'   \item "top"- The top argument used.
#' }
#' @export
#'
#' @examples
#' # No examples yet

qc_dendrogram <- function(data, groups = NULL, sampleLabels = NULL,
                           top = 500, stdize = TRUE,
                           clust.metric = "euclidean", clust.method = "complete",
                           cex.names = 1, title = NULL, legend = TRUE) {

  # Check arguments
  clust.metric <- rlang::arg_match(
    arg = clust.metric,
    values = c(
      "pearson", "sqrt pearson", "spearman", "absolute pearson",
      "uncentered correlation", "weird", "cosine", "euclidean",
      "maximum", "manhattan", "canberra", "binary", "minkowski"
    ),
    multiple = FALSE
  )

  clust.method <- rlang::arg_match(
    arg = clust.method,
    values = c(
      "ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
      "median", "centroid"
    ), multiple = FALSE
  )

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = groups)

  # Get top variable proteins
  o <- order(rowVars(as.matrix(data)), decreasing = TRUE)
  data <- data[o, ]
  top <- ifelse(nrow(data) >= top, top, nrow(data))
  data <- data[1:top, ]

  ## center/scale rows (proteins) mean=0;stdev=1
  if (stdize) {
    data <- t(scale(x = t(data), center = TRUE, scale = TRUE))
  }



  # Get original graphics pars
  orig_par <- graphics::par(no.readonly = TRUE)
  # Schedule them to be reset when we exit the function
  on.exit(graphics::par(orig_par), add = TRUE)

  b_mar <- left_margin(x = sampleLabels)
  r_mar <- right_margin(x = groups)
  if (!legend) r_mar <- 1



  # Set up graphics params
  graphics::par(xpd = T,
                mar = c(b_mar, 6, 5, r_mar),
                oma = c(0.5, 0, 0.5, r_mar/2 + 3))

  d <- ClassDiscovery::distanceMatrix(dataset = data, metric = clust.metric)
  hc <- stats::hclust(d, method = clust.method)
  ClassDiscovery::plotColoredClusters(hc,
                                      labs = sampleLabels, cols = colorGroup2(groups)[groups],
                                      lwd = 1.5, las = 2, cex.axis = 1.2,
                                      xlab = "", ylab = "",
                                      font = 1, cex = cex.names,
                                      line = -0.6)
  title(main = ifelse(is.null(title), "", title), font.main = 1, line = 2.5, cex = 1.5)
  if (legend) {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
           bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
           box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
           col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1
    )
  }

  data2 <- list(hc = hc, dist = d,
                corr.data = data,
                clust.method = clust.method, clust.metric = clust.metric,
                stdize = stdize, top = top)

  invisible(data2)
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

