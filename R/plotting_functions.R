
plotCorHM <- function(data,
                      groups, batch = NULL, sampleLabels = NULL,
                      dir = ".", save = FALSE) {
  groups <- make_factor(groups)
  if (is.null(batch)) {
    batch <- c(rep("1", ncol(data)))
  }
  batch <- make_factor(as.character(batch), prefix = NULL)
  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }

  if (length(groups) < 100) {
    width <- round(0.0871 * length(groups)^2 + 24.375 * length(groups) + 473.02, 0)
    fontsize <- 12
  }

  if (length(groups) >= 100) {
    width <- round(0.0035 * length(groups)^2 + 10.035 * length(groups) + 146.15, 0)
    width <- round(width + width * 0.3, 0)
    fontsize <- 10
  }

  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  if (save) {
    grDevices::png(
      filename = file.path(dir, "CorrHeatmap.png"), units = "px",
      width = width, height = width, pointsize = 15
    )
  }
  cor_mat <- stats::cor(data, use = "pairwise.complete.obs", method = "pearson")

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

  hm_corr <- ComplexHeatmap::Heatmap(cor_mat,
    name = "Pearson correlation", border = TRUE,
    col = circlize::colorRamp2(seq(min(cor_mat), 1, ((1 - min(cor_mat)) / 7)),
      colors = RColorBrewer::brewer.pal(8, "Blues"),
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
  graphics::par(mar = c(10, 10, 10, 10))
  ComplexHeatmap::draw(hm_corr, heatmap_legend_side = "top", ht_gap = grid::unit(2, "cm"))
  if (save) grDevices::dev.off()

  return(invisible(cor_mat))
}



plotBoxplot <- function(data,
                        groups = NULL, sampleLabels = NULL,
                        title = NULL, legend = TRUE,
                        dir = ".", save = FALSE) {
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)

  ## remove rows with an NA
  ## TODO: DO WE WANT TO DO THIS??
  data <- data[!apply(is.na(data), 1, any), ]

  ## plot margins
  x2 <- left_margin(x = sampleLabels)
  x3 <- right_margin(x = groups)

  width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)

  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  if (save) {
    grDevices::png(
      filename = file.path(dir, "BoxPlot.png"), units = "px",
      width = width, height = 750, pointsize = 15
    )
  }

  # Reorder group affiliations
  # Since groups is an ordered factor (see make_factor()),
  # this puts them in the order of their levels, then sorts by name.
  # Does not necessarily sort alphabetically.
  group_order <- sort.int(groups, index.return = T)$ix

  # Then, reorder the input data and groups
  data <- data[,group_order]
  sampleLabels <- sampleLabels[group_order]
  groups <- groups[group_order]

  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(x2, 6, 3, x3), graphics::par(oma = c(0.5, 0, 0.5, x3 / 2 + 3)))
  box <- graphics::boxplot(data,
    col = colorGroup2(groups)[groups], names = sampleLabels,
    notch = FALSE, horizontal = FALSE, outline = FALSE,
    las = 2, cex.axis = 1, cex.labs = 1, cex = 1
  )
  graphics::mtext(side = 2, text = "Normalized Intensity", font = 1, line = 3, cex = 1.2)
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.1)
  if (!legend) {
    x3 <- 1
  } else {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
      bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
      box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
      col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1
    )
  }
  if (save) grDevices::dev.off()

  return(invisible(box))
}

plotViolin <- function(data,
                       groups = NULL, sampleLabels = NULL,
                       title = NULL, legend = TRUE,
                       dir = ".", save = FALSE) {
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)

  ## remove rows with an NA
  data <- data[!apply(is.na(data), 1, any), ]

  ## convert data data.frame to a named list
  plotData <- list()
  for (k in 1:ncol(data)) {
    plotData[[k]] <- data[, k]
  }
  names(plotData) <- colnames(data)

  ## plot margins
  x2 <- left_margin(x = sampleLabels)
  x3 <- right_margin(x = groups)

  width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)

  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  if (save) {
    grDevices::png(
      filename = file.path(dir, "ViolinPlot.png"), units = "px",
      width = width, height = 750, pointsize = 15
    )
  }

  # Reorder group affiliations
  # Since groups is an ordered factor (see make_factor()),
  # this puts them in the order of their levels, then sorts by name.
  # Does not necessarily sort alphabetically.
  group_order <- sort.int(groups, index.return = T)$ix

  # Then, reorder the input data and groups
  plotData <- plotData[group_order]
  sampleLabels <- sampleLabels[group_order]
  groups <- groups[group_order]

  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(x2, 6, 3, x3), graphics::par(oma = c(0.5, 0, 0.5, x3 / 2 + 3)))
  vio <- vioplot::vioplot(plotData,
    col = colorGroup2(groups)[groups], names = sampleLabels,
    main = "", ylab = "", xlab = "",
    las = 2, font = 1, cex.axis = 1, cex.labs = 1, cex = 1
  )
  graphics::mtext(side = 2, text = "Density", line = 3, cex = 1.2)
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.1)
  if (!legend) {
    x3 <- 1
  } else {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
      bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
      box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
      col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1
    )
  }
  if (save) grDevices::dev.off()

  return(invisible(vio))
} ## VIOLIN



## data=normalized intensities, or zscore intensities. if stdize=TRUE then it scales each row
## (each protein) to have a mean of zero and standard deviation = 1 if data matrix is already
## scaled then set stdize=FALSE, dims= vector of PC to plot.
plotPCA <- function(data, groups = NULL, sampleLabels = NULL, title = NULL, top = 500, stdize = TRUE,
                    dims = c(1, 2), cex.dot = 2, xlim = NULL, ylim = NULL, legend = TRUE, dir = ".", save = FALSE) {
  if (is.null(sampleLabels)) {
    sampleLabels <- as.character(colnames(data))
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = as.character(groups), prefix = NULL)


  ## TOP VARIABLE PROTEINS
  ## remove rows with an NA and get top most variable proteins
  data <- data[!apply(is.na(data), 1, any), ]
  o <- order(matrixStats::rowVars(as.matrix(data)), decreasing = TRUE)
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
  ## % variance explained by each PC
  # eigs <- pca$sdev^2; eigs[1] / sum(eigs)
  # summary(pca)

  ## xlim
  max.char <- max(nchar(sampleLabels))
  cex.char <- ifelse(max.char <= 10, 1, 0.8)

  if (is.null(xlim)) {
    # max.char <- max(nchar(sampleLabels));max.char
    # cex.char <- ifelse(max.char <= 10,1,0.9);cex.char
    xmin <- abs(min(pca$x[, dims[1]])) + 0.10 * abs(min(pca$x[, dims[1]]))

    # xmax <- max(pca$x[,dims[1]]) + (0.05 * (chr * cex.char)) * max(pca$x[,dims[1]]);print(xmax)
    offset <- 0.6 * ((max.char / 10) * cex.char)
    if (ncol(data) > 50) {
      offset <- 0.05
    }
    xmax <- max(pca$x[, dims[1]]) + max(pca$x[, dims[1]]) * offset
    xlim <- c(-xmin, xmax)
  }

  ## ylim
  if (is.null(ylim)) {
    ymin <- abs(min(range(pca$x[, dims[2]]))) + 0.10 * abs(min(range(pca$x[, dims[2]])))
    ymax <- max(range(pca$x[, dims[2]])) + 0.10 * max(range(pca$x[, dims[2]]))
    ylim <- c(-ymin, ymax)
  }

  ## plot margins
  x2 <- left_margin(x = sampleLabels)
  x3 <- right_margin(x = groups)

  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  if (save) {
    grDevices::png(
      filename = file.path(dir, "PCAplot.png"), units = "px",
      width = 750, height = 650, pointsize = 15
    )
  }

  op <- graphics::par(no.readonly = TRUE)
  graphics::par(mar = c(x2, 6, 3, x3), graphics::par(oma = c(0.5, 0, 0.5, x3 / 2 + 3)))
  plot(
    x = pca$x[, dims[1]], y = pca$x[, dims[2]],
    pch = 21, bg = colorGroup2(groups)[groups], col = colorGroup2(groups)[groups],
    cex = cex.dot, lwd = 1, las = 1, cex.axis = 1.2, cex.lab = 1.3,
    xlim = xlim, ylim = ylim,
    xlab = paste0("PC ", dims[1], " (", round(summary(pca)$importance["Proportion of Variance", dims[1]] * 100, 2), " %)"),
    ylab = paste0("PC ", dims[2], " (", round(summary(pca)$importance["Proportion of Variance", dims[2]] * 100, 2), " %)")
  )
  title(main = ifelse(is.null(title), "", title), font.main = 1, cex.main = 1.3, line = 1.2)
  graphics::grid()
  if (ncol(data) <= 50) {
    graphics::text(
      labels = sampleLabels, x = pca$x[, dims[1]], y = pca$x[, dims[2]],
      cex = cex.char, adj = ifelse(max.char < 3, -0.6, ifelse(max.char < 5, -0.4, -0.2)),
      col = colorGroup2(groups)[as.character(groups)]
    )
  }

  if (!legend) {
    x3 <- 1
  } else {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
      bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
      box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
      col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1
    )
  }
  if (save) grDevices::dev.off()

  data2 <- list(pca = pca, dat = data)
  return(invisible(data2))
}


## normalized intensities input (cols=samples, rows=proteins.
## NAs removed, then top variable proteins identified.
## if stdize is true the rows of data are centered and scaled to 0 and 1 z-score,
## distance calc. then clustering.
plotDendrogram <- function(data, groups = NULL, sampleLabels = NULL,
                           top = 500, stdize = TRUE,
                           clust.metric = "euclidean", clust.meth = "complete",
                           cex.names = 1, xlim = NULL, title = NULL, legend = TRUE,
                           dir = ".", save = FALSE) {
  clust.metric <- match.arg(
    arg = clust.metric,
    choices = c(
      "pearson", "sqrt pearson", "spearman", "absolute pearson",
      "uncentered correlation", "weird", "cosine", "euclidean",
      "maximum", "manhattan", "canberra", "binary", "minkowski"
    ),
    several.ok = FALSE
  )

  clust.meth <- match.arg(
    arg = clust.meth,
    choices = c(
      "ward.D", "ward.D2", "single", "complete", "average", "mcquitty",
      "median", "centroid"
    ), several.ok = FALSE
  )

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }
  if (is.null(groups)) {
    groups <- c(rep("group", ncol(data)))
  }
  groups <- make_factor(x = groups)

  ## remove rows with an NA,get top variable proteins
  data <- data[!apply(is.na(data), 1, any), ]
  o <- order(matrixStats::rowVars(data), decreasing = TRUE)
  data <- data[o, ]
  top <- ifelse(nrow(data) >= top, top, nrow(data))
  data <- data[1:top, ]

  ## center and scale rows (protein) (mean=0, std=1)
  if (stdize) {
    data <- t(scale(t(data)))
  }

  x2 <- left_margin(x = sampleLabels)
  x3 <- right_margin(x = groups)

  width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)


  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
  }
  if (save) {
    grDevices::png(
      filename = file.path(dir, "Dendrogram.png"), units = "px",
      width = width, height = 650, pointsize = 15
    )
  }


  op <- graphics::par(no.readonly = TRUE)
  graphics::par(xpd = T, mar = graphics::par()$mar + c(3, 0, 2, 5))
  d <- ClassDiscovery::distanceMatrix(dataset = data, metric = clust.metric)
  hc <- stats::hclust(d, method = clust.meth)
  ClassDiscovery::plotColoredClusters(hc,
    labs = sampleLabels, cols = colorGroup2(groups)[groups],
    lwd = 1.5, las = 2, cex.axis = 1.2, xlab = "", ylab = "", font = 1, cex = cex.names,
    line = -0.6
  )
  title(main = ifelse(is.null(title), "", title), font.main = 1, line = 2.5, cex = 1.5)
  if (legend) {
    legend(graphics::par("usr")[2], graphics::par("usr")[4],
      bty = "n", xpd = NA, legend = levels(groups), pch = 22, cex = 1,
      box.col = "transparent", box.lwd = 1, bg = "transparent", pt.bg = colorGroup2(groups),
      col = colorGroup2(groups), pt.cex = 1.2, pt.lwd = 0, inset = 0.02, horiz = F, ncol = 1
    )
  }
  graphics::par(mar = c(5, 4, 4, 2) + 0.1)

  if (save) grDevices::dev.off()

  data2 <- list(hc = hc, d = d, dat = data, clust.meth = clust.meth, clust.metric = clust.metric, stdize = stdize, top = top)

  return(invisible(data2))
}
