make_qc_report <- function(normList,
                           norm.meth = "cycloess",
                           groups = NULL,
                           batch = NULL,
                           sampleLabels = NULL,
                           stdize = TRUE,
                           top = 500,
                           dims = c(1, 2),
                           cex.dot = 2,
                           clust.metric = "euclidean",
                           clust.meth = "complete",
                           cex.names = 1,
                           xlim = NULL,
                           ylim = NULL,
                           legend = TRUE,
                           enrich = c("protein", "phospho"),
                           dir = NULL,
                           file = NULL,
                           save = FALSE,
                           keep.png = FALSE) {

  param <- stats <- list()

  norm.meth <- match.arg(
    arg = norm.meth,
    choices = unique(c(
      names(normList), "log2", "median", "mean", "vsn", "quantile",
      "cycloess", "rlr", "gi"
    )), several.ok = FALSE
  )


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


  enrich <- match.arg(enrich, choices = c("protein", "phospho"), several.ok = FALSE)



  ## NORMALIZED INTENSITY DATA
  ## list object or dataframe/matrix of norm. intensities
  if (any(class(normList) == "list")) {
    data <- normList[[norm.meth]]
  }
  if (any(class(normList) %in% "data.frame")) {
    data <- as.matrix(normList)
  }
  if (any(class(normList) %in% "matrix")) {
    data <- as.matrix(normList)
  }

  top <- ifelse(top > nrow(data), nrow(data), top)


  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }

  if (is.null(groups)) {
    groups <- rep("group", ncol(data))
  }
  groups <- make_factor(x = as.character(groups))

  if (!is.null(batch)) {
    batch <- make_factor(as.character(batch))
  }

  if (!save) {
    keep.png <- FALSE
  }


  ## CREATE OUTPUT DIRECTORIES, QC REPORT FILE NAME
  if (save) {

    ## CREATE QC OUTPUT DIRECTORY
    ## if use enrich type to create QC output directory
    if (is.null(dir)) {
      if (enrich == "protein") {
        dir <- file.path("protein_analysis", "01_quality_control")
      }
      if (enrich == "phospho") {
        dir <- file.path("phospho_analysis", "01_quality_control")
      }
    }
    if (!is.null(dir)) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
      }
    }

    ## FILENAME NOT NULL
    if (!is.null(file)) {
      if (file_ext(file) != "pdf") {
        stop(
          "\nError! Invalid output file type...\nincorrect: file = '", file, "'",
          "\ncorrect:   file = 'QC_Report.pdf'"
        )
      }
      pngdir <- gsub(".pdf", "", file)
      if (file.exists(file.path(dir, file))) {
        file <- make_new_filename(x = file, dir = dir)
        no <- sub(".pdf", "", sub(".*_", "", file))
        pngdir <- gsub(".pdf", "", file)
      }
      if (!file.exists(file.path(dir, file))) {
        pngdir <- gsub(".pdf", "", file)
      }
    }

    if (is.null(file)) {
      file <- "QC_Report.pdf"
      no <- ""
      if (file.exists(file.path(dir, file))) {
        file <- make_new_filename(x = file, dir = dir)
        no <- sub(".pdf", "", sub(".*_", "", file))
        pngdir <- gsub(".pdf", "", file)
      }
      if (!file.exists(file.path(dir, file))) {
        pngdir <- gsub(".pdf", "", file)
      }
    }

    print(paste("QC output directory:", dir))
    print(paste("QC report file: ", file))
    print(paste("png directory: ", pngdir))
  } ## SAVE == TRUE


  ## -----------------
  ##  CREATE PLOTS
  ## -----------------
  ## BOX PLOTS
  box <- plotBoxplot(data = data,
                     groups = groups,
                     sampleLabels = sampleLabels,
                     title = "Box Plot",
                     legend = legend,
                     dir = dir,
                     save = save)
  if (!is.null(batch)) {
    if (save) {
      width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
      png(filename = file.path(dir, "BoxPlot2.png"),
          units = "px",
          width = width,
          height = 750,
          pointsize = 15)
    }
    box2 <- plotBoxplot(data = data,
                        groups = batch,
                        sampleLabels = sampleLabels,
                        title = "Box Plot",
                        legend = legend,
                        dir = dir,
                        save = FALSE)
    if (save == TRUE) {
      dev.off()
    }
  }

  ## VIOLIN PLOTS
  vio <- plotViolin(data,
                    groups = groups,
                    sampleLabels,
                    title = "Violin Plot",
                    legend,
                    dir,
                    save)
  if (!is.null(batch)) {
    width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
    if (save) {
      png(filename = file.path(dir, "violinPlot2.png"),
          units = "px",
          width = width,
          height = 750,
          pointsize = 15)
    }
    vio2 <- plotViolin(data = data,
                       groups = batch,
                       sampleLabels = sampleLabels,
                       title = "Violin Plot",
                       legend = legend,
                       dir = dir,
                       save = FALSE)
    if (save) {
      dev.off()
    }
  }

  ## PCA PLOTS
  pca <- plotPCA(
    data = data,
    groups = groups,
    sampleLabels = sampleLabels,
    title = "PCA Plot",
    top = top,
    stdize = stdize,
    dims = dims,
    cex.dot = cex.dot,
    xlim = xlim,
    ylim = ylim,
    legend = legend,
    dir = dir,
    save = save
  )
  if (!is.null(batch)) {
    if (save == TRUE) {
      png(filename = file.path(dir, "PCAplot2.png"), units = "px", width = 750, height = 650, pointsize = 15)
    }
    pca2 <- plotPCA(
      data = data,
      groups = batch,
      sampleLabels = sampleLabels,
      title = "PCA Plot",
      top = top,
      stdize = stdize,
      dims = dims,
      cex.dot = cex.dot,
      xlim = xlim,
      ylim = xlim,
      legend = legend,
      dir = dir,
      save = FALSE
    )
    if (save) {
      dev.off()
    }
  }

  ## CLUSTER DENDROGRAMS
  dendro <- plotDendrogram(
    data = data, groups = groups, sampleLabels = sampleLabels, top = top, stdize = stdize,
    clust.metric = clust.metric, clust.meth = clust.meth,
    cex.names = 1, xlim = NULL, title = "Cluster Dendrogram", legend = legend, dir = dir, save = save
  )
  if (!is.null(batch)) {
    if (save) {
      width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
      print(width)
      png(filename = file.path(dir, "Dendrogram2.png"), units = "px", width = width, height = 650, pointsize = 15)
    }
    dendro2 <- plotDendrogram(
      data = data, groups = batch, sampleLabels = sampleLabels, top = top, stdize = stdize,
      clust.metric = clust.metric, clust.meth = clust.meth,
      cex.names = 1, xlim = NULL, title = "Cluster Dendrogram", legend = legend, dir = dir, save = FALSE
    )
    if (save) {
      dev.off()
    }
  }

  ## SAMPLE CORRELATION HEATMAP (PEARSON)
  corhm <- plotCorHM(data = data,
                     groups = groups, batch = batch, sampleLabels = sampleLabels,
                     dir = dir, save = save)


  ## PLOTS LIST OBJECT
  if (is.null(batch)) {
    plotList <- list(
      box = box, vio = vio, pca = pca, dendro = dendro, corhm = corhm, norm.meth = norm.meth,
      dir = ifelse(save, dir, NULL),
      file = ifelse(save, file, NULL)
    )
  }
  if (!is.null(batch)) {
    plotList <- list(
      box = box, box2 = box2, vio = vio, vio2 = vio2, pca = pca, pca2 = pca2,
      dendro = dendro, dendro2 = dendro2, corhm = corhm, norm.meth = norm.meth,
      dir = ifelse(save, dir, NULL),
      file = ifelse(save, file, NULL)
    )
  }


  ## -----------------------------
  ##  SAVE QC REPORT PDF FILE
  ## -----------------------------
  if (save) {

    ##  MAKE PDF FILE
    pdf(file.path(dir, file), paper = "USr", pagecentre = TRUE, pointsize = 15, width = 12, height = 8)

    ## QC_REPORT.PDF (NO BATCH)
    if (is.null(batch)) {
      files <- c("BoxPlot.png", "ViolinPlot.png", "PCAplot.png", "Dendrogram.png", "CorrHeatmap.png")
      pnglist <- paste0(paste0(file.path(dir), "/"), files)
      pnglist
      thePlots <- lapply(1:length(pnglist), function(i) {
        grid::rasterGrob(png::readPNG(pnglist[i], native = F))
      })

      if (ncol(data) <= 50) {
        do.call(gridExtra::grid.arrange, c(thePlots[1:4], ncol = 2))
      }
      do.call(gridExtra::grid.arrange, c(thePlots[1], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[2], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[3], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[4], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[5], ncol = 1))
    } ## BATCH IS NULL


    ## QC_REPORT.PDF (BATCH INCLUDED)
    if (!is.null(batch)) {

      ## QC_REPORT.PDF
      files <- c(
        "BoxPlot.png", "BoxPlot2.png", "ViolinPlot.png", "ViolinPlot2.png",
        "PCAplot.png", "PCAplot2.png", "Dendrogram.png", "Dendrogram2.png", "CorrHeatmap.png"
      )
      pnglist <- paste0(paste0(file.path(dir), "/"), files)
      pnglist
      thePlots <- lapply(1:length(pnglist), function(i) {
        grid::rasterGrob(png::readPNG(pnglist[i], native = F))
      })

      ## QC REPORT PDF FILE
      if (ncol(data) <= 50) {
        do.call(gridExtra::grid.arrange, c(thePlots[1:4], ncol = 2))
        do.call(gridExtra::grid.arrange, c(thePlots[5:8], ncol = 2))
      }
      do.call(gridExtra::grid.arrange, c(thePlots[1], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[2], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[3], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[4], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[5], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[6], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[7], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[8], ncol = 1))
      do.call(gridExtra::grid.arrange, c(thePlots[9], ncol = 1))
    } ## BATCH IS NOT NULL

    dev.off()

    ## REMOVE PNG FILES
    if (!keep.png) {
      unlink(pnglist)
      print(file.path(dir, file))
      print("png files removed...")
    } else {
      if (!dir.exists(file.path(dir, pngdir))) {
        dir.create(file.path(dir, pngdir), recursive = TRUE)
      }
      lapply(files, function(x) {
        file.copy(from = file.path(dir, x), to = file.path(dir, pngdir, x))
        file.remove(file.path(dir, x))
      })
      print(paste("png files moved to :", file.path(dir, pngdir)))
    }
  }


  param[["norm.meth"]] <- norm.meth
  param[["batch"]] <- ifelse(is.null(batch), "NULL", paste(unique(batch), collapse = ", "))
  param[["stdize"]] <- stdize
  param[["top"]] <- top
  param[["dims"]] <- paste(dims, collapse = ", ")
  param[["clust.metric"]] <- clust.metric
  param[["clust.meth"]] <- clust.meth
  if (save) {
    param[["dir"]] <- dir
  }
  if (save) {
    param[["file"]] <- file
  }
  if (keep.png) {
    param[["png.dir"]] <- pngdir
  }
  param[["enrich"]] <- enrich

  stats[["no_samples"]] <- ncol(data)
  stats[["no_groups"]] <- length(unique(groups))
  stats[["no_batches"]] <- ifelse(is.null(batch), 0, length(unique(batch)))
  stats[["tot_no_rows"]] <- nrow(data)
  stats[["top_no_rows"]] <- top

  logs <- make_log(param, stats, title = "QC REPORT", save = TRUE)


  data2 <- list(plots = plotList, param = logs$param, stats = logs$stats)

  print("QC report files created. Success!!")

  return(invisible(data2))
}
