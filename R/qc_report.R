make_qc_report <- function(normList,
                           norm.method = "cycloess",
                           groups = NULL,
                           batch = NULL,
                           sampleLabels = NULL,
                           stdize = TRUE,
                           top = 500,
                           dims = c(1, 2),
                           cex.dot = 2,
                           clust.metric = "euclidean",
                           clust.method = "complete",
                           cex.names = 1,
                           xlim = NULL,
                           ylim = NULL,
                           legend = TRUE,
                           enrich = c("protein", "phospho"),
                           dir = NULL,
                           file = NULL,
                           save = FALSE,
                           keep.png = FALSE) {


  cli::cli_rule()

  #################
  ## Check args  ##
  #################
  # enrich: possible values encoded in function def
  enrich <- rlang::arg_match(enrich)

  # Here, possibilities are defined by use, and whatever
  # is available already in the normList.
  norm.method <- rlang::arg_match(
    arg = norm.method,
    values = unique(c(
      names(normList), "log2", "median", "mean", "vsn", "quantile",
      "cycloess", "rlr", "gi"
    )), multiple = FALSE
  )

  # here, the values are presumably from some distance function?
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



  ##################
  ## Set defaults ##
  ##################

  # TODO: see if this is needed?
  # If we're not saving, override keep.png
  if (!save) {
    keep.png <- FALSE
  }

  # Sort out some defaults if arguments are not supplied
  # Inform/alert user for groups, as this is semi-serious
  if (is.null(groups)) {
    groups <- rep("group", ncol(normList[[1]]))
    cli::cli_inform(cli::col_yellow("{.arg groups} argument is empty. Considering all samples/columns in {.arg normList} as one group."))
  }

  # Set default dir if not provided
  if (is.null(dir)) {
    out_dir <- file.path(paste0(enrich, "_analysis"), "01_quality_control")
    if (save) { # but only print alert if we're saving
      cli::cli_inform(cli::col_yellow("{.arg dir} argument is empty. Setting output directory to: {.path {out_dir}}"))
    }
  } else {
    out_dir <- dir
  }
  # Set default report name if not provided
  if (is.null(file)) {
    file <- "QC_Report.pdf"
    if (save) { # but only print alert if we're saving
      cli::cli_inform(cli::col_yellow("{.arg file} argument is empty. Saving report to: {.path {out_dir}/{file}}"))
    }
  }

  # Set defaults silently for sampleLabels, expected to just grab them
  # from the normList
  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(normList[[1]])
  }

  # This code doesn't seem to matter:
  # so far in testing, works the same when commented or uncommented?
  groups <- make_factor(x = as.character(groups))
  if (!is.null(batch)) {
    batch <- make_factor(as.character(batch))
  }


  ###########################
  ## If save, set up files ##
  ###########################
  if (save) {
    # Make directory if it doesn't exist
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = T)
    }

    # Check that the filename is a pdf
    if (tools::file_ext(file) != "pdf") {
      cli::cli_abort(c("Invalid filename for report",
                       "X" = "{.arg file} must end in {.path .pdf}, not {.path {tools::file_ext(file)}}"))
    }

    # Check if filename already exists, remake it if so
    if (file.exists(file.path(out_dir, file))) {
      new_file <- make_new_filename(x = file, dir = out_dir)
      cli::cli_inform("{.path {file}} already exists. Renaming as {.path {new_file}} instead.")
      file <- new_file
    }

    # Notify user of where tmp png files will be
    cli::cli_inform("Temporarily saving {.path .png} files in {.path {out_dir}}")
  }

  ## NORMALIZED INTENSITY DATA
  # TODO: figure out what exactly is going on here.
  # Isn't the data always a normList?
  ## list object or dataframe/matrix of norm. intensities
  if (any(class(normList) == "list")) {
    data <- normList[[norm.method]]
  }
  if (any(class(normList) %in% "data.frame")) {
    data <- as.matrix(normList)
  }
  if (any(class(normList) %in% "matrix")) {
    data <- as.matrix(normList)
  }


  # TODO: NEED TO FIGURE OUT MORE WHAT TOP DOES BEFORE I MESS WITH IT
  top <- ifelse(top > nrow(data), nrow(data), top)


  ## -----------------
  ##  CREATE PLOTS
  ## -----------------
  ## BOX PLOTS
  box <- plotBoxplot(data = data,
                     groups = groups,
                     sampleLabels = sampleLabels,
                     title = "Box Plot",
                     legend = legend,
                     dir = out_dir,
                     save = save)
  if (!is.null(batch)) {
    if (save) {
      width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
      grDevices::png(filename = file.path(out_dir, "BoxPlot2.png"),
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
                        dir = out_dir,
                        save = FALSE)
    if (save == TRUE) {
      grDevices::dev.off()
    }
  }

  ## VIOLIN PLOTS
  vio <- plotViolin(data,
                    groups = groups,
                    sampleLabels,
                    title = "Violin Plot",
                    legend,
                    out_dir,
                    save)
  if (!is.null(batch)) {
    width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
    if (save) {
      grDevices::png(filename = file.path(out_dir, "ViolinPlot2.png"),
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
                       dir = out_dir,
                       save = FALSE)
    if (save) {
      grDevices::dev.off()
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
    dir = out_dir,
    save = save
  )
  if (!is.null(batch)) {
    if (save == TRUE) {
      grDevices::png(filename = file.path(out_dir, "PCAplot2.png"), units = "px", width = 750, height = 650, pointsize = 15)
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
      dir = out_dir,
      save = FALSE
    )
    if (save) {
      grDevices::dev.off()
    }
  }

  ## CLUSTER DENDROGRAMS
  dendro <- plotDendrogram(
    data = data, groups = groups, sampleLabels = sampleLabels, top = top, stdize = stdize,
    clust.metric = clust.metric, clust.meth = clust.method,
    cex.names = 1, xlim = NULL, title = "Cluster Dendrogram", legend = legend, dir = out_dir, save = save
  )
  if (!is.null(batch)) {
    if (save) {
      width <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
      print(width)
      grDevices::png(filename = file.path(out_dir, "Dendrogram2.png"), units = "px", width = width, height = 650, pointsize = 15)
    }
    dendro2 <- plotDendrogram(
      data = data, groups = batch, sampleLabels = sampleLabels, top = top, stdize = stdize,
      clust.metric = clust.metric, clust.meth = clust.method,
      cex.names = 1, xlim = NULL, title = "Cluster Dendrogram", legend = legend, dir = out_dir, save = FALSE
    )
    if (save) {
      grDevices::dev.off()
    }
  }

  ## SAMPLE CORRELATION HEATMAP (PEARSON)
  corhm <- plotCorHM(data = data,
                     groups = groups, batch = batch, sampleLabels = sampleLabels,
                     dir = out_dir, save = save)


  ## PLOTS LIST OBJECT
  if (is.null(batch)) {
    plotList <- list(
      box = box, vio = vio, pca = pca, dendro = dendro, corhm = corhm, norm.meth = norm.method,
      dir = ifelse(save, out_dir, NULL),
      file = ifelse(save, file, NULL)
    )
  }
  if (!is.null(batch)) {
    plotList <- list(
      box = box, box2 = box2, vio = vio, vio2 = vio2, pca = pca, pca2 = pca2,
      dendro = dendro, dendro2 = dendro2, corhm = corhm, norm.meth = norm.method,
      dir = ifelse(save, out_dir, NULL),
      file = ifelse(save, file, NULL)
    )
  }


  ## -----------------------------
  ##  SAVE QC REPORT PDF FILE
  ## -----------------------------
  if (save) {

    ##  MAKE PDF FILE
    grDevices::pdf(file.path(out_dir, file), paper = "USr", pagecentre = TRUE, pointsize = 15, width = 12, height = 8)

    ## QC_REPORT.PDF (NO BATCH)
    if (is.null(batch)) {
      files <- c("BoxPlot.png", "ViolinPlot.png", "PCAplot.png", "Dendrogram.png", "CorrHeatmap.png")
      pnglist <- paste0(paste0(file.path(out_dir), "/"), files)
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
      pnglist <- paste0(paste0(file.path(out_dir), "/"), files)
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

    grDevices::dev.off()

    ## REMOVE PNG FILES
    if (!keep.png) {
      unlink(pnglist)
      print(file.path(out_dir, file))
      print("png files removed...")
    } else {
      if (!dir.exists(file.path(out_dir, pngdir))) {
        dir.create(file.path(out_dir, pngdir), recursive = TRUE)
      }
      lapply(files, function(x) {
        file.copy(from = file.path(out_dir, x), to = file.path(out_dir, pngdir, x))
        file.remove(file.path(out_dir, x))
      })
      print(paste("png files moved to :", file.path(out_dir, pngdir)))
    }
  }


  param <- stats <- list()
  param[["norm.method"]] <- norm.method
  param[["batch"]] <- ifelse(is.null(batch), "NULL", paste(unique(batch), collapse = ", "))
  param[["stdize"]] <- stdize
  param[["top"]] <- top
  param[["dims"]] <- paste(dims, collapse = ", ")
  param[["clust.metric"]] <- clust.metric
  param[["clust.method"]] <- clust.method
  if (save) {
    param[["dir"]] <- out_dir
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
