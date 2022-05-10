#' Create quality control report
#'
#' Creates, and optionally saves as a PDF report, a bunch of plots which give
#' information on the distribution, clustering, and correlation of protein intensities
#' across samples for the chosen normalization method. By default, the PCA and
#' clustering analyses are performed on the top 500 most variable proteins. This
#' function is a wrapper that calls many subfunctions: \itemize{
#'   \item Makes boxplots of per-sample intensities
#'     with \code{\link{qc_boxplot}}.
#'   \item Makes violin plots of per-sample intensities
#'     with \code{\link{qc_violin_plot}}.
#'   \item Performs and plots a PCA with \code{\link{qc_pca_plot}}.
#'   \item Does hierarchical clustering and plots a dendrogram with
#'     \code{\link{qc_dendrogram}}.
#'   \item Plots a correlation heatmap with \code{\link{qc_corr_hm}}.
#' } See the documentation of these subfunctions for more info.
#'
#'
#' @inheritParams make_proteinorm_report
#' @param norm.method The normalization method for which to perform the QC analysis.
#'   Should be a name of one of the datasets in normList.
#' @param file The file name of the report to be saved. Must end in .pdf. Will
#'   default to "QC_Report.pdf" if no filename is provided.
#' @param legend Include legends in the plots within the report? Default is TRUE.
#' @param top.proteins The number of most variable proteins to use for the analysis.
#'   Default is 500.
#' @param stdize Should input data be standardized to a mean of 0 and std.dev of
#'   1? If input data are not yet standardized, should be TRUE. Default is TRUE.
#' @param pca.axes  A numeric vector of length 2 which lists the PC axes to plot.
#'   Default is c(1,2), to plot the first two principal components.
#' @param pca.xlim Optional. Custom x-axis limits of the PCA plot. By default,
#'   the min is 10% below the min PC score on the x-axis, and the max is above
#'   the max PC score on the x-axis, with some padding for sample labels.
#' @param pca.ylim Optional. Custom y-axis limits of the PCA plot. Default is
#'   10% above and below the min/max PC score on the y-axis.
#' @param pca.dot Character expansion factor for the size of the points in the
#'   PCA plot. Default is 2.
#' @inheritParams qc_dendrogram
#' @param clust.label.cex Character expansion factor for sample labels in the
#'   dendrogram. Default is 1.
#'
#' @return Invisibly returns a list with three slots: \enumerate{
#'   \item "plots"- A large list, where each element in the list is the returned
#'     object from the corresponding plotting function for that type of plot.
#'   \item "stats"- A dataframe with statistics on the data.
#'   \item "param"- A dataframe giving the parameters used for making the report.
#' }
#'
#' @export
#'
#' @examples
#' # No examples yet
make_qc_report <- function(normList,
                           groups = NULL, batch = NULL, sampleLabels = NULL,
                           norm.method = NULL, enrich = c("protein", "phospho"),
                           dir = NULL, file = NULL, save = TRUE,
                           overwrite = FALSE, keep.png = FALSE,
                           legend = TRUE,
                           top.proteins = 500, # Number of top variable proteins to include in PCA and dendrogram
                           stdize = TRUE,
                           pca.axes = c(1, 2),
                           pca.xlim = NULL, pca.ylim = NULL, pca.dot = 2,
                           clust.metric = "euclidean", clust.method = "complete",
                           clust.label.cex = 1
                           ) {

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
    validate_filename(file, allowed_exts = c("pdf"))

    # Check if filename already exists, remake it if so
    if (file.exists(file.path(out_dir, file))) {
      if (overwrite) {
        cli::cli_inform("{.path {file}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
      } else {
        cli::cli_abort(c("{.path {file}} already exists in {.path {out_dir}}",
                     "!" = "and {.arg overwrite} == {.val {overwrite}}",
                     "i" = "Give {.arg file} a unique name or set {.arg overwrite} to {.val TRUE}"))
      }
    }

    # Notify user of where tmp png files will be
    cli::cli_inform("Temporarily saving {.path .png} files in {.path {out_dir}}")
  }

  ## NORMALIZED INTENSITY DATA
  # TODO: figure out what exactly is going on here.
  # Isn't the data always a list?
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

  ##################
  ## Create plots ##
  ##################

  # Set up widths
  # Will make individual for now, so we can fine tune if needed
  if (save) {
    width_box <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
    width_vio <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
    width_pca <- 750
    width_dendro <- round(0.0191 * length(groups)^2 + 12.082 * length(groups) + 671.75, 0)
    if (length(groups) < 100) {
      width_cor <- round(0.0871 * length(groups)^2 + 24.375 * length(groups) + 473.02, 0)
    } else {
      width_cor <- round((0.0035 * length(groups)^2 + 10.035 * length(groups) + 146.15)*1.3, 0)
    }
  }

  ####
  ## BOXPLOT(S)
  ####

  if (save) {
    # Open Graphics device for main plot
    grDevices::png(
      filename = file.path(out_dir, "BoxPlot.png"), units = "px",
      width = width_box, height = 750, pointsize = 15
    )
  }
  # Group boxplot
  box <- qc_boxplot(data = data,
                    groups = groups,
                    sampleLabels = sampleLabels,
                    title = "Grouped by group",
                    legend = legend)
  if (save) grDevices::dev.off()

  # If plotting by batch
  if (!is.null(batch)) {
    if (save) {
      grDevices::png(
        filename = file.path(out_dir, "BoxPlot2.png"), units = "px",
        width = width_box, height = 750, pointsize = 15
      )
    }
    box2 <- qc_boxplot(data = data,
                       groups = batch,
                       sampleLabels = sampleLabels,
                       title = "Grouped by batch",
                       legend = legend)
    if (save) grDevices::dev.off()
  }

  ####
  ## VIOLIN PLOT(S)
  ####
  if (save) {
    grDevices::png(
      filename = file.path(out_dir, "ViolinPlot.png"), units = "px",
      width = width_vio, height = 750, pointsize = 15)
  }
  # Main group plot
  vio <- qc_violin_plot(data,
                        groups = groups,
                        sampleLabels,
                        title = "Grouped by group",
                        legend)
  if (save) grDevices::dev.off()

  if (!is.null(batch)) {
    if (save) {
      grDevices::png(filename = file.path(out_dir, "ViolinPlot2.png"),
          units = "px", width = width_vio, height = 750, pointsize = 15)
    }
    vio2 <- qc_violin_plot(data = data,
                           groups = batch,
                           sampleLabels = sampleLabels,
                           title = "Grouped by batch",
                           legend = legend)
    if (save) grDevices::dev.off()
  }


  ####
  ## PCA PLOT(S)
  ####
  if (save) {
    grDevices::png(filename = file.path(out_dir, "PCAplot.png"), units = "px",
      width = width_pca, height = 650, pointsize = 15)
  }
  pca <- qc_pca_plot(
    data = data, groups = groups, sampleLabels = sampleLabels,
    title = "PCA, colored by group",
    top = top.proteins, stdize = stdize, dims = pca.axes,
    cex.dot = pca.dot, xlim = pca.xlim, ylim = pca.ylim,
    legend = legend)
  if (save) grDevices::dev.off()


  if (!is.null(batch)) {
    if (save) {
      grDevices::png(filename = file.path(out_dir, "PCAplot2.png"), units = "px",
                     width = width_pca, height = 650, pointsize = 15)
    }
    pca2 <- qc_pca_plot(
      data = data, groups = batch, sampleLabels = sampleLabels,
      title = "PCA, colored by batch",
      top = top.proteins, stdize = stdize, dims = pca.axes,
      cex.dot = pca.dot, xlim = pca.xlim, ylim = pca.ylim,
      legend = legend)
    if (save) grDevices::dev.off()
  }

  ####
  ## DENDROGRAM(S)
  ####
  if (save) {
    grDevices::png(filename = file.path(out_dir, "Dendrogram.png"), units = "px",
      width = width_dendro, height = 650, pointsize = 15)
  }
  dendro <- qc_dendrogram(
    data = data, groups = groups, sampleLabels = sampleLabels,
    top = top.proteins, stdize = stdize,
    clust.metric = clust.metric, clust.method = clust.method,
    cex.names = clust.label.cex,
    title = paste0("Method: ", clust.method, ". Metric: ", clust.metric, "\nColored by group"),
    legend = legend
  )
  if (save) grDevices::dev.off()
  if (!is.null(batch)) {
    if (save) {
      grDevices::png(filename = file.path(out_dir, "Dendrogram2.png"), units = "px",
                     width = width_dendro, height = 650, pointsize = 15)
    }
    dendro2 <- qc_dendrogram(
      data = data, groups = batch, sampleLabels = sampleLabels,
      top = top.proteins, stdize = stdize,
      clust.metric = clust.metric, clust.method = clust.method,
      cex.names = clust.label.cex,
      title = paste0("Method: ", clust.method, ". Metric: ", clust.metric, "\nColored by batch"),
      legend = legend
      )
    if (save) grDevices::dev.off()
  }

  ####
  ## HEATMAP(S)
  ####
  if (save) {
    grDevices::png(filename = file.path(out_dir, "CorrHeatmap.png"), units = "px",
      width = width_cor, height = width_cor, pointsize = 15)
  }
  corhm <- qc_corr_hm(data = data,
                     groups = groups, batch = batch, sampleLabels = sampleLabels)
  if (save) grDevices::dev.off()


  ##################
  ## Save report  ##
  ##################
  if (save) {

    cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")

    ##  MAKE PDF FILE
    grDevices::pdf(file.path(out_dir, file), paper = "USr", pagecentre = TRUE, pointsize = 15, width = 12, height = 8)

    if (is.null(batch)) {
      files <- c("BoxPlot.png", "ViolinPlot.png", "PCAplot.png", "Dendrogram.png", "CorrHeatmap.png")
    } else {
      files <- c("BoxPlot.png", "BoxPlot2.png", "ViolinPlot.png",
                 "ViolinPlot2.png", "PCAplot.png", "PCAplot2.png",
                 "Dendrogram.png", "Dendrogram2.png", "CorrHeatmap.png")
    }

    # Collect pngs as plot objects
    pnglist <- paste0(paste0(file.path(out_dir), "/"), files)
    thePlots <- lapply(1:length(pnglist), function(i) {
      grid::rasterGrob(png::readPNG(pnglist[i], native = F))
    })

    if (ncol(data) <= 50) { # When project is small enough
        # Put first 4 plots onto the first page
        gridExtra::grid.arrange(grobs = thePlots[1:4], ncol = 2)
        # And plot the next 4 plots on another page when batch
        if (!is.null(batch)) {
          gridExtra::grid.arrange(grobs = thePlots[5:8], ncol = 2)
        }
    }

    # Print all plots on their own page
    lapply(X = thePlots, FUN = gridExtra::grid.arrange, ncol = 1)
    grDevices::dev.off()

    # Deal with the .png files
    if (!keep.png) { # Delete
      unlink(pnglist)
      cli::cli_inform("Temporary {.path .png} files removed from {.path {out_dir}} ")
    } else { # Move to a new png directory
      # Setup dir name
      pngdir <- gsub(".pdf", "_pngs", file)
      # Create if it doesn't exist
      if (!dir.exists(file.path(out_dir, pngdir))) {
        dir.create(file.path(out_dir, pngdir), recursive = TRUE)
      }
      # Move files
      lapply(pnglist, function(x) {
        file.copy(from = x, to = file.path(out_dir, pngdir, basename(x)))
        file.remove(x)
      })
    }
  }

  ##############
  ## Logging  ##
  ##############
  param <- stats <- list()
  param[["norm.method"]] <- norm.method
  param[["batch"]] <- ifelse(is.null(batch), "NULL", paste(unique(batch), collapse = ", "))
  param[["stdize"]] <- stdize
  param[["top.proteins"]] <- top.proteins
  param[["dims"]] <- paste(pca.axes, collapse = ", ")
  param[["clust.metric"]] <- clust.metric
  param[["clust.method"]] <- clust.method
  param[["dir"]] <- ifelse(save, out_dir, NA)
  param[["file"]] <- ifelse(save, file, NA)
  param[["png.dir"]] <- ifelse(keep.png, pngdir, NA)
  param[["enrich"]] <- enrich

  stats[["num_samples"]] <- ncol(data)
  stats[["num_groups"]] <- length(unique(groups))
  stats[["num_batches"]] <- ifelse(is.null(batch), 0, length(unique(batch)))
  stats[["total_num_rows"]] <- nrow(data)
  stats[["top_num_rows"]] <- top.proteins

  logs <- make_log(param, stats, title = "QC REPORT", save = TRUE)


  ###############
  ## Returning ##
  ###############
  if (is.null(batch)) {
    plotList <- list(
      box = box, violin = vio, pca = pca, dendro = dendro, corr_hm = corhm
    )
  } else {
    plotList <- list(
      box_group = box, box_batch = box2,
      violin_group = vio, violin_batch = vio2,
      pca_group = pca, pca_batch = pca2,
      dendro_group = dendro, dendro_batch = dendro2,
      corr_hm = corhm
    )
  }

  output_data <- list(plots = plotList, stats = logs$stats, param = logs$param)

  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(output_data)
}
