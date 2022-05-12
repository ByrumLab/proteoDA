#' Create proteinorm report
#'
#' Creates, and optionally saves as a PDF report, a bunch of plots which give
#' information about the performance of different normalization metrics and
#' amount of missing data. A wrapper that calls many subfunctions: \itemize{
#'   \item Makes boxplots
#'     of some normalization metrics with \code{\link{proteinormMetricBoxplot}}.
#'   \item Makes density plots of the log2ratio between groups with
#'     \code{\link{plotLogRatioDensity}}.
#'   \item Plots total sample intensity with \code{\link{plotTotInten}}.
#'   \item  Makes missing value heatmaps with \code{\link{plotHeatmapsForReport}}
#' } Each of these functions has further subfunctions, see their documentation.
#'
#'
#' @param normList A list of normalized data matrices. Generally, the "normList"
#'   slot of the list that is output by \code{\link{process_data}}.
#' @param groups Optional, a vector describing what experimental or treatment
#'   group each sample belongs to. Generally, the "group" column in the targets
#'   data frame output by \code{\link{process_data}}. If not supplied, will
#'   warn the user and treat all samples as belonging to the same group.
#' @param batch Optional, a vector describing what batch each sample belongs to.
#'   If not supplied, will treat all samples as being from the same batch.
#' @param sampleLabels Optional, a set of sample labels to use, in the same order as
#'   the columns in the normList. If not supplied, defaults to using the column
#'   names in the normList.
#' @param enrich What type of analysis is this for? Options are
#'   "protein" and "phospho". Only used for making output dir name when one isn't
#'   supplied.
#' @param save Should the report be saved as a PDF? Default is TRUE. If FALSE,
#'   will print all plots to the R or RStudio graphics device.
#' @param overwrite Should report file be overwritten if it already exists?
#'   Default is FALSE.
#' @param dir The directory in which to save the report. If not provided,
#'   will default to either "protein_analysis/01_quality_control" or
#'   "phospho_analysis/01_quality_control" within the current working
#'   directory, depending on the enrich argument.
#' @param file The file name of the report to be saved. Must end in .pdf. Will
#'   default to "proteiNorm_Report.pdf" if no filename is provided.
#' @param keep.png Keep the individual .png files of each plot? Default is FALSE.
#' @param showAllProteins For missing data heatmaps, show all proteins (including
#'  ones with no missing data)? Default is FALSE.
#'
#' @return Invisibly, returns a list with nine slots:
#'   \enumerate{
#'     \item "pcv"- A list, equal in length to the input normList, where each
#'       element of the list gives the results of the \code{\link{PCV}} function
#'       for a given normaliztion method.
#'     \item "pmad"- A list, equal in length to the input normList, where each
#'       element of the list gives the results of the \code{\link{PMAD}} function
#'       for a given normaliztion method.
#'     \item "pev"- A list, equal in length to the input normList, where each
#'       element of the list gives the results of the \code{\link{PEV}} function
#'       for a given normaliztion method.
#'     \item "cor"- A list, equal in length to the input normList, where each
#'       element of the list gives the results of the \code{\link{COR}} function
#'       for a given normaliztion method.
#'     \item "lograt"- A list, equal in length to the input normList, where each
#'       element of the list is the density object of the log2ratio density for a
#'       given normalization method. See \code{\link{plotLogRatioDensity}} and
#'       \code{\link{log2ratio}}.
#'     \item "nahm"- A list of length 4, where the first element is the missing data
#'       matrix and the next three are the ComplexHeatmap plot objects.
#'     \item "totint"- A list, equal in length to the input normList, where each
#'       element of the list gives the total intensity for each sample for a
#'       given normaliztation method.
#'     \item "dir"- The directory for the report to be saved in.
#'     \item "file"- The file for the report to be saved in.
#'  }
#'
#'
#' @export
#'
#' @seealso \code{\link{proteinormMetricBoxplot}},
#'   \code{\link{plotLogRatioDensity}}, \code{\link{plotTotInten}},
#'   \code{\link{plotHeatmapsForReport}},
#'   \code{\link{norm_metrics}},
#'
#' @examples
#' # No examples yet
#'
make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   batch = NULL,
                                   sampleLabels = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   save = TRUE,
                                   overwrite = FALSE,
                                   dir = NULL,
                                   file = NULL,
                                   keep.png = FALSE,
                                   showAllProteins = FALSE) {

  cli::cli_rule()

  #################################
  ## Check args and set defaults ##
  #################################
  enrich <- rlang::arg_match(enrich)

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
    file <- "proteiNorm_Report.pdf"
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

  ####################
  ## Make the plots ##
  ####################

  # Boxplots of metrics
  metrics_to_boxplot <- c("PCV", "PMAD", "PEV", "COR")
  boxplot_plotting_res <- lapply(X = metrics_to_boxplot,
         FUN = proteinormMetricBoxplot,
         normList = normList,
         groups = groups,
         batch = batch,
         dir = out_dir,
         save = save)
  names(boxplot_plotting_res) <- stringr::str_to_lower(metrics_to_boxplot)

  # Log ratio plots
  lograt <- plotLogRatioDensity(normList = normList,
                                groups = groups,
                                batch = batch,
                                zoom=FALSE,
                                legend = TRUE,
                                inset = 0.02,
                                dir = out_dir,
                                save = save)
  lograt2 <- plotLogRatioDensity(normList = normList,
                                 groups = groups,
                                 batch = batch,
                                 zoom = TRUE,
                                 legend = TRUE,
                                 inset = 0.02,
                                 dir = out_dir,
                                 save = save)
  # Intensity
  totint <- plotTotInten(normList = normList,
                         groups = groups,
                         batch = batch,
                         sampleLabels = sampleLabels,
                         dir = out_dir,
                         save = save)

  ## Heatmap(s)
  nahm   <- plotHeatmapsForReport(data = normList[[1]],
                                  groups = groups,
                                  batch = batch,
                                  sampleLabels = sampleLabels,
                                  showAllProteins = showAllProteins,
                                  dir = out_dir,
                                  save = save)


  ###########################
  ## If save, make report  ##
  ## and deal with .pngs   ##
  ###########################
  if (save) {

    cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")

    grDevices::pdf(file.path(out_dir,file),
        paper = "USr",
        pagecentre = TRUE,
        pointsize = 10,
        width = 12,
        height = 8)
    # TODO: clean this up to work variably depending on which PNGs are produced?
    png_names <- c("PCVplot.png", "PMADplot.png", "PEVplot.png", "CORplot.png", "Log2RatioPlot.png",
               "Log2RatioPlot-zoom.png", "NaHMplot.png", "NaHMplot_clust.png", "NaHMplot_group.png",
               "NaHMplot_batch.png", "TotIntenPlot.png")
    png_paths <- paste0(paste0(file.path(out_dir), "/"), png_names)

    # Grab the plots
    thePlots <- lapply(1:length(png_paths), function(i) {grid::rasterGrob(png::readPNG(png_paths[i], native=F))})

    # Assemble on PDF
    # Make the first page, with the 6 main plots on one page
    gridExtra::grid.arrange(grobs = thePlots[1:6], ncol = 3)
    # Then print all plots on their own page
    lapply(X = thePlots, FUN = gridExtra::grid.arrange, ncol = 1)
    grDevices::dev.off()

    # Deal with the .png files
    if (!keep.png) { # Delete
      unlink(png_paths)
      cli::cli_inform("Temporary {.path .png} files removed from {.path {out_dir}} ")
    } else { # Move to a new png directory
      # Setup dir name
      pngdir <- gsub(".pdf", "_pngs", file)
      # Create if it doesn't exist
      if (!dir.exists(file.path(out_dir, pngdir))) {
        dir.create(file.path(out_dir, pngdir), recursive = TRUE)
      }
      # Move files
      lapply(png_paths, function(x) {
        file.copy(from = x, to = file.path(out_dir, pngdir, basename(x)))
        file.remove(x)
      })

      cli::cli_inform("{.path .png} files saved to {.path {file.path(out_dir, pngdir)}}")
    }

  }

  # Return data from the reports invisibly.
  # TODO: reconsider whether we return anything at all?
  output_data <-  c(boxplot_plotting_res,
                    list(lograt = lograt,
                         nahm = nahm,
                         totint = totint,
                         dir = out_dir,
                         file = file))
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  return(invisible(output_data))
}




#' Make boxplot of a proteinorm metric
#'
#' Makes, and optionally saves, a boxplot showing the selected normalization
#' metric across samples. See \code{\link{norm_metrics}} for info on
#' the available metrics.
#'
#' @inheritParams make_proteinorm_report
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param metric The normalization metric to calculate and plot. Can be "PCV",
#'   "PMAD", "PEV", or "COR". See \code{\link{norm_metrics}}.
#' @param dir The directory in which to save the plot, if saving. Default is the
#'   current working directory.
#' @param save Should the plot be saved (as a .png)? Default is FALSE.
#'
#' @return A list, giving the results of applying the given normalization metric
#'   to the input normList.
#' @export
#'
#' @seealso \code{\link{norm_metrics}}
#'
#' @examples
#' # No examples yet
#'
proteinormMetricBoxplot <- function(normList,
                                    metric = c("PCV", "PMAD", "PEV", "COR"),
                                    groups,
                                    batch = NULL,
                                    dir = ".",
                                    save = FALSE) {

  # check args
  metric <- rlang::arg_match(metric)

  if (is.null(batch)) {
    batch <- c(rep("1", ncol(normList[[1]])))
  }

  # coerce group and batch to factor???
  # Not sure we need to do this
  # these are either (1) already a factor,
  # and/or (2) the functions of make use of these don't require them to be factors?
  #groups <- make_factor(as.character(groups))
  #batch <- make_factor(as.character(batch), prefix = NULL)
  # So far, seems to work just fine without any coercion happening here.

  # Set up the titles, filenames, etc. that we'll use for each metric
  plot_labelling <- data.frame(type = c("PCV", "PMAD", "PEV", "COR"),
                               basefile = paste0(c("PCV", "PMAD", "PEV", "COR"), "plot.png"),
                               main = c("PCV", "PMAD", "PEV", "COR"),
                               yaxis = c("Pooled Coefficient of Variation",
                                         "Median Absolute Deviation",
                                         "Pooled Estimate of Variance",
                                         "Intragroup Correlation"))
  type <- NULL
  labels <- base::subset(plot_labelling, type == metric)


  # If saving, make sure directory exists then open the plotting device
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    grDevices::png(filename = file.path(dir, labels$basefile),
        units = "px",
        width = 650,
        height = 650,
        pointsize = 15)
  }


  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- graphics::par()$mar
  on.exit(graphics::par(mar = old_mar), add = TRUE)

  # Get plot data by applying our metric function across
  # the input normlist
  plotData <- base::lapply(normList, FUN = metric, groups = groups)

  # Make the plot
  graphics::par(mar = c(8, 6, 4, 3))
  # main plot
  graphics::boxplot(x = plotData,
          main = labels$main,
          las = 2,
          col = binfcolors[1:length(normList)],
          boxlwd = 1,
          yaxt = "n",
          xaxt = "n",
          cex.main = 1.5)
  # Y axis
  graphics::axis(side = 2, cex.axis = 1.2, las = 2)
  # X axis
  graphics::axis(side = 1,
       at = seq_along(names(normList)),
       labels = names(normList),
       cex.axis = 1.3,
       las = 2)
  # Y axis text
  graphics::mtext(side = 2,
        text = labels$yaxis,
        line = 4.5,
        cex = ifelse(save, 1.2, 0.9))
  # Points
  graphics::points(rep(seq_along(normList),
             each = length(plotData[[1]])),
         unlist(plotData),
         pch = "*",
         cex = 1)

  if (save) grDevices::dev.off()

  return(invisible(plotData))
}



#' Make a plot of the density of the log2ratio
#'
#' Makes, and optionally saves, a plot showing the density of the log2ratio.
#' Uses \code{\link{log2ratio}} for calculating the ratio.
#'
#' @inheritParams make_proteinorm_report
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param batch Optional. A character or factor vector, listing the batch(es)
#'   the samples belong to.
#' @param zoom Should the plot cover the full range of log2ratios, or zoom
#'   in around 0? Default is FALSE.
#' @param legend Include a legend in the plot? Default is TRUE.
#' @param inset Passed to \code{\link[graphics:legend]{graphics::legend}}, the
#'   inset distance from the margin as a fraction of the plot region. Default is
#'   0.02.
#' @param dir The directory in which to save the plot, if saving. Default is the
#'   current working directory.
#' @param save Should the plot be saved (as a .png)? Default is FALSE.
#'
#' @return A list, equal in length to the input normList, where each
#'       element of the list is the density object of the log2ratio density for a
#'       given normalization method. See \code{\link{log2ratio}}.
#'
#' @export
#'
#' @seealso \code{\link{log2ratio}}
#'
#' @examples
#' # No examples yet
#'
plotLogRatioDensity <- function(normList,
                                groups,
                                batch = NULL,
                                zoom = FALSE,
                                legend = TRUE,
                                inset = 0.02,
                                dir = ".",
                                save = FALSE) {
  # Prep args
  # Again, not sure we need to coerce to factor here
  groups <- make_factor(as.character(groups))


  if (is.null(batch)) {
    batch <- c(rep("1", ncol(normList[[1]])))
  }
  batch <- make_factor(as.character(batch), prefix = NULL)


  # Calculate the log2ratios for each element of the Normlist
  plotData <- lapply(normList, FUN = log2ratio, groups = groups)

  # Set up plotting area limits
  maxY <- max(unlist(base::lapply(plotData, FUN=function(x) max(stats::density(x, na.rm = T)$y))))
  minY=0

  if (zoom) {
    minX <- -0.3
    maxX <- 0.3
    maxY <- maxY + (0.2*maxY)
    minY <- maxY - (0.5*maxY)
  } else {
    minX <- 0.5 * min(unlist(min(stats::density(plotData[["vsn"]], na.rm = T)$x)))
    maxX <- 0.5 * max(unlist(max(stats::density(plotData[["vsn"]], na.rm = T)$x)))
  }

  # If saving, set up dir and filename
  # and open png
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    filename <- file.path(dir, paste0("Log2RatioPlot", ifelse(zoom,"-zoom.png", ".png")))
    grDevices::png(filename = filename,
        units = "px",
        width = 650,
        height = 650,
        pointsize = 15)
  }

  # Collect current par() options that we're going to change,
  # and set them back on exit
  ### FOR REASON THAT ARE SEMI-UNCLEAR TO ME,
  ### ALL THIS PAR STUFF HAS TO COME AFTER THE PNG MAKING.
  ### PAR WIL MAKE A NEW GRAPHICS DEVICE IF THERE ISN"T ALREADY ONE OPEN
  ### SO, WHEN IT IS IN FRONT, IT OPENS A GRAPHICS DEVICE.
  ### THEN, PNG OPENS ANOTHER GRAPHICS DEVICE. AND CLOSING ON EXIT
  ### WILL ONLY CLOSE ONE OF THEM.
  old_mar <- graphics::par()$mar
  on.exit(graphics::par(mar = old_mar), add = TRUE)

  if (legend) {
    graphics::par(mar=c(5,5,4,3))
  } else {
    graphics::par(mar=c(5,5,3,4))
  }
  # Initialize empty plot
  base::plot(NA, las = 1,
       xlim = c(minX, maxX),
       ylim = c(minY, maxY),
       xlab = "Log2 ratio",
       ylab = "Density",
       main = "Log2-ratio",
       cex.main = 1.5,
       cex.axis = 1.2,
       cex.lab = 1.3)
  graphics::abline(v = 0, lwd = 2, lty = 3, col = "grey")


  # Plot each density line
  densityList <- list()
  for (method in names(normList)) {
    graphics::lines(stats::density(plotData[[method]], na.rm = T),
          col = binfcolors[which(names(plotData) %in% method)],
          lwd = 3)

    densityList[[method]] <- stats::density(plotData[[method]], na.rm = T)
  }

  # Plot legend if including
  if (legend) {
    graphics::legend("topright",
           inset = c(inset, 0),
           names(plotData),
           bty = "n",
           xpd = TRUE,
           box.col = "transparent",
           box.lwd = 0,
           bg = "transparent",
           border = "transparent",
           col = "transparent",
           pch = 22,
           pt.bg = binfcolors[1:length(plotData)],
           pt.cex = 1.5,
           cex = 1,
           horiz = FALSE,
           ncol = 1)
  }
  if (save) grDevices::dev.off()

  return(invisible(densityList))
}



#' Make plots of total intensity for each sample
#'
#' Makes, and optionally saves, a set of plot showing the total intensity for
#' each sample across the normalization methods.
#'
#' @inheritParams proteinormMetricBoxplot
#' @param sampleLabels Optional, a set of sample labels to use. If not supplied,
#'   defaults to using the column names in the normList.
#' @param dir The directory in which to save the plot, if saving. Default is the
#'   current working directory.
#' @param save Should the plot be saved (as a .png)? Default is FALSE.
#'
#' @return A list, equal in length to the input normList, where each
#'       element of the list gives total intensity for each sample for a
#'       given normalization method.
#'
#' @export
#'
#' @examples
#' # No examples yet
#'
plotTotInten <- function(normList,
                         groups,
                         batch = NULL,
                         sampleLabels = NULL,
                         dir = ".",
                         save = FALSE) {
  # Prep args
  # Again, not sure we need to coerce to factor here
  groups <- make_factor(groups)
  if (is.null(batch)) {
    batch <- c(rep("1",ncol(normList[[1]])))
  }
  batch <- make_factor(as.character(batch), prefix = NULL)
  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(normList[[1]])
  }

  # Set up plotting area, variably by number of samples
  if (length(groups) < 100) {
    width <- round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02, 0)
    height <- 800
    ncols <- 3
  } else {
    width <- round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15, 0)
    height <- 2400
    ncols <- 1
  }


  # If saving, set up files
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    grDevices::png(filename = file.path(dir, "TotIntenPlot.png"),
        units = "px",
        width = width,
        height = height,
        pointsize = 15)
  }

  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- graphics::par()$mar
  old_oma <- graphics::par()$oma
  on.exit(graphics::par(mar = old_mar), add = TRUE)
  on.exit(graphics::par(oma = old_oma), add = TRUE)
  on.exit(graphics::layout(matrix(1)), add = TRUE)

  # Set up plotting parameters
  # Have to do this twice, can't have the pars above the png making
  if (length(groups) < 100) {
    graphics::par(oma = c(2, 1, 1, 1),
                  mar = c(8, 5, 5, 2))
  } else {
    graphics::par(oma = c(1, 5, 5, 5),
                  mar = c(8, 2, 2, 2))
  }


  # Reorder group affiliations
  # Since groups is an ordered factor (see make_factor()),
  # this puts them in the order of their levels, then sorts by name.
  # Does not necessarily sort alphabetically.
  group_order <- sort.int(groups, index.return = T)$ix

  # Then, reorder the sample labels and groups
  sampleLabels <- sampleLabels[group_order]
  groups <- groups[group_order]
  # Cols in data reordered below


  # Make a plot for each element of normList
  graphics::layout(matrix(1:9, ncol = ncols, byrow = TRUE))
  barList <- NULL
  for (i in names(normList)) {
    barList[[i]] <- colSums(normList[[i]], na.rm = T)
    graphics::barplot(barList[[i]][group_order],
            main = "",
            las = 2,
            yaxt = "n",
            cex.main = 1.5,
            cex.lab = 1.2,
            col = colorGroup2(groups)[groups],
            names.arg = sampleLabels)
    graphics::title(main = i, font.main = 1, cex.main = 1.5, line = 2)
    graphics::axis(side = 2, cex.axis = 1.2, las = 2)
    if (i == "VSN") {
      graphics::mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
    }
  }
  names(barList) <- names(normList)


  if (save) grDevices::dev.off()

  return(invisible(barList))
}



#' Create a missing values heatmap
#'
#' Makes a ComplexHeatmap object showing a heatmap of missing values in the
#' input data.
#'
#' @inheritParams plotTotInten
#' @param missing An input matrix describing the missing data. See internals
#'   of \code{\link{plotHeatmapsForReport}}.
#' @param column_sort How should the columns of the heatmap be sorted? Options
#'   are: "cluster"- sort by similarity in missing values, "group"- sort samples
#'   by group, "batch"- sort samples by batch.
#' @param groupColors A vector of colors to use for coloring samples by group
#' @param batchColors A vector of colors to use for coloring samples by batch.
#' @param legend Include a legend in the plot? Default is FALSE.
#'
#' @return A ComplexHeatmap object
#' @export
#'
#' @examples
#' # No examples yet.
#'
missingValueHeatmap <- function(missing,
                                groups,
                                batch = NULL,
                                column_sort = c("cluster", "group", "batch"),
                                groupColors,
                                batchColors,
                                sampleLabels = NULL,
                                legend = FALSE) {
  # Check args
  column_sort <- rlang::arg_match(column_sort)

  # Deal with batch and sample names if null
  if (is.null(batch)) {
    batch <- c(rep("1", ncol(missing)))
  }

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(missing)
  }

  # Set up the column ordering for the different options
  if (column_sort == "cluster") {
    # If clustering just order things as they are
    cluster <- T
    order <- seq_along(groups)
    title <- "Cluster"
  } else if (column_sort == "group") {
    # If ordering by group, sort by group then batch
    cluster <- F
    order <- order(groups, batch)
    title <- "Sorted By Group"
  } else if (column_sort == "batch") {
    # If ordering by batch, sort by bathc then group
    cluster <- F
    order <- order(batch, groups)
    title <- "Sorted By Batch"
  } else {
    cli::cli_abort("Invalid value for {.arg column_sort}: cannot be {.val {column_sort}}")
  }

  # Then order the groups/batches/etc
  sorted_groups <- groups[order]
  sorted_batches <- batch[order]
  sorted_labels <- sampleLabels[order]
  ordered_data <- missing[,order]

  # Set up heatmap annotation
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(
    Sample = sorted_groups,
    Batch = sorted_batches,
    col = list(Sample = groupColors,
               Batch = batchColors),
    annotation_legend_param = list(Sample = list(title = "Group",
                                                 at = unique(sorted_groups),
                                                 labels = paste("", unique(sorted_groups))),
                                   Batch = list(title = "Batch",
                                                at = unique(sorted_batches),
                                                labels = paste("Batch", unique(sorted_batches)))),
    show_legend = legend
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
    show_heatmap_legend = legend,
    top_annotation = ColAnn,
    cluster_columns = cluster,
    column_labels = sorted_labels
  )

  # Return the heatmap
  return(hm_clust)
}


#' Make a set of missing value heatmaps
#'
#' Makes, and optionally saves, a set of missing value heatmaps using
#' \code{\link{missingValueHeatmap}}.
#'
#' @inheritParams plotTotInten
#' @param data A normalized data matrix. Generally, the "log2" slot of the normList
#'   slot in the list output by \code{\link{process_data}}.
#' @param dir The directory in which to save the plot, if saving. Default is the
#'   current working directory.
#' @param showAllProteins For missing data heatmaps, show all proteins (including
#'  ones with no missing data)? Default is FALSE.
#' @param save Should the plot be saved (as a .png)? Default is FALSE.
#'
#' @return A list of length 4, where the first element is the missing data
#'       matrix and the next three are the ComplexHeatmap plot objects sorted
#'       by cluster, group, and batch, respectively.
#' @export
#'
#' @examples
#' # No examples yet.
#'
plotHeatmapsForReport <- function(data,
                                  groups,
                                  batch = NULL,
                                  sampleLabels = NULL,
                                  dir = ".",
                                  showAllProteins = FALSE,
                                  save = FALSE) {
  # Process arguments
  if (is.null(batch)) {
    batch <- c(rep("1",ncol(data)))
  }

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive=TRUE)
  }

  # Set up colors
  batchCol <- colorBatch(batch)
  groupCol <- colorGroup2(groups)

  #Prepare data
  missing <- !is.na(data)
  if (!showAllProteins) {
    complete = apply(missing, 1, all)
    completeNA = apply(!missing, 1, all)
    missing <- missing[!complete & !completeNA,]
  }



  # Set up width for
  if (length(groups) <100 ) {
    width <- round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02, 0)
    width_single <- width/200
    width_together <- width/72
  } else {
    width <- round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15, 0)
    width_single <- width/90
    width_together <- width/30
  }

  # If we're saving, make individual plots of each type for merging
  # with slightly different formatting than the plots alone.
  # Will only do this if saving (trying to display the three plots all together
  # in the Rstudio window is too small, no point in doing if save == F)
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot.png"),
        units="in",
        width = width_together,
        height = 8,
        res = 100,
        pointsize = 8)

    # Make each plot individually, plotting together
    hm_clust <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "cluster",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = F)
    hm_group <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "group",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = F)
    hm_batch <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "batch",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = T)

    # Giving these objects unique names in the
    # name slot, to avoid a warning. Doesn't seem to affect
    # the plot at all
    hm_clust@name <- "a"
    hm_group@name <- "b"
    hm_batch@name <- "c"

    ComplexHeatmap::draw(hm_clust + hm_group + hm_batch,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         ht_gap = grid::unit(2, "cm"),
                         column_title = "Missing Values")
    grDevices::dev.off()
  }


  # Then, do individual plots.
  # Cluster samples
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_clust.png"),
        units = "in",
        width = width_single,
        height = 8,
        res = 100,
        pointsize = 8)
  }
  hm_clust <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "cluster",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  # Will draw these whether we're saving or not
  ComplexHeatmap::draw(hm_clust)
  if (save) grDevices::dev.off()


  # Sample by group
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_group.png"),
        units = "in",
        width = width_single,
        height = 8,
        res = 100,
        pointsize = 8)
  }

  hm_group <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "group",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  ComplexHeatmap::draw(hm_group)
  if (save) grDevices::dev.off()


  # Sample by batch
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_batch.png"),
        units = "in",
        width = width_single,
        height = 8,
        res = 100,
        pointsize = 8)
  }

  hm_batch <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "batch",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  ComplexHeatmap::draw(hm_batch)
  if (save) grDevices::dev.off()



  data2 <- list(missing = missing,
                hm_clust = hm_clust,
                hm_group = hm_group,
                hm_batch = hm_batch)

  return(invisible(data2))
}

# TODO: Could do further decomposition of these functions, which are not at the same
# level. Right now, we have a sort of generic heatmapper that returns the plot object.
# And it has a helper function that does the plot saving around it.
# Our other functions (boxplots, etc) are fully contained: they make and save the plot
# and return data, not plot objects. would prefer to make them all like the heatmaps,
# but not a super-high priority at the moment. Might necessitate a switch to ggplot,
# which I would probably prefer anyway.

