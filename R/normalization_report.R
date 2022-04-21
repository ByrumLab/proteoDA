make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   batch = NULL,
                                   sampleLabels = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   save = TRUE,
                                   dir = NULL,
                                   file = NULL,
                                   keep.png = FALSE,
                                   legend = TRUE) {


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
  nahm   <- plotNaHM(normList = normList,
                     groups = groups,
                     batch = batch,
                     sampleLabels = sampleLabels,
                     dir = out_dir,
                     save = save)


  ###########################
  ## If save, make report  ##
  ## and deal with .pngs   ##
  ###########################
  if (save) {

    cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")

    pdf(file.path(out_dir,file),
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
    do.call(gridExtra::grid.arrange, c(thePlots[1:6], ncol=3))
    # Then, put all plots onto individual pages
    for (i in seq_along(thePlots)) {
      do.call(gridExtra::grid.arrange, c(thePlots[i], ncol=1))
    }

    # Close device
    dev.off()

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
                    list(cor = cor,
                         lograt = lograt,
                         nahm = nahm,
                         totint = totint,
                         dir = out_dir,
                         file = file))
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  return(invisible(output_data))
}


### New generic proteinorm boxplot function.

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

  # Upon exiting the function, close any open graphics devices
  # Need to possibly be a little careful with this:
  # The idea is that, if there's an error, we want
  # to close the graphics device that has been created
  # within this function. This function closes all graphics devices
  # including any from higher up in the stack.
  # There probably shouldn't be any higher up in the stack, but I think
  # at the moment this is not true for the normalization report code.
  on.exit(graphics.off(), add = TRUE)

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
  labels <- base::subset(plot_labelling, type == metric)


  # If saving, make sure directory exists then open the plotting device
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    png(filename = file.path(dir, labels$basefile),
        units = "px",
        width = 650,
        height = 650,
        pointsize = 15)
  }

  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- par()$mar
  on.exit(par(mar = old_mar), add = TRUE)

  # Get plot data by applying our metric function across
  # the input normlist
  plotData <- base::lapply(normList, FUN = metric, groups = groups)

  # Make the plot
  par(mar = c(8, 6, 4, 3))
  # main plot
  boxplot(x = plotData,
          main = labels$main,
          las = 2,
          col = binfcolors[1:length(normList)],
          boxlwd = 1,
          yaxt = "n",
          xaxt = "n",
          cex.main = 1.5)
  # Y axis
  axis(side = 2, cex.axis = 1.2, las = 2)
  # X axis
  axis(side = 1,
       at = seq_along(names(normList)),
       labels = names(normList),
       cex.axis = 1.3,
       las = 2)
  # Y axis text
  mtext(side = 2,
        text = labels$yaxis,
        line = 4.5,
        cex = ifelse(save, 1.2, 0.9))
  # Points
  points(rep(seq_along(normList),
             each = length(plotData[[1]])),
         unlist(plotData),
         pch = "*",
         cex = 1)

  # If saving, close the device
  if (save) dev.off()

  return(invisible(plotData))
}



### UPdated plotting of log Ratio density
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

  # Close any graphics devices if there's an error below before we're done plotting
  on.exit(graphics.off(), add = TRUE)

  # Calculate the log2ratios for each element of the Normlist
  plotData <- lapply(normList, FUN = log2ratio, groups = groups)

  # Set up plotting area limits
  maxY <- max(unlist(base::lapply(plotData, FUN=function(x) max(density(x, na.rm = T)$y))))
  minY=0

  if (zoom) {
    minX <- -0.3
    maxX <- 0.3
    maxY <- maxY + (0.2*maxY)
    minY <- maxY - (0.5*maxY)
  } else {
    minX <- 0.5 * min(unlist(min(density(plotData[["vsn"]], na.rm = T)$x)))
    maxX <- 0.5 * max(unlist(max(density(plotData[["vsn"]], na.rm = T)$x)))
  }

  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- par()$mar
  on.exit(par(mar = old_mar), add = TRUE)
  if (legend) {
    par(mar=c(5,5,4,3))
  } else {
    par(mar=c(5,5,3,4))
  }

  # If saving, set up dir and filename
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    filename <- file.path(dir, paste0("Log2RatioPlot", ifelse(zoom,"-zoom.png", ".png")))
    png(filename = filename,
        units = "px",
        width = 650,
        height = 650,
        pointsize = 15)
  }

  # Initialize empty plot
  plot(NA, las = 1,
       xlim = c(minX, maxX),
       ylim = c(minY, maxY),
       xlab = "Log2 ratio",
       ylab = "Density",
       main = "Log2-ratio",
       cex.main = 1.5,
       cex.axis = 1.2,
       cex.lab = 1.3)
  abline(v = 0, lwd = 2, lty = 3, col = "grey")

  # Plot each density line
  densityList <- list()
  for (method in names(normList)) {
    lines(density(plotData[[method]], na.rm = T),
          col = binfcolors[which(names(plotData) %in% method)],
          lwd = 3)

    densityList[[method]] <- density(plotData[[method]], na.rm = T)
  }

  # Plot legend if including
  if (legend) {
    legend("topright",
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

  if (save) dev.off()

  return(invisible(densityList))

}



### Updated total intensity plot
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

  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- par()$mar
  old_oma <- par()$oma
  on.exit(par(mar = old_mar), add = TRUE)
  on.exit(par(oma = old_oma), add = TRUE)

  # Close any graphics devices if there's an error below before we're done plotting
  on.exit(graphics.off(), add = TRUE)

  # Set up plotting area, variably by number of samples
    if (length(groups) < 100) {
    width <- round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02, 0)
    height <- 800
    ncols <- 3
    par(oma = c(2, 1, 1, 1),
        mar = c(8, 5, 5, 2))
    } else {
      width <- round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15, 0)
      height <- 2400
      ncols <- 1
      par(oma = c(1, 5, 5, 5),
          mar = c(8, 2, 2, 2))
  }

  # If saving, set up files
  if (save) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
    }
    png(filename = file.path(dir, "TotIntenPlot.png"),
        units = "px",
        width = width,
        height = height,
        pointsize = 15)
  }


  # Make a plot for each element of normList
  layout(matrix(1:9, ncol = ncols, byrow = TRUE))
  barList <- NULL
  for (i in names(normList)) {
    barList[[i]] <- colSums(normList[[i]], na.rm = T)
    barplot(barList[[i]],
            main = "",
            las = 2,
            yaxt = "n",
            cex.main = 1.5,
            cex.lab = 1.2,
            col = colorGroup2(groups)[groups],
            names.arg = sampleLabels)
    title(main = i, font.main = 1, cex.main = 1.5, line = 2)
    axis(side = 2, cex.axis = 1.2, las = 2)
    if (i == "VSN") {
      mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
    }
  }
  names(barList) <- names(normList)
  if (save) dev.off()

  return(invisible(barList))
}

