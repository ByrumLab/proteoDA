make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   batch = NULL,
                                   sampleLabels = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   save = TRUE, # seems like default should be TRUE, not false
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
  lograt <- plotLogRatio(normList = normList,
                         groups = groups,
                         batch = batch,
                         sampleLabels = sampleLabels,
                         zoom=FALSE,
                         legend = TRUE,
                         inset = 0.02,
                         dir = out_dir,
                         save = save)
  lograt2 <- plotLogRatio(normList = normList,
                          groups = groups,
                          batch = batch,
                          sampleLabels = sampleLabels,
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
    # Close the device if there's an error below before we're done plotting
    on.exit(dev.off(), add = TRUE)
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
