make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   batch = NULL,
                                   sampleLabels = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   save = TRUE,
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
  }


  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- par()$mar
  on.exit(par(mar = old_mar), add = TRUE)
  on.exit(dev.off(), add = TRUE)

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

  return(invisible(plotData))
}



### Updated plotting of log Ratio density
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

  # If saving, set up dir and filename
  # and open png
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

  # Collect current par() options that we're going to change,
  # and set them back on exit
  ### FOR REASON THAT ARE SEMI-UNCLEAR TO ME,
  ### ALL THIS PAR STUFF HAS FO COME AFTER THE PNG MAKING.
  ### PAR WIL MAKE A NEW GRAPHICS DEVICE IF THERE ISN"T ALREADY ONE OPEN
  ### SO, WHEN IT IS IN FRONT, IT OPENS A GRAPHICS DEVICE.
  ### THEN, PNG OPENS ANOTHER GRAPHICS DEVICE. AND CLOSING ON EXIT
  ### WILL ONLY CLOSE ONE OF THEM.
  old_mar <- par()$mar
  on.exit(par(mar = old_mar), add = TRUE)
  on.exit(dev.off(), add = TRUE)

  if (legend) {
    par(mar=c(5,5,4,3))
  } else {
    par(mar=c(5,5,3,4))
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
    png(filename = file.path(dir, "TotIntenPlot.png"),
        units = "px",
        width = width,
        height = height,
        pointsize = 15)
  }


  # Collect current par() options that we're going to change,
  # and set them back on exit
  old_mar <- par()$mar
  old_oma <- par()$oma
  on.exit(par(mar = old_mar), add = TRUE)
  on.exit(par(oma = old_oma), add = TRUE)
  on.exit(dev.off(), add = TRUE)

  # Set up plotting parameters
  # Have to do this twice, can't have the pars above the png making
  if (length(groups) < 100) {
    par(oma = c(2, 1, 1, 1),
        mar = c(8, 5, 5, 2))
  } else {
    par(oma = c(1, 5, 5, 5),
        mar = c(8, 2, 2, 2))
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

  return(invisible(barList))
}


### New generic heatmap function
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
    png(filename = file.path(dir, "NaHMplot.png"),
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
    print(class(hm_clust))

    ComplexHeatmap::draw(hm_clust + hm_group + hm_batch,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         ht_gap = grid::unit(2, "cm"),
                         column_title = "Missing Values")
    dev.off()

  }


  # Then, do individual plots.
  # Cluster samples
  if (save) {
    png(filename = file.path(dir, "NaHMplot_clust.png"),
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
  dev.off()


  # Sample by group
  if (save) {
    png(filename = file.path(dir, "NaHMplot_group.png"),
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
  dev.off()


  # Sample by batch
  if (save) {
    png(filename = file.path(dir, "NaHMplot_batch.png"),
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
  dev.off()



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
# but not a super-high priority at the moment.

