make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   batch = NULL,
                                   sampleLabels = NULL,
                                   legend = TRUE,
                                   enrich = c("protein", "phospho"),
                                   dir = NULL,
                                   file = NULL,
                                   save = FALSE,
                                   keep.png = FALSE) {

  enrich <- match.arg(enrich,
                      choices = c("protein", "phospho"),
                      several.ok = FALSE)

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(normList[[1]])
  }

  if (is.null(groups)) {
    groups <- rep("group", ncol(normList[[1]]))
  }

  groups <- make_factor(x = as.character(groups))

  if (!is.null(batch)) {
    batch <- make_factor(as.character(batch))
  }

  if (!save) {
    keep.png<-FALSE
    dir <- "."
  }


  ## CREATE QC OUTPUT DIRECTORY
  ## if use enrich type to create QC output directory
  if (save == TRUE) {

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


    print(paste("QC output directory:", dir))

    ## FILENAME DEFINED
    if (!is.null(file)) {
      if (tools::file_ext(file) != "pdf") {
        stop("\nError! Invalid output file type...\nincorrect: file = '",file,"'",
             "\ncorrect:   file = 'proteiNorm_Report.pdf'")
      }

      pngdir = gsub(".pdf", "", file)

      if (file.exists(file.path(dir, file))) {
        file <- make_new_filename(x = file, dir = dir)
        no <- sub(".pdf", "", sub(".*_", "", file))
        pngdir <- gsub(".pdf", "", file)
      }
      if (!file.exists(file.path(dir, file))) {
        pngdir = gsub(".pdf", "", file)
      }

    } ## FILE NOT NULL

    ## FILENAME NULL
    if (is.null(file)) {
      file <-"proteiNorm_Report.pdf"
      no <- ""

      if (file.exists(file.path(dir, file))) {
        file <- make_new_filename(x = file, dir = dir)
        no <- sub(".pdf", "", sub(".*_", "", file))
        pngdir <- gsub(".pdf", "", file)
      }
      if (!file.exists(file.path(dir, file))) {
        pngdir <- gsub(".pdf", "", file)
      }
    }  ## FILE NULL

    print(paste("QC output directory:", dir))
    print(paste("QC report file: ", file))
    print(paste("png directory: ", pngdir))

  } ## SAVE ==TRUE


  ## CREATE PLOTS

  metrics_to_boxplot <- c("PCV", "PMAD", "PEV", "COR")
  boxplot_plotting_res <- lapply(X = metrics_to_boxplot,
         FUN = proteinormMetricBoxplot,
         normList = normList,
         groups = groups,
         batch = batch,
         dir = dir,
         save = save)
  names(boxplot_plotting_res) <- stringr::str_to_lower(metrics_to_boxplot)

  ## Heatmap(s)
  nahm   <- plotNaHM(normList = normList,
                     groups = groups,
                     batch = batch,
                     sampleLabels = sampleLabels,
                     dir = dir,
                     save = save)
  lograt <- plotLogRatio(normList = normList,
                         groups = groups,
                         batch = batch,
                         sampleLabels = sampleLabels,
                         zoom=FALSE,
                         legend = TRUE,
                         inset = 0.02,
                         dir = dir,
                         save = save)
  plotLogRatio(normList = normList,
               groups = groups,
               batch = batch,
               sampleLabels = sampleLabels,
               zoom = TRUE,
               legend = TRUE,
               inset = 0.02,
               dir = dir,
               save = save)

  totint <- plotTotInten(normList = normList,
                         groups = groups,
                         batch = batch,
                         sampleLabels = sampleLabels,
                         dir = dir,
                         save = save)
  if (!save) {
    ComplexHeatmap::draw(nahm$hm_batch)
    ComplexHeatmap::draw(nahm$hm_clust)
    ComplexHeatmap::draw(nahm$hm_group)
    }

  data2 <- c(boxplot_plotting_res,
             list(cor = cor,
                  lograt = lograt,
                  nahm = nahm,
                  totint = totint))
  if (save) {
    data2 <-  c(boxplot_plotting_res,
                list(cor = cor,
                     lograt = lograt,
                     nahm = nahm,
                     totint = totint,
                     dir = dir,
                     file = file))
    }

  ##  SAVE PROTEINORM REPORT PDF
  if (save) {
    pdf(file.path(dir,file),
        paper = "USr",
        pagecentre = TRUE,
        pointsize = 10,
        width = 12,
        height = 8)

    files <- c("PCVplot.png", "PMADplot.png", "PEVplot.png", "CORplot.png", "Log2RatioPlot.png",
               "Log2RatioPlot-zoom.png", "NaHMplot.png", "NaHMplot_clust.png", "NaHMplot_group.png",
               "NaHMplot_batch.png", "TotIntenPlot.png")
    pnglist <- paste0(paste0(file.path(dir), "/"), files)
    thePlots <- lapply(1:length(pnglist), function(i) {grid::rasterGrob(png::readPNG(pnglist[i], native=F))})

    do.call(gridExtra::grid.arrange, c(thePlots[1:6], ncol=3))
    do.call(gridExtra::grid.arrange, c(thePlots[1], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[2], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[3], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[4], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[5], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[7], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[8], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[9], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[10], ncol=1))
    do.call(gridExtra::grid.arrange, c(thePlots[11], ncol=1))

    dev.off()

    ## REMOVE PNG FILES
    if (!keep.png) {
      unlink(pnglist)
      print(file.path(dir,file))
      print("png files removed...")
    } ## remove png files.

    ## KEEP PNG FILES
    if (keep.png) {
      if (!dir.exists(file.path(dir, pngdir))) {
        dir.create(file.path(dir, pngdir), recursive = TRUE)
      }
      lapply(files, function(x) {
        file.copy(from = file.path(dir, x), to = file.path(dir, pngdir, x))
        file.remove(file.path(dir, x))
      })
      print(paste("png files moved to :", file.path(dir, pngdir)))
    } ## KEEP

  } ## SAVE == TRUE

  return(invisible(data2))
}




### New generic proteinorm boxplot function.

proteinormMetricBoxplot <- function(normList,
                                    metric = c("PCV", "PMAD", "PEV", "COR"),
                                    groups,
                                    batch = NULL,
                                    dir = ".",
                                    save = FALSE) {

  # check args
  rlang::arg_match(metric)

  if (is.null(batch)) {
    batch <- c(rep("1", ncol(normList[[1]])))
  }

  # coerce group and batch to factor???
  # Not sure we need to do this
  # these are either (1) alredy a factor,
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
