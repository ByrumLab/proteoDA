## This file contains one function,
## which create the proteinorm report

## It uses functions in a variety of subfiles:
## normalization_metrics.R contains functions for numerically evaluating normalization methods
## normalization_plotting.R contains functions to plot these metrics

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



# a <- pn_plot_PCV(normList = norm_higgs$normList,
#                  grouping = norm_higgs$targets$group)
# b <- pn_plot_PMAD(normList = norm_higgs$normList,
#                   grouping = norm_higgs$targets$group)
# c <- pn_plot_PEV(normList = norm_higgs$normList,
#                   grouping = norm_higgs$targets$group)
# d <- pn_plot_COR(normList = norm_higgs$normList,
#                  grouping = norm_higgs$targets$group)
# e <- pn_plot_log2ratio(normList = norm_higgs$normList,
#                   grouping = norm_higgs$targets$group)
# f <- pn_plot_log2ratio(normList = norm_higgs$normList,
#                   grouping = norm_higgs$targets$group,
#                   zoom = T)
#
#
# a + b + c + d + e + f +
#   plot_layout(ncol = 3)
# ggsave(filename = "protein_analysis/01_quality_control/higgs_update.pdf", height = 8.7, width = 11, units = "in")
#



