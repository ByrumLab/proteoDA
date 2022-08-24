## This file contains one function,
## which create the proteinorm report

## It uses functions in a variety of files:
## normalization_metrics.R contains functions for numerically evaluating normalization methods
## normalization_plotting.R contains functions to plot these metrics

#' Create proteinorm report
#'
#' Creates and saves as a PDF report a variety of plots which give
#' information about the performance of different normalization metrics.
#'
#' @param normList A list of normalized data matrices. Generally, the "normList"
#'   slot of the list that is output by \code{\link{process_data}}.
#' @param groups Optional, a vector describing what experimental or treatment
#'   group each sample belongs to. Generally, the "group" column in the targets
#'   data frame output by \code{\link{process_data}}. If not supplied, will
#'   warn the user and treat all samples as belonging to the same group.
#' @param enrich What type of analysis is this for? Options are
#'   "protein" and "phospho". Only used for making output dir name when one isn't
#'   supplied.
#' @param overwrite Should report file be overwritten if it already exists?
#'   Default is FALSE.
#' @param dir The directory in which to save the report. If not provided,
#'   will default to either "protein_analysis/01_quality_control" or
#'   "phospho_analysis/01_quality_control" within the current working
#'   directory, depending on the enrich argument.
#' @param file The file name of the report to be saved. Must end in .pdf. Will
#'   default to "proteiNorm_Report.pdf" if no filename is provided.
#' @param suppress_zoom_legend Should the legend be removed from the zoomed
#'   log2ratio plot? Default is FALSE
#'
#' @return Invisibly, the filename of the created report
#'
#' @export
#'
#' @seealso \code{\link{norm_metrics}},
#'   \code{\link{eval_pn_metric_for_plot}},
#'   \code{\link{pn_plots_generic}},
#'   \code{\link{pn_plots}},
#'
#' @examples
#' # No examples yet
#'
make_proteinorm_report <- function(normList,
                                   groups = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   dir = NULL,
                                   file = NULL,
                                   overwrite = FALSE,
                                   suppress_zoom_legend = FALSE) {

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
    cli::cli_inform(cli::col_yellow("{.arg dir} argument is empty. Setting output directory to: {.path {out_dir}}"))
  } else {
    out_dir <- dir
  }
  # Set default report name if not provided
  if (is.null(file)) {
    file <- "proteiNorm_Report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg file} argument is empty. Saving report to: {.path {out_dir}/{file}}"))
  }


  ##################
  ## Set up files ##
  ##################

 # Make directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }

  # Check that the filename is a pdf
  validate_filename(file, allowed_exts = c("pdf"))

  # If file already exists, inform about overwriting or through an error
  if (file.exists(file.path(out_dir, file))) {
    if (overwrite) {
      cli::cli_inform("{.path {file}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
    } else {
      cli::cli_abort(c("{.path {file}} already exists in {.path {out_dir}}",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Give {.arg file} a unique name or set {.arg overwrite} to {.val TRUE}"))
    }
  }

  ####################
  ## Make the plots ##
  ####################

  a <- pn_plot_PCV(normList, groups)
  b <- pn_plot_PMAD(normList, groups)
  c <- pn_plot_PEV(normList, groups)
  d <- pn_plot_COR(normList, groups)
  e <- pn_plot_log2ratio(normList, groups)
  f <- pn_plot_log2ratio(normList, groups, zoom = T, legend = !suppress_zoom_legend)

  combined <- a + b + c + d + e + f +
    plot_layout(ncol = 3)


  ##############################
  ## Save plot, check, return ##
  ##############################

  cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")
  ggsave(combined, filename = file.path(out_dir, file), height = 8.5, width = 11, units = "in")

  if (!file.exists(file.path(out_dir, file))) {
    cli::cli_abort(c("Failed to create {.path {file.path(out_dir, file)}}"))
  }
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(file.path(out_dir, file))
}

