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
#' @param processed_data The output of the \code{\link{process_data}} function:
#'   a list object containing processed data and sample information..
#' @param grouping_column The name of column within the targets data frame which
#'   gives information on how to group samples for normalization. Must be supplied:
#'   some metrics can't be calculated for only one group.
#' @param enrich What type of analysis is this for? Options are
#'   "protein" and "phospho". Only used for making output dir name when one isn't
#'   supplied.
#' @param overwrite Should report file be overwritten if it already exists?
#'   Default is FALSE.
#' @param out_dir The directory in which to save the report. If not provided,
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
#' @importFrom ggplot2 ggsave
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
write_proteinorm_report <- function(processed_data,
                                   grouping_column = NULL,
                                   enrich = c("protein", "phospho"), #TODO: only used for making dir name...
                                   out_dir = NULL,
                                   file = NULL,
                                   overwrite = FALSE,
                                   suppress_zoom_legend = FALSE) {

  cli::cli_rule()

  #################################
  ## Check args and set defaults ##
  #################################

  # Check that processed_data has expected list structure
  if (!all(c("normList", "targets", "filt", "param", "stats") %in% names(processed_data))) {
    cli::cli_abort(c("{.arg processed_data} does not have expected structure:",
                     "i" = "Is it the object created by running {.code process_data()}?."))
  }

  enrich <- rlang::arg_match(enrich)

  # If provided, check that grouping column exists in the target dataframe
  # And set it
  if (is.null(grouping_column)) { # If no groups provided, abort
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  } else { # check that grouping column exists in the target dataframe
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))

    }
    if (grouping_column %notin% colnames(processed_data$targets)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in the targets dataframe of  in {.arg processed_data}",
                       "i" = "Check the column names with {.code colnames(processed_data$targets)}."))
    }
    # And set it
    groups <- as.character(processed_data$targets[,grouping_column])
    # And give error if there's only one group
    if (length(unique(groups)) < 2) {
      cli::cli_abort(c("Column {.arg {grouping_column}} does not contain at least two different groups",
                       "!" = "Cannot calculate all normaliztion metrics without at least two groups"))
    }
  }

  # Set default dir if not provided
  if (is.null(out_dir)) {
    out_dir <- file.path(paste0(enrich, "_analysis"), "01_quality_control")
    cli::cli_inform(cli::col_yellow("{.arg out_dir} argument is empty. Setting output directory to: {.path {out_dir}}"))
  }
  # Set default report name if not provided
  if (is.null(file)) {
    file <- "proteiNorm_Report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg file} argument is empty. Saving report to: {.path {out_dir}/{file}}"))
  }

  # Extract normList from input data
  normList <- processed_data$normList

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

  ###########################
  ## Make the plot objects ##
  ###########################

  # First page: all the combined plotting metrics
  a <- pn_plot_PCV(normList, groups)
  b <- pn_plot_PMAD(normList, groups)
  c <- pn_plot_PEV(normList, groups)
  d <- pn_plot_COR(normList, groups)
  e <- pn_plot_log2ratio(normList, groups)
  f <- pn_plot_log2ratio(normList, groups, zoom = T, legend = !suppress_zoom_legend)



  page_1 <- a + b + c + d + e + f +
    plot_layout(ncol = 3)

  # Second page: the faceted MD plots
  page_2 <- pn_plot_MD(normList, groups)

  ###############################
  ## Save plots, check, return ##
  ###############################

  # To save over multiple PDF pages
  # Make a list of plots
  # Need to convert page 1 from a patchwork object
  # into a Grob
  plots_list <-  list(patchwork::patchworkGrob(page_1), page_2)

  # Then save
  cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")
  ggsave(file.path(out_dir, file),
         plot = gridExtra::marrangeGrob(grobs = plots_list, nrow = 1, ncol = 1, top = NA),
         height = 8.5,
         width = 11,
         units = "in")

  if (!file.exists(file.path(out_dir, file))) {
    cli::cli_abort(c("Failed to create {.path {file.path(out_dir, file)}}"))
  }
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(file.path(out_dir, file))
}

