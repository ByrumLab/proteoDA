

## write_norm_report uses functions in multiple other files:
## normalization_metrics.R contains functions for numerically evaluating normalization methods
## normalization_plotting.R contains functions to plot these metrics

#' Create a normalization report
#'
#' Saves a PDF report containing a variety of plots which give
#' information about the performance of different normalization metrics.
#' The report is useful for choosing a normalization method to
#' use for downstream analysis.
#'
#' @param DAList A DAList.
#' @param grouping_column The name of the column in the metadata which
#'   gives information on how to group samples for normalization. Must be supplied:
#'   some metrics can't be calculated for only one group.
#' @param overwrite Should report file be overwritten if it already exists?
#'   Default is FALSE.
#' @param output_dir The directory in which to save the report. If not provided,
#'   will default to "protein_analysis/01_quality_control". If the directory does not
#'   exist, it will be created.
#' @param filename The file name of the report to be saved. Must end in .pdf. Will
#'   default to "proteiNorm_Report.pdf" if no file name is provided.
#' @param suppress_zoom_legend Should the legend be removed from the zoomed
#'   log2ratio plot? Default is FALSE.
#' @param use_ggrastr Should the \code{ggrastr} package be used to decrease
#'   file size? Default is FALSE. Requires installation of \code{ggrastr}.
#'
#' @return If report is created successfully, invisibly returns the input DAList.
#'
#' @importFrom ggplot2 ggsave
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Group samples according to group identities
#' # in the "treatment" column of the metadata
#' write_norm_report(DAList,
#'                         grouping_column = "treatment")
#'
#' # Change the default directory and file names
#' write_norm_report(DAList,
#'                         grouping_column = "treatment",
#'                         output_dir = "my/chosen/directory",
#'                         filename = "my_report.pdf")
#'
#' # Overwrite an existing report
#' write_norm_report(DAList,
#'                         grouping_column = "treatment",
#'                         overwrite = T)
#' }
#'
write_norm_report <- function(DAList,
                              grouping_column = NULL,
                              output_dir = NULL,
                              filename = NULL,
                              overwrite = FALSE,
                              suppress_zoom_legend = FALSE,
                              use_ggrastr = FALSE) {



  #################################
  ## Check args and set defaults ##
  #################################

  input_DAList <- validate_DAList(DAList)

  if (!is.null(DAList$tags$normalized)) {
    if (DAList$tags$normalized) {
      cli::cli_abort("Data in DAList are already normalized. Cannot evaluate normalization metrics")
    }
  }

  # If provided, check that grouping column exists in the target data frame
  # And set it
  if (is.null(grouping_column)) { # If no groups provided, abort
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  } else { # check that grouping column exists in the target data frame
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))

    }
    if (grouping_column %notin% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$metadata)}."))
    }
    # And set it
    groups <- as.character(DAList$metadata[,grouping_column])
    # And give error if there's only one group
    if (length(unique(groups)) < 2) {
      cli::cli_abort(c("Column {.arg {grouping_column}} does not contain at least two different groups",
                       "!" = "Cannot calculate all normaliztion metrics without at least two groups"))
    }
  }

  # Set default dir if not provided
  if (is.null(output_dir)) {
    output_dir <- file.path("protein_analysis", "01_quality_control")
    cli::cli_inform(cli::col_yellow("{.arg output_dir} argument is empty. Setting output directory to: {.path {output_dir}}"))
  }
  # Set default report name if not provided
  if (is.null(filename)) {
    filename <- "proteiNorm_Report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg filename} argument is empty. Saving report to: {.path {output_dir}/{filename}}"))
  }


  ###########################
  ## Do all normalizations ##
  ###########################
  cli::cli_inform("Starting normalizations")
  normList <- apply_all_normalizations(DAList$data)
  cli::cli_inform("Normalizations finished")

  ##################
  ## Set up files ##
  ##################

 # Make directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # Check that the filename is a pdf
  validate_filename(filename, allowed_exts = c("pdf"))

  # If file already exists, inform about overwriting or through an error
  if (file.exists(file.path(output_dir, filename))) {
    if (overwrite) {
      cli::cli_inform("{.path {filename}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
    } else {
      cli::cli_abort(c("{.path {filename}} already exists in {.path {output_dir}}",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Give {.arg filename} a unique name or set {.arg overwrite} to {.val TRUE}"))
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
    patchwork::plot_layout(ncol = 3)

  # Second page: the faceted MD plots
  page_2 <- pn_plot_MD(normList, groups, use_ggrastr)

  ###############################
  ## Save plots, check, return ##
  ###############################

  # To save over multiple PDF pages
  # Make a list of plots
  # Need to convert page 1 from a patchwork object
  # into a Grob
  plots_list <-  list(patchwork::patchworkGrob(page_1), page_2)

  # Then save
  cli::cli_inform("Saving report to: {.path {file.path(output_dir, filename)}}")
  ggsave(file.path(output_dir, filename),
         plot = gridExtra::marrangeGrob(grobs = plots_list, nrow = 1, ncol = 1, top = NA),
         height = 8.5,
         width = 11,
         units = "in",
         useDingbats = T)

  if (!file.exists(file.path(output_dir, filename))) {
    cli::cli_abort(c("Failed to create {.path {file.path(output_dir, filename)}}"))
  }

  invisible(input_DAList)
}


#' Apply all normalization methods to a set of raw data
#'
#' Takes in raw intensity data and applies 8 different normalization methods,
#' returning a list of the data normalized with each method. Used internally in
#' \code{\link{write_norm_report}}.
#'
#' @param data A data frame or matrix of raw data to be normalized. Rows are proteins and
#'   columns are raw intensity data.
#'
#' @return A list length 8, where each item in the list is a
#'       named matrix. Names give the normalization method, and the matrix
#'       contains the normalized data.
#'
#' @keywords internal
#'
apply_all_normalizations <- function(data) {

  normList <- NULL

  ## apply normalization methods using functions listed above.
  ## NOTE: most of the normalizations use log2(intensity) as input except
  ## VSN which normalizes using raw intensities with no log2 transformation.
  normList[["log2"]]     <- log2Norm(dat = data)
  normList[["median"]]   <- medianNorm(logDat = normList[["log2"]])
  normList[["mean"]]     <- meanNorm(logDat = normList[["log2"]])
  normList[["vsn"]]      <- vsnNorm(dat = data)
  normList[["quantile"]] <- quantileNorm(logDat = normList[["log2"]])
  normList[["cycloess"]] <- cycloessNorm(logDat = normList[["log2"]])
  normList[["rlr"]]      <- rlrNorm(logDat = normList[["log2"]])
  normList[["gi"]]       <- giNorm(dat = data)

  normList
}

