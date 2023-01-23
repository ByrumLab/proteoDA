#' Filter samples from a DAList
#'
#' This function is used to remove samples from a DAList, filtering using data
#' in the metadata data frame of the DAList. Samples which do not produce a
#' value of TRUE for the supplied condition are removed from the data, annotation,
#' and metadata slots of the DAList. If condition evaluates to NA, the function
#' will return an error.
#'
#' @param DAList A DAList object to be filtered.
#' @param condition An expression that returns a logical value, defined in terms
#'   of variables present in the metadata data frame of the supplied DAList.
#'   Samples are kept if the condition is TRUE for that sample.
#'
#'
#' @return A DAList, with samples that do not meet the condition are removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose the DAList$metadata data frame contains three columns:
#' # sample_ID = An alpha-numeric ID uniquely identifying a sample
#' # treatment = A character listing the treatment group a sample belongs to
#' # year = A numeric value listing the year a sample was collected
#'
#' # Remove a specific sample by ID
#' filtered <- filter_samples(DAList,
#'                            sample_ID != "abc123")
#'
#' # Remove any sample which contains the string "control" in the treatment:
#' filtered <- filter_samples(DAList,
#'                            grepl(pattern = "control",
#'                                  x = treatment))
#'
#' # Remove any sample from before 2010
#' filtered <- filter_samples(DAList,
#'                            year >= 2010)
#'
#'
#' # Filtering functions can be chained together
#' filtered <- DAList |>
#'   filter_samples(grepl(pattern = "control",
#'                        x = treatment)) |>
#'   filter_samples(year >= 2010)
#' }

filter_samples <- function(DAList, condition) {

  if (!(class(DAList) %in% c("DAList"))) {
    cli::cli_abort("{.arg DAList} must be a DAList object")
  }

  if (is.null(DAList$metadata)) {
    cli::cli_abort("{.arg DAList} does not contain metadata for filtering samples")
  }



  # get input metadata
  in_meta <- DAList$metadata

  # Add a col, keep, with the evaluation of the condition expression.
  condition_call <- substitute(condition)
  in_meta <- within(in_meta, {
    sample_to_be_kept <- eval(condition_call)
  })

  # Do subseting, keeping removed samples for stats.
  meta_removed <- subset(in_meta, subset = !sample_to_be_kept, select = c(-sample_to_be_kept))
  meta_kept <- subset(in_meta, subset = sample_to_be_kept, select = c(-sample_to_be_kept))

  # Check that sample numbers match
  if (nrow(meta_removed) + nrow(meta_kept) != nrow(in_meta)) {
    cli::cli_abort(c("Issue when filtering samples",
                     "!" = "Rows kept + rows removed != rows input",
                     "i" = "Is there an {.val NA} in the column you're filtering on?"))
  }

  # Update metadata samples
  DAList$metadata <- meta_kept
  # Update data, removing cols that are no longer present
  DAList$data <- DAList$data[, rownames(DAList$metadata)]
  validate_DAList(DAList)


  # print messages
  cli::cli_inform("Removed {.val {nrow(meta_removed)}} of the {.val {nrow(in_meta)}} {cli::qty(nrow(in_meta))} sample{?s} in DAList")
  cli::cli_inform("{cli::qty(nrow(meta_removed))} Sample{?s} removed:")
  message(paste0(utils::capture.output(meta_removed), collapse = "\n"))

  DAList
}
