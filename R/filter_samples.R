#' Filter samples from a DIAlist
#'
#' This function is used to remove samples from a DIAlist, filtering using data
#' in the metadata dataframe of the DIAlist. Samples which do not produce a
#' value of TRUE for the supplied condition are removed from the data, annotation,
#' and metadata slots of the DIAlist. If condition evaluates to NA, the function
#' will return an error.
#'
#' @param DIAlist A DIAlist object to be filtered.
#' @param condition An expression that returns a logical value, defined in terms
#'   of variables present in the metadata dataframe of the supplied DIAlist.
#'   Samples are kept if the condition is TRUE for that sample.
#'
#'
#' @return A DIAlist, with samples that do not meet the condition removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose the DIAlist$metadata data frame contains three columns:
#' # sample_ID = An alpha-numeric ID uniquely identifying a sample
#' # treatment = A character listing the treatment group a sample belongs to
#' # year = A numeric value listing the year a sample was collected
#'
#' # Remove a specific sample by ID
#' filtered <- filter_samples(DIAlist,
#'                            sample_ID != "abc123")
#'
#' # Remove any sample which contains the string "control" in the treatment:
#' filtered <- filter_samples(DIAlist,
#'                            grepl(pattern = "control",
#'                                  x = treatment))
#'
#' # Remove any sample from before 2010
#' filtered <- filter_samples(DIAlist,
#'                            year >= 2010)
#'
#'
#' # Filtering functions can be chained together
#' filtered <- DIAlist |>
#'   filter_samples(grepl(pattern = "control",
#'                        x = treatment)) |>
#'   filter_samples(year >= 2010)
#' }

filter_samples <- function(DIAlist, condition) {

  if (!(class(DIAlist) %in% c("DIAlist"))) {
    cli::cli_abort("{.arg DIAlist} must be a DIAlist object")
  }

  if (is.null(DIAlist$metadata)) {
    cli::cli_abort("{.arg DIAlist} does not contain metadata for filtering samples")
  }



  # get input metadata
  in_meta <- DIAlist$metadata

  # Add a col, keep, with the evaluation of the condition expression.
  condition_call <- substitute(condition)
  in_meta <- within(in_meta, {
    keep <- eval(condition_call)
  })

  # Do subseting, keeping removed samples for stats.
  meta_removed <- subset(in_meta, subset = !keep, select = c(-keep))
  meta_kept <- subset(in_meta, subset = keep, select = c(-keep))

  # Check that sample numbers match
  if (nrow(meta_removed) + nrow(meta_kept) != nrow(in_meta)) {
    cli::cli_abort(c("Issue when filtering samples",
                     "!" = "Rows kept + rows removed != rows input",
                     "i" = "Is there an {.val NA} in the column you're filtering on?"))
  }

  # Update metadata samples
  DIAlist$metadata <- meta_kept
  # Update data, removing cols that are no longer present
  DIAlist$data <- DIAlist$data[, rownames(DIAlist$metadata)]
  validate_DIAlist(DIAlist)


  # print messages
  cli::cli_inform("Removed {.val {nrow(meta_removed)}} of the {.val {nrow(in_meta)}} {cli::qty(nrow(in_meta))} sample{?s} in DIAlist")
  cli::cli_inform("{cli::qty(nrow(meta_removed))} Sample{?s} removed:")
  print(meta_removed)

  DIAlist
}
