#' Filter samples from a DIAlist
#'
#' Removes unwanted samples (rows) from the targets data frame This function is
#' basically \code{\link[base]{subset}} with extra functionality. It removes
#' rows if they contain certain values in a given column. By default, ignores
#' the case of strings in rm.val. Calls \code{\link{make_log}} as a subfunction.
#'
#' @param DIAlist The DIAlist to be filtered
#' @param filter_list A named list describing how to filter targets dataframe.
#'   The name for each element of the list gives a column to be filters, and
#'   each element of the list should be a character vector of strings to search for
#'   in the given column. See examples.
#' @param ignore.case Should the case of the strings in rm.vals be ignored? Is
#'   TRUE by default, such that case does not have to match. If FALSE, case must
#'   match for rows to be removed.
#'
#' @return A prototype of our new S3 list type.
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' DIAlist <- filter_samples(DIAlist
#'                           filter_list = list(group = "pool"
#'                                              sample = c("sampleA", "sampleB")))
#' }
#'

filter_samples <- function(DIAlist, filter_list, ignore.case = TRUE) {

  # check args
  if (!is.list(filter_list)) {
    cli::cli_abort(c("{.arg filter_list} must be a list"))
  }
  if (is.null(names(filter_list)) | length(names(filter_list)) == 0) {
    cli::cli_abort(c("{.arg filter_list} must be a named list"))
  }
  if (length(names(filter_list)) != length(filter_list)) {
    cli::cli_abort(c("Each list element in {.arg filter_list} must be named"))
  }

  cli::cli_rule()

  # Set up some vars to avoid R CMD check warning
  to_remove <- for_match <- remove_reason <- NULL

  # Mark for removal any rows that match/contain the supplied strings to exclude
  in_meta <- DIAlist$metadata
  in_meta$to_remove <- F



  for (column in names(filter_list)) {
    for (filter_string in unique(filter_list[[column]])) {
      if (ignore.case) {
        in_meta$to_remove <- stringr::str_to_lower(in_meta[,column]) %>%
          stringr::str_detect(stringr::str_to_lower(filter_string)) %>%
          ifelse(T, in_meta$to_remove)
      } else {
        in_meta$to_remove <- stringr::str_detect(in_meta[,column], filter_string) %>%
          ifelse(T, in_meta$to_remove)
      }
    }
  }

  # Do subsetting, keeping removed samples for stats.
  meta_removed <- subset(in_meta, subset = to_remove, select = c(-to_remove))
  meta_kept <- subset(in_meta, subset = !to_remove, select = c(-to_remove))

  # Check that sample numbers match
  # From testing, I think this would only happen if there's an NA in the column
  # being filtered on, such that the to_remove column then also has an NA in it.
  # Targets dataframes shouldn't really have NAs in them (not sure if that gets
  # checked earlier in the pipeline), so maybe good to throw an error here in case
  # that happens?
  if (nrow(meta_removed) + nrow(meta_kept) != nrow(in_meta)) {
    cli::cli_abort(c("Issue when subsetting samples",
                     "!" = "Rows kept + rows removed != rows input",
                     "i" = "Is there an {.val NA} in the column you're subsetting from?"))
  }

  # Update metadata samples
  DIAlist$metadata <- meta_kept
  # Update data, removing cols that are no longer present
  DIAlist$data <- DIAlist$data[, DIAlist$metadata$sampleIDs]
  validate_DIAlist(DIAlist)


  # print messages
  cli::cli_inform("Removed {nrow(meta_removed)} of the {nrow(in_meta)} {cli::qty(nrow(in_meta))} sample{?s} in DIAlist")
  cli::cli_inform("{cli::qty(nrow(meta_removed))} Sample{?s} removed:")
  print(meta_removed)
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  DIAlist
}
