#' Filter samples from targets file
#'
#' Removes unwanted samples (rows) from the targets data frame This function is
#' basically \code{\link[base]{subset}} with extra functionality. It removes
#' rows if they contain certain values in a given column. By default, ignores
#' the case of strings in rm.val. Calls \code{\link{make_log}} as a subfunction.
#'
#' @param targets The input targets dataframe to be filtered
#' @param filter_list A named list describing how to filter targets dataframe.
#'   The name for each element of the list gives a column to be filters, and
#'   each element of the list should be a character vector of strings to search for
#'   in the given column. See examples.
#' @param ignore.case Should the case of the strings in rm.vals be ignored? Is
#'   TRUE by default, such that case does not have to match. If FALSE, case must
#'   match for rows to be removed.
#'
#' @return A list with three elements:
#'   \enumerate{
#'     \item "targets", a dataframe of subsetted targets
#'     \item "stats", a dataframe with statistics on filtering
#'     \item "param", a dataframe giving the parameters used for filtering
#'   }
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' sub_targets <- subset_targets(targets = target,
#'                               filter_list = list(group = "pool"
#'                                                  sample = c("sampleA", "sampleB")))
#' }
#'

subset_targets <- function(targets, filter_list, ignore.case = TRUE) {

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
  targets$to_remove <- F
  for (column in names(filter_list)) {
    for (filter_string in unique(filter_list[[column]])) {
      if (ignore.case) {
        targets$to_remove <- stringr::str_to_lower(targets[,column]) %>%
          stringr::str_detect(stringr::str_to_lower(filter_string)) %>%
          ifelse(T, targets$to_remove)
      } else {
        targets$to_remove <- stringr::str_detect(targets[,column], filter_string) %>%
          ifelse(T, targets$to_remove)
      }
    }
  }

  # Do subsetting, keeping removed samples for stats.
  tar_removed <- subset(targets, subset = to_remove, select = c(-to_remove))
  tar_kept <- subset(targets, subset = !to_remove, select = c(-to_remove))

  # Check that sample numbers match
  # From testing, I think this would only happen if there's an NA in the column
  # being filtered on, such that the to_remove column then also has an NA in it.
  # Targets dataframes shouldn't really have NAs in them (not sure if that gets
  # checked earlier in the pipeline), so maybe good to throw an error here in case
  # that happens?
  if (nrow(tar_removed) + nrow(tar_kept) != nrow(targets)) {
    cli::cli_abort("Issue when subsetting targets",
                   "!" = "Rows kept + rows removed != rows input",
                   "i" = "Is there an {.val NA} in the column you're subsetting from?")
  }

  # Write info to logs:
  param <- stats <- list()
  # Param
  param[["filter_list"]] <- paste(lapply(names(filter_list), function (x) paste(x, "=", paste(filter_list[[x]], collapse = " or "))), collapse = "; ")
  param[["ignore.case"]] <- as.character(ignore.case)
  # Stats
  stats[["input_samples"]] <- nrow(targets)
  stats[["samples_removed"]] <- nrow(tar_removed)
  stats[["samples_kept"]] <- nrow(tar_kept)
  # Logfile
  title <- paste0("SUBSET TARGETS STATS (", paste(names(filter_list), collapse = ", "), ")")
  logs <- make_log(param = param, stats = stats, title = "test", save=TRUE)

  # print messages
  cli::cli_inform("Removed {nrow(tar_removed)} of the {nrow(targets)} {cli::qty(nrow(targets))} sample{?s}  from the targets file")
  cli::cli_inform("{cli::qty(nrow(tar_removed))} Sample{?s} removed:")
  print(tar_removed)
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  # Return
  list(targets = tar_kept, stats = logs$stats, param = logs$param)
}
