#' Filter samples from targets file
#'
#' Removes unwanted samples (rows) from the targets data frame This function is
#' basically \code{\link[base]{subset}} with extra functionality. It removes
#' rows if they contain certain values in a given column. Currently, can only
#' filter on one column at a time, but can provide multiple values to filter
#' against. By default, ignores the case of strings in rm.val. Calls
#' \code{\link{make_log}} as a subfunction.
#'
#' @param targets The input targets dataframe to be filtered
#' @param filter_column The column in which to check values.
#' @param rm.vals A character vector giving the values to search for in
#'   filter_column. Rows matching these values will be fitlered out.
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
#' @examples
#' # No examples yet
#'
subset_targets <- function(targets, filter_column, rm.vals, ignore.case = TRUE) {

  # Check args
  if (is.null(filter_column) | anyNA(filter_column)) {
    cli::cli_abort(c("{.arg filter_column} cannot be {.val {filter_column}}",
                     "i" = "Did you forget to specify the {.arg filter_column} argument?"))
  }
  if (length(filter_column) != 1) {
    cli::cli_abort(c("{.arg filter_column} must be of length 1, not {length(filter_column)}}"))
  }

  cli::cli_rule()

  # Set up some vars to avoid R CMD check warning
  to_remove <- for_match <- remove_reason <- NULL

  ## Set up removal patterns
  if (ignore.case) {
    removal_patterns <- stringr::str_to_lower(rm.vals)
    targets$for_match <- stringr::str_to_lower(targets[,filter_column])
  } else {
    removal_patterns <- rm.vals
    targets$for_match <- targets[,filter_column]
  }

  # Mark for removal any rows that match/contain the supplied strings to exclude
  targets$to_remove <- F
  targets$remove_reason <- NA
  for (pattern in removal_patterns) {
    targets$to_remove <- ifelse(stringr::str_detect(targets$for_match, pattern), T, targets$to_remove)
    targets$remove_reason <- ifelse(stringr::str_detect(targets$for_match, pattern), pattern, targets$remove_reason)
  }

  # Do subsetting, keeping removed samples for stats.
  tar_removed <- subset(targets, subset = to_remove, select = c(-for_match, -to_remove))
  tar_kept <- subset(targets, subset = !to_remove, select = c(-for_match, -to_remove, -remove_reason))

  # Check that sample numbers match
  # TODO: Not sure we even need this error message? I think it'd only happen with an NA
  # in the column being subsetted, which shouldnt happen and which would likely trigger
  # an error higher up??
  if (nrow(tar_removed) + nrow(tar_kept) != nrow(targets)) {
    cli::cli_abort("Issue when subsetting targets",
                   "!" = "Rows kept + rows removed != rows input",
                   "i" = "Is there an {.val NA} in the column you're subsetting from?")
  }

  # Write info to logs:
  param <- stats <- list()
  # Param
  param[["filter_column"]] <- filter_column
  param[["rm.vals"]] <- paste(unique(tar_removed$remove_reason), collapse=", ")
  param[["ignore.case"]] <- as.character(ignore.case)
  # Stats
  stats[["input_samples"]] <- nrow(targets)
  for (reason in unique(tar_removed$remove_reason)) {
    stats[[paste0(filter_column, " matched ", reason)]] <- nrow(subset(tar_removed, remove_reason == reason))
  }
  stats[["samples_removed"]] <- nrow(tar_removed)
  stats[["samples_kept"]] <- nrow(tar_kept)
  # Logfile
  title <- paste0("SUBSET TARGETS STATS (", filter_column, ")")
  logs <- make_log(param = param, stats = stats, title = title, save=TRUE)

  # print messages
  cli::cli_inform("Removed {nrow(tar_removed)} of the {nrow(targets)} sample{?s} {cli::qty(nrow(targets))} from the targets file")
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))


  # Return
  list(targets = tar_kept, stats = logs$stats, param = logs$param)
}
