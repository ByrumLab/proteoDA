#' Filter samples from a DIAlist
#'
#' NEED TO REWRITE
#'
#' @param DIAlist The DIAlist to be filtered
#' @param condition A logical expression indicating which samples to keep.
#'
#' @return A prototype of our new S3 list type.
#'
#' @export
#'
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' DIAlist <- filter_samples(DIAlist, group != "Pool")
#' }
#'

filter_samples <- function(DIAlist, condition) {

  if (!(class(DIAlist) %in% c("DIAlist"))) {
    cli::cli_abort("{.arg DIAlist} must be a DIAlist object")
  }

  if (is.null(DIAlist$metadata)) {
    cli::cli_abort("{.arg DIAlist} does not contain metadata for filtering samples")
  }

  cli::cli_rule()

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
