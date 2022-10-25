


#' Remove protein contaminant rows
#'
#' Not sure if this should be in internal function only, or if
#' it should stay public. Just calls the general function to
#' filter ptoreins based on annotations, but with some presets
#' that we use.
#'
#' @param DIAlist A DIAlist object to be filtered
#'
#' @return A DIAlist with contaminants filtered out
#'
#' @export
#'
#' @examples
#' # No examples yet
#'
filter_proteins_contaminants <- function(DIAlist) {

  # avoid R CMD check
  Protein.Name <- NULL

  out <- filter_proteins_by_annotation(DIAlist, !(stringr::str_detect(Protein.Name, "DECOY"))) %>%
    filter_proteins_by_annotation(!(stringr::str_detect(Protein.Name, "Group of")))


  num_proteins_in <- nrow(DIAlist$annotation)
  num_proteins_out <- nrow(out$annotation)

  cli::cli_inform("{num_proteins_in - num_proteins_out} contaminants removed")
  cli::cli_inform("{num_proteins_out} DIA protein entries retained")

  validate_DIAlist(out)
}


#' Remove proteins by annotations
#'
#' NEED TO REWRITE
#'
#' @param DIAlist A DIAlist object to be filtered
#' @param condition A logical expression indicating which samples to keep.
#'
#' @return A DIAlist with some proteins filtered out
#'
#' @export
#'
#' @examples
#' # No examples yet
#'
filter_proteins_by_annotation <- function(DIAlist, condition) {

  if (!(class(DIAlist) %in% c("DIAlist"))) {
    cli::cli_abort("{.arg DIAlist} must be a DIAlist object")
  }

  if (is.null(DIAlist$annotation)) {
    cli::cli_abort("{.arg DIAlist} does not contain annotation for filtering samples")
  }

  # get input annotation
  in_annot <- DIAlist$annotation

  # Add a col, keep, with the evaluation of the condition expression.
  condition_call <- substitute(condition)
  in_annot <- within(in_annot, {
    keep <- eval(condition_call)
  })

  # Do subseting, keeping removed samples for stats.
  annotation_removed <- subset(in_annot, subset = !keep, select = c(-keep))
  annotation_kept <- subset(in_annot, subset = keep, select = c(-keep))

  # Check that sample numbers match
  if (nrow(annotation_removed) + nrow(annotation_kept) != nrow(in_annot)) {
    cli::cli_abort(c("Issue when filtering proteins",
                     "!" = "Rows kept + rows removed != rows input",
                     "i" = "Is there an {.val NA} in the column you're filtering on?"))
  }

  # Update metadata samples
  DIAlist$annotation <- annotation_kept
  # Update data, removing proteins
  DIAlist$data <- DIAlist$data[rownames(DIAlist$annotation),]

  validate_DIAlist(DIAlist)
}






#' Filter protein data.
#'
#' A subfunction of \code{\link{process_data}}, this function removes samples that
#' are not present in the targets file and filters out proteins/rows that aren't
#' present in enough samples. For filtering, data rows (proteins/peptides) are
#' only retained if they have non-zero intensities in at least min.reps samples
#' and in at least min.grps groups. E.g., with min.reps = 2 and min.groups = 1,
#' a protein must be quantified in two samples within one group to be retained.
#' This function calls one subfunction, \code{\link{make_log}}.
#'
#' @inheritParams process_data
#'
#' @return A list with 6 elements:
#'   \enumerate{
#'     \item "data"- A dataframe of the samples and proteins that were retained
#'       for further analysis
#'     \item "targets" - A dataframe of targets. Should be nearly the same as
#'       the targets dataframe that was input, possibly with new rownames.
#'     \item "nonzeroSamplesPerGroup"- A data frame listing the number of samples
#'       with nonzero intensities for each protein (row) in each group (column).
#'     \item "rm.data"- A dataframe of the proteins/rows filtered out.
#'     \item "param"- A dataframe giving the parameters and arguments used for
#'       filtering
#'     \item "stats"- A dataframe giving statistics on filtering.
#'   }
#' @export
#'
#' @examples
#' # No examples yet
#'
#'

filter_data <- function(data,
                        targets,
                        min.reps = NULL,
                        min.grps = NULL,
                        group_col = "group") {

  # Check args
  ## Make sure targets are in the input data
  if (any(rownames(targets) %notin% colnames(data))) {
    missing_targets <- rownames(targets)[rownames(targets) %notin% colnames(data)]

    cli::cli_abort(c("Some targets are not present in input data.",
                     "x" = "Missing target{?s}: {missing_targets}",
                     "i" = "Check {.code colnames(data)}"))

  }
  # Make sure the group column is present in the targets data frame
  if (group_col %notin% colnames(targets)) {
    cli::cli_abort(c("Column {.arg {group_col}} not found in {.arg targets} data frame"))
  }
  # Make sure that min.reps and min.grps are specified, no more default values.
  if (is.null(min.reps)) {
    cli::cli_abort(c("You must specify a value for {.arg min.reps}",
                     "i" = "A good rule of thumb is to use 2/3 of the size of the sample group"))
  }
  if (is.null(min.grps)) {
    cli::cli_abort(c("You must specify a value for {.arg min.grps}",
                     "i" = "A good rule of thumb is that this should be at least 1",
                     "i" = "such that proteins must have the minimum number of replicates in at least one group to be analyzed"))
  }

  # Start processing.


  # Keep only data columns that match targets in the target file
  # Also reorders so they're the same.
  data_in_targets <- data[, rownames(targets)]
  # TODO: pretty sure this stopifnot check can be removed?
  # The line above should ensure that the below statement is always true, right?
  stopifnot(rownames(targets) == colnames(data_in_targets))

  num_removed_samples  <- ncol(data) - ncol(data_in_targets)
  removed_samples <- colnames(data)[colnames(data) %notin% colnames(data_in_targets)]


  cli::cli_inform("Removing {num_removed_samples} sample{?s} from the data matrix that {?is/are} not listed in the targets")
  cli::cli_inform("{cli::qty(num_removed_samples)} Sample{?s} removed:")
  cli::cli_inform("{.val {removed_samples}}")


  ## change data column names and targets row names
  ## to sample name i.e. sample column in targets
  if ("sample" %in% colnames(targets)) {
    # TODO: again, this check isn't needed, right?
    stopifnot(rownames(targets) == colnames(data_in_targets))
    rownames(targets) <- targets$sample
    colnames(data_in_targets)    <- rownames(targets)
    # TODO: And this one here? If the above two lines of code didn't fail, the below is true by definition?
    stopifnot(colnames(data_in_targets) == rownames(targets))
    cli::cli_inform(cli::col_yellow("\"sample\" column found in the input target file: renamed column names of {.arg data} and row names of {.arg targets} to sample names"))

  }

  ## TODO: do we need this to be a factor? Might be better to just do everything
  ## as a character vector, right?
  ## extract group column as character vector and make a factor.
  group_membership <- make_factor(as.character(targets[, group_col]))
  # And get just the different groups
  groups <- unique(group_membership)


  reps_per_group <- table(group_membership)      ##  number of samples in each group
  n_groups <- length(groups)   ##  number of groups

  ## If the min.reps set by the user exceeds the number of
  ## samples in any group, give an error
  group_meets_cutoff <- reps_per_group >= min.reps
  if (!all(group_meets_cutoff)) {

    # Prep for error message
    groups_below_cutoff <- names(group_meets_cutoff)[!group_meets_cutoff]
    reps_in_groups_below_cutoff <- reps_per_group[!group_meets_cutoff]

    cli::cli_abort(c("!" = "Some groups do not have the minimum number of replicates {.arg min.reps} = {.val {min.reps}}",
                     "x" = "Groups below threshold: {.val {groups_below_cutoff}}",
                     "i" = "Lower the {.arg min.reps} argument so that it is not greater",
                     "i" = "than the number of replicates in the smallest group, {.val {min(reps_per_group)}}"))
  }

  # Check whether enough groups satisfy the min_reps argument to be able
  # to meet the min.grps cutoff
  if (min.grps > sum(group_meets_cutoff)) {

    cli::cli_abort(c("!" = "Not enough groups have sufficient samples",
                     "x" = "Found {.val {cli::no(sum(group_meets_cutoff))}} group{?s} which meet the {.arg min.reps} = {.val {min.reps}} threshold,",
                     "x" = "but {.arg min.grps} = {.val {min.grps}}.",
                     "i" = "Lower either {.arg min.reps} or {.arg min.grps} so that enough groups meet the threshold"))
  }

  cli::cli_inform("Keeping only protein entries with intensity > 0 in at least {.val {min.reps}} sample{?s} {cli::qty(min.reps)} in at least {.val {min.grps}} group{?s} {cli::qty(min.grps)}")

  ## FILTERING
  ## zeros, if they exist, are replaced with NA
  tmpData <- data_in_targets
  tmpData[,][tmpData[,] == 0] <- NA

  ## calc. no samples in each group with intensities > 0
  # Build a results dataframe in advance, empty for now
  nonzero_samples_per_group_and_gene <- as.data.frame(
    matrix(
      data = NA,
      nrow = nrow(tmpData),
      ncol = length(groups),
      dimnames = list(rownames(tmpData),
                      as.character(groups))
    )
  )

  for (i in seq_along(groups)) { # For each group
    # Select only the cols from that group
    one_group_data <- tmpData[, group_membership == groups[i]]

    # count up the the number of samples in that group with
    # non-missing intensities for each protein/row
    # and add that to the results
    nonzero_samples_per_group_and_gene[, i] <- rowSums(!is.na(one_group_data))
  }


  # Then, find out which proteins/rows have the min number of reps in at least enough groups
  protein_meets_threshold <- apply(nonzero_samples_per_group_and_gene, 1, FUN = function(x) {sum(x >= min.reps) >= min.grps })
  # And filter
  kept_proteins <- tmpData[protein_meets_threshold, ]
  removed_proteins <- tmpData[!protein_meets_threshold, ]


  cli::cli_inform("Filtered {.val {nrow(removed_proteins)}} entr{?y/ies} {cli::qty(nrow(removed_proteins))} from the dataset leaving {.val {nrow(kept_proteins)}} entr{?y/ies} {cli::qty(nrow(kept_proteins))} for analysis")


  # Stats and params and logging
  param <- stats <- list()
  stats[["total_input_samples"]] <- ncol(data)
  stats[["num_removed_samples"]] <- num_removed_samples
  stats[["num_kept_samples"]] <- ncol(data_in_targets)
  stats[["total_input_rows"]] = nrow(data_in_targets)
  stats[["num_removed_rows"]] = nrow(removed_proteins)
  stats[["num_kept_rows"]] = nrow(kept_proteins)

  param[["group"]] <- group_col
  param[["groups"]] <- paste(groups, reps_per_group, sep = "=", collapse = ", ")
  param[["min.reps"]] <- min.reps
  param[["min.grps"]] <- min.grps

  logs <- make_log(param, stats, title="FILTERING X in Y", save = TRUE)

  # Output results
  list(data = kept_proteins,
       targets = targets,
       nonzeroSamplesPerGroup = nonzero_samples_per_group_and_gene,
       rm.data = removed_proteins,
       param = logs$param,
       stats = logs$stats)
}


