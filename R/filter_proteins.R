#' Remove protein contaminant rows
#'
#' Not sure if this should be a UAMS function only, or if
#' it should stay public. Just calls the general function to
#' filter proteins based on annotations, but with some presets
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

  cli::cli_inform("{.val {num_proteins_in - num_proteins_out}} contaminants removed")
  cli::cli_inform("{.val {num_proteins_out}} DIA protein entries retained")

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
  # Add tags to track filtering
  DIAlist$tags$filter_protein_by_annotation <- c(DIAlist$tags$filter_protein_by_annotation, list(condition = condition_call))

  validate_DIAlist(DIAlist)
}


#' Filter protein data.
#'
#' NEED TO REWRITE
#'
#' @param DIAlist A DIAlist to be filtered
#' @param min_reps The minimum number of replicates/samples within a group
#'   that need to have non-zero intensity for a given protein/peptide for that
#'   peptide to be considered as quantified within a group.
#' @param min_groups The minimum number of groups that must have at
#'   least min_reps non-zero samples for a given protein/peptide to be retained.
#' @param group_col Optional. The name of the column in the targets dataframe which
#'   lists the experimental group names. Default is "group"
#'
#' @return A filtered DIAlist
#' @export
#'
#' @examples
#' # No examples yet
#'
#'
filter_proteins_by_group <- function(DIAlist,
                                      min_reps = NULL,
                                      min_groups = NULL,
                                      group_col = "group") {

  validate_DIAlist(DIAlist)

  if (is.null(DIAlist$metadata)) {
    cli::cli_abort(c("The {.arg DIAlist} object must include metadata to filter proteins by group"))
  }

  # Make sure the group column is present in the metadata
  if (group_col %notin% colnames(DIAlist$metadata)) {
    cli::cli_abort(c("Column {.arg {group_col}} not found in metadata slot of {.arg DIAlist}"))
  }
  # Make sure that min_reps and min_groups are specified, no more default values.
  if (is.null(min_reps)) {
    cli::cli_abort(c("You must specify a value for {.arg min_reps}",
                     "i" = "A good rule of thumb is to use 2/3 of the size of the sample group"))
  }
  if (is.null(min_groups)) {
    cli::cli_abort(c("You must specify a value for {.arg min_groups}",
                     "i" = "A good rule of thumb is that this should be at least 1",
                     "i" = "such that proteins must have the minimum number of replicates in at least one group to be analyzed"))
  }

  ## extract group column as character vector
  group_membership <- as.character(DIAlist$metadata[, group_col])
  # And get just the different groups
  groups <- unique(group_membership)

  reps_per_group <- table(group_membership)
  n_groups <- length(groups)

  ## If the min_reps set by the user exceeds the number of
  ## samples in any group, give an error
  group_meets_cutoff <- reps_per_group >= min_reps
  if (!all(group_meets_cutoff)) {
    # Prep for error message
    groups_below_cutoff <- names(group_meets_cutoff)[!group_meets_cutoff]
    reps_in_groups_below_cutoff <- reps_per_group[!group_meets_cutoff]

    cli::cli_abort(c("!" = "Some groups do not have the minimum number of replicates {.arg min_reps} = {.val {min_reps}}",
                     "x" = "Groups below threshold: {.val {groups_below_cutoff}}",
                     "i" = "Lower the {.arg min_reps} argument so that it is not greater",
                     "i" = "than the number of replicates in the smallest group, {.val {min(reps_per_group)}}"))
  }

  # Check whether enough groups satisfy the min_reps argument to be able
  # to meet the min_groups cutoff
  if (min_groups > sum(group_meets_cutoff)) {

    cli::cli_abort(c("!" = "Not enough groups have sufficient samples",
                     "x" = "Found {.val {cli::no(sum(group_meets_cutoff))}} group{?s} which meet the {.arg min_reps} = {.val {min_reps}} threshold,",
                     "x" = "but {.arg min_groups} = {.val {min_groups}}.",
                     "i" = "Lower either {.arg min_reps} or {.arg min_groups} so that enough groups meet the threshold"))
  }

  cli::cli_inform("Keeping only protein entries with intensity > 0 in at least {.val {min_reps}} sample{?s} {cli::qty(min_reps)} in at least {.val {min_groups}} group{?s} {cli::qty(min_groups)}")

  ## FILTERING
  ## zeros, if they exist, are replaced with NA
  tmpData <- DIAlist$data
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
  protein_meets_threshold <- apply(nonzero_samples_per_group_and_gene, 1, FUN = function(x) {sum(x >= min_reps) >= min_groups })
  # And filter
  kept_proteins <- tmpData[protein_meets_threshold, ]
  removed_proteins <- tmpData[!protein_meets_threshold, ]

  cli::cli_inform("Filtered {.val {nrow(removed_proteins)}} entr{?y/ies} {cli::qty(nrow(removed_proteins))} from the dataset leaving {.val {nrow(kept_proteins)}} entr{?y/ies} {cli::qty(nrow(kept_proteins))} for analysis")

  out <- DIAlist
  # Update data and annotation
  out$data <- kept_proteins
  out$annotation <- out$annotation[rownames(out$data),]
  # add tags to track filtering
  out$tags$filter_protein_by_group <- c(out$tags$filter_protein_by_group, list(min_reps = min_reps, min_groups = min_groups, group_col = group_col))

  validate_DIAlist(out)
}
