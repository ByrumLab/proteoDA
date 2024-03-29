#' Remove proteins based on annotation data
#'
#' This function is used to remove proteins from a DAList, filtering using data
#' in the annotation data frame of the DAList. Proteins which do not produce a
#' value of TRUE for the supplied condition are removed from both the data and
#' annotation slots of the DAList. If condition evaluates to NA, the function
#' will return an error.
#'
#' @param DAList A DAList object to be filtered.
#' @param condition An expression that returns a logical value, defined in terms
#'   of variables present in the annotation data frame of the supplied DAList.
#'   Proteins are kept if the condition is TRUE for that protein.
#'
#' @return A DAList, with proteins that do not meet the condition removed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose the DAList$annotation data frame contains three columns:
#' # uniprot_id = An alpha-numeric ID uniquely identifying a protein
#' # protein_name = A character giving the protein name
#' # molecular_weight = A numeric value, giving the protein's MW in kDa.
#'
#' # Remove a specific protein by ID
#' filtered <- filter_proteins_by_annotation(DAList,
#'                                           uniprot_id != "abc123")
#'
#' # Remove any protein which contains "keratin" in the name:
#' filtered <- filter_proteins_by_annotation(DAList,
#'                                           !grepl(pattern = "keratin",
#'                                                  x = protein_name))
#'
#' # Remove any protein with molecular weight < 30 kDa
#' filtered <- filter_proteins_by_annotation(DAList,
#'                                           molecular_weight > 30)
#'
#'
#' # Filtering functions can be chained together
#' filtered <- DAList |>
#'   filter_proteins_by_annotation(!grepl(pattern = "keratin",
#'                                        x = protein_name)) |>
#'   filter_proteins_by_annotation(molecular_weight > 30)
#' }
#'
filter_proteins_by_annotation <- function(DAList, condition) {

  if (!(class(DAList) %in% c("DAList"))) {
    cli::cli_abort("{.arg DAList} must be a DAList object")
  }

  if (is.null(DAList$annotation)) {
    cli::cli_abort("{.arg DAList} does not contain annotation for filtering samples")
  }

  # get input annotation
  in_annot <- DAList$annotation

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
  DAList$annotation <- annotation_kept
  # Update data, removing proteins
  DAList$data <- DAList$data[rownames(DAList$annotation),]
  # Add tags to track filtering
  DAList$tags$filter_proteins_by_annotation <- c(DAList$tags$filter_proteins_by_annotation, list(condition = condition_call))

  validate_DAList(DAList)
}


#' Filter protein data by number of quantified samples in a group
#'
#' This function is used to remove proteins from a DAList, filtering out proteins
#' based on levels of missing values in the "data" data frame of the DAList. The
#' grouping_column must be a column in the metadata of the DAList, which lists the
#' group membership for each sample. The min_reps and min_groups arguments that determine
#' the number of replicates/samples per group (min_reps) and number of groups
#' (min_groups) in which a protein must have non-missing intensity values
#' in order to be retained. This function assumes that all missing values are encoded
#' as NA. See \code{\link{zero_to_missing}} and \code{\link{missing_to_zero}} for helper
#' functions to convert missing values to and from 0.
#'
#' @param DAList A DAList object to be filtered.
#' @param min_reps The minimum number of replicates/samples within a group
#'   that need to have a non-missing intensity value for a given protein in order for that
#'   protein to be considered as quantified within a group.
#' @param min_groups The minimum number of groups that must have at
#'   least min_reps non-zero samples for a given protein to be retained.
#' @param grouping_column The name of the column in the metadata which provides
#'   the group membership for each sample. Default is "group".
#'
#' @return A DAList, with proteins that are not present in sufficient samples
#'   and groups removed.
#' @export
#'
#' @examples
#'
#' \dontrun{
#'   # Suppose the DAList contains data from 20 samples across 4
#'   # experimental groups (5 samples per group), with the group membership
#'   # listed in a column named "group"
#'
#'   # Strict filtering:
#'   # no missing data
#'   # Proteins must be present in all samples in all groups
#'   filtered <- filter_proteins_by_group(DAList,
#'                                        min_reps = 5,
#'                                        min_groups = 4,
#'                                        grouping_column = "group")
#'   # Lax filtering:
#'   # protein must be present in at least one sample in each group
#'   filtered <- filter_proteins_by_group(DAList,
#'                                        min_reps = 1,
#'                                        min_groups = 4,
#'                                        grouping_column = "group")
#'
#'  # Filtering functions can be chained together
#'  filtered <- DAList |>
#'    filter_proteins_by_annotation(!grepl(pattern = "keratin",
#'                                         x = protein_name)) |>
#'    filter_proteins_by_group(min_reps = 1,
#'                             min_groups = 4,
#'                             grouping_column = "group")
#' }
#'
filter_proteins_by_group <- function(DAList,
                                     min_reps = NULL,
                                     min_groups = NULL,
                                     grouping_column = "group") {

  validate_DAList(DAList)

  if (is.null(DAList$metadata)) {
    cli::cli_abort(c("The {.arg DAList} object must include metadata to filter proteins by group"))
  }

  # Make sure the group column is present in the metadata
  if (grouping_column %notin% colnames(DAList$metadata)) {
    cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata slot of {.arg DAList}"))
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
  group_membership <- as.character(DAList$metadata[, grouping_column])
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

  cli::cli_inform("Keeping only protein entries with non-missing intensity in at least {.val {min_reps}} sample{?s} {cli::qty(min_reps)} in at least {.val {min_groups}} group{?s} {cli::qty(min_groups)}")

  ## FILTERING
  tmpData <- DAList$data

  ## calc. no samples in each group with non-missing intensity
  # Build a results data frame in advance, empty for now
  nonmissing_samples_per_group_and_gene <- as.data.frame(
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
    one_group_data <- tmpData[, group_membership == groups[i], drop = F]

    # count up the the number of samples in that group with
    # non-missing intensities for each protein/row
    # and add that to the results
    nonmissing_samples_per_group_and_gene[, i] <- rowSums(!is.na(one_group_data))
  }


  # Then, find out which proteins/rows have the min number of reps in at least enough groups
  protein_meets_threshold <- apply(nonmissing_samples_per_group_and_gene, 1, FUN = function(x) {sum(x >= min_reps) >= min_groups })
  # And filter
  kept_proteins <- tmpData[protein_meets_threshold, , drop = F]
  removed_proteins <- tmpData[!protein_meets_threshold, , drop = F]

  cli::cli_inform("Filtered {.val {nrow(removed_proteins)}} entr{?y/ies} {cli::qty(nrow(removed_proteins))} from the dataset leaving {.val {nrow(kept_proteins)}} entr{?y/ies} {cli::qty(nrow(kept_proteins))} for analysis")

  out <- DAList
  # Update data and annotation
  out$data <- kept_proteins
  out$annotation <- out$annotation[rownames(out$data), , drop = F]
  # add tags to track filtering
  out$tags$filter_proteins_by_group <- c(out$tags$filter_proteins_by_group, list(min_reps = min_reps, min_groups = min_groups, grouping_column = grouping_column))

  validate_DAList(out)
}


#' Filter protein data by proportion of quantified samples in group
#'
#' This function is used to remove proteins from a DAList, filtering out proteins
#' based on levels of missing data in the "data" data frame of the DAList. The
#' grouping_column must be a column in the metadata of the DAList which lists the
#' group membership for each sample. Proteins must have a non-missing intensity
#' value in at least min_prop of samples within each group in order to be retained.
#' When min_prop leads to a non-integer value for a given group, it is rounded up:
#' e.g., with 10 samples in a group and a min_prop of 0.75, a protein must be present in
#' at least 8 samples to be retained. This function assumes that all missing values are encoded
#' as NA. See \code{\link{zero_to_missing}} and \code{\link{missing_to_zero}} for helper
#' functions to convert missing values to and from 0.
#'
#' \code{filter_proteins_by_proportion()} is useful when sample sizes are not constant
#' across groups, allowing similar levels of missingness across groups even when sample sizes
#' are imbalanced.
#'
#' @inheritParams  filter_proteins_by_group
#' @param min_prop The minimum proportion of samples within a group
#'   in which a protein must be found for it to be retained.
#'   Must be a numeric value from 0 to 1.
#'
#' @return A DAList, with proteins that are not present in sufficient samples
#'   removed.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Suppose the DAList contains data from 40 samples across 2
#'   # experimental groups, with 15 samples in one group and 25 in the other.
#'   # and group membership listed in a column names "group"
#'
#'   # Strict filtering:
#'   # no missing data
#'   # Proteins must be present in all samples in all groups
#'   filtered <- filter_proteins_by_proportion(DAList,
#'                                             min_prop = 1,
#'                                             grouping_column = "group")
#'   # moderate filtering:
#'   # protein must be present in at 50% of samples within each group.
#'   # That is, 8 samples in the group of 15 and 13 samples in the group of 25.
#'   filtered <- filter_proteins_by_proportion(DAList,
#'                                             min_prop = 0.5,
#'                                             grouping_column = "group")
#'
#'   # Filtering functions can be chained together
#'   filtered <- DAList |>
#'     filter_proteins_by_annotation(!grepl(pattern = "keratin",
#'                                          x = protein_name)) |>
#'     filter_proteins_by_proportion(min_prop = 0.5,
#'                                   grouping_column = "group")
#' }
#'
#'
filter_proteins_by_proportion <- function(DAList,
                                          min_prop = NULL,
                                          grouping_column = "group") {
  validate_DAList(DAList)

  if (is.null(DAList$metadata)) {
    cli::cli_abort(c("The {.arg DAList} object must include metadata to filter proteins by group"))
  }

  # Make sure the group column is present in the metadata
  if (grouping_column %notin% colnames(DAList$metadata)) {
    cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata slot of {.arg DAList}"))
  }
  # Make sure that min_reps and min_groups are specified, no more default values.
  if (is.null(min_prop)) {
    cli::cli_abort(c("You must specify a value for {.arg min_prop}",
                     "i" = "A good rule of thumb is to use 0.66"))
  }

  if (min_prop < 0 | min_prop > 1) {
    cli::cli_abort(c("{.arg min_prop} must be from 0 and 1, not {.val {min_prop}}"))
  }

  cli::cli_inform("Keeping only protein entries with non-missing intensity in at least {.val {min_prop*100}}% of samples in each group")

  tmpData <- DAList$data

  # Prep for filtering by getting a vector of group memberships per sample
  # and calculating the threshold for each group
  group_membership <- as.character(DAList$metadata[, grouping_column])
  group_thresholds <-  ceiling(table(as.character(DAList$metadata[, grouping_column]))*min_prop)

  # Build a results data frame in advance, empty for now
  protein_passes_threshold_per_group <- as.data.frame(
    matrix(
      data = NA,
      nrow = nrow(tmpData),
      ncol = length(names(group_thresholds)),
      dimnames = list(rownames(tmpData),
                      as.character(names(group_thresholds)))
    )
  )
  # Loop over group thresholds, check if protein passes for that group
  for (group in names(group_thresholds)) {
    one_group_data <- tmpData[, group_membership == group]

    # count up the the number of samples in that group with
    # non-missing intensities for each protein/row
    # and check whether it is over the threshold
    protein_passes_threshold_per_group[, group] <- rowSums(!is.na(one_group_data)) >= group_thresholds[group]
  }

  # Keep only proteins where all groups are T
  protein_meets_threshold <- apply(X = protein_passes_threshold_per_group, MARGIN = 1, FUN = all)
  # And filter
  kept_proteins <- tmpData[protein_meets_threshold, ]
  removed_proteins <- tmpData[!protein_meets_threshold, ]

  cli::cli_inform("Filtered {.val {nrow(removed_proteins)}} entr{?y/ies} {cli::qty(nrow(removed_proteins))} from the dataset leaving {.val {nrow(kept_proteins)}} entr{?y/ies} {cli::qty(nrow(kept_proteins))} for analysis")

  out <- DAList
  # Update data and annotation
  out$data <- kept_proteins
  out$annotation <- out$annotation[rownames(out$data),, drop = F]
  # add tags to track filtering
  out$tags$filter_proteins_by_proportion <- c(out$tags$filter_proteins_by_proportion, list(min_prop = min_prop, grouping_column = grouping_column))

  validate_DAList(out)
}

