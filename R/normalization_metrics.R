
#' Metrics to evaluate normalization methods
#'
#' A set of functions that take in normalized sample data and a grouping factor
#' and calculate some metric of of variability, error, etc., that we can use
#' to evaluate normalization methods. The functions/metrics: \itemize{
#'   \item PCV- Calculated the pooled coefficient of variation (CV, standard
#'     deviation divided by the mean) for all proteins within a group. That is,
#'     the average (mean) of the per-protein CVs.
#'   \item PMAD- Calculate the pooled median absolute deviation. The MAD is the
#'     median of the absolute deviations of each sample from the group median.
#'     The pooled MAD is the group median of the per-protein PMADs.
#'   \item PEV- Calculates the pooled estimate of variance for each group. Uses
#'     the weighted average method to account for unequal sample sizes (see
#'     \url{https://en.wikipedia.org/wiki/Pooled_variance}).
#'   \item COR- Calculates all within-group pairwise correlations in intensity
#'     between samples.
#'   \item log2ratio- Calculates all possible between-group ratios of average
#'     log2-normalized intensity for each protein.
#' }
#'
#' @param data A numeric data frame, where each row is a protein and columns
#'   are densities. For most uses, these are probably normalized data.
#' @param groups A character or factor vector, listing the group(s) the samples
#'   belong to.
#'
#' @return For PCV, PMAD, and PEV: a named vector, with length equal to the
#'   number of groups, giving the metric for each group. For COR, a named vector
#'   giving all within-group pairwise correlations between samples. For log2ratio
#'   with keep_protein_ID = FALSE (the default), a non-named vector giving all
#'   possible between-group log2 ratios for each protein. For log2ratio with
#'   keep_protein_ID = TRUE, a data frame in which rows are proteins, columns are
#'   are the log2ratio for a pairwise comparisons, and rownames give protein ID.
#'
#' @name norm_metrics
#'


#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled coefficient of variation (average of per-protein CVs)
PCV <- function(data, groups) {
  PCV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    CVs <- rowSds(tempData, na.rm = FALSE)/rowMeans(tempData, na.rm = FALSE)
    PCV[group] <- mean(CVs, na.rm = T)
  }
  return(PCV)
}

#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled median absolute deviation (median of per-protein MADs)
PMAD <- function(data, groups) {
  PMAD <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    MAD <- rowMads(tempData, na.rm = F)
    PMAD[group] <- mean(MAD, na.rm = T)
  }
  return(PMAD)
}


#' @rdname norm_metrics
#' @keywords internal
#'
# Pooled estimate of variance
PEV <- function(data, groups) {
  PEV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])

    rowNonNACnt <- rowSums(!is.na(tempData)) - 1
    EV <- rowNonNACnt * rowVars(tempData, na.rm = FALSE)
    PEV[group] <- sum(EV, na.rm = TRUE)/sum(rowNonNACnt, na.rm = TRUE)
  }
  return(PEV)

}

#' @rdname norm_metrics
#' @keywords internal
#'
# Within-group pairwise correlations.
COR <- function(data, groups) {
  COR <- NULL
  for (group in unique(groups)) { # for each group
    corGroup <- NULL
    tempData <- as.matrix(data[,groups %in% group])
    if (ncol(tempData) == 1) { # If only one sample, cor with itself is 1
      corGroup <- 1
    }
    if (ncol(tempData) > 1) { # If multiple sample, get all pw cors
      corVals <- stats::cor(tempData,
                            use = "pairwise.complete.obs",
                            method = "pearson")
      corGroup <- corVals[lower.tri(corVals)]
    }
    COR[[group]] <- corGroup
  }
  return(unlist(COR))
}



#' @rdname norm_metrics
#' @param keep_protein_ID Should protein ID information be retained? Default is FALSE.
#' @keywords internal
#'
# log2ratio

log2ratio <- function(data, groups, keep_protein_ID = F) {

  # calculate mean intensity for each group
  # Original implementation by charity was a nested for loop that ended up calculating group
  # means multiple times
  # I a one-liner with aggregate, but it was slow as hell
  # Switch back to a loop, but now a non-nested one that uses the speedy rowMeans
  mean_prot_by_group <- as.data.frame(matrix(nrow = nrow(data),
                                             ncol = length(unique(groups)),
                                             dimnames = list(rownames(data),
                                                             sort(unique(groups)))))
  # Get group means in a loop
  for (group in sort(unique(groups))) {
    one_group_data <- data[, groups == group]
    mean_prot_by_group[, group] <- rowMeans(one_group_data, na.rm = T)
  }

  # Then, calculate all pairwise differences across cols/groups
  # Assuming these are all log2 normalized, the difference is the ratio
  # uses a custom function, all_pw_diffs (in utils), to apply this.
  # Need simplify = F and some ugly type conversion and transposition
  # To get this in the format we want (rows = proteins, cols = pw comparisons)
  # If simplify = T, cases with only 1 pw comparison are simplified to vectors.
  pw_diffs_per_gene <- as.data.frame(t(as.data.frame(apply(X = mean_prot_by_group, MARGIN = 1, FUN = all_pw_diffs, simplify = F))))


  # The result is a data frame, where each row is a gene and each col is the difference in one
  # pw comparison.
  # For the log2ratio plots, want to collapse protein info and just return
  # a vector like the other plots.
  # For the MD/MA plots, need to maintain protein info in a data frame
  if (keep_protein_ID) {
    output <- pw_diffs_per_gene
  } else {
    output <- unlist(pw_diffs_per_gene, use.names = F)
  }

  output
}

