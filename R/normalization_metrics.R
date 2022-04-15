##------------------------------
##  [02] PCV
##------------------------------
PCV <- function(data, groups) {
  PCV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    CVs <- genefilter::rowSds(tempData, na.rm = FALSE)/rowMeans(tempData, na.rm = FALSE)
    PCV[group] <- mean(CVs, na.rm = T)
  }
  return(PCV)
}


##------------------------------
##  [03] PMAD
PMAD <- function(data, groups) {
  PMAD <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])
    MAD <- matrixStats::rowMads(tempData, na.rm = FALSE)
    PMAD[group] <- mean(MAD, na.rm = T)
  }
  return(PMAD)
}


##------------------------------
##  [04] PEV
##------------------------------
PEV <- function(data, groups) {
  PEV <- NULL
  for (group in unique(groups)) {
    tempData <- as.matrix(data[, groups %in% group])

    rowNonNACnt <- rowSums(!is.na(tempData)) - 1
    EV <- rowNonNACnt * matrixStats::rowVars(tempData, na.rm = FALSE)
    PEV[group] <- sum(EV, na.rm = TRUE)/sum(rowNonNACnt, na.rm = TRUE)
  }
  return(PEV)

} ## PEV


##------------------------------
##  [05] COR
##------------------------------
COR <- function(data, groups) {
  COR <- NULL
  for (group in unique(groups)) {
    corGroup <- NULL
    tempData <- as.matrix(data[,groups %in% group])
    if (ncol(tempData) == 1) {
      corGroup <- 1
    }
    if (ncol(tempData) > 1) {
      corVals <- stats::cor(tempData,
                            use = "pairwise.complete.obs",
                            method = "pearson")
      for (index in seq_len(ncol(corVals) - 1)) {
        corGroup <- c(corGroup, corVals[index, -(seq_len(index)), drop="FALSE"])
      }
    }
    COR[[group]] <- corGroup
  }
  return(unlist(COR))
} ## COR
