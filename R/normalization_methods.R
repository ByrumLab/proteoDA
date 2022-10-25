
#' Apply all normalizations to a dataframe
#'
#' Takes in raw intensity data and applies 8 different normalization methods,
#' returning a list of the data normalized with each method. Used internally in
#' \code{\link{write_proteinorm_report}}.
#'
#' @param data A dataframe of raw data to be normalized. Rows are proteins and
#'   columns are raw intensity data.
#'
#' @return A list length 8, where each item in the list is a
#'       named dataframe. Names give the normalization method, and the dataframe
#'       gives the normalized intensity data
#'
#' @examples
#' # No examples yet
#'
apply_all_normalizations <- function(data) {


  normList <- NULL

  ## apply normalization methods using functions listed above.
  ## NOTE: most of the normalizations use log2(intensity) as input except
  ## VSN which normalizes using raw intensities with no log2 transformation.
  normList[["log2"]]     <- logNorm(dat = data)
  normList[["median"]]   <- medianNorm(logDat = normList[["log2"]])
  normList[["mean"]]     <- meanNorm(logDat = normList[["log2"]])
  normList[["vsn"]]      <- vsnNorm(dat = data)
  normList[["quantile"]] <- quantNorm(logDat = normList[["log2"]])
  normList[["cycloess"]] <- cycLoessNorm(logDat = normList[["log2"]])
  normList[["rlr"]]      <- rlrNorm(logDat = normList[["log2"]])
  normList[["gi"]]       <- giNorm(dat = data)

  normList
}



#' Normalization functions
#'
#' A set of functions that normalize sample data (in columns of a data frame
#' or matrix), according to various methods: \itemize{
#'   \item logNorm- A binary (base2) log transformation
#'   \item medianNorm- Divides by per-sample median, then multiplies by the
#'     average of the per-sample medians.
#'   \item meanNorm- Divides by per-sample mean, then multiplies by the
#'     average of the per-sample means
#'   \item vsnNorm-  Variance-stabilizing normalization (vsn),
#'     using the \code{\link[vsn:justvsn]{vsn::justvsn}} function.
#'   \item quantNorm-  Quantile normalization, using the
#'     \code{\link[preprocessCore:normalize.quantiles]{preprocessCore::normalize.quantiles}}
#'      function.
#'   \item cycLoessNorm- Cyclic Loess normalization, using the
#'     \code{\link[limma:normalizeCyclicLoess]{limma::normalizeCyclicLoess}}
#'     function.
#'   \item rlrNorm- Global linear regression normalization, inspired by the
#'     \code{NormalyzerDE::performGlobalRLRNormalization}
#'     function, but with a slightly different implementation.
#'   \item giNorm- Global intensity normalization, inspired by the
#'     \code{NormalyzerDE::globalIntensityNormalization}
#'     function, but with a slightly different implementation. In brief,
#'     normalizes each sample/column by dividing by the total intensity for the
#'     sample and multiplying by the median of per-sample intensities. Thus,
#'     like mean or median normalization, but using sum instead of
#' }
#'
#'
#' @param dat A numeric data frame, where each row is a protein and columns
#'   are the raw, unstandardized intensities
#' @param logDat A numeric data frame, where each row is a protein and columns
#'   are log2-transformed intensities.
#'
#' @return A matrix of normalized sample intensities
#'
#' @name norm_functions
#'
#' @examples
#' # No examples yet
#'
NULL
#> NULL



#' @rdname norm_functions
#' @export
#'
logNorm <- function(dat) {
  logInt <- log2(as.matrix(dat))
  logInt[is.infinite(as.matrix(logInt))] <- NA
  return(as.matrix(logInt))
}

#' @rdname norm_functions
#' @export
#'
medianNorm <- function(logDat) {
  # Find medians of each sample
  # Divide by median
  # Multiply by mean of medians
  sampleMed <- apply(logDat, 2, stats::median, na.rm = TRUE)
  meanMed <- mean(sampleMed, na.rm = TRUE)
  out <- t(t(logDat) / sampleMed)
  out <- out * meanMed
  return(as.matrix(out))
}

#' @rdname norm_functions
#' @export
#'
meanNorm <- function(logDat) {
  # Find means of each sample
  # Divide by mean
  # Multiply by mean of means
  sampleMean <- apply(logDat, 2, mean, na.rm = TRUE)
  meanMean <- mean(sampleMean, na.rm = TRUE)
  out <- t(t(logDat) / sampleMean)
  out <- out * meanMean
  return(as.matrix(out))
}

#' @rdname norm_functions
#' @export
#'
vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

#' @rdname norm_functions
#' @export
#'
quantNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy = TRUE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

#' @rdname norm_functions
#' @export
#'
cycLoessNorm <- function(logDat) {
  cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method = "fast")
  colnames(cycLoessNormed) <- colnames(logDat)
  row.names(cycLoessNormed) <- rownames(logDat)
  return(as.matrix(cycLoessNormed))
}

#' @rdname norm_functions
#' @export
#'
rlrNorm <- function(logDat) {
  # First, an internal function to do get the coefficients from rlm regression
  # for one column
  get_rlm_coeffs_and_slopes <- function(unnormalized_col, predictors) {
    MASS::rlm(unnormalized_col ~ predictors,
              na.action = stats::na.exclude)$coefficients
  }

  # Median intensities across samples are the predictors for the rlm
  row_medians <- rowMedians(logDat, na.rm = TRUE)

  # Apply the one-col function across the input matrix cols
  coefficients <- t(apply(X = logDat,
                        MARGIN = 2,
                        FUN = get_rlm_coeffs_and_slopes,
                        predictors = row_medians))
  stopifnot(nrow(t(logDat)) == nrow(coefficients))
  # Do the normalization: subtract the intercept from each point, then divide
  # by slope
  normalized_mat <- t((t(logDat) - coefficients[,1])/coefficients[,2])
  normalized_mat
}

#' @rdname norm_functions
#' @export
#'
giNorm <- function(dat) {
  # Make sure input data is matrix
  dat_mat <- as.matrix(dat)
  # Follow same pattern as mean and median, but with
  # sums
  col_sums <- colSums(dat_mat, na.rm = TRUE)
  median_col_sum <- stats::median(col_sums)

  normMatrix <- t(t(dat_mat)/col_sums) * median_col_sum
  # extra step to return log2 values, as with previous implementation
  normLog2Matrix <- log2(normMatrix)
  normLog2Matrix[is.infinite(normLog2Matrix)] <- NA
  normLog2Matrix
}

