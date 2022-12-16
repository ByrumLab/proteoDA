#' Normalize data in a DAList
#'
#' Normalizes data in a DAList, using one of the available normalization options (see Details).
#'
#' Available normalization methods: \itemize{
#'   \item log2- Binary (base2) log transformation.
#'   \item median- Divides each sample by the per-sample median, then multiplies all samples
#'     by the average of the per-sample medians.
#'   \item mean- Divides each sample by the per-sample mean, then multiplies all samples
#'     by the average of the per-sample means.
#'   \item vsn- Variance-stabilizing normalization (vsn),
#'     using the \code{\link[vsn:justvsn]{vsn::justvsn}} function.
#'   \item quantile- Quantile normalization, using the
#'     \code{\link[preprocessCore:normalize.quantiles]{preprocessCore::normalize.quantiles}}
#'      function.
#'   \item cycloess- Cyclic Loess normalization, using the
#'     \code{\link[limma:normalizeCyclicLoess]{limma::normalizeCyclicLoess}}
#'     function.
#'   \item rlr- Global linear regression normalization, inspired by the
#'     \code{NormalyzerDE::performGlobalRLRNormalization}
#'     function, but with a slightly different implementation.
#'   \item gi- Global intensity normalization, inspired by the
#'     \code{NormalyzerDE::globalIntensityNormalization}
#'     function, but with a slightly different implementation. Divides each sample by the
#'     total intensity for each sample, then multiplies all samples by the median of per-sample
#'     intensities.
#'  }
#'
#' @param DAList A DAList containing non-normalized data.
#' @param norm_method A normalization method to use. Options are "log2", "median",
#'   "mean", "vsn", "quantile", "cycloess", "rlr", and "gi". Default is "log2".
#'
#' @return A DAList with normalized data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' norm_data <- normalize_data(DAList, "log2")
#' }
#'
normalize_data <- function(DAList,
                           norm_method = c("log2", "median", "mean", "vsn", "quantile",
                                      "cycloess", "rlr", "gi")) {

  norm_method <- rlang::arg_match(norm_method)

  validate_DAList(DAList)

  if (!is.null(DAList$tags$normalized)) {
    if (DAList$tags$normalized) {
      cli::cli_abort("Data in DAList are already normalized.")
    }
  }

  # Some normalization methods work on raw data, some on
  # log2 transformed data
  if (norm_method %in% c("log2", "vsn", "gi")) {
    normalized_data <- do.call(what = paste0(norm_method, "Norm"),
                               args = list(dat = DAList$data))
  } else if (norm_method %in% c("median", "mean", "quantile", "cycloess", "rlr")) {
    log2_dat <- log2Norm(dat = DAList$data)
    normalized_data <- do.call(what = paste0(norm_method, "Norm"),
                               args = list(logDat = log2_dat))
  } else {
    cli::cli_abort("{.arg {norm_method}} is not a valid normalization method")
  }

  # Updata data
  DAList$data <- normalized_data
  DAList$tags$normalized <- T
  DAList$tags$norm_method <- norm_method

  validate_DAList(DAList)
}


#' Normalization functions
#'
#' A set of functions that normalize data (in columns of a data frame
#' or matrix), according to various methods: \itemize{
#'   \item log2Norm- Binary (base2) log transformation.
#'   \item medianNorm- Divides each sample by the per-sample median, then multiplies all samples
#'     by the average of the per-sample medians.
#'   \item meanNorm- Divides each sample by the per-sample mean, then multiplies all samples
#'     by the average of the per-sample means.
#'   \item vsnNorm- Variance-stabilizing normalization (vsn),
#'     using the \code{\link[vsn:justvsn]{vsn::justvsn}} function.
#'   \item quantileNorm- Quantile normalization, using the
#'     \code{\link[preprocessCore:normalize.quantiles]{preprocessCore::normalize.quantiles}}
#'      function.
#'   \item cycloessNorm- Cyclic Loess normalization, using the
#'     \code{\link[limma:normalizeCyclicLoess]{limma::normalizeCyclicLoess}}
#'     function.
#'   \item rlrNorm- Global linear regression normalization, inspired by the
#'     \code{NormalyzerDE::performGlobalRLRNormalization}
#'     function, but with a slightly different implementation.
#'   \item giNorm- Global intensity normalization, inspired by the
#'     \code{NormalyzerDE::globalIntensityNormalization}
#'     function, but with a slightly different implementation. Divides each sample by the
#'     total intensity for each sample, then multiplies all samples by the median of per-sample
#'     intensities.
#'  }
#'
#'
#' @param dat A numeric data frame, where each row is a protein and columns
#'   are the raw intensities
#' @param logDat A numeric data frame, where each row is a protein and columns
#'   are log2-transformed intensities.
#'
#' @return A matrix of normalized sample intensities
#'
#' @name norm_functions
#'
NULL



#' @rdname norm_functions
#' @keywords internal
#'
log2Norm <- function(dat) {
  logInt <- log2(as.matrix(dat))
  logInt[is.infinite(as.matrix(logInt))] <- NA
  return(as.matrix(logInt))
}

#' @rdname norm_functions
#' @keywords internal
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
#' @keywords internal
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
#' @keywords internal
#'
vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

#' @rdname norm_functions
#' @keywords internal
#'
quantileNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy = TRUE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

#' @rdname norm_functions
#' @keywords internal
#'
cycloessNorm <- function(logDat) {
  cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method = "fast")
  colnames(cycLoessNormed) <- colnames(logDat)
  row.names(cycLoessNormed) <- rownames(logDat)
  return(as.matrix(cycLoessNormed))
}

#' @rdname norm_functions
#' @keywords internal
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
#' @keywords internal
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

