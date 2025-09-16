#' Normalize data in a DAList
#'
#' Normalizes data in a DAList using one of several methods (see Details).
#'
#' **Available methods**
#' \itemize{
#'   \item \code{"log2"} — Binary (base-2) log transformation of raw intensities.
#'   \item \code{"median"} — Divide each sample by its median; rescale by the mean of sample medians.
#'   \item \code{"mean"} — Divide each sample by its mean; rescale by the mean of sample means.
#'   \item \code{"vsn"} — Variance-stabilizing normalization via \code{\link[vsn:justvsn]{vsn::justvsn}} (expects raw scale).
#'   \item \code{"quantile"} — Quantile normalization via \code{\link[preprocessCore:normalize.quantiles]{preprocessCore::normalize.quantiles}}.
#'   \item \code{"cycloess"} — Cyclic loess via \code{\link[limma:normalizeCyclicLoess]{limma::normalizeCyclicLoess}} (can be applied within groups).
#'   \item \code{"rlr"} — Robust linear regression (global) normalization.
#'   \item \code{"gi"} — Global intensity normalization (expects raw scale; returns log2).
#' }
#'
#' **Log scale behavior**
#' Methods \code{"median"}, \code{"mean"}, \code{"quantile"}, \code{"cycloess"}, and \code{"rlr"}
#' operate on log2 data. By default, this function will apply \code{log2} internally before those
#' methods. If your \code{DAList$data} are \emph{already} log2-transformed (e.g., after Perseus-style
#' imputation), set \code{input_is_log2 = TRUE} to avoid double logging.
#'
#' @param DAList A \code{DAList} containing non-normalized data (matrix/data frame in \code{DAList$data}).
#' @param norm_method Normalization method. One of \code{"log2"}, \code{"median"}, \code{"mean"},
#'   \code{"vsn"}, \code{"quantile"}, \code{"cycloess"}, \code{"rlr"}, or \code{"gi"}.
#' @param input_is_log2 Logical. If \code{TRUE}, indicates that \code{DAList$data} are already on the
#'   log2 scale for methods that expect log2 input (\code{"median"}, \code{"mean"}, \code{"quantile"},
#'   \code{"cycloess"}, \code{"rlr"}). Default \code{FALSE}.
#' @param ... Additional arguments forwarded to the underlying normalization function. For example,
#'   pass \code{groups=} to \code{cycloessNorm()} to normalize independently within groups.
#'
#' @return A \code{DAList} with \code{data} replaced by normalized values and tags updated.
#'
#' @examples
#' \dontrun{
#' # Example: log2 -> impute -> within-group cyclic loess (avoid double log)
#' DAList$data <- log2Norm(DAList$data)
#' DAList$data <- perseus_impute(DAList$data, seed = 1)
#' DAList <- normalize_data(
#'   DAList,
#'   norm_method   = "cycloess",
#'   input_is_log2 = TRUE,                 # <-- prevents double log2
#'   groups        = DAList$sample_metadata$group
#' )
#'
#' # Global cyclic loess on raw data (function will log2 internally)
#' DAList <- normalize_data(DAList, norm_method = "cycloess")
#'
#' # VSN works on raw (non-log) data
#' DAList <- normalize_data(DAList, norm_method = "vsn")
#' }
#'
#' @export
normalize_data <- function(DAList,
                           norm_method = c("log2","median","mean","vsn","quantile","cycloess","rlr","gi"),
                           input_is_log2 = FALSE,   # <— NEW
                           ...) {
  
  norm_method <- rlang::arg_match(norm_method)
  validate_DAList(DAList)
  
  if (!is.null(DAList$tags$normalized) && DAList$tags$normalized) {
    cli::cli_abort("Data in DAList are already normalized.")
  }
  
  if (norm_method %in% c("log2","vsn","gi")) {
    normalized_data <- do.call(
      what = paste0(norm_method, "Norm"),
      args = c(list(dat = DAList$data), list(...))
    )
  } else if (norm_method %in% c("median","mean","quantile","cycloess","rlr")) {
    log2_dat <- if (input_is_log2) as.matrix(DAList$data) else log2Norm(dat = DAList$data)
    normalized_data <- do.call(
      what = paste0(norm_method, "Norm"),
      args = c(list(logDat = log2_dat), list(...))
    )
  } else {
    cli::cli_abort("{.arg {norm_method}} is not a valid normalization method")
  }
  
  DAList$data <- normalized_data
  DAList$tags$normalized <- TRUE
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

##' @rdname norm_functions
##' @keywords internal
# cycloessNorm <- function(logDat) {
#   cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method = "fast")
#   colnames(cycLoessNormed) <- colnames(logDat)
#   row.names(cycLoessNormed) <- rownames(logDat)
#   return(as.matrix(cycLoessNormed))
# }

#' Cyclic loess normalization (optionally within groups)
#' Expects log2 data. If `groups` is given, normalizes each group independently.
#' @rdname norm_functions
#' @keywords internal
#'
cycloessNorm <- function(logDat, groups = NULL, method = "fast", ...) {
  logDat <- as.matrix(logDat)
  
  # Global (original behavior)
  if (is.null(groups)) {
    out <- limma::normalizeCyclicLoess(logDat, method = method, ...)
    colnames(out) <- colnames(logDat)
    rownames(out) <- rownames(logDat)
    return(out)
  }
  
  # Within-group normalization (independent)
  stopifnot(length(groups) == ncol(logDat))
  if (!is.factor(groups)) groups <- factor(groups, levels = unique(groups))
  
  out <- matrix(NA_real_, nrow = nrow(logDat), ncol = ncol(logDat),
                dimnames = dimnames(logDat))
  for (g in levels(groups)) {
    idx <- which(groups == g)
    sub <- logDat[, idx, drop = FALSE]
    normed <- limma::normalizeCyclicLoess(sub, method = method, ...)
    out[, idx] <- normed
  }
  out
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

#' Plot histograms highlighting imputed values (Perseus-style imputation)
#'
#' @param logDat_before  matrix/data.frame of log2 intensities with NAs (pre-imputation)
#' @param logDat_after   matrix/data.frame of log2 intensities after perseus_impute()
#' @param samples        character/integer vector of columns to plot (default: all)
#' @param bins           number of histogram bins
#' @param facet_ncol     number of columns in facet wrap
#' @param overlay        if TRUE, overlays observed + imputed in the same panel per sample;
#'                      if FALSE, uses position = "identity" stacked fills
#' @return ggplot object (faceted histograms across selected samples)
#' @details
#'   Values that were NA in `logDat_before` and non-NA in `logDat_after`
#'   are labeled as "Imputed". All others are "Observed".
plot_perseus_imputation <- function(logDat_before,
                                    logDat_after,
                                    samples = NULL,
                                    bins = 50,
                                    facet_ncol = 4,
                                    overlay = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install.packages('ggplot2') to use plot_perseus_imputation().")
  }
  stopifnot(nrow(logDat_before) == nrow(logDat_after),
            ncol(logDat_before) == ncol(logDat_after))
  
  X0 <- as.matrix(logDat_before)
  X1 <- as.matrix(logDat_after)
  
  # determine which columns to plot
  if (is.null(samples)) {
    cols <- seq_len(ncol(X0))
  } else if (is.character(samples)) {
    cols <- match(samples, colnames(X0))
    if (anyNA(cols)) stop("Some 'samples' not found in column names.")
  } else if (is.numeric(samples)) {
    cols <- samples
  } else {
    stop("'samples' must be NULL, character, or numeric indices.")
  }
  
  # build long data frame with "Observed" vs "Imputed" labels
  make_df <- function(j) {
    # imputed where before is NA but after is not NA
    was_na <- is.na(X0[, j])
    now_ok <- !is.na(X1[, j])
    imputed_flag <- was_na & now_ok
    
    # keep non-NA after-imputation values for plotting
    vals <- X1[!is.na(X1[, j]), j]
    status <- ifelse(imputed_flag[!is.na(X1[, j])], "Imputed", "Observed")
    data.frame(
      sample = rep(colnames(X1)[j], length(vals)),
      value = vals,
      status = factor(status, levels = c("Observed", "Imputed")),
      stringsAsFactors = FALSE
    )
  }
  dflist <- lapply(cols, make_df)
  DF <- do.call(rbind, dflist)
  
  library(ggplot2)
  p <- ggplot(DF, aes(x = value, fill = status)) +
    { if (overlay)
      geom_histogram(aes(y = after_stat(count)),
                     bins = bins, position = "identity", alpha = 0.6)
      else
        geom_histogram(bins = bins, position = "stack")
    } +
    facet_wrap(~ sample, scales = "free_y", ncol = facet_ncol) +
    labs(x = "log2 intensity", y = "Count",
         title = "Perseus-style imputation: observed vs imputed") +
    scale_fill_manual(values = c("Observed" = "#4C78A8", "Imputed" = "#E45756")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  p
}

#' Perseus-style MNAR imputation (left-censoring) on log2 data
#'
#' @param logDat   numeric matrix/data.frame (rows = features, cols = samples), on log2 scale
#' @param shift    numeric, how far to shift down the mean in SD units (default 1.8)
#' @param width    numeric, fraction of SD used as imputation SD (default 0.3)
#' @param robust   logical, use median/MAD (TRUE) or mean/SD (FALSE) per sample
#' @param min_obs_per_sample integer, minimum non-missing values needed in a sample to estimate stats (default 5).
#'                           If fewer, falls back to pooled stats across all samples.
#' @param seed     integer or NULL, set for reproducible draws
#' @return numeric matrix of same dimensions with NAs imputed
#' @details
#'   This implements the common Perseus approach:
#'   For each sample j, let μ_j and σ_j be location and scale of its observed (log2) values.
#'   Missing entries in sample j are drawn from Normal( μ_j - shift*σ_j, (width*σ_j)^2 ).
#'   Proteins missing in *all* samples remain NA (uninformative).
#'
#'   Typical defaults are shift ≈ 1.8 and width ≈ 0.3 on log2 data.
#'
#' @examples
#'   # X is log2-intensity matrix with NAs
#'   X_imp <- perseus_impute(X, shift=1.8, width=0.3, robust=TRUE, seed=123)
perseus_impute <- function(logDat,
                           shift = 1.8,
                           width = 0.3,
                           robust = TRUE,
                           min_obs_per_sample = 5,
                           seed = NULL) {
  # basic checks
  if (!is.matrix(logDat) && !is.data.frame(logDat)) {
    stop("logDat must be a matrix or data.frame of log2 intensities.")
  }
  X <- as.matrix(logDat)
  storage.mode(X) <- "double"
  
  if (!is.null(seed)) set.seed(seed)
  
  # helper functions
  est_loc <- function(x, robust) {
    if (robust) stats::median(x, na.rm = TRUE) else base::mean(x, na.rm = TRUE)
  }
  est_scale <- function(x, robust) {
    if (robust) stats::mad(x, center = stats::median(x, na.rm = TRUE), na.rm = TRUE) * 1.4826
    else stats::sd(x, na.rm = TRUE)
  }
  
  # pooled stats as fallback
  x_obs_all <- as.numeric(X[!is.na(X)])
  pooled_mu  <- est_loc(x_obs_all, robust)
  pooled_sd  <- est_scale(x_obs_all, robust)
  if (is.na(pooled_sd) || pooled_sd <= 0) pooled_sd <- 1.0  # guard
  
  n_col <- ncol(X)
  for (j in seq_len(n_col)) {
    xj <- X[, j]
    miss_idx <- which(is.na(xj))
    if (length(miss_idx) == 0) next
    
    obs_j <- xj[!is.na(xj)]
    if (length(obs_j) >= min_obs_per_sample) {
      mu_j <- est_loc(obs_j, robust)
      sd_j <- est_scale(obs_j, robust)
      if (is.na(sd_j) || sd_j <= 0) {
        mu_j <- pooled_mu; sd_j <- pooled_sd
      }
    } else {
      mu_j <- pooled_mu; sd_j <- pooled_sd
    }
    
    mu_imp <- mu_j - shift * sd_j
    sd_imp <- width * sd_j
    if (!is.finite(sd_imp) || sd_imp <= 0) sd_imp <- width * pooled_sd
    
    draws <- stats::rnorm(length(miss_idx), mean = mu_imp, sd = sd_imp)
    xj[miss_idx] <- draws
    X[, j] <- xj
  }
  
  # keep rows entirely missing as NA (uninformative across all samples)
  all_missing_rows <- which(rowSums(is.na(logDat)) == ncol(logDat))
  if (length(all_missing_rows) > 0) {
    X[all_missing_rows, ] <- NA_real_
  }
  
  dimnames(X) <- dimnames(logDat)
  X
}
