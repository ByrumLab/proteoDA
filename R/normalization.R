#' Normalize data in a DAList (contrast-aware)
#'
#' Normalizes data using one of several methods. If `DAList$data_per_contrast`
#' exists and is non-empty, normalization is applied **per-contrast** and the
#' results are written back to `DAList$data_per_contrast[[contrast]]`, leaving
#' `DAList$data` unchanged. Otherwise, it falls back to normalizing
#' `DAList$data` (legacy global behavior).
#'
#' **Available methods**
#' - "log2"       — Binary (base-2) log transformation of raw intensities.
#' - "median"     — Divide each sample by its median; rescale by the mean of sample medians.
#' - "mean"       — Divide each sample by its mean; rescale by the mean of sample means.
#' - "vsn"        — Variance-stabilizing normalization via `vsn::justvsn` (expects raw scale).
#' - "quantile"   — Quantile normalization via `preprocessCore::normalize.quantiles`.
#' - "cycloess"   — Cyclic loess via `limma::normalizeCyclicLoess` (can be applied within groups).
#' - "rlr"        — Robust linear regression (global) normalization.
#' - "gi"         — Global intensity normalization (expects raw scale; returns log2).
#'
#' **Log scale behavior**
#' Methods "median", "mean", "quantile", "cycloess", and "rlr" operate on log2 data.
#' By default, this function will apply `log2` internally before those methods.
#' If your inputs are already log2-transformed,
#' set `input_is_log2 = TRUE` to avoid double logging.
#'
#' @param DAList A `DAList` object.
#' @param norm_method One of: "log2","median","mean","vsn","quantile","cycloess","rlr","gi".
#' @param input_is_log2 Logical. If TRUE, indicates inputs are already on the log2 scale
#'   for methods that expect log2 input. Default FALSE.
#' @param contrasts Optional character vector to restrict which contrasts in
#'   `DAList$data_per_contrast` are normalized. Default = all available.
#' @param use_per_contrast_if_available Logical. Default TRUE. When TRUE and
#'   `DAList$data_per_contrast` is present/non-empty, normalize per-contrast and
#'   leave `DAList$data` unchanged.
#' @param ... Additional arguments forwarded to the underlying normalization function.
#'   For example, pass `groups=` to `cycloessNorm()` to normalize within groups;
#'   if supplied, `groups` will be subset/reordered to the columns present in each
#'   per-contrast matrix.
#'
#' @return The updated `DAList`. Tags are updated as follows:
#'   - Per-contrast path: `DAList$tags$per_contrast_normalized = TRUE`,
#'     `DAList$tags$norm_method_per_contrast` (named list), and
#'     `DAList$normalization_per_contrast[[contrast]]$diagnostics` (before/after
#'     per-sample summaries) are recorded. `DAList$data` is unchanged.
#'   - Global path: `DAList$tags$normalized = TRUE`, `DAList$tags$norm_method = norm_method`.
#'
#' @export
normalize_data <- function(DAList,
                           norm_method = c("log2","median","mean","vsn","quantile","cycloess","rlr","gi"),
                           input_is_log2 = FALSE,
                           contrasts = NULL,
                           use_per_contrast_if_available = TRUE,
                           ...) {
  norm_method <- rlang::arg_match(norm_method)
  validate_DAList(DAList)
  
  has_dpc <- use_per_contrast_if_available && !is.null(DAList$data_per_contrast) && length(DAList$data_per_contrast) > 0
  
  if (has_dpc) {
    # ---- Per-contrast normalization ----
    if (is.null(names(DAList$data_per_contrast))) {
      stop("DAList$data_per_contrast must be a *named* list of contrast matrices.")
    }
    target_contrasts <- if (is.null(contrasts)) names(DAList$data_per_contrast) else intersect(contrasts, names(DAList$data_per_contrast))
    if (!length(target_contrasts)) stop("No matching contrasts in DAList$data_per_contrast.")
    
    # Prepare containers for tags/diagnostics
    if (is.null(DAList$tags)) DAList$tags <- list()
    DAList$tags$per_contrast_normalized <- TRUE
    if (is.null(DAList$tags$norm_method_per_contrast)) DAList$tags$norm_method_per_contrast <- list()
    if (is.null(DAList$normalization_per_contrast)) DAList$normalization_per_contrast <- list()
    
    # Allow groups= to be supplied and subset per contrast
    dots <- list(...)
    has_groups_arg <- "groups" %in% names(dots)
    
    for (ct in target_contrasts) {
      dpc <- DAList$data_per_contrast[[ct]]
      
      # Extract a matrix to normalize and remember how to write it back
      extract_mat <- function(x) {
        if (is.matrix(x) || is.data.frame(x)) return(list(mat = as.matrix(x), write_back = function(m) { if (is.data.frame(x)) return(as.data.frame(m, check.names = FALSE)); m }))
        if (is.list(x) && "log2" %in% names(x) && (is.matrix(x$log2) || is.data.frame(x$log2))) {
          return(list(mat = as.matrix(x$log2), write_back = function(m) { x$log2 <- m; x }))
        }
        stop(sprintf("Unsupported structure for DAList$data_per_contrast[['%s']]. Expected matrix/data.frame or list with $log2.", ct))
      }
      wb <- extract_mat(dpc)
      X  <- wb$mat
      
      # Build method call arguments
      do_log_methods <- c("median","mean","quantile","cycloess","rlr")
      if (norm_method %in% c("log2","vsn","gi")) {
        call_args <- c(list(dat = X), dots)
        normalized <- do.call(paste0(norm_method, "Norm"), call_args)
      } else if (norm_method %in% do_log_methods) {
        # Optional groups subsetting
        local_dots <- dots
        if (has_groups_arg) {
          g <- dots$groups
          if (!is.null(colnames(X)) && !is.null(names(g))) {
            # Reorder by column names if named
            g <- g[colnames(X)]
          } else if (length(g) != ncol(X)) {
            stop(sprintf("Length of 'groups' (%d) must match ncol(X) (%d) for contrast '%s'.", length(g), ncol(X), ct))
          }
          local_dots$groups <- g
        }
        log2_X <- if (input_is_log2) X else log2Norm(dat = X)
        call_args <- c(list(logDat = log2_X), local_dots)
        normalized <- do.call(paste0(norm_method, "Norm"), call_args)
      } else {
        stop(sprintf("%s is not a valid normalization method", norm_method))
      }
      
      # Diagnostics (before/after per-sample summaries)
      diag <- list(
        before = list(
          sample_median = apply(if (exists("log2_X")) log2_X else (if (norm_method == "log2") log2Norm(dat = X) else X), 2, stats::median, na.rm = TRUE),
          sample_mean   = apply(if (exists("log2_X")) log2_X else (if (norm_method == "log2") log2Norm(dat = X) else X), 2, mean, na.rm = TRUE)
        ),
        after = list(
          sample_median = apply(as.matrix(normalized), 2, stats::median, na.rm = TRUE),
          sample_mean   = apply(as.matrix(normalized), 2, mean, na.rm = TRUE)
        ),
        method = norm_method,
        input_is_log2 = input_is_log2
      )
      
      # Write back preserving structure
      out_obj <- wb$write_back(normalized)
      DAList$data_per_contrast[[ct]] <- out_obj
      
      # Record tags/diagnostics
      DAList$tags$norm_method_per_contrast[[ct]] <- norm_method
      DAList$normalization_per_contrast[[ct]] <- list(diagnostics = diag)
    }
    
    validate_DAList(DAList)
    return(DAList)
  }
  
  # ---- Global (legacy) normalization ----
  if (!is.null(DAList$tags$normalized) && isTRUE(DAList$tags$normalized)) {
    cli::cli_abort("Global data in DAList$data are already normalized. To normalize per-contrast, set use_per_contrast_if_available = TRUE and populate $data_per_contrast.")
  }
  
  if (is.null(DAList$data)) cli::cli_abort("DAList$data is NULL; nothing to normalize.")
  
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
    cli::cli_abort(sprintf("%s is not a valid normalization method", norm_method))
  }
  
  DAList$data <- normalized_data
  DAList$tags$normalized <- TRUE
  DAList$tags$norm_method <- norm_method
  validate_DAList(DAList)
  DAList
}


# ---------------------------
# Normalization *kernels*
# ---------------------------

#' @rdname norm_functions
#' @keywords internal
log2Norm <- function(dat) {
  logInt <- log2(as.matrix(dat))
  logInt[is.infinite(as.matrix(logInt))] <- NA
  return(as.matrix(logInt))
}

#' @rdname norm_functions
#' @keywords internal
medianNorm <- function(logDat) {
  sampleMed <- apply(logDat, 2, stats::median, na.rm = TRUE)
  meanMed <- mean(sampleMed, na.rm = TRUE)
  out <- t(t(logDat) / sampleMed)
  out <- out * meanMed
  return(as.matrix(out))
}

#' @rdname norm_functions
#' @keywords internal
meanNorm <- function(logDat) {
  sampleMean <- apply(logDat, 2, mean, na.rm = TRUE)
  meanMean <- mean(sampleMean, na.rm = TRUE)
  out <- t(t(logDat) / sampleMean)
  out <- out * meanMean
  return(as.matrix(out))
}

#' @rdname norm_functions
#' @keywords internal
vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

#' @rdname norm_functions
#' @keywords internal
quantileNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy = TRUE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

#' Cyclic loess normalization (optionally within groups)
#' Expects log2 data. If `groups` is given, normalizes each group independently.
#' @rdname norm_functions
#' @keywords internal
cycloessNorm <- function(logDat, groups = NULL, method = "fast", ...) {
  logDat <- as.matrix(logDat)
  
  if (is.null(groups)) {
    out <- limma::normalizeCyclicLoess(logDat, method = method, ...)
    colnames(out) <- colnames(logDat)
    rownames(out) <- rownames(logDat)
    return(out)
  }
  
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
rlrNorm <- function(logDat) {
  get_rlm_coeffs_and_slopes <- function(unnormalized_col, predictors) {
    MASS::rlm(unnormalized_col ~ predictors,
              na.action = stats::na.exclude)$coefficients
  }
  row_medians <- rowMedians(logDat, na.rm = TRUE)
  coefficients <- t(apply(X = logDat,
                          MARGIN = 2,
                          FUN = get_rlm_coeffs_and_slopes,
                          predictors = row_medians))
  stopifnot(nrow(t(logDat)) == nrow(coefficients))
  normalized_mat <- t((t(logDat) - coefficients[,1]) / coefficients[,2])
  normalized_mat
}

#' @rdname norm_functions
#' @keywords internal
giNorm <- function(dat) {
  dat_mat <- as.matrix(dat)
  col_sums <- colSums(dat_mat, na.rm = TRUE)
  median_col_sum <- stats::median(col_sums)
  normMatrix <- t(t(dat_mat) / col_sums) * median_col_sum
  normLog2Matrix <- log2(normMatrix)
  normLog2Matrix[is.infinite(normLog2Matrix)] <- NA
  normLog2Matrix
}



#Add @export
#Document that if you passed a DAList to perseus_impute(..., save_before_after=TRUE), you can supply logDat_before / logDat_after from the stored slots.

#' Plot histograms highlighting imputed values (Perseus-style imputation)
#'
#' If you ran \code{perseus_impute(DAList, save_before_after=TRUE)}, you can pass:
#'   \code{logDat_before = DAList$data_per_contrast[[contrast]]$imputation$before_log2}
#'   \code{logDat_after  = DAList$data_per_contrast[[contrast]]$imputation$after_log2}
#' or use the global \code{DAList$imputation$before_log2/after_log2} when no per-contrast data exists.
#'
#' @param logDat_before  matrix/data.frame of log2 intensities with NAs (pre-imputation)
#' @param logDat_after   matrix/data.frame of log2 intensities after \code{perseus_impute()}
#' @param samples        character/integer vector of columns to plot (default: all)
#' @param bins           number of histogram bins
#' @param facet_ncol     number of columns in facet wrap
#' @param overlay        if TRUE, overlay observed+imputed per sample panel; else stacked
#' @return ggplot object (faceted histograms across selected samples)
#' @export
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
  
  make_df <- function(j) {
    was_na <- is.na(X0[, j])
    now_ok <- !is.na(X1[, j])
    imputed_flag <- was_na & now_ok
    
    vals <- X1[!is.na(X1[, j]), j]
    status <- ifelse(imputed_flag[!is.na(X1[, j])], "Imputed", "Observed")
    data.frame(
      sample = rep(colnames(X1)[j], length(vals)),
      value = vals,
      status = factor(status, levels = c("Observed", "Imputed")),
      stringsAsFactors = FALSE
    )
  }
  DF <- do.call(rbind, lapply(cols, make_df))
  
  ggplot2::ggplot(DF, ggplot2::aes(x = value, fill = status)) +
    { if (overlay)
      ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(count)),
                              bins = bins, position = "identity", alpha = 0.6)
      else
        ggplot2::geom_histogram(bins = bins, position = "stack")
    } +
    ggplot2::facet_wrap(~ sample, scales = "free_y", ncol = facet_ncol) +
    ggplot2::labs(x = "log2 intensity", y = "Count",
                  title = "Perseus-style imputation: observed vs imputed") +
    ggplot2::scale_fill_manual(values = c("Observed" = "#4C78A8", "Imputed" = "#E45756")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")
}

# Where it writes imputed data (per contrast): DAList$data_per_contrast[[contrast]]$log2
# This matches your updated plotting/limma flow that prefers per-contrast data when present.
# Fallback: If data_per_contrast is missing/empty, it imputes DAList$data.
# Matrix behavior unchanged: You can still call perseus_impute(X, ...) and get an imputed matrix back.
# Exported: The @export tag ensures the function is available after installation.

# Updated perseus_impute() (adds saving + mask)
# New args: save_before_after = FALSE and store_mask = TRUE
# Matrix input: returns the imputed matrix; if store_mask=TRUE, adds attr(result, "imputed_mask") (logical matrix of imputed positions).
# DAList input:
#   If data_per_contrast present, imputes per-contrast and writes to DAList$data_per_contrast[[contrast]]$log2.
#  If save_before_after=TRUE, also stores:
#   DAList$data_per_contrast[[contrast]]$imputation$before_log2
#   DAList$data_per_contrast[[contrast]]$imputation$after_log2
#   DAList$data_per_contrast[[contrast]]$imputation$imputed_mask
# If no data_per_contrast, imputes DAList$data and, if save_before_after=TRUE, saves:
#   DAList$imputation$before_log2
#   DAList$imputation$after_log2
#   DAList$imputation$imputed_mask
# Also adds @export for both functions.

#' Perseus-style MNAR imputation (left-censoring) on log2 data
#'
#' Applies the common Perseus imputation where, per sample, missing values
#' are drawn from a narrow normal distribution shifted down from the observed
#' distribution. Works on a plain log2 matrix/data.frame *or* on a DAList:
#' - If `DAList$data_per_contrast` exists and is non-empty, imputes each
#'   contrast's matrix and stores the result in
#'   `DAList$data_per_contrast[[contrast]]`.
#' - Otherwise, imputes the global `DAList$data` matrix in-place.
#'
#' Optionally stores before/after matrices and an imputed mask for diagnostics.
#'
#' @param x       A numeric log2 matrix/data.frame (rows=features, cols=samples),
#'                or a DAList with `$data` and optionally `$data_per_contrast`.
#' @param shift   numeric, how far to shift down the mean in SD units (default 1.8)
#' @param width   numeric, fraction of SD used as imputation SD (default 0.3)
#' @param robust  logical, use median/MAD (TRUE) or mean/SD (FALSE) per sample
#' @param min_obs_per_sample Integer. Minimum number of observed (non-missing)
#'   log2 values required in a sample before estimating its own mean (mu_j) and
#'   standard deviation (sd_j) for imputation. If a sample has fewer than this
#'   many observed values, the function falls back to using the pooled
#'   distribution across all samples. Default = 5.
#' @param seed    integer or NULL, set for reproducible draws
#' @param save_before_after logical, store full before/after matrices and mask
#'                into the DAList (or return attrs for matrix mode). Default FALSE.
#' @param store_mask logical, store/log an imputed logical mask. Default TRUE.
#'
#' @return If `x` is a matrix/data.frame, returns the imputed matrix (same dims);
#'         when `store_mask=TRUE`, attaches `attr(result, "imputed_mask")`.
#'         If `x` is a DAList, returns the updated DAList with imputed matrices
#'         written to the appropriate slots; when `save_before_after=TRUE`,
#'         persists `before_log2`, `after_log2`, and `imputed_mask`.
#'
#' @details
#' For each sample j, let μ_j and σ_j be location and scale of its observed (log2) values.
#' Missing entries in sample j are drawn from Normal( μ_j - shift*σ_j, (width*σ_j)^2 ).
#' Features missing in *all* samples remain NA (uninformative).
#'
#' This function assumes that the input data are on the **log2** scale. As a
#' safety check, if the maximum absolute value in a matrix exceeds 1e3, a
#' warning is issued suggesting that the data may not be log2-transformed.
#'
#' @examples
#' \dontrun{
#' # Matrix use
#' X_imp <- perseus_impute(X, seed=123)
#' mask  <- attr(X_imp, "imputed_mask")
#'
#' # DAList use (per-contrast if present)
#' results <- perseus_impute(results, save_before_after = TRUE, seed=1)
#' }
#'
#' @export
perseus_impute <- function(x,
                           shift = 1.8,
                           width = 0.3,
                           robust = TRUE,
                           min_obs_per_sample = 5,
                           seed = NULL,
                           save_before_after = FALSE,
                           store_mask = TRUE) {
  
  .perseus_impute_matrix <- function(logDat,
                                     shift, width, robust,
                                     min_obs_per_sample, seed,
                                     store_mask) {
    if (!is.matrix(logDat) && !is.data.frame(logDat)) {
      stop("logDat must be a matrix or data.frame of log2 intensities.")
    }
    X0 <- as.matrix(logDat)
    storage.mode(X0) <- "double"
    
    ## ---- NEW: scale sanity check -----------------------------------------
    max_abs <- suppressWarnings(max(abs(X0), na.rm = TRUE))
    if (is.finite(max_abs) && max_abs > 1e3) {
      warning("perseus_impute(): values appear very large (max abs > 1e3). ",
              "Are you sure these are log2 intensities?")
    }
    ## ----------------------------------------------------------------------
    
    if (!is.null(seed)) set.seed(seed)
    
    est_loc <- function(v, robust) {
      if (robust) stats::median(v, na.rm=TRUE) else base::mean(v, na.rm=TRUE)
    }
    est_scale <- function(v, robust) {
      if (robust) {
        stats::mad(v, center = stats::median(v, na.rm=TRUE), na.rm=TRUE) * 1.4826
      } else {
        stats::sd(v, na.rm=TRUE)
      }
    }
    
    x_obs_all <- as.numeric(X0[!is.na(X0)])
    pooled_mu <- est_loc(x_obs_all, robust)
    pooled_sd <- est_scale(x_obs_all, robust)
    if (is.na(pooled_sd) || pooled_sd <= 0) pooled_sd <- 1.0
    
    X1 <- X0
    imputed_mask <- matrix(
      FALSE,
      nrow = nrow(X0),
      ncol = ncol(X0),
      dimnames = dimnames(X0)
    )
    
    for (j in seq_len(ncol(X0))) {
      xj <- X1[, j]
      miss_idx <- which(is.na(xj))
      if (length(miss_idx) == 0) next
      
      obs_j <- X0[!is.na(X0[, j]), j]
      if (length(obs_j) >= min_obs_per_sample) {
        mu_j <- est_loc(obs_j, robust)
        sd_j <- est_scale(obs_j, robust)
        if (is.na(sd_j) || sd_j <= 0) {
          mu_j <- pooled_mu
          sd_j <- pooled_sd
        }
      } else {
        mu_j <- pooled_mu
        sd_j <- pooled_sd
      }
      
      mu_imp <- mu_j - shift * sd_j
      sd_imp <- width * sd_j
      if (!is.finite(sd_imp) || sd_imp <= 0) sd_imp <- width * pooled_sd
      
      draws <- stats::rnorm(length(miss_idx), mean = mu_imp, sd = sd_imp)
      xj[miss_idx] <- draws
      X1[, j] <- xj
      imputed_mask[miss_idx, j] <- TRUE
    }
    
    # rows entirely missing in input remain NA
    all_missing_rows <- which(rowSums(is.na(X0)) == ncol(X0))
    if (length(all_missing_rows) > 0) {
      X1[all_missing_rows, ] <- NA_real_
      imputed_mask[all_missing_rows, ] <- FALSE
    }
    
    dimnames(X1) <- dimnames(logDat)
    
    if (store_mask) attr(X1, "imputed_mask") <- imputed_mask
    X1
  }
  
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  # ---- Matrix path --------------------------------------------------------
  if (is.matrix(x) || is.data.frame(x)) {
    X1 <- .perseus_impute_matrix(
      x, shift, width, robust, min_obs_per_sample, seed, store_mask
    )
    if (save_before_after) {
      attr(X1, "before_log2") <- as.matrix(x)
    }
    return(X1)
  }
  
  # ---- DAList path --------------------------------------------------------
  if (is.list(x) && ("data" %in% names(x) || "data_per_contrast" %in% names(x))) {
    has_dpc <- !is.null(x$data_per_contrast) && length(x$data_per_contrast) > 0
    
    if (has_dpc) {
      # ensure a parallel diagnostics container that won't coerce data types
      if (is.null(x$imputation_per_contrast)) x$imputation_per_contrast <- list()
      
      for (contrast in names(x$data_per_contrast)) {
        dpc <- x$data_per_contrast[[contrast]]
        
        # Source is the matrix already on log2 scale
        if (is.matrix(dpc) || is.data.frame(dpc)) {
          source_mat  <- dpc
          orig_was_df <- is.data.frame(dpc)
        } else if (is.list(dpc) && "log2" %in% names(dpc) &&
                   (is.matrix(dpc$log2) || is.data.frame(dpc$log2))) {
          # (Support the alternative layout too, just in case)
          source_mat  <- dpc$log2
          orig_was_df <- is.data.frame(dpc$log2)
        } else if (!is.null(x$data)) {
          source_mat  <- x$data
          orig_was_df <- is.data.frame(x$data)
        } else {
          stop(sprintf("No usable matrix found for contrast '%s'.", contrast))
        }
        
        before <- as.matrix(source_mat)
        after  <- .perseus_impute_matrix(
          before, shift, width, robust, min_obs_per_sample, seed, store_mask
        )
        mask   <- attr(after, "imputed_mask"); attr(after, "imputed_mask") <- NULL

        # Write back in-place, preserving data.frame/matrix type
        x$data_per_contrast[[contrast]] <-
          if (orig_was_df) as.data.frame(after, check.names = FALSE) else after

        # Save diagnostics in a separate parallel list (does NOT change the data type)
        if (save_before_after) {
          x$imputation_per_contrast[[contrast]] <- list(
            before_log2  = before,
            after_log2   = as.matrix(x$data_per_contrast[[contrast]]),
            imputed_mask = if (store_mask) mask else NULL
          )
        }
      }
      
      x$imputation <- modifyList(
        x$imputation %||% list(),
        list(
          method             = "Perseus",
          shift              = shift,
          width              = width,
          robust             = robust,
          min_obs_per_sample = min_obs_per_sample,
          seed               = seed,
          scope              = "per-contrast"
        )
      )
      return(x)
    }
    
    # Global fallback
    if (is.null(x$data)) stop("DAList$data not found. Cannot impute.")
    before <- x$data
    after  <- .perseus_impute_matrix(
      x$data, shift, width, robust, min_obs_per_sample, seed, store_mask
    )
    mask   <- attr(after, "imputed_mask"); attr(after, "imputed_mask") <- NULL
    x$data <- after
    
    x$imputation <- modifyList(
      x$imputation %||% list(),
      list(
        method             = "Perseus",
        shift              = shift,
        width              = width,
        robust             = robust,
        min_obs_per_sample = min_obs_per_sample,
        seed               = seed,
        scope              = "global"
      )
    )
    if (save_before_after) {
      x$imputation$before_log2 <- before
      x$imputation$after_log2  <- after
      if (store_mask) x$imputation$imputed_mask <- mask
    }
    return(x)
  }
  
  stop("x must be a log2 matrix/data.frame or a DAList with $data and/or $data_per_contrast.")
}

#' Write Perseus imputation histograms for each contrast
#'
#' Uses DAList$imputation_per_contrast[[contrast]]$before_log2 / after_log2
#' (created by \code{perseus_impute(..., save_before_after = TRUE)}) to generate
#' diagnostic plots for each contrast.
#'
#' Each plot is a faceted histogram per sample, comparing observed vs imputed
#' log2 intensities. The helper \code{plot_perseus_imputation()} is used
#' internally.
#'
#' @param DAList       a DAList with \code{imputation_per_contrast} diagnostics.
#' @param out_dir      directory to save plots (created if missing). If NULL, do
#'                     not save plots to disk.
#' @param contrasts    character vector of contrast names to plot; default = all
#'                     available contrasts in \code{imputation_per_contrast}.
#' @param samples      subset of sample columns to plot (names or indices);
#'                     default = all samples.
#' @param bins         number of histogram bins.
#' @param facet_ncol   number of columns in facet layout.
#' @param overlay      logical; TRUE = overlay observed & imputed in same panel,
#'                     FALSE = stacked.
#' @param width,height plot size (inches) for saving.
#' @param dpi          dpi for saving.
#' @param device       graphics device for \code{ggplot2::ggsave}, e.g. "png".
#'
#' @return A named list of \code{ggplot} objects (invisible if you only inspect
#'         the written files).
#'
#' @export
write_perseus_imputation_plots <- function(
    DAList,
    out_dir      = "Imputation_Plots",
    contrasts    = NULL,
    samples      = NULL,
    bins         = 50,
    facet_ncol   = 4,
    overlay      = TRUE,
    width        = 7,
    height       = 5,
    dpi          = 300,
    device       = "png"
) {
  if (is.null(DAList$imputation_per_contrast) ||
      length(DAList$imputation_per_contrast) == 0) {
    stop("No per-contrast diagnostics found in DAList$imputation_per_contrast. ",
         "Run perseus_impute(..., save_before_after = TRUE) first.")
  }
  
  available <- names(DAList$imputation_per_contrast)
  if (is.null(available) || length(available) == 0) {
    stop("DAList$imputation_per_contrast has no named contrasts.")
  }
  
  if (is.null(contrasts)) {
    contrasts <- available
  } else {
    missing <- setdiff(contrasts, available)
    if (length(missing)) {
      warning("Skipping missing contrasts: ", paste(missing, collapse = ", "))
      contrasts <- intersect(contrasts, available)
    }
    if (!length(contrasts)) stop("No valid contrasts left after filtering.")
  }
  
  # create out_dir if saving
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  }
  
  plots <- list()
  
  for (ct in contrasts) {
    diag <- DAList$imputation_per_contrast[[ct]]
    if (is.null(diag$before_log2) || is.null(diag$after_log2)) {
      warning("Contrast '", ct, "' lacks before/after matrices; skipping.")
      next
    }
    
    p <- plot_perseus_imputation(
      logDat_before = diag$before_log2,
      logDat_after  = diag$after_log2,
      samples       = samples,
      bins          = bins,
      facet_ncol    = facet_ncol,
      overlay       = overlay
    ) + ggplot2::ggtitle(sprintf("Perseus imputation: %s", ct))
    
    plots[[ct]] <- p
    
    # save if requested
    if (!is.null(out_dir)) {
      safe <- gsub("[^A-Za-z0-9._-]+", "_", ct)
      fn   <- file.path(out_dir,
                        sprintf("perseus_imputation_%s.%s", safe, device))
      ggplot2::ggsave(
        filename = fn,
        plot     = p,
        width    = width,
        height   = height,
        dpi      = dpi,
        device   = device
      )
    }
  }
  
  invisible(plots)
}
