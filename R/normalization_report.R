## write_norm_report uses functions in multiple other files:
## normalization_metrics.R contains functions for numerically evaluating normalization methods
## normalization_plotting.R contains functions to plot these metrics

#' Create a normalization report (contrast-aware)
#'
#' Saves a PDF report containing plots/metrics to compare normalization methods.
#' If `DAList$data_per_contrast` exists, evaluates normalization **for each contrast**
#' using the matrices in that slot (leaving `DAList$data` untouched). Otherwise,
#' evaluates on `DAList$data` (legacy behavior).
#'
#' @param DAList A DAList.
#' @param grouping_column Name of the metadata column that gives sample groups
#'   (used by within-group cyclic loess). Must contain at least two groups.
#' @param output_dir Directory to save the report (created if missing). Defaults to cwd.
#' @param filename File name of the report (must end with .pdf). Default "normalization_report.pdf".
#' @param overwrite Overwrite if the report already exists? Default FALSE.
#' @param suppress_zoom_legend Remove the legend from the zoomed log2ratio plot? Default FALSE.
#' @param use_ggrastr Use ggrastr in MD plots to reduce size (if installed)? Default FALSE.
#' @param input_is_log2 Logical. If TRUE, indicates per-contrast (or global) input matrices
#'   are already on the log2 scale for methods that expect log2 input. Default FALSE.
#' @param contrasts Optional character vector: restrict report to these contrasts when
#'   `data_per_contrast` is present. Default = all.
#'
#' @return Invisibly returns the input DAList.
#' @export
write_norm_report <- function(DAList,
                              grouping_column = NULL,
                              output_dir = NULL,
                              filename = NULL,
                              overwrite = FALSE,
                              suppress_zoom_legend = FALSE,
                              use_ggrastr = FALSE,
                              input_is_log2 = FALSE,
                              contrasts = NULL,
                              sample_id_col = NULL,          # NEW: explicitly specify the metadata column matching matrix column names
                              groups_override = NULL,       # NEW: pass a named vector of groups (names = sample IDs)
                              metrics_csv = NULL        ) {   # <- NEW   
  
  #################################
  ## Check args and set defaults ##
  #################################
  input_DAList <- validate_DAList(DAList)
  
  # Reset ggplot theme to avoid interference from global theme changes in other tests
  old_theme <- ggplot2::theme_get()
  on.exit(ggplot2::theme_set(old_theme), add = TRUE)
  ggplot2::theme_set(ggplot2::theme_gray())
  
  # If data are already normalized, stop early — this report is meant for raw data.
  #If DAList$tags$normalized == TRUE, write_norm_report() now throws a dedicated error before ever looking at grouping_column.
  #That matches the test expecting an “already normalized” error rather than grouping_column cannot be empty.
  if (!is.null(DAList$tags$normalized) && isTRUE(DAList$tags$normalized)) {
    cli::cli_abort(
      "Input data already normalized; {.fn write_norm_report} expects raw (unnormalized) data."
    )
  }
  
  if (!is.null(grouping_column)) {
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))
    }
    if (grouping_column %notin% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata of {.arg DAList}",
                       "i" = "Check with {.code colnames(DAList$metadata)}."))
    }
    groups_vec <- as.character(DAList$metadata[[grouping_column]])
    if (length(unique(groups_vec)) < 2) {
      cli::cli_abort(c("Column {.arg {grouping_column}} does not contain at least two different groups",
                       "!" = "Cannot calculate all normalization metrics without at least two groups"))
    }
  } else {
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  }
  
  # Output path defaults
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform("{.arg output_dir} argument is empty. Using current working directory:\n{.path {output_dir}}")
  }
  if (is.null(filename)) {
    filename <- "normalization_report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg filename} argument is empty. Saving report to: {.path {file.path(output_dir, filename)}}"))
  }
  
  # Ensure dir / filename
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  validate_filename(filename, allowed_exts = c("pdf"))
  report_path <- file.path(output_dir, filename)
  if (file.exists(report_path)) {
    if (overwrite) {
      cli::cli_inform("{.path {filename}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
    } else {
      cli::cli_abort(c("{.path {filename}} already exists in {.path {output_dir}}",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Give {.arg filename} a unique name or set {.arg overwrite} to {.val TRUE}"))
    }
  }
  
  ########################################
  ## Resolve matrices (per-contrast/global)
  ########################################
  has_dpc <- !is.null(DAList$data_per_contrast) && length(DAList$data_per_contrast) > 0
  
  if (has_dpc) {
    all_ct <- names(DAList$data_per_contrast)
    if (is.null(all_ct) || !length(all_ct)) {
      cli::cli_abort("DAList$data_per_contrast must be a *named* list of contrast matrices.")
    }
    target_ct <- if (is.null(contrasts)) all_ct else intersect(contrasts, all_ct)
    if (!length(target_ct)) cli::cli_abort("No matching contrasts in DAList$data_per_contrast.")
    
    # Build list of per-contrast matrices and aligned group vectors
    mats <- list(); groups_map <- list()
    # We assume there is a sample id column matching colnames of matrices. Try common ones.
    # Resolve sample identifier column
    if (!is.null(groups_override)) {
      if (is.null(names(groups_override))) {
        cli::cli_abort("'groups_override' must be a named vector with names matching sample (column) IDs.")
      }
      id_map <- groups_override
    } else {
      id_col <- sample_id_col
      if (is.null(id_col)) {
        # Auto-detect best matching metadata column to matrix colnames
        candidates <- colnames(DAList$metadata)
        all_cols <- unique(unlist(lapply(DAList$data_per_contrast[names(DAList$data_per_contrast)], function(x){
          if (is.matrix(x) || is.data.frame(x)) return(colnames(x))
          if (is.list(x) && "log2" %in% names(x)) return(colnames(x$log2))
          character(0)
        })))
        scores <- vapply(candidates, function(cc){
          vals <- as.character(DAList$metadata[[cc]])
          mean(vals %in% all_cols)
        }, numeric(1))
        if (length(scores) && max(scores, na.rm = TRUE) > 0) {
          id_col <- names(which.max(scores))
        } else {
          id_col <- intersect(c("sample_id","Raw.file","file","sample"), colnames(DAList$metadata))[1]
          if (is.na(id_col)) id_col <- colnames(DAList$metadata)[1]
        }
      }
      id_map <- setNames(DAList$metadata[[grouping_column]], DAList$metadata[[id_col]])
    }
    
    extract_mat <- function(x) {
      if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
      if (is.list(x) && "log2" %in% names(x) && (is.matrix(x$log2) || is.data.frame(x$log2))) return(as.matrix(x$log2))
      stop("Unsupported per-contrast structure. Expected matrix/data.frame or list with $log2.")
    }
    
    for (ct in target_ct) {
      X <- extract_mat(DAList$data_per_contrast[[ct]])
      # Align groups by matrix columns (names in id_map)
      if (!is.null(colnames(X))) {
        g <- id_map[colnames(X)]
      } else {
        cli::cli_abort(sprintf("Per-contrast matrix for '%s' lacks column names; cannot align groups.", ct))
      }
      if (anyNA(g)) {
        missing <- colnames(X)[is.na(g)]
        cli::cli_abort(sprintf("Group labels missing for samples in contrast '%s': %s", ct, paste(missing, collapse=", ")))
      }
      mats[[ct]] <- X
      groups_map[[ct]] <- as.character(g)
    }
    
    ###########################
    ## Compute metrics/plots  ##
    ###########################
    cli::cli_inform("Starting per-contrast normalizations for {length(target_ct)} contrasts…")
    
    plot_pages <- list()
    for (ct in target_ct) {
      normList <- apply_all_normalizations_contrast(mats[[ct]], 
                                                    input_is_log2 = input_is_log2, 
                                                    groups = groups_map[[ct]])
      
      a <- pn_plot_PCV(normList, groups_map[[ct]])
      b <- pn_plot_PMAD(normList, groups_map[[ct]])
      c <- pn_plot_PEV(normList, groups_map[[ct]])
      d <- pn_plot_COR(normList, groups_map[[ct]])
      e <- pn_plot_log2ratio(normList, groups_map[[ct]])
      f <- pn_plot_log2ratio(normList, groups_map[[ct]], zoom = TRUE, legend = !suppress_zoom_legend)
      page_1 <- a + b + c + d + e + f + patchwork::plot_layout(ncol = 3) + patchwork::plot_annotation(title = paste0("Normalization metrics — ", ct))
      
      page_2 <- pn_plot_MD(normList, groups_map[[ct]], use_ggrastr)
      plot_pages[[length(plot_pages) + 1]] <- list(page_1 = page_1, page_2 = page_2)
      
      if (!is.null(metrics_csv)) {
        metrics_tbl <- collect_norm_metrics(normList, groups_map[[ct]], contrast = ct)
        metrics_accum <- if (exists(".norm_metrics_accum", inherits = FALSE)) .norm_metrics_accum else NULL
        .norm_metrics_accum <- rbind(metrics_accum, metrics_tbl)
      }
      
    }
    
    ###############################
    ## Save plots, check, return ##
    ###############################
    grobs <- unlist(lapply(seq_along(plot_pages), function(i) {
      pg <- plot_pages[[i]]
      list(patchwork::patchworkGrob(pg$page_1), pg$page_2)
    }), recursive = FALSE)
    
    if (!is.null(metrics_csv)) {
      if (exists(".norm_metrics_accum", inherits = FALSE) && !is.null(.norm_metrics_accum)) {
        utils::write.csv(.norm_metrics_accum, metrics_csv, row.names = FALSE)
      } else {
        cli::cli_warn("No metrics accumulated; CSV not written.")
      }
    }
    
    cli::cli_inform("Saving report to: {.path {report_path}}")
    ggplot2::ggsave(report_path,
                    plot = gridExtra::marrangeGrob(grobs = grobs, nrow = 1, ncol = 1, top = NA),
                    height = 8.5, width = 11, units = "in", useDingbats = TRUE)
    
    if (!file.exists(report_path)) cli::cli_abort(c("Failed to create {.path {report_path}}"))
    return(invisible(input_DAList))
  }
  
  ############################
  ## Legacy global behavior ##
  ############################
  cli::cli_inform("Starting normalizations (global)")
  normList <- apply_all_normalizations_contrast(DAList$data, input_is_log2 = input_is_log2, groups = groups_vec)
  cli::cli_inform("Normalizations finished")
  
  if (!is.null(metrics_csv)) {
    metrics_tbl <- collect_norm_metrics(normList, groups_vec, contrast = "GLOBAL")
    utils::write.csv(metrics_tbl, metrics_csv, row.names = FALSE)
  }
  
  a <- pn_plot_PCV(normList, groups_vec)
  b <- pn_plot_PMAD(normList, groups_vec)
  c <- pn_plot_PEV(normList, groups_vec)
  d <- pn_plot_COR(normList, groups_vec)
  e <- pn_plot_log2ratio(normList, groups_vec)
  f <- pn_plot_log2ratio(normList, groups_vec, zoom = TRUE, legend = !suppress_zoom_legend)
  page_1 <- a + b + c + d + e + f + patchwork::plot_layout(ncol = 3)
  
  page_2 <- pn_plot_MD(normList, groups_vec, use_ggrastr)
  
  cli::cli_inform("Saving report to: {.path {report_path}}")
  ggplot2::ggsave(report_path,
                  plot = gridExtra::marrangeGrob(grobs = list(patchwork::patchworkGrob(page_1), page_2), nrow = 1, ncol = 1, top = NA),
                  height = 8.5, width = 11, units = "in", useDingbats = TRUE)
  
  if (!file.exists(report_path)) cli::cli_abort(c("Failed to create {.path {report_path}}"))
  invisible(input_DAList)
}


#' Apply all normalization methods to a matrix (contrast-aware helper)
#'
#' Takes in a matrix of intensities and applies 8 normalization methods, returning
#' a named list of matrices. `input_is_log2` controls whether to log2 internally
#' for methods that expect log2 input. `groups` is only used for `cycloess`.
#'
#' @param data Matrix/data.frame (rows = features, cols = samples).
#' @param input_is_log2 Logical, data are already on log2 scale for log2 methods?
#' @param groups Character vector of group labels (length = ncol(data)); used for cycloess.
#'
#' @return A named list of 8 matrices.
#' @keywords internal
apply_all_normalizations_contrast <- function(data, input_is_log2 = FALSE, groups = NULL) {
  X <- if (is.matrix(data)) data else as.matrix(data)
  
  normList <- list()
  
  # Prepare log2 for the methods that need it
  log2_X <- if (input_is_log2) X else log2Norm(dat = X)
  
  ## Apply methods
  normList[["log2"]]     <- log2_X
  normList[["median"]]   <- medianNorm(logDat = log2_X)
  normList[["mean"]]     <- meanNorm(logDat = log2_X)
  normList[["vsn"]]      <- vsnNorm(dat = X)
  normList[["quantile"]] <- quantileNorm(logDat = log2_X)
  normList[["cycloess"]] <- cycloessNorm(logDat = log2_X, groups = groups)
  normList[["rlr"]]      <- rlrNorm(logDat = log2_X)
  normList[["gi"]]       <- giNorm(dat = X)
  
  normList
}

#' Apply all normalization methods (legacy wrapper)
#'
#' Thin wrapper around `apply_all_normalizations_contrast()` to preserve
#' the original function name used in tests and older code.
#'
#' @param data Matrix/data.frame of intensities (features x samples).
#' @param input_is_log2 Logical, whether `data` are already on log2 scale.
#' @param groups Optional character vector of group labels, used for cycloess.
#'
#' @return A named list of matrices as returned by `apply_all_normalizations_contrast()`.
#' @keywords internal
apply_all_normalizations <- function(data, input_is_log2 = FALSE, groups = NULL) {
  apply_all_normalizations_contrast(
    data,
    input_is_log2 = input_is_log2,
    groups        = groups
  )
}

#' Write a post-normalization evaluation report for a single method
#'
#' Creates a PDF report with the same panels used in `write_norm_report()`,
#' but **using only the already-normalized values** from your `DAList`.
#' This is ideal for evaluating what a chosen normalization (e.g., cyclic loess)
#' actually did to your data.
#'
#' If `DAList$data_per_contrast` exists, a two-page block is produced **for each contrast**:
#' metrics grid (PCV, PMAD, PEV, COR, Log2 ratio and zoom) and faceted MD plots.
#' If no per-contrast data exist, a single two-page report is produced from `DAList$data`.
#'
#' @param DAList A DAList containing already-normalized matrices (per-contrast or global).
#' @param norm_label Character scalar used as the list name in plots (e.g., "cycloess").
#' @param grouping_column Name of the metadata column giving sample groups. Must contain ≥2 groups.
#' @param output_dir Directory to save the PDF (created if missing). Default: working directory.
#' @param filename File name for the PDF (must end with `.pdf`). Default: "norm_eval_report.pdf".
#' @param overwrite Overwrite existing file? Default FALSE.
#' @param suppress_zoom_legend Remove legend from the zoomed log2 ratio plot? Default FALSE.
#' @param use_ggrastr Use ggrastr in MD plots (if installed) to reduce file size? Default FALSE.
#' @param input_is_log2 Logical, whether the matrices are already on the log2 scale for metrics
#'   that expect log2 (most are visualization/relative; this flag is unused here but kept for symmetry).
#' @param contrasts Optional character vector: restrict to these contrasts when `data_per_contrast` is present.
#' @param sample_id_col Optional name of the metadata column whose values match matrix colnames.
#'   If `NULL`, the function will try to auto-detect a suitable column.
#' @param groups_override Optional named vector of group labels with names equal to **matrix column names**.
#'   If supplied, overrides `grouping_column`/`sample_id_col` matching.
#'
#' @return Invisibly returns `DAList`.
#' @export
write_norm_eval_report <- function(DAList,
                                   norm_label = "cycloess",
                                   grouping_column,
                                   output_dir = NULL,
                                   filename = NULL,
                                   overwrite = FALSE,
                                   suppress_zoom_legend = FALSE,
                                   use_ggrastr = FALSE,
                                   input_is_log2 = TRUE,
                                   contrasts = NULL,
                                   sample_id_col = NULL,
                                   groups_override = NULL,
                                   metrics_csv = NULL) {
  # Validate
  validate_DAList(DAList)
  if (missing(grouping_column) || length(grouping_column) != 1) {
    cli::cli_abort("Provide a single {.arg grouping_column} (metadata column name).")
  }
  if (grouping_column %notin% colnames(DAList$metadata)) {
    cli::cli_abort("Column {.arg {grouping_column}} not found in DAList$metadata.")
  }
  groups_vec_all <- as.character(DAList$metadata[[grouping_column]])
  if (length(unique(groups_vec_all)) < 2) {
    cli::cli_abort("{.arg grouping_column} must contain at least two groups.")
  }
  
  # Output setup
  if (is.null(output_dir)) output_dir <- getwd()
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (is.null(filename)) filename <- "norm_eval_report.pdf"
  validate_filename(filename, allowed_exts = c("pdf"))
  report_path <- file.path(output_dir, filename)
  if (file.exists(report_path) && !overwrite) {
    cli::cli_abort(c("{.path {filename}} already exists in {.path {output_dir}}",
                     "!" = "and {.arg overwrite} == {.val FALSE}",
                     "i" = "Use a new {.arg filename} or set {.arg overwrite} = TRUE"))
  }
  
  has_dpc <- !is.null(DAList$data_per_contrast) && length(DAList$data_per_contrast) > 0
  
  extract_mat <- function(x) {
    if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
    if (is.list(x) && "log2" %in% names(x) && (is.matrix(x$log2) || is.data.frame(x$log2))) return(as.matrix(x$log2))
    stop("Unsupported structure: expected matrix/data.frame or list with $log2.")
  }
  
  # Determine ID mapping for groups
  if (!is.null(groups_override)) {
    if (is.null(names(groups_override))) {
      cli::cli_abort("'groups_override' must be a named vector with names that match matrix colnames.")
    }
    id_map <- groups_override
  } else if (!is.null(sample_id_col)) {
    if (sample_id_col %notin% colnames(DAList$metadata)) {
      cli::cli_abort("sample_id_col '{sample_id_col}' not found in metadata.")
    }
    id_map <- setNames(DAList$metadata[[grouping_column]], DAList$metadata[[sample_id_col]])
  } else {
    # Try to auto-detect best metadata column that matches colnames across contrasts
    candidates <- colnames(DAList$metadata)
    all_cols <- if (has_dpc) unique(unlist(lapply(DAList$data_per_contrast, function(x){
      if (is.matrix(x) || is.data.frame(x)) return(colnames(x))
      if (is.list(x) && "log2" %in% names(x)) return(colnames(x$log2))
      character(0)
    }))) else colnames(DAList$data)
    scores <- vapply(candidates, function(cc){
      vals <- as.character(DAList$metadata[[cc]])
      mean(vals %in% all_cols)
    }, numeric(1))
    if (length(scores) && max(scores, na.rm = TRUE) > 0) {
      sample_id_col <- names(which.max(scores))
      id_map <- setNames(DAList$metadata[[grouping_column]], DAList$metadata[[sample_id_col]])
    } else {
      cli::cli_abort("Could not auto-detect a metadata ID column that matches matrix colnames. Provide {.arg sample_id_col} or {.arg groups_override}.")
    }
  }
  
  grobs <- list()
  
  if (has_dpc) {
    all_ct <- names(DAList$data_per_contrast)
    target_ct <- if (is.null(contrasts)) all_ct else intersect(contrasts, all_ct)
    if (!length(target_ct)) cli::cli_abort("No matching contrasts found in data_per_contrast.")
    
    for (ct in target_ct) {
      X <- extract_mat(DAList$data_per_contrast[[ct]])
      if (is.null(colnames(X))) cli::cli_abort("Per-contrast matrix for '{ct}' has no column names.")
      g <- id_map[colnames(X)]
      if (anyNA(g)) {
        missing <- colnames(X)[is.na(g)]
        cli::cli_abort("Group labels missing for samples in contrast '{ct}': {paste(missing, collapse=", ")}")
      }
      
      # Build a single-method normList for plotting helpers
     # normList <- list(); normList[[norm_label]] <- X
      
      # NEW — always provide a 'log2' entry for plotting helpers
      X <- as.matrix(X)
      if (is.null(dim(X))) X <- matrix(X, nrow = length(X), ncol = 1)  # extra safety if ever 1-D
      normList <- list(log2 = X)
      
      if (!is.null(metrics_csv)) {
        mt <- collect_norm_metrics(normList, g, contrast = ct)
        metrics_accum <- if (exists(".norm_eval_metrics_accum", inherits = FALSE)) .norm_eval_metrics_accum else NULL
        .norm_eval_metrics_accum <- rbind(metrics_accum, mt)
      }
      
      a <- pn_plot_PCV(normList, g)
      b <- pn_plot_PMAD(normList, g)
      c <- pn_plot_PEV(normList, g)
      d <- pn_plot_COR(normList, g)
      e <- pn_plot_log2ratio(normList, g)
      f <- pn_plot_log2ratio(normList, g, zoom = TRUE, legend = !suppress_zoom_legend)
      page_1 <- a + b + c + d + e + f + patchwork::plot_layout(ncol = 3) +
        patchwork::plot_annotation(title = paste0("Post-normalization metrics — ", ct, " (", norm_label, ")"))
      
      page_2 <- pn_plot_MD(normList, g, use_ggrastr)
      
      grobs <- c(grobs, list(patchwork::patchworkGrob(page_1), page_2))
    }
  } else {
    X <- as.matrix(DAList$data)
    if (is.null(colnames(X))) cli::cli_abort("Global matrix has no column names.")
    # Attempt to align groups
    if (!is.null(names(id_map))) {
      g <- id_map[colnames(X)]
      if (anyNA(g)) {
        cli::cli_abort("Group labels missing for some samples in global data.")
      }
    } else {
      cli::cli_abort("Need named group mapping for global data alignment.")
    }
    
   # normList <- list(); normList[[norm_label]] <- X
    # some plotting functions expect normList$log2 to exist
    # NEW — always provide a 'log2' entry for plotting helpers
    X <- as.matrix(X)
    if (is.null(dim(X))) X <- matrix(X, nrow = length(X), ncol = 1)  # extra safety if ever 1-D
    normList <- list(log2 = X)
    
    
    
    a <- pn_plot_PCV(normList, g)
    b <- pn_plot_PMAD(normList, g)
    c <- pn_plot_PEV(normList, g)
    d <- pn_plot_COR(normList, g)
    e <- pn_plot_log2ratio(normList, g)
    f <- pn_plot_log2ratio(normList, g, zoom = TRUE, legend = !suppress_zoom_legend)
    page_1 <- a + b + c + d + e + f + patchwork::plot_layout(ncol = 3)
    
    page_2 <- pn_plot_MD(normList, g, use_ggrastr)
    
    grobs <- list(patchwork::patchworkGrob(page_1), page_2)
  }
  
  if (!is.null(metrics_csv)) {
    if (exists(".norm_eval_metrics_accum", inherits = FALSE) && !is.null(.norm_eval_metrics_accum)) {
      utils::write.csv(.norm_eval_metrics_accum, metrics_csv, row.names = FALSE)
    } else {
      cli::cli_warn("No metrics accumulated; CSV not written for post-normalization report.")
    }
  }
  
  cli::cli_inform("Saving post-normalization report to: {.path {report_path}}")
  ggplot2::ggsave(report_path,
                  plot = gridExtra::marrangeGrob(grobs = grobs, nrow = 1, ncol = 1, top = NA),
                  height = 8.5, width = 11, units = "in", useDingbats = TRUE)
  
  if (!file.exists(report_path)) cli::cli_abort("Failed to create {.path {report_path}}")
  invisible(DAList)
}

# -------------------------------
# Metrics collection helper
# -------------------------------
#' Collect basic normalization metrics from a normList and groups
#'
#' Computes summary metrics for each method in a `normList`:
#' - PCV: pooled coefficient of variation (computed on linear scale)
#' - PMAD: pooled median absolute deviation (log2 scale)
#' - PEV: pooled estimate of variance (log2 scale)
#' - COR: mean pairwise Pearson correlation within groups (log2 scale)
#'
#' @param normList Named list of matrices (methods × samples). Names are method labels.
#' @param groups Character vector of group labels aligned to columns of matrices.
#' @param contrast Character scalar to label the contrast in the output.
#'
#' @return data.frame with columns: contrast, method, metric, value
#' @keywords internal
collect_norm_metrics <- function(normList, groups, contrast = NA_character_) {
  stopifnot(is.list(normList), length(groups) > 1)
  out <- list()
  for (m in names(normList)) {
    Xlog <- as.matrix(normList[[m]])
    storage.mode(Xlog) <- "double"
    
    # PCV on linear scale
    Xlin  <- 2^Xlog
    cv_vec <- apply(Xlin, 1, function(v) {
      mu <- mean(v, na.rm = TRUE); sd <- stats::sd(v, na.rm = TRUE)
      if (!is.finite(mu) || mu == 0) return(NA_real_)
      sd / mu
    })
    pcv  <- stats::median(cv_vec, na.rm = TRUE)
    
    # PMAD and PEV on log2 scale
    mad_vec <- apply(Xlog, 1, stats::mad, na.rm = TRUE)
    pmad    <- stats::median(mad_vec, na.rm = TRUE)
    
    var_vec <- apply(Xlog, 1, stats::var, na.rm = TRUE)
    pev     <- mean(var_vec, na.rm = TRUE)
    
    # COR: mean pairwise correlation within each group, then average over groups
    cor_by_group <- tapply(seq_along(groups), groups, function(idx) {
      if (length(idx) < 2) return(NA_real_)
      C <- suppressWarnings(stats::cor(Xlog[, idx, drop = FALSE],
                                       use = "pairwise.complete.obs"))
      mean(C[upper.tri(C)], na.rm = TRUE)
    })
    cor_mean <- mean(unlist(cor_by_group), na.rm = TRUE)
    
    out[[length(out) + 1]] <- data.frame(
      contrast = contrast,
      method   = m,
      metric   = c("PCV","PMAD","PEV","COR"),
      value    = c(pcv, pmad, pev, cor_mean),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
