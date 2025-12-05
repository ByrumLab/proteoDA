#' Compute Rolling Standard Deviations and Z-scores for logFC Comparisons
#'
#' This function computes a rolling (moving) standard deviation of log fold change (logFC)
#' values for each comparison in the `DAList$results` list. The SD is computed over sorted data
#' using mean intensity for sorting (without modifying the data). It also computes logFC
#' Z-scores using the rolling SDs and stores all results in `DAList$tags`.
#'
#' If `binsize = "auto"`, the function evaluates a set of candidate bin sizes and chooses
#' the one that minimizes the average coefficient of variation (CV) of the moving SDs.
#'
#' If `plot = TRUE`, a plot of the final moving SDs for each comparison will be generated
#' and saved to `QC_dir` if available.
#'
#' @param DAList A DAList-like object containing `data`, `results`, and `tags` slots. Must contain:
#'   - `DAList$data`: A matrix/data.frame with numeric expression values (log2 intensities).
#'   - Optionally: `DAList$data_per_contrast`, a list of per-contrast matrices/data.frames.
#'   - `DAList$results`: A named list of comparisons each containing a `logFC` vector.
#'   - Optionally: `DAList$QC_dir` as a string path for saving QC plots.
#' @param binsize Integer or `"auto"`. The number of consecutive proteins (rows) to use in each rolling window,
#'   or `"auto"` to determine the best size automatically.
#' @param binsize_range Numeric vector. Only used if `binsize = "auto"`. A vector of candidate bin sizes to test.
#' @param plot Logical. If `TRUE`, plots the final moving SDs for each comparison and saves to QC_dir.
#' @param contrasts_file Optional CSV file with contrasts (one per row, `name = formula` syntax) if
#'   contrast names are not already available in the DAList.
#'
#' @return The input `DAList` object with updated `tags`:
#'   - `DAList$tags$movingSDs`: A list of rolling SD vectors for each comparison.
#'   - `DAList$tags$logFC_z_scores`: A list of Z-score vectors (logFC divided by rolling SD).
#'
#' @export
compute_movingSD_zscores <- function(DAList,
                                     binsize = "auto",
                                     binsize_range = c(50, 100, 200, 400),
                                     plot = FALSE,
                                     contrasts_file = NULL) {
  
  rolling_sd <- function(x, window_size) {
    n_x <- length(x)
    end_idx <- n_x - window_size + 1
    if (end_idx < 1) return(rep(NA_real_, n_x))
    sds <- vapply(
      X   = seq_len(end_idx),
      FUN = function(i) stats::sd(x[i:(i + window_size - 1)], na.rm = TRUE),
      FUN.VALUE = numeric(1)
    )
    c(sds, rep(tail(sds, 1), n_x - end_idx))
  }
  
  # Helper: choose which expression matrix to use for a given contrast
  # Prefers DAList$data_per_contrast[[contrast]] if available; otherwise DAList$data
  get_expr_for_contrast <- function(DAList, contrast) {
    # Per-contrast matrix (e.g., filtered/imputed + normalized)
    if (!is.null(DAList$data_per_contrast) &&
        !is.null(DAList$data_per_contrast[[contrast]])) {
      return(DAList$data_per_contrast[[contrast]])
    }
    # Fallback: global data matrix
    if (!is.null(DAList$data)) {
      return(DAList$data)
    }
    stop("No expression data found for contrast '", contrast,
         "': neither DAList$data_per_contrast[[contrast]] nor DAList$data is available.")
  }
  
  # Determine contrast names
  if (!is.null(DAList$filtered_proteins_per_contrast)) {
    contrast_names <- names(DAList$filtered_proteins_per_contrast)
  } else if (!is.null(DAList$design$contrast_vector)) {
    contrast_names <- stringr::str_trim(
      stringr::str_split_fixed(DAList$design$contrast_vector, "=", 2)[, 1]
    )
    cli::cli_alert_info("No filtered_proteins_per_contrast found. Using all proteins for each contrast.")
  } else if (!is.null(contrasts_file)) {
    contrast_table <- utils::read.csv(contrasts_file, header = FALSE, stringsAsFactors = FALSE)
    contrast_vector <- contrast_table[[1]]
    contrast_names <- stringr::str_trim(
      stringr::str_split_fixed(contrast_vector, "=", 2)[, 1]
    )
    cli::cli_alert_info("Loaded contrasts from file: {.file {contrasts_file}}")
  } else {
    stop("Cannot determine contrast names: no filtered_proteins_per_contrast, contrast_vector, or contrasts_file provided.")
  }
  
  # Ensure tags lists exist
  if (is.null(DAList$tags$movingSDs))      DAList$tags$movingSDs      <- list()
  if (is.null(DAList$tags$logFC_z_scores)) DAList$tags$logFC_z_scores <- list()
  
  # Auto-select binsize using chosen expression matrix per contrast
  if (identical(binsize, "auto")) {
    message("Evaluating candidate bin sizes: ", paste(binsize_range, collapse = ", "))
    
    evaluate_binsize <- function(b_size) {
      cv_list <- c()
      
      for (contrast in contrast_names) {
        df_contrast <- DAList$results[[contrast]]
        if (is.null(df_contrast) || !"logFC" %in% names(df_contrast)) next
        
        prot_ids <- rownames(df_contrast)
        if (length(prot_ids) == 0) next
        
        expr_mat <- get_expr_for_contrast(DAList, contrast)
        # Intersect to avoid mismatches
        prot_ids_use <- intersect(prot_ids, rownames(expr_mat))
        if (length(prot_ids_use) < 2) next
        
        mean_intens <- rowMeans(expr_mat[prot_ids_use, , drop = FALSE], na.rm = TRUE)
        ord_idx <- order(mean_intens)
        logFC_sorted <- df_contrast$logFC[match(prot_ids_use, prot_ids)][ord_idx]
        
        sds <- rolling_sd(logFC_sorted, b_size)
        
        sds_mean <- mean(sds, na.rm = TRUE)
        sds_sd   <- stats::sd(sds, na.rm = TRUE)
        
        if (is.finite(sds_mean) && sds_mean != 0) {
          cv_list <- c(cv_list, sds_sd / sds_mean)
        }
      }
      
      cv_list <- cv_list[is.finite(cv_list)]
      if (length(cv_list) == 0) return(Inf)
      mean(cv_list)
    }
    
    cv_values <- vapply(binsize_range, evaluate_binsize, FUN.VALUE = numeric(1))
    binsize   <- binsize_range[which.min(cv_values)]
    message("Auto-selected binsize: ", binsize)
  }
  
  # Main loop over contrasts
  for (contrast in contrast_names) {
    df_contrast <- DAList$results[[contrast]]
    if (is.null(df_contrast) || !"logFC" %in% names(df_contrast)) next
    
    prot_ids <- rownames(df_contrast)
    
    if (length(prot_ids) == 0) {
      df_contrast$movingSD        <- numeric(0)
      df_contrast$logFC_z_scores  <- numeric(0)
      DAList$results[[contrast]]  <- df_contrast
      DAList$tags$movingSDs[[contrast]]       <- setNames(numeric(0), character(0))
      DAList$tags$logFC_z_scores[[contrast]]  <- setNames(numeric(0), character(0))
      next
    }
    
    expr_mat <- get_expr_for_contrast(DAList, contrast)
    
    # Keep only proteins present in the expression matrix
    prot_ids_use <- intersect(prot_ids, rownames(expr_mat))
    if (length(prot_ids_use) == 0) {
      cli::cli_alert_warning(
        "No overlapping proteins between results and expression matrix for contrast '{contrast}'. Skipping."
      )
      next
    }
    
    # Sort by mean intensity (per-contrast if available, otherwise global)
    mean_intens <- rowMeans(expr_mat[prot_ids_use, , drop = FALSE], na.rm = TRUE)
    ord_idx <- order(mean_intens)
    
    # Match logFC to the same subset/order
    logFC_use    <- df_contrast$logFC[match(prot_ids_use, prot_ids)]
    logFC_sorted <- logFC_use[ord_idx]
    
    # Rolling SD on sorted logFC
    moving_sds_sorted <- rolling_sd(logFC_sorted, binsize)
    
    # Reorder back to prot_ids_use order
    reorder_back <- order(ord_idx)
    moving_sds_prot <- moving_sds_sorted[reorder_back]
    names(moving_sds_prot) <- prot_ids_use
    
    # Now map back to full prot_ids (fill with NA where no expression)
    moving_sds_full <- rep(NA_real_, length(prot_ids))
    names(moving_sds_full) <- prot_ids
    moving_sds_full[match(prot_ids_use, prot_ids)] <- moving_sds_prot
    
    # Compute z-scores AFTER movingSD is final
    logFC_z_scores <- df_contrast$logFC / moving_sds_full
    names(logFC_z_scores) <- prot_ids
    
    # Store in data frame
    df_contrast$movingSD       <- moving_sds_full
    df_contrast$logFC_z_scores <- logFC_z_scores
    DAList$results[[contrast]] <- df_contrast
    
    # Also store in tags (as named vectors)
    DAList$tags$movingSDs[[contrast]]      <- moving_sds_full
    DAList$tags$logFC_z_scores[[contrast]] <- logFC_z_scores
    
    # Optional plot (uses sorted curve for smoother viz)
    # What this adds:
    #   mean_intens_sorted: mean log2 intensity for each protein, sorted in the same order as the rolling SD curve.
    # A second ggplot (p_int) with:
    #   x-axis = mean_log2_intensity (your normalized log2 intensities),
    #   y-axis = movingSD.
    # Optional saving of both plots to DAList$QC_dir as:
    #   movingSD_<contrast>_index.png
    #   movingSD_<contrast>_intensity.png
    # you’ll get:
    #   A “shape” view (vs rank index).
    #   A “biological scale” view (vs log2 intensity), more comparable to the MD plots
    if (plot) {
      # 1) Existing plot: moving SD vs rank index
      df_plot_idx <- data.frame(
        Index    = seq_along(moving_sds_sorted),
        movingSD = moving_sds_sorted
      )
      p_idx <- ggplot2::ggplot(df_plot_idx, ggplot2::aes(x = Index, y = movingSD)) +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(paste("Moving SD (rank index):", contrast,
                               "(binsize =", binsize, ")"))
      print(p_idx)
      
      # 2) New plot: moving SD vs mean log2 intensity (MD-style)
      mean_intens_sorted <- mean_intens[ord_idx]
      df_plot_int <- data.frame(
        mean_log2_intensity = mean_intens_sorted,
        movingSD            = moving_sds_sorted
      )
      p_int <- ggplot2::ggplot(df_plot_int,
                               ggplot2::aes(x = mean_log2_intensity, y = movingSD)) +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Mean log2 normalized intensity") +
        ggplot2::ggtitle(paste("Moving SD vs intensity:", contrast,
                               "(binsize =", binsize, ")"))
      print(p_int)
      
      # Optional: save both if QC_dir is present
      if (!is.null(DAList$QC_dir)) {
        dir.create(DAList$QC_dir, showWarnings = FALSE, recursive = TRUE)
        ggplot2::ggsave(
          filename = file.path(DAList$QC_dir,
                               paste0("movingSD_", contrast, "_index.png")),
          plot = p_idx, width = 7, height = 5, dpi = 300
        )
        ggplot2::ggsave(
          filename = file.path(DAList$QC_dir,
                               paste0("movingSD_", contrast, "_intensity.png")),
          plot = p_int, width = 7, height = 5, dpi = 300
        )
      }
    }
  }
  
  return(DAList)
}
