#' Compute Rolling Standard Deviations and Z-scores for logFC Comparisons
#'
#' This function computes a rolling (moving) standard deviation of log fold change (logFC)
#' values for each comparison in the `raw$results` list. The SD is computed over sorted data
#' using mean intensity for sorting (without modifying the data). It also computes logFC
#' Z-scores using the rolling SDs and stores all results in `raw$tags`.
#'
#' If `binsize = "auto"`, the function evaluates a set of candidate bin sizes and chooses
#' the one that minimizes the average coefficient of variation (CV) of the moving SDs.
#'
#' If `plot = TRUE`, a plot of the final moving SDs for each comparison will be generated
#' and saved to `QC_dir` if available.
#'
#' @param raw An S3 object containing `data`, `results`, and `tags` slots. Must contain:
#'   - `raw$data`: A data frame with numeric expression values and a `Log2Fold` column.
#'   - `raw$results`: A named list of comparisons (e.g., `treat_vs_control`) each containing a `logFC` vector.
#'   - Optionally: `raw$QC_dir` as a string path for saving QC plots.
#' @param binsize Integer or `"auto"`. The number of consecutive proteins (rows) to use in each rolling window,
#'   or `"auto"` to determine the best size automatically.
#' @param binsize_range Numeric vector. Only used if `binsize = "auto"`. A vector of candidate bin sizes to test.
#' @param plot Logical. If `TRUE`, plots the final moving SDs for each comparison and saves to QC_dir.
#'
#' @return The input `raw` object with updated `tags`:
#'   - `raw$tags$movingSDs`: A list of rolling SD vectors for each comparison.
#'   - `raw$tags$logFC_z_scores`: A list of Z-score vectors (logFC divided by rolling SD).
#'
#' @examples
#' raw <- compute_movingSD_zscores(raw, binsize = "auto", binsize_range = c(50, 100, 200))
#' raw <- compute_movingSD_zscores(raw, binsize = 100, plot = FALSE)
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
    sds <- sapply(seq_len(end_idx), function(i) sd(x[i:(i + window_size - 1)], na.rm = TRUE))
    c(sds, rep(tail(sds, 1), n_x - end_idx))
  }
  
  if (!is.null(DAList$filtered_proteins_per_contrast)) {
    contrast_names <- names(DAList$filtered_proteins_per_contrast)
  } else if (!is.null(DAList$design$contrast_vector)) {
    contrast_names <- stringr::str_trim(stringr::str_split_fixed(DAList$design$contrast_vector, "=", 2)[, 1])
    cli::cli_alert_info("No filtered_proteins_per_contrast found. Using all proteins for each contrast.")
  } else if (!is.null(contrasts_file)) {
    contrast_table <- utils::read.csv(contrasts_file, header = FALSE, stringsAsFactors = FALSE)
    contrast_vector <- contrast_table[[1]]
    contrast_names <- stringr::str_trim(stringr::str_split_fixed(contrast_vector, "=", 2)[, 1])
    cli::cli_alert_info("Loaded contrasts from file: {.file {contrasts_file}}")
  } else {
    stop("Cannot determine contrast names: no filtered_proteins_per_contrast, contrast_vector, or contrasts_file provided.")
  }
  
  # Auto-select binsize
  if (binsize == "auto") {
    message("Evaluating candidate bin sizes: ", paste(binsize_range, collapse = ", "))
    evaluate_binsize <- function(b_size) {
      cv_list <- c()
      for (contrast in contrast_names) {
        df_contrast <- DAList$results[[contrast]]
        if (!"logFC" %in% names(df_contrast)) next
        
        prot_ids <- rownames(df_contrast)
        mean_intens <- rowMeans(DAList$data[prot_ids, , drop = FALSE], na.rm = TRUE)
        ord_idx <- order(mean_intens)
        logFC_sorted <- df_contrast$logFC[ord_idx]
        sds <- rolling_sd(logFC_sorted, b_size)
        
        sds_mean <- mean(sds, na.rm = TRUE)
        sds_sd <- sd(sds, na.rm = TRUE)
        
        if (is.finite(sds_mean) && sds_mean != 0) {
          cv_list <- c(cv_list, sds_sd / sds_mean)
        }
      }
      
      cv_list <- cv_list[is.finite(cv_list)]
      if (length(cv_list) == 0) return(Inf)
      mean(cv_list)
    }
    
    cv_values <- sapply(binsize_range, evaluate_binsize)
    binsize <- binsize_range[which.min(cv_values)]
    message("Auto-selected binsize: ", binsize)
  }
  
  for (contrast in contrast_names) {
    df_contrast <- DAList$results[[contrast]]
    if (!"logFC" %in% names(df_contrast)) next
    
    prot_ids <- rownames(df_contrast)
    
    if (length(prot_ids) == 0) {
      df_contrast$movingSD <- numeric(0)
      df_contrast$logFC_z_scores <- numeric(0)
      DAList$results[[contrast]] <- df_contrast
      DAList$tags$movingSDs[[contrast]] <- setNames(numeric(0), character(0))
      DAList$tags$logFC_z_scores[[contrast]] <- setNames(numeric(0), character(0))
      next
    }
    
    # Sort by mean intensity
    mean_intens <- rowMeans(DAList$data[prot_ids, , drop = FALSE], na.rm = TRUE)
    ord_idx <- order(mean_intens)
    logFC_sorted <- df_contrast$logFC[ord_idx]
    
    # Rolling SD on sorted logFC
    moving_sds_sorted <- rolling_sd(logFC_sorted, binsize)
    
    # Reorder back to original order
    reorder_back <- order(ord_idx)
    moving_sds_original <- moving_sds_sorted[reorder_back]
    names(moving_sds_original) <- prot_ids
    
    # Compute z-scores AFTER movingSD is final
    logFC_z_scores <- df_contrast$logFC / moving_sds_original
    names(logFC_z_scores) <- prot_ids
    
    # Store in data frame
    df_contrast$movingSD <- moving_sds_original
    df_contrast$logFC_z_scores <- logFC_z_scores
    DAList$results[[contrast]] <- df_contrast
    
    # Also store in tags (as named vectors)
    DAList$tags$movingSDs[[contrast]] <- moving_sds_original
    DAList$tags$logFC_z_scores[[contrast]] <- logFC_z_scores
    
    # Optional plot
    if (plot) {
      df_plot <- data.frame(Index = seq_along(moving_sds_sorted), movingSD = moving_sds_sorted)
      p <- ggplot(df_plot, aes(x = Index, y = movingSD)) +
        geom_line() +
        theme_minimal() +
        ggtitle(paste("Moving SD:", contrast, "(binsize =", binsize, ")"))
      print(p)
    }
  }
  
  return(DAList)
}

