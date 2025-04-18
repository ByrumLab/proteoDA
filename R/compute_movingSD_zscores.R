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
                                        plot = FALSE) {
  library(ggplot2)
  library(dplyr)
  
  # Helper function: rolling standard deviation over a window (window_size)
  rolling_sd <- function(x, window_size) {
    n_x <- length(x)
    end_idx <- n_x - window_size + 1
    if (end_idx < 1) {
      # If the window is too large, just return a vector of NAs
      return(rep(NA_real_, n_x))
    }
    sds <- sapply(seq_len(end_idx), function(i) {
      sd(x[i:(i + window_size - 1)], na.rm = TRUE)
    })
    # Extend the last computed value to the trailing elements
    sds_full <- c(sds, rep(tail(sds, 1), n_x - end_idx))
    return(sds_full)
  }
  
  #------------------------
  # 1) Determine optimal binsize if binsize="auto"
  #------------------------
  if (binsize == "auto") {
    message("Evaluating candidate bin sizes: ", paste(binsize_range, collapse = ", "))
    
    evaluate_binsize <- function(b_size) {
      cv_list <- c()
      for (contrast in names(DAList$results)) {
        
        # Pull the data for that contrast
        df_contrast <- DAList$results[[contrast]]
        if (!"logFC" %in% names(df_contrast)) {
          next
        }
        
        # Get filtered protein IDs and their logFC
        logFC_filtered <- df_contrast$logFC
        prot_ids <- rownames(df_contrast)
        
        # Calculate mean intensity for these filtered proteins only
        mean_intens <- rowMeans(DAList$data[prot_ids, , drop = FALSE], na.rm = TRUE)
        
        # Sort by ascending mean intensity
        ord_idx <- order(mean_intens, decreasing = FALSE)
        logFC_sorted <- logFC_filtered[ord_idx]
        
        # Compute rolling SD
        sds <- rolling_sd(logFC_sorted, window_size = b_size)
        
        # Compute coefficient of variation (CV)
        sds_mean <- mean(sds, na.rm = TRUE)
        sds_sd   <- sd(sds, na.rm = TRUE)
        cv       <- sds_sd / sds_mean
        cv_list  <- c(cv_list, cv)
      }
      mean(cv_list, na.rm = TRUE)
    }
    
    cv_values <- sapply(binsize_range, evaluate_binsize)
    binsize <- binsize_range[which.min(cv_values)]
    message("Auto-selected binsize: ", binsize)
  }
  
  #------------------------
  # 2) Compute rollingSD and logFC_z_scores for each contrast (filtered proteins only)
  #------------------------
  for (contrast in names(DAList$results)) {
    
    df_contrast <- DAList$results[[contrast]]
    
    if (!"logFC" %in% names(df_contrast)) {
      warning("Contrast ", contrast, " has no 'logFC' column; skipping.")
      next
    }
    
    # Filtered protein IDs and their logFC
    prot_ids <- rownames(df_contrast)
    logFC_filtered <- df_contrast$logFC
    
    # Nothing to compute if no proteins
    if (length(prot_ids) == 0) {
      df_contrast$movingSD <- numeric(0)
      df_contrast$logFC_z_scores <- numeric(0)
      DAList$results[[contrast]] <- df_contrast
      next
    }
    
    # Sort subset by mean intensity
    mean_intens <- rowMeans(DAList$data[prot_ids, , drop = FALSE], na.rm = TRUE)
    ord_idx <- order(mean_intens, decreasing = FALSE)
    logFC_sorted <- logFC_filtered[ord_idx]
    
    # Rolling SD in sorted order
    moving_sds_subset <- rolling_sd(logFC_sorted, binsize)
    
    # Map it back to the original (filtered) order in df_contrast
    reorder_back <- order(ord_idx)
    moving_sds_original <- moving_sds_subset[reorder_back]
    
    # Compute Z-scores
    z_scores_subset <- logFC_filtered / moving_sds_original
    
    # Add these columns (filtered proteins only) to the existing df
    df_contrast$movingSD       <- moving_sds_original
    df_contrast$logFC_z_scores <- z_scores_subset
    
    # Replace DAList$results[[contrast]] with the updated df
    DAList$results[[contrast]] <- df_contrast
    
    # (Optional) Plot if requested
    if (plot) {
      df_plot <- data.frame(
        Index = seq_along(moving_sds_subset),
        movingSD = moving_sds_subset
      )
      p <- ggplot(df_plot, aes(x = Index, y = movingSD)) +
        geom_line() +
        theme_minimal() +
        ggtitle(paste("Moving SD:", contrast, "(binsize =", binsize, ")"))
      print(p)
      # Optionally, ggsave(...) if you want to save the plot
    }
  }
  
  return(DAList)
}
