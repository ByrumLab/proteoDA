#' Compute Rolling Standard Deviations and Z-scores for logFC Comparisons
#'
#' This function computes a rolling (moving) standard deviation of log fold change (logFC)
#' values for each comparison in the `results$results` list. The SD is computed over sorted data,
#' using mean intensity for sorting (without modifying the data). It also computes logFC
#' Z-scores using the rolling SDs and stores all results in `results$tags`.
#'
#' @param results An S3 object containing `data`, `results`, and `tags` slots. Must contain:
#'   - `results$data`: A data frame with numeric expression values and a `Log2Fold` column.
#'   - `results$results`: A named list of comparisons (e.g., `treat_vs_control`) each containing a `logFC` vector.
#' @param binsize Integer. The number of consecutive proteins (rows) to use in each rolling window.
#'
#' @return The input `results` object with updated `tags`:
#'   - `results$tags$movingSDs`: A list of rolling SD vectors for each comparison.
#'   - `results$tags$logFC_z_scores`: A list of Z-score vectors (logFC divided by rolling SD).
#'
#' @examples
#' results <- compute_moving_sd_and_zscores(results, binsize = 100)
#'
#' @export
compute_movingSD_zscores <- function(results, binsize) {
  # Compute temporary mean intensity for sorting
  mean_intensity <- rowMeans(results$data[, -ncol(results$data)], na.rm = TRUE)
  order_index <- order(mean_intensity, decreasing = FALSE)
  n <- nrow(results$data)
  
  # Initialize containers
  results$tags$movingSDs <- list()
  results$tags$logFC_z_scores <- list()
  
  # Loop through each comparison
  for (comp in names(results$results)) {
    # Extract and sort logFC values
    logFC <- results$results[[comp]]$logFC
    logFC_sorted <- logFC[order_index]
    
    # Rolling SD calculation
    end <- n - binsize + 1
    moving_sd <- sapply(1:end, function(i) {
      sd(logFC_sorted[i:(i + binsize - 1)], na.rm = TRUE)
    })
    
    # Pad remaining values
    moving_sd_remain <- rep(tail(moving_sd, 1), n - end)
    moving_sd_full <- c(moving_sd, moving_sd_remain)
    
    # Reorder to original order
    moving_sd_original_order <- numeric(n)
    moving_sd_original_order[order_index] <- moving_sd_full
    
    # Store SD and Z-score
    results$tags$movingSDs[[comp]] <- moving_sd_original_order
    results$tags$logFC_z_scores[[comp]] <- logFC / moving_sd_original_order
  }
  
  return(results)
}
