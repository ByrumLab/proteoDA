#' Analyze Within-Group Log2FC Standard Deviations
#'
#' Calculates pairwise log2 fold change (log2FC) standard deviations within each sample group,
#' generates diagnostic density plots, flags groups with outlier SDs, and exports a summary table.
#'
#' @param count_data A data frame or matrix of expression data (rows = features, columns = samples).
#' @param sample_metadata A data frame with at least two columns: `sample` (matching column names of `count_data`) 
#'   and a grouping column (e.g., "group") that specifies sample groupings.
#' @param group_column Character string specifying the column name in `sample_metadata` used to group samples. Default is `"group"`.
#' @param output_dir Directory to save plots and summary table. Created if it doesn't exist. Default is `"QC_dir"`.
#' @param outlier_method Method for identifying SD outliers: `"IQR"` (default) or `"z-score"`.
#'   - `"IQR"` flags values outside 1.5 × IQR from the 1st and 3rd quartiles.
#'   - `"z-score"` flags values with Z-scores exceeding `z_thresh`.
#' @param z_thresh Z-score threshold for flagging outliers. Only used when `outlier_method = "z-score"`. Default is `3`.
#'
#' @return Invisibly returns a data frame with per-group statistics:
#'   \item{Group}{Group name}
#'   \item{Mean_SD}{Mean of pairwise log2FC SDs}
#'   \item{SD_of_SDs}{Standard deviation of SDs across log2FC pairs}
#'   \item{N_Pairs}{Number of pairwise comparisons}
#'   \item{Outlier}{Logical flag indicating SD outliers based on the selected method}
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   analyze_within_group_SDs(my_data$data, my_data$metadata, outlier_method = "z-score", z_thresh = 2.5)
#' }
analyze_within_group_SDs <- function(count_data,
                                     sample_metadata,
                                     group_column = "group",
                                     output_dir = "QC_dir",
                                     outlier_method = c("IQR", "z-score"),
                                     z_thresh = 3) {
  outlier_method <- match.arg(outlier_method)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  group_split <- split(sample_metadata$sample, sample_metadata[[group_column]])
  all_group_means <- numeric()
  
  stats_list <- list()
  
  for (grp in names(group_split)) {
    samples <- group_split[[grp]]
    data <- count_data[, samples, drop = FALSE]
    n <- ncol(data)
    
    if (n < 2) next
    
    #log2fc <- data.frame()
    log2fc <- as.data.frame(matrix(nrow = nrow(data), ncol = 0))
    
    
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        log2fc[[paste0(samples[j], "_vs_", samples[i])]] <- data[[samples[j]]] - data[[samples[i]]]
      }
    }
    
    sd_vals <- apply(log2fc, 2, sd, na.rm = TRUE)
    sd_mean <- mean(sd_vals, na.rm = TRUE)
    sd_sd <- sd(sd_vals, na.rm = TRUE)
    
    stats_list[[grp]] <- list(
      Mean_SD = sd_mean,
      SD_of_SDs = sd_sd,
      N_Pairs = length(sd_vals)
    )
    all_group_means <- c(all_group_means, sd_mean)
    
    p <- ggplot() +
      stat_function(fun = dnorm, args = list(mean = 0, sd = sd_mean), n = 1000) +
      theme_bw() +
      labs(title = paste("SD Distribution -", grp), x = "log2FC mean", y = "density") +
      xlim(-3, 3) +
      annotate("text", x = sd_mean, y = 1, label = paste0("SD = ", round(sd_mean, 5)), hjust = 1.1, vjust = 1) +
      annotate("text", x = 2 * sd_mean, y = 1, label = paste0("2SD = ", round(2 * sd_mean, 5)), hjust = -0.1, vjust = 1) +
      geom_vline(xintercept = sd_mean, lty = 4, col = "blue", lwd = 0.4) +
      geom_vline(xintercept = 2 * sd_mean, lty = 4, col = "blue", lwd = 0.4)
    
    ggsave(filename = file.path(output_dir, paste0("SD_dist_", grp, ".pdf")), plot = p, width = 6, height = 4)
  }
  
  outlier_flags <- logical(length(stats_list))
  names(outlier_flags) <- names(stats_list)
  
  if (length(stats_list) >= 3) {
    if (outlier_method == "IQR") {
      Q1 <- quantile(all_group_means, 0.25)
      Q3 <- quantile(all_group_means, 0.75)
      IQR <- Q3 - Q1
      for (grp in names(stats_list)) {
        val <- stats_list[[grp]]$Mean_SD
        outlier_flags[grp] <- val < (Q1 - 1.5 * IQR) | val > (Q3 + 1.5 * IQR)
      }
    } else if (outlier_method == "z-score") {
      mu <- mean(all_group_means)
      sigma <- sd(all_group_means)
      for (grp in names(stats_list)) {
        z <- (stats_list[[grp]]$Mean_SD - mu) / sigma
        outlier_flags[grp] <- abs(z) > z_thresh
      }
    }
  }
  
  sd_summary <- do.call(rbind, lapply(names(stats_list), function(grp) {
    data.frame(
      Group = grp,
      Mean_SD = stats_list[[grp]]$Mean_SD,
      SD_of_SDs = stats_list[[grp]]$SD_of_SDs,
      N_Pairs = stats_list[[grp]]$N_Pairs,
      Outlier = outlier_flags[grp],
      stringsAsFactors = FALSE
    )
  }))
  
  write.csv(sd_summary, file = file.path(output_dir, "within_group_SD_summary.csv"), row.names = FALSE)
  
  invisible(sd_summary)
}
