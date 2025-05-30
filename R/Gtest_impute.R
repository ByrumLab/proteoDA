#' G-test for Detection Imbalance and Conditional Imputation
#'
#' Performs a G-test (likelihood ratio test) on per-protein detection counts across two groups
#' based on observed and expected detection values. Missing values are imputed with the group-wise
#' minimum if either (1) the p-value is below a defined threshold, or (2) one group has all missing values.
#'
#' @param data_matrix A numeric matrix of protein intensities (rows = proteins, columns = samples).
#' @param group_labels A character or factor vector specifying the group for each sample (e.g., "N" or "P").
#' @param p_threshold A numeric p-value threshold for applying imputation (default = 0.05).
#' @param return_imputed_flags Logical; if TRUE, return a logical vector indicating which proteins were imputed.
#'
#' @return A list with:
#' \describe{
#'   \item{p_values}{Named numeric vector of G-test p-values per protein.}
#'   \item{imputed_data}{Matrix of the same dimension as `data_matrix` with imputations applied.}
#'   \item{imputed_flags}{(Optional) Logical vector indicating which proteins were imputed.}
#' }
#'
#' @importFrom stats pchisq
#' @export

GTest_impute <- function(data_matrix, group_labels, p_threshold = 0.01, return_imputed_flags = TRUE) {
  if (length(group_labels) != ncol(data_matrix)) {
    stop("Length of group_labels must match number of columns in data_matrix")
  }
  
  groups <- unique(group_labels)
  if (length(groups) != 2) {
    stop("Exactly two groups are required")
  }
  
  n_N <- sum(group_labels == groups[1])
  n_P <- sum(group_labels == groups[2])
  total_samples <- n_N + n_P
  
  p_values <- numeric(nrow(data_matrix))
  names(p_values) <- rownames(data_matrix)
  imputed_data <- data_matrix
  imputed_flags <- logical(nrow(data_matrix))
  imputed_reasons <- rep(NA_character_, nrow(data_matrix))
  
  for (i in 1:nrow(data_matrix)) {
    row_vals <- data_matrix[i, ]
    group_N_vals <- row_vals[group_labels == groups[1]]
    group_P_vals <- row_vals[group_labels == groups[2]]
    
    observed_N <- sum(!is.na(group_N_vals))
    observed_P <- sum(!is.na(group_P_vals))
    total_detected <- observed_N + observed_P
    
    if (total_detected == 0) {
      p_values[i] <- NA
      next
    }
    
    expected_N <- (total_detected / total_samples) * n_N
    expected_P <- (total_detected / total_samples) * n_P
    
    G_N <- if (observed_N == 0 || expected_N == 0) 0 else 2 * observed_N * log(observed_N / expected_N)
    G_P <- if (observed_P == 0 || expected_P == 0) 0 else 2 * observed_P * log(observed_P / expected_P)
    
    total_G <- G_N + G_P
    p_val <- pchisq(total_G, df = 1, lower.tail = FALSE)
    p_values[i] <- p_val
    
    imputed_flags[i] <- FALSE
    
    if ((p_val < p_threshold) || all(is.na(group_N_vals)) || all(is.na(group_P_vals))) {
      reason <- if (all(is.na(group_N_vals)) || all(is.na(group_P_vals))) {
        "all_missing_in_group"
      } else if (p_val < p_threshold) {
        "p_value"
      } else {
        NA_character_
      }
      
      missing_idx <- which(is.na(row_vals))
      if (length(missing_idx) > 0) {
        min_val <- suppressWarnings(min(row_vals, na.rm = TRUE))
        if (is.infinite(min_val)) {
          global_min <- suppressWarnings(min(data_matrix, na.rm = TRUE))
          min_val <- if (is.infinite(global_min)) 0 else global_min
        }
        imputed_data[i, missing_idx] <- min_val
        imputed_flags[i] <- TRUE
        imputed_reasons[i] <- reason
      }
    }
  }
  
  result <- list(
    p_values = p_values,
    imputed_data = imputed_data
  )
  if (return_imputed_flags) {
    result$imputed_flags <- imputed_flags
    result$imputed_reasons <- imputed_reasons
  }
  return(result)
}

#' Run G-test-based imputation on DAList data or contrast-specific data
#'
#' @param DAList A DAList object after filtering.
#' @param contrast Optional character string specifying contrast name (for use with data_per_contrast). If NULL, loops through all.
#' @param grouping_column Name of column in metadata to define group membership (default = "group").
#' @param p_threshold P-value threshold for triggering imputation (default = 0.05).
#'
#' @return The same DAList with imputed data written into $data or $data_per_contrast.
#' @export

impute_missing_by_gtest <- function(DAList,
                                    contrast = NULL,
                                    grouping_column = "group",
                                    p_threshold = 0.01) {
  if (!is.null(contrast)) {
    if (length(contrast) == 1 && file.exists(contrast)) {
      contrast_table <- read.csv(contrast, header = FALSE, stringsAsFactors = FALSE)
      contrast_vector <- contrast_table[[1]]
      contrast_list <- trimws(sub("=.*", "", contrast_vector))
    } else {
      contrast_list <- contrast
    }
  } else if (!is.null(DAList$data_per_contrast)) {
    contrast_list <- names(DAList$data_per_contrast)
  } else {
    contrast_list <- NULL
  }
  
  if (!is.null(contrast_list)) {
    for (ctr in contrast_list) {
      contrast_groups <- unlist(strsplit(ctr, "_vs_"))
      sample_ids <- rownames(DAList$metadata)[DAList$metadata[[grouping_column]] %in% contrast_groups]
      sample_meta <- DAList$metadata[sample_ids, , drop = FALSE]
      
      use_filtered <- !is.null(DAList$filtered_proteins_per_contrast) && !is.null(DAList$filtered_proteins_per_contrast[[ctr]])
      if (use_filtered) {
        filtered_ids <- DAList$filtered_proteins_per_contrast[[ctr]]
        contrast_data_matrix <- DAList$data[filtered_ids, sample_ids, drop = FALSE]
      } else {
        filtered_ids <- rownames(DAList$data)
        contrast_data_matrix <- DAList$data[, sample_ids, drop = FALSE]
      }
      group_labels <- as.character(sample_meta[[grouping_column]])
      
      result <- GTest_impute(
        data_matrix = contrast_data_matrix,
        group_labels = group_labels,
        p_threshold = p_threshold,
        return_imputed_flags = TRUE
      )
      
      imputed_data_ordered <- result$imputed_data[match(filtered_ids, rownames(result$imputed_data)), , drop = FALSE]
      if (is.null(DAList$data_per_contrast)) DAList$data_per_contrast <- list()
      DAList$data_per_contrast[[ctr]] <- imputed_data_ordered
      
      if (!is.null(DAList$annotation)) {
        if (is.null(DAList$annotation_per_contrast)) DAList$annotation_per_contrast <- list()
        DAList$annotation_per_contrast[[ctr]] <- DAList$annotation[filtered_ids, , drop = FALSE]
      }
      
      imputed_sample_per_protein <- lapply(seq_len(nrow(contrast_data_matrix)), function(i) {
        original_row <- contrast_data_matrix[i, ]
        new_row <- result$imputed_data[i, ]
        imputed_idx <- which(is.na(original_row) & !is.na(new_row))
        colnames(contrast_data_matrix)[imputed_idx]
      })
      
      result_table <- data.frame(
        protein = rownames(contrast_data_matrix),
        GTest_p_value = result$p_values,
        imputed = result$imputed_flags,
        imputed_reason = result$imputed_reasons,
        imputed_samples = sapply(imputed_sample_per_protein, function(x) paste(x, collapse = ";")),
        stringsAsFactors = FALSE
      )
      
      summary_tag <- list(
        contrast = ctr,
        grouping_column = grouping_column,
        p_threshold = p_threshold,
        imputed_proteins = sum(result$imputed_flags)
      )
      
      result_table$contrast <- ctr
      result_table$grouping_column <- grouping_column
      result_table$p_threshold <- p_threshold
      result_table$imputed_proteins <- sum(result$imputed_flags)
      result_table$imputation_tag <- I(list(summary_tag))
      
      if (is.null(DAList$tags$gtest_imputation_results)) {
        DAList$tags$gtest_imputation_results <- list()
      }
      DAList$tags$gtest_imputation_results[[ctr]] <- result_table
      
      export_table <- result_table
      export_table$imputation_tag <- NULL
      utils::write.csv(export_table, file = paste0("gtest_imputation_", ctr, ".csv"), row.names = FALSE)
    }
  }
  
  return(DAList)
}

