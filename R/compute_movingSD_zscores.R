#' Compute Rolling Standard Deviations and Z-scores for logFC Comparisons
#'
#' This function computes a rolling (moving) standard deviation of log fold change (logFC)
#' values for each comparison in the \code{DAList$results} list. For each contrast, proteins
#' are ordered by mean log2 intensity (without modifying the underlying data matrix).
#' A moving standard deviation of logFC is computed along this ordering and used to derive
#' logFC Z-scores. Results are stored both in \code{DAList$results[[contrast]]} and
#' \code{DAList$tags}.
#'
#' If \code{binsize = "auto"}, the function derives candidate bin sizes from the protein
#' count, evaluates their stability using the coefficient of variation (CV) of the moving SD
#' curves across contrasts, and chooses a binsize as the smallest value whose CV is within
#' 5\% of the minimum CV. Candidate binsizes are constrained between a minimum number
#' of proteins (40) and a maximum of 20\% of the total protein count.
#'
#' If \code{plot = TRUE}, two plots are generated per contrast:
#' \enumerate{
#'   \item Moving SD vs rank index (proteins sorted by intensity).
#'   \item Moving SD vs mean log2 intensity (MD-style).
#' }
#' Plots are printed to the current device and, if \code{DAList$QC_dir} is set,
#' saved as PNG files there. When \code{binsize = "auto"}, a global
#' \code{movingSD_binsize_CV.png} (CV vs binsize) is also saved in \code{QC_dir}.
#'
#' @param DAList A DAList-like object containing at least:
#'   \itemize{
#'     \item \code{data}: A matrix/data frame with numeric expression values
#'       (log2 intensities; proteins in rows, samples in columns).
#'     \item \code{results}: A named list of per-contrast data frames, each
#'       containing a \code{logFC} column and rownames as protein IDs.
#'     \item Optionally, \code{data_per_contrast}: a named list of per-contrast
#'       matrices/data frames. When present, \code{data_per_contrast[[contrast]]}
#'       is used for that contrast instead of the global \code{data}.
#'     \item Optionally, \code{QC_dir}: a string path for saving QC plots.
#'     \item Optionally, \code{design$contrast_vector} or
#'       \code{filtered_proteins_per_contrast} to infer contrast names.
#'   }
#'
#' @param binsize Integer or \code{"auto"}. The number of consecutive proteins (rows)
#'   to use in each rolling window, or \code{"auto"} to determine the best size
#'   automatically from the actual protein count.
#'
#' @param plot Logical. If \code{TRUE}, plots the final moving SDs for each
#'   contrast and (if \code{DAList$QC_dir} is available) saves them to disk.
#'
#' @param contrasts_file Optional CSV file with contrasts (one per row, \code{"name = formula"}
#'   syntax) used only to determine contrast names if they are not already
#'   available in the \code{DAList}.
#'
#' @return The input \code{DAList} object with updated \code{tags}:
#'   \itemize{
#'     \item \code{DAList$tags$movingSDs[[contrast]]}: rolling SD vector for each contrast.
#'     \item \code{DAList$tags$logFC_z_scores[[contrast]]}: Z-score vector for each contrast.
#'     \item \code{DAList$results[[contrast]]$movingSD}: moving SD per protein.
#'     \item \code{DAList$results[[contrast]]$logFC_z_scores}: logFC Z-score per protein.
#'   }
#'
#' @export
compute_movingSD_zscores <- function(DAList,
                                     binsize = "auto",
                                     plot = FALSE,
                                     contrasts_file = NULL) {
  
  ## --- rolling SD helper -----------------------------------------------------
  rolling_sd <- function(x, window_size) {
    n_x    <- length(x)
    end_ix <- n_x - window_size + 1
    if (end_ix < 1) return(rep(NA_real_, n_x))
    sds <- vapply(
      X   = seq_len(end_ix),
      FUN = function(i) stats::sd(x[i:(i + window_size - 1)], na.rm = TRUE),
      FUN.VALUE = numeric(1)
    )
    c(sds, rep(tail(sds, 1), n_x - end_ix))
  }
  
  ## --- expression matrix per contrast ----------------------------------------
  # Prefer DAList$data_per_contrast[[contrast]] if available; otherwise DAList$data
  get_expr_for_contrast <- function(DAList, contrast) {
    if (!is.null(DAList$data_per_contrast) &&
        !is.null(DAList$data_per_contrast[[contrast]])) {
      return(DAList$data_per_contrast[[contrast]])
    }
    if (!is.null(DAList$data)) {
      return(DAList$data)
    }
    stop("No expression data found for contrast '", contrast,
         "': neither DAList$data_per_contrast[[contrast]] nor DAList$data is available.")
  }
  
  ## --- determine contrast names ----------------------------------------------
  if (!is.null(DAList$filtered_proteins_per_contrast)) {
    contrast_names <- names(DAList$filtered_proteins_per_contrast)
  } else if (!is.null(DAList$design$contrast_vector)) {
    contrast_names <- stringr::str_trim(
      stringr::str_split_fixed(DAList$design$contrast_vector, "=", 2)[, 1]
    )
    cli::cli_alert_info("No filtered_proteins_per_contrast found. Using all proteins for each contrast.")
  } else if (!is.null(contrasts_file)) {
    contrast_table  <- utils::read.csv(contrasts_file, header = FALSE, stringsAsFactors = FALSE)
    contrast_vector <- contrast_table[[1]]
    contrast_names  <- stringr::str_trim(
      stringr::str_split_fixed(contrast_vector, "=", 2)[, 1]
    )
    cli::cli_alert_info("Loaded contrasts from file: {.file {contrasts_file}}")
  } else {
    stop("Cannot determine contrast names: no filtered_proteins_per_contrast, contrast_vector, or contrasts_file provided.")
  }
  
  ## --- ensure tags lists exist ----------------------------------------------
  if (is.null(DAList$tags)) DAList$tags <- list()
  if (is.null(DAList$tags$movingSDs))      DAList$tags$movingSDs      <- list()
  if (is.null(DAList$tags$logFC_z_scores)) DAList$tags$logFC_z_scores <- list()
  
  ## --- auto binsize selection -----------------------------------------------
  if (identical(binsize, "auto")) {
    # Max number of proteins across contrasts (based on results tables)
    n_per_contrast <- vapply(
      contrast_names,
      FUN = function(ct) {
        df_ct <- DAList$results[[ct]]
        if (is.null(df_ct)) return(0L)
        length(rownames(df_ct))
      },
      FUN.VALUE = integer(1)
    )
    n_proteins_max <- max(n_per_contrast, na.rm = TRUE)
    if (n_proteins_max <= 0L) {
      stop("No proteins found in DAList$results for any contrast.")
    }
    
    # Candidate binsizes from fractions of total N, constrained between
    # a minimum of 40 proteins and a maximum of 20% of N.
    frac_candidates    <- c(0.01, 0.02, 0.05, 0.10, 0.15, 0.20)
    binsize_candidates <- unique(round(n_proteins_max * frac_candidates))
    
    min_binsize <- 40L
    max_binsize <- max(min_binsize, floor(n_proteins_max / 5))  # 20% cap
    
    binsize_candidates <- binsize_candidates[
      binsize_candidates >= min_binsize & binsize_candidates <= max_binsize
    ]
    binsize_candidates <- sort(unique(as.integer(binsize_candidates)))
    
    if (length(binsize_candidates) == 0L) {
      binsize_candidates <- as.integer(max(min_binsize, floor(n_proteins_max / 10)))
    }
    
    cli::cli_alert_info(
      "Auto-selecting binsize from candidates: {paste(binsize_candidates, collapse = ', ')} (n_proteins_max = {n_proteins_max})"
    )
    
    # Use the same working logic to evaluate CV at each candidate binsize
    evaluate_binsize <- function(b_size) {
      cv_list <- c()
      
      for (contrast in contrast_names) {
        df_contrast <- DAList$results[[contrast]]
        if (is.null(df_contrast) || !"logFC" %in% names(df_contrast)) next
        
        prot_ids <- rownames(df_contrast)
        if (length(prot_ids) == 0) next
        
        expr_mat <- get_expr_for_contrast(DAList, contrast)
        prot_ids_use <- intersect(prot_ids, rownames(expr_mat))
        if (length(prot_ids_use) < b_size) next
        
        mean_intens <- rowMeans(expr_mat[prot_ids_use, , drop = FALSE], na.rm = TRUE)
        ord_idx     <- order(mean_intens)
        
        logFC_use    <- df_contrast$logFC[match(prot_ids_use, prot_ids)]
        logFC_sorted <- logFC_use[ord_idx]
        
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
    
    cv_values <- vapply(binsize_candidates, evaluate_binsize, FUN.VALUE = numeric(1))
    
    finite_idx <- which(is.finite(cv_values))
    if (length(finite_idx) == 0L) {
      cli::cli_alert_warning(
        "All CV values are infinite/NA during binsize auto-selection. Using median candidate."
      )
      binsize <- as.integer(round(stats::median(binsize_candidates)))
    } else {
      min_cv <- min(cv_values[finite_idx])
      tol    <- 0.05  # 5% tolerance
      
      good_idx <- which(cv_values <= min_cv * (1 + tol))
      if (length(good_idx) == 0L) {
        chosen_idx <- which.min(cv_values)
      } else {
        chosen_idx <- min(good_idx)  # smallest binsize with near-optimal CV
      }
      binsize <- as.integer(binsize_candidates[chosen_idx])
      
      message("Auto-selected binsize: ", binsize,
              " (CV=", signif(cv_values[chosen_idx], 3),
              ", min CV=", signif(min_cv, 3), ")")
    }
    
    # CV vs binsize QC plot (only finite CVs)
    cv_df <- data.frame(
      binsize = binsize_candidates,
      CV      = as.numeric(cv_values)
    )
    cv_df_finite <- cv_df[is.finite(cv_df$CV), , drop = FALSE]
    
    if (nrow(cv_df_finite) > 0L) {
      p_cv <- ggplot2::ggplot(cv_df_finite,
                              ggplot2::aes(x = binsize, y = CV)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle("CV of moving SD vs binsize (auto-selection)") +
        ggplot2::xlab("Binsize (proteins)") +
        ggplot2::ylab("Mean CV across contrasts")
      
      if (plot) {
        print(p_cv)
      }
      if (!is.null(DAList$QC_dir)) {
        dir.create(DAList$QC_dir, showWarnings = FALSE, recursive = TRUE)
        ggplot2::ggsave(
          filename = file.path(DAList$QC_dir, "movingSD_binsize_CV.png"),
          plot = p_cv, width = 7, height = 5, dpi = 300
        )
      }
    } else {
      cli::cli_alert_warning(
        "Could not compute finite CVs for any candidate binsize; skipping CV vs binsize plot."
      )
    }
  } else {
    binsize <- as.integer(round(binsize))
  }
  
  ## --- main loop over contrasts (unchanged) ----------------------------------
  for (contrast in contrast_names) {
    df_contrast <- DAList$results[[contrast]]
    if (is.null(df_contrast) || !"logFC" %in% names(df_contrast)) next
    
    prot_ids <- rownames(df_contrast)
    if (length(prot_ids) == 0) {
      df_contrast$movingSD       <- numeric(0)
      df_contrast$logFC_z_scores <- numeric(0)
      DAList$results[[contrast]] <- df_contrast
      DAList$tags$movingSDs[[contrast]]      <- setNames(numeric(0), character(0))
      DAList$tags$logFC_z_scores[[contrast]] <- setNames(numeric(0), character(0))
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
    ord_idx     <- order(mean_intens)
    
    # Match logFC to the same subset/order
    logFC_use    <- df_contrast$logFC[match(prot_ids_use, prot_ids)]
    logFC_sorted <- logFC_use[ord_idx]
    
    # Rolling SD on sorted logFC
    moving_sds_sorted <- rolling_sd(logFC_sorted, binsize)
    
    # Reorder back to prot_ids_use order
    reorder_back    <- order(ord_idx)
    moving_sds_prot <- moving_sds_sorted[reorder_back]
    names(moving_sds_prot) <- prot_ids_use
    
    # Map back to full prot_ids (fill with NA where no expression)
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
    
    # Optional plots (same as working version)
    if (plot) {
      # 1) Moving SD vs rank index
      df_plot_idx <- data.frame(
        Index    = seq_along(moving_sds_sorted),
        movingSD = moving_sds_sorted
      )
      p_idx <- ggplot2::ggplot(df_plot_idx,
                               ggplot2::aes(x = Index, y = movingSD)) +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::ggtitle(paste("Moving SD (rank index):", contrast,
                               "(binsize =", binsize, ")"))
      print(p_idx)
      
      # 2) Moving SD vs mean log2 intensity (MD-style)
      mean_intens_sorted <- mean_intens[ord_idx]
      df_plot_int <- data.frame(
        mean_log2_intensity = mean_intens_sorted,
        movingSD            = moving_sds_sorted
      )
      p_int <- ggplot2::ggplot(df_plot_int,
                               ggplot2::aes(x = mean_log2_intensity,
                                            y = movingSD)) +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::xlab("Mean log2 normalized intensity") +
        ggplot2::ggtitle(paste("Moving SD vs intensity:", contrast,
                               "(binsize =", binsize, ")"))
      print(p_int)
      
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
