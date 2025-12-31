## ------------------------------------------------------------------
## Default columns and filtering thresholds for limma tables / Excel
## ------------------------------------------------------------------

# Columns to use in the main DA table sheets (if present in results/annotation)
#DA_table_cols <- c("uniprot_id","Accession.Number","Protein.Description")  # DIANN

# Columns that are purely statistical results (subset of above)
#stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSD", "logFC_z_scores", "sig.PVal", "sig.FDR")

# Default filter thresholds used for the 'filter info' summary in Excel
filt_min_reps   <- NA_integer_
filt_min_groups <- NA_integer_


#' Summarize the number of DA proteins in a contrast
#'
#' Internal function to summarize the number of DE genes/proteins for a given
#' contrast.
#'
#' @param contrast_name The name of the contrast to summarize.
#' @param contrast_res_list A list of per-contrast DA results.
#'
#' @return A data frame summarizing differential abundance for the given contrast.
#' @keywords internal
summarize_contrast_DA <- function(contrast_name, contrast_res_list) {
  tmp <- contrast_res_list[[contrast_name]][,c("sig.PVal", "sig.FDR")]
  data.frame(cbind(contrast = contrast_name,
                   type = c("down", "nonsig", "up"),
                   rbind(colSums(tmp == -1, na.rm = TRUE),
                         colSums(tmp == 0, na.rm = TRUE),
                         colSums(tmp == 1, na.rm = TRUE))))
}

#' Choose the best intensity matrix for export (prefer normalized per-contrast)
#' @keywords internal
get_export_data <- function(DAList, prefer = c("data_per_contrast", "data")) {
  prefer <- match.arg(prefer)
  
  # Determine canonical sample order
  sample_ids <- colnames(DAList$data)
  if (!is.null(DAList$metadata) && all(rownames(DAList$metadata) %in% sample_ids)) {
    # keep DAList$data column order as the canonical order (safer)
    sample_ids <- colnames(DAList$data)
  }
  
  # Determine full protein universe (best = annotation)
  all_proteins <- NULL
  if (!is.null(DAList$annotation) && !is.null(rownames(DAList$annotation))) {
    all_proteins <- rownames(DAList$annotation)
  }
  
  # Collect per-contrast matrices if available
  dpc <- DAList$data_per_contrast
  has_dpc <- !is.null(dpc) && length(dpc) > 0 && any(vapply(dpc, is.data.frame, logical(1)) | vapply(dpc, is.matrix, logical(1)))
  
  if (prefer == "data_per_contrast" && has_dpc) {
    # union proteins from dpc (+ fallback to DAList$data and annotation)
    dpc_rows <- unique(unlist(lapply(dpc, rownames)))
    if (is.null(all_proteins)) {
      all_proteins <- unique(c(rownames(DAList$data), dpc_rows))
    } else {
      all_proteins <- unique(c(all_proteins, dpc_rows))
    }
    
    # initialize export matrix (data.frame to preserve colnames + types)
    export <- as.data.frame(matrix(NA_real_, nrow = length(all_proteins), ncol = length(sample_ids)),
                            stringsAsFactors = FALSE)
    rownames(export) <- all_proteins
    colnames(export) <- sample_ids
    
    warned <- FALSE
    
    for (nm in names(dpc)) {
      x <- dpc[[nm]]
      if (is.null(x)) next
      x <- as.data.frame(x, stringsAsFactors = FALSE)
      
      # align columns to canonical samples (drop extras, keep order)
      common_samples <- intersect(sample_ids, colnames(x))
      if (length(common_samples) == 0) next
      x <- x[, common_samples, drop = FALSE]
      
      # restrict to proteins we track
      common_prots <- intersect(all_proteins, rownames(x))
      if (length(common_prots) == 0) next
      
      # Fill: only overwrite NAs; detect conflicts if non-NA differs
      cur <- export[common_prots, common_samples, drop = FALSE]
      new <- x[common_prots, common_samples, drop = FALSE]
      
      # conflict check: both non-NA and different
      if (!warned) {
        both <- !is.na(cur) & !is.na(new)
        if (any(both)) {
          dif <- abs(as.matrix(cur[both]) - as.matrix(new[both])) > 1e-8
          if (any(dif, na.rm = TRUE)) {
            warning("Normalized values in data_per_contrast differ across contrasts for some protein/sample cells. Using first encountered values.",
                    call. = FALSE)
            warned <- TRUE
          }
        }
      }
      
      # overwrite only where export is NA and new is non-NA
      fill_idx <- is.na(cur) & !is.na(new)
      if (any(fill_idx)) {
        cur[fill_idx] <- new[fill_idx]
        export[common_prots, common_samples] <- cur
      }
    }
    
    return(export)
  }
  
  # Fallback: just use DAList$data (raw or normalized depending on pipeline)
  x <- as.data.frame(DAList$data, stringsAsFactors = FALSE)
  # align to canonical sample order
  x <- x[, sample_ids, drop = FALSE]
  
  # if annotation exists, expand to full protein universe
  if (!is.null(all_proteins)) {
    out <- as.data.frame(matrix(NA_real_, nrow = length(all_proteins), ncol = length(sample_ids)),
                         stringsAsFactors = FALSE)
    rownames(out) <- all_proteins
    colnames(out) <- sample_ids
    common <- intersect(all_proteins, rownames(x))
    out[common, ] <- x[common, , drop = FALSE]
    return(out)
  }
  
  x
}

# now aware of contrast names in the design matrix

#' Write per-contrast results csvs
#'
#' Internal function used to write per-contrast statistical results to .csv files.
#'
#' @param annotation_df A data frame of annotation data for each gene/protein.
#' @param data A data frame of average expression data for each sample.
#' @param results_statlist A list of per-contrast DE results.
#' @param output_dir The directory in which to save the per-contrast csv files.
#' @param annotation_cols Optional character vector of annotation column names to include.
#' @param metadata Optional data.frame of sample metadata (e.g. DAList$metadata).
#' @param group_col Character. Name of the grouping column within metadata to match contrast samples.
#' @param stat_cols Optional character vector of statistical columns to include from the results_statlist.
#' @param per_contrast_tags Optional list like DAList$tags$per_contrast; if present, will
#'   be used to determine which sample groups/levels are involved in each contrast.
#'
#' @return A logical vector indicating whether each contrast file was successfully written.
#' @keywords internal
write_per_contrast_csvs <- function(annotation_df,
                                    data,
                                    results_statlist,
                                    output_dir,
                                    annotation_cols = NULL,
                                    metadata = NULL,
                                    group_col = "group",
                                    stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSD", "logFC_z_scores"),
                                    per_contrast_tags = NULL) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Keep only requested annotation columns (if provided)
  if (!is.null(annotation_cols)) {
    keep_annot <- intersect(annotation_cols, colnames(annotation_df))
    annotation_df <- annotation_df[, keep_annot, drop = FALSE]
  }
  
  # Excel-friendly gene_symbol trick (leave as-is if you already rely on it)
  if ("gene_symbol" %in% colnames(annotation_df)) {
    annotation_df$gene_symbol <- paste0('"=""', annotation_df$gene_symbol, '"""')
  }
  csv_quote_cols <- which(colnames(annotation_df) != "gene_symbol")
  
  # Build a lookup of sample -> group
  sample_groups <- rep("all", ncol(data))
  names(sample_groups) <- colnames(data)
  if (!is.null(metadata) && !is.null(group_col) && group_col %in% colnames(metadata)) {
    sample_groups <- metadata[[group_col]]
    names(sample_groups) <- rownames(metadata)
  }
  
  # One contrast at a time
  per_contrast_results <- lapply(names(results_statlist), function(contrast) {
    # Default: include all samples
    included_samples <- colnames(data)
    
    # 1) Try tags first (preferred)
    ci <- NULL
    if (!is.null(per_contrast_tags) &&
        !is.null(per_contrast_tags[[contrast]]) &&
        !is.null(per_contrast_tags[[contrast]]$contrast_info)) {
      ci <- per_contrast_tags[[contrast]]$contrast_info
    }
    
    if (!is.null(ci) &&
        !is.null(ci$group_col) &&
        length(ci$involved_levels %||% character(0)) > 0 &&
        !is.null(metadata) &&
        ci$group_col %in% colnames(metadata)) {
      gcol <- ci$group_col
      cand <- rownames(metadata)[metadata[[gcol]] %in% ci$involved_levels]
      if (length(cand) > 0) {
        included_samples <- cand
      } else {
        warning(sprintf("Tag metadata present but no matching samples for contrast: %s", contrast))
      }
    } else {
      # 2) Fallback to legacy A_vs_B parsing
      groups_in_contrast <- unlist(strsplit(contrast, split = "_vs_"))
      if (!is.null(metadata) && !is.null(group_col) && group_col %in% colnames(metadata)) {
        cand <- rownames(metadata)[metadata[[group_col]] %in% groups_in_contrast]
        if (length(cand) > 0) {
          included_samples <- cand
        } else if (length(groups_in_contrast) > 1) {
          warning(paste("No matching sample groups found for contrast:", contrast))
        }
      }
    }
    
    subset_data <- data[rownames(data), included_samples, drop = FALSE]
    
    # Subset stat results to selected columns
    contrast_stats <- results_statlist[[contrast]]
    # if (!is.null(stat_cols)) {
    #   available <- colnames(contrast_stats)
    #   keep_cols <- intersect(stat_cols, available)
    #   missing <- setdiff(stat_cols, available)
    #   if (length(missing) > 0) {
    #     warning(paste("Missing columns in contrast", contrast, ":", paste(missing, collapse = ", ")))
    #   }
    #   contrast_stats <- contrast_stats[, keep_cols, drop = FALSE]
    # }
    if (!is.null(stat_cols)) {
      available <- colnames(contrast_stats)
      keep_cols <- intersect(stat_cols, available)
      missing   <- setdiff(stat_cols, available)
      
      if (length(missing) > 0) {
        # These are optional / tag-derived and may not always be present
        optional_cols   <- c("movingSD", "movingSDs", "logFC_z_scores")
        non_optional    <- setdiff(missing, optional_cols)
        
        # Only warn if something *other* than the optional columns is missing
        if (length(non_optional) > 0) {
          cli::cli_inform(c(
            "Some requested statistic columns were not found and will be skipped.",
            "x" = "{.val {non_optional}}",
            "i" = "Optional columns {.val {intersect(missing, optional_cols)}} were also absent."
          ))
        }
      }
      
      contrast_stats <- contrast_stats[, keep_cols, drop = FALSE]
    }
    
    
    cbind(
      annotation_df[rownames(data), , drop = FALSE],
      subset_data,
      contrast_stats
    )
  })
  
  names(per_contrast_results) <- names(results_statlist)
  
  # Write files
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  lapply(seq_along(per_contrast_results), function(i) {
    utils::write.csv(per_contrast_results[[i]],
                     file = filenames[i],
                     row.names = FALSE,
                     quote = csv_quote_cols)
  })
  
  ok <- file.exists(filenames)
  names(ok) <- names(per_contrast_results)
  ok
}



#' Combine statistical results
#'
#' Internal function used to construct a data frame of combined results.
#'
#' @param annotation A data frame of annotation data for each gene/protein.
#' @param data A data frame of normalized intensity data for each sample.
#' @param statlist A list of per-contrast DE results.
#'
#' @return A data frame of the combined results.
#' @keywords internal
create_combined_results <- function(annotation,
                                    data,
                                    statlist) {
  row_order <- rownames(data)
  results_for_combine <- lapply(names(statlist), function(x) {
    tmp <- statlist[[x]][row_order, , drop = FALSE]
    colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
    tmp
  })
  
  cbind(annotation[row_order, , drop = FALSE],
        data[row_order, , drop = FALSE],
        results_for_combine)
}

#' Make hyperlinks for an Excel column
#'
#' Adds hyperlinks to a column in a data frame to be exported to an Excel file.
#'
#' @param data The data frame in which to add hyperlinks.
#' @param url.col The name of the column to add hyperlinks to.
#' @param url The URL that will be prepended to the info in the column.
#'
#' @return The original data frame, now with hyperlinks in the desired column.
#' @keywords internal
make_excel_hyperlinks <- function(data, url.col, url) {
  ids <- data[[url.col]]
  tmp <- is.na(ids)
  url2 <- paste0(url, ids)
  ids2 <- paste0('HYPERLINK("', url2, '", "', ids, '")')
  ids2[tmp] <- NA
  data[[url.col]] <- ids2
  class(data[[url.col]]) <- "formula"
  data
}

## new helper function to define stat columns in Excel file.
#' Build a cleaned statlist for export to Excel
#'
#' Constructs a per-contrast statistics list from a \code{DAList} object with
#' full protein coverage. For each contrast, the function assembles a data
#' frame with one row per protein in \code{DAList$data} and selected columns
#' from the limma results plus optional moving standard deviations and
#' z-scores stored in \code{DAList$tags}.
#'
#' @param DAList A DAList object containing:
#'   \itemize{
#'     \item \code{$data}: matrix/data frame of log-intensities (rows = proteins).
#'     \item \code{$results}: named list of per-contrast limma result tables.
#'     \item \code{$tags$movingSDs}: optional named list of moving SD vectors
#'           (one per contrast).
#'     \item \code{$tags$logFC_z_scores}: optional named list of z-score vectors
#'           (one per contrast).
#'   }
#' @param stat_cols Character vector of column names to include for each
#'   contrast. Defaults to:
#'   \code{c("logFC", "P.Value", "adj.P.Val", "sig.PVal", "sig.FDR",
#'   "movingSDs", "logFC_z_scores")}.
#'
#' @return A named list of data frames, one per contrast, each with one row per
#'   protein in \code{DAList$data} and columns given by \code{stat_cols}.
#' @export
#'
#' @examples
#' \dontrun{
#' statlist <- build_statlist(DAList)
#' names(statlist)
#' head(statlist[[1]])
#' }
build_statlist <- function(DAList,
                           stat_cols = c("logFC", "P.Value", "adj.P.Val",
                                         "sig.PVal", "sig.FDR",
                                         "movingSDs", "logFC_z_scores")) {
  
  all_proteins <- rownames(DAList$data)
  final_statlist <- list()
  
  for (contrast in names(DAList$results)) {
    stats_df <- data.frame(row.names = all_proteins)
    contrast_res <- DAList$results[[contrast]]
    msd <- DAList$tags$movingSDs[[contrast]]
    zs  <- DAList$tags$logFC_z_scores[[contrast]]
    
    for (col in stat_cols) {
      vec <- rep(NA_real_, length(all_proteins))
      names(vec) <- all_proteins
      
      if (col %in% colnames(contrast_res)) {
        vec[rownames(contrast_res)] <- contrast_res[[col]]
        
      } else if (col == "movingSDs" && !is.null(msd)) {
        vec[names(msd)] <- msd
        
      } else if (col == "logFC_z_scores" && !is.null(zs)) {
        vec[names(zs)] <- zs
      }
      
      stats_df[[col]] <- vec
    }
    
    final_statlist[[contrast]] <- stats_df
  }
  
  final_statlist
}

# Define available normalization methods
norm.methods <- c("Log2" = "log2",
                  "Median" = "median",
                  "Quantile" = "quantile",
                  "DIANN" = "dian_quan")

#' Write an .xlsx file of limma results
#'
#' Creates a nicely formatted Excel spreadsheet of differential expression results.
#'
#' @param filename The filename of the Excel spreadsheet to be saved.
#' @param statlist A list of the per-contrast statistical results.
#' @param annotation A data frame of annotation data for each protein.
#' @param data A data frame containing the average expression data for each sample.
#' @param norm.method The normalization method used in the model.
#' @param pval_thresh The p-value threshold for significance.
#' @param lfc_thresh The log-fold change threshold for significance.
#' @param add_filter Logical. Whether to add Excel filters to columns.
#' @param color_palette A vector of colors for contrast sections (optional).
#' @param annot_cols Character vector of annotation column names to include.
#'   If \code{NULL}, all columns are used.
#'
#' @return A list with filename and workbook object (invisibly).
#' @keywords internal
write_limma_excel <- function(filename, statlist, annotation, data, norm.method,
                              pval_thresh, lfc_thresh, add_filter, 
                              color_palette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                                "#0072B2", "#D55E00", "#CC79A7", "#999999"),
                              annot_cols = NULL,
                              # new optional params (override values from globals)
                              filt_min_reps = NULL,
                              filt_min_groups = NULL) {
  
  # DEBUG: show what was passed in for filtering params
  # message("DEBUG write_limma_excel: filt_min_reps = ", paste0(filt_min_reps), 
  #         "; filt_min_groups = ", paste0(filt_min_groups))
  # 
  
  # [Full body remains unchanged -- from your original version]
  # See previous uploads for the exact function implementation
  # This placeholder is meant to conserve response length -- paste full version here as needed
  # Maybe some argument processing
  normName <- names(norm.methods)[grep(norm.method,norm.methods)]
  row_order <- rownames(data)
  
  # initialize workbook
  wb<-openxlsx::createWorkbook()
  
  # Add Overview worksheet first
  openxlsx::addWorksheet(wb, sheetName = "Overview", tabColour = "#A9A9A9")
  
  # Section 1: Dataset summary
  overview_data <- data.frame(
    Metric = c("Number of proteins", "Number of samples", "Number of contrasts"),
    Value = c(nrow(annotation), ncol(data), length(statlist))
  )
  
  # ------------------ Section 2: Filtering parameters ------------------
  # ------------------ Section 2: Filtering parameters (defensive) ------------------
  # Prefer explicit args (filt_min_reps, filt_min_groups), then look for globals.
  filt_min_reps_val <- if (!is.null(filt_min_reps)) {
    filt_min_reps
  } else if (exists("filt_min_reps", inherits = TRUE)) {
    get("filt_min_reps", inherits = TRUE)
  } else {
    NA
  }
  
  filt_min_groups_val <- if (!is.null(filt_min_groups)) {
    filt_min_groups
  } else if (exists("filt_min_groups", inherits = TRUE)) {
    get("filt_min_groups", inherits = TRUE)
  } else {
    NA
  }
  
  # Defensive: ensure scalar (length 1) and convert NULL/zero-length -> NA
  if (is.null(filt_min_reps_val) || length(filt_min_reps_val) == 0) filt_min_reps_val <- NA
  if (is.null(filt_min_groups_val) || length(filt_min_groups_val) == 0) filt_min_groups_val <- NA
  
  # If vectors supplied, take the first element (keeps behavior predictable)
  if (length(filt_min_reps_val) > 1) filt_min_reps_val <- filt_min_reps_val[[1]]
  if (length(filt_min_groups_val) > 1) filt_min_groups_val <- filt_min_groups_val[[1]]
  
  # convert NA -> "" for clean Excel cells and coerce to character
  filt_min_reps_val   <- if (is.na(filt_min_reps_val))   "" else as.character(filt_min_reps_val)
  filt_min_groups_val <- if (is.na(filt_min_groups_val)) "" else as.character(filt_min_groups_val)
  
  filtering_data <- data.frame(
    Parameter = c("Minimum replicates per group", "Minimum number of groups with detection"),
    Value = c(filt_min_reps_val, filt_min_groups_val),
    stringsAsFactors = FALSE
  )
  # -------------------------------------------------------------------------------
  
  
  # ---------------------------------------------------------------------
  
  # Section 3: Statistical thresholds used
  criteria_data <- data.frame(
    Criterion = c("Adjusted p-value threshold (FDR)", "Log2 fold-change threshold", "Significance definition"),
    Value = c(pval_thresh, lfc_thresh, paste0(
      "Up: logFC >= ", lfc_thresh, " & adj.P.Val <= ", pval_thresh, "; ",
      "Down: logFC <= -", lfc_thresh, " & adj.P.Val <= ", pval_thresh
    ))
  )
  
  # Section 4: Contrast-level summary
  overview_summary <- do.call("rbind", lapply(names(statlist), function(contrast_name) {
    tmp <- statlist[[contrast_name]][, c("logFC", "adj.P.Val"), drop = FALSE]
    sig_up <- sum(tmp$logFC >= lfc_thresh & tmp$adj.P.Val <= pval_thresh, na.rm = TRUE)
    sig_down <- sum(tmp$logFC <= -lfc_thresh & tmp$adj.P.Val <= pval_thresh, na.rm = TRUE)
    nonsig <- sum(!(tmp$logFC >= lfc_thresh & tmp$adj.P.Val <= pval_thresh |
                      tmp$logFC <= -lfc_thresh & tmp$adj.P.Val <= pval_thresh), na.rm = TRUE)
    data.frame(
      Contrast = contrast_name,
      Upregulated = sig_up,
      Downregulated = sig_down,
      Not_Significant = nonsig
    )
    
  }))
  
  overview_unadj_summary <- do.call("rbind", lapply(names(statlist), function(contrast_name) {
    tmp <- statlist[[contrast_name]][, c("logFC", "P.Value"), drop = FALSE]
    sig_up <- sum(tmp$logFC >= lfc_thresh & tmp$P.Value <= pval_thresh, na.rm = TRUE)
    sig_down <- sum(tmp$logFC <= -lfc_thresh & tmp$P.Value <= pval_thresh, na.rm = TRUE)
    nonsig <- sum(!(tmp$logFC >= lfc_thresh & tmp$P.Value <= pval_thresh |
                      tmp$logFC <= -lfc_thresh & tmp$P.Value <= pval_thresh), na.rm = TRUE)
    data.frame(
      Contrast = contrast_name,
      Upregulated_UnadjP = sig_up,
      Downregulated_UnadjP = sig_down,
      Not_Significant_UnadjP = nonsig
    )
  }))
  
  
  # Write everything to the Overview sheet
  openxlsx::writeData(wb, sheet = "Overview", x = "Dataset Summary", startRow = 1, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = overview_data, startRow = 2, startCol = 1)
  
  openxlsx::writeData(wb, sheet = "Overview", x = "Filtering Criteria", startRow = 6, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = filtering_data, startRow = 7, startCol = 1)
  
  openxlsx::writeData(wb, sheet = "Overview", x = "Significance Criteria", startRow = 10, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = criteria_data, startRow = 11, startCol = 1)
  
  contrast_start_row <- 15
  openxlsx::writeData(wb, sheet = "Overview", x = "Contrast Summary", startRow = contrast_start_row, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = overview_summary, startRow = contrast_start_row + 1, startCol = 1, withFilter = TRUE)
  
  unadj_summary_start_row <- contrast_start_row + nrow(overview_summary) + 4
  openxlsx::writeData(wb, sheet = "Overview", x = "Unadjusted P-Value Contrast Summary", startRow = unadj_summary_start_row, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = overview_unadj_summary, startRow = unadj_summary_start_row + 1, startCol = 1, withFilter = TRUE)
  
  # Styling
  header_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12, halign = "left")
  openxlsx::addStyle(wb, "Overview", header_style, rows = c(1, 6, 10, contrast_start_row), cols = 1, gridExpand = TRUE)
  
  table_header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAD3")
  openxlsx::addStyle(wb, "Overview", table_header_style, rows = c(2, 7, 11, contrast_start_row + 1), cols = 1:3, gridExpand = TRUE)
  
  openxlsx::addStyle(wb, "Overview", header_style, rows = unadj_summary_start_row, cols = 1, gridExpand = TRUE)
  openxlsx::addStyle(wb, "Overview", table_header_style, rows = unadj_summary_start_row + 1, cols = 1:4, gridExpand = TRUE)
  
  # Optional: freeze pane below contrast summary headers
  openxlsx::freezePane(wb, sheet = "Overview", firstActiveRow = contrast_start_row + 2)
  
  
  # Set up annotation columns -----------------------------------------------
  newNames <- stringr::str_replace_all(colnames(annotation), pattern = "uniprot_id", replacement = "UniProt ID") |>
    stringr::str_remove_all("[._]")
  
  annot.title <- "Protein Annotation"
  data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Intensities")
  sheetName <- "Protein Results"
  
  # Initialize worksheet ----------------------------------------------------
  openxlsx::addWorksheet(wb, sheetName = sheetName, tabColour = "#FFE100")
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 3, heights = 30) # increase row height for the title row
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 4, heights = 20) # increase row height for the header row
  title_row <- 3
  
  # Add annotation columns -----------------------------------------------------
  # Subset and rename for hyperlinks 
  if (!is.null(annot_cols)) {
    missing_cols <- setdiff(annot_cols, colnames(annotation))
    if (length(missing_cols) > 0) {
      cli::cli_inform(c(
        "Some requested annotation columns were not found and will be skipped.",
        "x" = "{.val {missing_cols}}"
      ))
    }
    keep_cols <- intersect(annot_cols, colnames(annotation))
    annot <- annotation[row_order, keep_cols, drop = FALSE]
  } else {
    annot <- annotation[row_order, , drop = FALSE]
  }
  
  
  # Rename columns
  newNames <- stringr::str_replace_all(colnames(annot), pattern = "uniprot_id", replacement = "UniProt ID") |>
    stringr::str_remove_all("[._]")
  colnames(annot) <- newNames
  
  # Add UniProt hyperlinks only if column is present
  if ("UniProt ID" %in% colnames(annot)) {
    annot <- make_excel_hyperlinks(data = annot,
                                   url.col = "UniProt ID",
                                   url = "https://www.uniprot.org/uniprot/")
  }
  
  
  ## merge annotation columns
  annot_col_start <- 1
  annot_col_end <- ncol(annot)
  
  # Make title row
  openxlsx::mergeCells(wb, sheet = sheetName,
                       cols = annot_col_start:annot_col_end,
                       rows = title_row)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = annot.title,
                      startCol = annot_col_start, startRow = title_row,
                      colNames = FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill = "#636262",
                                      halign = "center", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "white", fontSize = 14,
                                      numFmt = "TEXT")
  openxlsx::addStyle(wb = wb, sheet = sheetName,
                     style = mergeStyle, rows = 3, cols = 1, stack = TRUE)
  
  ## write annotation data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#d8d8d8",
                                    halign = "left", valign = "center",
                                    textDecoration = "bold",
                                    fontColour = "black", fontSize=11,
                                    numFmt="TEXT")
  openxlsx::addStyle(wb, sheet = sheetName, style = colStyle,
                     rows = title_row + 1,
                     cols = annot_col_start:annot_col_end, stack = TRUE)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = annot,
                      startCol = annot_col_start,
                      startRow = title_row + 1,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders = "columns",
                      keepNA = TRUE,
                      na.string = "NA",
                      sep = "\t")
  
  ## add hyperlink style (blue font, underlined)
  hyperlinkStyle <- openxlsx::createStyle(fontColour = "#0000FF",
                                          halign = "left",
                                          valign = "center",
                                          textDecoration = "underline")
  cols <- grep("UniProt ID", colnames(annot))
  openxlsx::addStyle(wb, sheet = sheetName,
                     style = hyperlinkStyle, cols = cols,
                     rows = (title_row + 2):(nrow(annot) + 4),
                     gridExpand = TRUE, stack = TRUE)
  
  
  
  
  # Add data columns --------------------------------------------------------
  ## start/stop column positions for data
  data_col_start <- annot_col_end + 1
  data_col_end <- annot_col_end + ncol(data)
  
  ## Set up title row
  openxlsx::mergeCells(wb, sheet = sheetName,
                       cols = data_col_start:data_col_end,
                       rows = title_row)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = data.title,
                      startCol = data_col_start,
                      startRow = title_row, colNames = FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill = "#516285",
                                      halign = "center", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "white",
                                      fontSize = 14, numFmt = "TEXT")
  openxlsx::addStyle(wb, sheet = sheetName,
                     style = mergeStyle,
                     rows = title_row, cols = data_col_start, stack = TRUE)
  
  ## write data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#cdd4e1",
                                    halign = "left", valign = "center",
                                    textDecoration = "bold",
                                    fontColour = "black", fontSize = 11,
                                    numFmt = "TEXT")
  openxlsx::addStyle(wb, sheet = sheetName, style = colStyle,
                     rows = title_row + 1, cols = data_col_start:data_col_end, stack=TRUE)
  
  openxlsx::writeData(wb, sheet = sheetName, x = data[row_order, ,drop = F],
                      startCol = data_col_start,
                      startRow = title_row + 1,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders = "columns",
                      keepNA = TRUE, na.string = "NA",
                      sep="\t")
  
  
  # Add stats columns -------------------------------------------------------
  # colors for stat columns
  if (length(names(statlist)) > length(color_palette)) {
    numColors<- ceiling(length(names(statlist))/length(color_palette)) + 1
    colors2 <- rep(color_palette, numColors)
    lightcolors2 <- rep(color_palette, numColors)
  } else {
    if (length(names(statlist)) <= length(color_palette)) {
      colors2 <- color_palette
      lightcolors2 <- color_palette
    }
  }
  stopifnot(length(colors2) >= length(names(statlist)))
  
  
  ## write stat data to worksheet
  stat.start<- NULL
  stat.end  <- NULL
  
  for(i in base::seq_along(names(statlist))) {
    
    stats <- statlist[[i]]
    stats<-stats[row_order, , drop = F]
    comparison <- names(statlist)[i]
    title_color <- colors2[i]
    header_color <- lightcolors2[i]
    
    ## stats start
    if(is.null(stat.start)){
      stat.start = data_col_end + 1
    } else {
      if(!is.null(stat.start)){
        stat.start = stat.end + 1}
    }
    
    ## stats end
    if (is.null(stat.end)) {
      stat.end = data_col_end + ncol(stats)
    } else {
      if (!is.null(stat.end)) {
        stat.end = stat.end + ncol(stats)
      }
    }
    
    # set up title row
    openxlsx::mergeCells(wb, sheet = sheetName,
                         cols = stat.start:stat.end, rows = title_row)
    openxlsx::writeData(wb, sheet = sheetName,
                        x = comparison,
                        startCol = stat.start, startRow = title_row, colNames = FALSE)
    mergeStyle <- openxlsx::createStyle(fgFill = title_color, halign = "center",
                                        valign = "center", textDecoration = "bold",
                                        fontColour = "white", fontSize = 14, numFmt = "TEXT")
    openxlsx::addStyle(wb, sheet = sheetName, style = mergeStyle, rows = title_row, cols = stat.start, stack=TRUE)
    
    
    ## write stats to worksheet
    colStyle <- openxlsx::createStyle(fgFill = header_color,
                                      halign = "left", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "black", fontSize = 11, numFmt = "TEXT")
    openxlsx::addStyle(wb, sheet = sheetName, style = colStyle, rows = title_row + 1,
                       cols = stat.start:stat.end, stack = TRUE)
    openxlsx::writeData(wb, sheet = sheetName, x = stats,
                        startCol = stat.start,
                        startRow = title_row + 1,
                        colNames = TRUE,
                        rowNames = FALSE,
                        borders = "columns",
                        keepNA = TRUE, na.string = "NA",
                        sep = "\t")
    
    # Add conditional formatting to stats section
    posStyle  <- openxlsx::createStyle(fontColour="#990000", bgFill="#FFC7CE")
    negStyle  <- openxlsx::createStyle(fontColour="#006100", bgFill="#C6EFCE")
    normStyle <- openxlsx::createStyle(fontColour="#000000", bgFill="#FFFFFF")
    
    # fc.col    <- grep("logFC", colnames(stats)) + stat.start - 1
    # fc.rule1  <- paste0(">=", lfc_thresh)
    # fc.rule2  <- paste0("<=", -lfc_thresh)
    # fdr.col   <- grep("adj.P.Val", colnames(stats)) + stat.start - 1
    # fdr.rule  <- paste0("<=", pval_thresh)
    # pval.col  <- grep("P.Value", colnames(stats)) + stat.start - 1
    # pval.rule <- paste0("<=", pval_thresh)
    
    fc.idx   <- grep("^logFC$", colnames(stats))
    fdr.idx  <- grep("^adj.P.Val$", colnames(stats))
    pval.idx <- grep("^P.Value$", colnames(stats))
    
    fc.col    <- fc.idx + stat.start - 1
    fdr.col   <- fdr.idx + stat.start - 1
    pval.col  <- pval.idx + stat.start - 1
    
    fc.rule1  <- paste0(">=", lfc_thresh)
    fc.rule2  <- paste0("<=", -lfc_thresh)
    fdr.rule  <- paste0("<=", pval_thresh)
    pval.rule <- paste0("<=", pval_thresh)
    
    
    
    openxlsx::conditionalFormatting(wb, sheet = sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) + 4),
                                    type = "expression", rule = fc.rule1, style = posStyle)
    openxlsx::conditionalFormatting(wb = wb, sheet = sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) + 4),
                                    type="expression", rule = fc.rule2, style = negStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) +4),
                                    type="contains", rule="NA", style = normStyle)
    
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=fdr.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)
    
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=pval.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)
    
    
  }
  
  ## add global styles
  headerStyle <- openxlsx::createStyle(halign = "center", valign="center", #fontColour = "#000000",
                                       border = "TopBottomLeftRight", borderColour = c("black","black","black","black"),
                                       textDecoration = "bold")
  openxlsx::addStyle(wb, sheet = sheetName, style = headerStyle,
                     cols = 1:stat.end, rows = title_row,
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, sheet = sheetName, style = headerStyle,
                     cols = 1:stat.end, rows = title_row + 1,
                     gridExpand = TRUE, stack = TRUE)
  if (add_filter) {
    openxlsx::addFilter(wb, sheet = sheetName, cols = 1:stat.end, rows = title_row + 1)
  }
  
  
  # Save workbook
  openxlsx::saveWorkbook(wb = wb,
                         file = filename,
                         overwrite = T)
  
  invisible(list(file = filename,
                 wb = wb))
  }

# 
#' Sanitize Excel worksheet names
#'
#' Ensures worksheet names are compatible with Excel and \code{openxlsx}.
#' Names are:
#' \itemize{
#'   \item Truncated to \code{max_len} characters (default 31, Excel's limit).
#'   \item Cleaned of invalid characters (\code{: \\ / ? * [ ]}).
#'   \item Made unique by appending numeric suffixes (e.g. \code{"name"}, \code{"name_1"}).
#' }
#' If any names are modified, a warning is emitted listing the original and
#' sanitized names.
#'
#' @param names Character vector of proposed worksheet names.
#' @param max_len Integer maximum length for each name (default 31).
#'
#' @return A character vector of sanitized, unique worksheet names.
#' @keywords internal
sanitize_sheet_names <- function(names, max_len = 31L) {
  if (!is.character(names)) {
    names <- as.character(names)
  }
  if (length(names) == 0L) {
    return(names)
  }
  if (max_len < 1L) {
    stop("max_len must be at least 1.", call. = FALSE)
  }
  
  orig <- names
  
  # Replace invalid Excel characters: : \ / ? * [ ]
  names <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", names)
  
  # Replace empty or NA names with a default base
  names[is.na(names) | names == ""] <- "Sheet"
  
  # Truncate to max_len
  names <- substr(names, 1L, max_len)
  
  # Ensure uniqueness by adding suffixes as needed
  seen <- character()
  out  <- character(length(names))
  
  for (i in seq_along(names)) {
    base_name <- names[i]
    candidate <- base_name
    
    # If already seen, start appending suffixes
    if (candidate %in% seen) {
      j <- 1L
      repeat {
        suffix <- paste0("_", j)
        # Leave room for suffix within max_len
        base_trunc <- substr(base_name, 1L, max_len - nchar(suffix))
        candidate  <- paste0(base_trunc, suffix)
        if (!(candidate %in% seen)) break
        j <- j + 1L
      }
    }
    
    out[i] <- candidate
    seen   <- c(seen, candidate)
  }
  
  changed <- orig != out
  if (any(changed)) {
    warning(
      "Some worksheet names were adjusted to be Excel-compatible (length/characters/uniqueness):\n",
      paste0("  '", orig[changed], "' -> '", out[changed], "'", collapse = "\n"),
      call. = FALSE
    )
  }
  
  out
}

#' Write tables of limma results
#'
#' Generates CSV and Excel summary tables from a limma differential
#' abundance analysis stored in a \code{DAList}.
#'
#' @details
#' For each contrast, per-contrast CSV files are written using
#' \code{write_per_contrast_csvs()}, combining annotation, data, and
#' statistics. If \code{DAList$tags$per_contrast[[label]]$contrast_info}
#' contains a usable \code{group_col} and non-empty \code{involved_levels},
#' those samples are selected for the per-contrast CSV. Otherwise the older
#' \code{"_vs_"} label parser is used. If that also fails, all samples are
#'  included with a warning.
#'
#' The function then:
#' \itemize{
#'   \item writes a summary CSV of up/down/non-significant counts per contrast;
#'   \item writes per-contrast CSV files into \code{contrasts_subdir};
#'   \item writes a combined results CSV with all contrasts side by side;
#'   \item writes an Excel workbook using \code{write_limma_excel()}, and
#'         optionally adds each per-contrast CSV as a worksheet.
#' }
#'
#' The columns included in the per-contrast and Excel tables are determined
#' by internal settings such as \code{DA_table_cols} and \code{stat_cols},
#' which are defined in the package configuration rather than as function
#' arguments.
#'
#' @param DAList A DAList object with statistical results in \code{$results},
#'   annotation in \code{$annotation}, data in \code{$data}, metadata in
#'   \code{$metadata}, and optional tags in \code{$tags}.
#' @param output_dir Directory to output tables. Defaults to the current
#'   working directory if \code{NULL}.
#' @param overwrite Logical; overwrite existing files if they are found in
#'   \code{output_dir}? Default is \code{FALSE}.
#' @param contrasts_subdir Subdirectory (within \code{output_dir}) for
#'   per-contrast CSV files.
#' @param summary_csv Filename for the summary CSV (per-contrast counts).
#' @param combined_file_csv Filename for the combined results CSV.
#' @param spreadsheet_xlsx Filename for the Excel spreadsheet.
#' @param add_filter Logical; add Excel autofilters to columns in the main
#'   worksheet produced by \code{write_limma_excel()}.
#' @param color_palette Optional vector of colors used for conditional
#'   formatting in the Excel output.
#' @param add_contrast_sheets Logical; if \code{TRUE}, each per-contrast CSV
#'   is added as a separate worksheet to the Excel workbook.
#' @param statlist Optional list of per-contrast result tables to use instead
#'   of \code{DAList$results}. If \code{NULL}, \code{DAList$results} is used.
#'
#' @return Invisibly returns the (validated) input \code{DAList}.
#' @export
write_limma_tables <- function(DAList,
                               output_dir = NULL,
                               overwrite = FALSE,
                               contrasts_subdir = NULL,
                               summary_csv = NULL,
                               combined_file_csv = NULL,
                               spreadsheet_xlsx = NULL,
                               add_filter = TRUE,
                               color_palette = NULL,
                               add_contrast_sheets = TRUE,
                               statlist = NULL) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
  
  if (is.null(color_palette)) {
    color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  }
  
  # validate DAList unless statlist overrides
  input_DAList <- if (is.null(statlist)) validate_DAList(DAList) else DAList
  
  if (is.null(DAList$results) && is.null(statlist)) {
    cli::cli_abort(c("Input DAList does not have results",
                     "i" = "Run {.code DAList <- extract_DA_results(DAList, ~ formula)}"))
  }
  
  # Choose best intensities for export (prefer normalized per-contrast)
  export_data <- get_export_data(input_DAList, prefer = "data_per_contrast")
  
  # results to use
  results_statlist <- statlist %||% DAList$results
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform("{.arg output_dir} argument is empty.")
    cli::cli_inform("Setting output directory to current working directory:")
    cli::cli_inform("{.path {output_dir}}")
  }
  
  if (is.null(contrasts_subdir)) contrasts_subdir <- "per_contrast_results"
  if (is.null(summary_csv)) summary_csv <- "DA_summary.csv"
  if (is.null(combined_file_csv)) combined_file_csv <- "combined_results.csv"
  if (is.null(spreadsheet_xlsx)) spreadsheet_xlsx <- "results.xlsx"
  
  for (filename in c(summary_csv, combined_file_csv)) validate_filename(filename, allowed_exts = "csv")
  validate_filename(spreadsheet_xlsx, allowed_exts = "xlsx")
  
  expected_per_contrast_tables <- file.path(output_dir, contrasts_subdir,
                                            paste0(names(results_statlist), ".csv"))
  summary_output_file  <- file.path(output_dir, summary_csv)
  combined_output_file <- file.path(output_dir, combined_file_csv)
  excel_output_file    <- file.path(output_dir, spreadsheet_xlsx)
  
  expected_files <- c(expected_per_contrast_tables, summary_output_file,
                      combined_output_file, excel_output_file)
  
  if (any(file.exists(expected_files))) {
    if (!overwrite) {
      cli::cli_abort(c("Results files already exist",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Change {.arg output_dir} or set {.arg overwrite} to {.val TRUE}"))
    } else {
      cli::cli_inform("Results files already exist, and {.arg overwrite} == {.val {overwrite}}. Overwriting results files.")
      unlink(expected_files)
    }
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Summary
  cli::cli_inform("Writing DA summary table to {.path {summary_output_file}}")
  summary <- do.call("rbind", lapply(names(results_statlist), summarize_contrast_DA, results_statlist))
  if (!is.null(DAList$tags$DA_criteria)) {
    summary$pval_thresh  <- DAList$tags$DA_criteria$pval_thresh
    summary$lfc_thresh   <- DAList$tags$DA_criteria$lfc_thresh
    summary$p_adj_method <- DAList$tags$DA_criteria$adj_method
  }
  utils::write.csv(summary, file = summary_output_file, row.names = FALSE)
  if (!file.exists(summary_output_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_output_file}}"))
  }
  
  # Per-contrast CSVs
 # per_contrast_dir <- file.path(output_dir, contrasts_subdir)
 # if (!dir.exists(per_contrast_dir)) dir.create(per_contrast_dir, recursive = TRUE)
  
  per_contrast_dir <- file.path(output_dir, contrasts_subdir)
  
  if (overwrite && dir.exists(per_contrast_dir)) {
    unlink(per_contrast_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(per_contrast_dir, recursive = TRUE, showWarnings = FALSE)
  
  cli::cli_inform("Writing per-contrast results {.path .csv} files to {.path {per_contrast_dir}}")
  
  contrast_csv_success <- write_per_contrast_csvs(
    annotation_df = DAList$annotation,
    data = export_data,                     #DAList$data,
    results_statlist = results_statlist,
    output_dir = per_contrast_dir,
    annotation_cols = DA_table_cols,      # from your params file
    metadata = DAList$metadata,
    group_col = "group",
    stat_cols = stat_cols,                # from your params file
    per_contrast_tags = DAList$tags$per_contrast %||% NULL
  )
  
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success]
    cli::cli_abort(c("Failed to write {.path .csv} results for contrast(s):",
                     "!" = "{.val {failed}}"))
  }
  
  # Combined CSV
  cli::cli_inform("Writing combined results table to {.path {combined_output_file}}")
  combined_results <- create_combined_results(DAList$annotation,
                                              export_data,  # DAList$data, 
                                              results_statlist)
  if ("gene_symbol" %in% colnames(combined_results)) {
    combined_results$gene_symbol <- paste0('"=""', combined_results$gene_symbol, '"""')
  }
  csv_quote_cols <- which(colnames(combined_results) != "gene_symbol")
  utils::write.csv(combined_results, file = combined_output_file, row.names = FALSE, quote = csv_quote_cols)
  
  # Excel workbook
  cli::cli_inform("Writing combined results Excel spreadsheet to {.path {excel_output_file}}")
  write_limma_excel(
    filename = excel_output_file,
    statlist = results_statlist,
    annotation = DAList$annotation,
    data = export_data,                      #DAList$data,
    norm.method = DAList$tags$norm_method,
    pval_thresh = DAList$tags$DA_criteria$pval_thresh,
    lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
    add_filter = add_filter,
    color_palette = color_palette,
    annot_cols = DA_table_cols,
    # Pass the filtering params saved on the DAList tags (may be NULL)
    filt_min_reps  = if (!is.null(input_DAList$tags$filt_min_reps))  input_DAList$tags$filt_min_reps  else NULL,
    filt_min_groups= if (!is.null(input_DAList$tags$filt_min_groups)) input_DAList$tags$filt_min_groups else NULL
  )
  
  if (add_contrast_sheets) {
    cli::cli_inform("Adding per-contrast CSVs as worksheets to {.path {excel_output_file}}")
    wb <- openxlsx::loadWorkbook(excel_output_file)
    per_contrast_files <- list.files(per_contrast_dir, pattern = "\\.csv$", full.names = TRUE)
    
    # --- NEW: sanitize worksheet names for Excel before adding ---
    raw_sheet_names <- tools::file_path_sans_ext(basename(per_contrast_files))
    sheet_names     <- sanitize_sheet_names(raw_sheet_names)
    
    for (i in seq_along(per_contrast_files)) {
      csv_file   <- per_contrast_files[i]
      sheet_name <- sheet_names[i]
      df <- utils::read.csv(csv_file, check.names = FALSE)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = df, withFilter = TRUE)
    }
    openxlsx::saveWorkbook(wb, file = excel_output_file, overwrite = TRUE)
  }
  
  invisible(input_DAList)
}
