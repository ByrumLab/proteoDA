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

#' Write per-contrast results csvs
#'
#' Internal function used to write per-contrast statistical results to .csv files.
#'
#' @param annotation_df A data frame of annotation data for each gene/protein.
#' @param data A data frame of average expression data for each sample.
#' @param results_statlist A list of per-contrast DE results.
#' @param output_dir The directory in which to save the per-contrast csv files.
#' @param annotation_cols Optional character vector of annotation column names to include.
#' @param metadata Optional list of metadata (e.g. DAList$metadata) for filtering data columns by sample group.
#' @param group_col Character. Name of the grouping column within metadata to match contrast samples.
#' @param stat_cols Optional character vector of statistical columns to include from the results_statlist. If NULL, include all.
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
                                    stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores")) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (!is.null(annotation_cols)) {
    annotation_df <- annotation_df[, annotation_cols, drop = FALSE]
  }
  
  if ("gene_symbol" %in% colnames(annotation_df)) {
    annotation_df$gene_symbol <- paste0('"=""', annotation_df$gene_symbol, '"""')
  }
  
  csv_quote_cols <- which(colnames(annotation_df) != "gene_symbol")
  
  sample_groups <- rep("all", ncol(data))
  names(sample_groups) <- colnames(data)
  
  if (!is.null(metadata) && !is.null(metadata[[group_col]])) {
    sample_groups <- metadata[[group_col]]
    names(sample_groups) <- rownames(metadata)
  }
  
  per_contrast_results <- lapply(names(results_statlist), function(contrast) {
    groups_in_contrast <- unlist(strsplit(contrast, split = "_vs_"))
    
    included_samples <- names(sample_groups)[sample_groups %in% groups_in_contrast]
    
    subset_data <- if (length(included_samples) > 0) {
      data[rownames(data), included_samples, drop = FALSE]
    } else {
      warning(paste("No matching sample groups found for contrast:", contrast))
      data[rownames(data), , drop = FALSE]
    }
    
    # Subset stat results to selected columns
    contrast_stats <- results_statlist[[contrast]]
    if (!is.null(stat_cols)) {
      available <- colnames(contrast_stats)
      keep_cols <- intersect(stat_cols, available)
      missing <- setdiff(stat_cols, available)
      if (length(missing) > 0) {
        warning(paste("Missing columns in contrast", contrast, ":", paste(missing, collapse = ", ")))
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
  
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  lapply(seq_along(per_contrast_results), function(x) {
    utils::write.csv(per_contrast_results[[x]], file = filenames[x], row.names = FALSE, quote = csv_quote_cols)
  })
  
  write_success <- file.exists(filenames)
  names(write_success) <- names(per_contrast_results)
  write_success
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
#' Build a cleaned statlist for export to Excel (renaming z_scores)
#'
#' @param results A named list of data frames. Each element is a contrast result (e.g., DAList$results$contrast_name).
#' @param movingSDs A named list of numeric vectors for each contrast (e.g., DAList$tags$movingSDs).
#' @param z_scores A named list of numeric vectors for each contrast (e.g., DAList$tags$logFC_z_scores).
#' @param stat_cols Character vector of column names to include from the results. Defaults include additional stats from tags.
#'
#' @return A named list of filtered data frames per contrast.
#' @export
build_statlist <- function(results, movingSDs, z_scores,
                           stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores")) {
  contrasts <- names(results)
  statlist <- list()
  
  for (contrast in contrasts) {
    res <- results[[contrast]]
    msd <- movingSDs[[contrast]]
    z   <- z_scores[[contrast]]
    
    base_cols <- intersect(stat_cols, colnames(res))
    df <- res[, base_cols, drop = FALSE]
    
    # Add tag columns only if requested
    if ("movingSDs" %in% stat_cols) {
      df$movingSDs <- msd
    }
    if ("logFC_z_scores" %in% stat_cols) {
      df$logFC_z_scores <- z
    }
    
    statlist[[contrast]] <- df
  }
  
  return(statlist)
}






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
#' @param annot_cols Character vector of annotation column names to include. If NULL, all columns are used.

#'
#' @return A list with filename and workbook object (invisibly).
#' @keywords internal
#' 

# Define available normalization methods
norm.methods <- c("Log2" = "log2", "Median" = "median", "Quantile" = "quantile", "DIANN" = "dian_quan")

write_limma_excel <- function(filename, statlist, annotation, data, norm.method,
                              pval_thresh, lfc_thresh, add_filter, color_palette,
                              annot_cols = NULL) {
  # [Full body remains unchanged — from your original version]
  # See previous uploads for the exact function implementation
  # This placeholder is meant to conserve response length — paste full version here as needed
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
  
  # Section 2: Statistical thresholds used
  criteria_data <- data.frame(
    Criterion = c("Adjusted p-value threshold (FDR)", "Log2 fold-change threshold", "Significance definition"),
    Value = c(pval_thresh, lfc_thresh, paste0(
      "Up: logFC ≥ ", lfc_thresh, " & adj.P.Val ≤ ", pval_thresh, "; ",
      "Down: logFC ≤ -", lfc_thresh, " & adj.P.Val ≤ ", pval_thresh
    ))
  )
  
  # Section 3: Contrast-level summary
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
  
  # Write everything to the Overview sheet
  openxlsx::writeData(wb, sheet = "Overview", x = "Dataset Summary", startRow = 1, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = overview_data, startRow = 2, startCol = 1)
  
  openxlsx::writeData(wb, sheet = "Overview", x = "Significance Criteria", startRow = 6, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = criteria_data, startRow = 7, startCol = 1)
  
  openxlsx::writeData(wb, sheet = "Overview", x = "Contrast Summary", startRow = 11, startCol = 1)
  openxlsx::writeData(wb, sheet = "Overview", x = overview_summary, startRow = 12, startCol = 1, withFilter = TRUE)
  
  # Styling
  header_style <- openxlsx::createStyle(textDecoration = "bold", fontSize = 12, halign = "left")
  openxlsx::addStyle(wb, "Overview", header_style, rows = c(1, 6, 11), cols = 1, gridExpand = TRUE)
  
  table_header_style <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#D9EAD3")
  openxlsx::addStyle(wb, "Overview", table_header_style, rows = c(2, 7, 12), cols = 1:3, gridExpand = TRUE)
  
  # Optional: freeze pane below criteria
  openxlsx::freezePane(wb, sheet = "Overview", firstActiveRow = 13)
  
  
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
  
  # Subset and rename
  #annot <- annotation[row_order, , drop = F]
    # colnames(annot) <- newNames
  # annot <- make_excel_hyperlinks(data = annot,
  #                                url.col = "UniProt ID",
  #                                url = "https://www.uniprot.org/uniprot/")
  # 
  
  # Subset and rename for hyperlinks 
  if (!is.null(annot_cols)) {
    missing_cols <- setdiff(annot_cols, colnames(annotation))
    if (length(missing_cols) > 0) {
      stop(paste("The following annotation columns were not found:", paste(missing_cols, collapse = ", ")))
    }
    annot <- annotation[row_order, annot_cols, drop = FALSE]
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
    
    fc.col    <- grep("logFC", colnames(stats)) + stat.start - 1
    fc.rule1  <- paste0(">=", lfc_thresh)
    fc.rule2  <- paste0("<=", -lfc_thresh)
    fdr.col   <- grep("adj.P.Val", colnames(stats)) + stat.start - 1
    fdr.rule  <- paste0("<=", pval_thresh)
    pval.col  <- grep("P.Value", colnames(stats)) + stat.start - 1
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

#' Write tables of limma results
#'
#' Generates CSV and Excel summary tables from a limma differential analysis list.
#'
#' @param DAList A DAList object with statistical results.
#' @param output_dir Directory to output tables. Defaults to working directory.
#' @param overwrite Logical. Overwrite existing files?
#' @param contrasts_subdir Subdirectory for per-contrast CSV files.
#' @param summary_csv Filename for summary CSV.
#' @param combined_file_csv Filename for combined results CSV.
#' @param spreadsheet_xlsx Filename for Excel spreadsheet.
#' @param add_filter Logical. Add filters to Excel columns?
#' @param color_palette Optional color palette for Excel output.
#' @param add_contrast_sheets Logical. Whether to add each per-contrast CSV as a worksheet in the Excel file.
#'
#' @return Invisibly returns the input DAList.
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
                               add_contrast_sheets = TRUE) {
  
  if (is.null(color_palette)) {
    color_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
  }
  
  input_DAList <- validate_DAList(DAList)
  
  if (is.null(DAList$results)) {
    cli::cli_abort(c("Input DAList does not have results",
                     "i" = "Run {.code DAList <- extract_DA_results(DAList, ~ formula)}"))
  }
  
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
                                            paste0(names(DAList$results), ".csv"))
  summary_output_file <- file.path(output_dir, summary_csv)
  combined_output_file <- file.path(output_dir, combined_file_csv)
  excel_output_file <- file.path(output_dir, spreadsheet_xlsx)
  
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
  
  # Write summary
  cli::cli_inform("Writing DA summary table to {.path {summary_output_file}}")
  summary <- do.call("rbind", lapply(names(DAList$results), summarize_contrast_DA, DAList$results))
  summary$pval_thresh <- DAList$tags$DA_criteria$pval_thresh
  summary$lfc_thresh <- DAList$tags$DA_criteria$lfc_thresh
  summary$p_adj_method <- DAList$tags$DA_criteria$adj_method
  utils::write.csv(summary, file = summary_output_file, row.names = FALSE)
  
  if (!file.exists(summary_output_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_output_file}}"))
  }
  
  # Write per-contrast CSVs
  per_contrast_dir <- file.path(output_dir, contrasts_subdir)
  cli::cli_inform("Writing per-contrast results {.path .csv} files to {.path {per_contrast_dir}}")
  
  contrast_csv_success <- write_per_contrast_csvs(annotation_df =DAList$annotation,
                                                  data = DAList$data,
                                                  results_statlist = statlist,
                                                  output_dir = per_contrast_dir,
                                                  annotation_cols = c("uniprot_id","Genes","Accession.Number","Identified.Peptides","Protein.Description"),
                                                  metadata = DAList$metadata,
                                                  group_col = "group",
                                                  stat_cols = c("logFC", "P.Value", "adj.P.Val", "movingSDs", "logFC_z_scores")) 
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success]
    cli::cli_abort(c("Failed to write {.path .csv} results for contrast(s):",
                     "!" = "{.val {failed}}"))
  }
  
  # Write combined CSV
  cli::cli_inform("Writing combined results table to {.path {combined_output_file}}")
  combined_results <- create_combined_results(DAList$annotation, DAList$data, DAList$results)
  if ("gene_symbol" %in% colnames(combined_results)) {
    combined_results$gene_symbol <- paste0('"=""', combined_results$gene_symbol, '"""')
  }
  csv_quote_cols <- which(colnames(combined_results) != "gene_symbol")
  utils::write.csv(combined_results, file = combined_output_file, row.names = FALSE, quote = csv_quote_cols)
  
  # Write Excel spreadsheet
  cli::cli_inform("Writing combined results Excel spreadsheet to {.path {excel_output_file}}")
  write_limma_excel(filename = excel_output_file,
                    statlist = statlist,     # DAList$results
                    annotation = DAList$annotation,
                    data = DAList$data,
                    norm.method = DAList$tags$norm_method,
                    pval_thresh = DAList$tags$DA_criteria$pval_thresh,
                    lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
                    add_filter = add_filter,
                    color_palette = color_palette,
                    annot_cols = c("uniprot_id", "Genes","Identified.Peptides" ,"Accession.Number", "Protein.Description"))
  
  if (add_contrast_sheets) {
    cli::cli_inform("Adding per-contrast CSVs as worksheets to {.path {excel_output_file}}")
    wb <- openxlsx::loadWorkbook(excel_output_file)
    per_contrast_files <- list.files(per_contrast_dir, pattern = "\\.csv$", full.names = TRUE)
    for (csv_file in per_contrast_files) {
      sheet_name <- tools::file_path_sans_ext(basename(csv_file))
      df <- utils::read.csv(csv_file, check.names = FALSE)
      openxlsx::addWorksheet(wb, sheet_name)
      openxlsx::writeData(wb, sheet = sheet_name, x = df, withFilter = TRUE)
    }
    
    openxlsx::saveWorkbook(wb, file = excel_output_file, overwrite = TRUE)
  }
  
  invisible(input_DAList)
}

