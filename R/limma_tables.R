#' Write tables of limma results
#'
#' This function prepares and then saves a number of tables of results:
#' \enumerate{
#'   \item A .csv summarizing the number of direction of differentially
#'     expressed genes for each contrast (see \code{summary_csv}).
#'   \item An Excel spreadsheet of the intensity and
#'     statistical results for all contrasts/terms in the results slot of the
#'     DAList (see \code{spreadsheet_xlsx}).
#'   \item A results .csv file with the same data as the Excel spreadsheet,
#'     but in a flat .csv file withour formatting (see \code{combined_file_csv}).
#'   \item A folder of results .csv file which contain the results for each
#'     contrast.
#' }
#'
#' @param DAList A DAList object, with statistical results in the results slot.
#' @param output_dir The directory in which to output tables. If not specified,
#'   will default to the current working directory.
#' @param overwrite Should results files be overwritten? Default is FALSE.
#' @param contrasts_subdir The subdirectory within output_dir to write the per-contrast
#'   result .csv files. If not specified, will be "per_contrast_results".
#' @param summary_csv The filename of the csv file giving a summary of the number
#'   of DA genes per contrast. If not specified, will be "DA_summary.csv".
#' @param combined_file_csv The filename of the combined results csv file. If
#'   not specified, will be combined_results.csv.
#' @param spreadsheet_xlsx The filename of the Excel spreadsheet containing the
#'   results. If not specified, will be results.xlsx.
#' @param add_filter Should per-column filters be added to the spreadsheet?
#'   Default is TRUE. These can sometimes cause unstable spreadsheet files,
#'   try setting to FALSE if you're having issues with the Excel output.
#'
#'
#' @return Invisibly returns the input DAList.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Using defaults
#'   write_limma_tables(DAList)
#'
#'  # Customize output directory
#'  # and filenames
#'  write_limma_tables(DAList,
#'                     output_dir = "DA_results",
#'                     contrasts_subdir = "by_contrast",
#'                     summary_csv = "my_summary.csv",
#'                     combined_file_csv = "all_results.csv",
#'                     spreadsheet_xlsx = "all_results.xlsx")
#'
#'  # If the xlsx files have issues, try removing the
#'  # per-column filters in the spreadsheet
#'  write_limma_tables(DAList,
#'                     add_filter = F)
#' }
write_limma_tables <- function(DAList,
                               output_dir = NULL,
                               overwrite = F,
                               contrasts_subdir = NULL,
                               summary_csv = NULL,
                               combined_file_csv = NULL,
                               spreadsheet_xlsx = NULL,
                               add_filter = T) {

  # Check input arguments generally
  input_DAList <- validate_DAList(DAList)

  # Make sure there's results present,
  # tell user to set it first if not
  if (is.null(DAList$results)) {
    cli::cli_abort(c("Input DAList does not have results",
                     "i" = "Run {.code DAList <- extract_DA_results(DAList, ~ formula)}"))
  }


  # Assign defaults if not overridden above
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform("{.arg output_dir} argument is empty.")
    cli::cli_inform("Setting output directory to current working directory:")
    cli::cli_inform("{.path {output_dir}}")
  }

  if (is.null(contrasts_subdir)) {
    contrasts_subdir <- "per_contrast_results"
  }
  if (is.null(summary_csv)) {
    summary_csv <- "DA_summary.csv"
  }

  if (is.null(combined_file_csv)) {
    combined_file_csv <- "combined_results.csv"
  }

  if (is.null(spreadsheet_xlsx)) {
    spreadsheet_xlsx <- "results.xlsx"
  }

  # Validate filenames
  for (filename in c(summary_csv, combined_file_csv)) {
    validate_filename(filename, allowed_exts = "csv")
  }
  validate_filename(spreadsheet_xlsx, allowed_exts = "xlsx")


  # Setup -------------------------------------------------------------------

  # Check for existence of output files
  expected_per_contrast_tables <-  file.path(
    output_dir,
    contrasts_subdir,
    apply(
      X = expand.grid(names(DAList$results),
                      ".csv"),
      MARGIN = 1,
      FUN = paste0,
      collapse = ""
    )
  )
  summary_output_file <- file.path(output_dir, summary_csv)
  combined_output_file <- file.path(output_dir, combined_file_csv)
  excel_output_file <- file.path(output_dir, spreadsheet_xlsx)

  expected_files <- c(expected_per_contrast_tables,
                      summary_output_file,
                      combined_output_file,
                      excel_output_file)

  # If any files already exist
  if (any(file.exists(expected_files))) {
    if (!overwrite) {
      cli::cli_abort(c("Results files already exist",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Change {.arg output_dir} or set {.arg overwrite} to {.val TRUE}"))

    } else {
      cli::cli_inform("Results files already exist, and {.arg overwrite} == {.val {overwrite}}. Overwriting results files.")
      # Delete old results files so results don't become unsynced if there are any issues
      unlink(expected_files)
    }
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # Write summary CSV -------------------------------------------------------
  cli::cli_inform("Writing DA summary table to {.path {summary_output_file}}")

  summary <- do.call("rbind",
                     lapply(X = names(DAList$results),
                            FUN = summarize_contrast_DA,
                            contrast_res_list = DAList$results))
  summary$pval_thresh <- DAList$tags$DA_criteria$pval_thresh
  summary$lfc_thresh <- DAList$tags$DA_criteria$lfc_thresh
  summary$p_adj_method <- DAList$tags$DA_criteria$adj_method

  utils::write.csv(x = summary,
                   file = summary_output_file,
                   row.names = F)

  if (!file.exists(summary_output_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_output_file}}")) #nocov
  }


  # Write per-contrast csvs -------------------------------------------------
  per_contrast_dir <- file.path(output_dir, contrasts_subdir)
  cli::cli_inform("Writing per-contrast results {.path .csv}
                  {cli::qty(length(DAList$results))} file{?s} to {.path {per_contrast_dir}}")


  contrast_csv_success <- write_per_contrast_csvs(annotation_df = DAList$annotation,
                                                  data = DAList$data,
                                                  results_statlist = DAList$results,
                                                  output_dir = per_contrast_dir)
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success] #nocov start
    cli::cli_abort(c("Failed to write {.path .csv} results for
                     {cli::qty(sum(!contrast_csv_success))} contrast{?s}:",
                     "!" = "{.val {failed}")) #nocovend
  }


  # Write combined results csv ----------------------------------------------
  cli::cli_inform("Writing combined results table to {.path {combined_output_file}}")


  ## combine annotation, data, and per contrast results
  combined_results <- create_combined_results(annotation = DAList$annotation,
                                              data = DAList$data,
                                              statlist = DAList$results)

  ## write combined results csv
  combined_results_csv <- combined_results

  # format gene_symbol values so that excel does not convert them to date format
  # e.g. SEPT11, SEPT9, SEPT7, EPT2, SEPT6, MARC1, SEPT8, SEPT10, MARC2
  if ("gene_symbol" %in% colnames(combined_results_csv)) {
    combined_results_csv$gene_symbol <- paste0('"=""', combined_results_csv$gene_symbol,'"""')
  }

  csv_quote_cols <- which(colnames(combined_results_csv) != "gene_symbol")
  utils::write.csv(x = combined_results_csv,
                   file = combined_output_file,
                   row.names = F,
                   quote = csv_quote_cols)


  if (!file.exists(combined_output_file)) {
    cli::cli_abort(c("Failed to write combined results {.path .csv} to {.path {combined_output_file}}")) #nocov
  }


  # Write excel spreadsheet -------------------------------------------------
  cli::cli_inform("Writing combined results Excel spreadsheet to {.path {excel_output_file}}")



  write_limma_excel(filename = excel_output_file,
                    statlist = DAList$results,
                    annotation = DAList$annotation,
                    data = DAList$data,
                    norm.method = DAList$tags$norm_method,
                    pval_thresh = DAList$tags$DA_criteria$pval_thresh,
                    lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
                    add_filter = add_filter)

  if (!file.exists(excel_output_file)) {
    cli::cli_abort(c("Failed to write combined results Excel spreadsheet to {.path {excel_output_file}}")) #nocov
  }

  # If everything works, return input DAList
  invisible(input_DAList)
}



#' Summarize the number of DA proteins in a contrast
#'
#' Internal function to summarize the number of DE genes/proteins for a given
#' contrast.
#'
#' @param contrast_name The name of the contrast to summarize.
#' @param contrast_res_list A list of per-contrast DA results.
#'
#' @return A data frame summarizing differential abundance for the given contrast.
#'
#' @keywords internal
#'
summarize_contrast_DA <- function(contrast_name, contrast_res_list) {
  tmp <- contrast_res_list[[contrast_name]][,c("sig.PVal", "sig.FDR")]

  data.frame(cbind(contrast = contrast_name,
                   type = c("down", "nonsig", "up"),
                   rbind(colSums(tmp == -1, na.rm = T),
                         colSums(tmp == 0, na.rm = T),
                         colSums(tmp == 1, na.rm = T))))

}


#' Write per-contrast results csvs
#'
#' Internal function used to write per-contrast statistical results to a .csv
#' file.
#'
#' @param annotation_df A data frame of annotation data for each gene/protein.
#' @param data A data frame of average expression data for each sample.
#' @param results_statlist A list of per-contrast DE results.
#' @param output_dir The directory in which to save the per-contrast csv files.
#'
#' @keywords internal
#'
#' @return A logical vector indicting whether each contrast file was
#'   successfully written.
#'
write_per_contrast_csvs <- function(annotation_df,
                                    data,
                                    results_statlist,
                                    output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # format annotation_df GN column so that excel does not convert them to date format
  # e.g. SEPT11, SEPT9, SEPT7, EPT2, SEPT6, MARC1, SEPT8, SEPT10, MARC2
  annotation_df$gene_symbol <- paste0('"=""',annotation_df$gene_symbol,'"""')
  csv_quote_cols <- which(colnames(annotation_df)!="gene_symbol")

  # make list combining
  # annotation, intensity data, and statistical results
  # reordering everything by data rownames to match and subset
  # Note: technically, these could be in different orders across contrasts,
  # If they differed across contrasts.
  per_contrast_results <- lapply(X = results_statlist,
                                 FUN = function(x) cbind(annotation_df[rownames(data), , drop = F],
                                                         data[rownames(data), , drop = F],
                                                         x[rownames(data), ,drop = F]))
  # Make corresponding filenames
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  # write
  tmp <- lapply(seq_along(per_contrast_results),
                function(x)
                {
                  utils::write.csv(per_contrast_results[[x]],
                                   file = filenames[x],
                                   row.names = F, quote = csv_quote_cols)
                }
  )

  # check
  write_success <- file.exists(filenames)
  names(write_success) <- names(per_contrast_results)

  write_success
}



#' Combine statistical results
#'
#' Internal function used to construct a data frame of combined results
#'
#' @param annotation A data frame of annotation data for each gene/protein.
#' @param data A data frame of normalized intensity data for each sample.
#' @param statlist A list of per-contrast DE results.
#'
#' @return A data frame of the combined results.
#'
#' @keywords internal
#'
create_combined_results <- function(annotation,
                                    data,
                                    statlist) {
  # Get common row order from data
  row_order <- rownames(data)
  # Reorder rows and rename cols for each element of the results statlist
  results_for_combine <- lapply(X = names(statlist),
                                function(x) {
                                  tmp <- statlist[[x]][row_order, , drop = F]
                                  colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
                                  tmp
                                })

  ## combine, annotation, data and combined stat results
  combined_results <- cbind(annotation[row_order, , drop = F],
                            data[row_order, , drop = F],
                            results_for_combine)

  combined_results
}


#' Write an .xlsx file of limma results
#'
#' Internal function for creating a nicely formatted Excel spreadsheet of
#' differential expression results.
#'
#' @param filename The filename of the Excel spreadsheet to be saved.
#' @param statlist A list of the per-contrast statistical results.
#' @param annotation A data frame of annotation date for each protein.
#' @param data A data frame containing the average expression data for each sample.
#' @param norm.method The method that was used to normalize the data for the
#'   statistical model being output.
#' @param pval_thresh The p-value threshold that was used to determine significance.
#' @param lfc_thresh The logFC threshold that was used to determine significance.
#' @param add_filter Should per-column filters be added to the spreadsheet?
#'
#' @return Invisibly returns a list, where the first element is the filename
#'   of the saved Excel spreadsheet and the second element is the openxlsx
#'   workbook object.
#'
#' @keywords internal
#'
write_limma_excel <- function(filename, statlist, annotation, data, norm.method,
                              pval_thresh, lfc_thresh, add_filter) {


  # Maybe some argument processing
  normName <- names(norm.methods)[grep(norm.method,norm.methods)]
  row_order <- rownames(data)

  # initialize workbook
  wb<-openxlsx::createWorkbook()

  # Set up annotation columns -----------------------------------------------
  newNames <- stringr::str_replace_all(colnames(annotation), pattern = "uniprot_id", replacement = "UniProt ID") |>
    stringr::str_remove_all("[._]")

  annot.title <- "Protein Annotation"
  data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Exclusive Intensities")
  sheetName <- "Protein Results"

  # Initialize worksheet ----------------------------------------------------
  openxlsx::addWorksheet(wb, sheetName = sheetName, tabColour = "#FFE100")
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 3, heights = 30) # increase row height for the title row
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 4, heights = 20) # increase row height for the header row
  title_row <- 3

  # Add annotation columns -----------------------------------------------------

  # Subset and rename
  annot <- annotation[row_order, , drop = F]
  colnames(annot) <- newNames
  annot <- make_excel_hyperlinks(data = annot,
                                 url.col = "UniProt ID",
                                 url = "https://www.uniprot.org/uniprot/")

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
  if (length(names(statlist)) > length(binfcolors)) {
    numColors<- ceiling(length(names(statlist))/length(binfcolors)) + 1
    colors2 <- rep(binfcolors, numColors)
    lightcolors2 <- rep(lightbinfcolors, numColors)
  } else {
    if (length(names(statlist)) <= length(binfcolors)) {
      colors2 <- binfcolors
      lightcolors2 <- lightbinfcolors
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


#' Make hyperlinks for an excel column
#'
#' Internal functions to add hyperlinks to a column in a data frame to
#' be exported to an excel file. Implicitly, works only on one column.
#'
#' @param data The data frame in which to add hyperlinks.
#' @param url.col The name of the column to add hyperlinks to.
#' @param url The url that will be prepended to the info in the column.
#'
#' @return The original data frame, now with hyperlinks in the desired column.
#'
#' @keywords internal
#'
make_excel_hyperlinks <- function(data, url.col, url) {

  ids <- data[,url.col, drop = T]
  tmp <- is.na(ids)
  url2 <- paste0(url, ids)
  ids2 <- paste0("HYPERLINK(\"", url2, "\", \"", ids, "\")")
  ids2[tmp] <- NA
  data[,url.col] <- ids2
  class(data[,url.col]) <- "formula"

  return(data)
}
