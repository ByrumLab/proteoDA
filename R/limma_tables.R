write_limma_results <- function(model_results,
                                annotation,
                                ilab,
                                out_dir = file.path(ilab, paste0(enrich, "_analysis")),
                                enrich = c("protein", "phospho"),
                                contrasts_subdir = "per_contrast_results",
                                summary_csv = "DE_summary.csv",
                                combined_file_csv = "combined_results.csv",
                                BQ_csv = paste0(ilab, "_results_BQ.csv"),
                                spreadsheet_xlsx = NULL) {

  ##########################
  ## Setup and arg checks ##
  ##########################

  enrich <- rlang::arg_match(enrich)

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }

  # Add filename validation

  # Add checks of model result rownames in annotation rownames.

  statlist <- model_results$stats_by_contrast
  data <- model_results$data

  #######################
  ## Write summary csv ##
  #######################

  model_results <- results_reb

  summary <- lapply(X = names(statlist),
                    FUN = summarize_contrast_DE,
                    contrast_res_list = statlist) %>%
    do.call("rbind", .)
  summary$pval_thresh <- model_results$min.pval
  summary$lfc_thresh <- model_results$min.lfc
  summary$p_adj_method <- model_results$adj.method

  summary_ouput_file <- file.path(out_dir, summary_csv)

  utils::write.csv(x = summary,
                   file = summary_ouput_file,
                   row.names = F)

  if (!file.exists(summary_ouput_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_ouput_file}}"))
  }


  #############################
  ## Write per-contrast csvs ##
  #############################

  contrast_csv_success <- write_per_contrast_csvs(annotation_df = annotation,
                                                  data = data,
                                                  results_statlist = statlist,
                                                  output_dir = file.path(out_dir, contrasts_subdir))
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success]
    cli::cli_abort(c("Failed to write {.path .csv} results for {cli::qty(sum(!contrast_csv_success))} contrast{?s}:",
                     "!" = "{.val {failed}"))
  }

  ################################
  ## Write combined results csv ##
  ################################

  # Get common row order from the first element of the stat list
  row_order <- rownames(statlist[[1]])
  # Reorder rows and rename cols for each element of the results statlist
  results_for_combine <- lapply(X = names(statlist),
                                function(x) {
                                  tmp <- statlist[[x]][row_order, ,drop = F]
                                  colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
                                  tmp
                                })
  combined_output_file <- file.path(out_dir, combined_file_csv)
  combined_results <- cbind(annotation[row_order, ],
                            data[row_order, ],
                            results_for_combine)
  utils::write.csv(x = combined_results,
                   file = combined_output_file,
                   row.names = F)

  if (!file.exists(combined_output_file)) {
    cli::cli_abort(c("Failed to write combined results {.path .csv} to {.path {combined_output_file}}"))
  }

  #################################
  ## Write combined BiqQuery csv ##
  #################################

  BQ_output_file <- file.path(out_dir, BQ_csv)
  BQ_output <- combined_results
  BQ_output[is.na(BQ_output)] <- 0
  BQ_output[BQ_output == ""] <- 0

  # colnames can't begin with numbers?
  col_beings_number <- stringr::str_detect(colnames(BQ_output), "^[:digit:]")
  if (any(col_beings_number)) {
    colnames(BQ_output)[col_beings_number] <- paste0("X", colnames(BQ_output)[col_beings_number])
  }
  utils::write.csv(BQ_output, file = BQ_output_file, row.names = F)

  if (!file.exists(BQ_output_file)) {
    cli::cli_abort(c("Failed to write BigQuery results {.path .csv} to {.path {BQ_output_file}}"))
  }

  #############################
  ## Write excel spreadsheet ##
  #############################

  # TODO: deal with this.

  ############
  ## Finish ##
  ############

  # If everything works, return combined results
  invisible(combined_results)
}


summarize_contrast_DE <- function(contrast_name, contrast_res_list) {
  tmp <- contrast_res_list[[contrast_name]][,c("sig.PVal", "sig.FDR")]

  data.frame(cbind(contrast = contrast_name,
                   type = c("down", "nonsig", "up"),
                   rbind(colSums(tmp == -1, na.rm = T),
                         colSums(tmp == 0, na.rm = T),
                         colSums(tmp == 1, na.rm = T))))

}


write_per_contrast_csvs <- function(annotation_df,
                                    data,
                                    results_statlist,
                                    output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # make list combining
  # annotation, intensity data, and statistical results
  # reordering everything by rownames to match and subset
  # Note: technically, these could be in different orders across contrasts,
  # If they differed across contrasts.
  per_contrast_results <- lapply(X = results_statlist,
                                 FUN = function(x) cbind(annotation_df[rownames(x), ],
                                                         data[rownames(x), ],
                                                         x))
  # Make corresponding filenames
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  # write
  tmp <- lapply(seq_along(per_contrast_results),
                function(x)
                  {
                  utils::write.csv(per_contrast_results[[x]],
                                   file = filenames[x],
                                   row.names = F)
                  }
                )

  # check
  write_success <- file.exists(filenames)
  names(write_success) <- names(per_contrast_results)

  write_success
}

<<<<<<< HEAD








add_limma_results <- function(wb, sheetName=NULL, statList, annot, data, norm.method,
                              min.pval, min.lfc, pipe, enrich){

  ## merge column colors
  binfcolors <-c("#1F5EDC","#EE0010","#32CD32","#FF1493","#FF7F00",
                 "#A342FC","#00C8FF","#ADFF2F","#FFE100","#E36EF6","#009ACE","#996633")
  names(binfcolors)<-c("blueberry","cherry","apple","barbie","fanta",
                       "grape","ocean","mtndew","gold","orchid","aceblue","poop")

  ## column heading colors
  lightbinfcolors<-c("#c8d8f7","#ffdadd","#d0f3d0","#ffb1db","#ffe1c4",
                     "#dbb6fe","#c4f2ff","#e3ffb8","#fff8c4","#f1b8fb",
                     "#baeeff","#eddcca")
  names(lightbinfcolors)<-c("lightblueberry","lightcherry","lightapple","lightbarbie","lightfanta",
                            "lightgrape","lightocean","lightmtndew","lightgold","lightorchid",
                            "lightaceblue","lightpoop")


  normName <- names(norm.methods)[grep(norm.method,norm.methods)];normName

  if(pipe=="DIA" & enrich=="protein"){

    annotCols=c("id", "Protein.Name","Accession.Number","Molecular.Weight", "Protein.Group.Score",
                "Identified.Peptide.Count","Exclusivity",
                "UniprotID", "Gene_name","Description")
    newNames=c("id", "Protein Name","Accession Number","Molecular Weight", "Protein Group Score",
               "Identified Peptide Count","Exclusivity",
               "UniProt ID", "Gene Name","Description")
    annot.title <- "Protein Annotation"
    data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Intensities");data.title

    if(is.null(sheetName)){ sheetName <- "Protein Results" }
    ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
    if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
    print(sheetName)
  }

  if(pipe=="TMT" & enrich=="protein"){
    annotCols= c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",
                 "Gene_name","Description")
    newNames= c("id", "Fasta headers","Majority protein IDs", "Score", "UniProt ID",
                "Gene Name","Description")
    annot.title <- "Protein Annotation"
    data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
    # sheetName   <- "Protein Results"

    if(is.null(sheetName)){ sheetName <- "Protein Results" }
    ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
    if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
    print(sheetName)

  }


  if(pipe=="LF" & enrich=="protein"){
    annotCols= c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",
                 "Gene_name","Description")
    newNames= c("id", "Fasta headers","Majority protein IDs", "Score", "UniProt ID",
                "Gene Name","Description")
    annot.title <- "Protein Annotation"
    data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Intensities");data.title
    # sheetName   <- "Protein Results"

    if(is.null(sheetName)){ sheetName <- "Protein Results" }
    ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
    if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
    print(sheetName)
  }


  if(pipe=="phosphoTMT" & enrich=="protein"){

    annotCols= c("id", "Fasta.headers", "Majority.protein.IDs", "Phospho..STY..site.IDs", "Score",
                 "UniprotID","Gene_name","Description")
    newNames= c("id", "Fasta headers", "Majority protein IDs", "Phospho (STY) site IDs", "Score",
                "UniProt ID","Gene Name","Description")
    annot.title <- "Protein Annotation"
    data.title   <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
    # sheetName   <- "Protein Results"

    if(is.null(sheetName)){ sheetName <- "Protein Results" }
    ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
    if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
    print(sheetName)

  }

  if(pipe=="phosphoTMT" & enrich=="phospho"){

    annotCols <- c("id","Fasta.headers", "proGroupID","UniprotID","Gene_name",
                   "Description", "PEP","Score", "Localization.prob", "Phospho..STY..Probabilities",
                   "Flanking","phosAAPosition","Class")
    newNames <- c("id","Fasta headers", "Protein Group ID","UniProt ID","Gene Name",
                  "Description", "PEP","Score", "Localization prob", "Phospho (STY) Probabilities",
                  "Flanking","Site","Class")

    annot.title <- "Protein/Phospho Annotation"
    data.title  <- paste("Log2",ifelse(normName=="Log2","",normName),"Normalized Exclusive MS1 Intensities");data.title
    # sheetName   <- "Phospho Results"

    if(is.null(sheetName)){ sheetName <- "Phospho Results" }
    ifelse(nchar(sheetName)>30,substring(sheetName,1,30), sheetName)
    if(sheetName %in% names(wb)){ sheetName=paste("Sheet",(length(wb$worksheets)+1),sep="") };sheetName
    print(sheetName)

  }


  ## get subset of annotation columns and rename them
  annot <- annot[, annotCols];head(annot)
  colnames(annot) <- newNames;head(annot)
  annot <- make_excel_hyperlinks(data=annot, url.col="UniProt ID", url="https://www.uniprot.org/uniprot/")


  ## add worksheet to workbook
  openxlsx::addWorksheet(wb=wb, sheetName=sheetName, tabColour="#FFE100")
  openxlsx::setRowHeights(wb=wb, sheet=sheetName, rows = 3, heights = 30);## increase row height for the header
  openxlsx::setRowHeights(wb=wb, sheet=sheetName, rows = 4, heights = 20);## increase row height for the header


  ##--------------
  ##  ANNOTATION
  ##--------------
  ## merge annotation columns
  annot.start=1
  annot.end=ncol(annot);annot.end

  openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=annot.start:annot.end, rows=3)
  openxlsx::writeData(wb=wb,sheet=sheetName, x=annot.title, startCol=1, startRow=3, colNames=FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill="#636262", halign="center", valign="center", textDecoration="bold",
                                      fontColour="white", fontSize=14, numFmt="TEXT")
  openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=1, stack=TRUE)


  ## write annotation data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#d8d8d8", halign="left", valign="center", textDecoration="bold",
                                    fontColour="black", fontSize=11, numFmt="TEXT")
  openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=annot.start:annot.end, stack=TRUE)
  openxlsx::writeData(wb=wb, sheet=sheetName, x=annot,
                      startCol = annot.start,
                      startRow = 4,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders="columns",
                      keepNA=TRUE, na.string="NA",
                      sep="\t")

  ## add hyperlink style (blue font, underlined)
  hyperlinkStyle <- openxlsx::createStyle(fontColour="#0000FF",halign="left",valign="center",textDecoration="underline")
  cols <- grep("UniProt ID",colnames(annot));cols
  openxlsx::addStyle(wb=wb, sheet=sheetName, style=hyperlinkStyle, cols=cols, rows=5:(nrow(annot)+4),
                     gridExpand=TRUE,stack=TRUE)



  ##--------------
  ##  DATA
  ##--------------
  ## start/stop column positions for data
  data.start=annot.end+1;data.start
  data.end=annot.end+ncol(data);data.end

  ## data merged column
  openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=data.start:data.end, rows=3)
  openxlsx::writeData(wb=wb,sheet=sheetName, x=data.title, startCol=data.start, startRow=3, colNames=FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill = "#516285", halign="center", valign="center", textDecoration="bold",
                                      fontColour="white", fontSize=14, numFmt="TEXT")
  openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=data.start, stack=TRUE)

  ## write data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#cdd4e1", halign="left", valign="center", textDecoration="bold",
                                    fontColour="black", fontSize=11, numFmt="TEXT")
  openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=data.start:data.end, stack=TRUE)
  openxlsx::writeData(wb=wb, sheet=sheetName, x=data,
                      startCol = data.start,
                      startRow = 4,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders="columns",
                      keepNA=TRUE, na.string="NA",
                      sep="\t")


  ##--------------
  ##  STATS
  ##--------------
  ## colors for stat columns
  if(length(names(statList)) > length(binfcolors)){
    numColors<- ceiling(length(names(statList))/length(binfcolors))+1
    colors2 <- rep(binfcolors,numColors)
    lightcolors2 <- rep(lightbinfcolors,numColors)
  } else {
    if(length(names(statList))<=length(binfcolors)){
      colors2<-binfcolors
      lightcolors2<-lightbinfcolors
    }
  }
  stopifnot(length(colors2) > length(names(statList)))


  ## write stat data to worksheet
  stat.start=NULL;stat.end=NULL
  for(i in base::seq_along(names(statList))){
    stats <- statList[[i]]
    comparison <- names(statList)[i];comparison
    mycolor<-colors2[i];mycolor
    mylightcolor<-lightcolors2[i];mylightcolor

    ## stats start
    if(is.null(stat.start)){
      stat.start=data.end+1;stat.start
    } else {
      if(!is.null(stat.start)){
        stat.start=stat.end+1;stat.start }
    }

    ## stats end
    if(is.null(stat.end)){
      stat.end=data.end+ncol(stats);stat.end
    } else {
      if(!is.null(stat.end)){
        stat.end=stat.end+ncol(stats);stat.end }
    }


    openxlsx::mergeCells(wb=wb,sheet=sheetName, cols=stat.start:stat.end, rows=3)
    openxlsx::writeData(wb=wb,sheet=sheetName, x=comparison, startCol=stat.start, startRow=3, colNames=FALSE)
    mergeStyle <- openxlsx::createStyle(fgFill = mycolor, halign="center", valign="center", textDecoration="bold",
                                        fontColour="white", fontSize=14, numFmt="TEXT")
    openxlsx::addStyle(wb=wb,sheet=sheetName, style=mergeStyle, rows=3, cols=stat.start, stack=TRUE)


    ## write stats to worksheet
    colStyle <- openxlsx::createStyle(fgFill = mylightcolor, halign="left", valign="center", textDecoration="bold",
                                      fontColour="black", fontSize=11, numFmt="TEXT")
    openxlsx::addStyle(wb=wb,sheet=sheetName, style=colStyle, rows=4, cols=stat.start:stat.end, stack=TRUE)
    openxlsx::writeData(wb=wb, sheet=sheetName, x=stats,
                        startCol = stat.start,
                        startRow = 4,
                        colNames = TRUE,
                        rowNames = FALSE,
                        borders="columns",
                        keepNA=TRUE, na.string="NA",
                        sep="\t")


    posStyle  <- openxlsx::createStyle(fontColour="#990000", bgFill="#FFC7CE")
    negStyle  <- openxlsx::createStyle(fontColour="#006100", bgFill="#C6EFCE")
    normStyle <- openxlsx::createStyle(fontColour="#000000", bgFill="#FFFFFF")

    fc.col    <- grep("logFC", colnames(stats))+stat.start-1;fc.col
    fc.rule1  <- paste0(">=",min.lfc);fc.rule1
    fc.rule2  <- paste0("<=",-min.lfc);fc.rule2
    fdr.col   <- grep("adj.P.Val", colnames(stats))+stat.start-1; fdr.col
    fdr.rule  <- paste0("<=",min.pval);fdr.rule
    pval.col  <- grep("P.Value", colnames(stats))+stat.start-1;pval.col
    pval.rule <- paste0("<=",min.pval);pval.rule

    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=fc.rule1, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=fc.rule2, style=negStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fc.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=fdr.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=pval.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)


  } ## FOR LOOP


  ## add global styles
  headerStyle <- openxlsx::createStyle(halign = "center", valign="center", #fontColour = "#000000",
                                       border = "TopBottomLeftRight", borderColour = c("black","black","black","black"),
                                       textDecoration = "bold")
  openxlsx::addStyle(wb=wb, sheet=sheetName, style=headerStyle, cols=1:stat.end, rows=3, gridExpand=TRUE, stack=TRUE)
  openxlsx::addStyle(wb=wb, sheet=sheetName, style=headerStyle, cols=1:stat.end, rows=4, gridExpand=TRUE, stack=TRUE)
  openxlsx::addFilter(wb=wb, sheet=sheetName, cols=1:stat.end, rows=4)


} ## ADD LIMMA RESULTS



}

=======
>>>>>>> main
