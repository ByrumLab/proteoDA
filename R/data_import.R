
#' Loading Maxquant data
#'
#' A function to extract and load MaxQuant data. Works across multiple
#' data types (DIA, TMT, etc.), and does slightly different things for each one.
#' The \code{extract_data} function is the main function, and it calls a number
#' of other subfunctions to process and prepare the data. See the links to
#' subfunctions for more info. Returns a list of data from maxquant, along with
#' some stats and other info. Also has some side effects: creates a BigQuery
#' report and some log files.
#'
#'
#' Subfunctions (DIA): \code{\link{import_data}}, \code{\link{remove_contaminants}},
#' \code{\link{extract_protein_data}}, \code{\link{make_log}}.
#'
#' @param file The file of Maxquant data to be imported. Format and filetype
#'   vary depending on the type of data to be analyzed. If not supplied, R will
#'   ask you to pick a file from the file browser.
#' @param pipe Which analysis pipeline are you running? Options are "DIA", "TMT",
#'   "phosphoTMT", and "LF". Output varies slightly based on the pipeline you
#'   choose.
#' @param enrich Another aspect of pipeline? Options are "protein" and "phospho".
#'   Seems semi-redundant?: DIA, TMT, and LF seem to always be "protein" in the
#'   code. PhosphoTMT can be either protein or phospho, though that seems to vary
#'   across subfunctions a little?
#' @param sampleIDs Optional: a list of sampleIDs to analyze. Default is NULL, in which
#'   case all sample IDs are analyzed. These are determined from column names.
#' @param min.prob Optional: the minimum probability, for use in local probability filter of
#'   phosphoTMT data. Unused for other pipelines. Default is 0.75.

#'
#' @return A list of extracted data. The structure of the list depends on the
#'   pipeline and enrichment used. For DIA/protein, a list with 5 elements:
#'   \enumerate{
#'     \item "data", a dataframe with N rows and M cols, where N is proteins and M is
#'       individual sample intensities. Rownames in this data frame give a protein ID,
#'       in the form of UniprotID_Gene_name_id. These
#'       proteins have been partially filtered: contaminants removed.
#'     \item "annot". a dataframe with N rows and X cols, where N is proteins and the X
#'       cols are various protein annotation data. Rownames in this data frame give a protein ID,
#'       in the form of UniprotID_Gene_name_id.
#'     \item "param", a dataframe with 1 column giving the value of some of the
#'       parameters used to call the function, but also other info it parses from the
#'       file. Rownames give the info on what each value is.
#'     \item  "stats", a named list of statistics on sample #, contaminants
#'       removed, etc.
#'     \item a character vector of the ilab variable (also listed in the
#'     param dataframe).
#'   }
#' @export
#'
#' @examples
#' # No examples yet
#'
extract_data <- function(file = NULL,
                         pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                         enrich = c("protein", "phospho"),
                         sampleIDs = NULL,
                         min.prob = 0.75) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)


  ## if no file is input by user the open file dialog box is
  ## called to allow the user to navigate to and select the file.
  if (is.null(file)) {
    file <- file.choose()
  }

  cli::cli_rule()

  ## input parameters are added to param variable
  param <- stats <- list()
  param[["file"]] <- file
  param[["pipe"]] <- pipe
  param[["enrich"]] <- enrich
  param[["sampleIDs"]] <- ifelse(is.null(sampleIDs), "NULL", paste(sampleIDs, collapse = ", "))
  param[["min.prob"]] <- ifelse(enrich == "protein", "NULL", min.prob)


  ## DIA PROTEIN
  if (pipe == "DIA" & enrich == "protein") {


    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    ## the first 10 lines are skipped and last line removed. column 1 (.X)
    ## is renamed "id"
    DIADATA <- import_data(file = file, pipe = pipe, enrich = enrich)
    diaData <- DIADATA$data
    reqCols <- DIADATA$reqCols

    param[["file"]] <- DIADATA$file
    stats[["total_input_rows"]] <- DIADATA$num_input

    ## REMOVE IRRELEVANT COLUMNS
    ## irrelevant columns, columns 2-3 (visible,star) are removed from
    ## the sample report
    if (colnames(diaData)[1] != "id") {
      colnames(diaData)[1] <- "id"
    }
    remove <- colnames(diaData) %in% c("Visible", "Star")
    diaData <- diaData[, remove == FALSE]


    ## REPLACE MISSING VALUES WITH ZERO
    ## intensity data is not a numeric data.frame and requires
    ## additional processing numbers are character strings separated
    ## by commas (4,000,000), some cells contain text 'Missing Value'
    ## if the protein does not have detectable expression.
    ## replace 'Missing Value' with zeros, remove comma's, convert
    ## data in each column to numeric values
    diaData[, ][diaData[, ] == "Missing Value"] <- "0"
    replaceNumbers <- colnames(diaData) %in% reqCols
    replaceColums <- colnames(diaData)[replaceNumbers == FALSE]
    for (i in replaceColums) {
      diaData[, i] <- remove_commas(diaData[, i])
    }


    ## SAVE BIG QUERY INPUT FILE
    ## save a copy of samples report for upload into Big Query.
    ## This file cannot contain NA or blank values, column names
    ## cannot begin with a number, replace all NA and blank cells
    ## with zeros. If column names begin with a number add X to
    ## beginning of column name. change file name to ilab_Samples_Report_BQ.csv
    bqData <- data.frame(diaData)
    bqData[, ][is.na(bqData[, ])] <- 0
    bqData[, ][bqData[, ] == ""] <- 0

    ## if column names begin with a number append X to the start of each name.
    testColumNumber <- substr(colnames(bqData), 1, 1)
    if (length(grep("[[:digit:]]", testColumNumber)) > 0) {
      colnames(bqData) <- paste0("X", colnames(bqData))
    }
    bn <- basename(file)
    bn <- gsub("Samples Report of ", "", bn)
    ilab <- gsub(paste0(".", file_extension(bn)), "", bn)
    ## save DIA Big Query Input File
    filename <- paste0(ilab, "_Samples_Report_BQ.csv")
    utils::write.csv(bqData, file.path(".", filename), row.names = FALSE)

    cli::cli_inform("BigQuery samples report saved to {.file {filename}}")
    param[["ilab"]] <- ilab


    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data = diaData, pipe = pipe, enrich = enrich)
    stats[["num_contam_removed"]] <- qfilterData$num_contam

    ## EXTRACT PROTEIN DATA
    extProt <- extract_protein_data(
      data = qfilterData$data, sampleIDs = sampleIDs,
      pipe = pipe, enrich = enrich
    )

    stats[["total_input_samples"]] <- extProt$num_samples
    stats[["num_extracted_rows"]] <- extProt$num_extracted

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param = param, stats = stats, title = "EXTRACTED DATA", save = TRUE)

    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(
      data = extProt$rawData, annot = extProt$rawAnnot,
      param = logs$param, stats = logs$stats[c(3, 1, 2, 4), ], ilab = ilab
    )
    return(data2)
  } ## DIA PROTEIN


  ## TMT PROTEIN
  if ((pipe == "TMT" | pipe == "LF") & enrich == "protein") {

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file = file, pipe = pipe, enrich = enrich)
    stats[["total_input_rows"]] <- maxQuant$num_input

    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data = maxQuant$data, pipe = pipe, enrich = enrich)
    stats[["num_contam_removed"]] <- qfilterData$num_contam


    ## EXTRACT PROTEIN DATA AND ANNOTATION
    extProt <- extract_protein_data(
      data = qfilterData$data, sampleIDs = sampleIDs,
      pipe = pipe, enrich = enrich
    )

    stats[["total_input_samples"]] <- extProt$num_samples
    stats[["num_extracted_rows"]] <- extProt$num_extracted
    param[["pattern"]] <- extProt$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param = param, stats = stats, title = "EXTRACTED DATA", save = TRUE)

    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(
      data = extProt$rawData, annot = extProt$rawAnnot,
      param = logs$param, stats = logs$stats
    )
    return(data2)
  } ## TMT PROTEIN


  ## phosphoTMT PROTEIN
  if (pipe == "phosphoTMT" & enrich == "protein") {

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file = file, pipe = pipe, enrich = enrich)
    stats[["total_input_rows"]] <- maxQuant$num_input


    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data = maxQuant$data, pipe = pipe, enrich = enrich)
    stats[["num_contam_removed"]] <- qfilterData$num_contam


    ## EXTRACT PROTEIN DATA
    extProt <- extract_protein_data(
      data = qfilterData$data, sampleIDs = sampleIDs,
      pipe = pipe, enrich = enrich
    )

    stats[["total_input_samples"]] <- extProt$num_samples
    stats[["num_extracted_rows"]] <- extProt$num_extracted
    param[["pattern"]] <- extProt$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param = param, stats = stats, title = "EXTRACTED DATA", save = TRUE)


    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(
      data = extProt$rawData, annot = extProt$rawAnnot,
      param = logs$param, stats = logs$stats[c(3, 1, 2, 4), ]
    )
    return(data2)
  } ## phosphoTMT PROTEIN


  ## phosphoTMT PHOSPHO
  if (pipe == "phosphoTMT" & enrich == "phospho") {

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("phospho_analysis/01_quality_control")){
    #    dir.create("phospho_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file = file, pipe = pipe, enrich = enrich)
    stats[["total_input_rows"]] <- maxQuant$num_input

    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data = maxQuant$data, pipe = pipe, enrich = enrich)
    stats[["num_contam_removed"]] <- qfilterData$num_contam

    ## LOCAL PROBABILITY FILTER
    qfilterData <- local_prob_filter(
      data = qfilterData$data, min.prob = min.prob,
      pipe = pipe, enrich = enrich
    )
    param[["min.prob"]] <- qfilterData$min.prob
    stats[["num_localprob_removed"]] <- qfilterData$num_localprob_removed

    ## EXTRACT PHOSPHO DATA
    extPhos <- extract_phospho_data(
      data = qfilterData$data, sampleIDs = sampleIDs,
      pipe = pipe, enrich = enrich
    )
    stats[["total_input_samples"]] <- extPhos$classList$class_1$num_class_samples
    stats[["num_extracted_rows"]] <- extPhos$classList$class_1$num_class_phospho
    param[["pattern"]] <- extPhos$classList$class_1$idList$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param = param, stats = stats, title = "EXTRACTED DATA", save = TRUE)


    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data for class 1 phosphopeptides,
    ## corresponding extracted annotation, stats, and input parameters
    ## classList corresponds to data/annotation for each phosphopeptide
    ## class e.g. class___1, class___2, class___3

    data2 <- list(
      data = extPhos$classList$class_1$rawData, ## class 1 data
      annot = extPhos$classList$class_1$rawAnnot, ## class 1 annot
      classList = extPhos$classList, ## class_1$data/annot, class_2$data/annot etc.
      sampleIDs = extPhos$classList$class_1$idList$sampleIDs,
      pattern = extPhos$classList$class_1$idList$pattern,
      param = logs$param, stats = logs$stats[c(4, 1, 2, 3, 5), ]
    )
    return(data2)
  } ## phosphoTMT PHOSPHO
}



#' Import MaxQuant data
#'
#' First subfunction called by \code{\link{extract_data}}. Reads in the tabular
#' data and does some checks for required columns.
#'
#' @inheritParams extract_data
#'
#' @return A list. Items in list vary depending on the pipeline used. For
#'   DIA/protein, a list with 3 slots:
#'   \enumerate{
#'      \item "data": a data frame, where each row is a (raw) protein, the first
#'        few columns are protein ID data, and the last columns are the individual
#'        sample intensities (as character vectors).
#'     \item "num_input", a named integer listing the total number of raw proteins (ie,
#'        the number of rows of the data slot)
#'     \item a character vector of the required columns for this pipeline/enrichment
#'        combination.
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet


import_data <-function(file,
                       pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                       enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)


  ## check that the file is a csv, tsv, or text file
  filext <- file_extension(file)
  if (filext %notin% c("csv","txt","tsv")) {
    cli::cli_abort(c("Problem with input file",
                     "x" = "Input file must end in {.file .csv}, {.file .tsv}, or {.file .txt}"))
  }

  ## check that file exists
  if (!file.exists(file)) {
    cli::cli_abort(c("Cannot find file.",
                     "x" = "{.file {file}} does not exist",
                     "i" = "Did you specify the filepath correctly?"))
  }

  if(filext=="txt" | filext=="tsv"){ sep="\t" }
  if(filext=="csv"){ sep=","}

  ## vector of required column names based on analysis pipeline
  ## and type of enrichment. these columns are processed by other
  ## functions that extract annotation, remove contaminants, etc.
  ## values for each variable are defined according to R syntax
  ## rules. Values for each variable are listed at the top of the
  ## functions.R file.
  if(all(pipe=="DIA" & enrich=="protein")){
    reqCols<-unique(c(diaAnnotationColums,diaContamColums))
  }
  if(all(pipe=="TMT" & enrich=="protein")){
    reqCols<-c(tmtAnnotationColums, tmtContamColums)
  }
  if(all(pipe=="phosphoTMT" & enrich=="protein")){
    reqCols <- c(proteinAnnotationColums, proteinContamColums)
  }
  if(all(pipe=="phosphoTMT" & enrich=="phospho")){
    reqCols <- c(phosphoAnnotationColums, phosphoContamColums, phosphoLocalProbColum)
  }
  if(all(pipe=="LF" & enrich=="protein")){
    reqCols<-c(tmtAnnotationColums, tmtContamColums)
  }

  ## IMPORT DIA DATA
  if(pipe=="DIA"){

    ## sample report generated by scaffold DIA is imported as data.frame.
    ## Sample Reports generated by Scaffold DIA contain meta data in the
    ## first 10 lines, and text in the last row to mark the end of the file.
    ## Thus, when reading in the file the first 10 lines are skipped and the
    ## last row is removed. column name of column 1 (X. on import) is changed to id
    ## column 1 must be changed to "id" to allow extraction and other processing.
    ## Note: column names of the input file are defined according to R syntax rules.
    ## special characters, spaces, are converted to periods.
    data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE,
                            header=TRUE, check.names=TRUE, skip=10)
    if(data[nrow(data),1]=="END OF FILE"){ data <- data[-nrow(data),] }

    ## name of column 1 is changed to 'id'
    if(colnames(data)[1]=="X."){ colnames(data)[1] <- "id" }
    data <- data.frame(data)
    num_input <- nrow(data)

    ## checks to make sure the input file contains the required annotation/contamination
    ## columns. This ensures that the file imported correctly and helps to identify
    ## sampleIDs columns and extract targets info. later.
    if(all(reqCols %in% colnames(data))==TRUE){

      cli::cli_inform("{pipe} {enrich} file {.file {file}} imported")
      cli::cli_inform("Input file contains {num_input} {enrich} entries")

      data2<-list(data=data, num_input=num_input, reqCols=reqCols)
      return(data2)

    } else {
      missing_cols <- reqCols[reqCols %notin% colnames(data)]
      cli::cli_abort(c("Problem with columns in input file.",
                     "x" = "Required column(s) not present in {.file {file}}",
                     "x" = "Missing column{?s}: {missing_cols}"))

    }
  } ## DIA IMPORT


  ## IMPORT TMT/phosp DATA
  if(any(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF")){

    ## imports file as data.frame. Note: column names of the input file
    ## are defined according to R syntax rules. special characters, spaces,
    ## are converted to periods.
    data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE,
                            header=TRUE, check.names=TRUE)
    num_input <- nrow(data)

    ## checks that file contains required annotation/contaminant columns,
    ## the id column of MaxQuant output files is a unique key, if the values
    ## are non-numeric or duplicated the file may contain some rows with trash
    ## information that needs to be removed. if this is the case a msg is returned
    ## alerting the user that the file needs manual inspection. A data.frame
    ## of the input data is returned if it contains all the required columns
    ## and the id column contains numeric unique values (0-n)
    if(all(reqCols %in% colnames(data))==TRUE){
      if(all(!duplicated(data$id) & is.integer(data$id))){

        cli::cli_inform("{pipe} {enrich} file {.file {file}} imported")
        cli::cli_inform("Input file contains {num_input} {enrich} entries")

        data2 <- list(data=data, num_input=num_input)
        return(data2)

      } else {
        stop(message(paste("Error! 'id' column values. The id' column of",
                           "maxQuant output files should contain unique",
                           "integer values (e.g. 0-n). Inspect the file manually.",
                           "e.g. sort id column a -> z. Remove junk information",
                           "typically at the end of the file.")))
      }
    } else {
      stop(message(paste("Error! the input file is missing one or more required",
                         "annotation/contamination columns. For this pipeline",
                         "the file should include the following column names: ",
                         paste(reqCols, collapse=", "))))
    }

  } ## TMT IMPORT


}


#' Remove protein contaminant rows
#'
#' The second subfunction called by \code{\link{extract_data}}. Removes
#' contaminant rows from the data.
#'
#' @param data Raw Maxquant intensity data, from which contaminants need to be
#'   removed. Within the \code{\link{extract_data}} function, this is "diaData",
#'   which is the data object from the \code{\link{import_data}} function that has
#'   gone through additional processing within the body of the
#'   \code{\link{extract_data}} function
#' @inheritParams extract_data
#'
#' @return A list of data with contaminants removed, with three elements:
#'   \enumerate{
#'     \item "data", a data frame with the contaminants removed. Rows are
#'       proteins, columns include both protein info and individual intensities.
#'     \item "num_contam", an integer giving the number of contaminant rows removed
#'     \item "num_qfilter", an integer giving the number of retained rows
#'       (ie, nrow() of the "data" slot).
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet
#'

remove_contaminants <- function(data,
                                pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                                enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)


  ## vector of required contaminant column names based on
  ## analysis pipeline and type of enrichment. these columns
  ## are processed by other functions that remove contaminants, etc.
  ## values for each variable are defined according to R syntax
  ## rules. Values for each variable are listed at the top of the
  ## functions.R file.
  if(pipe=="TMT" & enrich=="protein"){ contamColums <- tmtContamColums }
  if(pipe=="phosphoTMT" & enrich=="protein"){ contamColums <- proteinContamColums }
  if(pipe=="phosphoTMT" & enrich=="phospho"){ contamColums <- phosphoContamColums }
  if(pipe=="DIA" & enrich=="protein"){ contamColums <- diaContamColums }
  if(pipe=="LF" & enrich=="protein"){ contamColums <- tmtContamColums }

  if(pipe=="DIA"){

    num_rows <- nrow(data)
    data <- data[!grepl("DECOY", data[, contamColums]),]
    data <- data[!grepl("Group of", data[, contamColums]),]

    num_contam  <- num_rows - nrow(data)
    num_qfilter <- nrow(data)

    cli::cli_inform("{num_contam} contaminantes removed")
    cli::cli_inform("{num_qfilter} {pipe} {enrich} entries retained")

    data2 <- list(data=data, num_contam=num_contam, num_qfilter=num_qfilter)
    return(data2)

  } ## DIA CONTAM

  if(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF"){

    num_rows <- nrow(data)
    for(x in contamColums){
      data[,x][is.na(data[,x])] <- ""
      remove<-data[,x]=="+"
      data <- data[remove==FALSE, ]
    }

    num_contam  <- num_rows - nrow(data)
    num_qfilter <- nrow(data)

    print(paste(num_contam, "contaminants removed. Success!!"))
    print(paste(num_qfilter, pipe, enrich, "entries retained ..."))

    data2 <- list(data=data, num_contam=num_contam, num_qfilter=num_qfilter)
    return(data2)

  } ## TMT/phosphoTMT/LF

}


#' Extract protein data
#'
#' Third subfunction called by \code{\link{extract_data}}. Processes the raw
#' protein annotation strings to extract standard gene names,
#' accession numbers, and Ids.
#'
#' @param data Data from which to extract protein data. In the
#'   \code{\link{extract_data}} function, this is \code{qfilterData$data}, which
#'   is the "data" slot of the result of the \code{\link{remove_contaminants}}
#'   output: the data frame of both protein and intensity data.
#' @inheritParams extract_data
#' @return A list. Items in list vary depending on the pipeline used. For
#'   DIA/protein, a list with 9 slots:
#'   \enumerate{
#'     \item "rawData": A dataframe, where each row is a proetin and each column
#'       is a sample. Not actually raw, I think(?): at this point in the pipeline,
#'       contaminants have been removed. Rownames in this data frame give a protein ID,
#'       in the form of UniprotID_Gene_name_id.
#'     \item "rawAnnot"-A dataframe, where each row is a protein and the columns
#'       are further protein annotation information. Rownames in this data frame give
#'       a protein ID, in the form of UniprotID_Gene_name_id.
#'     \item "sampleIDs"- A vector of the sample IDs for the data. Seems like the
#'       column names from the rawData slot. Redundant?
#'     \item "annotColumns"- The column names of the "rawAnnot" slot. Redundant?
#'     \item "num_samples"- The number of samples.
#'       Currently seems redundant: just ncol() of rawData slot. Maybe useful if we
#'       moved back to merging annotation and sample data?
#'     \item "num_extracted"- The number of proteins (nrow() of rawData slot).
#'       Redundant?
#'     \item "pipe"- The pipeline argument used.
#'     \item "enrich"- The enrich argument used.
#'     \item "pattern"- Unclear what this is at the moment. Comes from the
#'       results of the \code{\link{extract_sampleIDs}} function. For DIA data, though,
#'       there is no "pattern" slot in the result list from
#'       \code{\link{extract_sampleIDs}}. For TMT and other data, seems to be the the
#'       regex pattern that was used for getting sample IDs from columns?
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet


extract_protein_data <- function(data,
                                 sampleIDs=NULL,
                                 pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                                 enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)

  ## extract protein data and annotation from DIA experiment
  if(all(pipe=="DIA" & enrich=="protein")){ ## DIA/PROTEIN

    ## get sampleIDs defined by user
    if(!is.null(sampleIDs)){
      IDs <- extract_sampleIDs(colNames=colnames(data),sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
      pattern   <- IDs$pattern
      sampleIDs <- IDs$sampleIDs
    }

    ## if sampleIDs are not defined by user then use all column names
    ## of the input data in the search.
    if(is.null(sampleIDs)){
      IDs <- extract_sampleIDs(colNames=colnames(data), sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
      pattern   <- IDs$pattern
      sampleIDs <- IDs$sampleIDs
    }

    ## extract gene name, gene symbol, and uniprot id info. from the Fasta.header
    ## append to qfilterData add unique protein ID (proKey = uniprot_GN_id
    data$Accession.Number <- gsub("(.+?)(\\ .*)","\\1", stringr::str_extract(data$Protein.Name,"(?<=\\|)[^\\|]+(?=\\ )"))
    data$UniprotID    <- stringr::str_extract(data$Protein.Name, "(?<=\\|)[^\\|]+(?=\\|)")
    data$Gene_name    <- stringr::str_extract(data$Protein.Name, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
    data$Description  <- stringr::str_extract(data$Protein.Name, "(?<= )[^\\|]+(?= OS\\=)")

    if(pipe=="DIA"){ annotColums <- diaAnnotationColums }
    annotColums <- c(annotColums, "UniprotID", "Gene_name","Description")
    rownames(data) <- paste(data$UniprotID, data$Gene_name, data$id, sep="_")

    ## extract annotation columns and protein intensity data
    rawAnnot <- data[, annotColums]
    rawData  <- data[, sampleIDs]

    ## intensity data is not a numeric data.frame and requires additional processing
    ## numbers are character strings separated by commas (4,000,000), some cells
    ## contain text 'Missing Value' if the protein does not have detectable expression.
    ## replace 'Missing Value' with zeros, remove comma's, convert data in each column
    ## to numeric values
    rawData[,][rawData[,]=="Missing Value"] <- 0
    for(i in 1:ncol(rawData)){ rawData[,i] <- remove_commas(rawData[,i]) }

    num_samples   <- ncol(rawData)
    num_extracted <- nrow(rawData)

    cli::cli_inform(c("Intensity data for {num_extracted} {pipe} {enrich} entries and {num_samples} samples extracted"))
    cli::cli_rule()
    cli::cli_inform(c("v" = "Success!!"))
    # print(paste("intensity data for", num_extracted,
    #             pipe, enrich, "entries and ",num_samples, "samples extracted.",
    #             "Success!!"));cat("\n")

    data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                  annotColums=annotColums, num_samples=num_samples,
                  num_extracted=num_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
    return(data2)


  } ## DIA/PROTEIN


  ## after filtering contaminants append additional protein annotation info.,
  ## extract protein data, parse into annotation and sample protein intensities.
  if(all((pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF") & enrich=="protein")){ ## TMT/PROTEIN

    ## get sampleIDs if defined by user.
    if(!is.null(sampleIDs)){
      IDs <- extract_sampleIDs(colNames=colnames(data),sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)
      pattern   <- IDs$pattern
      sampleIDs <- IDs$sampleIDs
    }

    ## if sampleIDs are not defined by user then use all column names
    ## of the input data in the search.
    if(is.null(sampleIDs)){
      IDs <- extract_sampleIDs(colNames=colnames(data), sampleIDs=sampleIDs,pipe=pipe, enrich=enrich)
      pattern   <- IDs$pattern
      sampleIDs <- IDs$sampleIDs
    }

    if(pipe=="TMT" | pipe=="LF"){ annotColums <- tmtAnnotationColums }
    if(pipe=="phosphoTMT"){ annotColums <- proteinAnnotationColums }

    ## extract gene name, gene symbol, and uniprot id info. from the Fasta.header
    ## append to qfilterData add unique protein ID (proKey = uniprot_GN_id
    data$UniprotID   <- stringr::str_extract(data$Fasta.headers, "(?<=\\|)[^\\|]+(?=\\|)")
    data$Gene_name   <- stringr::str_extract(data$Fasta.headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
    data$Description <- stringr::str_extract(data$Fasta.headers, "(?<= )[^\\|]+(?= OS\\=)")
    annotColums <- c(annotColums, "UniprotID", "Gene_name","Description")

    ## set row names as a combination of uniprotid_gene.name_id (protein key)
    rownames(data) <- paste(data$UniprotID, data$Gene_name, data$id, sep="_")

    ## extract annotation columns and corrected protein intensity data
    ## intensity corrected columns) for the protein lysate samples.
    rawAnnot <- data[, annotColums]
    rawData  <- data[, sampleIDs]

    num_samples   <- ncol(rawData)
    num_extracted <- nrow(rawData)

    print(paste("corrected reporter intensity data for", num_extracted,
                pipe, enrich, "entries and ",num_samples, "samples extracted.",
                "Success!!"));cat("\n")

    data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                  annotColums=annotColums, num_samples=num_samples,
                  num_extracted=num_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
    return(data2)

  } ## TMT/PROTEIN

}


#' Make a log file
#'
#' Fourth subfunction called by \code{\link{extract_data}}. Used to write out
#' some logs during the data extraction process. Might want to rethink how this works.
#'
#' @param param The param object within the \code{\link{extract_data}} function
#'   body. Slightly unclear what the requirements for this are.
#' @param stats The stats object within the \code{\link{extract_data}} function
#'   body. Slightly unclear what the requirements for this are.
#' @param title A title for the printed parameter block. Default is "".
#' @param save TRUE/FALSE: should the parameters be save dot a text file?
#'
#' @return A list, with the param and stats objects. If SAVE = TRUE, has side
#'   effect of creating log files.
#'
#' @export
#'
#' @examples
#' # No examples yet

make_log <- function(param, stats, title="", save=TRUE){

  param<-data.frame(t(as.data.frame(t(param))))
  colnames(param)<-""

  stats<-data.frame(t(as.data.frame(t(stats))))
  colnames(stats)<-""

  if(save==TRUE){

    if(!dir.exists("logs")){dir.create("logs",recursive=TRUE)}

    if(all(length(param) > 0 & nrow(param)>0)){

      if(!file.exists(file.path("logs","parameters.log"))){
        sink(file="./logs/parameters.log",append=FALSE)
      } else { sink(file="./logs/parameters.log",append=TRUE) }

      # sink(file="./param.txt",append=TRUE)
      cat(paste0("\n##",paste(rep("-",40),collapse="")))
      cat(paste0("\n##  ",title,"\n"))
      cat(paste0("##",paste(rep("-",40),collapse="")))
      print(param); cat("\n\n")
      sink()
    }
    colnames(param)<-"Value"


    if(all(length(stats) > 0 & nrow(stats) > 0)){

      if(!file.exists(file.path("logs","processing.log"))){
        sink(file="./logs/processing.log",append=FALSE)
      } else { sink(file="./logs/processing.log",append=TRUE) }

      # sink(file="./log.txt",append=TRUE)
      cat(paste0("\n##",paste(rep("-",40),collapse="")))
      cat(paste0("\n##  ",title,"\n"))
      cat(paste0("##",paste(rep("-",40),collapse="")))
      print(stats); cat("\n\n")
      sink()
    }
    colnames(stats)<-"Value"
  } ## SAVE

  data2<-list(param=param, stats=stats)
  return(data2)

}

###### NEED TO DOCUMENT AT SOME POINT
###### NOT USED FOR DIA, THOUGH
local_prob_filter <- function(data,
                              min.prob,
                              pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                              enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)

  num_rows  <- nrow(data)
  min.prob <- ifelse(is.numeric(min.prob),min.prob,0.75)
  remove   <- data[, phosphoLocalProbColum] < min.prob
  data     <- data[remove==FALSE, ]

  num_localprob_removed <- num_rows - nrow(data)
  num_localprob_kept    <- nrow(data)

  print(paste(num_localprob_removed, pipe, "entries with localization",
              "probabilities < ", min.prob, "removed. Success!!"))
  print(paste(num_localprob_kept, pipe, "entries kept for further analysis..."))

  data2 <- list(data=data, min.prob=min.prob, num_localprob_removed=num_localprob_removed,
                num_localprob_kept=num_localprob_kept)
  return(data2)


}




#' Extract sample IDs
#'
#' A subfunction called by \code{\link{extract_protein_data}}. For DIA, Sort of does two
#' different things. In the "default" case, with sampleIDs=NULL, it takes in a
#' vector of column names (in our pipeline, these are the column names from
#' \code{qfilterData$data}, which is the "data" slot of the result
#' of the \code{\link{remove_contaminants}} output and is a data frame
#' of both protein and intensity data). Then, it extracts sample names from those
#' column names. If a vector of sampleIDs are supplied, it checks to see if those
#' IDs are present in the column names. If so, it returns them. If not, gives an
#' error.
#'
#' @param colNames A list of column names.
#' @inheritParams extract_data
#' @return A list, with a single slot. I think a vector of sample IDs.
#'
#' @export
#'
#' @examples
#' # No examples yet

extract_sampleIDs <- function(colNames,
                              sampleIDs=NULL,
                              pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                              enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)

  ## the first 11 columns of DIA sample reports contain
  ## 2 col. of junk & 9 annotation info. sample data should
  ## begin in column 12, with pool samples used to create the
  ## library in the first 3 columns (12-14) followed by the
  ## experimental samples.
  ## NOTE: the 2 columns of the file are removed when the file
  ## is imported. Thus samples will start in column 10.
  if (pipe == "DIA" & enrich == "protein") {

    if (!is.null(sampleIDs)) {

      # Maybe re-write this section so we don't have to do a separate check for the case when 0 are found?
      keep<-colNames%in%sampleIDs
      sampleIDs_out<-colNames[keep==TRUE]

      # If None of the sampleIDs are found
      if (identical(sampleIDs_out, character(0))) {
        cli::cli_abort(c("Provided sampleId{?s} not found in column names of input file {cli::qty(sampleIDs)}",
                         "x" = "missing sampleID{?s}: {.arg {sampleIDs}}"))
      }

      # If only some are found
      if (!all(sampleIDs %in% sampleIDs_out)) {
        missing <- sampleIDs[sampleIDs %notin% sampleIDs_out]

        cli::cli_abort(c("Provided sampleID{?s} not found in column names of input file {cli::qty(missing)}",
                         "x" = "missing sampleID{?s}: {.arg {missing}}"))
      }

    } ## NOT NULL

    if(is.null(sampleIDs)){
      ## columns containing sampleIDs should end with .mzML extension.
      keep <- grep("mzML", colNames)
      sampleIDs_out <- colNames[keep]

      ## if any annotation/contaminant column names are still included in list
      ## sampleIDs remove them.
      remove    <- sampleIDs_out %in% unique(diaAnnotationColums,diaContamColums)
      sampleIDs_out <- sampleIDs_out[remove==FALSE]

    } ## NULL


    data2 <- list(sampleIDs=sampleIDs_out)
    return(data2)

  } ## DIA

  ## label free looking for iBAQ.XXX columns. user should input ids
  if(pipe=="LF" & enrich=="protein"){

    if(!is.null(sampleIDs)){

      IDs<-sampleIDs
      pattern="NULL"
      keep<-colNames%in%IDs
      sampleIDs<-colNames[keep==TRUE]
      stopifnot(length(sampleIDs)==length(IDs))
      if(identical(sampleIDs, character(0))){
        stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
      }
    } ## NOT NULL

    if(is.null(sampleIDs)){

      pattern="iBAQ."
      sampleIDs <- grep(pattern, colNames, value=TRUE)
      sampleIDs <- sampleIDs[-grep("iBAQ.peptides",sampleIDs)]

      if(identical(sampleIDs, character(0))){
        stop(paste0("\nsampleIDs could not be identified based on pattern='",pattern,
                    "' sampleIDs==character(0)"))
      }

    } ## NULL

    data2<-list(sampleIDs=sampleIDs, pattern=pattern)
    return(data2)

  } ## LF

  ## TMT/phosphoTMT
  if(pipe=="TMT" | pipe=="phosphoTMT"){

    if(enrich=="protein"){

      if(!is.null(sampleIDs)){

        pattern="NULL"
        IDs<-sampleIDs

        keep<-colNames%in%IDs
        sampleIDs<-colNames[keep==TRUE]
        stopifnot(length(sampleIDs)==length(IDs))
        if(identical(sampleIDs, character(0))){
          stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
        }
      } ## NOT NULL


      ## MaxQuant TMT sample data columns of interest. X is the index
      ## of the TMT tag starting with 0.
      ## https://adinasarapu.github.io/posts/2020/01/blog-post-tmt/
      ## proteinGroups.txt		   "Reporter intensity corrected X <experiment name(s)>"
      ## peptides.txt			   "Reporter intensity corrected X <experiment name(s)>"
      ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>"

      ## get the corrected reporter intensity column names for each protein
      ## lysate sample
      ## "Reporter.intensity.corrected.X.<experiment name(s)> _Lysate"

      if(is.null(sampleIDs)){
        pattern <- "Reporter\\.intensity\\.corrected\\.[[:digit:]]+.*[L/l]ysate*"
        sampleIDs <- grep(pattern, colNames, value=TRUE)

        ## if pattern doesn't match any column names then reduce pattern to more
        ## generic format.
        ## "Reporter.intensity.corrected.X"
        if(identical(sampleIDs, character(0))){
          pattern = "Reporter\\.intensity\\.corrected\\.[[:digit:]]+"
          sampleIDs = grep(pattern, colNames, value=TRUE)
          cat("sampleIDs identified using generic pattern.")
        }

      }## NULL

      data2 <- list(sampleIDs=sampleIDs, pattern=pattern)
      return(data2)

    } ## PROTEIN

    if(enrich=="phospho"){

      if(!is.null(sampleIDs)){

        pattern="NULL"
        IDs<-sampleIDs
        keep<-colNames%in%IDs

        classIDs  <- colNames[keep==TRUE]
        stopifnot(length(classIDs)==length(IDs))
        classNums <- gsub("^.*_{3}","",classIDs) ## keep everything after ____ (class number)
        sampleIDs <- gsub("_{3}.*","",classIDs)  ## keep everything before ____ (should match targets)

        if(identical(sampleIDs, character(0))){
          stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
        }
      } ## NOT NULL


      ## MaxQuant TMT sample data columns of interest. X is the index
      ## of the 'tag', starting with 0. ___XX indicates peptides with
      ## single, double, triple phospho sites)
      ## better to analyze ___1, ___2, ___3, than summarized site
      ## intensities because in a biological system the same protein
      ## may have distinct functions when it is differentially phosphorylated.
      ## https://adinasarapu.github.io/posts/2020/01/blog-post-tmt/
      ## proteinGroups.txt		   "Reporter intensity corrected X <experiment name(s)>"
      ## peptides.txt			   "Reporter intensity corrected X <experiment name(s)>"
      ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>"
      ## Phospho (STY)Sites.txt	"Reporter intensity corrected X <experiment name(s)>___XX"

      ## get column names for the corrected reporter intensities for each
      ## phospho class
      ## "Reporter.intensity.corrected.X.<experiment name(s)> _Phospho___1"
      if(is.null(sampleIDs)){

        pattern   <- "Reporter\\.intensity\\.corrected\\.[[:digit:]]+.*[P/p]hos*"
        classIDs  <- grep(pattern, colNames, value = TRUE)
        classNums <- gsub("^.*_{3}","",classIDs) ## keep everything after ____ (class number)
        sampleIDs <- gsub("_{3}.*","",classIDs)  ## keep everything before ____ (should match targets)

        ## if pattern doesn't match any column names then reduce pattern
        ## to more generic format.
        ## "Reporter.intensity.corrected.X"
        if(identical(classIDs, character(0))){
          pattern   <- "Reporter\\.intensity\\.corrected\\.*[[:digit:]].*\\_{3}[[:digit:]]$"
          classIDs  <- grep(pattern, colNames, value = TRUE)
          classNums <- gsub("^.*_{3}","", classIDs) ## keep everything after ____
          sampleIDs <- unique(gsub("_{3}.*","",classIDs)) ## keep everything before ____
          cat("classIDs identified using generic pattern.")
        }


      } ## NULL


      data2 <- list(classIDs=classIDs, sampleIDs=sampleIDs, classNums=classNums,
                    pattern=pattern)
      return(data2)

    } ## PHOSPHO

  } ## TMT/PHOSPHOTMT


}


extract_phospho_data <- function(data, sampleIDs=NULL, pipe, enrich){

  if(pipe=="phosphoTMT" & enrich=="phospho"){

    annotColums <- phosphoAnnotationColums

    ## after filtering contaminants & low prob sites parse the Fasta.header
    ## column to obtain gene name (Description), gene symbol (Gene_name),
    ## UniprotID, and Flanking aa sequence (7aa-site-7aa), extract first id
    ## from protein.group.IDs, first a.a. position within protein from
    ## Positions.within.proteins, and create extract gene name,
    ## gene symbol, uniprot id info. from phospho(STY)sites Fasta_header and
    ## Flanking a.a.sequence around phospho sites from phospho(Group)STY)sites
    ## Sequence_window and add to the qfilterData data.frame
    data$UniprotID   <- stringr::str_extract(data$Fasta.headers, "(?<=\\|)[^\\|]+(?=\\|)")
    data$Gene_name   <- stringr::str_extract(data$Fasta.headers, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
    data$Description <- stringr::str_extract(data$Fasta.headers, "(?<= )[^\\|]+(?= OS\\=)")
    data$Flanking    <- gsub("\\;.*$", "", data$Sequence.window) %>% stringr::str_sub(9,23) %>% paste0("-p")

    ## get protein group ID (id col. proteinGroups) and site position for
    ## the first protein entry i.e. corresponds to Majority.protein.ID
    ## to get the first entry remove everything after the first semicolon.
    data$proGroupID     <- sub(";.*", "", data$Protein.group.IDs)
    data$phosPosition   <- sub(";.*", "", data$Positions.within.proteins)
    data$phosAAPosition <- paste(data$Amino.acid, data$phosPosition,sep="")

    annotColums <- c(annotColums, "proGroupID", "UniprotID", "Gene_name",
                     "Description", "Flanking", "phosPosition", "phosAAPosition")
    rownames(data) <- paste(data$UniprotID, data$Gene_name, data$phosAAPosition,sep="_")

    ## GET SAMPLE IDS
    ## use user defined sampleIDs to extract the corrected intensity data for
    ## single, double, and tripple phosphopeptides ___1, ___2, ____3.
    if(!is.null(sampleIDs)){
      table(sampleIDs %in% colnames(data))
      IDs <- extract_sampleIDs(colNames=sampleIDs, pipe=pipe, enrich=enrich)
    }

    ## use extract_sampleIDs() function to extract the corrected intensity data for
    ## single, double, and tripple phosphopeptides ___1, ___2, ____3.
    if(is.null(sampleIDs)){ IDs <- extract_sampleIDs(colNames=colnames(data),
                                                     pipe=pipe, enrich=enrich) }
    classIDs  <- IDs$classIDs; classIDs
    classNums <- IDs$classNums
    pattern   <- IDs$pattern

    ## extract annotation and corrected intensity data
    ## for single, double, tripple phosphopeptides
    rawAnnot <- data[, annotColums]; dim(rawAnnot)
    rawData  <- data[, classIDs]; dim(rawData)

    num_samples   <- ncol(rawData)
    num_extracted <- nrow(rawData)

    cli::cli_inform(c("Intensity data for {num_extracted} {pipe} {enrich} entries and {num_samples} samples extracted."))

    ## EXTRACT CLASS DATA
    classList <- vector("list", max(IDs$classNums))
    names(classList) <- c(paste("class_",1:max(IDs$classNums),sep=""))
    classData <- NULL; classAnnot <- NULL

    for(i in 1:length(classList)){ ## INDIV CLASSES
      ## get class == i
      keep      <- IDs$classNums==i; table(keep)
      classIDs  <- IDs$classIDs[keep==TRUE]   ## sampleID___1
      sampleIDs <- IDs$sampleIDs[keep==TRUE]  ## sampleID
      classNums <- IDs$classNums[keep]        ## 1
      stopifnot(all(classNums==i))

      ## INDIVIDUAL CLASS DATA
      classData <- rawData[, classIDs]  ##__1
      stopifnot(colnames(classData)==classIDs)
      colnames(classData) <- sampleIDs  ## rename minus __1

      ## INDIVIDUAL CLASS ANNOTATION
      classAnnot <- rawAnnot[, annotColums]
      classAnnot$Class <- rep(i, nrow(classAnnot))  ## Class
      annotColums2 <- c(annotColums, "Class")
      stopifnot(rownames(classData)==rownames(classAnnot))

      num_class_samples <- ncol(classData)
      num_class_phospho <- nrow(classData)

      print(paste("classData",i,":","intensity data for",num_class_phospho,
                  "class", i, pipe,enrich, "entries and", num_class_samples,
                  "samples extracted. Success!!",sep=" "))

      ## id list from extract_sampleIDs() for class i
      idList <- list(classIDs=classIDs, sampleIDs=sampleIDs, classNums=classNums,
                     annotColums=annotColums2, pattern=pattern)
      classList[[i]] <- list(rawData=classData, rawAnnot=classAnnot[rownames(classData), ],
                             idList=idList, num_class_samples=num_class_samples,
                             num_class_phospho=num_class_phospho)
    }## INDIV CLASSES

    data2 <- list(classList=classList)
    return(data2)

  }

}

