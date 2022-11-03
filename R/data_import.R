
#' Loading Maxquant data
#'
#' A function to extract and load DIA MaxQuant data. The \code{read_DIA_data}
#' function is the main function, and it calls a number
#' of other subfunctions to process and prepare the data. See the links to
#' subfunctions for more info. Returns a list of data from maxquant, along with
#' some stats and other info. Also has some side effects: creates some log files.
#'
#'
#' Subfunctions (DIA): \code{\link{read_maxquant_delim}}, \code{\link{remove_contaminants}},
#' \code{\link{extract_protein_data}}, \code{\link{make_log}}.
#'
#' @param input_file The file of Maxquant data to be imported. Format and filetype
#'   vary depending on the type of data to be analyzed. If not supplied, R will
#'   ask you to pick a file from the file browser.
#' @param sample_IDs Optional: a vector of sample IDs to analyze. Default is NULL, in which
#'   case all sample IDs are analyzed. These are determined from column names.
#'
#' @return A list of extracted data, with 5 elements:
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
read_DIA_data <- function(input_file = NULL,
                         sample_IDs = NULL) {


  ## if no file is input by user the open file dialog box is
  ## called to allow the user to navigate to and select the file.
  if (is.null(input_file)) {
    input_file <- file.choose()
  }

  cli::cli_rule()

  ## input parameters are added to param variable
  param <- stats <- list()
  param[["file"]] <- input_file
  param[["pipe"]] <- "DIA"
  param[["enrich"]] <- "protein"
  param[["sample_IDs"]] <- ifelse(is.null(sample_IDs), "NULL", paste(sample_IDs, collapse = ", "))
  param[["min.prob"]] <- "NULL"

  ## IMPORT DATA
  ## the first 10 lines are skipped and last line removed. column 1 (.X)
  ## is renamed "id"
  DIADATA <- read_maxquant_delim(input_file = input_file)
  diaData <- DIADATA$data
  reqCols <- DIADATA$reqCols

  param[["input_file"]] <- DIADATA$file
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
  replaceNumbers <- colnames(diaData) %in% c(reqCols, diaAnnotationColums)
  replaceColums <- colnames(diaData)[replaceNumbers == FALSE]
  for (i in replaceColums) {
    diaData[, i] <- remove_commas(diaData[, i])
  }

  bn <- basename(input_file)
  bn <- gsub("Samples Report of ", "", bn)
  ilab <- gsub(paste0(".", file_extension(bn)), "", bn)
  param[["ilab"]] <- ilab


  ## REMOVE CONTAMINANTS
  qfilterData <- remove_contaminants(data = diaData)
  stats[["num_contam_removed"]] <- qfilterData$num_contam

  ## EXTRACT PROTEIN DATA
  extProt <- extract_protein_data(
    data = qfilterData$data, sample_IDs = sample_IDs
  )

  stats[["total_input_samples"]] <- extProt$num_samples
  stats[["num_extracted_rows"]] <- extProt$num_extracted

  ## SAVE PARAM/STATS TO LOG FILE
  logs <- make_log(param = param, stats = stats, title = "EXTRACTED DATA", save = TRUE)

  ## return a list of data.frame containing the extracted quality
  ## filtered intensity data, corresponding extracted annotation,
  ## stats, and input parameters
  list(
    data = extProt$rawData, annot = extProt$rawAnnot,
    param = logs$param, stats = logs$stats[c(3, 1, 2, 4), ], ilab = ilab
  )
}



#' Import MaxQuant data
#'
#' First subfunction called by \code{\link{read_DIA_data}}. Reads in the tabular
#' data and does some checks for required columns.
#'
#' @inheritParams read_DIA_data
#'
#' @return A list with 3 slots:
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

read_maxquant_delim <-function(input_file) {
  ## check that the file is a csv, tsv, or text file
  filext <- stringr::str_to_lower(file_extension(input_file))
  if (filext %notin% c("csv","txt","tsv")) {
    cli::cli_abort(c("Problem with input file",
                     "x" = "Input file must end in {.file .csv}, {.file .tsv}, or {.file .txt}"))
  }

  ## check that file exists
  if (!file.exists(input_file)) {
    cli::cli_abort(c("Cannot find file.",
                     "x" = "{.file {input_file}} does not exist",
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
  reqCols <- unique(c(diaContamColums, "Protein.Name", "Accession.Number"))


  ## sample report generated by scaffold DIA is imported as data.frame.
  ## Sample Reports generated by Scaffold DIA contain meta data in the
  ## first 10 lines, and text in the last row to mark the end of the file.
  ## Thus, when reading in the file the first 10 lines are skipped and the
  ## last row is removed. column name of column 1 (X. on import) is changed to id
  ## column 1 must be changed to "id" to allow extraction and other processing.
  ## Note: column names of the input file are defined according to R syntax rules.
  ## special characters, spaces, are converted to periods.
  data <- utils::read.csv(file = input_file, sep = sep, stringsAsFactors = FALSE,
                          header = TRUE, check.names = TRUE, skip = 10)

  if (data[nrow(data),1] == "END OF FILE") {
    data <- data[-nrow(data),]
  }

  ## name of column 1 is changed to 'id'
  if (colnames(data)[1] == "X.") {colnames(data)[1] <- "id"}

  reqCols <- c(reqCols,"id")
  data <- data.frame(data)
  num_input <- nrow(data)

  ## checks to make sure the input file contains the required annotation/contamination
  ## columns. This ensures that the file imported correctly and helps to identify
  ## sampleIDs columns and extract targets info. later.
  if (!all(reqCols %in% colnames(data))) {
    missing_cols <- reqCols[reqCols %notin% colnames(data)]
    cli::cli_abort(c("Problem with columns in input file.",
                     "x" = "Required column(s) not present in {.file {input_file}}",
                     "x" = "Missing column{?s}: {missing_cols}"))
  } else {
    cli::cli_inform("DIA protein file {.file {input_file}} imported")
    cli::cli_inform("Input file contains {num_input} protein entries")

    output <- list(data = data,
                  num_input = num_input,
                  reqCols = reqCols)
  }

  output
}

#' Remove protein contaminant rows
#'
#' The second subfunction called by \code{\link{read_DIA_data}}. Removes
#' contaminant rows from the data.
#'
#' @param data Raw Maxquant intensity data, from which contaminants need to be
#'   removed. Within the \code{\link{read_DIA_data}} function, this is "diaData",
#'   which is the data object from the \code{\link{read_maxquant_delim}} function that has
#'   gone through additional processing within the body of the
#'   \code{\link{read_DIA_data}} function
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

remove_contaminants <- function(data) {
  ## vector of required contaminant column names
  contamColums <- diaContamColums

  num_rows <- nrow(data)
  data <- data[!grepl("DECOY", data[, contamColums]),]
  data <- data[!grepl("Group of", data[, contamColums]),]

  num_contam  <- num_rows - nrow(data)
  num_qfilter <- nrow(data)

  cli::cli_inform("{num_contam} contaminantes removed")
  cli::cli_inform("{num_qfilter} DIA protein entries retained")

  list(data = data,
       num_contam = num_contam,
       num_qfilter = num_qfilter)
}


#' Extract protein data
#'
#' Third subfunction called by \code{\link{read_DIA_data}}. Processes the raw
#' protein annotation strings to extract standard gene names,
#' accession numbers, and Ids.
#'
#' @param data Data from which to extract protein data. In the
#'   \code{\link{read_DIA_data}} function, this is \code{qfilterData$data}, which
#'   is the "data" slot of the result of the \code{\link{remove_contaminants}}
#'   output: the data frame of both protein and intensity data.
#' @inheritParams read_DIA_data
#' @return A list with 9 slots:
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
#'     \item "pipe"- The pipeline used, "DIA".
#'     \item "enrich"- The enrichment type, "protein".
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
                                 sample_IDs = NULL) {

  ## extract protein data and annotation from DIA experiment
  ## get sampleIDs defined by user
  IDs <- extract_sampleIDs(colNames = colnames(data), sample_IDs = sample_IDs)
  pattern   <- IDs$pattern
  sampleIDs <- IDs$sample_IDs


  ## extract gene name, gene symbol, and uniprot id info. from the Fasta.header
  ## append to qfilterData add unique protein ID (proKey = uniprot_GN_id
  data$Accession.Number <- gsub("(.+?)(\\ .*)","\\1", stringr::str_extract(data$Protein.Name,"(?<=\\|)[^\\|]+(?=\\ )"))
  data$UniprotID    <- stringr::str_extract(data$Protein.Name, "(?<=\\|)[^\\|]+(?=\\|)")
  data$Gene_name    <- stringr::str_extract(data$Protein.Name, "(?<=GN\\=)[^\\|]+(?= PE\\=)")
  data$Description  <- stringr::str_extract(data$Protein.Name, "(?<= )[^\\|]+(?= OS\\=)")

  annotColums <- diaAnnotationColums
  annotColums <- c(annotColums, "UniprotID", "Gene_name","Description")
  rownames(data) <- paste(data$UniprotID, data$Gene_name, data$id, sep="_")

  ##instead of requiring the data to include all annotation columns in the list
  ## only the columns that are present in the annotation list are extracted)
  # rawAnnot <- data[, annotColums]
  rawAnnot <- data[, colnames(data) %in% annotColums]
  rawData  <- data[, sampleIDs]

  ## intensity data is not a numeric data.frame and requires additional processing
  ## numbers are character strings separated by commas (4,000,000), some cells
  ## contain text 'Missing Value' if the protein does not have detectable expression.
  ## replace 'Missing Value' with zeros, remove comma's, convert data in each column
  ## to numeric values
  rawData[,][rawData[,]=="Missing Value"] <- NA
  for(i in 1:ncol(rawData)){ rawData[,i] <- remove_commas(rawData[,i]) }


  num_samples   <- ncol(rawData)
  num_extracted <- nrow(rawData)

  cli::cli_inform(c("Intensity data for {num_extracted} DIA protein entries and {num_samples} samples extracted"))
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success!!"))

  # Return results list
  list(
    rawData = rawData,
    rawAnnot = rawAnnot,
    sampleIDs = sampleIDs,
    annotColums = annotColums,
    num_samples = num_samples,
    num_extracted = num_extracted,
    pipe = "DIA",
    enrich = "protein",
    pattern = pattern
    )
}

#' Extract sample IDs
#'
#' A subfunction called by \code{\link{extract_protein_data}}. Does two
#' different things. In the "default" case, with sample_IDs = NULL, it takes in a
#' vector of column names (in our pipeline, these are the column names from
#' \code{qfilterData$data}, which is the "data" slot of the result
#' of the \code{\link{remove_contaminants}} output and is a data frame
#' of both protein and intensity data). Then, it extracts sample names from those
#' column names. If a vector of sample_IDs are supplied, it checks to see if those
#' IDs are present in the column names. If so, it returns them. If not, gives an
#' error.
#'
#' @param colNames A list of column names.
#' @inheritParams read_DIA_data
#' @return A list, with a single slot. I think a vector of sample IDs.
#'
#' @export
#'
#' @examples
#' # No examples yet
extract_sampleIDs <- function(colNames,
                              sample_IDs = NULL) {

  ## the first 11 columns of DIA sample reports contain
  ## 2 col. of junk & 9 annotation info. sample data should
  ## begin in column 12, with pool samples used to create the
  ## library in the first 3 columns (12-14) followed by the
  ## experimental samples.
  ## NOTE: the 2 columns of the file are removed when the file
  ## is imported. Thus samples will start in column 10.
  if (!is.null(sample_IDs)) {

    # Maybe re-write this section so we don't have to do a separate check for the case when 0 are found?
    keep <- colNames %in% sample_IDs
    sampleIDs_out <- colNames[keep == TRUE]

    # If None of the sampleIDs are found
    if (identical(sampleIDs_out, character(0))) {
      cli::cli_abort(c("Provided sampleId{?s} not found in column names of input file {cli::qty(sample_IDs)}",
                       "x" = "missing sampleID{?s}: {.arg {sample_IDs}}"))
    }

    # If only some are found
    if (!all(sample_IDs %in% sampleIDs_out)) {
      missing <- sample_IDs[sample_IDs %notin% sampleIDs_out]

      cli::cli_abort(c("Provided sampleID{?s} not found in column names of input file {cli::qty(missing)}",
                       "x" = "missing sampleID{?s}: {.arg {missing}}"))
    }

  } else {

    ## columns containing sampleIDs should end with .mzML extension.
    keep <- grep("mzML", colNames)
    sampleIDs_out <- colNames[keep]

    ## if any annotation/contaminant column names are still included in list
    ## sampleIDs remove them.
    remove    <- sampleIDs_out %in% unique(diaAnnotationColums,diaContamColums)
    sampleIDs_out <- sampleIDs_out[remove == FALSE]
  }

  list(sample_IDs = sampleIDs_out)
}

