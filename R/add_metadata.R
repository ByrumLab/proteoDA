
#' Add metadata to a DAList
#'
#' Add a dataframe of metadata to a DAList object. The metadata file defines the
#' sample labels, groups, and other factors (such as batch, gender, age, paired samples, etc.)
#' required for the analysis and design matrix.
#'
#' @param DAList A DAList to which metadata will be added.
#' @param metadata_file A metadata file defining sample information necessary for the statistical design.
#'
#' @return An S3 list type with added metadata.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data <- add_metadata(data, "path/to/metadata_file.csv")
#'
#' # An example metadata file contains columns: number, sample, batch, group. The number column
#' in the metadata must match the sample name in the data file. For instance, metadata file number "1"
#' corresponds to data column "Sample_01". This matching defines the other factors about the sample.
#' }
#'
add_metadata <- function(DAList,
                         metadata_file) {



  validate_DAList(DAList)

  if (!is.null(DAList$metadata)) {
    cli::cli_abort(c("!" = "input DAList already contains metadata"))

  }

  ## extract sample number from sample_IDs and use as data key
  number  <- get_dia_sample_number(sample_IDs = colnames(DAList$data))

  # Import the metadata
  # which has some formatting rules (see documentation)
  metadata <- import_meta(input_file = metadata_file)

  # Then, cross-check the sample numbers read in from the sample_IDs
  # against the sample numbers in the metadata
  # If they don't match, give some warnings and output non-integrated lists
  if ((length(number) > length(metadata$number)) | (any(number %notin% metadata$number))) {
    # more samples in data than in metadata

    missing_samples <- number[number %notin% metadata$number]

    cli::cli_abort(c("!" = "Not all samples supplied in the data are present in metadata file",
                    "!" = "sample number{?s} {missing_samples} not found in metadata file"))
  }

  if ((length(metadata$number) > length(number)) | (any(metadata$number %notin% number))){
    # more samples in metadata than in data columns

    missing_samples <- metadata$number[metadata$number %notin% number]


    cli::cli_abort(c("!" = "Not all samples in metadata file are present in the data",
                    "!" = "sample number{?s} {missing_samples} in metadata file not found in data"))
  }

  ## create basic targets info.
  targets <- data.frame(sampleIDs = colnames(DAList$data),
                        number = number,
                        row.names = colnames(DAList$data))

  ## use sample number info. to sort targets info. so that it is in the same order
  ## as the meta data file. combine targets and metadata into a single data.frame
  m <- match(metadata$number, number)
  targets <- targets[m, ]
  remove  <- colnames(metadata) %in% colnames(targets)
  targets <- cbind(targets, metadata[,remove==FALSE])


  ## change data column names and targets row names
  ## to sample name i.e. sample column in targets
  #TODO NEED TO FIGURE THIS OUT
  # NEEDED FOR OUR CURRENT DATA, SHOULD rethink for public
  if ("sample" %in% colnames(targets)) {
    rownames(targets) <- targets$sample
    colnames(DAList$data) <- rownames(targets)

    cli::cli_inform(cli::col_yellow("\"sample\" column found in the input metadata file: renamed column names of {.arg data} and row names of {.arg metadata} to sample names"))
  }

  DAList$metadata <- targets
  validate_DAList(DAList)
}


#' Import metadata file
#'
#' Imports a delimited metadata file, checking for some required columns.
#'
#' @param input_file A metadata file, giving info on the samples.
#'
#' @return A dataframe of the imported metadata.
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
import_meta <-function(input_file) {

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

  if (filext == "txt" | filext=="tsv") {sep= "\t"}
  if (filext == "csv") {sep=","}

  ## the metadata file should contain the following required columns
  reqCols <- diaMetaColums

  ## import file
  data <- utils::read.csv(file = input_file,
                          sep = sep,
                          stringsAsFactors = FALSE,
                          header = TRUE,
                          check.names = TRUE,
                          strip.white = T)


  if (!all(reqCols %in% colnames(data))) {
    missing_cols <- reqCols[reqCols %notin% colnames(data)]

    cli::cli_abort(c("Problem with columns in metadata file.",
                     "x" = "Required column(s) not present in {.file {file}}",
                     "x" = "Missing column{?s}: {missing_cols}"))

  } else {

    cli::cli_inform("metadata file successfully imported")
    return(data)

  }
}


#' Get DIA sample number from sample IDs
#'
#' Extracts the sample number from a vector of sample IDs
#'
#' @param sample_IDs A character vector of sample IDs, from which the targets
#'   file is constructed. Generally, these are the column names of the
#'   data slot in the object returned by \code{\link{read_DIA_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sample_IDs have a required format:
#'   "Pool" samples must include the regex pattern \code{[P/p]ool_[[:digit:]]} in
#'   the column name, and samples must contain the regex pattern
#'   \code{[S/s]ample_[[:digit:]]} in the column name.
#'
#' @return A numeric vector of the extracted sample numbers
#' @keywords internal
#'
#' @examples
#' # No examples yet
get_dia_sample_number <- function(sample_IDs) {

  ## scaffold DIA sample reports should contain
  ## Pool_1/Sample_1 as part of the sample_IDs
  ## e.g. "Balachandran_061721_Pool_3.mzML",
  ## e.g. "Balachandran_061721_Sample_01_CS1.mzML"
  pattern1 <- "[P/p]ool_[[:digit:]]"   ## POOL PATTERN (Pool_1)
  pattern2 <- "[S/s]ample_[[:digit:]]" ## SAMPLE PATTERN (Sample_1)
  pat1     <- grep(pattern1, sample_IDs)
  pat2     <- grep(pattern2, sample_IDs)
  pattern  <- c(pat1, pat2)

  if(!all(length(pattern) == length(sample_IDs))) { # Double-check this check

    supplied_ids <- cli::cli_vec(sample_IDs, list(vec_trunc = 5))
    cli::cli_abort(c("sample numbers could not be extracted from sample_IDs",
                     "x" = "Pool sample(s) must contain a string matching {.code {pattern1}}",
                     "!" = "For example: \"Pool_1\" or \"pool_01\"",
                     "x" = "Other sample(s) must contain a string matching {.code {pattern2}}",
                     "!" = "For example: \"Sample_1\" or \"sample_01\"",
                     "!" = "Supplied sample ids: {.val {supplied_ids}}"))
  } else {
    ## pool samples
    if (length(pat1) != 0) {
      number1 <- gsub("\\.", "_", sample_IDs[pat1])
      number1 <- gsub(".*[P/p]ool_", "", number1)  ## removes everything b4 Pool_
      number1 <- gsub("_.*", "", number1)  ## removes everything after number _mzML
      number1 <- paste0("P", as.numeric(number1)) ## adds capital P to number == P1, P2, P3
    } else {
      number1 <- NULL
    }

    ## samples
    number2 <- gsub("\\.","_", sample_IDs[pat2])
    number2 <- gsub(".*[S/s]ample_", "", number2) ## removes everything b4 Sample_
    number2 <- gsub("_.*", "", number2)   ## removes everthing after number _mzML
    number2 <- as.numeric(number2)      ## if 01,02 == 1,2
    number  <- c(number1, number2)       ## P1,P2,P3, 1,2,3,..., n

    cli::cli_inform("sample numbers extracted from DIA sample_IDs.")
    return(number)
  }
}


