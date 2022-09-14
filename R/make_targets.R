
#' Make a targets file
#'
#' Creates a dataframe of "targets". Output depends on whether a metadata file
#' is supplied, and whether the metadata and targets can be matched successfully
#' (see Value, below). Calls two subfunctions: \code{\link{import_meta}} and
#' \code{\link{get_dia_sample_number}}.
#'
#' @param sample_IDs A character vector of sample IDs, from which the targets
#'   file is constructed. Generally, these are the column names of the
#'   data slot in the object returned by \code{\link{read_DIA_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sample_IDs have a required format:
#'   "Pool" samples must include the regex pattern \code{[P/p]ool_[[:digit:]]} in
#'   the column name, and samples must contain the regex pattern
#'   \code{[S/s]ample_[[:digit:]]} in the column name.
#' @param input_file Optional: a metadata file, giving info on the samples.
#'
#' @return It depends:
#'   \itemize{
#'     \item If a metadata file is supplied and could be successfully matched
#'       against the input column names, a full targets dataframe is returned.
#'     \item If a metadata file is supplied but matching against column names
#'       was unsuccessful, returns, with a warning, a list with two slots:
#'       (1) the dataframe of assembled targets, and (2) the dataframe of
#'       imported metadata. These must be manually combined before further
#'       analysis.
#'     \item When a metadata file is not supplied, returns, with a warning, an
#'       incomplete targets dataframe.
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet
#'
make_targets <- function(input_file = NULL,
                         sample_IDs) {

  cli::cli_rule()

  ## extract sample number from sample_IDs and use as data key
  number  <- get_dia_sample_number(sample_IDs = sample_IDs)


  # If the user didn't supply a metadata file,
  # build a targets dataframe just from the sample IDs provided
  if (is.null(input_file)) {
    targets <- data.frame(sampleIDs = sample_IDs,
                          number = number,
                          pipe = rep("DIA", length(sample_IDs)),
                          enrichment= rep("protein", length(sample_IDs)),
                          row.names=sample_IDs)

    cli::cli_rule()
    cli::cli_warn(c("!" = "metadata file not provided",
                    "!" = "returning protein targets template for a DIA experiment",
                    "!" = "import metadata for these samples and add it to the returned target template dataframe"))
    output <- targets

  } else { # Otherwise if metadata supplied

    # Import the metadata
    # which has some formatting rules (see documentation)
    metadata <- import_meta(input_file = input_file,
                            sample_IDs = sample_IDs)

    # Then, cross-check the sample numbers read in from the sample_IDs
    # against the sample numbers in the metadata
    # If they don't match, give some warnings and output non-integrated lists
    if ((length(number) > length(metadata$number)) | (any(number %notin% metadata$number))) {
      # more samples in data than in metadata

      missing_samples <- number[number %notin% metadata$number]

      cli::cli_rule()
      cli::cli_warn(c("!" = "Not all samples supplied in {.arg samples_ID} are present in metadata file",
                      "!" = "sample number{?s} {missing_samples} not found in metadata file",
                      "!" = "returning basic target info and metadata as separate dataframes within a list",
                      "!" = "Combine them manually"))
      ## create basic targets info.
      targets <- data.frame(sampleIDs = sample_IDs,
                            number = number,
                            pipe = rep("DIA", length(sample_IDs)),
                            enrichment = rep("protein", length(sample_IDs)),
                            row.names = sample_IDs)

      output <- list(targets = targets,
                     metadata = metadata)
    } else if ((length(metadata$number) > length(number)) | (any(metadata$number %notin% number))){
      # more samples in metadata than in data columns

      missing_samples <- metadata$number[metadata$number %notin% number]

      cli::cli_rule()
      cli::cli_warn(c("!" = "Not all samples in metadata file are present in {.arg sample_IDs}",
                      "!" = "sample number{?s} {missing_samples} in metadata file not found in {.arg sample_IDs}",
                      "!" = "returning basic target info and metadata as separate dataframes within a list",
                      "!" = "Combine them manually"))
      ## create basic targets info.
      targets <- data.frame(sampleIDs = sample_IDs,
                            number = number,
                            pipe = rep("DIA", length(sample_IDs)),
                            enrichment = rep("protein", length(sample_IDs)),
                            row.names = sample_IDs)

      output <- list(targets = targets,
                     metadata = metadata)
    } else {
      # double check that everything is equal

      ## create basic targets info.
      targets <- data.frame(sampleIDs = sample_IDs,
                            number = number,
                            pipe = rep("DIA", length(sample_IDs)),
                            enrichment = rep("protein", length(sample_IDs)),
                            row.names = sample_IDs)

      ## use sample number info. to sort targets info. so that it is in the same order
      ## as the meta data file. combine targets and metadata into a single data.frame
      m <- match(metadata$number, number)
      targets <- targets[m, ]
      remove  <- colnames(metadata) %in% colnames(targets)
      targets <- cbind(targets, metadata[,remove==FALSE])


      cli::cli_inform(c("metadata and targets info combined"))
      cli::cli_rule()
      cli::cli_inform(c("v" = "Success"))

      output <- targets
    }
  }
  output
}


#' Import metadata file
#'
#' Imports a delimted metadata file, checking for some required columns.
#'
#' @inheritParams make_targets
#'
#' @return A dataframe of the imported metadata.
#' @export
#'
#' @examples
#' # No examples yet
#'
import_meta <-function(input_file,
                       sample_IDs) {

  ## check that the file is a csv, tsv, or text file
  filext <- file_extension(input_file)
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
#' @inheritParams make_targets
#'
#' @return A numeric vector of the extracted sample numbers
#' @export
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


