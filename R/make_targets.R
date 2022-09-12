
#' Make a targets file
#'
#' Creates a dataframe of "targets". Output depends on whether a metadata file
#' is supplied, and whether the metadata and targets can be matched successfully
#' (see Value, below). Calls two subfunctions: \code{\link{import_meta}} and
#' \code{\link{get_dia_sample_number}}.
#'
#' @param sample_IDs A character vector of sample IDs, from which the targets
#'   file is constructed. In the DIA pipeline, these are the column names of the
#'   data slot in the object returned by \code{\link{read_DIA_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sample_IDs have a required format:
#'   "Pool" samples must include the regex pattern \code{[P/p]ool_[[:digit:]]} in
#'   the column name, and samples must contain the regex pattern
#'   \code{[S/s]ample_[[:digit:]]} in the column name.
#' @param input_file Optional: a metadata file, giving info on the samples. Required
#'   format and columns for metadata file vary by pipeline.
#'   Need to document somewhere.
#'
#' @return It depends:
#'   \itemize{
#'     \item If a metadata file is supplied and could be successfully matched
#'       against the input column names, a full targets dataframe is returned,
#'       and a targets.csv file is created in the working directory.
#'     \item If a metadata file is supplied but matching against column names
#'       was unsuccessful, returns, with a warning, a list with two slots:
#'       (1) the dataframe of assembled targets, and (2) the dataframe of
#'       imported metadata. These must be manually combined before further
#'       analysis. No targets.csv file is created.
#'     \item When a metadata file is not supplied, returns, with a warning, an
#'       incomplete targets dataframe, and does not create a targets.csv file.
#'   }
#'
#' @export
#'
#' @examples
#' # No examples yet
make_targets <- function(input_file = NULL,
                         sample_IDs) {

  cli::cli_rule()

  ## METADATA FILE INPUT BY USER
  if (!is.null(file)) {

    ## DIA METADATA FILE
    ## sample number extracted from sample_IDs is matched to sample number
    ## in meta file. pool samples in metadata number column should be labeled
    ## P1, P2, P3. samples in metadata number column should be labeled 1,2,3,4,
    ## if combining the two is successful a file named targets.csv is saved to
    ## the project directory.

    ## import metadata file. use sample number as metadata key
    metadata <- import_meta(file = file,
                            sample_IDs = sample_IDs,
                            pipe = pipe)

    ## extract sample number from sample_IDs and use as data key
    number <- get_dia_sample_number(sample_IDs = sample_IDs)

    ## if the extracted sample numbers and metadata sample numbers match.
    ## create basic targets info. and combine with metadata. save targets
    ## to csv file and return combined metadata data.frame.
    if (all(number %in% metadata$number)) {
      # TODO change the check here to actual matching?? still returns true if number < metadata$number, which
      # necessitates the stopifnot check below.

      ## create basic targets info.
      targets <- data.frame(sample_IDs = sample_IDs,
                            number = number,
                            pipe = rep("DIA", length(sample_IDs)),
                            enrichment = rep("protein", length(sample_IDs)),
                            row.names = sample_IDs)


      ## use sample number info. to sort targets info. so that it is in the same order
      ## as the meta data file. combine targets and metadata into a single data.frame
      m <- match(metadata$number, number)
      stopifnot(number[m] == metadata$number) # TODO MORE INFORMATIVE MESSAGE HERE
      targets <- targets[m, ]
      remove  <- colnames(metadata) %in% colnames(targets)
      targets <- cbind(targets, metadata[,remove==FALSE])

      ## save targets as a csv file in the main project directory.
      filename <- "targets.csv"
      utils::write.csv(targets, file = filename, row.names = FALSE)

      cli::cli_inform(c("{pipe} targets file {.file {filename}} saved to current working directory"))
      cli::cli_inform(c("metadata and targets info combined"))
      cli::cli_rule()
      cli::cli_inform(c("v" = "Success"))

      output <- targets

    } else {
      ## sample numbers do not match. return basic extracted targets info. and
      ## imported metadata as separate data.frames. integrate manually.
      missing_samples <- number[number %notin% metadata$number]
      cli::cli_rule()
      cli::cli_warn(c("!" = "Not all samples found in column names are present in metadata",
                      "!" = "sample number{?s} {missing_samples} not found in metadata file",
                      "!" = "returning basic target info and metadata as separate dataframes within a list",
                      "!" = "Combine them manually"))
      ## create basic targets info.
      targets <- data.frame(sample_IDs = sample_IDs,
                            number = number,
                            pipe = rep("DIA", length(sample_IDs)),
                            enrichment = rep("protein", length(sample_IDs)),
                            row.names = sample_IDs)

      output <- list(targets = targets,
                     metadata = metadata)
    }

  } else {
    ## extract sample number from sample_IDs and use as data key
    number  <- get_dia_sample_number(sample_IDs = sample_IDs)

    targets <- data.frame(sample_IDs = sample_IDs,
                          number = number,
                          pipe = rep("DIA", length(sample_IDs)),
                          enrichment= rep("protein", length(sample_IDs)),
                          row.names=sample_IDs)

    cli::cli_rule()
    cli::cli_warn(c("!" = "metadata file not provided",
                    "!" = "returning protein targets template for a DIA experiment",
                    "!" = "import metadata for these samples and add it to the returned target template dataframe"))
    output <- targets
  }
  output
}


#' Import metadata file
#'
#' Imports a delimted metadata file, checking for some required columns.
#'
#' @param file The filename of the metadata file to import
#' @param sample_IDs A character vector of sample IDs, from which the targets
#'   file is constructed. In the DIA pipeline, these are the column names of the
#'   data slot in the object returned by \code{\link{read_DIA_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sample_IDs have a required format:
#'   "Pool" samples must include the regex pattern \code{[P/p]ool_[[:digit:]]} in
#'   the column name, and samples must contain the regex pattern
#'   \code{[S/s]ample_[[:digit:]]} in the column name.
#' @inheritParams make_targets
#'
#' @return A dataframe of the imported metadata.
#' @export
#'
#' @examples
#' # No examples yet
#'
import_meta <-function(file,
                       sample_IDs,
                       pipe = c("DIA", "TMT", "phosphoTMT")) {
  # Check args
  rlang::arg_match(pipe)


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

  if (filext == "txt" | filext=="tsv") {sep= "\t"}
  if (filext == "csv") {sep=","}

  ## the metadata file should contain the following required columns
  if (pipe == "phosphoTMT" | pipe=="TMT") {reqCols<- tmtMetaColums}
  if (pipe == "DIA") {reqCols <- diaMetaColums}
  if (pipe == "LF") {reqCols <- lfMetaColums}

  ## import file
  data <- utils::read.csv(file = file,
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


