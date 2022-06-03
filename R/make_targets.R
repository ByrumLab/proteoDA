
#' Make a targets file
#'
#' Creates a dataframe of "targets". Output depends on whether a metadata file
#' is supplied, and whether the metadata and targets can be matched successfully
#' (see Value, below). Calls two subfunctions: \code{\link{import_meta}} and
#' \code{\link{get_dia_sample_number}}.
#'
#' @param sampleIDs A character vector of sample IDs, from which the targets
#'   file is constructed. In the DIA pipeline, these are the column names of the
#'   data slot in the object returned by \code{\link{extract_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sampleIDs have a required format:
#'   "Pool" samples must include the regex pattern \code{[P/p]ool_[[:digit:]]} in
#'   the column name, and samples must contain the regex pattern
#'   \code{[S/s]ample_[[:digit:]]} in the column name.
#' @param pipe Which analysis pipeline are you running? Options are "DIA", "TMT",
#'   "phosphoTMT". Output varies slightly based on the pipeline you
#'   choose.
#' @param enrich Another aspect of pipeline? Options are "protein" and "phospho".
#'   Seems semi-redundant?: DIA, TMT, and LF seem to always be "protein" in the
#'   code. PhosphoTMT can be either protein or phospho, though that seems to vary
#'   across subfunctions a little?
#' @param file Optional: a metadata file, giving info on the samples. Required
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
make_targets <- function(file = NULL,
                         sampleIDs,
                         pipe = c("DIA", "TMT", "phosphoTMT"),
                         enrich = c("protein", "phospho")) {

  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)

  cli::cli_rule()

  ## METADATA FILE INPUT BY USER
  if (!is.null(file)) {

    ## DIA METADATA FILE
    ## sample number extracted from sampleIDs is matched to sample number
    ## in meta file. pool samples in metadata number column should be labeled
    ## P1, P2, P3. samples in metadata number column should be labeled 1,2,3,4,
    ## if combining the two is successful a file named targets.csv is saved to
    ## the project directory.
    if (pipe == "DIA") { ## DIA

      ## import metadata file. use sample number as metadata key
      metadata <- import_meta(file = file,
                              sampleIDs = sampleIDs,
                              pipe = pipe)

      ## extract sample number from sampleIDs and use as data key
      number <- get_dia_sample_number(sampleIDs = sampleIDs)

      ## if the extracted sample numbers and metadata sample numbers match.
      ## create basic targets info. and combine with metadata. save targets
      ## to csv file and return combined metadata data.frame.
      if (all(number %in% metadata$number)) {
        # TODO change the check here to actual matching?? still returns true if number < metadata$number, which
        # necessitates the stopifnot check below.

        ## create basic targets info.
        targets <- data.frame(sampleIDs = sampleIDs,
                              number = number,
                              pipe = rep(pipe, length(sampleIDs)),
                              enrichment = rep(enrich, length(sampleIDs)),
                              row.names = sampleIDs)


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

        return(targets)

      } else {
        ## sample numbers do not match. return basic extracted targets info. and
        ## imported metadata as separate data.frames. integrate manually.
        ## NOTE: make sure the sample numbers for the pool samples in the imported
        ## metadata file are P1, P2, P3 and samples are 1,2,3,4 ...

        missing_samples <- number[number %notin% metadata$number]


        cli::cli_rule()
        cli::cli_warn(c("!" = "Not all samples found in column names are present in metadata",
                        "!" = "sample number{?s} {missing_samples} not found in metadata file",
                        "!" = "returning basic target info and metadata as separate dataframes within a list",
                        "!" = "Combine them manually"))

        ## create basic targets info.
        targets <- data.frame(sampleIDs = sampleIDs,
                              number = number,
                              pipe = rep(pipe, length(sampleIDs)),
                              enrichment = rep(enrich, length(sampleIDs)),
                              row.names = sampleIDs)

        data2 <- list(targets=targets, metadata=metadata)

        return(data2)
      }
    } ## DIA


    ## TMT/phosphoTMT METADATA FILE
    ## batch number and channel indexes are extracted from sampleIDs
    ## this combined key is used to match the combined batch/tag/channel info.
    ## in the metadata file.remember the batch column refers to TMT batch number
    ## in the metadata file.
    if (pipe == "TMT" | pipe == "phosphoTMT") { ## TMT/phosphoTMT

      ## import meta data file. set row names as batch.channel key
      metadata <- import_meta(file = file,
                              sampleIDs = sampleIDs,
                              pipe = pipe)
      rownames(metadata) <- paste(metadata$batch, metadata$tag, sep = ".")

      ## list of TMT tags and TMT batches in metadata file
      unq_tags    <- unique(metadata$tag)   ## list of unique TMT tags
      unq_bats    <- unique(metadata$batch) ## list of unique TMT batches
      num_unq_tags <- length(unq_tags)    ## no of unique TMT tags
      num_unq_bats <- length(unq_bats)    ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs) == num_unq_tags * num_unq_bats)

      ## use no. unique tags to determine type of TMT kit used in experiment.
      # TODO?: switch to ifelse statements?
      if (all(num_unq_tags==6  & all(unq_tags %in% TMT6plex)))  {plex <- TMT6plex}
      if (all(num_unq_tags==10 & all(unq_tags %in% TMT10plex))) {plex <- TMT10plex}
      if (all(num_unq_tags==11 & all(unq_tags %in% TMT11plex))) {plex <- TMT11plex}
      if (all(num_unq_tags==16 & all(unq_tags %in% TMT16plex))) {plex <- TMT16plex}
      if (num_unq_tags %notin% c(6,10,11,16)) {
        return(message("Error! TMTplex could not be determined."))
      }

      ## use TMT tags to add channel index number to metadata
      ## create named vector
      ## combine each TMT batch with the known TMT tag list
      ## e.g. (1.126, 1.130C, 2.126,2.130C) like metadata rownames)
      ## add channel column to metadata
      metadata$channel <- rep(NA, nrow(metadata))
      for (i in seq_along(plex)) {
        m <- metadata$tag %in% plex[i] # TODO == ?
        metadata$channel[m == TRUE] <-names(plex)[i] #TODO Just m, not m == TRUE?
      }
      if (any(is.na(metadata$channel))) {
        stop(paste("Error! adding channel index to metadata failed.",
                   "NAs in channel column remain:",
                   paste(metadata$channel,collapse=", ")))
      }

      ##--------------------------
      ##    PROCESS SAMPLEIDS
      ##--------------------------
      print("processing sampleIDs ...")

      ## extract channel index number from sampleIDs
      channel <- gsub("Reporter\\.intensity\\.corrected\\.", "", sampleIDs)
      channel <- as.numeric(sub("\\..*", "", channel))


      ## GET BATCH NUMBERS
      b <- get_tmt_batches(sampleIDs = sampleIDs,
                           pipe = pipe,
                           enrich = enrich)
      print(paste0("batch info. for ", length(unique(b)), " TMT batches created. Success!!"))

      ## unique TMT tags and TMT batches extracted from sampleIDs
      unq_idx     <- unique(channel)    ## list of unique TMT tags
      unq_b       <- unique(b)            ## list of unique TMT batches
      num_unq_idx  <- length(unq_idx) ## no of unique TMT tags
      num_unq_b    <- length(unq_b)     ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs) == num_unq_idx * num_unq_b) # TODO more informative error message

      ## use no. unique tags to determine type of TMT kit used in experiment.
      # TODO: switch to ifelse?
      if (all(num_unq_idx == 6  & all(unq_idx %in% names(TMT6plex))))  {plex2 <- TMT6plex}
      if (all(num_unq_idx == 10 & all(unq_idx %in% names(TMT10plex)))) {plex2 <- TMT10plex}
      if (all(num_unq_idx == 11 & all(unq_idx %in% names(TMT11plex)))) {plex2 <- TMT11plex}
      if (all(num_unq_idx == 16 & all(unq_idx %in% names(TMT16plex)))) {plex2 <- TMT16plex}
      if (num_unq_idx %in% c(6,10,11,16) == FALSE) {
        return(stop(paste("Error! TMTplex tags could not be determined from channel",
                          "index extracted from sampleIDs. :(")))
        }

      ## channel indexes from sampleIDs used to obtain TMT tag info.
      tag <- plex2[names(plex2)[channel]]

      ## basic targets info. created setting rownames as batch.tag (1.131C)
      ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
      targets <- data.frame(sampleIDs = sampleIDs,
                            batch = b,
                            tag = tag,
                            channel = channel,
                            pipe = rep(pipe,length(sampleIDs)),
                            enrichment = rep(enrich, length(sampleIDs)),
                            row.names=NULL) #TODO? addign rownames here, if we're gonna use them
      rownames(targets) <- paste(targets$batch,targets$tag,sep=".")


      ## TMT COMBINED
      if (all(rownames(targets) %in% rownames(metadata))) {
        ## merge metadata with targets information removing duplicate columns in targets
        targets <- targets[rownames(metadata), ]
        remove  <- colnames(targets) %in% colnames(metadata)
        targets2<-as.data.frame(targets[, remove == FALSE])
        colnames(targets2)<-colnames(targets)[remove == FALSE]
        targets <- cbind(metadata,targets2)
        stopifnot(targets$channel==metadata$channel)
        rownames(targets) <- targets$sampleIDs

        ## save targets file to project directory
        if(pipe=="phosphoTMT" & enrich=="protein"){ filename <- "targets.pro.csv" }
        if(pipe=="phosphoTMT" & enrich=="phospho"){ filename <- "targets.phos.csv" }
        if(pipe=="TMT" & enrich=="protein"){ filename <- "targets.csv" }
        utils::write.csv(targets, file=filename, row.names=FALSE)
        print(paste("targets file: ", filename, "saved to project directory. Success!!"))

        ## return the combined targets info. and imported metadata
        cat("\n")
        cat("targets metadata created. Success!!")
        return(targets)

      } ## TMT COMBINED


      ## TMT NOT COMBINED
      ## batch.channel index info. do not match. return basic extracted targets info.
      ## and imported metadata as separate data.frames. integrate manually.
      if(all(rownames(targets) %in% rownames(metadata))==FALSE){
        message("Error! The extracted batch.channel info. does not match the metadata file.\n",
                "The following batch.channel info. does not match the metadata file:\n",
                paste(rownames(targets)[rownames(targets) %in% rownames(metadata)==F],collapse="\n "),"\n")
        message("Basic targets info. (e.g. sampleIDs/batch/tag/channel)\n",
                "and the imported metadata could not be combined.\n",
                "Returning metadata and basic targets info. as separate data.frames.\n",
                "Integration failed. Combine the two data.frames manually. :(")

        data2 <- list(targets=targets, metadata=metadata)
        return(data2)


      } ## TMT NOT COMBINED

    } ## TMT/PHOSPHOTMT


  } else {

    ## DIA (NO META DATA FILE)
    if (pipe == "DIA") { ## DIA NULL
      ## extract sample number from sampleIDs and use as data key
      number  <- get_dia_sample_number(sampleIDs = sampleIDs)

      targets <- data.frame(sampleIDs = sampleIDs,
                            number = number,
                            pipe = rep(pipe, length(sampleIDs)),
                            enrichment= rep(enrich, length(sampleIDs)),
                            row.names=sampleIDs)

      cli::cli_rule()
      cli::cli_warn(c("!" = "metadata file not provided",
                      "!" = "returning {enrich} targets template for {pipe} experiment",
                      "!" = "import metadata for these samples and add it to the returned target template dataframe"))

      return(targets)

    } ## DIA NULL


    ## TMT/phosphoTMT (NO METADATA FILE)
    if (pipe == "TMT" | pipe == "phosphoTMT") { ## TMT/PHOSPHO NULL

      ## extract channel index number from sampleIDs
      channel <- gsub("Reporter\\.intensity\\.corrected\\.","",sampleIDs); channel
      channel <- sub("\\..*","",channel)
      channel<-as.numeric(channel)
      ## extract batch number from sampleIDs
      b <- get_tmt_batches(sampleIDs=sampleIDs, pipe=pipe, enrich=enrich)

      ## list of TMT tags and TMT batches
      unq_idx     <- unique(channel)   ## list of unique TMT tags
      unq_b       <- unique(b)           ## list of unique TMT batches
      num_unq_idx  <- length(unq_idx)## no of unique TMT tags
      num_unq_b    <- length(unq_b)    ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs)==num_unq_idx * num_unq_b)

      ## use no. unique tags to determine type of TMT kit used in experiment.
      if(all(num_unq_idx==6  & all(unq_idx %in% names(TMT6plex)))){ plex2  <- TMT6plex }
      if(all(num_unq_idx==10 & all(unq_idx %in% names(TMT10plex)))){ plex2 <- TMT10plex }
      if(all(num_unq_idx==11 & all(unq_idx %in% names(TMT11plex)))){ plex2 <- TMT11plex }
      if(all(num_unq_idx==16 & all(unq_idx %in% names(TMT16plex)))){ plex2 <- TMT16plex }
      plex2
      ## channel indexes from sampleIDs used to obtain TMT tag info.
      if(!is.null(plex2)){ tag <- plex2[names(plex2)[channel]] }
      if(is.null(plex2)){
        tag=rep("", length(sampleIDs))
        message(paste("Warning! TMTplex tags could not be determined from channel",
                      "index values extracted from sampleIDs. :("))
      }

      ## create tagets basic info. data.frame.
      ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
      targets <- data.frame(sampleIDs=sampleIDs, batch=b, tag=tag, channel=channel,
                            pipe=rep(pipe,length(sampleIDs)),
                            enrichment=rep(enrich, length(sampleIDs)), row.names=NULL)
      rownames(targets) <- targets$sampleIDs
      print(paste(enrich, "basic targets info. for", pipe, "experiment created. Success!!",
                  "Import metadata and combine the two data.frames."))

      return(targets)

    } ## TMT/PHOSPHO NULL


  } ## NULL
}


#' Import metadata file
#'
#' Imports a delimted metadata file, checking for some required columns.
#'
#' @param file The filename of the metadata file to import
#' @param sampleIDs A character vector of sample IDs, from which the targets
#'   file is constructed. In the DIA pipeline, these are the column names of the
#'   data slot in the object returned by \code{\link{extract_data}}, e.g.,
#'   \code{colnames(ext$data)}. Currently, these sampleIDs have a required format:
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
                       sampleIDs,
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
get_dia_sample_number <- function(sampleIDs) {

  ## scaffold DIA sample reports should contain
  ## Pool_1/Sample_1 as part of the sampleIDs
  ## e.g. "Balachandran_061721_Pool_3.mzML",
  ## e.g. "Balachandran_061721_Sample_01_CS1.mzML"
  pattern1 <- "[P/p]ool_[[:digit:]]"   ## POOL PATTERN (Pool_1)
  pattern2 <- "[S/s]ample_[[:digit:]]" ## SAMPLE PATTERN (Sample_1)
  pat1     <- grep(pattern1, sampleIDs)
  pat2     <- grep(pattern2, sampleIDs)
  pattern  <- c(pat1, pat2)

  if(!all(length(pattern) == length(sampleIDs))) { # Double-check this check

    supplied_ids <- cli::cli_vec(sampleIDs, list(vec_trunc = 5))
    cli::cli_abort(c("sample numbers could not be extracted from sampleIDs",
                     "x" = "Pool sample(s) must contain a string matching {.code {pattern1}}",
                     "!" = "For example: \"Pool_1\" or \"pool_01\"",
                     "x" = "Other sample(s) must contain a string matching {.code {pattern2}}",
                     "!" = "For example: \"Sample_1\" or \"sample_01\"",
                     "!" = "Supplied sample ids: {.val {supplied_ids}}"))
  } else {
    ## pool samples
    if (length(pat1) != 0) {
      number1 <- gsub("\\.", "_", sampleIDs[pat1])
      number1 <- gsub(".*[P/p]ool_", "", number1)  ## removes everything b4 Pool_
      number1 <- gsub("_.*", "", number1)  ## removes everything after number _mzML
      number1 <- paste0("P", as.numeric(number1)) ## adds capital P to number == P1, P2, P3
    } else {
      number1 <- NULL
    }

    ## samples
    number2 <- gsub("\\.","_", sampleIDs[pat2])
    number2 <- gsub(".*[S/s]ample_", "", number2) ## removes everything b4 Sample_
    number2 <- gsub("_.*", "", number2)   ## removes everthing after number _mzML
    number2 <- as.numeric(number2)      ## if 01,02 == 1,2
    number  <- c(number1, number2)       ## P1,P2,P3, 1,2,3,..., n

    cli::cli_inform("sample numbers extracted from DIA sampleIDs.")
    return(number)
  }
}


## TO BE DOCUMENTED STILL
get_tmt_batches <- function(sampleIDs,
                            pipe = c("DIA","TMT","phosphoTMT"),
                            enrich = c("protein","phospho")) {
  # Check args
  rlang::arg_match(pipe)
  rlang::arg_match(enrich)

  ## remove Reporter.intensity.corrected.XX text from sampleIDs
  pattern <- gsub("Reporter\\.intensity\\.corrected\\.[[:digit:]]+", "", sampleIDs)

  ## determine if remaining experiment info. text contains a number
  ## if a number is present we assume this refers to different batches
  hasNumber <- grep("[[:digit:]]", pattern)

  ## if experiment info. part of sampleID contains a number
  ## then isolate number by removing periods, underscores, and letters
  ## from the text
  if (all(length(hasNumber) == length(sampleIDs))) {

    b <- gsub("\\.","",pattern)   ## removes period
    b <- gsub("_","",b)         ## removes underscore
    b <- gsub("[A-Z]","",b)      ## removes capital letters
    b <- gsub("[a-z]","",b)    ## removes lowercase letters

    ## check that b is now a number and free of other
    ## characters and symbols by checking against a vector of 1:100
    if (all(b %in% 1:100)) {
      b <- as.numeric(b)
      } else {
      stop(paste("b is either not a number or is a number > 100. :(",
                 "unique batch values include:", paste(unique(b), collapse=", ")))
      }
  } else {
    ## if experiment text does not contain a number then use unique
    ## values of the text to assign batches assuming that the first
    ## unique value corresponds to experiment 1
    b=pattern
    for (i in 1:length(unique(b))) {
      b[b %in% unique(b)[i]] <- i
    }

    if (all(b %in% 1:100)) {
      b <- as.numeric(b)
      } else {
      stop(paste("b is either not a number or is a number > 100. :(",
                 "unique batch values include:", paste(unique(b), collapse=", ")))
        }
  }
  cat("\n");print(b);cat("\n");print(unique(b))

  return(b)
}
