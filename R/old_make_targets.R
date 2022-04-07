

make_targets <- function(file=NULL, sampleIDs, pipe=c("DIA","TMT","phosphoTMT"),
                         enrich=c("protein","phospho")){

  ## METADATA FILE INPUT BY USER
  if(!is.null(file)){

    ## DIA METADATA FILE
    ## sample number extracted from sampleIDs is matched to sample number
    ## in meta file. pool samples in metadata number column should be labeled
    ## P1, P2, P3. samples in metadata number column should be labeled 1,2,3,4,
    ## if combining the two is successful a file named targets.csv is saved to
    ## the project directory.
    if(pipe=="DIA"){ ## DIA

      ## import metadata file. use sample number as metadata key
      metadata <- import_meta(file=file, sampleIDs=sampleIDs, pipe=pipe)

      ## extract sample number from sampleIDs and use as data key
      number <- get_dia_sample_number(sampleIDs=sampleIDs);number

      ## if the extracted sample numbers and metadata sample numbers match.
      ## create basic targets info. and combine with metadata. save targets
      ## to csv file and return combined metadata data.frame.
      if(all(number %in% metadata$number)){
        ## create basic targets info.
        targets <- data.frame(sampleIDs=sampleIDs, number=number,
                              pipe=rep(pipe,length(sampleIDs)),
                              enrichment=rep(enrich,length(sampleIDs)))
        rownames(targets) <- targets$sampleIDs;head(targets)

        ## use sample number info. to sort targets info. so that it is in the same order
        ## as the meta data file. combine targets and metadata into a single data.frame
        m <- match(metadata$number, number);m
        stopifnot(number[m]==metadata$number)
        targets <- targets[m, ];head(targets)
        remove  <- colnames(metadata)%in%colnames(targets);remove
        targets <- cbind(targets, metadata[,remove==FALSE]);head(targets)

        ## save targets as a csv file in the main project directory.
        filename <- "targets.csv"
        utils::write.csv(targets, file=filename, row.names=FALSE)
        print(paste(pipe, "targets file: ", filename, "saved to project directory. Success!!"))
        print("metadata and targets info. combined. Success!!")
        return(targets)

      } ## COMBINED DIA

      ## sample numbers do not match. return basic extracted targets info. and
      ## imported metadata as separate data.frames. integrate manually.
      ## NOTE: make sure the sample numbers for the pool samples in the imported
      ## metadata file are P1, P2, P3 and samples are 1,2,3,4 ...
      if(all(number %in% metadata$number)==FALSE){
        message("Error! The extracted sample numbers do not match the metadata file.\n",
                "The following sample numbers are not present in the metadata file:",
                paste(number[number %in% metadata$number==F],collapse=", "),"\n")
        message("Basic targets info. (e.g. sampleIDs/number) extracted from the input\n",
                " sampleIDs and the imported metadata could not be combined.\n",
                "Returning metadata and basic targets info. as separate data.frames.\n",
                "Integration failed. Combine manually :(")

        ## create basic targets info.
        targets <- data.frame(sampleIDs=sampleIDs, number=number,
                              pipe=rep(pipe,length(sampleIDs)),
                              enrichment=rep(enrich,length(sampleIDs)))
        rownames(targets) <- targets$sampleIDs;head(targets)

        data2 <- list(targets=targets, metadata=metadata)
        return(data2)

      } ## NOT COMBINED DIA


      # } ## PATTERN 1

    } ## DIA


    ## TMT/phosphoTMT METADATA FILE
    ## batch number and channel indexes are extracted from sampleIDs
    ## this combined key is used to match the combined batch/tag/channel info.
    ## in the metadata file.remember the batch column refers to TMT batch number
    ## in the metadata file.
    if(pipe=="TMT" | pipe=="phosphoTMT"){ ## TMT/phosphoTMT

      ## import meta data file. set row names as batch.channel key
      metadata <- import_meta(file=file, sampleIDs=sampleIDs, pipe=pipe);metadata
      rownames(metadata) <- paste(metadata$batch, metadata$tag,sep=".");metadata

      ## list of TMT tags and TMT batches in metadata file
      unq_tags    <- unique(metadata$tag);unq_tags   ## list of unique TMT tags
      unq_bats    <- unique(metadata$batch);unq_bats ## list of unique TMT batches
      no_unq_tags <- length(unq_tags);no_unq_tags    ## no of unique TMT tags
      no_unq_bats <- length(unq_bats);no_unq_bats    ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs)==no_unq_tags * no_unq_bats)

      ## use no. unique tags to determine type of TMT kit used in experiment.
      if(all(no_unq_tags==6  & all(unq_tags %in% TMT6plex))){ plex  <- TMT6plex }
      if(all(no_unq_tags==10 & all(unq_tags %in% TMT10plex))){ plex <- TMT10plex }
      if(all(no_unq_tags==11 & all(unq_tags %in% TMT11plex))){ plex <- TMT11plex }
      if(all(no_unq_tags==16 & all(unq_tags %in% TMT16plex))){ plex <- TMT16plex }
      if(no_unq_tags %in% c(6,10,11,16)==FALSE){
        return(message("Error! TMTplex could not be determined."))
      }; plex

      ## use TMT tags to add channel index number to metadata
      ## create named vector
      ## combine each TMT batch with the known TMT tag list
      ## e.g. (1.126, 1.130C, 2.126,2.130C) like metadata rownames)
      ## add channel column to metadata
      metadata$channel <- rep(NA,nrow(metadata))
      for(i in base::seq_along(plex)){
        m <- metadata$tag %in% plex[i]
        metadata$channel[m==TRUE] <-names(plex)[i]
      }
      if(any(is.na(metadata$channel))==TRUE){
        stop(paste("Error! adding channel index to metadata failed.",
                   "NAs in channel column remain:",
                   paste(metadata$channel,collapse=", ")))
      }



      ##--------------------------
      ##    PROCESS SAMPLEIDS
      ##--------------------------
      print("processing sampleIDs ...")

      ## extract channel index number from sampleIDs
      channel <- gsub("Reporter\\.intensity\\.corrected\\.","",sampleIDs); channel
      channel <- as.numeric(sub("\\..*","",channel));channel


      ## GET BATCH NUMBERS
      b <- get_tmt_batches(sampleIDs=sampleIDs, pipe=pipe, enrich=enrich);b
      print(paste0("batch info. for ", length(unique(b)), " TMT batches created. Success!!"))

      ## unique TMT tags and TMT batches extracted from sampleIDs
      unq_idx     <- unique(channel);unq_idx    ## list of unique TMT tags
      unq_b       <- unique(b);unq_b            ## list of unique TMT batches
      no_unq_idx  <- length(unq_idx);no_unq_idx ## no of unique TMT tags
      no_unq_b    <- length(unq_b);no_unq_b     ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs)==no_unq_idx * no_unq_b)

      ## use no. unique tags to determine type of TMT kit used in experiment.
      if(all(no_unq_idx==6  & all(unq_idx %in% names(TMT6plex)))){  plex2  <- TMT6plex }
      if(all(no_unq_idx==10 & all(unq_idx %in% names(TMT10plex)))){ plex2 <- TMT10plex }
      if(all(no_unq_idx==11 & all(unq_idx %in% names(TMT11plex)))){ plex2 <- TMT11plex }
      if(all(no_unq_idx==16 & all(unq_idx %in% names(TMT16plex)))){ plex2 <- TMT16plex }
      if(no_unq_idx %in% c(6,10,11,16)==FALSE){
        return(stop(paste("Error! TMTplex tags could not be determined from channel",
                          "index extracted from sampleIDs. :("))) }
      plex2

      ## channel indexes from sampleIDs used to obtain TMT tag info.
      tag <- plex2[names(plex2)[channel]];tag

      ## basic targets info. created setting rownames as batch.tag (1.131C)
      ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
      targets <- data.frame(sampleIDs=sampleIDs, batch=b, tag=tag, channel=channel,
                            pipe=rep(pipe,length(sampleIDs)),
                            enrichment=rep(enrich, length(sampleIDs)), row.names=NULL)
      rownames(targets) <- paste(targets$batch,targets$tag,sep=".");head(targets)


      ## TMT COMBINED
      if(all(rownames(targets) %in% rownames(metadata))){
        ## merge metadata with targets information removing duplicate columns in targets
        targets <- targets[rownames(metadata), ]
        remove  <- colnames(targets) %in% colnames(metadata);table(remove)
        targets2<-as.data.frame(targets[,remove==FALSE])
        colnames(targets2)<-colnames(targets)[remove==FALSE]; head(targets2)
        targets <- cbind(metadata,targets2);head(targets)
        stopifnot(targets$channel==metadata$channel)
        rownames(targets) <- targets$sampleIDs;head(targets)

        ## save targets file to project directory
        if(pipe=="phosphoTMT" & enrich=="protein"){ filename <- "targets.pro.csv" }
        if(pipe=="phosphoTMT" & enrich=="phospho"){ filename <- "targets.phos.csv" }
        if(pipe=="TMT" & enrich=="protein"){ filename <- "targets.csv" }
        utils::write.csv(targets, file=filename, row.names=FALSE)
        print(paste("targets file: ", filename, "saved to project directory. Success!!"))

        ## return the combined targets info. and imported metadata
        cat("\n");print("targets metadata created. Success!!")
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


  } ## NOT NULL


  ## NO METADATA FILE
  if(is.null(file)){

    ## DIA (NO META DATA FILE)
    if(pipe=="DIA"){ ## DIA NULL
      ## extract sample number from sampleIDs and use as data key
      number  <- get_dia_sample_number(sampleIDs=sampleIDs);number
      targets <-data.frame(sampleIDs=sampleIDs, number=number,
                           pipe=rep(pipe,length(sampleIDs)),
                           enrichment=rep(enrich,length(sampleIDs)), row.names=NULL)
      rownames(targets)<-targets$sampleIDs;head(targets)
      print(paste(enrich, "targets template for", pipe, "experiment created. Success!!",
                  "Import metadata and combine."))

      return(targets)

    } ## DIA NULL


    ## TMT/phosphoTMT (NO METADATA FILE)
    if(pipe=="TMT" | pipe=="phosphoTMT"){ ## TMT/PHOSPHO NULL

      ## extract channel index number from sampleIDs
      channel <- gsub("Reporter\\.intensity\\.corrected\\.","",sampleIDs); channel
      channel <- sub("\\..*","",channel);channel
      channel<-as.numeric(channel);channel
      ## extract batch number from sampleIDs
      b <- get_tmt_batches(sampleIDs=sampleIDs, pipe=pipe, enrich=enrich);b

      ## list of TMT tags and TMT batches
      unq_idx     <- unique(channel);unq_idx   ## list of unique TMT tags
      unq_b       <- unique(b);unq_b           ## list of unique TMT batches
      no_unq_idx  <- length(unq_idx);no_unq_idx## no of unique TMT tags
      no_unq_b    <- length(unq_b);no_unq_b    ## no of unique batches
      ## stop if each batch does not contain same no samples
      stopifnot(length(sampleIDs)==no_unq_idx * no_unq_b)

      ## use no. unique tags to determine type of TMT kit used in experiment.
      if(all(no_unq_idx==6  & all(unq_idx %in% names(TMT6plex)))){ plex2  <- TMT6plex }
      if(all(no_unq_idx==10 & all(unq_idx %in% names(TMT10plex)))){ plex2 <- TMT10plex }
      if(all(no_unq_idx==11 & all(unq_idx %in% names(TMT11plex)))){ plex2 <- TMT11plex }
      if(all(no_unq_idx==16 & all(unq_idx %in% names(TMT16plex)))){ plex2 <- TMT16plex }
      plex2
      ## channel indexes from sampleIDs used to obtain TMT tag info.
      if(!is.null(plex2)){ tag <- plex2[names(plex2)[channel]];tag }
      if(is.null(plex2)){
        tag=rep("", length(sampleIDs))
        message(paste("Warning! TMTplex tags could not be determined from channel",
                      "index values extracted from sampleIDs. :("))
      }

      ## create tagets basic info. data.frame.
      ## https://stackoverflow.com/questions/23534066/cbind-warnings-row-names-were-found-from-a-short-variable-and-have-been-discar
      targets <- data.frame(sampleIDs=sampleIDs, batch=b, tag=tag, channel=channel,
                            pipe=rep(pipe,length(sampleIDs)),
                            enrichment=rep(enrich, length(sampleIDs)), row.names=NULL);head(targets)
      rownames(targets) <- targets$sampleIDs;targets;dim(targets)
      print(paste(enrich, "basic targets info. for", pipe, "experiment created. Success!!",
                  "Import metadata and combine the two data.frames."))

      return(targets)

    } ## TMT/PHOSPHO NULL


  } ## NULL
}


import_meta <-function(file, sampleIDs, pipe){

  ## check that file exists in specified location.
  ## if file exits and is csv/tsv/txt
  file <-gsub("\\./","",file)
  file.dir <- dirname(file);file.dir
  file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),
                                        several.ok=FALSE))
  filext <- tools::file_ext(file);filext
  filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=FALSE)
  if(filext=="txt" | filext=="tsv"){ sep="\t" }
  if(filext=="csv"){ sep=","}

  ## the metadata file should contain the following required columns
  if(pipe=="phosphoTMT" | pipe=="TMT"){ reqCols<- tmtMetaColums }
  if(pipe=="DIA"){ reqCols <- diaMetaColums }
  if(pipe=="LF"){ reqCols <- lfMetaColums }

  ## import file
  data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE,
                          header=TRUE, check.names=TRUE)

  if(all(reqCols %in% colnames(data))==TRUE){
    print("metadata file imported. Success!!")
    return(data)
  }

  if(all(reqCols %in% colnames(data))==FALSE){
    stop(message(paste("Error! For",pipe,"experiments the sample metadata file",
                       "should include the following required columns:",
                       paste(reqCols,collapse=", "), ". The following columns are",
                       "missing from the metadata file: ",
                       paste(reqCols[reqCols %in% colnames(data)==FALSE],
                             collapse=", "))))
  }

}


get_dia_sample_number <- function(sampleIDs){

  ## scaffold DIA sample reports should contain
  ## Pool_1/Sample_1 as part of the sampleIDs
  ## e.g. "Balachandran_061721_Pool_3.mzML",
  ## e.g. "Balachandran_061721_Sample_01_CS1.mzML"
  pattern1 <- "[P/p]ool_[[:digit:]]"   ## POOL PATTERN (Pool_1)
  pattern2 <- "[S/s]ample_[[:digit:]]" ## SAMPLE PATTERN (Sample_1)
  pat1     <- grep(pattern1,sampleIDs)
  pat2     <- grep(pattern2,sampleIDs)
  pattern  <- c(pat1, pat2);pattern

  if(all(length(pattern)==length(sampleIDs))==TRUE){

    ## pool samples
    number1 <- gsub("\\.","_", sampleIDs[pat1]);number1
    number1 <- gsub(".*[P/p]ool_","",number1)  ## removes everything b4 Pool_
    number1 <- gsub("_.*","",number1);number1  ## removes everything after number _mzML
    number1 <- paste0("P",as.numeric(number1)) ## adds capital P to number == P1, P2, P3

    ## samples
    number2 <- gsub("\\.","_",sampleIDs[pat2])
    number2 <- gsub(".*[S/s]ample_","",number2) ## removes everything b4 Sample_
    number2 <- gsub("_.*","",number2);number2   ## removes everthing after number _mzML
    number2 <- as.numeric(number2);number2      ## if 01,02 == 1,2
    number  <- c(number1, number2);number       ## P1,P2,P3, 1,2,3,..., n

    print("sample numbers extracted from DIA sampleIDs. Success!!")
    return(number)
  }

  if(all(length(pattern)!=length(sampleIDs))==FALSE){
    stop(message(paste("Error! Sample numbers could not be extracted",
                       "from the input sampleIDs. DIA sampleIDs should",
                       "contain the following patterns:", pattern1,
                       "(pool samples e.g. P/pool_1/P/pool_01) & ",
                       pattern2,"(samples e.g. Sample_1/Sample_01 (not S1)).")))
  }

}
