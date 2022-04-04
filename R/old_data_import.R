
extract_data <- function(file=NULL, sampleIDs=NULL, min.prob=0.75,
                         pipe=c("DIA","TMT","phosphoTMT","LF"),
                         enrich=c("protein","phospho")){

  ## if no file is input by user the open file dialog box is
  ## called to allow the user to navigate to and select the file.
  if(is.null(file)){ file=file.choose() }

  ## input parameters are added to param variable
  param <- stats <- list()
  param[["file"]]      <- file
  param[["pipe"]]      <- pipe
  param[["enrich"]]    <- enrich
  param[["sampleIDs"]] <- ifelse(is.null(sampleIDs),"NULL", paste(sampleIDs,collapse=", "))
  param[["min.prob"]]  <- ifelse(enrich=="protein", "NULL",min.prob)


  ## DIA PROTEIN
  if(pipe=="DIA" & enrich=="protein"){

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    ## the first 10 lines are skipped and last line removed. column 1 (.X)
    ## is renamed "id"
    DIADATA <- import_data(file=file, pipe=pipe, enrich=enrich)
    diaData <- DIADATA$data
    reqCols <- DIADATA$reqCols

    param[["file"]] <- DIADATA$file
    stats[["total_input_rows"]] <- DIADATA$no_input

    ## REMOVE IRRELEVANT COLUMNS
    ## irrelevant columns, columns 2-3 (visible,star) are removed from
    ## the sample report
    if(colnames(diaData)[1] != "id"){ colnames(diaData)[1] <-"id" }
    remove  <- colnames(diaData) %in% c("Visible", "Star");remove
    diaData <- diaData[, remove==FALSE]


    ## REPLACE MISSING VALUES WITH ZERO
    ## intensity data is not a numeric data.frame and requires
    ## additional processing numbers are character strings separated
    ## by commas (4,000,000), some cells contain text 'Missing Value'
    ## if the protein does not have detectable expression.
    ## replace 'Missing Value' with zeros, remove comma's, convert
    ## data in each column to numeric values
    diaData[,][diaData[,]=="Missing Value"] <- "0"
    replaceCommas  <- function(x){ x <- as.numeric(gsub("\\,", "", x)) }
    replaceNumbers <- colnames(diaData) %in% reqCols
    replaceColums  <- colnames(diaData)[replaceNumbers==FALSE];replaceColums
    for(i in replaceColums){ diaData[,i] <- replaceCommas(diaData[,i]) }


    ## SAVE BIG QUERY INPUT FILE
    ## save a copy of samples report for upload into Big Query.
    ## This file cannot contain NA or blank values, column names
    ## cannot begin with a number, replace all NA and blank cells
    ## with zeros. If column names begin with a number add X to
    ## beginning of column name. change file name to ilab_Samples_Report_BQ.csv
    bqData <- data.frame(diaData)
    bqData[,][is.na(bqData[,])] <- 0
    bqData[,][bqData[,]==""]    <- 0

    ## if column names begin with a number append X to the start of each name.
    testColumNumber <-substr(colnames(bqData),1,1);testColumNumber
    if(length(grep("[[:digit:]]",testColumNumber)) > 0){
      colnames(bqData) <- paste0("X",colnames(bqData))
    }
    bn   <- basename(file);bn
    bn   <- gsub("Samples Report of ","",bn);bn
    ilab <- gsub(paste0(".",tools::file_ext(bn)),"",bn);ilab
    ## save DIA Big Query Input File
    filename <- paste0(ilab,"_Samples_Report_BQ.csv")
    utils::write.csv(bqData, file.path(".",filename), row.names=FALSE)
    print(paste("Big Querry Samples Report (",filename,
                ") saved to project directory. Success!!"))
    param[["ilab"]] <- ilab


    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data=diaData, pipe=pipe, enrich=enrich)
    stats[["no_contam_removed"]] <- qfilterData$no_contam

    ## EXTRACT PROTEIN DATA
    extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs,
                                    pipe=pipe, enrich=enrich)

    stats[["total_input_samples"]] <- extProt$no_samples
    stats[["no_extracted_rows"]]   <- extProt$no_extracted

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)

    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot,
                  param=logs$param, stats=logs$stats[c(3,1,2,4),], ilab=ilab)
    return(data2)


  } ## DIA PROTEIN


  ## TMT PROTEIN
  if((pipe=="TMT" | pipe=="LF") & enrich=="protein"){

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
    stats[["total_input_rows"]] <- maxQuant$no_input

    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
    stats[["no_contam_removed"]] <- qfilterData$no_contam


    ## EXTRACT PROTEIN DATA AND ANNOTATION
    extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs,
                                    pipe=pipe, enrich=enrich)

    stats[["total_input_samples"]] <- extProt$no_samples
    stats[["no_extracted_rows"]]   <- extProt$no_extracted
    param[["pattern"]] <- extProt$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)

    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot,
                  param=logs$param,stats=logs$stats)
    return(data2)

  } ## TMT PROTEIN


  ## phosphoTMT PROTEIN
  if(pipe=="phosphoTMT" & enrich=="protein"){

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("protein_analysis/01_quality_control")){
    #    dir.create("protein_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
    stats[["total_input_rows"]] <- maxQuant$no_input


    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
    stats[["no_contam_removed"]] <- qfilterData$no_contam


    ## EXTRACT PROTEIN DATA
    extProt <- extract_protein_data(data=qfilterData$data, sampleIDs=sampleIDs,
                                    pipe=pipe, enrich=enrich)

    stats[["total_input_samples"]] <- extProt$no_samples
    stats[["no_extracted_rows"]] <- extProt$no_extracted
    param[["pattern"]] <- extProt$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)


    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data, corresponding extracted annotation,
    ## stats, and input parameters
    data2 <- list(data=extProt$rawData, annot=extProt$rawAnnot,
                  param=logs$param, stats=logs$stats[c(3,1,2,4),])
    return(data2)

  } ## phosphoTMT PROTEIN


  ## phosphoTMT PHOSPHO
  if(pipe=="phosphoTMT" & enrich=="phospho"){

    ## CREATE QUALITY CONTROL OUTPUT DIRECTORY
    # if(!dir.exists("phospho_analysis/01_quality_control")){
    #    dir.create("phospho_analysis/01_quality_control", recursive=TRUE)
    #    print("quality control directory created...")
    # }

    ## IMPORT DATA
    maxQuant <- import_data(file=file, pipe=pipe, enrich=enrich)
    stats[["total_input_rows"]] <- maxQuant$no_input

    ## REMOVE CONTAMINANTS
    qfilterData <- remove_contaminants(data=maxQuant$data, pipe=pipe, enrich=enrich)
    stats[["no_contam_removed"]] <- qfilterData$no_contam

    ## LOCAL PROBABILITY FILTER
    qfilterData <- local_prob_filter(data=qfilterData$data, min.prob=min.prob,
                                     pipe=pipe, enrich=enrich)
    param[["min.prob"]] <- qfilterData$min.prob
    stats[["no_localprob_removed"]] <- qfilterData$no_localprob_removed

    ## EXTRACT PHOSPHO DATA
    extPhos <- extract_phospho_data(data=qfilterData$data, sampleIDs=sampleIDs,
                                    pipe=pipe, enrich=enrich)
    stats[["total_input_samples"]] <- extPhos$classList$class_1$no_class_samples
    stats[["no_extracted_rows"]] <- extPhos$classList$class_1$no_class_phospho
    param[["pattern"]] <- extPhos$classList$class_1$idList$pattern

    ## SAVE PARAM/STATS TO LOG FILE
    logs <- make_log(param=param, stats=stats, title="EXTRACTED DATA",save=TRUE)


    ## return a list of data.frame containing the extracted quality
    ## filtered intensity data for class 1 phosphopeptides,
    ## corresponding extracted annotation, stats, and input parameters
    ## classList corresponds to data/annotation for each phosphopeptide
    ## class e.g. class___1, class___2, class___3

    data2 <- list(data=extPhos$classList$class_1$rawData,   ## class 1 data
                  annot=extPhos$classList$class_1$rawAnnot, ## class 1 annot
                  classList=extPhos$classList, ## class_1$data/annot, class_2$data/annot etc.
                  sampleIDs=extPhos$classList$class_1$idList$sampleIDs,
                  pattern=extPhos$classList$class_1$idList$pattern,
                  param=logs$param,stats=logs$stats[c(4,1,2,3,5), ])
    return(data2)


  } ## phosphoTMT PHOSPHO


}





import_data <-function(file, pipe, enrich) {

  ## check that file exists in specified location
  ## and is csv/tsv/txt
  file <-gsub("\\./","",file);file
  file.dir <- dirname(file);file.dir
  file <- file.path(file.dir, match.arg(basename(file), list.files(file.dir),
                                        several.ok=FALSE));file
  filext <- tools::file_ext(file);filext
  filext <- match.arg(filext, c("csv","txt","tsv"),several.ok=FALSE)
  if(filext=="txt" | filext=="tsv"){ sep="\t" }
  if(filext=="csv"){ sep=","}

  ## vector of required column names based on analysis pipeline
  ## and type of enrichment. these columns are processed by other
  ## functions that extract annotation, remove contaminants, etc.
  ## values for each variable are defined according to R syntax
  ## rules. Values for each variable are listed at the top of the
  ## functions.R file.
  if(all(pipe=="DIA" & enrich=="protein")){
    reqCols<-c(diaAnnotationColums,diaContamColums)
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
    if(data[nrow(data),1]=="END OF FILE"){ data <- data[-nrow(data),] };data[nrow(data),1]

    ## name of column 1 is changed to 'id'
    if(colnames(data)[1]=="X."){ colnames(data)[1] <- "id" }
    data <- data.frame(data)
    no_input <- nrow(data)

    ## checks to make sure the input file contains the required annotation/contamination
    ## columns. This ensures that the file imported correctly and helps to identify
    ## sampleIDs columns and extract targets info. later.
    if(all(reqCols %in% colnames(data))==TRUE){

      print(paste(pipe, enrich, "file imported. Success!!"))
      print(paste("the file contains",no_input, enrich, "entries"));cat("\n")

      data2<-list(data=data, no_input=no_input, reqCols=reqCols)
      return(data2)

    } else { stop(paste0("Error! The input file may be missing one or more",
                         "required annotation/contamination \ncolumns or ",
                         "id not import properly. \n\nFor DIA experiments ",
                         "the input sample report file should include the following",
                         "column names:\n", paste(reqCols, collapse="', '"), "\n",
                         " The annotation/sample data should also begin in row 11 ",
                         " of the report."))
    }
  } ## DIA IMPORT


  ## IMPORT TMT/phosp DATA
  if(any(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF")){

    ## imports file as data.frame. Note: column names of the input file
    ## are defined according to R syntax rules. special characters, spaces,
    ## are converted to periods.
    data <- utils::read.csv(file=file, sep=sep, stringsAsFactors=FALSE,
                            header=TRUE, check.names=TRUE)
    no_input <- nrow(data)

    ## checks that file contains required annotation/contaminant columns,
    ## the id column of MaxQuant output files is a unique key, if the values
    ## are non-numeric or duplicated the file may contain some rows with trash
    ## information that needs to be removed. if this is the case a msg is returned
    ## alerting the user that the file needs manual inspection. A data.frame
    ## of the input data is returned if it contains all the required columns
    ## and the id column contains numeric unique values (0-n)
    if(all(reqCols %in% colnames(data))==TRUE){
      if(all(!duplicated(data$id) & is.integer(data$id))){

        print(paste(pipe, enrich, "file imported. Success!!"))
        print(paste("the file contains",no_input, enrich, "entries"));cat("\n")

        data2 <- list(data=data, no_input=no_input)
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



remove_contaminants <- function(data, pipe, enrich){

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

    no_rows <- nrow(data);no_rows
    data <- data[!grepl("DECOY", data[, contamColums]),];nrow(data)
    data <- data[!grepl("Group of", data[, contamColums]),];nrow(data)

    no_contam  <- no_rows - nrow(data)
    no_qfilter <- nrow(data)

    print(paste(no_contam, "contaminants removed. Success!!"))
    print(paste(no_qfilter, pipe, enrich, "entries retained ..."))

    data2 <- list(data=data, no_contam=no_contam, no_qfilter=no_qfilter)
    return(data2)

  } ## DIA CONTAM

  if(pipe=="TMT" | pipe=="phosphoTMT" | pipe=="LF"){

    no_rows <- nrow(data)
    for(x in contamColums){
      data[,x][is.na(data[,x])] <- ""
      remove<-data[,x]=="+"
      data <- data[remove==FALSE, ]
    }

    no_contam  <- no_rows - nrow(data)
    no_qfilter <- nrow(data)

    print(paste(no_contam, "contaminants removed. Success!!"))
    print(paste(no_qfilter, pipe, enrich, "entries retained ..."))

    data2 <- list(data=data, no_contam=no_contam, no_qfilter=no_qfilter)
    return(data2)

  } ## TMT/phosphoTMT/LF

}


extract_protein_data <- function(data, sampleIDs=NULL, pipe, enrich){

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
    data$Accession.Number <- gsub("(.+?)(\\ .*)","\\1",stringr::str_extract(data$Protein.Name,"(?<=\\|)[^\\|]+(?=\\ )"))
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
    replaceCommas <- function(x){ x<-as.numeric(gsub("\\,", "", x)) }
    for(i in 1:ncol(rawData)){ rawData[,i] <- replaceCommas(rawData[,i]) }

    no_samples   <- ncol(rawData)
    no_extracted <- nrow(rawData)

    print(paste("intensity data for", no_extracted,
                pipe, enrich, "entries and ",no_samples, "samples extracted.",
                "Success!!"));cat("\n")

    data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                  annotColums=annotColums, no_samples=no_samples,
                  no_extracted=no_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
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

    no_samples   <- ncol(rawData)
    no_extracted <- nrow(rawData)

    print(paste("corrected reporter intensity data for", no_extracted,
                pipe, enrich, "entries and ",no_samples, "samples extracted.",
                "Success!!"));cat("\n")

    data2 <- list(rawData=rawData, rawAnnot=rawAnnot, sampleIDs=sampleIDs,
                  annotColums=annotColums, no_samples=no_samples,
                  no_extracted=no_extracted, pipe=pipe, enrich=enrich, pattern=pattern)
    return(data2)

  } ## TMT/PROTEIN

}



make_log <- function(param, stats, title="", save=TRUE){

  param<-data.frame(t(as.data.frame(t(param))))
  colnames(param)<-"";param

  stats<-data.frame(t(as.data.frame(t(stats))))
  colnames(stats)<-"";stats

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
    colnames(param)<-"Value";param


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
    colnames(stats)<-"Value";stats
  } ## SAVE

  data2<-list(param=param, stats=stats)
  return(data2)

}



local_prob_filter <- function(data, min.prob, pipe, enrich){

  no_rows  <- nrow(data)
  min.prob <- ifelse(is.numeric(min.prob),min.prob,0.75)
  remove   <- data[, phosphoLocalProbColum] < min.prob
  data     <- data[remove==FALSE, ]

  no_localprob_removed <- no_rows - nrow(data)
  no_localprob_kept    <- nrow(data)

  print(paste(no_localprob_removed, pipe, "entries with localization",
              "probabilities < ", min.prob, "removed. Success!!"))
  print(paste(no_localprob_kept, pipe, "entries kept for further analysis..."))

  data2 <- list(data=data, min.prob=min.prob, no_localprob_removed=no_localprob_removed,
                no_localprob_kept=no_localprob_kept)
  return(data2)


}

extract_sampleIDs <- function(colNames, sampleIDs=NULL, pipe, enrich){

  ## the first 11 columns of DIA sample reports contain
  ## 2 col. of junk & 9 annotation info. sample data should
  ## begin in column 12, with pool samples used to create the
  ## library in the first 3 columns (12-14) followed by the
  ## experimental samples.
  ## NOTE: the 2 columns of the file are removed when the file
  ## is imported. Thus samples will start in column 10.
  if(pipe=="DIA" & enrich=="protein"){

    if(!is.null(sampleIDs)){

      IDs<-sampleIDs
      keep<-colNames%in%IDs;keep;table(keep)
      sampleIDs<-colNames[keep==TRUE]
      stopifnot(length(sampleIDs)==length(IDs));sampleIDs
      if(identical(sampleIDs, character(0))){
        stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
      }
    } ## NOT NULL

    if(is.null(sampleIDs)){
      ## columns containing sampleIDs should end with .mzML extension.
      keep <- grep("mzML", colNames);keep
      sampleIDs <- colNames[keep];sampleIDs

      ## if any annotation/contaminant column names are still included in list
      ## sampleIDs remove them.
      remove    <- sampleIDs %in% unique(diaAnnotationColums,diaContamColums);remove
      sampleIDs <- sampleIDs[remove==FALSE];sampleIDs

    } ## NULL


    data2 <- list(sampleIDs=sampleIDs)
    return(data2)

  } ## DIA

  ## label free looking for iBAQ.XXX columns. user should input ids
  if(pipe=="LF" & enrich=="protein"){

    if(!is.null(sampleIDs)){

      IDs<-sampleIDs
      pattern="NULL"
      keep<-colNames%in%IDs;keep;table(keep)
      sampleIDs<-colNames[keep==TRUE]
      stopifnot(length(sampleIDs)==length(IDs));sampleIDs
      if(identical(sampleIDs, character(0))){
        stop(paste0("\nsampleIDs could not be identified based on input\n
                           sampleIDs='", paste(IDs,collapse=","),"\nsampleIDs==character(0)"))
      }
    } ## NOT NULL

    if(is.null(sampleIDs)){

      pattern="iBAQ."
      sampleIDs <- grep(pattern, colNames, value=TRUE)
      sampleIDs <- sampleIDs[-grep("iBAQ.peptides",sampleIDs)];sampleIDs

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

        keep<-colNames%in%IDs;keep;table(keep)
        sampleIDs<-colNames[keep==TRUE]
        stopifnot(length(sampleIDs)==length(IDs));sampleIDs
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
          sampleIDs = grep(pattern, colNames, value=TRUE);sampleIDs
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
        keep<-colNames%in%IDs;keep;table(keep)

        classIDs  <- colNames[keep==TRUE]
        stopifnot(length(classIDs)==length(IDs));classIDs
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
          classIDs  <- grep(pattern, colNames, value = TRUE);classIDs
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
    data$Flanking    <- gsub("\\;.*$", "", data$Sequence.window) %>% str_sub(9,23) %>% paste0("-p")

    ## get protein group ID (id col. proteinGroups) and site position for
    ## the first protein entry i.e. corresponds to Majority.protein.ID
    ## to get the first entry remove everything after the first semicolon.
    data$proGroupID     <- sub(";.*", "", data$Protein.group.IDs)
    data$phosPosition   <- sub(";.*", "", data$Positions.within.proteins)
    data$phosAAPosition <- paste(data$Amino.acid, data$phosPosition,sep="")

    annotColums <- c(annotColums, "proGroupID", "UniprotID", "Gene_name",
                     "Description", "Flanking", "phosPosition", "phosAAPosition");annotColums
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
    classNums <- IDs$classNums;classNums
    pattern   <- IDs$pattern

    ## extract annotation and corrected intensity data
    ## for single, double, tripple phosphopeptides
    rawAnnot <- data[, annotColums]; dim(rawAnnot)
    rawData  <- data[, classIDs]; dim(rawData)

    no_samples   <- ncol(rawData)
    no_extracted <- nrow(rawData)

    print(paste("intensity data for", no_extracted,pipe, enrich,
                "entries and", no_samples, "samples extracted. Success!!"))
    cat("\n")

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

      no_class_samples <- ncol(classData)
      no_class_phospho <- nrow(classData)

      print(paste("classData",i,":","intensity data for",no_class_phospho,
                  "class", i, pipe,enrich, "entries and", no_class_samples,
                  "samples extracted. Success!!",sep=" "))

      ## id list from extract_sampleIDs() for class i
      idList <- list(classIDs=classIDs, sampleIDs=sampleIDs, classNums=classNums,
                     annotColums=annotColums2, pattern=pattern)
      classList[[i]] <- list(rawData=classData, rawAnnot=classAnnot[rownames(classData), ],
                             idList=idList, no_class_samples=no_class_samples,
                             no_class_phospho=no_class_phospho)
    }## INDIV CLASSES

    data2 <- list(classList=classList)
    return(data2)

  }

}

