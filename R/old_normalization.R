

normalize_data <- function(data, targets){

  if(ncol(data) > nrow(targets)){
    warning("Warning! data < targets. the data set will be subsetted to match",
            "the samples listed in the targets file prior to normalization...")
  }
  if(any(rownames(targets)%in%colnames(data))==FALSE){
    stop("some row names in targets do not match the column names in data...")
  }
  data <- data[,rownames(targets)]
  stopifnot(colnames(data)==rownames(targets))

  ## change target row names to sample name
  ## change column names of data to target row names
  ## i.e. sample name
  if("sample"%in%colnames(targets)){
    print("column names of data and row names of targets converted to targets sample names...");cat("\n")
    rownames(targets) <- targets$sample
    colnames(data)<-rownames(targets)
    stopifnot(colnames(data)==rownames(targets))
  }

  ## create named list object to hold output dataframes of normalized data.
  ## names= normalization methods
  normList <- vector("list", length(norm.methods))
  names(normList) <- norm.methods

  ## apply normalization methods using functions listed above.
  ## NOTE: most of the normalizations use log2(intensity) as input except
  ## VSN which normalizes using raw intensities with no log2 transformation.
  normList[["log2"]]     <- logNorm(dat=data)
  normList[["median"]]   <- medianNorm(logDat=normList[["log2"]])
  normList[["mean"]]     <- meanNorm(logDat=normList[["log2"]])
  normList[["vsn"]]      <- vsnNorm(dat=data)
  normList[["quantile"]] <- quantNorm(logDat=normList[["log2"]])
  normList[["cycloess"]] <- cycLoessNorm(logDat=normList[["log2"]])
  normList[["rlr"]]      <- rlrNorm(logDat=normList[["log2"]])
  normList[["gi"]]       <- giNorm(logDat=normList[["log2"]])

  data2 <- list(normList=normList, targets=targets)
  print(paste("Data normalized... Success!!"))
  return(data2)

}



logNorm <- function(dat) {
  logInt <- log2(dat)
  #logInt <- replace(is.infinite(logInt), NA)
  logInt[is.infinite(as.matrix(logInt))] <- NA
  return(as.matrix(logInt))
}

medianNorm <- function(logDat) {
  # Find medians of each sample
  # Divide by median
  # Multiply by mean of medians
  sampleMed <- apply(logDat, 2, median, na.rm=TRUE)
  meanMed <- mean(sampleMed, na.rm=TRUE)
  out <- t(t(logDat) / sampleMed)
  out <- out * meanMed
  return(as.matrix(out))
}

meanNorm <- function(logDat) {
  # Find means of each sample
  # Divide by mean
  # Multiply by mean of means
  sampleMean <- apply(logDat, 2, mean, na.rm=TRUE)
  meanMean <- mean(sampleMean, na.rm=TRUE)
  out <- t(t(logDat) / sampleMean)
  out <- out * meanMean
  return(as.matrix(out))
}

vsnNorm <- function(dat) {
  vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
  colnames(vsnNormed) <- colnames(dat)
  row.names(vsnNormed) <- rownames(dat)
  return(as.matrix(vsnNormed))
}

quantNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy=FALSE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

cycLoessNorm <- function(logDat) {
  cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method="fast")
  colnames(cycLoessNormed) <- colnames(logDat)
  row.names(cycLoessNormed) <- rownames(logDat)
  return(as.matrix(cycLoessNormed))
}

rlrNorm <- function(logDat) {
  rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(logDat), noLogTransform=TRUE)
  colnames(rlrNormed) <- colnames(logDat)
  row.names(rlrNormed) <- rownames(logDat)
  return(as.matrix(rlrNormed))
}

giNorm <- function(logDat) {
  giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform=TRUE)
  colnames(giNormed) <- colnames(logDat)
  row.names(giNormed) <- rownames(logDat)
  return(as.matrix(giNormed))
}

