
process_data <- function(data, targets, group="group", min.reps=3, min.grps=1){

  ## conditional filter applied to the data
  filt <- filter_data(data=data, targets=targets, group=group,
                      min.reps=min.reps, min.grps=min.grps)

  ## apply 8 norm. methods to the filtered data set.
  ## output returned includes list object containing df for each norm. tech.
  ## filtered unnormalized dataset (data) is also returned
  norm <- normalize_data(data=filt$data, targets=filt$targets)

  param<-filt$param;colnames(param)<-NULL
  stats<-filt$stats;colnames(stats)<-NULL
  logs<-make_log(param=as.list(unlist(param)), stats=as.list(unlist(stats)),
                 title="DATA PROCESSING",save=TRUE)

  print("data processing complete... returning filtered and normalized intensity data...")
  data2<-list(normList=norm$normList, targets=norm$targets, filt=filt,
              param=logs$param,stats=logs$stats)

  return(data2)


}


filter_data <- function(data, targets, group="group", min.reps=3, min.grps=1){

  param <- stats <- list()

  stats[["total_input_samples"]] <- ncol(data)

  ## reorder/subset data matrix to match the row order in targets
  if(any(rownames(targets)%in%colnames(data))==FALSE){
    stop("Erorr! some row names of targets do not match any column names in data.")
  }
  data <- data[, rownames(targets)]
  stopifnot(rownames(targets)==colnames(data))

  stats[["no_removed_samples"]] <- stats[["total_input_samples"]] - ncol(data)
  stats[["no_filtered_samples"]] <- ncol(data)
  stats[["total_input_rows"]]=nrow(data)

  print(paste(stats[["no_removed_samples"]], "samples were removed from the data matrix..."))

  ## change data column names and targets row names
  ## to sample name i.e. sample column in targets
  if("sample" %in% colnames(targets)){
    stopifnot(rownames(targets)==colnames(data))
    print("column names of data and row names of targets converted to targets sample names...");cat("\n")
    rownames(targets) <- targets$sample
    colnames(data)    <- rownames(targets)
    stopifnot(colnames(data)==rownames(targets))
  }

  ## check that the group name provided is a valid column name in targets.
  ## extract group column as character vector and make a factor.
  if(group%in%colnames(targets)==FALSE){
    stop("Error! group should be a column name in targets (e.g. group='group')") }
  groups <- make_factor(as.character(targets[, group]))

  param[["group"]] <- group
  param[["groups"]] <- paste(paste0(names(table(groups)),"=",table(groups)),collapse=", ")

  nreps <- table(groups); nreps              ##  number of samples in each group
  ngrps <- length(unique(groups));ngrps      ##  number of groups
  # print(group.names);print(nreps);print(ngrps)



  ## if min.reps is NULL set min.reps to the smallest group size
  ## if min.grps is NULL set the min.grps to 1
  if(is.null(min.reps)){ min.reps <- min(nreps);min.reps }
  if(is.null(min.grps)){ min.grps <- 1;  min.grps }

  ## if input no min.reps exceeds the max no. reps / group
  ## min.reps value is lowered to the smallest samle group.
  repsCutoff <- min.reps <= nreps;repsCutoff
  if(all(repsCutoff)==FALSE){
    warning(paste0("The min.reps threshold min.reps = ", min.reps, " exceeds the max. ",
                   "number of replicates per group: ", paste(names(nreps),nreps,sep="=",collapse=", "),
                   ". \nThe min.reps ",
                   "threshold was lowered to equal the smallest sample group.\n",
                   "min.reps = ", min(nreps)))
    min.reps <- min(nreps); min.reps
    repsCutoff <- min.reps <= nreps;repsCutoff
  }

  ## repsCutoff is used to determine how many groups meet the min.reps criteria.
  ## if min.grps parameter is greater than grpsCutoff, then the threshold is lowered
  ## to one (i.e. )
  grpsCutoff <- sum(as.numeric(repsCutoff));grpsCutoff
  if(min.grps > grpsCutoff){
    message(paste0("Based on the min.reps parameter ",min.reps,". The min.grps threshold : ",
                   min.grps," exceeds what is allowed by the dataset: ", grpsCutoff,".\n",
                   "The min.grps threshold has been lowered to one ",
                   "group. \nmin.grps = ", 1))
    min.grps <- 1
  }

  param[["min.reps"]] <- min.reps
  param[["min.grps"]] <- min.grps
  print(paste0("extracting entries with intensity > 0 in at least ", min.reps,
               " of the samples in ", min.grps," or more groups..."))

  ## FILTERING
  ## zeros are replaced with NA
  tmpData <- data; dim(tmpData)
  tmpData[,][tmpData[,] == 0] <- NA; head(tmpData)

  ## calc. no samples in each group with intensities > 0
  noSamplesPerGroup <- NULL
  for(x in unique(groups)){
    keep<-groups %in% x;table(keep)
    grpData<- data.frame(tmpData[,keep]);head(grpData)
    ## no samples in each group with intensities > 0
    noSamplesPerGroup <- cbind(noSamplesPerGroup, apply(grpData, 1, FUN = function(x){sum(!is.na(x))}))
  }
  noSamplesPerGroup<-data.frame(noSamplesPerGroup)
  colnames(noSamplesPerGroup) <- unique(groups)
  rownames(noSamplesPerGroup) <- rownames(grpData)
  head(noSamplesPerGroup);dim(noSamplesPerGroup)

  ## min X samples in at least Y groups
  aboveCutoff <- apply(noSamplesPerGroup, 1 , FUN=function(x){ sum(x >= min.reps) >= min.grps });table(aboveCutoff)
  noSamplesPerGroup$aboveCutoff <-aboveCutoff;head(noSamplesPerGroup)
  filterData <- tmpData[aboveCutoff == TRUE, ]; dim(filterData)
  removeData <- tmpData[aboveCutoff == FALSE, ];dim(removeData)

  stats[["no_removed_rows"]]=nrow(removeData)
  stats[["no_filtered_rows"]]=nrow(filterData)

  ## filtering processing stats
  print(paste("A total of",nrow(removeData), "entries were removed from the data set leaving",
              nrow(filterData),"entries for further analysis. Success!!"))

  logs<-make_log(param,stats,title="FILTERING X in Y",save=TRUE)

  data2 <- list(data=filterData, targets=targets, noSamplesPerGroup=noSamplesPerGroup,
                rm.data=removeData, param=logs$param, stats=logs$stats)
  return(data2)


}


make_factor <- function(x, prefix=NULL){
  if(is.numeric(x)){
    prefix <- ifelse(is.null(prefix), "X", prefix)
    x <- paste(prefix, x, sep="")
  }
  x <- factor(x, levels=ordered(unique(x)))
  return(x)
}




