subset_targets <- function(targets, factor, rm.vals){

  lst <- c(); rm <- 0;
  param <- stats <- list()

  if(any(is.null(factor) | length(factor)!=1)){
    stop(paste("Error! factor is either NULL or contains multiple values",
               "filtering cannot be performed.","'factor' should refer to",
               " 1 column name in targets."))
  }

  no_input <- nrow(targets)
  stats[[factor]] <- lst
  stats[["total_input_samples"]] <- no_input

  ## REMOVE ROWS
  if(!is.null(factor)){
    for(x in rm.vals){
      no_rows<-nrow(targets);no_rows
      if(any(targets[,factor] %in% x)){
        remove <- targets[,factor] %in% x
        targets <- targets[remove==FALSE,]
        lst <- c(lst,x)
        rm_rows <- no_rows-nrow(targets)
        stats[[paste0(factor," = ",x)]] <- rm_rows
        rm<-rm + rm_rows
        message(paste0("[",factor," == ",x,"]", " : ", rm_rows," samples removed",
                       " : ", nrow(targets), " samples kept"))
      }}

  } ## RM ROWS

  stopifnot(no_input == rm + nrow(targets))
  param[["factor"]] <-factor
  param[["rm.vals"]] <- paste(lst,collapse=", ")
  stats[["total_samples_removed"]] <- rm
  stats[["total_samples_kept"]] <- nrow(targets)

  title<-paste0("SUBSET TARGETS STATS (",factor,")");title
  logs <-make_log(param=param, stats=stats,title=title,save=TRUE)


  message("metadata for ",rm," of the ",no_input," samples were removed from the targets file. Success!!\n")
  data2 <- list(targets=targets, stats=logs$stats,param=logs$param)
  return(data2)


}
