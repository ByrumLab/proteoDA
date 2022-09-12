#' Make a log file
#'
#' Function used to write out logs across various top-level functions
#'
#' @param param The param object within the \code{\link{read_DIA_data}} function
#'   body. Slightly unclear what the requirements for this are.
#' @param stats The stats object within the \code{\link{read_DIA_data}} function
#'   body. Slightly unclear what the requirements for this are.
#' @param title A title for the printed parameter block. Default is "".
#' @param save TRUE/FALSE: should the parameters be save dot a text file?
#'
#' @return A list, with the param and stats objects. If SAVE = TRUE, has side
#'   effect of creating log files.
#'
#' @export
#'
#' @examples
#' # No examples yet

make_log <- function(param, stats, title="", save=TRUE){

  param<-data.frame(t(as.data.frame(t(param))))
  colnames(param)<-""

  stats<-data.frame(t(as.data.frame(t(stats))))
  colnames(stats)<-""

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
    colnames(param)<-"Value"


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
    colnames(stats)<-"Value"
  } ## SAVE

  data2<-list(param=param, stats=stats)
  return(data2)

}
