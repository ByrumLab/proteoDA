
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Removes commas in character representation of a number
#'
#' @param x Character string of a number, with commas in it.
#'
#' @return Integer
#'
#' @keywords internal
#'
#'
remove_commas <- function(x){ x<-as.numeric(gsub("\\,", "", x)) }


#' %notin% operator, opposite of %in%
#'
#'
#' @inheritParams base::"%in%"
#' @return A logical vector, indicating if a match was NOT located for each
#' element of of the search list.
#'
#' @keywords internal
#'
#'
#'
`%notin%` <- function(x,table) !`%in%`(x,table)



#' Convert a column to a factor
#'
#' Converts a vector/column to a factor. For numeric inputs, adds an "X"
#' in front of the number by default, or can supply your own prefix to use
#' instead.
#'
#' @param x Input vctor to be converted to a factor
#' @param prefix Optional: a prefix to put in front of numbers when
#'   converting to a factor. Default is "X"
#'
#' @return A factor vector
#'
#' @export
#'
#' @examples
#' # No examples yet
#'

# TODO: check out forcats for this? how often do we use the prefix thing?
make_factor <- function(x, prefix = "X") {
  if (is.numeric(x)) {
    x <- paste(prefix, x, sep = "")
  }
  x <- factor(x, levels = ordered(unique(x)))
  return(x)
}




make_new_filename <- function(x,dir="."){
  kk=0
  # library(xfun,quietly=TRUE)
  filelist<-c(list.files(path=dir));filelist
  ext <- tools::file_ext(x);ext
  xx=x;xx
  success=FALSE
  while(success==FALSE){
    kk=kk+1
    # xx=sub("_[^_]+$", "", xx);xx# run the function again
    xx <- paste0(gsub(paste0(".",ext),"",xx),"_",kk,".",ext);xx
    success=(file.path(dir,xx)%in%file.path(dir,filelist)==FALSE);success
    if(success==FALSE){
      xx=sub("_[^_]+$", "", xx);xx# run the function again
      xx=x
    };xx
  };xx
  if(success==TRUE){return(xx)}
}
