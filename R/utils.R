

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
