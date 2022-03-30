
#' Do something useless
#'
#' A function that doesn't do anything (besides serving as an example about how to test and document functions).
#'
#' @param ... Any argument you want, named or unnamed.
#'
#' @return invisible(NULL)
#' @export
#'
#' @examples
#' x <- useless_function(1234)
#'
useless_function <- function(...) {
  message("This function doesn't do anything")
  return(invisible(NULL))
}
