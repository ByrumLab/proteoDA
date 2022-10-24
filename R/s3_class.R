#' DIAlist internal constructor
#'
#' Internal function for constructing a DIAlist object with no/minimal checks.
#'
#' @param x An object to be converted into our DIAlist type
#'
#' @return A prototype of our new S3 list type (currently just a list).
#'
#' @examples
#' # No examples yet
#'
new_DIAlist <- function(x = list()) {
  stopifnot(is.list(x))

  structure(x, class = "DIAlist")
}


#' DIAlist internal constructor
#'
#' Internal function for constructing a DIAlist object with no/minimal checks.
#'
#' @param x An object to be converted into our DIAlist type
#'
#' @return A prototype of our new S3 list type (currently just a list).
#'
#' @examples
#' # No examples yet
#'
validate_DIAlist <- function(x) {

  # This will get complicated
  # TODO: maybe separate out into separate functions to check each slot??
  # Kinda depends on dependencies between slots...

  # Data must be a matrix or data frame
  # TODO: decide on this. Should it just be a data frame?
  # Is matrix allowed?
  if (!any(c(is.data.frame(x$data), is.matrix(x$data)))) {
    cli::cli_abort("The {.arg data} slot of a DIAlist must be a dataframe or matrix")
  }

  # Data must be numeric
  if (!all(apply(x$data, 2, is.numeric))) {
    cli::cli_abort("The {.arg data} slot of a DIAlist must contain only numeric data")
  }


  # CHECKS FOR TARGETS/METADATA

  # If all checks pass, return input
  x
}




