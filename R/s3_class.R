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

  # data and annotation slots should have same number of rows
  if (nrow(x$data) != nrow(x$annotation)) {
    cli::cli_abort("The {.arg data} slot and {.arg annotation} slots of a DIAlist must have the same number of rows")
  }

  # Data and annotation should have matching rownames
  if (!(all(rownames(x$data) == rownames(x$annotation)))) {
    cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DIAlist must match")
  }


  # CHECKS FOR TARGETS/METADATA
  # Maybe need more??
  if (!is.null(x$metadata)) {
    if (nrow(x$metadata) != ncol(x$data)) {
      cli::cli_abort("The number of samples in the metadata ({nrow(x$metadata)}) do not match the number of samples in the data ({ncol(x$data)})")
    }
    if (any(colnames(x$data) != rownames(x$metadata))) {
      cli::cli_abort("The row names of the metadata do not match the column names of the data")
    }
  }

  # Checks for design:
  # If design matrix exists, needs to have same # of rows as metadata
  # also, rownames of design matrix need to equal colnames of data

  # If there's a random effect blocking factor, make sure its in the metadata?
  # if you have a formula, must have a matrix and vice versa?
  # And, if you have a random effect, must have a matrix and formula?


  # Checks for contrasts
  # If you have a contrasts, must have a design
  # if you have a contrasts_vec, must have contrasts_matrix and vice versa?

  # Checks for eBayes fit
  # if not null, then the first item should be an MArrayLM. If random effect,
  # should have a non-null correlation term. If no random effect,
  # should have a null correlation term.


  # Checks for results
  if (!is.null(x$results)) {

    # if not null, pval.thresh, lfc.thresh, and adj.method should be set in tags??
    # if not null, length should match either ncol(design_matrix) or ncol(contrasts_matrix)
    # if not null, nrow and rownames for each element of the results should
    # match the data (which matches the annotation, as we check above. )


    # # data and annotation slots should have same number of rows
    # if (nrow(x$data) != nrow(x$annotation)) {
    #   cli::cli_abort("The {.arg data} slot and {.arg annotation} slots of a DIAlist must have the same number of rows")
    # }
    #
    # # Data and annotation should have matching rownames
    # if (!(all(rownames(x$data) == rownames(x$annotation)))) {
    #   cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DIAlist must match")
    # }

  }







  # If all checks pass, return input
  x
}





