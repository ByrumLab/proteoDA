#' DAList internal constructor
#'
#' Internal function for constructing a DAList object with no/minimal checks.
#'
#' @param x An object to be converted into our DAList type
#'
#' @return A DAList object derived from x
#'
#'
new_DAList <- function(x = list()) {
  stopifnot(is.list(x))

  structure(x, class = "DAList")
}


#' Create a DAList
#'
#' A function to create a DAList from existing data.
#'
#' @param data
#' @param annotation
#' @param metadata
#' @param design
#' @param eBayes_fit
#' @param results
#' @param tags
#'
#' @return A DAList object
#'
#' @export
#'
#'
DAList <- function(data,
                   annotation,
                   metadata,
                   design = NULL,
                   eBayes_fit = NULL,
                   results = NULL,
                   tags = NULL) {

  # Do some testing or checking?
  # In particular, try to re-order elements, rename data columns
  # or sample IDs, anything like that?


  # User internal constructor to make into DAList
  out <- new_DAList(
    x = list(data = data,
             annotation = annotation,
             metadata = metadata,
             design = design,
             eBayes_fit = eBayes_fit,
             results = results,
             tags = tags)
    )


  # Check validity
  validate_DAList(out)

}



#' DAList validator
#'
#' Internal function for validating a DAList object
#'
#' @param x An object to be tested if it is a valid DAList.
#'
#' @return If x was a valid DAList, returns x
#'
#' @keywords internal
#'
validate_DAList <- function(x) {

  ## Check overall structure

  # Should be 7 elements with the correct names
  slots <- c("data", "annotation", "metadata", "design", "eBayes_fit", "results", "tags")

  # Check for proper number and order
  if (!identical(names(x), slots)) {
    # Possible problems: missing slots, extra slots, or wrong order.
    # Technically, may have multiple of these problems, but users can solve one at a time

    # May be missing slots
    missing_slots <- slots[slots %notin% names(x)]
    if (length(missing_slots) > 0) {
      cli::cli_abort("The DAList is missing the following slot{?s}: {missing_slots}")
    }
    # May have extra slots
    extra_slots <- names(x)[names(x) %notin% slots]
    if (length(extra_slots) > 0) {
      cli::cli_abort("The DAList contains the following extra slot{?s}: {extra_slots}")
    }

    # Otherwise, out of order.
    # Reorder with warning
    cli::cli_alert("Slots in DAList out of order. Reordering.")
    x <- x[slots]
  }

  ## Check each element

  # Data
  if (!any(c(is.data.frame(x$data), is.matrix(x$data)))) {
    cli::cli_abort("The {.arg data} slot of a DAList must be a dataframe or matrix")
  }

  # Data must be numeric
  if (!all(apply(x$data, 2, is.numeric))) {
    cli::cli_abort("The {.arg data} slot of a DAList must contain only numeric data")
  }

  # Data and annotation should match
  # data and annotation slots should have same number of rows
  if (nrow(x$data) != nrow(x$annotation)) {
    cli::cli_abort("The {.arg data} slot and {.arg annotation} slots of a DAList must have the same number of rows")
  }

  # Data and annotation should have matching rownames
  if (!(all(rownames(x$data) == rownames(x$annotation)))) {
    cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DAList must match")
  }

  # TODO Annotation-specific checks??


  # Metadata checks
  # TODO add more?
  if (!is.null(x$metadata)) {
    if (nrow(x$metadata) != ncol(x$data)) {
      cli::cli_abort("The number of samples in the metadata ({nrow(x$metadata)}) do not match the number of samples in the data ({ncol(x$data)})")
    }
    if (any(colnames(x$data) != rownames(x$metadata))) {
      cli::cli_abort("The row names of the metadata do not match the column names of the data")
    }
  }

  # TODO Checks for design:
  # If design matrix exists, needs to have same # of rows as metadata
  # also, rownames of design matrix need to equal colnames of data

  # If there's a random effect blocking factor, make sure its in the metadata?
  # if you have a formula, must have a matrix and vice versa?
  # And, if you have a random effect, must have a matrix and formula?


  # TODO Checks for contrasts
  # If you have a contrasts, must have a design
  # if you have a contrasts_vec, must have contrasts_matrix and vice versa?

  # TODO Checks for eBayes fit
  # if not null, then the first item should be an MArrayLM. If random effect,
  # should have a non-null correlation term. If no random effect,
  # should have a null correlation term.


  # TODO Checks for results
  if (!is.null(x$results)) {

    # if not null, pval_thresh, lfc_thresh, and adj_method should be set in tags??
    # if not null, length should match either ncol(design_matrix) or ncol(contrasts_matrix)
    # if not null, nrow and rownames for each element of the results should
    # match the data (which matches the annotation, as we check above. )


    # # data and annotation slots should have same number of rows
    # if (nrow(x$data) != nrow(x$annotation)) {
    #   cli::cli_abort("The {.arg data} slot and {.arg annotation} slots of a DAList must have the same number of rows")
    # }
    #
    # # Data and annotation should have matching rownames
    # if (!(all(rownames(x$data) == rownames(x$annotation)))) {
    #   cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DAList must match")
    # }

  }

  # TODO Any tags checks?
  # Not sure...

  # If all checks pass, return input
  x
}





