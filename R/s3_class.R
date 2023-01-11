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
#' A function to create a DAList from existing data. Need to add in notes and further
#' documentation as we decide on exactly what requirements we have.
#' In particular, should note that: data and annotation need to have the same number of
#' rows and should already be in order. Columns in data should be same as number of rows in metadata.
#' TODO: need to expand on documentation as we finalize things. Also need to decide if we
#' should include all the other slots (design through tags). Also need to add examples.
#'
#' @param data A dataframe or matrix containing protein intensity data for each sample. Rows are proteins, columns are samples.
#' @param annotation A dataframe containing protein annotation data. Must contain a column named uniprot_id containing unique entries for each row.
#' @param metadata A dataframe containing metadata on each sample. REQUIRED COLUMNS?
#' @param design Optional, possibly unneeded. But a list of design information.
#' @param eBayes_fit Optional, possibly unneeded. An ebayes fit/MAarray object from limma ebayes
#' @param results Optional, possibly unneeded. list of table of stats.
#' @param tags Optional, possibly unneeded. List of tags.
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


  # Check for a uniprot_id column in annotation
  if ("uniprot_id" %notin% colnames(annotation)) {
    cli::cli_abort("The annotation data must contain a column named \"uniprot_id\"")
  }

  # Check that the uniprot_id column contains only unique info
  if (any(duplicated(annotation$uniprot_id))) {
    cli::cli_abort(c("Entries in the {.val uniprot_id} column of the annotation are not unique.",
                     "i" = "Find duplicate items with {.fun base::duplicated} or {.fun base::anyDuplicated} and make them unqiue"))
  }

  # Assign uniprot_id column as rownames for data and annotation
  rownames(data) <- annotation$uniprot_id
  rownames(annotation) <- annotation$uniprot_id


  # Use internal constructor to make into DAList
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


  # Annotation-specific checks
  # Check for a uniprot_id column in annotation
  if ("uniprot_id" %notin% colnames(x$annotation)) {
    cli::cli_abort("The annotation data must contain a column named \"uniprot_id\"")
  }

  # Check that the uniprot_id column contains only unique info
  if (any(duplicated(x$annotation$uniprot_id))) {
    cli::cli_abort(c("Entries in the {.val uniprot_id} column of the annotation are not unique.",
                     "i" = "Find duplicate items with {.fun base::duplicated} or {.fun base::anyDuplicated} and make them unqiue"))
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

  # rownames of data and annotation should equal the uniprot_id column in the annotation
  if (!(all(rownames(x$data) == x$annotation$uniprot_id))) {
    cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DAList must the {.val uniprot_id} column of the annotation.")
  }

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





