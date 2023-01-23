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
#' Creates a DAList from existing raw data. DALists are S3 objects of class
#' DAList and contain 7 slots which hold data and results for a quantitative
#' proteomics project (see Details for a full description of each slot). The
#' first three slots (data, annotation, and metadata) are required to create a
#' DAList. The other slots can be supplied, but are generally added via their
#' respective analysis functions. *NB*- you must ensure that the data,
#' annotation, and metadata are in the proper order. To help with this,
#' \code{DAList} requires the row names of the metadata to match the column
#' names of the data. However, DAList cannot know if the rows in the data and
#' the annotation are in the proper order: double-check that they are!
#'
#' A DAList contains 7 slots: \enumerate{
#'
#'   \item data- A matrix or data frame containing protein intensity data.
#'     Rows are proteins, columns are samples. This array should contain numeric
#'     information only, and the number of rows (proteins) should match the
#'     number of rows in the annotation data frame. The columns (samples) in the
#'     data should be in the same order as the rows of samples in the metadata
#'     data frame, with column names that match the row names of the metadata.
#'   \item annotation- A data frame of protein information, containing the same
#'     number of rows as the data slot and in the same order. This dataframe
#'     must contain a column named "uniprot_id" which contains a unique value
#'     for each row/protein (ideally, the UniProt ID). The values in this column
#'     will be used to uniquely identify proteins in the DAList object.
#'     Additional columns containing other information (gene names,
#'     descriptions, molecular weights, etc) can be supplied.
#'   \item metadata- A data frame containing information on sample metadata,
#'     one row per sample. Should contain the same number of rows (samples) as
#'     there are columns in the data dataframe, in the same order, with rownames
#'     that match the column names of the data.
#'   \item design- A list containing information on the statistical design used
#'     for analyzing differential abundance. Set using \code{\link{add_design}}
#'     and \code{\link{add_design}}.
#'   \item eBayes_fit- A list containing the model fit object. This is added
#'     using \code{\link{fit_limma_model}}, which uses functions from the limma
#'     package to perform differential abundance analysis. Specifically, this
#'     slot should be an MAarrayLM object output from
#'     \code{\link[limma:eBayes]{limma::eBayes}}, which is used within
#'     \code{\link{fit_limma_model}}.
#'   \item results- A list of data frames, in which each data frame contains
#'     differential abundance results for a given statistical term or contrast.
#'     These tables are derived from the output of
#'     \code{\link[limma:decideTests]{limma::decideTests}} and
#'     \code{\link[limma:topTable]{limma::topTable}}.
#'   \item tags- A list containing flags, parameter values, and other settings
#'     for the analysis. Generally, should not be set or modified by hand, but
#'     will be set by various pipeline functions.
#' }
#'
#' @param data A data frame or matrix containing protein intensity data for
#'   each sample. Rows are proteins, columns are samples.
#' @param annotation A data frame containing protein annotation data.
#'   Must contain a column named uniprot_id containing unique entries for each
#'   row.
#' @param metadata A data frame containing metadata on each sample.
#' @param design Optional: a list of design information, see Details.
#' @param eBayes_fit Optional: a model fit object, see Details.
#' @param results Optional: a list of data frame(s) summarizing differential
#'   abundance statistics, see Details.
#' @param tags Optional: a list containing flags, parameter values, and other
#'   settings, see Details.
#'
#' @return A DAList object
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Prepare data
#' # In most cases, these would be imported from external
#' # files and processed.
#' # Here, an example with 4 samples and 3 proteins.
#'
#' # Raw numeric data
#' raw_data <- data.frame(sampleA = c(10000, 15000, 12500),
#'                        sampleB = c(20000, 30000, 25000),
#'                        sampleC = c(18000, 23000, 20500),
#'                        sampleD = c(36000, 46000, 41000))
#'
#' # Annotation data
#' # Must be in same order as numeric values in data
#' # and requires one column named uniprot_id containing unique values
#' protein_info <- data.frame(uniprot_id = c("A0A023ZSD2", "M9NH73", "H9AZR4"),
#'                            gene = c("Mcf", "shh", "WntA"))
#'
#' # sample metadata
#' # should be in same order as samples in the data,
#' # with row names that match the column names of the data
#' sample_info <- data.frame(sample_id = c("sampleA", "sampleB",
#'                                         "sampleC", "sampleD"),
#'                           group = c("control",   "control",
#'                                     "treatment", "treatment"))
#' rownames(sample_info) <- colnames(raw_data)
#'
#'
#' # Assemble the DAList from existing raw data.
#' raw <- DAList(data = raw_data,
#'               annotation = protein_info,
#'               metadata = sample_info)
#' }
DAList <- function(data,
                   annotation,
                   metadata,
                   design = NULL,
                   eBayes_fit = NULL,
                   results = NULL,
                   tags = NULL) {

  # Check for a uniprot_id column in annotation
  if ("uniprot_id" %notin% colnames(annotation)) {
    cli::cli_abort(
      "The annotation data must contain a column named \"uniprot_id\""
    )
  }

  # Check that the uniprot_id column contains only unique info
  if (any(duplicated(annotation$uniprot_id))) {
    cli::cli_abort(c("Entries in the {.val uniprot_id} column of the annotation are not unique.",
                     "i" = "Find duplicate items with {.fun base::duplicated} or {.fun base::anyDuplicated} and make them unqiue"))
  }

  if (nrow(data) != nrow(annotation)) {
    cli::cli_abort("The {.arg data} and {.arg annotation} data frames must have the same number of rows")
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
    cli::cli_abort("The {.arg data} slot of a DAList must be a data frame or matrix")
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
    cli::cli_abort("Rownames for the {.arg data} and {.arg annotation} slots of the DAList must equal the {.val uniprot_id} column of the annotation.")
  }

  # Metadata checks
  if (!is.null(x$metadata)) {
    if (nrow(x$metadata) != ncol(x$data)) {
      cli::cli_abort("The number of samples in the metadata ({nrow(x$metadata)}) do not match the number of samples in the data ({ncol(x$data)})")
    }
    if (any(colnames(x$data) != rownames(x$metadata))) {
      cli::cli_abort("The row names of the metadata do not match the column names of the data")
    }
  }

  # Checks for design and contrasts:
  if (!is.null(x$design)) {

    # if present, should be a list with at least two slots: design_formula and design_matrix
    if (length(x$design) < 2) {
      cli::cli_abort(c("The object in the design slot should be a list with length of at least 2:",
                       "Current length: {length(x$design)}"))
    }

    # two required items when design is present
    if (!(all(c("design_matrix", "design_formula") %in% names(x$design)))) {
      cli::cli_abort(c("The list in the design slot must contain at least a design_formula and a design_matrix"))
    }

    # design_matrix should be an R design matrix, so it should have an "assign" attribute.
    if ("assign" %notin% names(attributes(x$design$design_matrix))) {
      cli::cli_abort(c("The design_matrix in the design slot should have an \"assign\" attribute"))
    }

    # there are 5 possible names for the elements of design
    design_possible_names <- c("design_formula", "design_matrix", "random_factor", "contrast_matrix", "contrast_vector")
    if (any(names(x$design) %notin% design_possible_names)) {
      problem_names <- names(x$design)[names(x$design) %notin% design_possible_names]
      cli::cli_abort(c("The list in the design slot may contain the following elements:",
                       "{.val {design_possible_names}}",
                       "Problem {cli::qty(problem_names)} item{?s} in design list: {.val {problem_names}}"))
    }

    # Same number of rows as metadata
    if (nrow(x$design$design_matrix) != nrow(x$metadata)) {
      cli::cli_abort("The number of samples in the design matrix ({nrow(x$design$design_matrix)}) does not match the number of samples in the metadata ({nrow(x$metadata)})")
    }

    # rownames equal colnames of data
    if (any(colnames(x$data) != rownames(x$design$design_matrix))) {
      cli::cli_abort("The row names of the design matrix do not match the column names of the data")
    }

    # If there's a random effect, must be present in the metadata
    if (!is.null(x$design$random_factor)) {
      if (x$design$random_factor %notin% colnames(x$metadata)) {
        cli::cli_abort("The random factor term is not present in the metadata")
      }
    }

    # Checks for contrasts
    # If any contrast stuff is present,
    if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) > 0) {
      # Must have both the vector and the matrix
      if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) != 2) {
        cli::cli_abort(c("if contrast information is present in the design slot, must have both a contrast_vector and contrast_matrix"))
      }
      # And rows in the contrast matrix should equal cols in the design matrix
      if (!(all(rownames(x$design$contrast_matrix) == colnames(x$design$design_matrix)))) {
        cli::cli_abort(c("The rows in the contrast matrix do not match the column in the design matrix"))
      }
    }
  }

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

  # For the moment, no tag checks.
  # Can revisit if needed.

  # If all checks pass, return input
  x
}





