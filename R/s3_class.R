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
#'     number of rows as the data slot and in the same order. This data frame
#'     must contain a column named "uniprot_id", which contains a unique value
#'     for each row/protein (ideally, the UniProtID). The values in this column
#'     will be used to uniquely identify proteins in the DAList object.
#'     Additional columns containing other information (gene names,
#'     descriptions, molecular weights, etc) can be supplied.
#'   \item metadata- A data frame containing information on sample metadata,
#'     one row per sample. It should contain the same number of rows (samples) as
#'     there are columns in the data data frame, in the same order, with row names
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
#'   Must contain a column named "uniprot_id" containing unique entries for each
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

  # Check if data, annotation, and metadata
  # are tibble, convert to data frames if so.
  if (any(c("tbl_df", "tbl") %in% class(data))) {
    cli::cli_warn(
      "Input data is a tibble, converting to data frame"
    )
    data <- as.data.frame(data)
  }

  if (any(c("tbl_df", "tbl") %in% class(annotation))) {
    cli::cli_warn(
      "Input annotation is a tibble, converting to data frame"
    )
    annotation <- as.data.frame(annotation)
  }

  if (any(c("tbl_df", "tbl") %in% class(metadata))) {
    cli::cli_warn(
      "Input metadata is a tibble, converting to data frame"
    )
    metadata <- as.data.frame(metadata)
  }

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

  # Validate and return the object
  out <- validate_DAList(out)
  return(out)
}



#' DAList validator
#'
#' Internal function for validating a DAList object.
#' Tolerant of arbitrary result labels (custom contrast names).
#'
#' @keywords internal
validate_DAList <- function(x) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  required_slots <- c("data", "annotation", "metadata", "design", "eBayes_fit", "results", "tags")
  optional_slots <- c("filtered_proteins_per_contrast", "data_per_contrast", "annotation_per_contrast")
  allowed_slots  <- c(required_slots, optional_slots)
  
  missing_required <- setdiff(required_slots, names(x))
  if (length(missing_required) > 0) {
    cli::cli_abort("The DAList is missing the following required slot{?s}: {missing_required}")
  }
  
  # keep canonical order; allow extra slots
  final_order <- c(required_slots, intersect(optional_slots, names(x)))
  x <- x[final_order]
  class(x) <- "DAList"
  
  # ---- Basic shape checks ----
  if (any(c("tbl_df", "tbl") %in% class(x$data)))
    cli::cli_abort("Data slot is a tibble, must be a matrix or data frame")
  if (any(c("tbl_df", "tbl") %in% class(x$annotation)))
    cli::cli_abort("Annotation slot is a tibble, must be a matrix or data frame")
  if (any(c("tbl_df", "tbl") %in% class(x$metadata)))
    cli::cli_abort("Metadata slot is a tibble, must be a matrix or data frame")
  
  if (!any(c(is.data.frame(x$data), is.matrix(x$data))))
    cli::cli_abort("The 'data' slot must be a data frame or matrix")
  if (!all(apply(x$data, 2, is.numeric)))
    cli::cli_abort("The 'data' slot must contain only numeric data")
  
  if (!"uniprot_id" %in% colnames(x$annotation))
    cli::cli_abort("The annotation data must contain a column named 'uniprot_id'")
  if (any(duplicated(x$annotation$uniprot_id)))
    cli::cli_abort("The 'uniprot_id' values in annotation must be unique")
  if (nrow(x$data) != nrow(x$annotation))
    cli::cli_abort("'data' and 'annotation' must have the same number of rows")
  if (!all(rownames(x$data) == rownames(x$annotation)))
    cli::cli_abort("Rownames for 'data' and 'annotation' must match")
  if (!all(rownames(x$data) == x$annotation$uniprot_id))
    cli::cli_abort("Rownames for 'data'/'annotation' must equal the 'uniprot_id' column")
  
  if (nrow(x$metadata) != ncol(x$data))
    cli::cli_abort("The number of samples in 'metadata' must match columns in 'data'")
  if (any(colnames(x$data) != rownames(x$metadata)))
    cli::cli_abort("Row names of 'metadata' must match column names of 'data'")
  
  # ---- Design checks ----
  if (!is.null(x$design)) {
    need <- c("design_matrix", "design_formula")
    if (length(x$design) < 2)
      cli::cli_abort("'design' must be a list of at least 2 items: design_formula and design_matrix")
    if (!all(need %in% names(x$design)))
      cli::cli_abort("'design' must include 'design_matrix' and 'design_formula'")
    
    if (!"assign" %in% names(attributes(x$design$design_matrix)))
      cli::cli_abort("'design_matrix' must have an 'assign' attribute")
    
    if (nrow(x$design$design_matrix) != nrow(x$metadata))
      cli::cli_abort("Rows in design_matrix must match metadata")
    if (any(rownames(x$design$design_matrix) != rownames(x$metadata)))
      cli::cli_abort("Row names of design_matrix must match row names of metadata")
    
    if (!is.null(x$design$random_factor)) {
      if (!x$design$random_factor %in% colnames(x$metadata))
        cli::cli_abort("Random factor not found in metadata")
    }
    if (sum(c("contrast_matrix", "contrast_vector") %in% names(x$design)) == 1) {
      cli::cli_abort("If contrast information is included, both 'contrast_matrix' and 'contrast_vector' must be present")
    }
    if ("contrast_matrix" %in% names(x$design)) {
      if (!all(rownames(x$design$contrast_matrix) == colnames(x$design$design_matrix)))
        cli::cli_abort("Rows in contrast_matrix must match columns in design_matrix")
    }
  }
  
  # ---- eBayes fit checks ----
  if (!is.null(x$eBayes_fit)) {
    if (is.list(x$eBayes_fit) && !"MArrayLM" %in% class(x$eBayes_fit)) {
      all_marrays <- all(vapply(x$eBayes_fit, function(obj) "MArrayLM" %in% class(obj), logical(1)))
      if (!all_marrays) cli::cli_abort("All elements of eBayes_fit list must be of class 'MArrayLM'")
    } else if (!is.null(x$eBayes_fit) && !"MArrayLM" %in% class(x$eBayes_fit)) {
      cli::cli_abort("eBayes_fit must be of class 'MArrayLM' or a list of them")
    }
  }
  
  # ---- Results checks (lenient on labels) ----
  if (!is.null(x$results)) {
    term_names <- if (!is.null(x$design$contrast_matrix)) {
      colnames(x$design$contrast_matrix)
    } else {
      colnames(x$design$design_matrix)
    }
    
    expected <- if (!is.null(x$tags$extract_intercept) && isTRUE(x$tags$extract_intercept)) {
      term_names
    } else {
      term_names[!grepl("ntercept", term_names, ignore.case = TRUE)]
    }
    
    res_names <- names(x$results)
    
    if (!identical(res_names, expected)) {
      unmapped <- setdiff(res_names, expected)
      if (length(unmapped) > 0) {
        dcols <- colnames(x$design$design_matrix)
        still_bad <- character(0)
        for (lab in unmapped) {
          ci <- tryCatch(x$tags$per_contrast[[lab]]$contrast_info, error = function(e) NULL)
          ok <- FALSE
          if (!is.null(ci) && !is.null(ci$contrast_expression)) {
            expr <- as.character(ci$contrast_expression)[1]
            tokens <- unique(unlist(regmatches(expr, gregexpr("[A-Za-z0-9_.:]+", expr))))
            if (length(intersect(tokens, dcols)) > 0) ok <- TRUE
          }
          if (!ok) still_bad <- c(still_bad, lab)
        }
        if (length(still_bad) > 0) {
          cli::cli_warn("Results names do not match design/contrast terms and lack usable tags: {still_bad}. Proceeding.")
        }
      } else {
        cli::cli_warn("Results names differ from expected design/contrast terms; proceeding with arbitrary labels.")
      }
    }
    
    # Per-contrast row checks: must be subset of data rows (order may differ)
    for (term in res_names) {
      result_rows <- rownames(x$results[[term]])
      data_rows   <- rownames(x$data)
      
      if (!all(result_rows %in% data_rows)) {
        cli::cli_abort("Some rows in results for contrast '{term}' are not present in 'data'")
      }
      
      if ("filtered_proteins_per_contrast" %in% names(x) && !is.null(x$filtered_proteins_per_contrast[[term]])) {
        expected_subset <- x$filtered_proteins_per_contrast[[term]]
        if (!all(result_rows %in% expected_subset)) {
          cli::cli_warn("Some rows in results for contrast '{term}' do not match the filtered protein list")
        }
      }
    }
  }
  
  if ("filtered_proteins_per_contrast" %in% names(x)) {
    if (!is.list(x$filtered_proteins_per_contrast))
      cli::cli_abort("'filtered_proteins_per_contrast' must be a list")
  }
  
  x
}
