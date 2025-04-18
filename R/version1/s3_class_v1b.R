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
#' Internal function for validating a DAList object
#'
#' @param x An object to be tested if it is a valid DAList.
#'
#' @return If x was a valid DAList, returns x
#'
#' @keywords internal
#'
validate_DAList <- function(x) {
  
  ## Define required vs optional slots
  required_slots <- c("data", "annotation", "metadata", 
                      "design", "eBayes_fit", "results", "tags")
  optional_slots <- c("filtered_proteins_per_contrast")
  allowed_slots <- c(required_slots, optional_slots)
  
  # Check that all required slots are present
  missing_required <- setdiff(required_slots, names(x))
  if (length(missing_required) > 0) {
    cli::cli_abort("The DAList is missing the following required slot{?s}: {missing_required}")
  }
  
  # Check for any extraneous slots not in allowed_slots
  extra_slots <- setdiff(names(x), allowed_slots)
  if (length(extra_slots) > 0) {
    cli::cli_abort("The DAList contains unknown slot{?s}: {extra_slots}")
  }
  
  # Reorder the slots so the required ones come first (in the usual order),
  # followed by any optional slot (if present).
  final_order <- c(required_slots, intersect(optional_slots, names(x)))
  x <- x[final_order]
  class(x) <- "DAList"
  ## --- BEGIN EXISTING VALIDATIONS FOR THE 7 CORE SLOTS ---
  
  # Make sure data, annotation, metadata not tibbles
  if (any(c("tbl_df", "tbl") %in% class(x$data))) {
    cli::cli_abort("Data slot is a tibble, must be a matrix or data frame")
  }
  if (any(c("tbl_df", "tbl") %in% class(x$annotation))) {
    cli::cli_abort("Annotation slot is a tibble, must be a matrix or data frame")
  }
  if (any(c("tbl_df", "tbl") %in% class(x$metadata))) {
    cli::cli_abort("Metadata slot is a tibble, must be a matrix or data frame")
  }
  
  # Data must be data frame or matrix with only numeric columns
  if (!any(c(is.data.frame(x$data), is.matrix(x$data)))) {
    cli::cli_abort("The 'data' slot of a DAList must be a data frame or matrix")
  }
  if (!all(apply(x$data, 2, is.numeric))) {
    cli::cli_abort("The 'data' slot must contain only numeric data")
  }
  
  # Annotation checks
  if ("uniprot_id" %notin% colnames(x$annotation)) {
    cli::cli_abort('The annotation data must contain a column named "uniprot_id"')
  }
  if (any(duplicated(x$annotation$uniprot_id))) {
    cli::cli_abort(c(
      "Entries in the 'uniprot_id' column of the annotation are not unique.",
      "i" = "Check duplicates with duplicated(...) or anyDuplicated(...)."
    ))
  }
  if (nrow(x$data) != nrow(x$annotation)) {
    cli::cli_abort("The 'data' and 'annotation' slots must have the same number of rows")
  }
  if (!all(rownames(x$data) == rownames(x$annotation))) {
    cli::cli_abort("Rownames for 'data' and 'annotation' must match")
  }
  if (!all(rownames(x$data) == x$annotation$uniprot_id)) {
    cli::cli_abort("Rownames for 'data'/'annotation' must equal the 'uniprot_id' column.")
  }
  
  # Metadata checks
  if (nrow(x$metadata) != ncol(x$data)) {
    cli::cli_abort(
      "The number of samples in 'metadata' ({nrow(x$metadata)}) must match the columns in 'data' ({ncol(x$data)})"
    )
  }
  if (any(colnames(x$data) != rownames(x$metadata))) {
    cli::cli_abort("Row names of 'metadata' must match the column names of 'data'")
  }
  
  # Design checks
  if (!is.null(x$design)) {
    if (length(x$design) < 2) {
      cli::cli_abort(
        "The 'design' slot should be a list of length >= 2: (design_formula, design_matrix, ...)"
      )
    }
    if (!all(c("design_matrix", "design_formula") %in% names(x$design))) {
      cli::cli_abort(
        "The 'design' slot must have at least 'design_formula' and 'design_matrix' elements"
      )
    }
    if ("assign" %notin% names(attributes(x$design$design_matrix))) {
      cli::cli_abort("The design_matrix should have an 'assign' attribute")
    }
    
    design_possible_names <- c("design_formula", "design_matrix", 
                               "random_factor", "contrast_matrix", "contrast_vector")
    if (any(names(x$design) %notin% design_possible_names)) {
      problem_names <- setdiff(names(x$design), design_possible_names)
      cli::cli_abort(c(
        "The 'design' list may only contain elements:",
        "{.val {design_possible_names}}",
        "Problem elements: {.val {problem_names}}"
      ))
    }
    
    # Check design_matrix rows match metadata
    if (nrow(x$design$design_matrix) != nrow(x$metadata)) {
      cli::cli_abort(
        "Number of samples in design_matrix ({nrow(x$design$design_matrix)}) does not match metadata ({nrow(x$metadata)})"
      )
    }
    if (any(colnames(x$data) != rownames(x$design$design_matrix))) {
      cli::cli_abort("Row names of design_matrix must match the column names of 'data'")
    }
    
    # If there's a random factor, it must be in metadata
    if (!is.null(x$design$random_factor)) {
      if (x$design$random_factor %notin% colnames(x$metadata)) {
        cli::cli_abort("The random factor term is not present in the metadata")
      }
    }
    
    # If contrasts are present, must have both contrast_matrix & contrast_vector
    if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) > 0) {
      if (sum(c("contrast_vector", "contrast_matrix") %in% names(x$design)) != 2) {
        cli::cli_abort(
          "If contrast info is present in 'design', must have both contrast_vector and contrast_matrix"
        )
      }
      if (!all(rownames(x$design$contrast_matrix) == colnames(x$design$design_matrix))) {
        cli::cli_abort("Rows in contrast_matrix do not match columns in design_matrix")
      }
    }
  }
  
  # eBayes_fit checks
  # In validate_DAList(), replace the block that checks eBayes_fit with:
  
  # eBayes_fit checks
  if (!is.null(x$eBayes_fit)) {
    
    if (is.list(x$eBayes_fit) && !inherits(x$eBayes_fit, "MArrayLM")) {
      # eBayes_fit is a *list*, not a single MArrayLM
      # => Verify that *each* element is an MArrayLM
      all_marrays <- all(vapply(
        X   = x$eBayes_fit,
        FUN = function(obj) "MArrayLM" %in% class(obj),
        FUN.VALUE = logical(1)
      ))
      if (!all_marrays) {
        cli::cli_abort(
          "If 'eBayes_fit' is a list, each element must be an object of class 'MArrayLM' (from limma)."
        )
      }
    } else {
      # eBayes_fit is not NULL, not a list => must be a single MArrayLM
      if ("MArrayLM" %notin% class(x$eBayes_fit)) {
        cli::cli_abort(
          "The eBayes_fit slot must be either a single 'MArrayLM' or a list of 'MArrayLM' objects."
        )
      }
    }
    
    # Now handle correlation checks if using a random factor
    if (!is.null(x$design$random_factor)) {
      # If there's a random factor, each MArrayLM (or the single MArrayLM) must have $correlation
      if (is.list(x$eBayes_fit)) {
        for (obj in x$eBayes_fit) {
          if ("correlation" %notin% names(obj)) {
            cli::cli_abort(
              "Random factor in 'design', but one or more MArrayLM objects in 'eBayes_fit' have no 'correlation' entry."
            )
          }
        }
      } else {
        if ("correlation" %notin% names(x$eBayes_fit)) {
          cli::cli_abort(
            "Random factor in 'design', but 'eBayes_fit' has no 'correlation' entry."
          )
        }
      }
    } else {
      # If no random factor is used, ensure the MArrayLM(s) do not have correlation
      if (is.list(x$eBayes_fit)) {
        for (obj in x$eBayes_fit) {
          if ("correlation" %in% names(obj)) {
            cli::cli_abort(
              "Model fit has a 'correlation' term, but 'design' has no random factor. Mismatch."
            )
          }
        }
      } else {
        if ("correlation" %in% names(x$eBayes_fit)) {
          cli::cli_abort(
            "Model fit has a 'correlation' term, but 'design' has no random factor. Mismatch."
          )
        }
      }
    }
  }
  
  
  # results checks
  if (!is.null(x$results)) {
    if (!is.null(x$design$contrast_matrix)) {
      expected_terms <- colnames(x$design$contrast_matrix)
      terms_from <- "contrast matrix"
    } else {
      expected_terms <- colnames(x$design$design_matrix)
      terms_from <- "design matrix"
    }
    
    # If user extracted intercept, results names match all terms
    # else skip the intercept
    if (isTRUE(x$tags$extract_intercept)) {
      if (!all(names(x$results) == expected_terms)) {
        cli::cli_abort(c(
          "Discrepancy between statistical results and design:",
          "{length(expected_terms)} Term name{?s} in {terms_from} do not match names in 'results'."
        ))
      }
    } else {
      non_intercept <- expected_terms[!grepl("ntercept", expected_terms)]
      if (!all(names(x$results) == non_intercept)) {
        cli::cli_abort(c(
          "Discrepancy between statistical results and design:",
          "Non-intercept terms in {terms_from} do not match 'results' names."
        ))
      }
    }
  }
    
    # Each results df must match rownames of data
  #   for (term in names(x$results)) {
  #     if (!all(rownames(x$results[[term]]) == rownames(x$data))) {
  #       cli::cli_abort(c(
  #         "Discrepancy between statistical results and data:",
  #         "Rownames of the {.arg {term}} results do not match rownames of 'data'."
  #       ))
  #     }
  #   }
  # }
    
    # Allow results to be a subset if filtered_proteins_per_contrast is present
    for (term in names(x$results)) {
      result_rows <- rownames(x$results[[term]])
      data_rows   <- rownames(x$data)
      
      if ("filtered_proteins_per_contrast" %in% names(x)) {
        expected_subset <- x$filtered_proteins_per_contrast[[term]]
        if (!all(result_rows %in% data_rows)) {
          cli::cli_abort(c(
            "Discrepancy between results and data for contrast {.val {term}}:",
            "Some rownames in results are not present in 'data'."
          ))
        }
        if (!all(result_rows %in% expected_subset)) {
          cli::cli_warn(c(
            "Some rows in the results for contrast {.val {term}}",
            "do not match the filtered protein list. Double-check filtering."
          ))
        }
      } else {
        if (!all(result_rows == data_rows)) {
          cli::cli_abort(c(
            "Discrepancy between statistical results and data:",
            "Rownames of the {.arg {term}} results do not match rownames of 'data'."
          ))
        }
      }
    }
    
  ## --- END CORE VALIDATIONS ---
  
  ## OPTIONAL: Validate 'filtered_proteins_per_contrast' if present
  if ("filtered_proteins_per_contrast" %in% names(x)) {
    if (!is.list(x$filtered_proteins_per_contrast)) {
      cli::cli_abort(
        "'filtered_proteins_per_contrast' must be a list (e.g., a list of retained IDs per contrast)"
      )
    }
    # Add any additional checks desired:
    # e.g. are they all character vectors of valid rownames?
  }
  
  # If all checks passed, explicitly return x
  #return(x)
  class(x) <- "DAList"
  return(x)
  
}





