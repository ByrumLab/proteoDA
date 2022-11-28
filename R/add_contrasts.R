#' Prepare limma contrasts matrix
#'
#' Create the contrasts matrix, for use in limma.
#'
#' @param DIAlist A DIAlist. Must have a non-empty statistical design.
#' @param contrasts_vector A vector of contrasts
#' @param contrasts_file The path to the contrasts file listing the desired contrasts.
#'   Must be a .csv, .tsv, or .txt file.
#'
#' @return A DIAlist with added contrasts
#' @export
#'
#' @examples
#' # No examples yet
#'
add_contrasts <- function(DIAlist,
                          contrasts_vector = NULL,
                          contrasts_file = NULL) {

  # Check input arguments generally
  validate_DIAlist(DIAlist)

  # Make sure there's a design matrix present already,
  # tell user to set it first if not
  if (is.null(DIAlist$design)) {
    cli::cli_abort(c("Input DIAlist does not have a statistical design",
                     "i" = "Run {.code DIAlist <- add_design(DIAlist, ~ formula)} before adding contrasts"))
  }

  # Must have either contrasts vec or contrasts file, not both.
  if (sum(!is.null(contrasts_vector), !is.null(contrasts_file)) != 1) {
    cli::cli_abort(c("Must supply either {.arg contrasts_vector} OR {.arg contrasts_file}"))
  }

  # When contrasts_file is supplied, process it into a contrasts_vector
  if (!is.null(contrasts_file)) {
    # check that the file is a csv, tsv, or text file
    filext <- stringr::str_to_lower(file_extension(contrasts_file))
    if (filext %notin% c("csv","txt","tsv")) {
      cli::cli_abort(c("Problem with {.arg contrasts_file}",
                       "x" = "Input file must end in {.file .csv}, {.file .tsv}, or {.file .txt}"))
    }

    # check that file exists
    if (!file.exists(contrasts_file)) {
      cli::cli_abort(c("Cannot find contrasts file.",
                       "x" = "{.file {contrasts_file}} does not exist",
                       "i" = "Did you specify the filepath correctly?"))
    }

    if (filext == "txt" | filext=="tsv") {sep= "\t"}
    if (filext == "csv") {sep=","}


    # imports contrast file
    contrast_table <- utils::read.delim(file = contrasts_file,
                                      sep = sep,
                                      stringsAsFactors = F,
                                      header = F)

    # Check that there's only one column
    if (ncol(contrast_table) != 1) {
      cli::cli_abort(c("Imported contrasts file contains more than one column",
                       "i" = "Check the formatting of your contrasts file"))
    }
    contrasts_vector <- contrast_table[,1, drop = T]
  }

  contrasts_vector <- validate_contrasts(contrasts_vector)
  contrast_vec_squish <- stringr::str_remove_all(contrasts_vector, " ")

  # extract groups included in each contrast
  contrast_groups <- extract_contrast_groups(contrast.vec = contrast_vec_squish)

  # if the groups defined in the contrast file do no match the design
  # then error
  if (!all(unique(unlist(contrast_groups)) %in% colnames(DIAlist$design$design_matrix))) {
    invalidGroups <- unique(unlist(contrast_groups))[unique(unlist(contrast_groups)) %notin% colnames(DIAlist$design$design_matrix)]

    cli::cli_abort(c("Some groups present in the contrasts are not present in the design matrix",
                     "!" = "{cli::qty(length(invalidGroups))} Contrast group{?s} {.val {invalidGroups}} not found in the colnames of the design matrix",
                     "i" = "Check for typos in the contrasts",
                     "i" = "Make sure all groups for which you define contrasts are included in the design matrix.",
                     "i" = "Use {.code colnames(DIAlist$design$design_matrix)} to see the names of the design matrix.",
                     "i" = "If using an intercpet model, you may need to reorder/relevel your factors."))
  }

  # Get contrasts with limma
  contrasts <- limma::makeContrasts(contrasts = contrast_vec_squish, levels = DIAlist$design$design_matrix)
  colnames(contrasts) <- stringr::str_remove(colnames(contrasts), "=.*")

  if (!is.null(DIAlist$design$contrast_matrix) | !is.null(DIAlist$design$contrast_vector)) {
    cli::cli_inform("DIAlist already contains contrasts. Overwriting.")
    # Get rid of any old stuff
    DIAlist$design$contrast_matrix <- NULL
    DIAlist$design$contrast_vector <- NULL
  }

  # Set everything here
  DIAlist$design$contrast_matrix <- contrasts
  DIAlist$design$contrast_vector <- contrasts_vector


  validate_DIAlist(DIAlist)
}


#' Validate contrasts format
#'
#' Internal function that checks the format of a contrasts vector
#' and provides informative error messages about issues.
#'
#' @param contrast_vector A
#'
#' @return If all checks pass, invisibly returns contrast vector.
#'
#' @examples
#' # No examples yet

validate_contrasts <- function(contrast_vector) {
  # Check that the contrast definition has an equals sign in it
  if (any(stringr::str_detect(contrast_vector, "=", negate = T))) {
    problem_rows <- which(stringr::str_detect(contrast_vector, "=", negate = T))
    cli::cli_abort(c("Formatting issue in contrasts",
                     "!" = "{cli::qty(length(problem_rows))} Contrasts{?s} {.val {problem_rows}} {?do/does} not contain an {.val =}"))
  }

  # Check that the contrast definition only has one equals sign
  if (any(stringr::str_count(contrast_vector, "=") > 1)) {
    problem_rows <- which(stringr::str_count(contrast_vector, "=") > 1)
    cli::cli_abort(c("Formatting issue in contrast file",
                     "!" = "{cli::qty(length(problem_rows))} Contrasts{?s} {.val {problem_rows}} {?contain/contains} more than one {.val =}"))
  }

  # If we pass all checks, return contrasts invisibly
  invisible(contrast_vector)
}


#' Extract groups from a vector of contrasts
#'
#' @param contrast.vec A vector of contrasts, most likely from reading in
#'   a contrasts text file
#'
#' @return A named list of length equal to the number of contrasts. Each element
#'   in the list is a character vector of the groups involved in that contrast,
#'   and the names of the list give the contrast names.
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
extract_contrast_groups <- function(contrast.vec) {

  # Remove blank spaces
  input <- stringr::str_squish(as.character(contrast.vec))

  # Split at equal sign, into name of comparison
  # and definittion of comparison.
  input_split <- matrix(unlist(stringr::str_split(input, "=")), ncol = 2, byrow = T)
  contrast_names <- input_split[,1]
  contrasts <- input_split[,2]

  # Basically replicating a lot of the stuff that
  # is already there, but with stringr. Though I have no idea what some of this is for
  processed_contrasts <-  contrasts %>%
    stringr::str_remove_all(" ") %>% # remove blank spaces
    stringr::str_remove_all("\\/[[:digit:]]+") %>% ## removes division by one digit number (/2)
    stringr::str_remove_all("\\-[[:digit:]]+") %>% ## removes subtraction by one digit number (-2)
    stringr::str_remove_all("\\+[[:digit:]]+") %>% ## removes addition of one digit number (+2)
    stringr::str_remove_all("\\*[[:digit:]]+") %>% ## removes multiplication by one digit number (*2)
    stringr::str_replace_all("\\+", "-") %>%       ## replaces + with - (+ -> -)
    stringr::str_remove_all("[[()]]") %>%          ## remove parentheses
    stringr::str_replace_all("\\/", "-")           ## replaces / with -

  # Split on "-", to get just the groups.
  tmp_split <- strsplit(processed_contrasts, "-")
  output_list <- NULL
  for (i in 1:length(tmp_split)) {
    output_list[[i]] <- unique(tmp_split[[i]])
  }
  names(output_list) <- contrast_names
  output_list
}
