#' Prepare limma contrasts matrix
#'
#' Create the contrasts matrix, for use in a limma model. This function utilizes the
#' function \link[limma:makeContrasts]{limma::makeContrasts} with a user provided file or vector of a list of comparisons.
#' Note: The label on the plots is defined by what is written in the contrast statement prior to the equal sign.
#'
#'
#' @param DAList A DAList. Must have a non-empty statistical design.
#' @param contrasts_vector A vector of contrasts.
#' @param contrasts_file The path to the contrasts file listing the desired contrasts.
#'   Must be a .csv, .tsv, or .txt file.
#'
#' @return A DAList with added contrasts associated with the limma design
#' @export
#'
#' @examples
#' \dontrun{
#' # An example of a .csv file with three comparisons
#' data -> add_contrasts(data, contrasts_file = "path/to/file.csv")
#'
#' # file info
#' Treatment1_vs_Control= Treatment1 - Control
#' Treatment2_vs_Control= Treatment2 - Control
#' Treatment2_vs_Treatment1= Treatment2 - Treatment1
#'
#' # An example using a vector
#' data <- add_contrasts(contrasts_vector = c("Treat1_vs_Control=Treat1-Control",
#'      "Treat2_vs_Control=Treat2-Control"))
#' }
#'
add_contrasts <- function(DAList,
                          contrasts_vector = NULL,
                          contrasts_file = NULL) {

  # Check input arguments generally
  validate_DAList(DAList)

  # Make sure there's a design matrix present already,
  # tell user to set it first if not
  if (is.null(DAList$design)) {
    cli::cli_abort(c("Input DAList does not have a statistical design",
                     "i" = "Run {.code DAList <- add_design(DAList, ~ formula)} before adding contrasts"))
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
  contrast_groups <- unique(unlist(extract_contrast_groups(contrast.vec = contrast_vec_squish)))

  # if the groups defined in the contrast file do no match the design
  # then error
  if (!all(contrast_groups %in% colnames(DAList$design$design_matrix))) {
    invalidGroups <- contrast_groups[contrast_groups %notin% colnames(DAList$design$design_matrix)]

    cli::cli_abort(c("Some groups present in the contrasts are not present in the design matrix",
                     "!" = "{cli::qty(length(invalidGroups))} Contrast group{?s} {.val {invalidGroups}} not found in the colnames of the design matrix",
                     "i" = "Check for typos in the contrasts",
                     "i" = "Make sure all groups for which you define contrasts are included in the design matrix.",
                     "i" = "Use {.code colnames(DAList$design$design_matrix)} to see the names of the design matrix.",
                     "i" = "If using an intercept model, you may need to reorder/relevel your factors."))
  }

  # Get contrasts with limma
  contrasts <- limma::makeContrasts(contrasts = contrast_vec_squish, levels = DAList$design$design_matrix)
  colnames(contrasts) <- stringr::str_remove(colnames(contrasts), "=.*")

  if (!is.null(DAList$design$contrast_matrix) | !is.null(DAList$design$contrast_vector)) {
    cli::cli_inform("DAList already contains contrasts. Overwriting.")
    # Get rid of any old stuff
    DAList$design$contrast_matrix <- NULL
    DAList$design$contrast_vector <- NULL
  }

  # Delete any existing model fits and results, which may no longer match the
  # contrasts
  if (!is.null(DAList$eBayes_fit)) {
    cli::cli_inform("DAList contains a model fit, deleting.")
    # Get rid of any old stuff
    DAList["eBayes_fit"] <- list(NULL)
  }
  if (!is.null(DAList$results)) {
    cli::cli_inform("DAList contains a statistical results, deleting.")
    # Get rid of any old stuff
    DAList["results"] <- list(NULL)
  }

  # Set everything here
  DAList$design$contrast_matrix <- contrasts
  DAList$design$contrast_vector <- contrasts_vector


  validate_DAList(DAList)
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
    tmp <- as.numeric(length(problem_rows))
    cli::cli_abort(c("Formatting issue in contrasts",
                     "!" = "{cli::qty(tmp)} Contrast{?s} {.val {problem_rows}} {cli::qty(tmp)} {?does/do} not contain an {.val =}"))
  }

  # Check that the contrast definition only has one equals sign
  if (any(stringr::str_count(contrast_vector, "=") > 1)) {
    problem_rows <- which(stringr::str_count(contrast_vector, "=") > 1)
    tmp <- as.numeric(length(problem_rows))
    cli::cli_abort(c("Formatting issue in contrast file",
                     "!" = "{cli::qty(tmp)} Contrast{?s} {.val {problem_rows}} {cli::qty(tmp)} {?contains/contain} more than one {.val =}"))
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
  processed_contrasts <-  contrasts |>
    stringr::str_remove_all(" ") |> # remove blank spaces
    stringr::str_remove_all("\\/[[:digit:]]+") |> ## removes division by number (/2, or /2000)
    stringr::str_remove_all("\\*[[:digit:]]+") |> ## removes multiplication by number (*2, or *2000)
    stringr::str_replace_all("\\+", "-") |>       ## replaces + with - (+ -> -)
    stringr::str_remove_all("[[()]]") |>          ## remove parentheses
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
