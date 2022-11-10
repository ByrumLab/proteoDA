#' Prepare limma model design matrix
#'
#' Need to rewrite
#'
#' @param DIAlist The DIAlist to be filtered
#' @param design_formula NEED TO UPDATE WITH SOME SPECIFICS ON WHATS ALLOWED
#'
#' @return A DIAlist object, with a design
#' @export
#'
#' @examples
#' # No examples yet
#'
#'
#'
add_design <- function(DIAlist,
                       design_formula = NULL) {

  # Check input arguments generally
  validate_DIAlist(DIAlist)

  # validate and parse formula
  formula <- validate_formula(design_formula)

  # Check that the terms in the user-supplied formula are present in
  # the metadata
  if (!all(all.vars(formula) %in% colnames(DIAlist$metadata))) {
    problem_terms <- all.vars(formula)[all.vars(formula) %notin% colnames(DIAlist$metadata)]

    cli::cli_abort(c("{cli::qty(problem_terms)} Term{?s} in {.arg design_formula} not found in metadata.",
                     "x" = "Missing term{?s}: {problem_terms}"))
  }

  formula_terms <- attributes(terms(formula))$term.labels
  formula_has_random <- any(stringr::str_detect(formula_terms, "\\|"))

  # After all checks, strip out any random effects into a fixed-only formula
  if (formula_has_random) {
    form_fixed_only <- drop.terms(terms(formula),
                                  dropx = which(stringr::str_detect(formula_terms, "\\|")))
    elements <- formula_terms[stringr::str_detect(formula_terms, "\\|")] %>%
      stringr::str_split_fixed(pattern = "\\|", n = 2)

    random_factor <- stringr::str_remove_all(elements[[2]], " ")
  } else {
    form_fixed_only <- formula
  }

  # Make the model matrix
  design_matrix <- stats::model.matrix(form_fixed_only, data = DIAlist$metadata)
  extratext <- rownames(attr(terms(form_fixed_only), which="factors"))
  colnames(design_matrix) <- stringr::str_replace_all(colnames(design_matrix), "\\(Intercept\\)", "Intercept")
  colnames(design_matrix) <- stringr::str_replace_all(colnames(design_matrix), "\\:", ".")
  colnames(design_matrix) <- stringr::str_remove_all(colnames(design_matrix), paste(extratext, collapse="|"))

  if (!is.null(DIAlist$design)) {
    cli::cli_inform("DIAlist already contains a statistical design. Overwriting.")
    # Get rid of any old stuff
    DIAlist$design <- NULL
  }

  DIAlist$design <- list(design_formula = paste0(as.character(formula), collapse = ""),
                         design_matrix = design_matrix)

  # If theres a random factor, add it in
  if (formula_has_random) {
    DIAlist$design$random_factor <- random_factor
  }

  validate_DIAlist(DIAlist)
}


#' Validate a statistical formula
#'
#' An internal function to validate
#'
#' @param design_formula A string or expression which represents a statistical formula
#'
#' @return If all checks pass, a formula object of the design_formula
#'
#' @examples
#' # No examples yet
#'
#'

validate_formula <- function(design_formula) {

  # First, try to coerce to formula
  formula <- try(formula(eval(parse(text = design_formula))), silent = T)
  if (inherits(formula, "try-error")) {
    cli::cli_abort(c("Could not parse supplied design_formula",
                     "!" = "Did you forget the {.val ~} or include a response variable?"))
  }

  # Then, check that the formula is only the RHS
  if (attributes(terms(formula))$response == 1) {
    cli::cli_abort(c("{.arg design_formula} must be the right-hand side of a statistical formula.",
                     "i" = "It must be a string or expression that starts with {.val ~}"))
  }

  # Extract terms
  terms <- attributes(terms(formula))$term.labels

  # Multiple terms within parentheses not allowed
  if (stringr::str_count(paste0(as.character(formula), collapse = ""), "\\(") > 1) {
    cli::cli_abort(c("{.arg design_formula} can only contain 1 random factor term",
                     "i" = "Do not include multiple sets of terms within parentheses"))
  }
  # Multiple terms with | not allowed, no || allowed
  if (stringr::str_count(paste0(as.character(formula), collapse = ""), "\\|") > 1) {
    cli::cli_abort(c("{.arg design_formula} can only contain 1 random factor term",
                     "i" = "Do not include multiple terms with the | character",
                     "i" = "Do not include any terms with ||"))
  }

  # All | terms, if they exist, must be enclosed in parentheses.
  if (stringr::str_detect(paste0(as.character(formula), collapse = ""), "\\|")) {
    if (!stringr::str_detect(paste0(as.character(formula), collapse = ""), "\\(.*\\|.*\\)")) {
      cli::cli_abort(c("In {.arg design_formula}, random effects terms must be enclosed in parentheses.",
                       "i" = "The | character can only be used within parentheses."))
    }
  }

  if (any(stringr::str_detect(terms, "\\|"))) {
    # Grab the random term
    random_term <- terms[stringr::str_detect(terms, "\\|")]

    # Get the randomg group term and the grouping
    elements <- stringr::str_split_fixed(random_term, pattern = "\\|", n = 2)
    rand_group_term <- stringr::str_remove_all(elements[[1]], " ")
    rand_groups <- stringr::str_remove_all(elements[[2]], " ")

    # Random group term has to equal "1"
    if (!identical(rand_group_term, "1")) {
      cli::cli_abort(c("Random effects can only influence the intercept, and not other terms",
                       "i" = "The random effect must be specified as (1 | random_factor)"))
    }

    # Random groups has to be only one factor
    if (stringr::str_count(rand_groups, "\\w+") > 1) {
      cli::cli_abort(c("Only one grouping/blocking factor can be specified as a random effect.",
                       "i" = "The random effect must be specified as (1 | random_factor)",
                       "i" = "Multiple effects like (1 | batch + group) aren't allowed"))
    }

    # Random group cannot be the same as one of the other terms (right?)
    if (rand_groups %in% terms[!stringr::str_detect(terms, "\\|")]) {
      cli::cli_abort(c("Random factors cannot also be fixed factors."))
    }
  }

  # If all checks pass, return the input formula as a formula
  formula
}
