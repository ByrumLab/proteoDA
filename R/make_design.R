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

  # After all checks, strip out any random effects into a fixed-only formula
  # STILL NEED TO IMPLEMENT
  # Som epossible code ideas: update(mod, formula=drop.terms(mod$terms, 2, keep.response=TRUE)  )
  form_fixed_only <- formula

  # Make the model matrix
  design_matrix <- stats::model.matrix(form_fixed_only, data = DIAlist$metadata)
  colnames(design_matrix) <- stringr::str_replace_all(colnames(design_matrix), "\\(Intercept\\)", "Intercept")
  colnames(design_matrix) <- stringr::str_replace_all(colnames(design_matrix), "\\:", ".")

  # then, process fixed effects
  # Or maybe lower down..

  # Do I need to do anything to the metadata?
  # I don't think so???

  if (!is.null(DIAlist$design)) {
    cli::cli_inform("DIAlist already contains a statistical design. Overwriting.")
  }

  DIAlist$design <- list(design_formula = paste0(as.character(formula), collapse = ""),
                         design_matrix = design_matrix)

  # Add fixed effect info as needed.
  # unlike before, won't mess with the metadata at all
  # will just add in the design list a slot for the
  # random factor

  validate_DIAlist(DIAlist)

  #
  #     tar <- targets[ , rownames(attr(formulaobject, which = "factors")), drop = F]
  #     formulaobject <- stats::terms(eval(parse(text = design_formula)), data = tar)
  #     design <- stats::model.matrix(eval(parse(text = design_formula)), data = tar)
  #     extratext <- rownames(attr(formulaobject, which="factors"))
  #
  #     # Fix colnames of design matrix to be compatible with limma
  #     colnames(design) <- stringr::str_remove_all(colnames(design), paste(extratext, collapse="|"))
  #     colnames(design) <- stringr::str_replace_all(colnames(design), "\\(Intercept\\)", "Intercept")
  #     colnames(design) <- stringr::str_replace_all(colnames(design), "\\:", ".")
  #
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
    random_term <- terms[stringr::str_detect(terms, "\\|")]
    print("Random term:")
    print(random_term)
  }

  ## NEED TO IMPLEMENT A BUNCH OF CHECKS FOR APPROPRIATENESS
  # can't have more than one term within the intercept part of the random effect
  # can't have nested terms in the random effect: right side of random effect should be 1 word, no special characters
  # can't have random slopes, only random intercepts.: left side of random effect with | must be 1


  formula
}
