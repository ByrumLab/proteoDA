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


# design_formula <- "~ 0 + group"
# design_formula <- "group"
# design_formula <- "~ group + (1 | replicate)"
# design_formula <- "~ group + (replicate | replicate)"
# design_formula <- "~ group*sample + XXX + (replicate | replicate)"

add_design <- function(DIAlist,
                       design_formula = NULL) {

  # Check input arguments generally
  validate_DIAlist(DIAlist)

  ## NEED TO FIX THIS CHECK
  # DOESN"T WORK PROPERLY FOR STRINGS WITHOUT ~
  if ("formula" %notin% class(eval(parse(text=design_formula)))) {
    cli::cli_abort(c("Could not process {.arg design_formula} as a formula"))
  }

  formula <- formula(eval(parse(text=design_formula)))

  # Check that the terms in the user-supplied formula are present in
  # the metadata
  if (!all(all.vars(formula) %in% colnames(DIAlist$metadata))) {
    problem_terms <- all.vars(formula)[all.vars(formula) %notin% colnames(DIAlist$metadata)]

    cli::cli_abort(c("{cli::qty(problem_terms)} Term{?s} in design formula not found in metadata.",
                     "x" = "Missing term{?s}: {problem_terms}"))
  }

  ## NEED TO IMPLEMENT A BUNCH OF CHECKS FOR APPROPRIATENESS
  # checks: can't have more than 1 random effect
  # can't have more than one term within the intercept part of the random effect
  # can't have nested terms in the random effect
  # can't have random slopes, only random intercepts.

  # After all checks, strip out any random effects into a fixed-only formula
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
