#' Prepare limma model design matrix
#'
#' Creates the model design matrix, targets data frame, and design formula for use in
#' limma. Either group_column or design_formula must be specified. If design_formula is
#' specified all other column parameters are ignored. Currently, if you want to use a
#' random factor model with paired_column, you cannot specify your model via design_formula.
#'
#' @param targets A targets dataframe. In the pipeline, usually the "targets" slot
#'   of the list output by \code{\link{process_data}}.
#' @param group_column The name of the column in the targets dataframe that
#'   describes the main groups to be compared in the limma model.
#' @param factor_columns Optional. The name of the column(s) in the targets dataframe
#'   that describes additional statistical factors for analysis. Default is NULL.
#' @param paired_column Optional. The name of the column in the targets dataframe
#'   that lists paired samples, for use in limma's paired mixed effect model.
#' @param design_formula Optional. A character string that can be coerced into a terms
#'   object that is appropriate for input as the object for stats:model.matrix(). When
#'   designformula is specified group_columns, factor_columns and paired_column is ignored.
#'
#' @return A list with three slots \enumerate{
#'   \item "design"- The model design matrix for the specified factors.
#'   \item "targets"- An updated targets data frame including only the grouping,
#'     factor, and paired columns specified.
#'   \item "designformula"- The model design formula for the specified factors.
#' }
#' @export
#'
#' @examples
#' # No examples yet
#'
make_design <- function(targets,
                        group_column = NULL,
                        factor_columns = NULL,
                        paired_column = NULL,
                        design_formula = NULL) {
  cli::cli_rule()

  ## if input values are all column names in targets
  if (is.null(design_formula)) {
    if (is.null(group_column)) {
      cli::cli_abort(c("Either {.arg design_formula} or {.arg group_column} must be specified"))
    }
    pass <- c(group_column, paired_column, factor_columns) %in% colnames(targets)

    if (!all(pass)) {
        # Stop if all input values are not column names in targets
        invalidCols <- c(group_column, paired_column, factor_columns)[pass == FALSE]
        cli::cli_abort(c("Some input column names not found in the targets data frame",
                        "x" = "{cli::qty(length(invalidCols))} Column{?s} {.val {invalidCols}} not present in {.arg targets}} "))
    }

    # Reorder columns for target
    tar <- targets[ , c(group_column, paired_column, factor_columns), drop = F]

    ## check that each factor contains 2 or more levels
    enough_levels <- apply(tar, 2, FUN = function(x) {length(unique(x)) >= 2})
    if (!all(enough_levels)) {
        invalidLevels <- names(enough_levels)[!enough_levels]
        cli::cli_abort(c("Columns must have at least two unique levels to create a design matrix",
                        "x" = "{cli::qty(length(invalidLevels))} Column{?s} {.val {invalidLevels}} {?has/have} 1 or fewer levels."))
    }

    ## make group a named list and create means model no intercept formula
    groupformula <- paste0("~0+", group_column)

    ## if paired factor supplied then design matrix and targets file for
    ## mixed effects model created. i.e. paired column in targets renamed
    ## "paired"
    if (!is.null(paired_column)) {
        ## make paired a named list. if paired.type is paired then
        ## add to formula.
        p.col <- grep(paired_column, colnames(tar))
        stopifnot(colnames(tar)[p.col] == paired_column) #TODO: DOES THIS EVER FAIL?
        # IT DOES IF YOU SUPPLY MULTIPLE COLUMNS HERE. BUT SHOULD WE BE CHECKING THAT ABOVE ANYWAY?
        colnames(tar)[p.col] <- "paired" ## rename column
        cli::cli_inform("Creating design matrix and targets for paired sample design using limma's mixed effect model")
        cli::cli_inform("Renamed input {.arg paired_column} {.val {paired_column}} to {.val paired}")
        paired_column <- NULL ## so paired is not included in design formula
    }
    # From my reading of the code, we always want pairedformula to be NULL
    # And it is never included in the design formula and matrix, only the
    # target??
    # TODO: double check on this, and streamline it out if that's true.
    pairedformula <- NULL


    ## named lists of factors. add to formula
    if (!is.null(factor_columns)) {
        facs <- as.list.data.frame(tar[ , factor_columns, drop = F])
    } else {
        facs <- NULL
    }
    additiveformula <- names(facs)

    # Coerce target columns into factors
    # adding the colname as a prefix to numeric vectors
    for (i in base::seq_along(tar)) {
        if (is.numeric(tar[,i]) | any(substr(tar[,i],1,1) %in% c(0:9))) {
        tar[,i] <- make_factor(tar[,i], prefix = colnames(tar)[i])
        } else {
        tar[,i] <- make_factor(tar[,i])
        }
    }

    ## DESIGN FORMULA
    # Start it with the group formula
    designformula <- groupformula
    # Add paired section if needed?
    if (length(pairedformula) > 0) { # TODO: DOES THIS EVER GET CALLED???
        # SEEMS THAT ABOVE WE SET pairedformula == NULL whether we're including
        # a pairing column or not??
        designformula <- paste(designformula, pairedformula, sep = "+")
    }
    # Add additive/factor sections if needed
    if (length(additiveformula) > 0) {
        additiveformula <- paste(additiveformula, collapse = "+")
        designformula <- paste(designformula, additiveformula, sep = "+")
    }
    ## create design matrix
    # TODO: explicitly set up contrasts here?
    # right now, as far as I can tell, it relies on the
    # environment options, which could lead to some very weird behavior if
    # they were ever changed accidentally.
    design <- stats::model.matrix(eval(parse(text = designformula)), data = tar)
    # TODO: Did some testing by putting, e.g., replicate as a factor
    # The design matrix I got out used different encoding: left off replicate 1
    # The group part of the matrix used a full encoding, while the non-group
    # part used use encoding where one level was all 0s.

    ## change column names of design.
    desCols <- levels(tar[,group_column])
    for (x in c(paired_column, factor_columns)) {
        desCols <- c(desCols, levels(tar[,x])[-1])
        # TODO: why drop the first level of the paired or factor cols here??
        # Rather, I get why we do it: to make what the design matrix makes above
        # But why does the design matrix do that above? I guess it defaults to treatment
        # contrasts for the later terms of the model formula?
    }
    colnames(design) <- desCols

    output <- list(design = design,
                   targets = tar,
                   designformula = designformula)

  } else { # When supplying design_formula

      formulaobject <- stats::terms(eval(parse(text=design_formula)))
      if (class(formulaobject)[2] != "formula") {
          cli::cli_abort(c("String provided in {.arg design_formula} could not be coerced to a formula"))
      }
      if (!is.null(group_column)) {
          cli::cli_inform(c("{.arg group_column} is ignored when {.arg design_formula} is provided"))
      }
      ## if paired factor supplied then design matrix and targets file for
      ## mixed effects model created. i.e. paired column in targets renamed
      ## "paired"
      if (!is.null(paired_column)) {
          ## make paired a named list. if paired.type is paired then
          ## add to formula.
          p.col <- grep(paired_column, colnames(targets))
          stopifnot(colnames(targets)[p.col] == paired_column) #TODO: DOES THIS EVER FAIL?
          # IT DOES IF YOU SUPPLY MULTIPLE COLUMNS HERE. BUT SHOULD WE BE CHECKING THAT ABOVE ANYWAY?
          colnames(targets)[p.col] <- "paired" ## rename column
          cli::cli_inform("Creating design matrix and targets for paired sample design using limma's mixed effect model")
          cli::cli_inform("Renamed input {.arg paired_column} {.val {paired_column}} to {.val paired}")
          paired_column <- NULL ## so paired is not included in design formula
      }
      # From my reading of the code, we always want pairedformula to be NULL
      # And it is never included in the design formula and matrix, only the
      # target??
      # TODO: double check on this, and streamline it out if that's true.
      pairedformula <- NULL
      if (!is.null(factor_columns)) {
          cli::cli_inform(c("{.arg factor_columns} is ignored when {.arg design_formula} is provided"))
      }
      if (attr(formulaobject, which="intercept") == 0) {
          cli::cli_inform(c("A no-intercept model has been specified"))
      } else{
          cli::cli_inform(c("An intercept model has been specified"))
      }

      tar <- targets[ , rownames(attr(formulaobject, which = "factors")), drop = F]
      formulaobject <- stats::terms(eval(parse(text = design_formula)), data = tar)
      design <- stats::model.matrix(eval(parse(text = design_formula)), data = tar)
      extratext <- rownames(attr(formulaobject, which="factors"))

      # Fix colnames of design matrix to be compatible with limma
      colnames(design) <- stringr::str_remove_all(colnames(design), paste(extratext, collapse="|"))
      colnames(design) <- stringr::str_replace_all(colnames(design), "\\(Intercept\\)", "Intercept")
      colnames(design) <- stringr::str_replace_all(colnames(design), "\\:", ".")

      output <- list(design = design,
                     targets = tar,
                     designformula = design_formula)
    }

  cli::cli_inform("Design matrix and targets created")
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  output

}
