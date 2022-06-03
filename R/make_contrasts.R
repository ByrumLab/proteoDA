#' Prepare limma contrasts matrix
#'
#' Create the contrasts matrix, for use in limma. Calls
#' \code{\link{extract_contrast_groups}} \code{\link{validate_contrasts}}
#' and \code{\link{make_log}} as subfunctions.
#'
#' @param file The path to the contrasts file listing the desired contrasts.
#'   Must be a .csv, .tsv, or .txt file.
#' @param design The limma design matrix for the model. In the pipeline, this
#'   is usually the "design" slot of the list returned by
#'   \code{\link{make_design}}.
#'
#' @return A list with 4 slots: \enumerate{
#'   \item "contrasts"- The contrast matrix, for use in limma
#'   \item "contrast.vec"- The vector of contrasts, as read from the input
#'     contrasts file
#'   \item "param"- A dataframe of parametes used for calling the function
#'   \item "stats"- A dataframe giving stats (the total number of contrasts).
#' }
#' @export
#'
#' @examples
#' # No examples yet
#'
make_contrasts <- function(file = NULL,
                           design) {

  param <- stats <- list()

  # TODO: I think we should take away this functionality
  # Script isn't trackable/reproducible if we're selecting files
  # through clicking in Rstudio
  if (is.null(file)) {
    file = file.choose()
  }

  ## check that the file is a csv, tsv, or text file
  filext <- file_extension(file)
  if (filext %notin% c("csv","txt","tsv")) {
    cli::cli_abort(c("Problem with input file",
                     "x" = "Input file must end in {.file .csv}, {.file .tsv}, or {.file .txt}"))
  }

  ## check that file exists
  if (!file.exists(file)) {
    cli::cli_abort(c("Cannot find file.",
                     "x" = "{.file {file}} does not exist",
                     "i" = "Did you specify the filepath correctly?"))
  }

  if (filext == "txt" | filext=="tsv") {sep= "\t"}
  if (filext == "csv") {sep=","}

  cli::cli_rule()
  ## imports contrast file
  contrast.tbl <- utils::read.delim(file = file,
                                    sep = sep,
                                    stringsAsFactors = F,
                                    header = F)
  cli::cli_inform("Contrast file imported, checking format...")
  validate_contrasts(contrast.tbl)
  cli::cli_inform("Contrast file passed formatting checks")

  contrast.vec <- contrast.tbl[,1, drop = T]
  contrast.vec <-gsub(" ", "", contrast.vec)

  ## extract groups included in each contrast
  contrast.grps <- extract_contrast_groups(contrast.vec = contrast.vec)

  ## if the groups defined in the contrast file do no match the design
  ## then return the imported contrast.vec.
  if (!all(unique(unlist(contrast.grps)) %in% colnames(design))) {
    invalidGroups <- unique(unlist(contrast.grps))[unique(unlist(contrast.grps)) %notin% colnames(design)]

    cli::cli_abort(c("Some groups present in contrast file are not present in the design matrix",
                     "!" = "{cli::qty(length(invalidGroups))} Contrast group{?s} {.val {invalidGroups}} not found in the colnames of the design matrix",
                     "i" = "Check for typos of group names in the contrast file,",
                     "i" = "And make sure all groups for which you define contrasts are included in the statistical design."))
  }
  ## define contrasts
  contrasts <- limma::makeContrasts(contrasts = contrast.vec, levels = design)
  # rename cols to just the name of the contrasts, not the full vector
  colnames(contrasts) <- gsub("=.*", "", colnames(contrasts))


  #TODO: add in a check on rank of matrix? (see the F1000 paper Charity sent)

  # Write logs
  param[["file"]] <- file
  stats[["num_contrasts"]] <- ncol(contrasts)
  logs <- make_log(param = param,
                   stats = stats,
                   title = "CONTRASTS",
                   save = TRUE)

  cli::cli_inform("Contrasts imported successfully")
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  # return data
  list(contrasts = contrasts,
       contrast.vec = contrast.vec,
       param = logs$param,
       stats = logs$stats)
}


#' Validate contrasts format
#'
#' Internal function that checks the format of the imported contrasts
#' and provides informative error messages about issues.
#'
#' @param contrast.tbl The imported data frame of contrasts.
#'
#' @return If all checks pass, invisible(TRUE).
#'
#' @examples
#' # No examples yet
validate_contrasts <- function(contrast.tbl) {

  # Check that there's only one column
  if (ncol(contrast.tbl) != 1) {
    cli::cli_abort(c("Imported contrasts file contains more than one column",
                     "i" = "Check the formatting of your contrasts file"))
  }

  # Check that the contrast definition has an equals sign in it
  if (any(stringr::str_detect(contrast.tbl[,1], "=", negate = T))) {
    problem_rows <- which(stringr::str_detect(contrast.tbl[,1], "=", negate = T))
    cli::cli_abort(c("Formatting issue in contrast file",
                    "!" = "{cli::qty(length(problem_rows))} Line{?s} {.val {problem_rows}} {?do/does} not contain an {.val =}"))
  }

  # Check that the contrast definition only has one equals sign
  if (any(stringr::str_count(contrast.tbl[,1], "=") > 1)) {
    problem_rows <- which(stringr::str_count(contrast.tbl[,1], "=") > 1)
    cli::cli_abort(c("Formatting issue in contrast file",
                     "!" = "{cli::qty(length(problem_rows))} Line{?s} {.val {problem_rows}} {?contain/contains} more than one {.val =}"))
  }


  # Check that the names portion of the contrast definition includes a "_vs_"
  contrast_names <- vapply(contrast.tbl[,1],
                           function(x) {
                             stringr::str_split(x, "=")[[1]][1]
                             },
                           FUN.VALUE = "text",
                           USE.NAMES = F)
  if (any(stringr::str_detect(contrast_names, "_vs_", negate = T))) {
    problem_rows <- which(stringr::str_detect(contrast_names, "_vs_", negate = T))
    cli::cli_abort(c("Formatting issue in contrast file",
                     "!" = "The name portion of the contrast (before the \"=\")",
                     "!" = "in {cli::qty(length(problem_rows))} line{?s} {.val {problem_rows}} does not contain {.val _vs_}"))
  }

  # If we pass all checks, return true invisibly
  invisible(TRUE)
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
#' @export
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
    stringr::str_remove_all(" ") %>% # remove blank spaces TODO: probably unneeded
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
