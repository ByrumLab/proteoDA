#' Prepare limma contrasts matrix
#'
#' Create the contrasts matrix, for use in limma. Calls
#' \code{\link{extract_contrast_groups}} and \code{\link{make_log}} as
#' subfunctions.
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
  filext <- tools::file_ext(file)
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

  ## imports contrast file
  ## TODO: From Duah's experience that one time, when she had a minus sign instead
  # of an equals sign: might be good
  # to do more checking of the format of this file
  # frankly, there shouldn't be any more than one columns:
  # should probably through an error if there is.
  contrast.vec <- utils::read.delim(file = file,
                                    sep = sep,
                                    stringsAsFactors = F,
                                    header = F)
  contrast.vec <- contrast.vec[,1]
  contrast.vec <-gsub(" ", "", contrast.vec)

  ## extract groups included in each contrast
  contrast.grps <- extract_contrast_groups(contrast.vec = contrast.vec)

  ## if the groups defined in the contrast file do no match the design
  ## then return the imported contrast.vec.
  if (!all(unique(unlist(contrast.grps)) %in% colnames(design))) {
    invalidGroups <- unique(unlist(contrast.grps))[unique(unlist(contrast.grps)) %notin% colnames(design)]

    cli::cli_inform(cli::col_yellow("Some groups present in contrast file do not match design matrix"))
    cli::cli_inform("{cli::qty(length(invalidGroups))} Contrast group{?s} {.val {invalidGroups}} not found in design matrix")
    cli::cli_inform("Returning the imported contrast.vec for additional processing")
    cli::cli_inform("Define contrasts manually using the {.fn limma::makeContrasts} function, e.g.:")
    cli::cli_inform(c("{.code contrasts <- makeContrasts(contrasts=contrast.vec, levels=design)}",
                      "{.code colnames(contrasts) <- gsub('=.*','',colnames(contrasts))}"))
    return(contrast.vec)
  }
  ## define contrasts
  contrasts <- limma::makeContrasts(contrasts = contrast.vec, levels = design)
  # rename cols to just the name of the contrasts, not the full vector
  colnames(contrasts) <- gsub("=.*", "", colnames(contrasts))

  # Write logs
  param[["file"]] <- file
  stats[["no_contrasts"]] <- ncol(contrasts)
  logs <- make_log(param = param,
                   stats = stats,
                   title = "CONTRASTS",
                   save = TRUE)

  # return data
  list(contrasts = contrasts,
       contrast.vec = contrast.vec,
       param = logs$param,
       stats = logs$stats)
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
  # TODO: Issue with all the single-digit operation removals is that they
  # won't remove more digits, if there's an operation by a two-digit number
  # e.g., /22 will get turned into just "2".
  # Seems suboptimal, but I'm not really sure what the point of this section is
  # I've changed the regexes so that they will remove operations by any # of digits.
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
  cli::cli_inform("Contrast file imported successfully.")
  output_list
}
