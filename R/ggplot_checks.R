
#' Validate Metadata Structure
#'
#' This internal function checks the validity of a metadata \code{data.frame}, ensuring
#' it meets structural requirements such as unique column names, non-empty row names,
#' and appropriate formatting.
#'
#' @param x A \code{data.frame} containing sample metadata. Each row represents a sample.
#' @param verbose A logical value indicating whether to check if row and column names exceed
#'   a character length threshold. Default is \code{FALSE}.
#'
#' @return Invisibly returns the validated input \code{data.frame}.
#'
#' @details
#' - Ensures \code{x} is a \code{data.frame} and not empty.
#' - Checks that column names are present, unique, and free from \code{NA} or blank values.
#' - Ensures row names are not missing, blank, or duplicated.
#' - If \code{verbose = TRUE}, warns if row or column names exceed 15 characters.
#'
#' @keywords internal
#' @rdname check_metadata
#'
#' @importFrom cli cli_abort cli_warn
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#' metadata <- data.frame(sample = c("C1", "C2", "T1", "T2"),
#'                        sample_id = c("C100", "C202", "T303", "T100"),
#'                        condition = c("Control", "Control", "Treatment",
#'                                      "Treatment"))
#'
#' metadata2 <- data.frame(sample = c("C1", "C2", "T1", "T2"),
#'                         sample = c("C100", "C202", "T303", "T100"),
#'                         condition = c("Control", "Control", "Treatment",
#'                                       "Treatment"),
#'                         check.names = FALSE)
#'
#' metadata3 <- data.frame(long_sample_column_name = c("C1", "C2", "T1", "T2"),
#'                         sample_id = c("C100", "C202", "T303", "T100"),
#'                         condition = c("Control", "Control", "Treatment",
#'                                       "Treatment"),
#'                         row.names = c("C100", "C202", "T303",
#'                                       "T100_long_row_name"))
#'
#' metadata4 <- data.frame(sample_id = c("C1", "C2", "T1", "T2"),
#'                         sample = c("C100", "C202", "T303", "T100"),
#'                         condition = c("Control", "Control", "Treatment",
#'                                       "Treatment"),
#'                         row.names = c("C1", "C2", "T2", ""), check.rows = FALSE)
#'
#'
#' metadata  <- check_metadata(x = metadata, verbose = TRUE)
#' metadata2 <- check_metadata(x = metadata2, verbose = TRUE)
#' metadata3 <- check_metadata(x = metadata3, verbose = TRUE)
#' metadata4 <- check_metadata(x = metadata4, verbose = TRUE)
#'
#'  }
#' }
#'
check_metadata <- function(x,
                           verbose = FALSE) {


 ## check metadata is data.frame
 if(!is.data.frame(x)) {
  cli::cli_abort("The {.arg metadata} argument is not a data.frame.")
 }


 ## check metadata is not an empty data.frame
 is_empty <- dim(x) == 0L
 if(any(is_empty)){
  cli::cli_abort(c("The metadata argument appears to be an empty data.frame."))
 }


 ## check column names are not absent (ie. NULL)
 null_colms <- is.null(colnames(x))
 if (any(null_colms)) {
  cli::cli_abort("The metadata data.frame does not have any column names
                   (i.e. NULL).")
 }


 ## check column names do not include NAs, blanks, or duplicates
 na_colms <- is.na(colnames(x))
 if(any(na_colms)) {
  cli::cli_abort("Column names of the metadata data.frame
                   should not include missing values.")
 }

 ## check if any values are empty strings i.e. blank
 blank_colms <- nchar(colnames(x)) == 0L
 if(any(blank_colms)) {
  cli::cli_abort("Column names of the metadata data.frame
                   should not include blank values.")
 }

 ## check column names are unique
 dup_column_names <- check_dup(x = colnames(x))
 if(!is.null(dup_column_names)) {
  cli::cli_abort(c("The metadata data.frame should contain unique
                     column names.",
                   "x" = "The following column names are duplicated:",
                   "{.val {dup_column_names}}"))
 }


 ## check row names are not absent (ie. NULL)
 null_rows <- is.null(rownames(x))
 if (any(null_rows)) {
  cli::cli_abort("The metadata data.frame does not have any row names
                   (i.e. NULL).")
 }


 ## check row names do not include NAs, blanks, or duplicates
 ## NEED TO INCLUDE WHEN NA, NA.1, NA.2 ... AND ALL VALUES ARE NA
 ## rowSums(is.na(metadata)) == ncol(metadata)
 na_rows <- is.na(rownames(x))
 if(any(na_rows)){
  cli::cli_abort("Row names of the metadata data.frame should not include
                   missing values")
 }


 ## check if any values are blank
 blank_rows <- nchar(row.names(x)) == 0L
 if(any(blank_rows)) {
  cli::cli_abort("Row names of the metadata data.frame should not
                   include blank values.")
 }

 ## check verbose argument
 if(!is.logical(x = verbose)){
  cli::cli_abort("The verbose argument should be either TRUE or FALSE")
 }


 ## if verbose = TRUE warn if any column or row names are
 ## longer than 15 characters in length.
 if(verbose){

  ## warn if column names are > 15 characters in length
  col_char_thresh <- 15
  long_column_names <- check_long(x = colnames(x),
                                  char_thresh = col_char_thresh)
  if(!is.null(long_column_names)) {
   cli::cli_warn(c("Some column names of the metadata data.frame are
                    very long (> {.val {col_char_thresh}} characters).",
                   "i" = "You may consider shortening the following column
                    names to avoid", "sub-optimal visualizations downstream:",
                   "{.val {long_column_names}}"))
  }


  ## warn if row names are > 15 characters in length
  row_char_thresh <- 15
  long_row_names <- check_long(x = row.names(x),
                               char_thresh = row_char_thresh)
  if(!is.null(long_row_names)) {
   if(length(long_row_names) > 20L) {
    long_row_names <- long_row_names[1:20]
   }
   cli::cli_warn(c("Some row names of the metadata data.frame
                    are very long (> {.val {row_char_thresh}} characters).",
                   "i" = "You may consider shortening the following row
                    names to avoid", "sub-optimal visualizations downstream
                    (max. of 20 displayed below):",
                   "{.val {long_row_names}}"))
  }
 }

 invisible(x)

}



#' Metadata Label Column Validator
#'
#' Internal function for validating the sample label column in a metadata
#' `data.frame` for a `DGEList` object.
#'
#' @param metadata A `data.frame` containing sample metadata, where rows
#'   represent samples.
#' @param label_column The name or index of the column in `metadata` that contains
#'   labels used to identify each sample. The function checks if `label_column`:
#'   - Exists in `metadata`.
#'   - Contains non-missing and non-blank values.
#'   - Has unique values.
#'   - Adheres to R syntax rules (a warning is issued if not).
#'   Additionally, if `verbose = TRUE`, a warning is displayed when label values
#'   exceed 15 characters to encourage better plot formatting.
#' @param verbose Logical. If `TRUE`, checks whether the label values exceed
#'   15 characters. Default: `FALSE`.
#'
#' @rdname check_label_column
#' @keywords internal
#'
#' @return The validated column name of `label_column`, if all checks pass.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'
#' metadata <- data.frame(
#'   sample = c("C1", "C2", "T1", "T2"),
#'   sample_id = c("C100", "C202", "T303", "T100"),
#'   condition = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' metadata2 <- data.frame(
#'   sample = c("C1", "C2", "T1", "T2"),
#'   sample = c("C100", "C202", "T303", "T100"),
#'   condition = c("Control", "Control", "Treatment", "Treatment"),
#'   check.names = FALSE
#' )
#'
#' metadata3 <- data.frame(
#'   sample = c("C2", "C2", "T1", "T2"),
#'   sample_id = c("C100", "C202", NA, "T100"),
#'   sample_name = c("sample_1", "sample_2", "", "sample_4"),
#'   condition = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' metadata4 <- data.frame(
#'   long_sample_column_name = c("C1", "C2", "T1", "T2"),
#'   sample_id = c("C100", "C202", "chemo_treatment_303", "chemo_treatment_100"),
#'   condition = c("Control", "Control", "Treatment", "Treatment"),
#'   row.names = c("C100", "C202", "T303", "T100_long_row_name")
#' )
#'
#' check_label_column(metadata = metadata, label_column = "sample_id")
#' check_label_column(metadata = metadata, label_column = 2)
#' check_label_column(metadata = metadata, label_column = "replicates")
#' check_label_column(metadata = metadata2, label_column = "sample")
#' check_label_column(metadata = metadata3, label_column = 1)
#' check_label_column(metadata = metadata3, label_column = 2)
#' check_label_column(metadata = metadata3, label_column = "sample_name")
#' check_label_column(metadata = metadata4, label_column = "sample_id")
#'
#' }
#' }
#'
#'
check_label_column <- function(metadata,
                               label_column,
                               verbose = FALSE) {


 ## check metadata is a data.frame
 metadata <- check_metadata(x = metadata, verbose = FALSE)


 ## check label_column is one value
 is_single <- length(label_column) == 1L
 if(!is_single) {
  cli::cli_abort(c("x" = "Length of {.arg label_column} does not equal 1.",
                   "i" = "Only specify one column name or column
                     index for the {.arg label_column} argument."))
 }


 ## check label_column is character or numeric value
 if(all((!check_string(label_column)),
        (!check_num(label_column)))) {
  cli::cli_abort(c("x"="{.arg label_column}: {.val {label_column}}
                     is not a character string or numeric value.",
                   "i" = "The {.arg label_column} argument should be a
                     column number between {.val {1}} and
                     {.val {ncol(metadata)}} or", "one of the following
                     column names:", "{.val {colnames(metadata)}}."))
 }


 ## check label_column is column index
 if(check_num(label_column)) {
  is_index <- label_column %in% 1:ncol(metadata)
  if(!is_index) {
   cli::cli_abort(c("x" = "{.arg label_column}: {.val {label_column}} is
                       not a valid column index.",
                    "i" = "The {.arg label_column} argument should be
                       a column number between {.val {1}} and
                       {.val {ncol(metadata)}} or", "one of the following
                       column names:", "{.val {colnames(metadata)}}."))
  }


  ## if index get column name
  label_column <- colnames(metadata)[label_column]

 } ## number


 ## check label_column is a unique column name
 if(check_string(label_column)) {
  is_name <- label_column %in% colnames(metadata)
  if(!is_name) {
   cli::cli_abort(c("x"="The {.arg label_column} argument
                       {.val {label_column}} is not a column name in
                       {.arg metadata}.",
                    "i" = "Available columns include:",
                    "{.val {colnames(metadata)}}."))
  }


  ## check only once occurance of label_column in metadata
  unique_name <- grep(paste0("^",label_column,"$"), colnames(metadata))
  if(length(unique_name) != 1L) {
   cli::cli_abort(c("x" = "The {.arg label_column} argument
                       {.val {label_column}} is not a unique column name
                       in {.arg metadata}."))
  }
 }


 ## get vector of label values in label_column
 label_column <- as.character(label_column)
 label_values <- metadata[, label_column]


 ## check if any values are NA
 is_na <- is.na(label_values)
 if(any(is_na)) {
  cli::cli_abort(c("x" = "The {.val {label_column}} {.arg label_column}
                     in {.arg metadata} cannot have missing values."))
 }


 ## check if any values are blank
 is_blank <- nchar(label_values) == 0L
 if(any(is_blank)) {
  cli::cli_abort(c("x" = "The {.val {label_column}} {.arg label_column}
                     in {.arg metadata} cannot have blank values."))
 }


 ## check if any values are duplicated
 dup_values <- check_dup(x = label_values)
 if(!is.null(dup_values)) {
  num_dup_values <- length(dup_values)
  if(length(dup_values) > 10L) {
   dup_values <- dup_values[1:10]
  }
  cli::cli_abort(c("x" = "The {.val {label_column}} {.arg label_column}
                     in {.arg metadata} cannot have duplicate values.",
                   "x" = "A total of {.val {num_dup_values}} duplicates
                     were identified.",
                   "i" = "Duplicate values include (max. of 10 displayed below):",
                   "{.val {dup_values}}"))
 }


 ## warn if values do not follow R syntax rules
 is_correct_syntax <- check_syntax(x = label_values)
 if(!is_correct_syntax){
  cli::cli_inform(c("x" = "The values in the {.val {label_column}}
                     {.arg label_column} do not follow R syntax rules.",
                    "This may cause problems when matching between data.frames.",
                    "The functions {.fun base::make.names} and
                     {.fun base::make.unique} can be used to create",
                    "syntactically valid names."))
 }


 if(verbose){

  ## warn if label values are > 15 characters in length
  char_thresh <- 15
  long_label_values <- check_long(x = label_values,
                                  char_thresh = char_thresh)
  if(!is.null(long_label_values)) {
   if(length(long_label_values) > 20L) {
    long_label_values <- long_label_values[1:20]
   }
   cli::cli_warn(c("i" = "Some values in the {.val {label_column}}
                    {.arg label_column} of {.arg metadata} are","very long
                    (> {.val {char_thresh}} characters).",
                   "i" = "You may consider shortening the following labels
                    to avoid", "sub-optimal visualizations downstream
                    (max. of 20 displayed below):",
                   "{.val {long_label_values}}"))
  }
 }

 label_column

}



#' Metadata Grouping Column Validator
#'
#' Internal function for validating the grouping column in a metadata `data.frame`.
#'
#' @param metadata A `data.frame` containing sample metadata, where rows represent samples.
#' @param grouping_column The name or index of the column in `metadata` that provides
#'   information on how to group the samples. The function checks if `grouping_column`:
#'   - Exists in `metadata`.
#'   - Contains non-missing and non-blank values.
#'   - Has unique values (if required for analysis).
#'   Additionally, if `verbose = TRUE`, a warning is displayed when grouping values
#'   exceed 15 characters to encourage better formatting.
#' @param verbose Logical. If `TRUE`, checks whether values in `grouping_column`
#'   exceed 15 characters. Default: `FALSE`.
#'
#' @rdname check_grouping_column
#' @keywords internal
#'
#' @return The validated column name of `grouping_column`, if all checks pass.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'
#' metadata <- data.frame(
#'   sample = c("C1", "C2", "T1", "T2"),
#'   sample_id = c("C100", "C202", "T303", "T100"),
#'   condition = c("Control", "Control", "Treatment", "Treatment")
#' )
#'
#' metadata2 <- data.frame(
#'   sample = c("C1", "C2", "T1", "T2"),
#'   sample = c("C100", "C202", "T303", "T100"),
#'   condition = c("Control", "Control", "Treatment", "Treatment"),
#'   check.names = FALSE
#' )
#'
#' metadata3 <- data.frame(
#'   sample = c("C2", "C2", "T1", "T2"),
#'   sample_id = c("C100", "C202", NA, "T100"),
#'   sample_name = c("sample_1", "sample_2", "", "sample_4"),
#'   condition = c("Baseline_Saline_Control", "Baseline_Saline_Control", "Treatment", "Treatment"),
#'   treatment = c("", "Control", "Chemo", "Chemo"),
#'   group = c("Control", NA, "Treatment", ""),
#'   chemo = c("Control", "Saline", "Treatment", "Treatment"),
#'   check.names = FALSE
#' )
#'
#' check_grouping_column(metadata = metadata, grouping_column = "condition")
#' check_grouping_column(metadata = metadata, grouping_column = "replicates")
#' check_grouping_column(metadata = metadata, grouping_column = 7)
#' check_grouping_column(metadata = metadata2, grouping_column = 3)
#' check_grouping_column(metadata = metadata2, grouping_column = "sample")
#' check_grouping_column(metadata = metadata3, grouping_column = "condition")
#' check_grouping_column(metadata = metadata3, grouping_column = 5)
#' check_grouping_column(metadata = metadata3, grouping_column = "group")
#' check_grouping_column(metadata = metadata3, grouping_column = 7)
#'
#' }
#' }
#'

check_grouping_column <- function(metadata,
                                  grouping_column,
                                  verbose = FALSE) {


 ## chceck metadata is data.frame
 metadata <- check_metadata(x = metadata, verbose = FALSE)


 ## check grouping_column is not empty i.e. NULL
 if(is.null(grouping_column)){
  cli::cli_abort(c("x"="{.arg grouping_column} cannot be empty.",
                   "i" = "The {.arg grouping_column} argument should be a
                     column number between {.val {1}} and
                     {.val {ncol(metadata)}} or", "one of the following
                     column names:", "{.val {colnames(metadata)}}."))
 }


 ## check grouping_column is one value
 if(length(grouping_column) != 1L) {
  cli::cli_abort(c("x"="Length of {.arg grouping_column} does not equal 1.",
                   "i" = "Only specify one column name or column
                     index for the","{.arg grouping_column} argument."))
 }


 ## check grouping_column is character or numeric value
 # if(all(!check_string(x = grouping_column),
 #        !check_int(x = grouping_column))) {
 if(all(!is.character(grouping_column),
        !is.numeric(grouping_column))) {
  cli::cli_abort(c("x"="{.arg grouping_column}: {.val {grouping_column}}
                     is not a character string or numeric value.",
                   "i" = "The {.arg grouping_column} argument should be a
                     column number between {.val {1}} and
                     {.val {ncol(metadata)}} or", "one of the following
                     column names:", "{.val {colnames(metadata)}}."))
 }


 ## check if grouping_column is a column index
 if(is.numeric(grouping_column)) {
  is_index <- grouping_column %in% 1:ncol(metadata)
  if(!is_index) {
   cli::cli_abort(c("x" = "{.arg grouping_column}: {.val {grouping_column}}
                       is not a valid column index.",
                    "i" = "The {.arg grouping_column} argument should be
                       a column number between {.val {1}} and
                       {.val {ncol(metadata)}} or", "one of the following
                       column names:", "{.val {colnames(metadata)}}."))
  }


  ## if index get column name
  grouping_column <- colnames(metadata)[grouping_column]

 } ## number


 ## check grouping_column is a unique column name
 if(is.character(grouping_column)) {
  is_name <- grouping_column %in% colnames(metadata)
  if(!is_name) {
   cli::cli_abort(c("x"="The {.arg grouping_column} argument
                       {.val {grouping_column}} is not a column name in
                       {.arg metadata}.",
                    "i" = "Available columns include:",
                    "{.val {colnames(metadata)}}."))
  }


  ## check only once occurance of grouping_column in metadata
  unique_name <- grep(paste0("^",grouping_column,"$"), colnames(metadata))
  if(length(unique_name) != 1L) {
   cli::cli_abort(c("x" = "The {.arg grouping_column} argument
                       {.val {grouping_column}} is not a unique column name
                       in {.arg metadata}."))
  }
 }


 ## get vector of values in grouping_column
 grouping_column <- as.character(grouping_column)
 group_labels    <- as.character(metadata[, grouping_column])


 ## check if any values are NA
 is_na <- is.na(group_labels)
 if(any(is_na)) {
  cli::cli_abort(c("x" = "The {.val {grouping_column}} {.arg grouping_column}
                     in {.arg metadata} cannot have missing values."))
 }


 ## check if any values are blank
 is_blank <- nchar(group_labels) == 0L
 if(any(is_blank)) {
  cli::cli_abort(c("x" = "The {.val {grouping_column}} {.arg grouping_column}
                     in {.arg metadata} cannot have empty strings."))
 }


 ## check groups have at least 2 replicates
 num_per_group <- table(group_labels)
 is_dup <- num_per_group > 1L
 if(any(!is_dup)) {
  no_dup_values <- names(num_per_group)[!is_dup]
  cli::cli_warn(c("i" = "The {.val {grouping_column}} {.arg grouping_column}
                    in {.arg metadata} contains {.val {length(no_dup_values)}}
                    groups","that do not have at least 2 replicates.
                    Downstream quality control","and statistical analyses
                    typically require at least 2 replicates per group.",
                  "Single value groups include: {.val {no_dup_values}}"))
 }


 if(verbose){

  ## warn if grouping values are > 15 characters in length
  char_thresh <- 15
  long_group_labels <- check_long(x = group_labels,
                                  char_thresh = char_thresh)
  if(!is.null(long_group_labels)) {
   long_group_labels <- unique(long_group_labels)
   if(length(long_group_labels) > 10L) {
    long_group_labels <- long_group_labels[1:10]
   }
   cli::cli_warn(c("i" = "Some values in the {.val {grouping_column}}
                    {.arg grouping_column} of {.arg metadata} are",
                   "very long (> {.val {char_thresh}} characters).",
                   "i" = "You may consider shortening the following
                    group names to avoid", "sub-optimal visualizations
                    downstream (max. of 10 displayed below):",
                   "{.val {long_group_labels}}"))
  }

 }

 grouping_column

}


################################################################################
################################################################################
##                                                                            ##
##                            SUPPORT FUNCTIONS                               ##
##                                                                            ##
################################################################################
################################################################################



#' Check if an object is a data.frame or matrix
#'
#' Internal helper function to determine whether an object is a `data.frame` or `matrix`.
#'
#' @param x An object to be tested (e.g., counts, metadata, annotation).
#'
#' @return Logical. Returns `TRUE` if `x` is a `data.frame` or `matrix`, otherwise returns `FALSE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_data(data.frame(a = 1:3, b = 4:6)) # TRUE
#'   check_data(matrix(1:9, nrow = 3))        # TRUE
#'   check_data(list(a = 1, b = 2))           # FALSE
#' }
#' }
check_data <- function(x) {
  any(is.data.frame(x), is.matrix(x))
}


#' Check if a file exists
#'
#' Internal helper function to verify if an input string corresponds to an existing file.
#'
#' @param x A character string representing a file path.
#'
#' @return Logical. Returns `TRUE` if `x` is a valid file path and the file exists.
#' Returns `FALSE` if `x` is not a valid string or the file does not exist.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_file("path/to/existing_file.txt") # TRUE if file exists
#'   check_file("nonexistent_file.txt")      # FALSE
#' }
#' }

check_file <- function(x) {
  if (is.character(x) && length(x) == 1L && nchar(x) > 0L) {
    file.exists(x)
  } else {
    FALSE
  }
}

# check_file <-  function(x){
#
#  ## check is character string of length 1 that is not empty
#  if(all(is.character(x),
#         (length(x) == 1L),
#         (nchar(x) > 0L))){
#   ## if is string check if is a file that exists
#   ## if file exists return TRUE if not return FALSE
#   is_file <- file.exists(x)
#  } else {
#   ## if not string return false
#   is_file <- FALSE
#  }
#  is_file
#
# }



#' Check if a vector follows R syntax rules
#'
#' Internal helper function to verify whether values in a vector conform to R's variable naming rules.
#'
#' @param x A vector of values to check.
#'
#' @return Logical. Returns `TRUE` if all values in `x` follow R syntax rules, otherwise returns `FALSE`.
#'
#' @details
#' This function checks whether values in `x` comply with R's syntactic name rules using `make.names()`.
#' - A valid name consists of letters, numbers, dots (`.`), or underscores (`_`).
#' - It must start with a letter or a dot not followed by a number.
#' - See: `base::make.names()` and `base::make.unique()`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_syntax(c("C1", "C2", "T1", "T2"))       # TRUE
#'   check_syntax(c("C 100", "C-202", "T 303"))    # FALSE
#'   check_syntax(c("Var_1", "Var.2", ".Valid"))   # TRUE
#'   check_syntax(c(1, 2, 3, 4))                   # TRUE (coerced to strings)
#' }
#' }
check_syntax <- function(x) {
  if (!is.vector(x) || length(x) < 1L || is.null(x)) {
    cli::cli_abort("The {.arg x} argument must be a vector of length >= 1.")
  }

  x_test <- base::make.names(names = x, unique = FALSE, allow_ = TRUE)
  identical(x, x_test)
}


#' Identify Positive Integer Values
#'
#' Internal helper function to check if a value is a positive whole number (including zero).
#'
#' @param x A single value to be tested. Must be of length 1L.
#'
#' @return Logical. Returns `TRUE` if `x` is a positive integer (including 0), otherwise returns `FALSE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_int(5)      # TRUE
#'   check_int(0)      # TRUE
#'   check_int(-3)     # FALSE
#'   check_int(4.5)    # FALSE
#'   check_int("ABC")  # FALSE
#' }
#' }
check_int <- function(x) {
  if (length(x) != 1L) {
    cli::cli_abort("{.arg x} must have length 1")
  }
  grepl("(^[1-9][0-9]*$|^0$)", x)
}



#' Identify Positive Numeric Values
#'
#' Internal helper function to verify if `x` is a valid positive numeric value (≥ 0).
#'
#' @param x A single value to be tested.
#'
#' @return Logical. Returns `TRUE` if `x` is a valid positive numeric value (≥ 0).
#' Returns `FALSE` for non-numeric inputs, characters, multiple values, empty strings, white spaces,
#' `NA`, `NULL`, `Inf`, or negative numbers.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_num(1.4323)   # TRUE
#'   check_num(0)        # TRUE
#'   check_num(-45)      # FALSE
#'   check_num(Inf)      # FALSE
#'   check_num("snap")   # FALSE
#'   check_num("  ")     # FALSE
#'   check_num(NA)       # FALSE
#'   check_num(1e6)      # TRUE
#' }
#' }
check_num <- function(x) {
  all(is.numeric(x),
      !is.character(x),
      !(grepl("^\\s*$", x)),  # Checks for blank/whitespace/tab/newlines
      length(x) == 1L,
      nchar(x) > 0L,  # Ensures x is not an empty string
      !is.na(x),
      !is.null(x),
      is.finite(x),
      !is.infinite(x),
      x >= 0)
}


#' Validate Single Character String
#'
#' Internal helper function to check whether `x` is a valid character string of length 1.
#'
#' @param x A value to be checked.
#'
#' @return Logical. Returns `TRUE` if `x` is a character string of length 1.
#' Returns `FALSE` if `x` is numeric, a vector (length > 1L), an empty string, white spaces,
#' `NA`, `NULL`, `Inf`, etc.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_string("ProjectID_122522") # TRUE
#'   check_string(123)                # FALSE
#'   check_string(NA)                 # FALSE
#'   check_string(NULL)               # FALSE
#'   check_string(c("oh", "snap"))    # FALSE
#'   check_string("     ")            # FALSE
#'   check_string("")                 # FALSE
#'   check_string(Inf)                # FALSE
#' }
#' }
check_string <- function(x) {
  all(!is.numeric(x),
      is.character(x),
      length(x) == 1L,
      nchar(x) > 0L,
      !(grepl("^\\s*$", x)),  # Checks for empty/whitespace-only strings
      !is.na(x),
      !is.null(x))
}



#' Check if a Value is a TRUE or FALSE Logical of Length 1
#'
#' Internal helper function to verify if `x` is a logical value (TRUE or FALSE) of length 1.
#'
#' @param x A value to be checked.
#'
#' @return Logical. Returns `TRUE` if `x` is a logical value (`TRUE` or `FALSE`) of length 1.
#' Returns `FALSE` if `x` is any of the following: numeric, character string, length greater than 1,
#' empty string, a series of blank spaces, `NA`, `NULL`, or non-logical values.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check_logical(TRUE)    # TRUE
#'   check_logical(FALSE)   # TRUE
#'   check_logical(1)       # FALSE
#'   check_logical("TRUE")  # FALSE
#'   check_logical(c(TRUE, FALSE)) # FALSE
#'   check_logical(NA)      # FALSE
#' }
#' }
check_logical <- function(x) {
  is_logical <- all(is.logical(x),
                    length(x) == 1L,
                    nchar(x) > 0L,
                    !(grepl("^\\s*$", x)),
                    !is.character(x),
                    !is.na(x),
                    !is.null(x),
                    !is.numeric(x))
  is_logical
}


#' Check if a Vector Includes a Set of Reference Values
#'
#' This function checks if a reference vector (`ref`) is a subset of the input vector (`x`).
#' The function ensures that `x` has the same or greater length than `ref` and that no duplicates exist in `ref`.
#' Both `x` and `ref` must follow R's syntax rules for valid variable names (e.g., column names and row names).
#'
#' See \code{\link[base]{make.names}} for details on syntactic names.
#'
#' @param x A vector of values to be checked. The length of `x` must be equal to or greater than `ref` with no duplicates.
#'   This is typically a metadata or annotation column.
#' @param ref A vector of allowed reference values. All values in `ref` must be unique, and it must contain no `NA` or empty values.
#'   This is typically column names or row names in a counts matrix.
#'
#' @return Logical. Returns `TRUE` if all values in `ref` are found in `x`, otherwise returns `FALSE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # ref is a subset of x
#'   x <- paste0("gene_", 1:10)
#'   ref <- paste0("gene_", 3:7)
#'   check_vec(x = x, ref = ref)
#'
#'   # x does not include all ref values
#'   x <- paste0("gene_", 1:5)
#'   ref <- paste0("gene_", 4:6)
#'   check_vec(x = x, ref = ref)
#'
#'   # x does not include any ref values
#'   x <- paste0("gene_", 1:5)
#'   ref <- paste0("gene_", 20:22)
#'   check_vec(x = x, ref = ref)
#'
#'   # x contains duplicates, NA, empty strings, or blank spaces
#'   ref <- c("gene_", "NA", " ")
#'   check_vec(x = x, ref = ref)
#'
#' ## x is < ref
#' ## x contains duplicates
#' ## x contains NA
#' ## ref includes duplicates
#' ## ref includes NA
#' ## ref includes empty strings
#' ## ref includes blank spaces
#'
#' ## ref contains values equal to blank spaces
#'
#'
#' samples <- data.frame(sample_id = c("C100","C202","T303","T100"),
#'                       sample    = c("C1","C2","T1","T2"),
#'                       group     = c("Con","Con","Treat","Treat"))
#'
#' counts <- data.frame(C100 = c(0,10,21,3,4),
#'                      C202 = c(0,0,0,2,1),
#'                      T100 = c(15,33,21,55,42),
#'                      row.names = c(paste0("gene_",1:5)))
#'
#' genes <- data.frame(gene_id   = c(paste0("gene_",1:5),"NA.","NA..1",""),
#'                     gene_name = LETTERS[1:8])
#'
#'
#' is_id <- unlist(lapply(samples,
#'                        function(x) {
#'                          check_vec(x = x, ref = colnames(counts))
#'                      }))
#'
#' is_id <- sapply(as.list.data.frame(samples),
#'                 function(x) {
#'                   check_vec(x = x, ref = colnames(counts))
#'                 })
#'
#' id_column <- colnames(samples)[is_id][1]
#'
#' sort(samples[, id_column]); sort(colnames(counts))
#'
#'
#' ## find id column in annotation matching counts row names
#' is_id <- unlist(lapply(genes,
#'                        function(x) {
#'                          check_vec(x = x, ref = row.names(counts))
#'                      }))
#'
#'
#'  }
#' }
check_vec <- function(x, ref){

 ## check arguments

 ## check if ref or x are factors if so change to character
 if(is.factor(ref)) { ref <- as.character(ref) }
 if(is.factor(x)) { x <- as.character(x) }

 ## check ref is a vector of length >= 1
 ref_vec <- any(!is.vector(ref), !(length(ref) >= 1L))

 ## check if ref has any missing values
 ref_na <- any(is.na(ref))

 ## check if ref has any empty strings
 ref_blank <- any(nchar(ref) == 0L)

 ## check if ref has values equal to blank spaces
 ref_spaces <- any(grepl("^\\s*$", ref))

 ## check if ref has duplicate values
 ref_dup <- !is.null(check_dup(x = ref))

 ## check if x is a vector of length >= 1
 x_vec <- any(!is.vector(x), !(length(x) >= 1L))

 ## check if x length is < ref
 x_short <- length(x) < length(ref)

 if(any(ref_vec,
        ref_na,
        ref_blank,
        ref_spaces,
        ref_dup,
        x_vec,
        x_short
        # x_dup
 )) {

  is_match <- FALSE

 } else {

  ## check x contains all values in ref
  is_match <- all(ref %in% x)
  if(is_match){

   ## if all ref values are in x
   ## get subset of x that match the reference
   m <- match(ref, x)
   x_in_ref <- x[m]
   x_not_in_ref <- x[!x %in% x_in_ref]

   ## then check if any values in x subset are duplicated
   x_dup <- length(check_dup(x = x_in_ref)) > 0L

   if(!x_dup){
    is_match <- TRUE
   } else {
    is_match <- FALSE
   }
  }
 }

 is_match

}


#' Identify ID Column
#'
#' This function identifies the first column in a data frame `x` that contains
#' all the values specified in the `ref` vector. The `ref` values should typically
#' represent sample or gene IDs, and the function checks for the presence of these
#' reference values in the columns of `x`.
#'
#' @param x A data frame. The columns in `x` will be queried against the reference
#'   values provided in `ref`. `x` is typically a data frame of metadata or annotation.
#' @param ref A character vector of reference values. The values in `ref` must be
#'   unique, and should not contain any `NA` or blank values. The values in `ref` are
#'   typically the column or row names of a counts matrix (e.g., sample IDs or gene IDs).
#'
#' @return The function returns the name of the first column in `x` that contains
#'   all the reference values. If no matching column is found, the function will return an error.
#'
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#' metadata <- data.frame(sample = c("C1","C2","T1","T2"),
#'                        sample_id = c("C 100","C 202","T 303","T 100"),
#'                        condition = c("Control","Control","Treatment",
#'                                      "Treatment"))
#'
#' counts <- data.frame(C.100 = c(0,10,21,3,4),
#'                      C.202 = c(0,0,0,2,1),
#'                      T.100 = c(15,33,21,55,42),
#'                      T.303 = c(12,22,54,33,33),
#'                      row.names = paste0("gene_",1:5))
#'
#' annotation <- data.frame(gene_id = c("gene_2","gene_3","gene_4",
#'                                      "gene_1","gene_5"),
#'                          gene_name = c("SNP","CRC","POP","RCE","CSP"),
#'                          description = c("Snap","crackle","pop","rice",
#'                                          "crispies"))
#'
#' ## rsyntax error
#' find_id_column(x = metadata, ref = colnames(counts))
#' find_id_column(x = annotation, ref = row.names(counts))
#'
#'   }
#' }
#'
find_id_column <- function(x, ref) {

 ## check ref arguments

 ## check ref is a vector of at least 1
 if (any(!is.vector(ref),
         !(length(ref) >= 1L))) {
  cli::cli_abort(c("The {.arg ref} argument must be a vector of
                     length >= 1."))
 }


 ## check if there are any missing values
 is_na <- is.na(ref)
 if (any(is_na)) {
  cli::cli_abort("The {.arg ref} argument cannot include missing values.")
 }

 ## check if any values are blank
 is_blank <- nchar(ref) == 0L
 if (any(is_blank)) {
  cli::cli_abort("The {.arg ref} argument cannot contain blank values.")
 }

 ## check ref for duplicate values
 dup_ref <- check_dup(x = ref)
 if (!is.null(dup_ref)) {
  if (length(dup_ref) > 20L) {
   dup_ref <- dup_ref[1:20]
  }
  cli::cli_abort(c("The {.arg ref} argument cannot contain duplicate values.",
                   "x" = "The following duplicate values were identified
                     (max. of 20 displayed below): {.val {dup_ref}}"))
 }



 ## check data argument

 ## check data is a data.frame
 if (!check_data(x = x)) {
  cli::cli_abort("The {.arg x} argument is not a data.frame or matrix.")
 }


 ## check data for duplicate column names
 dup_columns <- check_dup(x = colnames(x))
 if (!is.null(dup_columns)) {
  cli::cli_abort(c("The {.arg x} data.frame cannot contain duplicate
                     column names.",
                   "x" = "The following duplicate column names were
                     identified {.val {dup_columns}}."))
 }


 ## check number of rows in data is less than the number of ref
 if (nrow(x) < length(ref)) {
  cli::cli_abort("The number of rows in {.arg x} (n = {.val {nrow(x)}})
                   is less than the number of {.arg ref} values
                   (n = {.val {length(ref)}}).")
 }

 ## identify data column containing all ref
 check_id_column <- unlist(
  lapply(x,
         function(xx) {
          check_vec(x = xx, ref = ref)
         })
 )

 ## get matching column names
 matching_id_column <- colnames(x)[check_id_column]

 ## error if no column identified
 if (length(matching_id_column) == 0L) {
  cli::cli_abort(c("A column in the {.arg x} data.frame with values
                     matching the values in {.arg ref} could not
                     be identified."))
 }

 ## return column name of first matching id column
 matching_id_column[1]

}



#' Identify Duplicate Values
#'
#' This function identifies duplicate values in a vector. It can be used to find duplicate
#' column names, row names, or values in any particular column of a data frame.
#'
#' @param x A vector of values to be tested. This is typically a vector of column names,
#'   row names, or values in a particular column of a data frame.
#'
#' @return The function returns a vector of duplicate values. If no duplicates are found,
#'   it returns `NULL`. This allows for custom error handling if desired.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'   x <- c("A", "A", "B", "C", "D", "D", "D", "E")
#'   y <- c(1, 2, 3, "A", "B", "C")
#' ## returns vector of duplicates identified
#'   dup_present <- check_dup(x = x)
#' ## returns NULL if all unique
#'   dup_not_present <- check_dup(x = y)
#'
#' }
#' }

check_dup <- function(x) {

 ## check arguments: vector with at least 1 value and is not NULL
 if(any(!is.vector(x),
        !(length(x) >= 1L),
        is.null(x))) {
  cli::cli_abort(c("The {.arg x} argument must be a vector of values
                     with length >= 1."))
 }

 ## returns a vector of duplicates identified or NULL if all are unique
 is_dup <- duplicated(x)
 if(any(is_dup)) {
  dup_vals <- x[is_dup]
 } else {
  dup_vals <- NULL
 }
 unique(dup_vals)

}



#' Identify Values Above a Specified Character Length
#'
#' This function identifies values in a vector that exceed a specified maximum number
#' of characters. This is often used to identify long sample names or group labels that
#' could cause issues in visualizations or presentations.
#'
#' @param x A vector of values. Typically column names, row names, or metadata values.
#' @param char_thresh A non-negative numeric value specifying the maximum allowed
#'   character length. Values exceeding this length are returned.
#' @param na.rm Logical; if `TRUE`, `NA` values are removed before checking. If `FALSE`,
#'   `NA` and `NaN` may appear in the output.
#'
#' @return A vector of values in `x` whose character lengths exceed `char_thresh`.
#'   Returns `NULL` if no values meet the criterion.
#'
#' @details
#' Behavior under different combinations of `x` contents and `na.rm`:
#'
#' * **No NAs in `x`**  
#'   - `na.rm = TRUE` or `FALSE`, and **all values are below** the threshold → returns `NULL`.  
#'   - `na.rm = TRUE` or `FALSE`, and **some values exceed** the threshold → returns those values.
#'
#' * **NAs present in `x`**  
#'   - `na.rm = TRUE`: NAs are removed; behavior is based only on non-NA values.  
#'   - `na.rm = FALSE`:  
#'     - If all non-NA values are below threshold → returns `c(NA, NaN)` if present.  
#'     - If some exceed threshold → returns the long values and any `NA`/`NaN` values.
#'
#' @examples
#' \dontrun{
#' x <- c("Con_01", "Con_02", NA, "TAC_Treatment_01", "", Inf, NaN, "TAC_Treat_02")
#' y <- c("Con_01", "Con_02", "Treat_01", "Treat_101")
#' z <- c("Con_01", NA, "TAC_Treatment_01", "TAC_Treat_02", "", Inf)
#'
#' # identify values with more than 15 characters
#' check_long(x = x, char_thresh = 15, na.rm = FALSE)
#'
#' # identify values with more than 8 characters
#' check_long(x = y, char_thresh = 8)
#'
#' # identify values with more than 6 characters
#' check_long(x = y, char_thresh = 6)
#' }

check_long <- function(x, char_thresh = 15, na.rm = FALSE) {

 ## check x argument is NULL
 if(is.null(x)){
  cli::cli_abort("{.arg x} cannot be NULL.")
 }

 ## check x argument is NA or a vector of NAs
 if(all(is.na(x))){
  cli::cli_abort("{.arg x} cannot be a vector of NAs.")
 }

 ## check na.rm argument is TRUE or FALSE
 if(!check_logical(na.rm)){
  cli::cli_abort("{.arg na.rm} must be one of TRUE or FALSE.")
 }

 ## if na.rm = TRUE remove NA NaN values from x
 if(na.rm){
  x <- x[!is.na(x)]
 }

 ## check x is a vector of values with length >= 1
 if(any(!is.vector(x),
        !(length(x) >= 1L))) {
  cli::cli_abort(c("The {.arg x} argument must be a vector with
                     length >= 1."))
 }

 ## check char_thresh argument is a single value
 if (!(length(char_thresh) == 1L)){
  cli::cli_abort(c("The {.arg char_thresh} must be an integer value of length 1
                     that is greater >= 0."))
 }

 ## check char_thresh argument is a number >= 0
 if (!check_int(x = char_thresh)){
  cli::cli_abort(c("The {.arg char_thresh} must be an integer value that is
                     greater >= 0."))
 }

 is_long <- nchar(x) > char_thresh
 is_long <- is_long | is.na(is_long)
 if(any(is_long)) {
  long_vals <- x[is_long]
 } else {
  long_vals <- NULL
 }

 long_vals

}


