#' Metadata validator.
#'
#' Internal function for validating the structure of an input metadata
#' data.frame.
#'
#' @param x a data.frame of sample metadata. Rows are samples.
#' @param verbose logical. Should the row names be checked for
#'   character length? Default: FALSE.
#'
#' @rdname check_metadata
#'
#' @keywords internal

#' @return If all checks pass the function returns the input metadata as a
#'   data.frame (invisible).
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
#' @export
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






#' Metadata label_column validator.
#'
#' Internal function for validating metadata sample label_column information
#' for a DGEList object.
#'
#' @param metadata a data.frame of sample metadata. Rows are samples.
#' @param label_column The name of column or column index within
#'   the metadata data.frame which contains labels to use to identify each sample.
#'   To facilitate good plot formatting, when verbose = TRUE a warning message
#'   is invoked if sample label values have more than 15 characters.
#'   The function will give an error if the /code{label_column} is not a column
#'   in /code{metadata}, label values have missing values, blanks, or are not
#'   unique. A warning is displayed if the label values do not follow R syntax
#'   rules.
#' @param verbose logical. Should the values of the label_column be checked for
#'   character length? Default: False
#'
#' @rdname check_label_column
#'
#' @keywords internal
#'
#' @return If all checks pass the function returns the /code{label_column}
#'   column name.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#' metadata <- data.frame(sample = c("C1","C2","T1","T2"),
#'                        sample_id = c("C100","C202","T303","T100"),
#'                        condition = c("Control","Control","Treatment",
#'                                      "Treatment"))
#'
#' metadata2 <- data.frame(sample = c("C1","C2","T1","T2"),
#'                         sample = c("C100","C202","T303","T100"),
#'                         condition = c("Control","Control","Treatment",
#'                                       "Treatment"), check.names = FALSE)
#'
#' metadata3 <- data.frame(sample = c("C2","C2","T1","T2"),
#'                         sample_id = c("C100","C202",NA,"T100"),
#'                         sample_name = c("sample_1","sample_2","","sample_4"),
#'                         condition = c("Control","Control","Treatment",
#'                                       "Treatment"))
#'
#' metadata4 <- data.frame(long_sample_column_name = c("C1","C2","T1","T2"),
#'                         sample_id = c("C100","C202","chemo_treatment_303",
#'                         "chemo_treatment_100"),
#'                         condition = c("Control","Control","Treatment",
#'                                       "Treatment"),
#'                         row.names = c("C100","C202","T303","T100_long_row_name"))
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
#'   }
#' }
#'
#' @export
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



#' metadata group column validator.
#'
#' Internal function for validating metadata grouping column information
#'
#' @param metadata a data.frame of sample metadata. Rows are samples.
#' @param grouping_column The name of column or column index within
#'   the metadata data.frame which gives information on how to group the samples.
#' @param verbose logical. Should the grouping_column be checked for
#'   character length? Default: False
#'
#' @rdname check_grouping_column
#'
#' @keywords internal
#'
#' @return If all checks pass the function returns the /code{grouping_column}
#'   column name.
#'
#' @examples
#'
#' \dontrun{
#'
#' metadata <- data.frame(sample = c("C1","C2","T1","T2"),
#'                        sample_id = c("C100","C202","T303","T100"),
#'                        condition = c("Control","Control","Treatment","Treatment"))
#'
#' metadata2 <- data.frame(sample = c("C1","C2","T1","T2"),
#'                         sample = c("C100","C202","T303","T100"),
#'                         condition = c("Control","Control","Treatment","Treatment"),
#'                         check.names = FALSE)
#'
#' metadata3 <- data.frame(sample = c("C2","C2","T1","T2"),
#'                         sample_id = c("C100","C202",NA,"T100"),
#'                         sample_name = c("sample_1","sample_2","","sample_4"),
#'                         condition = c("Baseline_Saline_Control",
#'                                       "Baseline_Saline_Control",
#'                                       "Treatment","Treatment"),
#'                         treatment = c("","Control","Chemo","Chemo"),
#'                         group = c("Control",NA,"Treatment",""),
#'                         chemo = c("Control","Saline","Treatment","Treatment"),
#'                         check.names = FALSE)
#'
#'
#' check_grouping_column(metadata = metadata,  grouping_column = "condition")
#' check_grouping_column(metadata = metadata,  grouping_column = "replicates")
#' check_grouping_column(metadata = metadata,  grouping_column = 7)
#' check_grouping_column(metadata = metadata2, grouping_column = 3)
#' check_grouping_column(metadata = metadata2, grouping_column = "sample")
#' check_grouping_column(metadata = metadata3, grouping_column = "condition")
#' check_grouping_column(metadata = metadata3, grouping_column = 5)
#' check_grouping_column(metadata = metadata3, grouping_column = "group")
#' check_grouping_column(metadata = metadata4, grouping_column = 7)
#'
#' }
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



#' check if x is a data.frame or matrix
#'
#' for internal use. helper function. NEED TO UPDATE
#'
#' @param x an object to be tested. (e.g. counts, metadata, annotation)
#'
#' @return The function returns TRUE if x is a data.frame or matrix
#' If x is not a data.frame or matrix the function returns FALSE.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#'
check_data <- function(x){
 
 if(any(is.data.frame(x),
        is.matrix(x))){
  is_data <- TRUE
 } else {
  is_data <- FALSE
 }
 is_data
 
}



#' check if x is a file that exists
#'
#' for internal use. helper function. NEED TO UPDATE
#'
#'
#' @param x an object to be tested.
#'
#' @return if x is a character string that points to
#' a file that exists the function returns TRUE. if x is not a character string
#' of length 1 that is not empty or is not a valid file name the function
#' returns FALSE
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#' # no examples yet
#'  }
#' }
#'
check_file <-  function(x){
 
 ## check is character string of length 1 that is not empty
 if(all(is.character(x),
        (length(x) == 1L),
        (nchar(x) > 0L))){
  ## if is string check if is a file that exists
  ## if file exists return TRUE if not return FALSE
  is_file <- file.exists(x)
 } else {
  ## if not string return false
  is_file <- FALSE
 }
 is_file
 
}



#' check if a vector of values follow R syntax rules.
#'
#' for internal use. helper function. NEED TO UPDATE
#'
#' NEED TO UPDATE
#' The input vector is compared to a test vector created using
#' base::make.names() function (duplicate values allowed).
#' Note: A syntactically valid name consists of letters, numbers and the dot or
#' underline characters and starts with a letter or the dot not followed
#' by a number. see: base::make.names and base::make.unique.
#'
#' @param x vector of values to check.
#'
#' @return If the input and test vectors are identical the function returns TRUE
#'   and FALSE otherwise.
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' check_syntax(x = c("C1","C2","T1","T2"))
#' check_syntax(x = c("C 100","C 202","T 303","T 100"))
#' check_syntax(x = c("Treat-A","Treat-A","Treat-B","Treat-B"))
#' check_syntax(x = c(1, 2, 1, 2))
#'
#' metadata <- data.frame(sample = c("C1","C2","T1","T2"),
#'                        sample_id = c("C 100","C 202","T 303","T 100"),
#'                        condition = c("Treat-A","Treat-A","Treat-B","Treat-B"),
#'                        batch = c(1, 2, 1, 2))
#'
#' sapply(metadata, function(x){ check_syntax(x = x)})
#' apply(X = metadata, MARGIN = 2, FUN = check_syntax)
#'
#'
#' tmp <- lapply(metadata, function(i){
#'   ## if all values are numbers then return original values
#'   if(all(sapply(i, FUN = check_num))){
#'     values <- i
#'     cli::cli_inform(c("v"="all numeric. no check. returning
#'                                     ORIGINAL values."))
#'   } else {
#'     ## if any values are non-numeric check R syntax
#'     ## if vector does not follow R syntax return converted values
#'     if(!check_syntax(x = i)){
#'       values <- make.names(names = i, unique = FALSE, allow_ = TRUE)
#'       cli::cli_inform(c("x"="failed. returning NEW values."))
#'     } else {
#'       ## if any values are non-numeric check R syntax
#'       ## if vector follows R syntax return original values
#'       values <- i
#'       cli::cli_inform(c("v"="passed. returning ORIGINAL values."))
#'     }
#'   }
#'   values
#' })
#' metadata_syntax <- data.frame(do.call("cbind", tmp))
#'
#'   }
#' }
#'
check_syntax <- function(x){
 
 ## check arguments: vector with at least 1 value and is not NULL
 if (any((!is.vector(x)),
         (!(length(x) >= 1L)),
         is.null(x))) {
  cli::cli_abort(c("The {.arg x} argument must be a vector of values
                     with length >= 1."))
 }
 
 ## test vector: values converted to R syntax using make.names function
 ## duplicate values allowed
 x_test <- base::make.names(names = x, unique = FALSE, allow_ = TRUE)
 
 ## if test vector and input vector are the same then
 ## input vector follows R syntax rules.
 if (identical(x, x_test)){
  is_syntax <- TRUE
 } else {
  is_syntax <- FALSE
 }
 is_syntax
 
}



#' function identifies positive integer values
#'
#' this function does some stuff NEED TO UPDATW
#'
#' @param x description x a value to be tested. x must be single value of
#' length = 1L The function identifies positive whole numbers including 0
#' https://www3.ntu.edu.sg/home/ehchua/programming/howto/Regexe.html
#'
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' df <- data.frame(c1 = c(1, "ABC", Inf),
#'                  c2 = c(NA, 4.3, 6))
#'
#' sapply(unlist(df), check_int)
#'
#' mat <- matrix(c(1, NA, 2, Inf, -2, 10),
#'               nrow = 3, ncol = 2, byrow = TRUE,
#'               dimnames = list(c(paste0("gene_",1:3)), c("C1","C2")))
#'
#' sapply(c(mat), check_int)
#'
#' c(mat)[sapply(c(mat), check_int) == FALSE]
#'
#' }
#'}
#'
#'
check_int <- function(x){
 
 ## check arg
 # Check if an argument is a single numeric value
 if(any(!(length(x) == 1L))) {
  cli::cli_abort("{.arg x} must have length 1")
 }
 is_int <- grepl("(^[1-9][0-9]*$|^0$)", x)
 # test_int(x = x, na.ok = FALSE, lower = 0, null.ok = FALSE)
 is_int
 
 # is_integerlike
 # Numerical tolerance used to check whether a double or complex can be
 # converted to an integer. Default is sqrt(.Machine$double.eps).
 # abs(x - round(x, 0)) > sqrt(.Machine$double.eps)
 
}




#' function identifies positive numeric values
#'
#' this function is used to check that the norm.factors are valid numeric values
#'
#' @param x a value to be tested. x must be value of length = 1L
#'
#' @return the function returns TRUE if x is a positive numeric value
#'   greater than or equal to zero. The function returns FALSE if x is
#'   any of the following: non-numeric, a character string, vector with
#'   multiple values (length > 1L), an empty string (blank), white space
#'   e.g. blank spaces, tabs, new lines), NA, NULL, Inf, or a number less
#'    than zero i.e. negative number etc.
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' check_num(x = "snap")
#' check_num(x = 1.4323)
#' check_num(x = c(12,130))
#' check_num(x = NA)
#' check_num(x = -45)
#' check_num(x = Inf)
#' check_num(x = 0)
#' check_num(x = "    ")
#' check_num(x = " \n")
#' check_num(x = "")
#' check_num(x = 1e6)
#' check_num(x = 1/12)
#'
#'  }
#' }
#'
#'
check_num <- function(x){
 is_num <- all(is.numeric(x),
               !is.character(x),
               !(grepl("^\\s*$", x)), ## blank/whitespace/tab/newlines
               length(x) == 1L,
               nchar(x) > 0L, ## not empty string
               !is.na(x),
               !is.null(x),
               is.finite(x),
               !is.infinite(x),
               x >= 0)
 is_num
 
}




#' checks that x is a character string of length 1.
#'
#' for internal use. helper function. NEED TO UPDATE
#'
#'
#' @param x a value to be checked.
#'
#'
#' @returns the function returns TRUE if x is a character string of length 1.
#' the function returns FALSE if x is any of the following numeric,
#' has length > 1L, is an empty string, series of blank spaces,
#' NA, NULL, Inf etc.
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' check_string(x = "ProjectID_122522")
#' check_string(x = 123)
#' check_string(x = NA)
#' check_string(x = NULL)
#' check_string(x = c("oh","snap"))
#' check_string(x = TRUE)
#' check_string(x = "     ")
#' check_string(x = "")
#' check_string(x = Inf)
#'
#'  }
#'}
#'
#'
check_string <- function(x){
 
 ## check that project_id is a character string of length 1.
 ## project_id cannot be numeric, an emtpy string, vector of values,
 ##  NA or NULL
 ## "^\\s*$" asks for 0 or more (*) spaces (\\s) between
 ## beginning (^) and end ($) of string.
 is_string <- all(!is.numeric(x),
                  is.character(x),
                  length(x) == 1L,
                  nchar(x) > 0L,
                  !(grepl("^\\s*$", x)),
                  !is.na(x),
                  !is.null(x))
 is_string
 
}



#' checks that x is a TRUE or FALSE logical of length 1.
#'
#' for internal use. helper function. NEED TO UPDATE
#'
#'
#' @param x a value to be checked.
#'
#'
#' @returns the function returns TRUE if x is a TRUE OR FALSE logical of length 1.
#' the function returns FALSE if x is any of the following numeric,
#' has length > 1L, is an empty string, a series of blank spaces,
#' NA, NULL, or a character string.
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' check_logical(x = "ProjectID_122522")
#' check_logical(x = 123)
#' check_logical(x = 1.2342)
#' check_logical(x = NA)
#' check_logical(x = NULL)
#' check_logical(x = c("oh","snap"))
#' check_logical(x = TRUE)
#' check_logical(x = FALSE)
#' check_logical(x = c(TRUE,TRUE))
#' check_logical(x = "TRUE")
#' check_logical(x = "     ")
#' check_logical(x = "")
#' check_logical(x = Inf)
#'
#' }
#'}
#'
#'
check_logical <- function(x){
 
 ## "^\\s*$" asks for 0 or more (*) spaces (\\s) between
 ## beginning (^) and end ($) of string.
 is_logical <- all(is.logical(x),
                   (length(x) == 1L),
                   (nchar(x) > 0L),
                   !(grepl("^\\s*$", x)),
                   !is.character(x),
                   !is.na(x),
                   !is.null(x),
                   !is.numeric(x))
 is_logical
 
}




#' Check if a vector x includes a set of reference values.
#'
#' General purpose function that checks if a reference vector is a
#' subset of a given input vector x.
#'
#' x should have the same or longer length as the number of values
#' in ref with no duplicated values. The values in x and ref must
#' follow R's syntax rules for variable names (i.e. column names
#' and row names).
#'
#' See \code{\link{base::make.names}} for naming requirements.
#'
#' @param x a vector of values to be checked. The length of x
#'   must be equal to or greater than the ref vector with
#'   no duplicate values matching the reference. This parameter
#'   is typically the values from a particular metadata
#'   or annotation column.
#' @param ref a vector of allowed values. All values in the ref
#'   set must be unique with no NA or empty values. This parameter
#'   is typically the column names or row names of a counts matrix.
#'
#' @rdname check_vec
#'
#'
#' @return The function returns TRUE if the check is successful
#'   and FALSE if the check is unsuccessful.
#'
#' @examples
#'
#' \dontrun{
#' if(interactive()){
#'
#' ## ref is a subset of x
#' x   <- paste0("gene_", 1:10)
#' ref <- paste0("gene_", 3:7)
#' check_vec(x = x, ref = ref)
#'
#' ## x does not include all ref values
#' x   <- paste0("gene_", 1:5)
#' ref <- paste0("gene_", 4:6)
#' check_vec(x = x, ref = ref)
#'
#' ## x does not include any ref values
#' x   <- paste0("gene_", 1:5)
#' ref <- paste0("gene_", 20:22)
#' check_vec(x = x, ref = ref)
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





#' Identify id column.
#'
#' NEED TO UPDATE
#'
#' @param x a data.frame. The values in each column will be queried against
#' the reference. x is typically a data.frame of metadata or annotation.
#' @param ref a character vector of reference values. The values in ref must be
#' unique with no NAs or blanks. ref values are typically the column names or
#' row names of a counts matrix (sample ids / gene ids) with x defined as
#' metadata or annotation.
#'
#' @return The function returns the first column name in `x` that contains all
#' reference values.
#'
#' @examples
#'
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






#' Identify duplicate values.
#'
#' General purpose function that identifies duplicate values in a vector.
#'
#' @param x a vector of values to be tested. This parameter is typically a vector of
#' column names, row names, or vector of values in a particular column.
#'
#' @return The function returns a vector of duplicate values.
#' If no duplicates are present the function returns NULL to allow
#' follow-up custom error messaging.

#'
#' @examples
#' \dontrun{
#'   if(interactive()){
#'
#' x <- c("A", "A", "B", "C", "D", "D", "D","E")
#' y <- c(1, 2, 3, "A", "B", "C")
#'
#' ## returns vector of duplicates identified
#' dup_present <- check_dup(x = x)
#' dup_present
#'
#' ## returns NULL if all unique
#' dup_not_present <- check_dup(x = y)
#' dup_not_present
#'
#'   }
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





#' Identify values above a specified character length.
#'
#' General purpose function that identifies values in a vector that exceed a
#' maximum number of characters in length. This function is used to identify
#' long sample names and group labels that may lead to sub-optimal visualizations
#'
#' @param x a vector of values. This parameter is typically a vector of
#'   column names, row names, or vector of values in a particular metadata column.
#' @param char_thresh numeric value greater than or equal to 0.
#'   (numeric; >= 1).
#'   which values in x are long.
#' @param na.rm logical, Should NA's be omitted?. Default: FALSE
#'
#'
#' @return The function returns a vector of values in x that exceed the
#'   maximum number of characters threshold. If all values are below the
#'   char_thresh the function returns NULL. If \emph{na.rm = FALSE}
#'   any NA and NaN values in x will be included in the returned output.
#'
#'   ## x(no NAs) + na.rm = TRUE/FALSE + all below ==> NULL
#'   ## x(NAs) + na.rm = TRUE + all below ==> NULL
#'   ## x(NAs) + na.rm = FALSE + all below ==> c(NA, NaN)
#'
#'   ## x(no NAs) + na.rm = TRUE/FALSE + above == > c(values)
#'   ## x(NAs) + na.rm = TRUE + above ==> c(values)
#'   ## x(NAs) + na.rm = FALSE + above ==> c(values, NA, NaN)
#'
#'
#' @examples
#'
#' \dontrun{
#'
#'  x  <- c("Con_01", "Con_02", NA, "TAC_Treatment_01", "", Inf, NaN, "TAC_Treat_02")
#'  y  <- c("Con_01", "Con_02", "Treat_01", "Treat_101")
#'  z  <- c("Con_01", NA, "TAC_Treatment_01", "TAC_Treat_02", "", Inf)
#'
#'  ## identify values with more than 15 characters
#'  check_long(x = x, char_thresh = 15, na.rm = FALSE)
#'
#'  ## identify values with more than 8 characters
#'  check_long(x = y, char_thresh = 8)
#'
#'  ## identify values with more than 6 characters
#'  check_long(x = y, char_thresh = 6)
#'
#' }
#'
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


