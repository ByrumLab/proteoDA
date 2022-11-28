
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Removes commas in character representation of a number
#'
#' @param x Character string of a number, with commas in it.
#'
#' @return Integer
#'
#' @keywords internal
#'
#'
remove_commas <- function(x) {
  as.numeric(gsub("\\,", "", x))
}


#' %notin% operator, opposite of %in%
#'
#'
#' @inheritParams base::"%in%"
#' @return A logical vector, indicating if a match was NOT located for each
#' element of of the search list.
#'
#' @keywords internal
#'
#'
#'
`%notin%` <- function(x,table) !`%in%`(x,table)



#' Convert a column to a factor
#'
#' Converts a vector/column to a factor. For numeric inputs, adds an "X"
#' in front of the number by default, or can supply your own prefix to use
#' instead.
#'
#' @param x Input vector to be converted to a factor
#' @param prefix Optional: a prefix to put in front of numbers when
#'   converting to a factor. Default is "X"
#'
#' @return A factor vector
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
make_factor <- function(x, prefix = "X") {
  if (is.numeric(x)) {
    x <- paste(prefix, x, sep = "")
  }
  x <- factor(x, levels = ordered(unique(x)))
  return(x)
}



#' All pairwise differences of a vector
#'
#' Taken from \url{https://stackoverflow.com/questions/48445003/compute-all-pairwise-differences-of-elements-in-a-vector},
#' uses a nice R matrix algebra trick, taking the outer product of the vectors
#' (in this case, the difference instead of the product). Then, just return the
#' lower triangle, so we don't double-count pw comparisons
#'
#' @param vector A numeric vector of items for which you want all pairwise
#'   differences
#'
#' @return A vector of all pairwise differences between the elements in the vector
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
all_pw_diffs <- function(vector) {
  all_diffs <- outer(vector, vector, "-")
  all_diffs[lower.tri(all_diffs)]
}


#' Get colors for groups, function 2
#'
#' Used to get colors for the groups in our missing value heatmaps.
#'
#' @param group A vector of group names.
#'
#' @return A vector of colors for each unique group
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
#'
colorGroup <- function(group) {
  blueberry <- "#1F5EDC";  cherry <- "#EE0010"
  apple     <- "#32CD32";  barbie <- "#FF1493"
  fanta     <- "#FF7F00";  grape  <- "#A342FC"
  ocean     <- "#00C8FF";  mtndew <- "#ADFF2F"
  gold      <- "#FFE100";  orchid <- "#E36EF6"
  aceblue   <- "#009ACE";  poop   <- "#996633"

  binfcolors <-c(blueberry, cherry, apple, barbie, fanta,
                 grape, ocean, mtndew, gold, orchid, aceblue, poop)
  names(binfcolors)<-c("blueberry", "cherry", "apple", "barbie", "fanta",
                       "grape", "ocean", "mtndew", "gold", "orchid", "aceblue", "poop")

  ## group=6-12
  if(length(unique(group)) > 5){
    groupCol <- binfcolors[1:length(unique(group))]
    names(groupCol) <- unique(group)
  }
  if(length(unique(group)) > 12) { # Rare
    if (!requireNamespace("grDevices", quietly = TRUE)) {
      cli::cli_abort(c("Package \"grDevices\" must be installed to make plots for more than 12 groups"))
    }
    groupCol <- c(grDevices::rainbow(length(unique(group))))
    names(groupCol) <- unique(group)
  }
  if(length(unique(group)) <= 5) {
    groupCol<-binfcolors[1:length(unique(group))]
    names(groupCol) <- unique(group)
  }
  if(length(unique(group)) == 3) {
    groupCol<-c(binfcolors[c("blueberry","barbie","apple")])
    names(groupCol) <- unique(group)
  }
  if(length(unique(group)) == 2) {
    groupCol<-binfcolors[c("aceblue","cherry")]
    names(groupCol) <- unique(group)
  }
  return(groupCol)
}

#' Validate filename
#'
#' Internal function used to check constraints on user-defined filenames. By
#' default checkes whether there are spaces in the filename, and whether the
#' extension is in the list of allowed extensions.
#'
#' @param filename File name to check
#' @param allowed_exts Character vector of allowable extensions, with no period
#' @param check_space Check for whitespace in file name?
#'
#' @return If filename checks out, returns invisible TRUE. Otherwise, throws
#' error.
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
validate_filename <- function(filename, allowed_exts, check_space = T) {

  # Check spaces
  if (check_space) {
    if (stringr::str_detect(filename, " ")) {
      cli::cli_abort(c("Invalid file name supplied: {.path {filename}}",
                       "x" = "Supplied file name {.val {filename}} contains spaces, which are not allowed",
                       "i" = "Edit the supplied file name to remove spaces"))
    }
  }

  # check extension:
  if (file_extension(filename) %notin% allowed_exts) {
    cli::cli_abort(c("Invalid file name supplied: {.path {filename}}",
                     "x" = "File name cannot end in {.path {file_extension(filename)}}",
                     "i" = "Permitted file {cli::qty(length(allowed_exts))} extension{?s}: {.val {allowed_exts}}"))

  }

  invisible(TRUE)
}

#' Get the file extension from a file
#'
#' @param filepath The filepath from which to extract the extension
#'
#' @return A string of the file extension. When no extension was detected, an empty string.
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet
file_extension <- function(filepath) {
  ext <- stringr::str_extract(filepath, "\\.[0-9A-Za-z]+$") %>%
    stringr::str_remove_all("\\.")
  if (is.na(ext)) ext <- ""
  ext
}

#' Calculate per-row variances of a numeric array
#'
#' Have done testing, gives same results as the genefilter::rowVars() function.
#'
#' @param x The array for which to calculate per-row variances
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row variances.
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.
rowVars <- function(x, ...) {
  # Get variance of each row by vectorized subtraction of rowmeans
  # This is much faster than, e.g., applying the built-in var() function over the rows of a matrix
  rowSums((x - rowMeans(x, ...))^2, ...)/(rowSums(!is.na(x)) - 1)
}


#' Calculate per-row standard deviations of a numeric array
#'
#' Have done testing, gives same results as the genefilter::rowSds() function.
#'
#' @param x The array for which to calculate per-row Sds
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row standard deviations
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.
rowSds <- function(x, ...) {
  sqrt(rowVars(x, ...))
}

#' Calculate per-row medians of a numeric array
#'
#' Have done testing, gives same results as the matrixStats::rowMedians() function
#' it replaces, though it is much slower (the matrixStats version uses C code).
#'
#' @param x The array for which to calculate per-row medians
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row medians
#'
#' @keywords internal
#'
#' @examples
#' # No examples yet.
rowMedians <- function(x, ...) { # Much slower than the matrixStats version, which uses C
  # but probably worth losing the dependency.
  apply(x, MARGIN = 1, FUN = stats::median, ...)
}

#' Calculate per-row MADs of a numeric array
#'
#' Have done testing, gives same results as the matrixStats::rowMads() function
#' it replaces, though it is much slower (the matrixStats version uses C code).
#'
#' @param x The array for which to calculate per-row MADs
#' @param ... Additional arguments to be passed to internal functions.
#'   Meant for na.rm.
#'
#' @keywords internal
#'
#' @return A numeric vector of appropriate length, named if input was named, with
#'   per-row MADs
#'
#' @examples
#' # No examples yet.
rowMads <- function(x, ...) {
  apply(x, MARGIN = 1, FUN = stats::mad, ...)
}


#' Internal function that checks that the rownames of set of dataframes have matching rownames in a reference set
#'
#' @param obj A list of dataframes to evaluate
#' @param ref_rows A vector of row names to which the row names of each obj dataframe will be compared
#'
#' @keywords internal
#'
#' @return A logical, with TRUE indicating all dataframes have matching rows in the reference set, and
#' FALSE indicates that the one or more dataframes do not have matching rownames with the reference
check_rows_in <- function(obj=list(), ref_rows=c()){
  all(unlist(lapply(obj, function(x){
    all(rownames(x) %in% ref_rows)
  })))
}

