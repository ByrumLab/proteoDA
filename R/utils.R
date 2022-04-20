
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
remove_commas <- function(x){ x<-as.numeric(gsub("\\,", "", x)) }


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
#' @export
#'
#' @examples
#' # No examples yet
#'

# TODO: check out forcats for this? how often do we use the prefix thing?
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
#' @param vector
#'
#' @return A vector of all pairwise differences between the elements in the vector
#' @export
#'
#' @examples
#' # No examples yet
#'
all_pw_diffs <- function(vector) {
  all_diffs <- outer(vector, vector, "-")
  all_diffs[lower.tri(all_diffs)]
}



#' Make new filename to avoid overwriting
#'
#' Takes in a filename and a directory for a file that already exists. Tries to
#' make a new filename by adding "_01" between the filename and extension, and
#' continues to do so, increasing the number, until it finds a filename that does
#' not yet exist. Returns that filename.
#'
#' @param x Non-unique filename we're trying to replace.
#' @param dir directory of the file we're trying to create
#'
#' @return A new, unique filename
#'
#' @export
#'
#' @examples
#' # No examples yet
#'
make_new_filename <- function(x, dir) {
  # Parse input filename and get current file list
  ext <- tools::file_ext(x)
  base_name <- stringr::str_remove(x, paste0(".", ext))
  current_files <- list.files(path=dir)

  # Loop setup
  attempt <- 0
  success <- FALSE
  # Start looping
  while (!success) {
    if (attempt >= 50) {
      cli::cli_abort("Could not create new unique filename to replace {.path {x}} after {.val {attempt}} tries. ")
    }
    attempt <- attempt + 1

    suffix <- stringr::str_pad(as.character(attempt), width = 2, side = "left", pad = "0")
    new_name <- paste0(base_name, "_", suffix, ".", ext)
    success <- new_name %notin% current_files
  }

  new_name
}

