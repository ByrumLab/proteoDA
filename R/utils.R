
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
#' @param vector A numeric vector of items for which you want all pairwise
#'   differences
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



#' Get colors for batches
#'
#' Used to get colors for the batches in our missing value heatmaps.
#'
#' @param batch A vector of batch names.
#'
#' @return A vector of colors for each unique batch
#' @export
#'
#' @examples
#' # No examples yet
#'
colorBatch <- function(batch){

  batchCol <- unlist(ifelse(length(unique(batch)) == 1, grDevices::rainbow(1, start=0.5),
                            list(grDevices::rainbow(length(unique(batch))))))
  names(batchCol) <- unique(batch)

  return(batchCol)

}


#' Get colors for groups, function 2
#'
#' Used to get colors for the groups in our missing value heatmaps.
#'
#' @param group A vector of group names.
#'
#' @return A vector of colors for each unique group
#' @export
#'
#' @examples
#' # No examples yet
#'
colorGroup2 <- function(group) {
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
    if(length(unique(group)) > 12) {
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


##------------------------------
##  PLOT RIGHT MARGIN
##------------------------------
right_margin <- function(x) {
  maxchar <- max(nchar(as.character(x))) + 10
  right <- maxchar/5
  return(right)
}


##------------------------------
##  PLOT LEFT MARGIN
##------------------------------
left_margin <- function(x) {
  maxchar <- max(nchar(as.character(x))) + 30
  left <- maxchar/5
  return(left)
}


##------------------------------
##  PLOT HEIGHT
##------------------------------
plot_height <- function(x){
  height <- (600/20)*length(x)
  height <- ifelse(height > 1000, 1000, height)
  return(height)
}

