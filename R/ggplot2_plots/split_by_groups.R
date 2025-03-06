

#' Split a data.frame or matrix according to group membership.
#'
#' Split a data.frame or matrix according group membership, creating a list
#' where each element consists of a data.frame/matrix for a particular
#' experimental group.
#'
#' @param data a data.frame or matrix. For our pipeline this is typically a counts
#'   matrix or matrix of counts-per-million values.
#' @param groups a vector or factor giving the experimental group/condition for
#'   each sample (i.e. column) in data. The length of groups must equal the
#'   number of columns in data.
#'
#' @rdname split_by_groups
#'
#' @keywords internal
#'
#' @return The function returns a list in which each element is a data.frame
#'   for an individual group.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#'  counts <- matrix(rpois(4*5, lambda = 5), nrow = 4,ncol = 5)
#'  rownames(counts) <- paste0("gene_", 1:4)
#'  colnames(counts) <- paste0("Sample_", 1:5)
#'  y <- read_rnaseq_data(counts = counts, project_id = "LastName_Date")
#'  y <- validate_DGEList(y = y)
#'
#'  tmp <- split_by_groups(data   = y$counts,
#'                         groups = c("control","treat","treat","control","treat"))
#'
#'  }
#' }
#'
split_by_groups <- function(data,
                            groups){
 
 
 ## check data argument
 if (!check_data(x = data)) {
  cli::cli_abort("{.arg data} must be a data.frame or matrix.")
 }
 
 
 ## check groups argument
 if (length(groups) != ncol(data)) {
  cli::cli_abort("{.arg groups} is not the same length as the number of
                     columns in data")
 }
 
 
 if (any(is.na(groups),
         (nchar(as.character(groups)) == 0L))) {
  cli::cli_abort("{.arg groups} cannot contain missing values or
                   empty strings.")
 }
 
 
 group <- unique(groups)
 group_indx <- lapply(group, function(x){ which(groups %in% x) })
 split_data <- lapply(group_indx, function(x){
  if(length(x) == 1L){
   df <- data.frame(data[, x])
   colnames(df) <- colnames(data)[x]
  } else {
   df <- data.frame(data[, x])
  }
  df
 })
 names(split_data) <- group
 
 
 split_data
 
}

