

#' Rename samples and/or reorganize DGEList by group levels
#'
#' internal function. NEED TO UPDATE DECRIPTION
#'
#' @param DAList a DAList object with metadata
#' @param label_column Optional. The name of column or column number within the
#'   metadata data frame which contains labels to use for each sample.
#'   When not supplied, defaults to using the column names of the data.
#'   Default: NULL
#' @param grouping_column Optional. The name of the column or column number in
#'   the metadata which gives information on how to group the samples.
#'    Default: NULL
#' @param sort logical. Should the samples in the DAList object be sorted
#'   according to the group levels of the metadata grouping_column? If
#'   grouping_column is NULL sort is ignored (ie set as FALSE)
#'   Default: FALSE
#'
#' @export
#'
#' @return returns a DAList object with updated sample names if label_column
#'   is defined and updated sample order based on group levels of input
#'   grouping_column if supplied and sort is TRUE.
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'
#'  data <- matrix(rnbinom(5 * 8, mu = 5, size = 10), 5, 8)
#'  rownames(data) <- paste0("protein_", 1:5);
#'  colnames(data) <- paste0("S000_", 1:8)
#'
#'  metadata <- data.frame(sample = paste0("S",1:8),
#'                         sample_id = paste0("S000_",1:8),
#'                         condition = c("S","S","P","P","A","A","N","N"),
#'                         batch = c(1, 2, 1, 2, 1, 2, 1, 2),
#'                         row.names = paste0("S000_",1:8))
#'
#' annotation <- data.frame(uniprot_id = paste0("protein_", 1:5),
#'                          check.names = FALSE, check.rows = FALSE,
#'                          fix.empty.names = FALSE, stringsAsFactors = FALSE,
#'                          row.names = paste0("protein_", 1:5))
#'
#' d <- proteoDA:::DAList(data = data, metadata = metadata, annotation = annotation)
#' d <- proteoDA:::validate_DAList(x = d)

#'  ## rename samples
#'  d2 <- rename_samples(DAList = d, label_column = "sample",
#'                       grouping_column = NULL, sort = FALSE)
#'
#'
#'  ## group levels = unique values i.e. make_factor()
#'  d2$metadata$condition <- as.character(d2$metadata$condition)
#'  unique(d2$metadata$condition)
#'  d3 <- rename_samples(DAList = d2, label_column = NULL, grouping_column = "condition",
#'                       sort = TRUE)
#'
#'  ## group levels = group levels i.e. group order defined by user
#'  d3$metadata$condition <- factor(d3$metadata$condition, levels = c("S","N","A","P"))
#'  levels(d3$metadata$condition)
#'  d4 <- rename_samples(DAList = d3, label_column = NULL, grouping_column = "condition",
#'                       sort = TRUE)
#'
#'
#'  }
#' }
#'
rename_samples <- function(DAList,
                           label_column = NULL,
                           grouping_column = NULL,
                           sort = FALSE){

  ## validate DAList object i.e. data, metadata, annotation
  ## have same row/column names and in same order etc.
  DAList <- proteoDA:::validate_DAList(x = DAList)


  ## CHANGE SAMPLE NAMES
  if(!is.null(label_column)){

    ## check label column is a column name and that values are unique
    ## no NAs, blanks, etc.
    label_column <- check_label_column(metadata     = DAList$metadata,
                                       label_column = label_column)

    ## set column names of counts and rownames of samples slot
    ## as label_column values
    rownames(DAList$metadata) <- DAList$metadata[, label_column]
    colnames(DAList$data)  <- rownames(DAList$metadata)

    cli::cli_inform("{.val {label_column}} labels used to rename samples.")
    cli::cli_rule()

  }


  ## SORT DGELIST BY LEVELS OF GROUPING COLUMN
  ## reorder samples according to group levels of grouping_column in metadata slot
  if(!is.null(grouping_column)){

    ## check grouping_column is column name and does not contain NA, blanks,
    ## Inf etc values
    grouping_column <- check_grouping_column(metadata = DAList$metadata,
                                             grouping_column = grouping_column)


    ## if not already make grouping column a factor
    ## group levels will be used to sort samples
    if(is.factor(DAList$metadata[, grouping_column])){
      DAList$metadata[, grouping_column] <- droplevels(x = DAList$metadata[, grouping_column])
    } else {
      DAList$metadata[, grouping_column] <- proteoDA:::make_factor(x = DAList$metadata[, grouping_column])
    }


    ## sort DGEList according to group levels
    if(sort){

      ## sort samples and counts according to group levels
      o <- order(ordered(DAList$metadata[, grouping_column], levels = levels(DAList$metadata[, grouping_column])))
      DAList$metadata <- DAList$metadata[o, ]
      DAList$data  <- DAList$data[, rownames(DAList$metadata)]

      cli::cli_inform(c("sort = {.val {sort}}.",
                        "The DAList object was sorted according to
                        {.val {grouping_column}}
                        group levels:","{.val {levels(DAList$metadata[, grouping_column])}}"))
      cli::cli_rule()

    } else {
      cli::cli_inform(c("sort = {.val {sort}}."," Group level sorting was not
                        applied to the DAList object"))
      cli::cli_rule()
    }

  }

 DAList <- proteoDA:::validate_DAList(x = DAList)

  DAList


}
