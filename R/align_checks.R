#' Align and validate sample correspondence between data and metadata
#'
#' @description
#' In many quantitative proteomics workflows, mass spectrometry sample files are
#' named using numeric suffixes (e.g. `Sample1`, `Sample2`, ..., `Sample10`).
#' When read into R or other languages, these can be automatically sorted
#' lexicographically (`Sample1`, `Sample10`, `Sample11`, `Sample2`, ...),
#' causing incorrect sample-to-group mapping and misaligned metadata.
#'
#' This helper ensures that the \code{data} and \code{metadata} objects are
#' correctly matched and ordered before creating a \code{DAList}. It automatically
#' fixes ordering issues, groups replicates of the same condition together, and
#' performs natural ("human") sorting of sample names containing numbers.
#'
#' @param data A numeric matrix or data frame of protein intensities,
#'   where rows are proteins and columns are samples.
#' @param metadata A data frame of sample metadata, one row per sample.
#' @param sample_col Optional name of the column in \code{metadata}
#'   that contains the sample identifiers. If NULL, row names are used.
#' @param group_col Optional column name in \code{metadata} used to group
#'   replicates together (e.g. "group", "condition").
#' @param prefer_group_blocks Logical; if TRUE (default) and \code{group_col}
#'   is supplied, both \code{data} and \code{metadata} will be reordered to
#'   keep replicates of the same group contiguous and sorted naturally.
#'   If FALSE, only the metadata will be reordered to match data columns.
#' @param strict Logical; if TRUE, stops with an error if sample orders differ.
#'   If FALSE (default), automatically reorders when sets match.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{data}{The reordered numeric intensity data.}
#'   \item{metadata}{The reordered sample metadata.}
#'   \item{changes}{A data frame summarizing the old vs new order of samples.}
#' }
#'
#' @details
#' The function performs three major checks:
#' \enumerate{
#'   \item Ensures the same sample IDs appear in both \code{data} and \code{metadata}.
#'   \item Detects and removes any ordering mismatches.
#'   \item (Optionally) groups replicates of the same experimental condition together
#'         using natural sorting.
#' }
#'
#' If sample names differ only by order, the function reorders them safely.
#' If samples are missing or duplicated, it stops with a descriptive error.
#'
#' @examples
#' \dontrun{
#' # Example data with misordered samples
#' dat <- data.frame(
#'   Sample1  = c(10, 20, 30),
#'   Sample10 = c(11, 21, 31),
#'   Sample2  = c(12, 22, 32)
#' )
#' rownames(dat) <- paste0("Prot", 1:3)
#'
#' meta <- data.frame(
#'   sample = c("Sample1", "Sample2", "Sample10"),
#'   group  = c("A", "A", "B")
#' )
#' rownames(meta) <- meta$sample
#'
#' # Align automatically, grouping by condition
#' aligned <- align_data_and_metadata(
#'   data = dat,
#'   metadata = meta,
#'   sample_col = "sample",
#'   group_col = "group",
#'   prefer_group_blocks = TRUE
#' )
#'
#' # Check results
#' aligned$changes
#' head(aligned$data)
#' head(aligned$metadata)
#' }
#'
#' @export
align_data_and_metadata <- function(data,
                                    metadata,
                                    sample_col = NULL,
                                    group_col = NULL,
                                    prefer_group_blocks = TRUE,
                                    strict = FALSE) {
  .natural_order <- function(x) {
    num <- suppressWarnings(as.integer(sub(".*?(\\d+)$", "\\1", x)))
    prefix <- sub("(.*?)(\\d+)$", "\\1", x)
    num[is.na(num)] <- Inf
    prefix[is.na(prefix)] <- x[is.na(prefix)]
    order(prefix, num, x, method = "radix")
  }
  
  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("`data` must be numeric in all columns (one column per sample).")
  }
  
  md <- metadata
  if (!is.null(sample_col)) {
    if (!sample_col %in% colnames(md)) {
      stop(sprintf("`sample_col` '%s' not found in metadata.", sample_col))
    }
    sample_ids <- as.character(md[[sample_col]])
  } else {
    if (is.null(rownames(md)) || any(rownames(md) == "")) {
      stop("Metadata has no rownames; set rownames or provide `sample_col`.")
    }
    sample_ids <- rownames(md)
  }
  
  if (any(duplicated(colnames(data))))
    stop("Duplicate sample IDs in `data` columns: ", paste(unique(colnames(data)[duplicated(colnames(data))]), collapse = ", "))
  if (any(duplicated(sample_ids)))
    stop("Duplicate sample IDs in `metadata`: ", paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", "))
  
  dat_ids <- colnames(data)
  miss_in_md  <- setdiff(dat_ids, sample_ids)
  miss_in_dat <- setdiff(sample_ids, dat_ids)
  
  if (length(miss_in_md) > 0 || length(miss_in_dat) > 0) {
    msg <- c()
    if (length(miss_in_md)  > 0) msg <- c(msg, paste0("Missing in metadata: ", paste(miss_in_md,  collapse = ", ")))
    if (length(miss_in_dat) > 0) msg <- c(msg, paste0("Missing in data: ",     paste(miss_in_dat, collapse = ", ")))
    stop(paste(c("Sample ID mismatch between data and metadata.", msg), collapse = "\n - "))
  }
  
  new_data <- data
  new_md   <- md
  
  if (prefer_group_blocks && !is.null(group_col) && group_col %in% colnames(md)) {
    group <- md[[group_col]]
    group_key <- as.character(group)
    ord_groups <- order(group_key, method = "radix")
    grouped_ids <- split(seq_along(sample_ids)[ord_groups], group_key[ord_groups])
    final_idx <- unlist(lapply(grouped_ids, function(ix) ix[.natural_order(sample_ids[ix])]), use.names = FALSE)
    target_ids <- sample_ids[final_idx]
    if (!identical(sample_ids, target_ids)) {
      new_md <- md[match(target_ids, sample_ids), , drop = FALSE]
      rownames(new_md) <- if (is.null(sample_col)) target_ids else rownames(md)[match(target_ids, sample_ids)]
    }
    if (!identical(colnames(data), target_ids)) {
      new_data <- data[, target_ids, drop = FALSE]
    }
  } else {
    if (!identical(dat_ids, sample_ids)) {
      new_md <- md[match(dat_ids, sample_ids), , drop = FALSE]
      rownames(new_md) <- if (is.null(sample_col)) dat_ids else rownames(md)[match(dat_ids, sample_ids)]
    }
  }
  
  changes <- data.frame(
    position = seq_along(colnames(new_data)),
    sample   = colnames(new_data),
    old_pos_in_metadata = match(colnames(new_data), sample_ids),
    old_pos_in_data     = match(colnames(new_data), dat_ids),
    stringsAsFactors = FALSE
  )
  
  if (strict && (!identical(dat_ids, colnames(new_data)) || !identical(sample_ids, rownames(new_md)))) {
    stop("Order mismatch detected. Set strict=FALSE to auto-align, or fix input files and rerun.")
  }
  
  list(data = new_data, metadata = new_md, changes = changes)
}


#' Check data–metadata alignment consistency
#'
#' @description
#' A lightweight check to verify that the columns in a data matrix and
#' the rows in the corresponding metadata frame refer to the same samples.
#' Typically used after importing data but before building a \code{DAList}.
#'
#' @param data A numeric matrix or data frame of intensity values.
#' @param metadata A data frame of sample metadata with row names matching
#'   the sample IDs in \code{data}.
#'
#' @return Invisibly returns TRUE if the match is valid.
#' Stops with an informative error if there are missing or mismatched samples.
#'
#' @examples
#' \dontrun{
#' dat <- data.frame(S1 = 1:3, S2 = 4:6)
#' meta <- data.frame(group = c("A", "B"), row.names = c("S1", "S2"))
#' check_data_metadata_match(dat, meta)  # passes
#'
#' bad_meta <- data.frame(group = c("A", "B"), row.names = c("S1", "S3"))
#' check_data_metadata_match(dat, bad_meta)  # stops with descriptive error
#' }
#'
#' @export
check_data_metadata_match <- function(data, metadata) {
  md_ids <- if (!is.null(rownames(metadata))) rownames(metadata) else {
    stop("Metadata must have rownames, or use align_data_and_metadata() with `sample_col`.")
  }
  if (!setequal(colnames(data), md_ids)) {
    miss_in_md  <- setdiff(colnames(data), md_ids)
    miss_in_dat <- setdiff(md_ids, colnames(data))
    msg <- c()
    if (length(miss_in_md))  msg <- c(msg, paste0("Missing in metadata: ", paste(miss_in_md,  collapse = ", ")))
    if (length(miss_in_dat)) msg <- c(msg, paste0("Missing in data: ",     paste(miss_in_dat, collapse = ", ")))
    stop(paste(c("Sample ID mismatch between data and metadata.", msg), collapse = "\n - "))
  }
  invisible(TRUE)
}
