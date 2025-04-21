#' Make an HTMLwidget for the interactive report
#'
#' This internal function is called within the .Rmd report template as
#' part of the user-facing \code{\link{write_limma_plots}} function. It takes in
#' data on statistical results, annotation, and sample intensities (partially
#' assembled in \code{\link{write_limma_plots}} and then passed on to the .RMD
#' environment), does some further processing, packages it into the form needed
#' for our interactive report, and then outputs an HTMLwidget for use in the HTML
#' file generated from the .Rmd.
#'
#' @param DAList A DAList object.
#' @param contrast The contrast name.
#' @param display.columns A vector of columns to display in the output table.
#' @param grouping_column The metadata column indicating sample groups.
#' @param status.cols A vector of colors to use for downregulated, nonDE,
#'   and upregulated proteins. Must be of length 3.
#' @param sample.cols A vector of colors for each sample. Should have the same
#'   length as groups.
#' @param height The height of the interactive report objects, in pixels.
#' @param width The width of the interactive report objects, in pixels.
#'
#' @return An HTMLwidget containing our interactive plots and tables
#' @keywords internal

uams_glimmaXY <- function(DAList,
                          contrast,
                          display.columns,
                          grouping_column,
                          status.cols,
                          sample.cols,
                          width,
                          height) {
  
  model_data <- prep_plot_model_data(DAList$results, contrast)
  
  if (!is.null(DAList$data_per_contrast) && contrast %in% names(DAList$data_per_contrast)) {
    counts <- DAList$data_per_contrast[[contrast]]
  } else {
    counts <- DAList$data[rownames(model_data), , drop = FALSE]
  }
  
  if (!is.null(DAList$annotation_per_contrast) && contrast %in% names(DAList$annotation_per_contrast)) {
    anno <- DAList$annotation_per_contrast[[contrast]]
  } else {
    anno <- DAList$annotation[rownames(model_data), , drop = FALSE]
  }
  
  anno <- anno[rownames(model_data), , drop = FALSE]
  anno$p <- round(model_data$P.Value, 4)
  anno$adjusted_p <- round(model_data$adj.P.Val, 4)
  
  if (is.null(model_data$internal_title_column)) {
    model_data$internal_title_column <- rownames(model_data)
  }
  
  groups_in_contrast <- unlist(strsplit(contrast, split = "_vs_"))
  
  sample_groups <- DAList$metadata[[grouping_column]]
  names(sample_groups) <- rownames(DAList$metadata)
  included_samples <- names(sample_groups)[sample_groups %in% groups_in_contrast]
  
  counts <- counts[, colnames(counts) %in% included_samples, drop = FALSE]
  grouping_vector <- sample_groups[colnames(counts)]
  if (is.factor(grouping_vector)) grouping_vector <- as.character(grouping_vector)
  
  if (length(status.cols) != 3) {
    stop("status.cols must have exactly 3 elements for [downreg, notDE, upreg]")
  }
  if (ncol(counts) != length(grouping_vector)) {
    stop("Length of group vector must match number of columns in counts")
  }
  if (nrow(model_data) == 0) stop("model_data is empty")
  if (!"logFC" %in% colnames(model_data)) stop("logFC missing in model_data")
  if (!"P.Value" %in% colnames(model_data)) stop("P.Value missing in model_data")
  if (!"adj.P.Val" %in% colnames(model_data)) stop("adj.P.Val missing in model_data")
  if (!"uniprot_id" %in% colnames(anno)) stop("uniprot_id missing in annotation")
  
  table_for_widget <- cbind(model_data, anno)
  table_for_widget$internal_id_for_brushing <- rownames(model_data)
  display.columns <- c("internal_id_for_brushing", display.columns)
  
  table_for_widget <- data.frame(lapply(table_for_widget, function(col) {
    if (is.numeric(col)) round(col, 4) else col
  }))
  table_for_widget <- data.frame(index = 0:(nrow(table_for_widget) - 1), table_for_widget)
  
  req_cols <- c("p", "adjusted_p", "uniprot_id",
                "sig.pval.fct", "sig.FDR.fct",
                "negLog10rawP", "negLog10adjP",
                "logFC", "average_intensity",
                "internal_id_for_brushing", "internal_title_column", "index")
  output_cols <- unique(c(req_cols, display.columns))
  
  table_for_widget <- table_for_widget[, colnames(table_for_widget) %in% output_cols, drop = FALSE]
  counts <- data.frame(round(counts, 4))
  
  groups_df <- data.frame(group = grouping_vector, sample = colnames(counts))
  
  xData <- list(
    data = list(
      table = table_for_widget,
      cols = display.columns,
      counts = counts,
      groups = groups_df,
      expCols = colnames(groups_df),
      numUniqueGroups = length(unique(groups_df$group)),
      statusColours = status.cols,
      sampleColours = if (is.null(sample.cols)) {-1} else {sample.cols},
      values = as.matrix(counts[rownames(model_data), , drop = FALSE]),
      title = model_data$internal_title_column
    )
  )
  
  htmlwidgets::createWidget(
    name = "glimmaXY",
    xData,
    package = "proteoDAstjude",
    width = width,
    height = height,
    elementId = NULL,
    sizingPolicy = htmlwidgets::sizingPolicy(
      defaultWidth = width,
      defaultHeight = height,
      browser.fill = TRUE,
      viewer.suppress = TRUE
    )
  )
}
