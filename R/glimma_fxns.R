#' Make an HTMLwidget for the interactive report
#'
#' This internal function is called within the .Rmd report template as
#' part of the user-facing \code{\link{write_limma_plots}} function. It takes in
#' data on statistical results, annotation, and sample intensities (partially
#' assembled in \code{\link{write_limma_plots}} and then passed on the to .RMD
#' environment), does some further processing, packages
#' it into the form needed for our interactive report, and then outputs an
#' HTMLwidget for use in the HTML file generated from the .Rmd.
#'
#' @param DAList A DAList object.
#' @param model_data The output from the \code{\link{prep_plot_model_data}}
#'   function, containing statistical results for a single contrast.
#' @param contrast The contrast name.
#' @param display.columns A vector of columns to display in the output table.
#' @param grouping_column The metadata column indicating sample groups.
#' @param status.cols A vector of colors to use for down regulated, nonDE,
#'   and upregulated proteins. Must be of length 3.
#' @param sample.cols A vector of colors for each sample. Should have the same
#'   length as groups.
#' @param height The height of the interactive report objects, in pixels.
#' @param width The width of the interactive report objects, in pixels.
#'
#' @return An HTMLwidget containing our interactive plots and tables
#'
#' @keywords internal
#'
uams_glimmaXY <- function(DAList,
                          model_data,
                          contrast,
                          display.columns,
                          grouping_column,
                          status.cols,
                          sample.cols,
                          width,
                          height) {
  
  # Pull counts and annotation conditionally
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
  
  # Align annotation to model_data
  anno <- anno[rownames(model_data), , drop = FALSE]
  
  # Add p-values
  anno$p <- round(model_data$P.Value, 4)
  anno$adjusted_p <- round(model_data$adj.P.Val, 4)
  
  # Add internal title column if not already present
  if (is.null(model_data$internal_title_column)) {
    model_data$internal_title_column <- rownames(model_data)
  }
  
  # Ensure group vector is aligned with counts
  grouping_vector <- DAList$metadata[colnames(counts), grouping_column, drop = TRUE]
  if (is.factor(grouping_vector)) grouping_vector <- as.character(grouping_vector)
  
  if (length(status.cols) != 3) {
    stop("status.cols must have exactly 3 elements for [downreg, notDE, upreg]")
  }
  if (ncol(counts) != length(grouping_vector)) {
    stop("Length of group vector must match number of columns in counts")
  }
  
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
  
  groups_df <- data.frame(group = grouping_vector,
                          sample = colnames(counts))
  
  xData <- list(
    data = list(
      table = table_for_widget,
      cols = display.columns,
      counts = counts,
      groups = groups_df,
      expCols = colnames(groups_df),
      numUniqueGroups = length(unique(groups_df$group)),
      statusColours = status.cols,
      sampleColours = if (is.null(sample.cols)) {-1} else {sample.cols}
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