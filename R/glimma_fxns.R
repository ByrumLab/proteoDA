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
  
  # Ensure annotation is aligned with model_data
  anno <- anno[rownames(model_data), , drop = FALSE]
  
  # Add p-values to annotation
  anno$p <- round(model_data$P.Value, 4)
  anno$adjusted_p <- round(model_data$adj.P.Val, 4)
  
  # Add optional stats
  optional_stats <- c("movingSDs", "logFC_z_scores", "sig.PVal", "sig.FDR")
  for (stat_col in optional_stats) {
    if (!stat_col %in% colnames(model_data) && !is.null(DAList$tags[[stat_col]][[contrast]])) {
      model_data[[stat_col]] <- NA_real_
      stat_vec <- DAList$tags[[stat_col]][[contrast]]
      model_data[names(stat_vec), stat_col] <- stat_vec
    }
  }
  anno$movingSDs <- model_data$movingSDs
  anno$logFC_z_scores <- model_data$logFC_z_scores
  anno$sig.PVal <- model_data$sig.PVal
  anno$sig.FDR <- model_data$sig.FDR
  
  # Assign fallback title if none is provided
  if (is.null(model_data$internal_title_column) ||
      any(is.na(model_data$internal_title_column)) ||
      length(model_data$internal_title_column) != nrow(model_data)) {
    model_data$internal_title_column <- rownames(model_data)
  }
  model_data$internal_title_column <- unname(as.character(model_data$internal_title_column))
  
  # Add name column for JS compatibility
  model_data$name <- model_data$internal_title_column
  
  # Get groups for contrast
  groups_in_contrast <- unlist(strsplit(contrast, split = "_vs_"))
  sample_groups <- DAList$metadata[[grouping_column]]
  names(sample_groups) <- rownames(DAList$metadata)
  included_samples <- names(sample_groups)[sample_groups %in% groups_in_contrast]
  
  counts <- counts[, colnames(counts) %in% included_samples, drop = FALSE]
  grouping_vector <- sample_groups[colnames(counts)]
  if (is.factor(grouping_vector)) grouping_vector <- as.character(grouping_vector)
  
  if (length(status.cols) != 3) stop("status.cols must have exactly 3 elements [down, nonDE, up]")
  if (ncol(counts) != length(grouping_vector)) stop("Group vector length must match counts columns")
  if (nrow(model_data) == 0) stop("model_data is empty")
  if (!"logFC" %in% colnames(model_data)) stop("Missing 'logFC' in model_data")
  if (!"P.Value" %in% colnames(model_data)) stop("Missing 'P.Value' in model_data")
  if (!"adj.P.Val" %in% colnames(model_data)) stop("Missing 'adj.P.Val' in model_data")
  if (!"uniprot_id" %in% colnames(anno)) stop("Missing 'uniprot_id' in annotation")
  
  # Create table for widget
  table_for_widget <- cbind(model_data, anno)
  table_for_widget$internal_id_for_brushing <- rownames(model_data)
  display.columns <- unique(c("internal_id_for_brushing", "name", display.columns))
  
  table_for_widget <- data.frame(lapply(table_for_widget, function(col) {
    if (is.numeric(col)) round(col, 4) else col
  }))
  table_for_widget <- data.frame(index = 0:(nrow(table_for_widget) - 1), table_for_widget)
  table_for_widget$name <- model_data$internal_title_column  # ensure name is available to JS
  
  req_cols <- c("p", "adjusted_p", "uniprot_id",
                "sig.pval.fct", "sig.FDR.fct",
                "negLog10rawP", "negLog10adjP",
                "logFC", "average_intensity",
                "movingSDs", "logFC_z_scores", "sig.PVal", "sig.FDR",
                "internal_id_for_brushing", "internal_title_column", "index", "name")
  output_cols <- unique(c(req_cols, display.columns))
  
  # Filter to expected columns
  table_for_widget <- table_for_widget[, colnames(table_for_widget) %in% output_cols, drop = FALSE]
  
  # Build group/sample structure first
  counts <- data.frame(round(counts, 4))
  groups_df <- data.frame(group = grouping_vector, sample = colnames(counts))
  
  # Extract values matrix
  values <- counts[rownames(model_data), , drop = FALSE]
  rownames(values) <- model_data$internal_title_column
  rownames(counts) <- model_data$internal_title_column
  values <- as.matrix(values)
  
  # Print debug info
  cat("\nDEBUG INFO:\n")
  print(str(list(
    counts_dim = dim(counts),
    values_dim = dim(values),
    table_rows = nrow(table_for_widget),
    table_cols = colnames(table_for_widget),
    display_columns = display.columns,
    num_samples = length(grouping_vector),
    title_preview = head(model_data$internal_title_column)
  )))
  
  # Create interactive widget
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
      values = values,
      title = setNames(as.character(model_data$internal_title_column), table_for_widget$index)
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
