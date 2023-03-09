uams_glimmaXY <- function(model_data,
                          counts,
                          groups,
                          status,
                          anno,
                          display.columns,
                          status.cols,
                          sample.cols,
                          width,
                          height) {

  # Some possible mild argument checking?
  # Though user will never really call this
  if (length(status.cols) != 3) {
    stop("status.cols\n arg must have exactly 3 elements for [downreg, notDE, upreg]")
  }
  if (is.null(groups)) {
    groups <- factor("group")
  }
  else {
    if (ncol(counts) != length(groups))
      stop("Length of groups must be equal to the number of columns in counts.\n")
  }

  # The model table should be the output from
  # prep_plot_model_data.
  table_for_widget <- cbind(model_data, anno)
  table_for_widget$internal_id_for_brushing <- rownames(model_data)
  display.columns <- c("internal_id_for_brushing", display.columns)


  # Round all numeric values to 4 digits, to make everything look nicer in tables
  table_for_widget <- data.frame(lapply(table_for_widget, FUN = function(col) {
    if (is.numeric(col)) {
      output <- round(col, digits = 4)
    } else {
      output <- col
    }
    output
  }))

  table_for_widget <- data.frame(index = 0:(nrow(table_for_widget) - 1), table_for_widget)

  # To reduce final size, only output cols we need
  # Required cols for plotting
  req_cols <- c("p", "adjusted_p", "uniprot_id",
                "sig.pval.fct", "sig.FDR.fct",
                "negLog10rawP", "negLog10adjP",
                "logFC", "average_intensity",
                "internal_id_for_brushing", "internal_title_column", "index")
  # And the user-requested display columns
  output_cols <- unique(c(req_cols, display.columns))

  table_for_widget <- table_for_widget[,colnames(table_for_widget) %in% output_cols, drop = F]


  # Counts is a matrix, can just round it
  counts <- data.frame(round(counts, digits = 4))

  groups_df <- data.frame(group = groups,
                          sample = colnames(counts))

  # Build the data to pass to the htmlwidget
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

  # make the widget
  widget <- htmlwidgets::createWidget(
    name = "glimmaXY",
    xData,
    package = "proteoDA",
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

  widget
}
