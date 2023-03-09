uams_glimmaXY <- function(model_data,
                          counts,
                          groups,
                          status,
                          anno,
                          display.columns,
                          status.cols,
                          sample.cols,
                          main = "",
                          width = 920,
                          height = 920) {

  # The model table should be the output from
  # prep_plot_model_data.
  table <- model_data
  table$internal_id_for_brushing <- rownames(model_data)
  display.columns <- c("internal_id_for_brushing", display.columns)

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


  # Round all numeric values to 4 digits, to make everything look nicer in tables
  table <- data.frame(lapply(table, FUN = function(col) {
    if (is.numeric(col)) {
      output <- round(col, digits = 4)
    } else {
      output <- col
    }
    output
  }))

  anno <- data.frame(lapply(anno, FUN = function(col) {
    if (is.numeric(col)) {
      output <- round(col, digits = 4)
    } else {
      output <- col
    }
    output
  }))
  # Counts is a matrix, can just round it
  counts <- round(counts, digits = 4)


  # Build the data to pass to the htmlwidget
  xData <- UAMS_buildXYData(table, status, main, display.columns,
                            anno, counts, status.cols, sample.cols, groups)

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


UAMS_buildXYData <- function(table, status, main, display.columns, anno, counts,
                             status.cols, sample.cols, groups) {

  counts <- data.frame(counts)
  level <- levels(groups)
  groups <- data.frame(group = groups)
  groups <- cbind(groups, sample = colnames(counts))

  if (!is.null(anno)) {
    table <- cbind(table, anno)
  }

  table <- data.frame(index = 0:(nrow(table) - 1), table)

  xData <- list(
    data = list(
      table = table,
      cols = display.columns,
      counts = counts,
      groups = groups,
      levels = level,
      expCols = colnames(groups),
      numUniqueGroups = length(unique(groups$group)),
      annoCols = if (is.null(anno)) {-1} else {colnames(anno)},
      statusColours = status.cols,
      sampleColours = if (is.null(sample.cols)) {-1} else {sample.cols},
      samples = colnames(counts),
      title = main
    )
  )
  xData
}
