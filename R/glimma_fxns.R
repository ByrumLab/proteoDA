uams_glimmaXY <- function(model_data,
                          counts,
                          groups,
                          status,
                          anno,
                          display.columns,
                          status.cols,
                          sample.cols,
                          main = "",
                          html = NULL,
                          width = 920,
                          height = 920) {

  # The model table should be the output from
  # prep_plot_model_data.
  # Originally these were rounded to 4 digits, but maybe don't need to do that?
  # could if we really need to save space
  table <- model_data
  table$gene <- rownames(model_data)

  xData <- UAMS_buildXYData(table, status, main, display.columns,
                            anno, counts, status.cols, sample.cols, groups)
  return(UAMS_glimmaXYWidget(xData, width, height, html))
}


UAMS_glimmaXYWidget <- function (xData, width, height, html) {
  widget <- htmlwidgets::createWidget(name = "glimmaXY", xData, package = "proteoDA",
                                      width = width, height = height, elementId = NULL,
                                      sizingPolicy = htmlwidgets::sizingPolicy(defaultWidth = width,
                                                                               defaultHeight = height, browser.fill = TRUE, viewer.suppress = TRUE))
  if (is.null(html)) {
    return(widget)
  }
  else {
    message("Saving widget...")
    htmlwidgets::saveWidget(widget, file = html)
    message(html, " generated.")
  }
}


UAMS_buildXYData <- function(table, status, main, display.columns, anno, counts,
                             status.cols, sample.cols, groups) {

  counts <- data.frame(counts)
  if (is.null(groups)) {
    groups <- factor("group")
  }
  else {
    if (ncol(counts) != length(groups))
      stop("Length of groups must be equal to the number of columns in counts.\n")
  }

  level <- levels(groups)
  groups <- data.frame(group = groups)
  groups <- cbind(groups, sample = colnames(counts))

  status <- sapply(status, function(x) {
    switch(as.character(x), `-1` = "downReg", `0` = "nonDE",
           `1` = "upReg")
  })
  if (length(status) != nrow(table)) {
    stop("Status vector\n     must have the same number of genes as the main arguments.")
  }


  table <- cbind(table, status = as.vector(status))
  if (!is.null(anno)) {
    table <- cbind(table, anno)
  }
  if (!("gene" %in% display.columns)) {
    display.columns <- c("gene", display.columns)
  }
  table <- data.frame(index = 0:(nrow(table) - 1), table)
  if (length(status.cols) != 3) {
    stop("status.cols\n arg must have exactly 3 elements for [downreg, notDE, upreg]")
  }
  xData <- list(data = list(table = table,
                            cols = display.columns, counts = counts, groups = groups,
                            levels = level, expCols = colnames(groups), annoCols = if (is.null(anno)) {
                              -1
                            } else {
                              colnames(anno)
                            }, statusColours = status.cols, sampleColours = if (is.null(sample.cols)) {
                              -1
                            } else {
                              sample.cols
                            }, samples = colnames(counts), title = main))
  return(xData)
}
