uams_glimmaXY <- function (x, y, xlab = "x", ylab = "y", dge = NULL, counts = dge$counts,
          groups = dge$samples$group, status = rep(0, length(x)), anno = NULL,
          display.columns = NULL, status.cols = c("#1052bd", "silver",
                                                  "#cc212f"), sample.cols = NULL, main = "XY Plot", html = NULL,
          width = 920, height = 920)
{
  if (length(x) != length(y))
    stop("Error: x and y args must have the same length.")
  table <- data.frame(signif(x, digits = 4), signif(y, digits = 4))
  colnames(table) <- c(xlab, ylab)
  if (!is.null(counts)) {
    table <- cbind(gene = rownames(counts), table)
  }
  else if (!is.null(rownames(x))) {
    table <- cbind(gene = rownames(x), table)
  }
  else if (!is.null(rownames(y))) {
    table <- cbind(gene = rownames(y), table)
  }
  else {
    table <- cbind(gene = seq_along(x), table)
  }
  xData <- UAMS_buildXYData(table, status, main, display.columns,
                       anno, counts, xlab, ylab, status.cols, sample.cols, groups)
  return(UAMS_glimmaXYWidget(xData, width, height, html))
}


UAMS_glimmaXYWidget <- function (xData, width, height, html)
{
  widget <- htmlwidgets::createWidget(name = "glimmaXY", xData, package = "proteomicsDIA",
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


UAMS_buildXYData <- function (table, status, main, display.columns, anno, counts,
          xlab, ylab, status.cols, sample.cols, groups)
{
  if (is.null(counts)) {
    counts <- -1
    level <- NULL
  }
  else {
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
  }
  status <- sapply(status, function(x) {
    switch(as.character(x), `-1` = "downReg", `0` = "nonDE",
           `1` = "upReg")
  })
  if (length(status) != nrow(table))
    stop("Status vector\n     must have the same number of genes as the main arguments.")
  table <- cbind(table, status = as.vector(status))
  if (!is.null(anno)) {
    table <- cbind(table, anno)
  }
  if (is.null(display.columns)) {
    display.columns <- colnames(table)
  }
  else {
    if (!(xlab %in% display.columns))
      display.columns <- c(display.columns, xlab)
    if (!(ylab %in% display.columns))
      display.columns <- c(display.columns, ylab)
    if (!("gene" %in% display.columns))
      display.columns <- c("gene", display.columns)
  }
  table <- data.frame(index = 0:(nrow(table) - 1), table)
  if (length(status.cols) != 3)
    stop("status.cols\n          arg must have exactly 3 elements for [downreg, notDE, upreg]")
  xData <- list(data = list(x = xlab, y = ylab, table = table,
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
