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
  xData <- Glimma:::buildXYData(table, status, main, display.columns,
                       anno, counts, xlab, ylab, status.cols, sample.cols, groups,
                       "none")
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
