#' Evaluate a normalization metric
#'
#' Applies a metric (see \code{\link{norm_metrics}}) across a list of normalized
#' data matrices and outputs the results in a data frame for plotting.
#'
#' @param normList A named list of normalized data matrices. Generally, the "normList"
#'   slot of the list that is output by \code{\link{process_data}}.
#' @param grouping A character or factor vector, listing the group(s) the samples
#'   belong to.
#' @param metric The normalization metric to calculate. Can be "PCV",
#'   "PMAD", "PEV", or "COR". See \code{\link{norm_metrics}}.
#' @return A data frame containing the selected metric,
#'   ready to be used for plotting.
#'
#' @export
#'
#' @seealso \code{\link{norm_metrics}}
#'
#' @examples
#' # No examples yet
#'
eval_pn_metric_for_plot <- function(normList,
                                    grouping,
                                    metric = c("PCV", "PMAD", "PEV", "COR", "log2ratio")) {

  # check args
  metric <- rlang::arg_match(metric)

  plotData <- NULL
  for (norm_method in names(normList)) {
    metric_results <- do.call(what = metric, args= list(data = as.data.frame(normList[norm_method]), groups = grouping))
    if (metric == "log2ratio") { #log2ratio does not have the per-group info the other metrics do, only per-normalization method
      x <- data.frame(method = norm_method,
                      value = metric_results, row.names = NULL)
    } else {
      x <- data.frame(method = norm_method,
                      group = names(metric_results),
                      value = metric_results, row.names = NULL)
    }
    plotData <- rbind(plotData, x)
  }

  plotData$method <- factor(plotData$method, levels = names(normList))

  return(plotData)
}


#' Generic plotting functions for normalization report
#'
#' An internal set of functions that take in processed data ready for plotting and
#' create generic plots (e.g., means, violin plots, line plots). These generic plots
#' are then passed to metric-specific plotting functions (see \code{\link{pn_plots}})
#' to create the final plot objects.
#'
#' @param plotData A dataframe of data to be plotted, created with the
#'   \code{\link{eval_pn_metric_for_plot}} function.
#' @return A ggplot object of the plot.
#'
#' @name pn_plots_generic
#'
#' @examples
#' # No examples yet
#'

#' @rdname pn_plots_generic
#'
pn_mean_plot <- function(plotData) {

  means <- aggregate(value ~ method, data = plotData, mean)
  sds <- aggregate(value ~ method, data = plotData, sd)
  ns <- aggregate(value ~ method, data = plotData, FUN = length)

  summary <- data.frame(method = means$method,
                        mean = means$value,
                        se = sds$value/sqrt(ns$value))

  # make the plot
  result <- summary  %>%
    ggplot(aes(col = method)) +
    geom_pointrange(aes(x = method, y = mean,
                        ymin = mean - se,
                        ymax = mean + se),
                    pch = 18,
                    size = 1.15) +
    scale_color_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank())

  return(result)
}

#' @rdname pn_plots_generic
#'
pn_violin_plot <- function(plotData) {
  # make the plot
  result <- plotData  %>%
    ggplot(aes(fill = method)) +
    geom_violin(aes(x = method, y = value),
                draw_quantiles = c(0.5),
                col = "black") +
    scale_fill_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          panel.grid = element_blank())

  return(result)
}

#' @rdname pn_plots_generic
#'
pn_density_plot <- function(plotData) {
  # make the plot
  result <- plotData  %>%
    ggplot(aes(color = method, group = method)) +
    geom_vline(aes(xintercept = 0), col = "grey80") +
    geom_line(aes(x = value),
              stat = "density",
              na.rm = T,
              trim = T) + # need trim = T for the density calculation to happen per-group
    scale_color_manual(values = unname(binfcolors)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank())

  return(result)
}


#' Functions for plotting normalization metrics
#'
#' A set of functions that take in normalized sample data and a grouping factor
#' and calculate some metric of of variability, error, etc., that we can use
#' to evaluate normalization methods. See \code{\link{norm_metrics}} for explanations
#' of each metric.
#'
#' @inheritParams eval_pn_metric_for_plot
#' @param zoom Should the plot cover the full range of log2ratios, or zoom
#'   in around 0? Default is FALSE.
#' @param legend Should the plot include the legend? Default is TRUE.
#' @return A ggplot object of the plot.
#'
#' @name pn_plots
#'
#' @examples
#' # No examples yet
#'

#' @rdname pn_plots
#' @export
#'
pn_plot_PCV <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PCV") %>%
    pn_mean_plot(.) +
    ylab("Pooled Coefficient of Variation") +
    ggtitle("PCV")
}

#' @rdname pn_plots
#' @export
#'
pn_plot_PMAD <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PMAD") %>%
    pn_mean_plot(.) +
    ylab("Median Absolute Deviation") +
    ggtitle("PMAD")
}

#' @rdname pn_plots
#' @export
#'
pn_plot_PEV <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "PEV") %>%
    pn_mean_plot(.) +
    ylab("Pooled Estimate of Variance") +
    ggtitle("PEV")
}

#' @rdname pn_plots
#' @export
#'
pn_plot_COR <- function(normList, grouping) {
  eval_pn_metric_for_plot(normList,
                          grouping,
                          metric = "COR") %>%
    pn_violin_plot(.) +
    ylab("Intragroup Correlation") +
    ggtitle("COR")
}

#' @rdname pn_plots
#' @export
#'
pn_plot_log2ratio <- function(normList, grouping, zoom = F, legend = T) {

  # Get plot data
  plotData <- eval_pn_metric_for_plot(normList,
                                      grouping,
                                      metric = "log2ratio")

  # Find max density
  # and min/max of ratios
  # for zooming later
  maxY <- max(aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$y))$value)
  min_xlim <- 0.5 * min(aggregate(value ~ method, data = plotData, FUN = function(x) min(stats::density(x, na.rm = T)$x))$value)
  max_xlim <- 0.5 * max(aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$x))$value)

  # Build base plot
  base <- plotData %>%
    pn_density_plot(.) +
    xlab("") +
    ylab("Density") +
    ggtitle("Log2 ratio")

  # Do zooming
  if (zoom) {
    base <- base +
      coord_cartesian(xlim = c(-0.3, 0.3),
                      ylim = c(maxY - (0.5*maxY), maxY + (0.2*maxY)))
  } else {
    base <- base +
      coord_cartesian(xlim = c(min_xlim, max_xlim))
  }

  # Output baseplot with or without legend
  if (legend) {
    result <- base +
      theme(legend.justification = c(1,1), # sets upper-right corner as locating point
            legend.position = c(0.99, 0.99), # puts locating point in upper-right corner
            legend.title = element_blank(),
            legend.text = element_text(size = 8),
            legend.key.size = unit(0.35, "cm"),
            legend.background = element_rect(fill = NULL))
  } else {
    result <- base +
      theme(legend.position = "none")
  }

  result
}











# Make plots of total intensity for each sample
#
# Makes, and optionally saves, a set of plot showing the total intensity for
# each sample across the normalization methods.
#
# @inheritParams proteinormMetricBoxplot
# @param sampleLabels Optional, a set of sample labels to use. If not supplied,
#   defaults to using the column names in the normList.
# @param dir The directory in which to save the plot, if saving. Default is the
#   current working directory.
# @param save Should the plot be saved (as a .png)? Default is FALSE.
#
# @return A list, equal in length to the input normList, where each
#       element of the list gives total intensity for each sample for a
#       given normalization method.
#
# @export
#
# @examples
# # No examples yet

# plotTotInten <- function(normList,
#                          groups,
#                          batch = NULL,
#                          sampleLabels = NULL,
#                          dir = ".",
#                          save = FALSE) {
#   # Prep args
#   # Again, not sure we need to coerce to factor here
#   groups <- make_factor(groups)
#   if (is.null(batch)) {
#     batch <- c(rep("1",ncol(normList[[1]])))
#   }
#   batch <- make_factor(as.character(batch), prefix = NULL)
#   if (is.null(sampleLabels)) {
#     sampleLabels <- colnames(normList[[1]])
#   }
#
#   # Set up plotting area, variably by number of samples
#   if (length(groups) < 100) {
#     width <- round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02, 0)
#     height <- 800
#     ncols <- 3
#   } else {
#     width <- round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15, 0)
#     height <- 2400
#     ncols <- 1
#   }
#
#
#   # If saving, set up files
#   if (save) {
#     if (!dir.exists(dir)) {
#       dir.create(dir, recursive = TRUE)
#     }
#     grDevices::png(filename = file.path(dir, "TotIntenPlot.png"),
#                    units = "px",
#                    width = width,
#                    height = height,
#                    pointsize = 15)
#   }
#
#   # Collect current par() options that we're going to change,
#   # and set them back on exit
#   old_mar <- graphics::par()$mar
#   old_oma <- graphics::par()$oma
#   on.exit(graphics::par(mar = old_mar), add = TRUE)
#   on.exit(graphics::par(oma = old_oma), add = TRUE)
#   on.exit(graphics::layout(matrix(1)), add = TRUE)
#
#   # Set up plotting parameters
#   # Have to do this twice, can't have the pars above the png making
#   if (length(groups) < 100) {
#     graphics::par(oma = c(2, 1, 1, 1),
#                   mar = c(8, 5, 5, 2))
#   } else {
#     graphics::par(oma = c(1, 5, 5, 5),
#                   mar = c(8, 2, 2, 2))
#   }
#
#
#   # Reorder group affiliations
#   # Since groups is an ordered factor (see make_factor()),
#   # this puts them in the order of their levels, then sorts by name.
#   # Does not necessarily sort alphabetically.
#   group_order <- sort.int(groups, index.return = T)$ix
#
#   # Then, reorder the sample labels and groups
#   sampleLabels <- sampleLabels[group_order]
#   groups <- groups[group_order]
#   # Cols in data reordered below
#
#
#   # Make a plot for each element of normList
#   graphics::layout(matrix(1:9, ncol = ncols, byrow = TRUE))
#   barList <- NULL
#   for (i in names(normList)) {
#     barList[[i]] <- colSums(normList[[i]], na.rm = T)
#     graphics::barplot(barList[[i]][group_order],
#                       main = "",
#                       las = 2,
#                       yaxt = "n",
#                       cex.main = 1.5,
#                       cex.lab = 1.2,
#                       col = colorGroup2(groups)[groups],
#                       names.arg = sampleLabels)
#     graphics::title(main = i, font.main = 1, cex.main = 1.5, line = 2)
#     graphics::axis(side = 2, cex.axis = 1.2, las = 2)
#     if (i == "VSN") {
#       graphics::mtext(side = 2, text = "Total Intensity", line = 6, cex = 1.5)
#     }
#   }
#   names(barList) <- names(normList)
#
#
#   if (save) grDevices::dev.off()
#
#   return(invisible(barList))
# }
