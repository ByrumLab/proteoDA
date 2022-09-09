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
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggplot2 geom_pointrange geom_violin geom_vline geom_line geom_point geom_hline geom_smooth
#' @importFrom ggplot2 scale_fill_manual scale_color_manual
#' @importFrom ggplot2 theme_bw theme element_text element_blank element_rect
#' @importFrom ggplot2 xlab ylab ggtitle coord_cartesian unit facet_wrap
#'
#' @examples
#' # No examples yet
#'
#' @rdname pn_plots_generic
#'
pn_mean_plot <- function(plotData) {

  means <- stats::aggregate(value ~ method, data = plotData, mean)
  sds <- stats::aggregate(value ~ method, data = plotData, stats::sd)
  ns <- stats::aggregate(value ~ method, data = plotData, FUN = length)

  summary <- data.frame(method = means$method,
                        mean = means$value,
                        se = sds$value/sqrt(ns$value))

  # make the plot
  result <- summary  %>%
    ggplot(aes(col = .data$method)) +
    geom_pointrange(aes(x = .data$method, y = .data$mean,
                        ymin = .data$mean - .data$se,
                        ymax = .data$mean + .data$se),
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
    ggplot(aes(fill = .data$method)) +
    geom_violin(aes(x = .data$method, y = .data$value),
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
    ggplot(aes(color = .data$method, group = .data$method)) +
    geom_vline(aes(xintercept = 0), col = "grey80") +
    geom_line(aes(x = .data$value),
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
    pn_mean_plot() +
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
    pn_mean_plot() +
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
    pn_mean_plot() +
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
    pn_violin_plot() +
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
  maxY <- max(stats::aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$y))$value)
  min_xlim <- 0.5 * min(stats::aggregate(value ~ method, data = plotData, FUN = function(x) min(stats::density(x, na.rm = T)$x))$value)
  max_xlim <- 0.5 * max(stats::aggregate(value ~ method, data = plotData, FUN = function(x) max(stats::density(x, na.rm = T)$x))$value)

  # Build base plot
  base <- plotData %>%
    pn_density_plot() +
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


#' @rdname pn_plots
#' @export
#'
pn_plot_MD <- function(normList, grouping) {

  # Assemble data for plotting
  log2ratios <- NULL
  for (norm_method in names(normList)) {
    metric_results <- do.call(what = "log2ratio", args= list(data = as.data.frame(normList[norm_method]),
                                                             groups = grouping, keep_protein_ID = T))

    long <- utils::stack(metric_results)
    long$protein  <- rep(rownames(metric_results), times = ncol(metric_results))
    long$method <- norm_method

    log2ratios <- rbind(log2ratios, long)
  }

  # Merge together the log2ratios for each pairwise comparison
  # and their average intensities
  plotData <- merge(log2ratios,
                    data.frame(mean_intensity = rowMeans(normList$log2, na.rm = T),
                               protein = rownames(normList$log2)))
  plotData$method <- factor(plotData$method, levels = names(normList))

  # make plot
  plotData %>%
  ggplot(aes(x = .data$mean_intensity, y = .data$values)) +
    geom_point(alpha = 0.3, na.rm = T) +
    geom_hline(yintercept = 0, color = "dodgerblue3") +
    geom_smooth(method = "lm", formula = "y ~ x", na.rm = T, color = "darkorange") +
    facet_wrap("method", ncol = 4) +
    ylab("logFC") +
    xlab("mean intensity (log2)") +
    theme_bw() +
    theme(panel.grid = element_blank())
}


