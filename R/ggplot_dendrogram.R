#' Dendrogram for QC report
#'
#' Performs and then plots a hierarchical clustering analysis of sample intensities
#' across samples. Sample labels are colored according to the `groups` argument,
#' with colors determined by the \code{\link{colorGroup}} function. By default,
#' uses only the 500 most variable proteins for the analysis.
#'
#' @param data A numeric matrix or data frame with samples as columns and features (e.g. proteins) as rows.
#' @param groups A vector indicating sample groupings for color coding in the dendrogram.
#' @param sample_labels A vector of sample labels. If NULL, column names of `data` are used.
#' @param text.sizes A numeric vector specifying text sizes for title, leaf labels, and legend text.
#' @param point.size A numeric value specifying the size of the node circles.
#' @param top The number of the most variable proteins to use in the analysis.
#' @param standardize A logical value indicating whether to standardize data (mean = 0, SD=1).
#' @param dist_metric A character string specifying the distance metric (e.g., "euclidean", "manhattan").
#' @param clust_method A character string specifying the clustering method (e.g., "complete", "ward.D2").
#' @param colors A vector of colors corresponding to group levels.
#' @param legend.position A character string specifying legend position (e.g., "right", "bottom").
#' @param title A character string specifying the plot title. 
#' @param subtitle A character string specifying the plot subtitle. It can be NULL.
#' @param show.plot A logical value indicating whether to display the plot. 
#'
#' @return A ggplot list containing
#' \itemize{
#'   \item \code{p} - a ggplot dendrogram object.
#'   \item \code{hc} - The hierarchical clustering object
#'   \item \code{data_na} - the processed data used in the clustering
#'   \item \code{sample_group_info} - A data frame containing sample labels and groups
#'   \item \code{param} - A list of parameters used in the function
#'}   
#'   
#'
#' @importFrom ggplot2 aes theme element_text margin
#' @keywords internal
#' @export
#'
#' @examples
#' \dontrun{
#' den <- qc_dendrogram(data = results$data,
#'                  groups        = results$metadata$group,
#'                  sample_labels = results$metadata$sample,
#'                  top           = nrow(results$data),
#'                  standardize   = TRUE,
#'                  dist_metric   = "euclidean",
#'                  clust_method  = "complete",
#'                  colors        = all_colors2$group,
#'                  point.size    = 3,
#'                  text.sizes    = c(14,3,9),
#'                  legend.position = "right",
#'                  title         = "",
#'                  show.plot     = FALSE)
#'
#' }

qc_dendrogram <- function(data,
                          groups = NULL,           # vector
                          sample_labels = NULL,    # vector
                          text.sizes = c(14, 3, 12),  # title, leaf cex, legend cex
                          point.size = 3,          # size of node circle
                          top = 500,
                          standardize = TRUE,
                          dist_metric = "euclidean",
                          clust_method = "complete",
                          colors = NULL,           # vector length = levels(groups)
                          legend.position = "right",
                          title = NULL,
                          subtitle = NULL,
                          show.plot = TRUE) {
  
  # Require ggtree only if this function is used (ggtree in Suggests)
  if (!requireNamespace("ggtree", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg ggtree} is required for {.fn qc_dendrogram}. Please install it."
    )
  }
  
  # Check arguments
  dist_metric <- rlang::arg_match(
    arg    = dist_metric,
    values = c("euclidean", "maximum", "manhattan",
               "canberra", "binary", "minkowski"),
    multiple = FALSE
  )
  
  clust_method <- rlang::arg_match(
    arg    = clust_method,
    values = c("ward.D", "ward.D2", "single",
               "complete", "average", "mcquitty",
               "median", "centroid"),
    multiple = FALSE
  )
  
  if (is.null(title))    title    <- ""
  if (is.null(subtitle)) subtitle <- ""
  
  if (is.null(groups)) {
    groups <- rep(1, ncol(data))
  }
  groups <- make_factor(as.character(groups), prefix = NULL)
  stopifnot(length(groups) == ncol(data))
  
  # get plot colors 1 per group level
  if (is.null(colors)) {
    colors <- colorGroup(levels(groups))
  }
  if (!all.equal(length(colors), length(levels(groups)))) {
    stop("length of colors and levels of groups do not match.")
  }
  names(colors) <- levels(groups)
  
  # labels: default to column names
  if (is.null(sample_labels)) {
    sample_labels <- colnames(data)
  }
  if (!all(length(sample_labels) == ncol(data))) {
    stop("vector of sample_labels < the number of data columns.")
  }
  sample_labels <- make_factor(as.character(sample_labels))
  
  # metadata for plotting
  sample_group_info <- data.frame(
    ind    = colnames(data),
    groups = groups,
    labels = sample_labels
  )
  
  # top variable proteins
  data <- data[!apply(is.na(data), 1, any), ]
  data <- data[order(rowVars(as.matrix(data)), decreasing = TRUE), ]
  
  if (is.null(top)) top <- nrow(data)
  top2 <- ifelse(nrow(data) >= top, top, nrow(data))
  data <- data[1:top2, ]
  
  # center/scale rows (proteins)
  if (standardize) {
    data <- t(scale(x = t(data), center = TRUE, scale = TRUE))
  }
  
  # clustering
  hc <- stats::hclust(stats::dist(t(data), method = dist_metric),
                      method = clust_method)
  
  # dendrogram plot
  p <- ggtree::ggtree(hc)
  p <- ggtree::`%<+%`(p, sample_group_info) +
    ggtree::layout_dendrogram() +
    ggtree::geom_tippoint(aes(color = .data$groups), size = point.size) +
    ggtree::geom_tiplab(
      aes(color = .data$groups, label = .data$labels),
      size = text.sizes[2],
      angle = 90,
      hjust = 1,
      offset = -0.3,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = colors[groups],
      limits = levels(groups),
      name   = NULL
    ) +
    ggtree::theme_dendrogram(plot.margin = margin(6, 6, 110, 6)) +
    theme(
      plot.title   = element_text(hjust = 0, size = text.sizes[1]),
      legend.position = legend.position,
      legend.margin   = margin(30, 0, 0, 0),
      legend.text     = element_text(size = text.sizes[3])
    ) +
    ggtitle(label = title, subtitle = subtitle)
  
  # increase point size in legend
  p <- p + guides(colour = guide_legend(override.aes = list(size = point.size)))
  
  if (show.plot) {
    plot(p)
  }
  
  invisible(list(
    p    = p,
    hc   = hc,
    data_na = data,
    sample_group_info = sample_group_info,
    param = list(
      groups       = groups,
      top          = top,
      dist_metric  = dist_metric,
      clust_method = clust_method
    )
  ))
}

#' Dendrogram for QC Report with Subgroups
#' 
#' Generates dendrograms for different subgroups within a dataset. 
#' Each dendrogram highlights one subgroup while keeping others in grayscale. 
#' 
#' @param DAList A list containing data and metadata. 
#' @param grouping_column A character string specifying the sample metadata column to use for subgrouping (e.g., "group","batch", etc)
#' @param label_column An optional character string specifying the sample metadata column for sample labels. 
#' @param text.sizes A numeric vector specifying the text sizes for the title, leaf labels, and legend text. 
#' @param point.size A numeric vector specifying the size of the node circles. 
#' @param group_color A character string specifying the highlight color for the subgroup. 
#' @param top The number of the most variable proteins to use in the analysis. 
#' @param standardize A logical value indicating whether to standardize data (e.g., mean = 0, SD= 1).
#' @param dist_metric A character string specifying the distance metric (e.g., "euclidean", "manhattan").
#' @param clust_method A character string specifying the clustering method (e.g, "complete", "ward.D2")
#' @param legend.position A character string specifying the legend position (e.g., "right", "left", "top", "bottom").
#' @param show.plot A logical value indicating whether to display the plot. 
#' 
#' @return A list of dendrogram plots for each subgroup. 
#'
#' @importFrom ggplot2 aes theme element_text margin
#' @export
#'
#' @examples
#' \dontrun{
#' # An example of plot
#' den-> qc_dendrogram_subgroups(DAList = results, 
#'                                grouping_column = 

#' den <- qc_dendrogram(DAList           = results,
#'                     grouping_column  = results$metadata$group,
#'                     label_column     = results$metadata$sample,
#'                     top           = nrow(results$data),
#'                     standardize   = TRUE,
#'                     dist_metric   = "euclidean",
#'                     clust_method  = "complete",
#'                     colors        = all_colors2$group,
#'                     point.size    = 3,
#'                     text.sizes    = c(14,3,9),
#'                     legend.position = "right",
#'                     title         = "",
#'                     show.plot     = FALSE)
#'
#' }


qc_dendrogram_subgroups <- function(DAList,
                                    grouping_column,
                                    label_column = NULL,
                                    text.sizes = c(12,4,12),
                                    point.size = 3,
                                    group_color = NULL,
                                    top = 500,
                                    standardize = TRUE,
                                    dist_metric = "euclidean",
                                    clust_method = "complete",
                                    legend.position = "right",
                                    show.plot = TRUE){


  ## validate DAList
  DAList <- validate_DAList(x = DAList)

  ## check grouping column is a column in metadata. if not a factor make a factor
  ## then define groups vector
  grouping_column <- check_grouping_column(metadata = DAList$metadata, grouping_column = grouping_column)
  groups <- DAList$metadata[, grouping_column]

  ## make groups a factor and drop empty levels
  if(!is.factor(groups)){
    groups <- make_factor(x = as.character(groups))
  }
  groups <- droplevels(x = groups)


  ## if label column is not defined use rownames of metadata. if supplied check label column
  ## is single value, all values are unique and not NA, or blank
  ## then define sample_labels vector
  if(is.null(label_column)){
    sample_labels <- rownames(DAList$metadata)
  } else {
    label_column <- check_label_column(metadata = DAList$metadata, label_column = label_column)
    sample_labels <- DAList$metadata[, label_column]
  }


  ## if group color is not defined set as blueberry from binfcolors
  ## if defined check that a single value is provided if not error
  if(is.null(group_color)){
    group_color <- binfcolors[1]
  }
  if((length(group_color) != 1L)){
    cli::cli_abort("The group_color argument should be one color.")
  }


  out <- lapply(levels(groups), function(x){

    ## define gray color palette for all groups then set group color for
    ## subgroup of interest that will be highlighted in plot
    plot_colors <- rep("gray", length(unique(groups)))
    names(plot_colors) <- levels(groups)

    idx <- which(names(plot_colors) %in% x)
    plot_colors[idx] <- group_color

    title <- paste("Colored by:", grouping_column, "-", x)
    subtitle <- paste("cluster method:", clust_method,
                      "/  distance metric:", dist_metric)

    qc_dendrogram(data = DAList$data,
                  groups = groups,
                  sample_labels = sample_labels,
                  top = top,
                  standardize = standardize,
                  dist_metric = dist_metric,
                  clust_method = clust_method,
                  colors = plot_colors,
                  point.size = point.size,
                  text.sizes = text.sizes,
                  legend.position = legend.position,
                  title = title,
                  subtitle = subtitle,
                  show.plot = show.plot)

  })
  names(out) <- levels(groups)

  invisible(out)
  # p <- lapply(out, function(x){
  #   x$p
  # })
  #
  # p

}



# theme(axis.text.x  = element_text(angle = 45, vjust = 0.9, hjust = 1, size=text.sizes[3]),
#       axis.title.x = element_blank(),
#       axis.title.y = element_text(size=text.sizes[2]),
#       axis.text.y  = element_text(size=text.sizes[3]),
#       plot.title   = element_text(size=text.sizes[1])) +
#
#   # geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25) +
#   theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm")) +
#   theme(legend.position = legend.position, legend.text = element_text(size=text.sizes[4]), legend.margin = margin(90,0.2,0.2,0.2)) +
#   theme(legend.key.size = unit(0.5, 'cm')) +
#   ggtitle(title)

