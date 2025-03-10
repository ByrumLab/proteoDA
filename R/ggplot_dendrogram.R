


#' Dendrogram for QC report
#'
#' Performs and then plots a hierarchical clustering analysis of sample intensities
#' across samples. Sample labels are colored according to the "groups" argument,
#' with colors determined by the \code{\link{colorGroup}} function. By default,
#' uses only the 500 most variable proteins for the analysis.
#'
#' @inheritParams qc_mds_plot
#' @inheritParams write_qc_report
#'
#' @return A ggplot object of the plot.
#' @keywords internal
#'
#' @importFrom ggplot2 aes theme element_text margin
#' @importFrom ggtree %<+% ggtree layout_dendrogram geom_tippoint geom_tiplab theme_dendrogram
#'
qc_dendrogram <- function(data,
                          groups = NULL, ## vector
                          sample_labels = NULL, ## vector
                          text.sizes = c(14,3,12), ## title, leaf cex, legend cex
                          point.size = 3, ## size of node circcle
                          top = 500,
                          standardize = TRUE,
                          dist_metric = "euclidean",
                          clust_method = "complete",
                          colors = NULL, ## vector length = levels gorups
                          legend.position = "right",
                          title = NULL,
                          subtitle = NULL,
                          show.plot = TRUE) {

  require(ggtree)
  # Check arguments
  dist_metric <- rlang::arg_match(
    arg = dist_metric,
    values = c(
      "euclidean", "maximum", "manhattan",
      "canberra", "binary", "minkowski"
    ),
    multiple = FALSE
  )

  clust_method <- rlang::arg_match(
    arg = clust_method,
    values = c(
      "ward.D", "ward.D2", "single",
      "complete", "average", "mcquitty",
      "median", "centroid"
    ), multiple = FALSE
  )



  if (is.null(title)) { title <- "" }
  if (is.null(subtitle)) { subtitle <- "" }

  if(is.null(groups)) { groups <- rep(1, ncol(data)) }
  groups <- proteoDA:::make_factor(as.character(groups), prefix=NULL)
  stopifnot(length(groups)==ncol(data))

  ## get plot colors 1 per group levels i.e. unique groups
  if(is.null(colors)){ colors <- proteoDA:::colorGroup(levels(groups)) }
  if(!all.equal(length(colors), length(levels(groups)))){
    stop("length of colors and levels of groups do not match.")
  }
  names(colors) <- levels(groups)



  # if (is.null(sample_labels)) {
  #   sample_labels <- colnames(data)
  # } else {
  #   colnames(data) <- sample_labels
  # }
  #   # Get sample info for colors
  # sample_group_info <- data.frame(text = sample_labels, group = groups)
  #


  ## labels null set as col names; change col names to labels
  if(is.null(sample_labels)){
    sample_labels <- colnames(data)
  }
  if(!all(length(sample_labels) == ncol(data))){
    stop("vector of samples_labels < the number of data columns.")
  }
  sample_labels <- proteoDA:::make_factor(as.character(sample_labels))
  ## meta data
  sample_group_info <- data.frame(ind = colnames(data), groups = groups,
                                  labels = sample_labels)

  ## top variable genes
  data <- data[!apply(is.na(data), 1, any), ]
  data <- data[order(proteoDA:::rowVars(as.matrix(data)), decreasing = TRUE), ]

  if(is.null(top)) { top <- nrow(data) }
  top2 <- ifelse(nrow(data) >= top, top, nrow(data))
  data <- data[1:top2, ]


  ## center/scale rows (proteins) mean=0;stdev=1
  if (standardize) {
    data <- t(scale(x = t(data), center = TRUE, scale = TRUE))
  }

  # Do clustering
  hc <- stats::hclust(stats::dist(t(data), method = dist_metric), method = clust_method)


  # Make plot
  p <- ggtree::ggtree(hc) %<+% sample_group_info +
        ggtree::layout_dendrogram() +
    ggtree::geom_tippoint(aes(color = .data$groups), size = point.size) +
    ggtree::geom_tiplab(aes(color = .data$groups, label = .data$labels),
                        size =text.sizes[2], ## cex leaf text
                        angle=90, hjust = 1, offset = -0.3, show.legend=FALSE) +
    scale_color_manual(values = colors[groups],#proteoDA:::colorGroup(groups),
                       limits = levels(groups), name = NULL) +
    ggtree::theme_dendrogram(plot.margin = margin(6,6,110,6)) +
    theme(plot.title = element_text(hjust = 0, size = text.sizes[1])) +
    theme(legend.position = legend.position,
          legend.margin = margin(30,0,0,0),
          legend.text = element_text(size = text.sizes[3])) +
    # theme(legend.margin = margin(0.2,0.2,60,0.2))+#text.sizes[4])) +
    ggtitle(label = title, subtitle = subtitle)



  ## increase point size of dots in legend
  p <- p +   guides(colour = guide_legend(override.aes = list(size = point.size)))
  if(show.plot){  plot(p) }


 invisible(list(p = p, hc = hc, data_na=data, sample_group_info=sample_group_info,
       param = list(groups=groups, top = top, dist_metric=dist_metric,
                    clust_method=clust_method)))

}






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
  DAList <- proteoDA:::validate_DAList(x = DAList)

  ## check grouping column is a column in metadata. if not a factor make a factor
  ## then define groups vector
  grouping_column <- check_grouping_column(metadata = DAList$metadata, grouping_column = grouping_column)
  groups <- DAList$metadata[, grouping_column]

  ## make groups a factor and drop empty levels
  if(!is.factor(groups)){
    groups <- proteoDA:::make_factor(x = as.character(groups))
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
    group_color <- proteoDA:::binfcolors[1]
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

