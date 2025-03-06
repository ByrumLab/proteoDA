
#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param DAList PARAM_DESCRIPTION
#' @param grouping_column PARAM_DESCRIPTION
#' @param label_column PARAM_DESCRIPTION, Default: NULL
#' @param group_color PARAM_DESCRIPTION, Default: NULL
#' @param top PARAM_DESCRIPTION, Default: 500
#' @param standardize PARAM_DESCRIPTION, Default: TRUE
#' @param pca_axes PARAM_DESCRIPTION, Default: c(1, 2)
#' @param text.sizes PARAM_DESCRIPTION, Default: c(14, 14, 12, 10)
#' @param point.size PARAM_DESCRIPTION, Default: 3
#' @param label.size PARAM_DESCRIPTION, Default: 4
#' @param max.labels PARAM_DESCRIPTION, Default: 100
#' @param legend.position PARAM_DESCRIPTION, Default: 'right'
#' @param show.plot PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[proteoDA]{validate_DAList}}, \code{\link[proteoDA]{make_factor}}, \code{\link[proteoDA]{character(0)}}
#'  \code{\link[cli]{cli_abort}}
#' @rdname qc_pca_plot7_subgroups
#' @export
#' @importFrom proteoDA validate_DAList make_factor binfcolors
#' @importFrom cli cli_abort
qc_pca_plot7_subgroups <- function(DAList,
                                   grouping_column,
                                   label_column = NULL,
                                   group_color = NULL,
                                   top = 500,
                                   standardize = TRUE,
                                   pca_axes = c(1,2),
                                   text.sizes = c(14,14,12,10),
                                   point.size = 3,
                                   label.size = 4,
                                   max.labels = 100,
                                   legend.position = "right",
                                   show.plot=TRUE,
                                   max_pc = NULL){


 DAList <- proteoDA:::validate_DAList(x = DAList)

 grouping_column <- check_grouping_column(metadata = DAList$metadata,
                                          grouping_column = grouping_column)
 groups <- DAList$metadata[, grouping_column]

 if(!is.factor(groups)){
  groups <- proteoDA:::make_factor(x = groups)
 } else {
  groups <- droplevels(x = groups)
 }

 if(is.null(label_column)){
  sample_labels <- colnames(DAList$data)
 } else {
  label_column <- check_label_column(metadata = DAList$metadata,
                                     label_column = label_column)
  sample_labels <- DAList$metadata[, label_column]
 }

 if(is.null(group_color)){
  group_color <- proteoDA:::binfcolors[1]
 }

 if((length(group_color) != 1L)){
  cli::cli_abort("The group_color argument should be one color.")
 }


 out <- lapply(levels(groups), function(x){

  plot_colors <- rep("gray85", length(unique(groups)))
  names(plot_colors) <- levels(groups)

  idx <- which(names(plot_colors) %in% x)
  plot_colors[idx] <- group_color

  sample_labels2 <- rep("", length(sample_labels))
  sample_labels2[which(groups %in% x)] <- sample_labels[which(groups %in% x)]


  qc_pca_plot7(data            = DAList$data,
               groups          = groups,
               sample_labels   = sample_labels2,
               text.sizes      = text.sizes,
               point.size      = point.size,
               label.size      = label.size,
               max.labels      = max.labels,
               top             = top,
               standardize     = standardize,
               pca_axes        = pca_axes,
               colors          = plot_colors,
               title           = paste0("colored by ", grouping_column, ": ", x),
               legend.position = legend.position,
               show.plot       = show.plot,
               max_pc          = max_pc)


  })
 names(out) <- levels(groups)

 out

}











#' @title FUNCTION_TITLE
#' @description FUNCTION_DESCRIPTION
#' @param data A data frame of intensity data, likely normalized. Rows should
#' be proteins and columns should be samples.
#' @param groups A character or factor vector, listing the group(s) the samples
#' belong to., Default: NULL
#' @param sample_labels Optional, a vector of sample labels to use. If not
#' supplied, defaults to using the column names in the data, Default: NULL
#' @param text.sizes PARAM_DESCRIPTION, Default: c(14, 14, 12, 10)
#' @param point.size PARAM_DESCRIPTION, Default: 3
#' @param label.size PARAM_DESCRIPTION, Default: 4
#' @param max.labels PARAM_DESCRIPTION, Default: 100
#' @param top The number of most variable proteins to use for the PCA and
#' dendrogram clustering. Default is 500., Default: 500
#' @param standardize Should input data be standardized to a mean of 0 and
#' std.dev of 1 before performing PCA and dendrogram clustering? If input
#' data are not yet standardized, should be TRUE. , Default: TRUE
#' @param pca_axes A numeric vector of length 2 which lists the PC axes to plot.
#' Default is c(1,2), to plot the first two principal components., Default: c(1, 2)
#' @param colors PARAM_DESCRIPTION, Default: NULL
#' @param title PARAM_DESCRIPTION, Default: NULL
#' @param legend.position PARAM_DESCRIPTION, Default: 'right'
#' @param show.plot PARAM_DESCRIPTION, Default: TRUE
#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[rlang]{arg_match}}
#'  \code{\link[cli]{cli_abort}}
#'  \code{\link[proteoDA]{make_factor}}, \code{\link[proteoDA]{colorGroup}}, \code{\link[proteoDA]{rowVars}}
#'  \code{\link[stats]{prcomp}}
#'  \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{aes}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#'  \code{\link[ggplotify]{as.ggplot}}
#' @rdname qc_pca_plot7
#' @export
#' @importFrom rlang arg_match
#' @importFrom cli cli_abort
#' @importFrom proteoDA make_factor colorGroup rowVars
#' @importFrom stats prcomp
#' @importFrom ggplot2 ggplot aes
#' @importFrom ggrepel geom_label_repel
#' @importFrom ggplotify as.ggplot
qc_pca_plot7 <- function (data, ## matrix normalized data
                          groups = NULL,
                          sample_labels = NULL,
                          text.sizes = c(14,14,12,10), # title, axis.title, axis.label, legend
                          point.size = 3, ## size of points on graph
                          label.size = 4, ## size of point labels
                          max.labels = 100, ## if more than 50 samples in data don't show point labels
                          top = 500, ## top most variable proteins
                          standardize = TRUE, # z-score
                          pca_axes = c(1, 2),
                          colors = NULL, ## unique groups
                          title  = NULL,
                          max_pc = NULL, ## number PC to display in scree plot
                          legend.position = "right",
                          show.plot = TRUE) {



 ## if legend.position = none set show.legend = FALSE
 ## if legend.position = left,right,top,bottom set show.legend=TRUE
 legend.position <- rlang::arg_match(arg = legend.position,
                                     values= c("none", "top", "right", "bottom","left"),
                                     multiple = FALSE)
 show.legend <- ifelse(legend.position == "none", FALSE, TRUE)


 if(legend.position == "right"){
  # margin order t, r, b, l
  legend.margin = c(0, 0.5, 0, 0)
 }

 if(legend.position == "bottom"){
  # margin order t, r, b, l
  legend.margin = c(0.2, 0, 0.5, 0)
 }

 if(legend.position == "none"){
  legend.margin = c(0, 0, 0, 0)
 }

 if (is.null(title)) {
  title <- ""
 }

 if(is.null(groups)) {
  groups <- rep(1, ncol(data))
 }

 if(length(groups) != ncol(data)){
  cli::cli_abort("groups should be the same length as the number of columns in data")
 }

 if(!is.factor(groups)){
  groups <- proteoDA:::make_factor(as.character(groups), prefix=NULL)
 }
 groups <- droplevels(x = groups)


 ## get plot colors 1 per group levels i.e. unique groups
 if(is.null(colors)){
  colors <- proteoDA:::colorGroup(levels(groups))
 }

 if(!all(length(colors) == length(levels(groups)))){
  cli::cli_abort("length of colors and levels of groups do not match.")
 }
 names(colors) <- levels(groups)

 if(is.null(max.labels)){
  max.labels=ncol(data)
 }

 ## labels null set as col names; change col names to labels
 if(is.null(sample_labels)){
  sample_labels = colnames(data)
 }
 sample_labels <- proteoDA:::make_factor(as.character(sample_labels))

 ## meta data
 sample_group_info <- data.frame(ind = colnames(data),
                                 groups = groups,
                                 labels=sample_labels)


 ## top variable genes
 dat <- data[!apply(is.na(data), 1, any), ]
 dat <- dat[order(proteoDA:::rowVars(as.matrix(dat)), decreasing = TRUE), ]

 if(is.null(top)) { top <- nrow(dat) }
 top2 <- ifelse(nrow(dat) >= top, top, nrow(dat))
 dat <- dat[1:top2, ]


 ## center and scale
 ## col=samples; rows=genes
 if (standardize) {
  dat <- t(scale(x = t(dat), center = TRUE, scale = TRUE))
 }

 ## perform PCA analysis
 pca <- stats::prcomp(t(dat), scale = FALSE)
 pca$summary <- summary(pca)$importance

 ## df of components to plot x/y axis
 pca_data <- as.data.frame(pca$x[, pca_axes])
 pca_data$ind <- rownames(pca_data)
 colnames(pca_data) <- c("x", "y", "ind")

 #                  x          y     ind
 # HDLEC-1  -90.49757  25.548628 HDLEC-1
 # HDLEC-2  -86.00977   7.705991 HDLEC-2

 ## convert to long format
 plot_data <- merge(pca_data, sample_group_info, by = "ind");head(plot_data)

 #       ind          x           y group   labels
 # 1 HDLEC-1  0.2327423 -0.03698871  HDLEC   norm
 # 2 HDLEC-2  0.2484879 -0.07642031  HDLEC   norm


########################
##    PCA PLOT
########################
  p <-  ggplot2::ggplot(plot_data , ggplot2::aes(x = .data$x,
                                                 y = .data$y,
                                                 color = .data$groups,
                                                 label = .data$labels))

 # p <-  ggplot2::ggplot(data =  plot_data, aes(x = x, y = y, color = groups, label = labels)) +

 p <- p +  geom_point(show.legend = show.legend, size = point.size) +
           scale_color_manual(values = colors[groups], name = NULL) +

  ## size of points in legend
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  theme_gray() +

  ## axis labels / axis text
  theme(axis.title.x = element_text(size = text.sizes[2]),
        axis.text.x  = element_text(size = text.sizes[3], hjust = 0.5,
                                    margin = margin(t = 0.3, b = 0.15, unit = "cm")),
        axis.title.y = element_text(size = text.sizes[2]),
        axis.text.y  = element_text(size = text.sizes[3], vjust = 0.5,
                                    margin = margin(r = 0.3, l = 0.15, unit = "cm"))) +

  ## axis ticks, line width, lenght and color
  theme(axis.ticks.y = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
        axis.ticks.length.y = unit(0.25, "cm"),
        axis.ticks.x = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
        axis.ticks.length.x = unit(0.25, "cm")) +


  ## grid lines and border around plot area
  theme(panel.grid.major = element_line(linetype = "dotted", colour = "white", linewidth = 0.3),
        panel.grid.minor = element_line(linetype = "dotted", colour = "white", linewidth = 0.3),
        # panel.background = element_rect(color="gray80"),
        panel.border     = element_rect(linetype = "solid", fill = NA, colour = "gray90", linewidth = 0.3)) +

  theme(plot.title = element_text(size = text.sizes[1], face = "bold", hjust = 0)) +
  # geom_text(aes(label=labels), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(plot.margin     = unit(c(0.2,0.2,0.2,0.2), "cm")) +

  theme(legend.position = legend.position,
        legend.text     = element_text(size = text.sizes[4], vjust = 0.5)) +
  theme(legend.key.size = unit(0.5, 'cm'),
        legend.key = element_rect(fill = NA)) +
  # margin order t, r, b, l
  # theme(legend.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
  # theme(legend.margin = margin(c(0,0,0,5), unit = "cm")) +
  theme(legend.margin = margin(c(legend.margin), unit = "cm")) +

  # theme(plot.title = element_text(hjust = 0)) +
  xlab(paste0("PC", pca_axes[1], " (", round(pca$summary["Proportion of Variance", pca_axes[1]] * 100, 2), " %)")) +
  ylab(paste0("PC", pca_axes[2], " (", round(pca$summary["Proportion of Variance", pca_axes[2]] * 100, 2), " %)")) +
  # xlim(c(-15,20)) +
  # ylim(c(-10,10)) +
  ggtitle(title)
 p


 if (nrow(plot_data) < max.labels) {
  p <- p + ggrepel::geom_label_repel(show.legend = F,
                                     label.size = NA,
                                     fill = NA,
                                     size = label.size,
                                     max.time = 5,
                                     max.iter=1e9,
                                     max.overlaps=10,
                                     force=1,
                                     force_pull=0.5,
                                     label.padding=0.25,
                                     point.padding=0.25)
 }
 p <- ggplotify::as.ggplot(p)

 if(show.plot){  plot(p) }



 ####################
 ## SCREE PLOT
 ####################
 scree <- qc_pca_scree_plot7(pca = pca,
                             max_pc = max_pc,
                             text.sizes = text.sizes[1:3])

 # if(show.plot){ plot(scree$p) }


 invisible(list(p=p, pca=pca, plot_data = plot_data,
                sample_group_info=sample_group_info,
                param = list(top=top, standardize=standardize),
                scree = scree))

}







 ## pca = qc_pca_plot7() pca prcomp object in output list
 ## perform PCA analysis
 ## pca <- stats::prcomp(t(dat), scale = FALSE)
 ## pca$summary <- summary(pca)$importance
qc_pca_scree_plot7 <- function(pca,
                               max_pc = NULL,
                               text.sizes = c(16, 14, 14)){

  if(!is(pca, "prcomp")){
    cli::cli_abort("prcomp argument must be of class prcomp. see qc_pca_plot7 pca output.")
  }

  ## check text.sizes argument
  if(length(text.sizes) != 3L){
    cli::cli_abort(c("the text.sizes argument must be a vector of length 3
                   with positive numeric values (title, axis title, axis text)."))
  }

  is_cex <- sapply(text.sizes, check_num)
  if(!all(is_cex)){
    cli::cli_abort("the text.sizes argument must be a vector of length 3
                 with positive numeric values (title, axis title, axis text)")
  }


  ## max_pc is null set as total number of PC
  if(is.null(max_pc)){
    max_pc <- ncol(pca$summary)
  }

  ## check that max_pc is number >= 1
  if(!check_num(x = max_pc)){
    cli::cli_abort(c("{.arg max_pc} argument is invalid.",
                     "the max number of PCs displayed must be a number between
                  {.val {1}} and {.val {nrow(data)}}."))
  }

  ## if max_pc is larger than total number PCs set as total number PCs
  num_pc <- ifelse(max_pc > ncol(pca$summary), ncol(pca$summary), max_pc)



  ## get data for plotting
  data <- data.frame(PC = c(1:ncol(pca$summary)),
                     PropVar = c(pca$summary["Proportion of Variance", ]),
                     check.rows = F,
                     check.names = F,
                     fix.empty.names = F,
                     stringsAsFactors = F)
  data$PC <- proteoDA:::make_factor(as.character(data$PC))

  head(data)
  #     PC PropVar
  # PC1  1 0.38527
  # PC2  2 0.18913
  # PC3  3 0.14051


  ## subset data for number of desired PCs
  scree_data <- data[1:num_pc, ]

  # Basic line plot with points
  p <- ggplot2::ggplot(data=scree_data, ggplot2::aes(x = .data$PC, y = .data$PropVar, group = 1)) +
    geom_line(linetype = "solid", color = "black", linewidth = 1.2) +
    geom_point(size = 3) +

    ## axis labels / axis text
    theme(axis.title.x = element_text(size = text.sizes[2]),
          axis.text.x  = element_text(size = text.sizes[3], hjust = 0.5,
                                      margin = margin(t = 0.3, b = 0.15, unit = "cm")),
          axis.title.y = element_text(size = text.sizes[2]),
          axis.text.y  = element_text(size = text.sizes[3], vjust = 0.5,
                                      margin = margin(r = 0.3, l = 0.15, unit = "cm"))) +

    ## axis ticks, line width, lenght and color
    theme(axis.ticks.y = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.y = unit(0.25, "cm"),
          axis.ticks.x = element_line(linetype = "solid", colour = "black", linewidth = unit(0.4, 'cm')),
          axis.ticks.length.x = unit(0.25, "cm")) +


    ## grid lines and border around plot area
    theme(panel.grid.major = element_line(linetype = "dotted", colour = "white", linewidth = 0.4),
          panel.grid.minor = element_line(linetype = "dotted", colour = "white", linewidth = 0.4),
          # panel.background = element_rect(color="gray80"),
          panel.border     = element_rect(linetype = "solid", fill = NA, colour = "gray90",
                                          linewidth = 0.5)) +

    theme(plot.title = element_text(size = text.sizes[1], face = "bold", hjust = 0)) +
    # geom_text(aes(label=labels), position=position_dodge(width=0.9), vjust=-0.25) +
    theme(plot.margin     = unit(c(0.2,0.2,0.2,0.2), "cm")) +

    # theme(legend.position = legend.position,
    #       legend.text     = element_text(size = text.sizes[4], vjust = 0.5)) +
    # theme(legend.key.size = unit(0.5, 'cm'),
    #       legend.key = element_rect(fill = NA)) +
    # # margin order t, r, b, l
    # # theme(legend.margin = margin(t = 0, r = 0.5, b = 0, l = 0, unit = "cm")) +
    # # theme(legend.margin = margin(c(0,0,0,5), unit = "cm")) +
    # theme(legend.margin = margin(c(legend.margin), unit = "cm")) +

    # theme(plot.title = element_text(hjust = 0)) +
  xlab(paste0("Principle Component (PC)")) +
    ylab(paste0("Proportion of Variance")) +
    # xlim(c(-15,20)) +
    # ylim(c(-10,10)) +
    ggtitle(paste0("Scree Plot: Showing PC1 - PC",num_pc, " of ", nrow(data), " PCs"))
  p

  out <- list(p = p, scree_data = data, max_pc = max_pc)

  out


}
