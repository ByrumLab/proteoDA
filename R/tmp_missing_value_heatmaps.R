# A temporary location for the missing value heatmaps as I figure out what to do
# with them and which report to put them in.


#' Create a missing values heatmap
#'
#' Makes a ComplexHeatmap object showing a heatmap of missing values in the
#' input data.
#'
#' @inheritParams plotTotInten
#' @param missing An input matrix describing the missing data. See internals
#'   of \code{\link{plotHeatmapsForReport}}.
#' @param column_sort How should the columns of the heatmap be sorted? Options
#'   are: "cluster"- sort by similarity in missing values, "group"- sort samples
#'   by group, "batch"- sort samples by batch.
#' @param groupColors A vector of colors to use for coloring samples by group
#' @param batchColors A vector of colors to use for coloring samples by batch.
#' @param legend Include a legend in the plot? Default is FALSE.
#'
#' @return A ComplexHeatmap object
#' @export
#'
#' @examples
#' # No examples yet.
#'
missingValueHeatmap <- function(missing,
                                groups,
                                batch = NULL,
                                column_sort = c("cluster", "group", "batch"),
                                groupColors,
                                batchColors,
                                sampleLabels = NULL,
                                legend = FALSE) {
  # Check args
  column_sort <- rlang::arg_match(column_sort)

  # Deal with batch and sample names if null
  if (is.null(batch)) {
    batch <- c(rep("1", ncol(missing)))
  }

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(missing)
  }

  # Set up the column ordering for the different options
  if (column_sort == "cluster") {
    # If clustering just order things as they are
    cluster <- T
    order <- seq_along(groups)
    title <- "Cluster"
  } else if (column_sort == "group") {
    # If ordering by group, sort by group then batch
    cluster <- F
    order <- order(groups, batch)
    title <- "Sorted By Group"
  } else if (column_sort == "batch") {
    # If ordering by batch, sort by bathc then group
    cluster <- F
    order <- order(batch, groups)
    title <- "Sorted By Batch"
  } else {
    cli::cli_abort("Invalid value for {.arg column_sort}: cannot be {.val {column_sort}}")
  }

  # Then order the groups/batches/etc
  sorted_groups <- groups[order]
  sorted_batches <- batch[order]
  sorted_labels <- sampleLabels[order]
  ordered_data <- missing[,order]

  # Set up heatmap annotation
  ColAnn <- ComplexHeatmap::HeatmapAnnotation(
    Sample = sorted_groups,
    Batch = sorted_batches,
    col = list(Sample = groupColors,
               Batch = batchColors),
    annotation_legend_param = list(Sample = list(title = "Group",
                                                 at = unique(sorted_groups),
                                                 labels = paste("", unique(sorted_groups))),
                                   Batch = list(title = "Batch",
                                                at = unique(sorted_batches),
                                                labels = paste("Batch", unique(sorted_batches)))),
    show_legend = legend
  )
  # Plot heatmap
  hm_clust <- ComplexHeatmap::Heatmap(
    ordered_data + 0,
    col = c("white", "black"),
    column_names_side = "top",
    column_title = title,
    show_row_names = FALSE,
    show_column_names = TRUE,
    name = "Status",
    column_names_gp = grid::gpar(fontsize=7),
    heatmap_legend_param = list(at = c(0, 1),
                                labels = c("Missing", "Valid")),
    show_heatmap_legend = legend,
    top_annotation = ColAnn,
    cluster_columns = cluster,
    column_labels = sorted_labels
  )

  # Return the heatmap
  return(hm_clust)
}


#' Make a set of missing value heatmaps
#'
#' Makes, and optionally saves, a set of missing value heatmaps using
#' \code{\link{missingValueHeatmap}}.
#'
#' @inheritParams plotTotInten
#' @param data A normalized data matrix. Generally, the "log2" slot of the normList
#'   slot in the list output by \code{\link{process_data}}.
#' @param dir The directory in which to save the plot, if saving. Default is the
#'   current working directory.
#' @param showAllProteins For missing data heatmaps, show all proteins (including
#'  ones with no missing data)? Default is FALSE.
#' @param save Should the plot be saved (as a .png)? Default is FALSE.
#'
#' @return A list of length 4, where the first element is the missing data
#'       matrix and the next three are the ComplexHeatmap plot objects sorted
#'       by cluster, group, and batch, respectively.
#' @export
#'
#' @examples
#' # No examples yet.
#'
plotHeatmapsForReport <- function(data,
                                  groups,
                                  batch = NULL,
                                  sampleLabels = NULL,
                                  dir = ".",
                                  showAllProteins = FALSE,
                                  save = FALSE) {
  # Process arguments
  if (is.null(batch)) {
    batch <- c(rep("1",ncol(data)))
  }

  if (is.null(sampleLabels)) {
    sampleLabels <- colnames(data)
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive=TRUE)
  }

  # Set up colors
  batchCol <- colorBatch(batch)
  groupCol <- colorGroup2(groups)

  #Prepare data
  missing <- !is.na(data)
  if (!showAllProteins) {
    complete = apply(missing, 1, all)
    completeNA = apply(!missing, 1, all)
    missing <- missing[!complete & !completeNA,]
  }



  # Set up width for
  if (length(groups) <100 ) {
    width <- round(0.0871*length(groups)^2 + 24.375*length(groups) + 473.02, 0)
    width_single <- width/200
    width_together <- width/72
  } else {
    width <- round(0.0035*length(groups)^2 + 10.035*length(groups) + 146.15, 0)
    width_single <- width/90
    width_together <- width/30
  }

  # If we're saving, make individual plots of each type for merging
  # with slightly different formatting than the plots alone.
  # Will only do this if saving (trying to display the three plots all together
  # in the Rstudio window is too small, no point in doing if save == F)
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot.png"),
                   units="in",
                   width = width_together,
                   height = 8,
                   res = 100,
                   pointsize = 8)

    # Make each plot individually, plotting together
    hm_clust <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "cluster",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = F)
    hm_group <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "group",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = F)
    hm_batch <- missingValueHeatmap(missing = missing,
                                    groups = groups,
                                    batch = batch,
                                    column_sort = "batch",
                                    sampleLabels = sampleLabels,
                                    groupColors = colorGroup2(groups),
                                    batchColors = colorBatch(batch),
                                    legend = T)

    # Giving these objects unique names in the
    # name slot, to avoid a warning. Doesn't seem to affect
    # the plot at all
    hm_clust@name <- "a"
    hm_group@name <- "b"
    hm_batch@name <- "c"

    ComplexHeatmap::draw(hm_clust + hm_group + hm_batch,
                         heatmap_legend_side = "right",
                         annotation_legend_side = "right",
                         ht_gap = grid::unit(2, "cm"),
                         column_title = "Missing Values")
    grDevices::dev.off()
  }


  # Then, do individual plots.
  # Cluster samples
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_clust.png"),
                   units = "in",
                   width = width_single,
                   height = 8,
                   res = 100,
                   pointsize = 8)
  }
  hm_clust <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "cluster",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  # Will draw these whether we're saving or not
  ComplexHeatmap::draw(hm_clust)
  if (save) grDevices::dev.off()


  # Sample by group
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_group.png"),
                   units = "in",
                   width = width_single,
                   height = 8,
                   res = 100,
                   pointsize = 8)
  }

  hm_group <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "group",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  ComplexHeatmap::draw(hm_group)
  if (save) grDevices::dev.off()


  # Sample by batch
  if (save) {
    grDevices::png(filename = file.path(dir, "NaHMplot_batch.png"),
                   units = "in",
                   width = width_single,
                   height = 8,
                   res = 100,
                   pointsize = 8)
  }

  hm_batch <- missingValueHeatmap(missing = missing,
                                  groups = groups,
                                  batch = batch,
                                  column_sort = "batch",
                                  sampleLabels = sampleLabels,
                                  groupColors = colorGroup2(groups),
                                  batchColors = colorBatch(batch),
                                  legend = T)
  ComplexHeatmap::draw(hm_batch)
  if (save) grDevices::dev.off()



  data2 <- list(missing = missing,
                hm_clust = hm_clust,
                hm_group = hm_group,
                hm_batch = hm_batch)

  return(invisible(data2))
}
