#' Create a quality control report
#'
#' Saves a PDF report containing a variety of plots which provide information on the
#' distribution, clustering, and correlation of protein intensities
#' across samples. See arguments for options for customizing the report.
#'
#' @inheritParams write_norm_report
#' @param color_column The name of the column in the metadata which
#'   gives information on how to color samples in plots within the report. If not
#'   supplied, all samples will be the same color.
#' @param label_column Optional. The name of column within the targets data frame
#'   which contains labels to use for plotting figures. When not supplied,
#'   defaults to using the column names of the data in processed_data. To ensure
#'   good plot formatting, sample labels will be truncated to 20 characters, and
#'   the function will give an error if the sample labels are not unique.
#' @param filename The file name of the report to be saved. Must end in .pdf. Will
#'   default to "QC_Report.pdf" if no filename is provided.
#' @param top_proteins The number of most variable proteins to use
#'  for the PCA and dendrogram clustering. Default is 500.
#' @param standardize Should input data be standardized to a mean of 0 and std.dev of
#'   1 before performing PCA and dendrogram clustering? If input data are
#'   not yet standardized, should be TRUE. Default is TRUE.
#' @param pca_axes  A numeric vector of length 2 which lists the PC axes to plot.
#'   Default is c(1,2), to plot the first two principal components.
#' @param dist_metric The metric used to define distance for dendrogram clustering.
#'   Default is "euclidean". See \code{\link[stats:dist]{stats::dist}} for options.
#' @param clust_method The agglomeration method to use for dendrogram clustering.
#'   Default is "complete", See \code{\link[stats:hclust]{stats::hclust}} for options.
#' @param show_all_proteins Should all proteins be shown in missing value heatmap,
#'  of only those with missing data? Default is F (only those with missing data).
#'
#' @return If report is created successfully, invisibly returns the input DAList.
#'
#' @export
#'
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' \dontrun{
#' # Color samples according to group identities
#' # in the "treatment" column of the metadata
#' write_qc_report(DAList,
#'                 color_column = "treatment")
#'
#' # Change the default directory and file names
#' write_qc_report(DAList,
#'                 color_column = "treatment",
#'                 output_dir = "my/chosen/directory",
#'                 filename = "my_report.pdf")
#'
#' # Overwrite an existing report
#' write_qc_report(DAList,
#'                 color_column = "treatment",
#'                 overwrite = T)
#'
#' # Customize PCA and clustering plots
#' write_qc_report(DAList,
#'                 color_column = "treatment",
#'                 top_proteins = 1000,
#'                 pca_aces = c(2,3),
#'                 dist_metric = "manhattan",
#'                 clust_method = "average")
#'
#' }

write_qc_report <- function(DAList,
                            color_column = NULL,
                            label_column = NULL,
                            output_dir = NULL,
                            filename = NULL,
                            overwrite = FALSE,
                            top_proteins = 500,
                            standardize = TRUE,
                            pca_axes = c(1,2),
                            dist_metric = "euclidean",
                            clust_method = "complete",
                            show_all_proteins = F) {

  #################################
  ## Check args and set defaults ##
  #################################

  input_DAList <- validate_DAList(DAList)

  if (!is.null(DAList$tags$normalized)) {
    if (!DAList$tags$normalized) {
      cli::cli_warn("Data in DAList are not normalized. Writing QC report for raw data.")
      normalized <- F
    } else {
      normalized <- T
    }
  } else {
    cli::cli_warn("Data in DAList are not normalized. Writing QC report for raw data.")
    normalized <- F
  }

  # If provided, check that grouping column exists in the target data frame
  # And set it
  if (!is.null(color_column)) {
    if (length(color_column) != 1) {
      cli::cli_abort(c("Length of {.arg color_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg color_column}"))

    }
    if (color_column %notin% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {color_column}} not found in the metadata of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$metadata)}."))
    }
    groups <- as.character(DAList$metadata[,color_column])
  } else { # If no groups provided, set them but warn user
    groups <- rep("group", ncol(DAList$data))
    cli::cli_inform(cli::col_yellow("{.arg groups} argument is empty. Considering all samples in {.arg DAList} as one group."))
  }

  # If provided, check that label column is present in the metadata
  # and set it
  if (!is.null(label_column)) {
    if (length(label_column) != 1) {
      cli::cli_abort(c("Length of {.arg label_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg label_column}"))

    }
    if (label_column %notin% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {label_column}} not found in the metadata of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$metadata)}."))
    }
    sample_labels <- as.character(DAList$metadata[,label_column])
    if (any(duplicated(sample_labels))) {
      cli::cli_abort(c("Sample labels in {label_column} column are not unique"))
    }
  } else { # Use colnames of the data if not provided
    sample_labels <- colnames(DAList$data)
  }

  # If necessary, truncate sample labels
  if (any(stringr::str_length(sample_labels) > 20)) {
    sample_labels <- stringr::str_trunc(sample_labels, width = 20, side = "right")
    num_too_long <- sum(stringr::str_length(sample_labels) > 20, na.rm = T)
    cli::cli_inform(c("{cli::qty(num_too_long)} {?A/Some} sample label{?s} longer than 20 characters",
                      "i" = "Truncating sample labels to 20 characters for plotting"))

    if (any(duplicated(sample_labels))) {
      cli::cli_abort(c("Sample labels are not unique after truncation",
                       "i" = "Supply a different column to use as sample labels with {.arg label_column}"))
    }
  }

  ############
  ## SET UP ##
  ############

  # Set default dir if not provided
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform("{.arg output_dir} argument is empty.")
    cli::cli_inform("Setting output directory to current working directory:")
    cli::cli_inform("{.path {output_dir}}")
  }

  # Set default report name if not provided
  if (is.null(filename)) {
    filename <- "QC_Report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg filename} argument is empty. Saving report to: {.path {output_dir}/{filename}}"))
  }

  # Make directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # Check that the filename is a pdf
  validate_filename(filename, allowed_exts = c("pdf"))

  # Check if filename already exists
  if (file.exists(file.path(output_dir, filename))) {
    if (overwrite) {
      cli::cli_inform("{.path {filename}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
    } else {
      cli::cli_abort(c("{.path {filename}} already exists in {.path {output_dir}}",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Give {.arg filename} a unique name or set {.arg overwrite} to {.val TRUE}"))
    }
  }

  # Select data from chosen normalization method for further use
  norm_data <- DAList$data

  #######################
  ## MAKE PLOT OBJECTS ##
  #######################

  # Violin plot
  violin_plot <- qc_violin_plot(data = norm_data,
                                groups = groups,
                                sample_labels = sample_labels) +
    ggtitle(paste0("Colored by ", color_column))

  # Change label if data aren't normalized
  if (!normalized) {
    violin_plot <- violin_plot +
      ylab("Intensity")
  }



  # PCA
  pca_plot <- qc_pca_plot(data = norm_data,
                          groups = groups,
                          sample_labels = sample_labels,
                          top_proteins = top_proteins,
                          standardize = standardize,
                          pca_axes = pca_axes) +
    ggtitle(paste0("PCA, colored by ", color_column))


  # dendrogram
  dendro_plot <- qc_dendro_plot(data = norm_data,
                                groups = groups,
                                sample_labels = sample_labels,
                                top_proteins = top_proteins,
                                standardize = standardize,
                                dist_metric = dist_metric,
                                clust_method = clust_method) +
    ggtitle(paste0("Cluster method: ", clust_method, "\n",
                   "Distance metric: ", dist_metric, "\n",
                   "Colored by: ", color_column))

  # correlation heatmap
  correlation_heatmap <- qc_corr_hm(data = norm_data,
                                    groups = groups,
                                    sample_labels = sample_labels)

  # missing data heatmap - cluster by similarity
  miss_heatmap_cluster <- qc_missing_hm(data = norm_data,
                                        groups = groups,
                                        sample_labels = sample_labels,
                                        column_sort = "cluster",
                                        group_var_name = color_column,
                                        show_all_proteins = show_all_proteins)

  # missing value heatmap <- cluster by grouping column
  miss_heatmap_groups <- qc_missing_hm(data = norm_data,
                                       groups = groups,
                                       sample_labels = sample_labels,
                                       column_sort = "group",
                                       group_var_name = color_column,
                                       show_all_proteins = show_all_proteins)

  ###############################
  ## Save plots, check, return ##
  ###############################

  # When > 50 samples,
  # save plots individually on each page
  if (ncol(norm_data) > 50) {
    plots_list <- list(pca_plot,
                       dendro_plot,
                       correlation_heatmap,
                       miss_heatmap_cluster,
                       miss_heatmap_groups)
  } else { # when <50 samples, group first plots together on 1 page
    joint_plot <- patchwork::wrap_plots(violin_plot, pca_plot, dendro_plot) &
      theme(plot.margin = margin(0,0,0,0))

    plots_list <-  list(patchwork::patchworkGrob(joint_plot),
                        correlation_heatmap,
                        miss_heatmap_cluster,
                        miss_heatmap_groups)
  }

  # Will adjust PDF plot sizes based on the number of samples
  # This also occurs in the correlation heatmap functions, where size has to be
  # specified at time of object creation (not of plotting).
  # If you make changes here, make corresponding changes in the
  # qc_corr_hm and qc_missing_hm functions
  if (ncol(norm_data) > 50) {
    # With > 50 samples, sample labels
    # are in 10 pt font. 72pt per inch
    # Dynamically size, with a minimum
    height <- max((10/72)*ncol(norm_data), 14)
    width <- max((10/72)*ncol(norm_data), 14)
  } else {
    height <- 10
    width <- 30
  }

  # Then save
  cli::cli_inform("Saving report to: {.path {file.path(output_dir, filename)}}")
  ggsave(file.path(output_dir, filename),
         plot = gridExtra::marrangeGrob(grobs = plots_list, nrow = 1, ncol = 1, top = NA),
         height = height,
         width = width,
         units = "in")

  if (!file.exists(file.path(output_dir, filename))) {
    cli::cli_abort(c("Failed to create {.path {file.path(output_dir, filename)}}")) #nocov
  }

  invisible(input_DAList)
}
