#' Create quality control report
#'
#' Creates, and optionally saves as a PDF report, a bunch of plots which give
#' information on the distribution, clustering, and correlation of protein intensities
#' across samples for the chosen normalization method. By default, the PCA and
#' clustering analyses are performed on the top 500 most variable proteins. This
#' function is a wrapper that calls many subfunctions: \itemize{
#'   \item Makes violin plots of per-sample intensities
#'     with \code{\link{qc_violin_plot}}.
#'   \item Performs and plots a PCA with \code{\link{qc_pca_plot}}.
#'   \item Does hierarchical clustering and plots a dendrogram with
#'     \code{\link{qc_dendro_plot}}.
#'   \item Plots a correlation heatmap with \code{\link{qc_corr_hm}}.
#' } See the documentation of these subfunctions for more info.
#'
#'
#' @inheritParams write_proteinorm_report
#' @param chosen_norm_method The normalization method for which to perform the QC analysis.
#'   Should be a name of one of the datasets in normList.
#' @param label_column Optional. The name of column within the targets data frame
#'   which contains labels to use for plotting figures. When not supplied,
#'   defaults to using the column names of the data in processed_data.
#' @param file The file name of the report to be saved. Must end in .pdf. Will
#'   default to "QC_Report.pdf" if no filename is provided.
#' @param top_proteins The number of most variable proteins to use for the analysis.
#'   Default is 500.
#' @param standardize Should input data be standardized to a mean of 0 and std.dev of
#'   1? If input data are not yet standardized, should be TRUE. Default is TRUE.
#' @param pca_axes  A numeric vector of length 2 which lists the PC axes to plot.
#'   Default is c(1,2), to plot the first two principal components.
#' @inheritParams qc_dendro_plot
#' @inheritParams qc_missing_hm
#'
#' @return Invisibly returns a list with three slots: \enumerate{
#'   \item "plots"- A large list, where each element in the list is the returned
#'     object from the corresponding plotting function for that type of plot.
#'   \item "stats"- A dataframe with statistics on the data.
#'   \item "param"- A dataframe giving the parameters used for making the report.
#' }
#'
#' @export
#'
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' # No examples yet
#'

write_qc_report <- function(processed_data,
                            chosen_norm_method = NULL,
                            grouping_column = NULL,
                            label_column = NULL,
                            enrich = c("protein", "phospho"),
                            out_dir = NULL,
                            file = NULL,
                            overwrite = FALSE,
                            top_proteins = 500,
                            standardize = TRUE,
                            pca_axes = c(1,2),
                            dist_metric = "euclidean",
                            clust_method = "complete",
                            show_all_proteins = F) {

  cli::cli_rule()

  #################
  ## Check args  ##
  #################

  # Check that processed_data has expected list structure
  if (!all(c("normList", "targets", "filt", "param", "stats") %in% names(processed_data))) {
    cli::cli_abort(c("{.arg processed_data} does not have expected structure:",
                     "i" = "Is it the object created by running {.code process_data()}?."))
  }

  # enrich: possible values encoded in function def
  enrich <- rlang::arg_match(enrich)

  # chosen_norm_method
  chosen_norm_method <- rlang::arg_match(
    arg = chosen_norm_method,
    values = unique(c(
      names(processed_data$normList), "log2", "median", "mean", "vsn", "quantile",
      "cycloess", "rlr", "gi"
    )), multiple = FALSE
  )

  # If provided, check that grouping column exists in the target dataframe
  # And set it
  if (!is.null(grouping_column)) {
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))

    }
    if (grouping_column %notin% colnames(processed_data$targets)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in the targets dataframe of  in {.arg processed_data}",
                       "i" = "Check the column names with {.code colnames(processed_data$targets)}."))
    }
    groups <- as.character(processed_data$targets[,grouping_column])
  } else { # If no groups provided, set them but warn user
    groups <- rep("group", ncol(processed_data$normList[[chosen_norm_method]]))
    cli::cli_inform(cli::col_yellow("{.arg groups} argument is empty. Considering all samples/columns in {.arg processed_data} as one group."))
  }

  # If provided, check that label column is present in the target dataframe
  # and set it
  if (!is.null(label_column)) {
    if (length(label_column) != 1) {
      cli::cli_abort(c("Length of {.arg label_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg label_column}"))

    }
    if (label_column %notin% colnames(processed_data$targets)) {
      cli::cli_abort(c("Column {.arg {label_column}} not found in the targets dataframe of  in {.arg processed_data}",
                       "i" = "Check the column names with {.code colnames(processed_data$targets)}."))
    }
    sample_labels <- processed_data$targets[,label_column]
  } else { # Use colnames of the normlist if not provided
    sample_labels <- colnames(processed_data$normList[[chosen_norm_method]])
  }

  ############
  ## SET UP ##
  ############

  # Set default dir if not provided
  if (is.null(out_dir)) {
    out_dir <- file.path(paste0(enrich, "_analysis"), "01_quality_control")
    cli::cli_inform(cli::col_yellow("{.arg dir} argument is empty. Setting output directory to: {.path {out_dir}}"))
  }

  # Set default report name if not provided
  if (is.null(file)) {
    file <- "QC_Report.pdf"
    cli::cli_inform(cli::col_yellow("{.arg file} argument is empty. Saving report to: {.path {out_dir}/{file}}"))
  }

  # Make directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = T)
  }

  # Check that the filename is a pdf
  validate_filename(file, allowed_exts = c("pdf"))

  # Check if filename already exists
  if (file.exists(file.path(out_dir, file))) {
    if (overwrite) {
      cli::cli_inform("{.path {file}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting.")
    } else {
      cli::cli_abort(c("{.path {file}} already exists in {.path {out_dir}}",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Give {.arg file} a unique name or set {.arg overwrite} to {.val TRUE}"))
    }
  }

  # Select data from chosen normalization method for further use
  norm_data <- processed_data$normList[[chosen_norm_method]]

  #######################
  ## MAKE PLOT OBJECTS ##
  #######################

  # Violin plot
  violin_plot <- qc_violin_plot(data = norm_data,
                           groups = groups,
                           sample_labels = sample_labels) +
    ggtitle(paste0("Grouped by ", grouping_column))


  # PCA
  pca_plot <- qc_pca_plot(data = norm_data,
                          groups = groups,
                          sample_labels = sample_labels,
                          top_proteins = top_proteins,
                          standardize = standardize,
                          pca_axes = pca_axes) +
    ggtitle(paste0("PCA, colored by ", grouping_column))


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
                   "Colored by: ", grouping_column))

  # correlation heatmap
  correlation_heatmap <- qc_corr_hm(data = norm_data,
                                    groups = groups,
                                    sample_labels = sample_labels)

  # missing data heatmap - cluster by similarity
  miss_heatmap_cluster <- qc_missing_hm(data = norm_data,
                                        groups = groups,
                                        sample_labels = sample_labels,
                                        column_sort = "cluster",
                                        group_var_name = grouping_column,
                                        show_all_proteins = show_all_proteins)

  # missing value heatmap <- cluster by grouping column
  miss_heatmap_groups <- qc_missing_hm(data = norm_data,
                                        groups = groups,
                                        sample_labels = sample_labels,
                                        column_sort = "group",
                                        group_var_name = grouping_column,
                                        show_all_proteins = show_all_proteins)

  ###############################
  ## Save plots, check, return ##
  ###############################

  # When > 50 samples,
  # save plots individually on each page
  if (ncol(norm_data) > 50) {
    plots_list <- list(violin_plot,
                       pca_plot,
                       dendro_plot,
                       correlation_heatmap,
                       miss_heatmap_cluster,
                       miss_heatmap_groups)
  } else { # when <50 samples, group first plots together on 1 page
    plots_list <-  list(patchwork::patchworkGrob(violin_plot + pca_plot + dendro_plot),
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
    height <- 17
    width <- 17
  } else {
    height <- 8.5
    width <- 22
  }

  # Then save
  cli::cli_inform("Saving report to: {.path {file.path(out_dir, file)}}")
  ggsave(file.path(out_dir, file),
         plot = gridExtra::marrangeGrob(grobs = plots_list, nrow = 1, ncol = 1, top = NA),
         height = height,
         width = width,
         units = "in")

  if (!file.exists(file.path(out_dir, file))) {
    cli::cli_abort(c("Failed to create {.path {file.path(out_dir, file)}}"))
  }
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(file.path(out_dir, file))
}
