#' Make interactive reports on differential abundance
#'
#' Creates and saves interactive HTML reports summarizing
#' differential abundance analyses for each contrast in the results slot of the DAList.
#' Creates one HTML report for each contrast. Also creates a subfolder containing static .pdf
#' versions of all interactive plots. Currently, overwrites any previous reports
#' or other files with the same name.
#'
#'
#' Users can modify some aspects of the report output. First, users can modify
#' the width and height of the interactive report with the height and width
#' arguments, specified in pixels.
#'
#' Users can also control what data are displayed in the interactive table of
#' the report via the table_columns argument. Users can supply a vector of
#' column names from the annotation data frame to include, these columns will be
#' displayed in the order provided (column names may be changed slightly if the
#' original names cause issues with the Javascript code used for plotting). By
#' default, only the uniprot_id column is displayed.
#'
#' Finally, users can use the title_column argument to change the title of the
#' protein intensity (the plots on the right of the report). By default, the
#' values in the uniprot_id column of the annotation are used. To avoid plotting
#' issues, the values in the user-provided column will be truncated to 15
#' characters for use in the title: these values must remain unique
#' after truncation. If a title_column is supplied, it will be added to the
#' table_columns so that it is displayed in the table as well.
#'
#'
#' @param DAList A DAList object, with statistical results.
#' @param grouping_column The name of the column in the metadata which
#'   gives information on how to group samples for the interactive
#'   abundance plot.
#' @param table_columns Optional: the name of the column(s) in the annotation
#'   data frame that will be included in the interactive table of statistical
#'   results in the report. These columns will also be displayed in the tooltips
#'   of the volcano and MD plots. By default, only the uniprot_id column is
#'   displayed. See Details for more
#'   information.
#' @param title_column Optional: the name of a column in the annotation data
#'   frame from which to take values to use as the title for the protein
#'   intensity plots in the report. If not supplied, will use the
#'   info in the uniprot_id. If supplied, the title column will also be
#'   displayed (without truncation) in the results table. See Details for more
#'   information.
#' @param output_dir The directory in which to create the reports and save the
#'   plot files. If not specified, will default to the current working directory.
#' @param tmp_subdir The subdirectory within the output directory in which to
#'   store temporary files. Deleted by default. Default is "tmp".
#' @param height The height of the interactive report objects, in pixels.
#'   Default is 1000.
#' @param width The width of the interactive report objects, in pixels.
#'   Default is 1000.
#'
#' @return Invisibly returns the input DAList.
#'
#' @importFrom ggplot2 ggsave
#' @export
#'
#' @examples
#' \dontrun{
#'   # Using defaults
#'   write_limma_plots(DAList,
#'                     grouping_column = "treatment")
#'
#'   # Adjust size of report
#'   write_limma_plots(DAList,
#'                     grouping_column = "treatment",
#'                     height = 1500,
#'                     width = 1500)
#'
#'  # Customize output directory
#'  write_limma_plots(DAList,
#'                    grouping_column = "treatment",
#'                    output_dir = "DA_results")
#'
#'  # Add titles to the intensity plot,
#'  # using the values in the protein_id column
#'  # in the annotation of the DAList,
#'  # and display additional columns in the table
#'  write_limma_plots(DAList,
#'                    grouping_column = "treatment",
#'                    table_columns = "description",
#'                    title_column = "protein_id")
#' }
#'
write_limma_plots <- function(DAList = NULL,
                              grouping_column = NULL,
                              table_columns = c("uniprot_id"),
                              title_column = NULL,
                              output_dir = NULL,
                              tmp_subdir = "tmp",
                              height = 1000,
                              width = 1000) {

  # TODO:
  # Add overwriting checks?

  # Check input arguments generally
  input_DAList <- validate_DAList(DAList)

  # Make sure there's a design matrix present already,
  # tell user to set it first if not
  if (is.null(DAList$results)) {
    cli::cli_abort(c("Input DAList does not have a results design",
                     "i" = "Run {.code DAList <- extract_DA_results(DAList, ~ formula)}"))
  }

  # If provided, check that grouping column exists in the target data frame
  # And set it
  if (is.null(grouping_column)) { # If no groups provided, abort
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  } else { # check that grouping column exists in the target data frame
    if (length(grouping_column) != 1) {
      cli::cli_abort(c("Length of {.arg grouping_column} does not equal 1",
                       "i" = "Only specify one column name for {.arg grouping_column}"))

    }
    if (grouping_column %notin% colnames(DAList$metadata)) {
      cli::cli_abort(c("Column {.arg {grouping_column}} not found in metadata of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$metadata)}."))
    }
    # And set it
    groups <- as.character(DAList$metadata[,grouping_column])
  }

  # check that the table columns exist in the annotation
  if (!all(table_columns %in% colnames(DAList$annotation))) {
    problem_columns <- table_columns[which(table_columns %notin% colnames(DAList$annotation))]
    cli::cli_abort(c("Column{?s} {problem_columns} not found in annotation of {.arg DAList}",
                     "i" = "Check the column names with {.code colnames(DAList$annotation)}."))
  }

  # If title column is supplied
  if (!is.null(title_column)) {
    # Check that it exists in the annotation
    if (!title_column %in% colnames(DAList$annotation)) {
      cli::cli_abort(c("Column {.arg {title_column}} not found in annotation of {.arg DAList}",
                       "i" = "Check the column names with {.code colnames(DAList$annotation)}."))
    }
    # Grab the possible titles and truncate them
    temp_keys <- stringr::str_trunc(DAList$annotation[,title_column], width = 20, side = "right")

    # check that they're not duplicated
    if (any(duplicated(temp_keys))) {
      cli::cli_abort(c("values in {.arg {title_column}} were not unique after truncating to 15 characters.",
                       "i" = "Use a different column with unique values."))
    }

    # And add title column to the display table
    table_columns <- unique(c(table_columns, title_column))
  }

  # check height and width
  if (any(!is.numeric(height),
          !(height > 0))) {
    cli::cli_abort(c("{.arg height} must be a numeric value greater > 0."))
  }

  if (any(!is.numeric(width),
          !(width > 0))) {
    cli::cli_abort(c("{.arg width} must be a numeric value greater > 0."))
  }

  # defauly output dir if not provided
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform("{.arg output_dir} argument is empty.")
    cli::cli_inform("Setting output directory to current working directory:")
    cli::cli_inform("{.path {output_dir}}")
  }

  # TODO: add output filename validation once we decide on directory format
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  # Capture original wd, and setup function to return to original wd
  # upon error or function exit if it has changed
  old_wd <- getwd()
  on.exit(expr = {
    if (getwd() != old_wd) {
      cli::cli_inform("Returning working directory to {.path {old_wd}}")
      setwd(old_wd)
    }
  }, add = T)

  # If output directory isn't current wd, change wd so we can copy
  # files into the right spot.
  if (output_dir != old_wd) {
    cli::cli_inform("Setting working directory to output directory:")
    cli::cli_inform("{.path {file.path(old_wd, output_dir)}}")
    setwd(output_dir)
  }

  # Copy templates over
  file.copy(from = system.file("report_templates/glimma_xy_plot.Rmd",
                               package = "proteoDA"),
            to = "plot_template.Rmd", overwrite = T)
  file.copy(from = system.file("report_templates/limma_report_per_contrast.Rmd",
                               package = "proteoDA"),
            to = "report_template.Rmd", overwrite = T)

  # once we create the files, ensure they're deleted if there's an error below
  on.exit(expr = {
    cli::cli_inform("Removing temporary files from {.path {output_dir}}")
    unlink(c("logo_higherres.png", "plot_template.Rmd", "report_template.Rmd", tmp_subdir), recursive = T, expand = F)
  }, add = T, after = F)

  # Set up static plot folder
  if (!dir.exists("static_plots")) {
    dir.create("static_plots")
  }

  # Prep to loop over contrasts

  contrast_count <- 1
  num_contrasts <- length(names(DAList$results))

  # Loop over contrasts, making static plots and reports for each
  for (contrast in names(DAList$results)) {

    # Prep data
    data <- prep_plot_model_data(DAList$results, contrast)
    counts <- DAList$data[rownames(data),]
    counts[which(is.na(counts))] <- -9 # reassign missing to -9, so we can filter out later when plotting in Vega
    anno <- DAList$annotation[rownames(data), ]

    if (any(stringr::str_detect(table_columns, "\\."))) {
      # If any of the column names in the table columns have periods
      # the javascript will be mad. Replace with underscores.
      tbl_cols <- which(colnames(anno) %in% table_columns)
      colnames(anno)[tbl_cols] <- stringr::str_replace_all(colnames(anno)[tbl_cols],
                                                          "\\.",
                                                          "_")
      table_columns <- stringr::str_replace_all(table_columns,
                                                "\\.",
                                                "_")
    }

    # Change unique ids if title column was specified
    # checks were done above to ensure these are OK.
    if (!is.null(title_column)) {
        rownames(counts) <- temp_keys
        rownames(anno) <- temp_keys
        rownames(data) <- temp_keys
    }
    cli::cli_inform("Writing report for contrast {contrast_count} of {num_contrasts}: {.val {contrast}}")
    # make and save static plots
    for (type in c("raw", "adjusted")) {
      volcano <- static_volcano_plot(data,
                                     lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
                                     pval_thresh = DAList$tags$DA_criteria$pval_thresh,
                                     contrast = contrast, pval_type = type)
      MD <- static_MD_plot(data,
                           lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
                           contrast = contrast, pval_type = type)
      ggsave(filename = file.path("static_plots", paste0(paste(contrast, "volcano", type, "pval", sep = "-"), ".pdf")),
             plot = volcano,
             height = 6,
             width = 7,
             units = "in")
      ggsave(filename = file.path("static_plots", paste0(paste(contrast, "MD", type, "pval", sep = "-"), ".pdf")),
             plot = MD,
             height = 6,
             width = 7,
             units = "in")
      rm(volcano, MD)
    }


    pval_hist <- static_pval_histogram(data = data, contrast = contrast)
    ggsave(filename = file.path("static_plots", paste0(contrast, "_pval-hist.pdf")),
           plot = pval_hist,
           height = 6,
           width = 11,
           units = "in")
    rm(pval_hist)


    # Appease R CMD check
    # Right now, all the knitr function calls are in my .Rmd templates
    # so R CMD check thinks I'm importing knitr for no reason.
    # Making a dummy function call here, might be a better way to get around this
    x <- knitr::rand_seed
    rm(x)

    # make and save report
    rmarkdown::render("report_template.Rmd",
                      knit_root_dir = getwd(),
                      intermediates_dir = tmp_subdir,
                      output_file = paste0(contrast, "_DA_report.html"),
                      quiet = T)
    contrast_count <- contrast_count + 1
  }

  # TODO: add checks that files exist?
  invisible(input_DAList)
}


#' Prepare per-contrast model data for plotting
#'
#' Internal function used to prepare a results data frame for both static and interactive
#' plots in reports.
#'
#' @param model_results The results slot of a DAList object.
#' @param contrast The name of the contrast for which to prep the model data.
#'
#' @return A data frame of model results for the given contrast.
#'
#' @keywords internal
#'
prep_plot_model_data <- function(model_results, contrast) {
  # Get just the contrast we're interested in,
  # rename cols,
  # convert missing values in sig cols to 0
  # and add factor columns for static plots
  data <- model_results[[contrast]]
  data$`P value` <- data$P.Value
  data$`Adjusted P value` <- data$adj.P.Val
  data$negLog10rawP <- -log(data$P.Value, 10)
  data$negLog10adjP <- -log(data$adj.P.Val, 10)
  data$sig.PVal <- ifelse(is.na(data$sig.PVal), 0, data$sig.PVal)
  data$sig.FDR <- ifelse(is.na(data$sig.FDR), 0, data$sig.FDR)
  data$sig.pval.fct <- factor(x = data$sig.PVal,
                              levels = c(-1, 0, 1),
                              labels = c("downReg", "nonDE", "upReg"))
  data$sig.FDR.fct <- factor(x = data$sig.FDR,
                             levels = c(-1, 0, 1),
                             labels = c("downReg", "nonDE", "upReg"))
  data
}

#' Make a DE Volcano plot.
#'
#' Internal function for plotting static versions of volcano plots.
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param lfc_thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). LogFC are base 2.
#' @param pval_thresh The p-value threshold used to determine significance
#'   (significant when p < pval_thresh).
#' @param contrast The contrast being plotted. Used for generating the plot title.
#' @param pval_type The type of p-value to plot. Can be "raw" or "adjusted".
#'
#' @importFrom ggplot2 ggplot geom_vline geom_hline geom_point aes xlab ylab scale_color_manual scale_alpha_manual theme_bw ggtitle scale_y_continuous
#'
#' @return A ggplot object.
#'
#' @keywords internal
#'
static_volcano_plot <- function(data, lfc_thresh, pval_thresh, contrast, pval_type) {

  base <- ggplot(data = data) +
    geom_vline(xintercept = lfc_thresh*c(-1,1),
               linetype = "dashed",
               color = "grey50") +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed",
               color = "grey50") +
    xlab("log FC") +
    scale_color_manual(values = c("downReg" = "#00bfff",
                                  "nonDE" = "#858585",
                                  "upReg" = "#ff3030"),
                       name = "DE status") +
    scale_alpha_manual(values = c("downReg" = 1,
                                  "nonDE" = 0.6,
                                  "upReg" = 1),
                       name = "DE status") +
    theme_bw() +
    ggtitle(paste0(stringr::str_replace(contrast, "_vs_", " vs "), ", ", pval_type, " p-value"))

  if (pval_type == "raw") {
    final <- base +
      geom_point(aes(x = .data$logFC,
                     y = -log10(.data$P.Value),
                     color = .data$sig.pval.fct,
                     alpha = .data$sig.pval.fct), na.rm = T) +
      scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(-log10(data$P.Value), na.rm = T)) + 1, by = 1)) +
      ylab("-log10(P)")
  } else if (pval_type == "adjusted") {
    final <- base +
      geom_point(aes(x = .data$logFC,
                     y = -log10(.data$adj.P.Val),
                     color = .data$sig.FDR.fct,
                     alpha = .data$sig.FDR.fct), na.rm = T) +
      scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(-log10(data$adj.P.Val), na.rm = T)) + 1, by = 1)) +
      ylab("-log10(adjP)")
  } else{
    stop("invalid value for pval_type")
  }
  final
}

#' Make a DE MD plot.
#'
#' Internal function for plotting static versions of MD plots.
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param lfc_thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). LogFC are base 2.
#' @param contrast The contrast being plotted. Used for generating the plot title.
#' @param pval_type The type of p-value to plot. Can be "raw" or "adjusted".
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot geom_hline geom_point aes xlab ylab theme_bw ggtitle scale_color_manual scale_alpha_manual
#' @keywords internal
#'
static_MD_plot <- function(data, lfc_thresh, contrast, pval_type) {

  base <- ggplot(data = data) +
    geom_hline(yintercept = lfc_thresh*c(-1,1),
               linetype = "dashed",
               color = "grey50") +
    xlab("average normalized intensity") +
    scale_color_manual(values = c("downReg" = "#00bfff",
                                  "nonDE" = "#858585",
                                  "upReg" = "#ff3030"),
                       name = "DE status") +
    scale_alpha_manual(values = c("downReg" = 1,
                                  "nonDE" = 0.6,
                                  "upReg" = 1),
                       name = "DE status") +
    ylab("log FC") +
    theme_bw() +
    ggtitle(paste0(stringr::str_replace(contrast, "_vs_", " vs "), ", ", pval_type, " p-value"))

  if (pval_type == "raw") {
    final <- base +
      geom_point(aes(x = .data$AveExpr,
                     y = .data$logFC,
                     color = .data$sig.pval.fct,
                     alpha = .data$sig.pval.fct), na.rm = T)
  } else if (pval_type == "adjusted") {
    final <- base +
      geom_point(aes(x = .data$AveExpr,
                     y = .data$logFC,
                     color = .data$sig.FDR.fct,
                     alpha = .data$sig.FDR.fct), na.rm = T)
  } else{
    stop("invalid value for pval_type")
  }
  final
}

#' Make a p-value histogram plot
#'
#' Internal function for plotting p-value histograms of raw- and adjusted p-values.
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param contrast The contrast being plotted. Used for generating the plot title.
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot geom_histogram aes xlim theme_bw xlab theme element_rect
#'
#' @keywords internal
#'
static_pval_histogram <- function(data, contrast) {
  output <- ggplot(data) +
    geom_histogram(aes(x = .data$P.Value), binwidth = 0.025, color = "black", na.rm = T) +
    xlim(c(-0.05,1.05)) +
    theme_bw() +
    xlab("raw P value") +
    theme(panel.border = element_rect(fill = NA, color = "grey30")) +
    ggplot(data) +
    geom_histogram(aes(x = .data$adj.P.Val), binwidth = 0.025, color = "black", na.rm = T) +
    xlim(c(-0.05,1.05)) +
    theme_bw() +
    xlab("adjusted P value") +
    theme(panel.border = element_rect(fill = NA, color = "grey30")) +
    patchwork::plot_annotation(
      title = stringr::str_replace(contrast, "_vs_", " vs ")
    )
  output
}
