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
#' @param overwrite Should results files be overwritten? Default is FALSE.
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
                              overwrite = FALSE,
                              height = 1000,
                              width = 1000,
                              control_proteins = NULL,
                              highlight_by = "uniprot_id",
                              image_formats = c("pdf", "png")) {
  
  # Preserve optional slots (so we don't drop them on return)
  optional_slots <- c("filtered_proteins_per_contrast")
  optional_preserved <- DAList[intersect(names(DAList), optional_slots)]
  
  # Validate core DAList structure
  DAList <- validate_DAList(DAList)
  
  if (is.null(DAList$results)) {
    cli::cli_abort(c(
      "Input DAList does not have results",
      "i" = "Run {.code DAList <- extract_DA_results(DAList, ~ formula)}"
    ))
  }
  
  ## ---- Argument checks -----------------------------------------------------
  
  # grouping_column
  if (is.null(grouping_column)) {
    cli::cli_abort("{.arg grouping_column} cannot be empty")
  }
  if (length(grouping_column) != 1) {
    cli::cli_abort("{.arg grouping_column} must have length 1")
  }
  if (grouping_column %notin% colnames(DAList$metadata)) {
    cli::cli_abort("Grouping column {.val {grouping_column}} not found in metadata")
  }
  
  # table_columns must exist in annotation
  if (!all(table_columns %in% colnames(DAList$annotation))) {
    missing_cols <- table_columns[!table_columns %in% colnames(DAList$annotation)]
    cli::cli_abort("Table column(s) {.val {missing_cols}} not found in annotation")
  }
  
  # title_column: only check arguments & ensure it is included in table_columns.
  # Do NOT touch `data` or `anno` here.
  if (!is.null(title_column)) {
    if (length(title_column) != 1) {
      cli::cli_abort("{.arg title_column} must have length 1")
    }
    if (title_column %notin% colnames(DAList$annotation)) {
      cli::cli_abort("Title column {.val {title_column}} not found in annotation")
    }
    # Guarantee the title column is also in the interactive table
    table_columns <- unique(c(table_columns, title_column))
  }
  
  # height / width
  if (!is.numeric(height) || height <= 0) {
    cli::cli_abort("{.arg height} must be numeric > 0")
  }
  if (!is.numeric(width) || width <= 0) {
    cli::cli_abort("{.arg width} must be numeric > 0")
  }
  
  # output_dir
  if (is.null(output_dir)) {
    output_dir <- getwd()
    cli::cli_inform(
      "No {.arg output_dir} provided; using current working directory: {.path {output_dir}}"
    )
  }
  
  ## ---- Expected outputs & overwrite check ----------------------------------
  
  stopifnot(!is.null(DAList$results), length(DAList$results) > 0)
  contrast_names <- names(DAList$results)
  
  # Volcano + MD plots (raw/adjusted) in each format
  expected_xy_plots <- unlist(lapply(image_formats, function(ext) {
    file.path(
      output_dir, "static_plots",
      apply(
        X = expand.grid(
          contrast_names,
          c("volcano", "MD"),
          c("raw", "adjusted"),
          paste0("pval.", ext)
        ),
        MARGIN = 1,
        FUN = paste,
        collapse = "-"
      )
    )
  }))
  
  # p-value histograms
  expected_histograms_plots <- unlist(lapply(image_formats, function(ext) {
    file.path(
      output_dir, "static_plots",
      paste0(contrast_names, "-pval-hist.", ext)
    )
  }))
  
  # HTML interactive reports
  expected_reports <- file.path(
    output_dir,
    paste0(contrast_names, "_DA_report.html")
  )
  
  expected_results <- unique(c(
    expected_histograms_plots,
    expected_xy_plots,
    expected_reports
  ))
  
  existing <- file.exists(expected_results)
  if (any(existing, na.rm = TRUE)) {
    if (!overwrite) {
      cli::cli_abort(c(
        "Results files already exist",
        "!" = "and {.arg overwrite} == {.val {overwrite}}",
        "i" = "Change {.arg output_dir} or set {.arg overwrite} to {.val TRUE}"
      ))
    } else {
      cli::cli_inform(
        "Results files already exist, and {.arg overwrite} == {.val {overwrite}}. Overwriting results files."
      )
      unlink(expected_results[existing], recursive = FALSE, expand = FALSE)
    }
  }
  
  ## ---- Working directory + template setup ----------------------------------
  
  old_wd <- getwd()
  on.exit(
    expr = {
      if (getwd() != old_wd) {
        cli::cli_inform("Returning working directory to {.path {old_wd}}")
        setwd(old_wd)
      }
    },
    add = TRUE
  )
  
  # Ensure output/static_plots dirs exist
  if (!dir.exists(file.path(output_dir, "static_plots"))) {
    dir.create(file.path(output_dir, "static_plots"), recursive = TRUE)
  }
  
  # Change wd if output_dir is different
  if (output_dir != old_wd) {
    cli::cli_inform("Setting working directory to output directory:")
    cli::cli_inform("{.path {file.path(old_wd, output_dir)}}")
    setwd(output_dir)
  }
  
  # Choose template package (uams_internal flag) and copy Rmd
  if (!is.null(DAList$tags$uams_internal)) {
    template_package <- "proteoDAuams" # nocov
  } else {
    template_package <- "proteoDAstjude"
  }
  
  file.copy(
    from = system.file("report_templates/limma_report_per_contrast.Rmd",
                       package = template_package),
    to   = "report_template.Rmd",
    overwrite = TRUE
  )
  template_file <- "report_template.Rmd"
  
  # Clean up temp files on exit
  on.exit(
    expr = {
      cli::cli_inform("Removing temporary files from {.path {output_dir}}")
      unlink(c("CPM_Hz.png", "report_template.Rmd", tmp_subdir),
             recursive = TRUE, expand = FALSE)
    },
    add = TRUE,
    after = FALSE
  )
  
  ## ---- Shared annotation prep ----------------------------------------------
  
  shared_anno <- DAList$annotation
  
  # Deal with dots in table column names (JS compatibility)
  if (any(stringr::str_detect(table_columns, "\\."))) {
    tbl_cols <- which(colnames(shared_anno) %in% table_columns)
    colnames(shared_anno)[tbl_cols] <- stringr::str_replace_all(
      colnames(shared_anno)[tbl_cols],
      "\\.",
      "_"
    )
    internal_table_columns <- stringr::str_replace_all(table_columns, "\\.", "_")
  } else {
    internal_table_columns <- table_columns
  }
  
  ## ---- Loop over contrasts -------------------------------------------------
  
  contrast_count <- 1
  num_contrasts <- length(contrast_names)
  
  for (contrast in contrast_names) {
    
    # Per-contrast results
    data <- prep_plot_model_data(DAList$results, contrast)
    rownames(data) <- rownames(DAList$results[[contrast]])
    
    cols_to_display <- c(
      internal_table_columns,
      "average_intensity", "logFC", "p", "adjusted_p"
    )
    
    # Use per-contrast filtered data if available
    if (!is.null(DAList$data_per_contrast) &&
        contrast %in% names(DAList$data_per_contrast)) {
      counts <- DAList$data_per_contrast[[contrast]]
      groups <- DAList$metadata[colnames(counts), grouping_column, drop = TRUE]
    } else {
      counts <- DAList$data[rownames(data), , drop = FALSE]
      groups <- DAList$metadata[[grouping_column]]
    }
    
    # Align counts rows with data rows
    counts <- counts[rownames(data), , drop = FALSE]
    counts[is.na(counts)] <- -9
    
    # Annotation: per-contrast if present, else shared
    if (!is.null(DAList$annotation_per_contrast) &&
        contrast %in% names(DAList$annotation_per_contrast)) {
      anno <- DAList$annotation_per_contrast[[contrast]]
    } else {
      anno <- DAList$annotation[rownames(data), , drop = FALSE]
    }
    
    # Sanitize per-contrast annotation colnames (dots -> underscores)
    if (any(stringr::str_detect(colnames(anno), "\\."))) {
      colnames(anno) <- stringr::str_replace_all(colnames(anno), "\\.", "_")
    }
    
    # Ensure row alignment
    anno <- anno[rownames(data), , drop = FALSE]
    
    # Add p, adjusted_p, etc. for table/tooltip use
    anno$p                 <- round(data$P.Value,    digits = 4)
    anno$adjusted_p        <- round(data$adj.P.Val,  digits = 4)
    anno$logFC             <- data$logFC
    anno$movingSD          <- data$movingSD
    anno$logFC_zscore      <- data$logFC_z_scores
    anno$average_intensity <- data$average_intensity
    
    # ---- Title column per contrast (now that `anno` & `data` exist) ----
    if (!is.null(title_column)) {
      internal_title_col <- stringr::str_replace_all(title_column, "\\.", "_")
      if (internal_title_col %in% colnames(anno)) {
        aligned_titles <- anno[rownames(data), internal_title_col, drop = TRUE]
        data$internal_title_column <- stringr::str_trunc(
          aligned_titles,
          width   = 20,
          side    = "right",
          ellipsis = "..."
        )
      } else {
        cli::cli_abort(
          "title_column {.val {title_column}} not found in annotation for contrast {.val {contrast}}"
        )
      }
    } else {
      data$internal_title_column <- rownames(data)
    }
    
    cli::cli_inform(
      "Writing report for contrast {contrast_count} of {num_contrasts}: {.val {contrast}}"
    )
    
    ## ---- Static plots: volcano + MD ----------------------------------------
    
    for (type in c("raw", "adjusted")) {
      volcano <- static_volcano_plot(
        data,
        lfc_thresh  = DAList$tags$DA_criteria$lfc_thresh,
        pval_thresh = DAList$tags$DA_criteria$pval_thresh,
        contrast    = contrast,
        pval_type   = type,
        control_proteins = control_proteins,
        anno        = anno,
        highlight_by = highlight_by
      )
      
      MD <- static_MD_plot(
        data,
        lfc_thresh = DAList$tags$DA_criteria$lfc_thresh,
        contrast   = contrast,
        pval_type  = type
      )
      
      for (ext in image_formats) {
        out_volcano <- file.path(
          "static_plots",
          paste0(paste(contrast, "volcano", type, "pval", sep = "-"), ".", ext)
        )
        out_md <- file.path(
          "static_plots",
          paste0(paste(contrast, "MD", type, "pval", sep = "-"), ".", ext)
        )
        
        if (ext == "png") {
          ggplot2::ggsave(
            filename = out_volcano, plot = volcano,
            height = 6, width = 7, units = "in",
            dpi = 300, bg = "white"
          )
          ggplot2::ggsave(
            filename = out_md, plot = MD,
            height = 6, width = 7, units = "in",
            dpi = 300, bg = "white"
          )
        } else {
          ggplot2::ggsave(
            filename = out_volcano, plot = volcano,
            height = 6, width = 7, units = "in",
            bg = "white"
          )
          ggplot2::ggsave(
            filename = out_md, plot = MD,
            height = 6, width = 7, units = "in",
            bg = "white"
          )
        }
      }
      
      rm(volcano, MD)
    }
    
    ## ---- Static plots: p-value histograms ----------------------------------
    
    pval_hist <- static_pval_histogram(data = data, contrast = contrast)
    for (ext in image_formats) {
      out_hist <- file.path("static_plots", paste0(contrast, "-pval-hist.", ext))
      if (ext == "png") {
        ggplot2::ggsave(
          filename = out_hist, plot = pval_hist,
          height = 6, width = 11, units = "in",
          dpi = 300, bg = "white"
        )
      } else {
        ggplot2::ggsave(
          filename = out_hist, plot = pval_hist,
          height = 6, width = 11, units = "in",
          bg = "white"
        )
      }
    }
    rm(pval_hist)
    
    # Dummy knitr call to appease R CMD check
    x <- knitr::rand_seed
    rm(x)
    
    ## ---- Model equation label ----------------------------------------------
    
    design_label <- NULL
    if (!is.null(DAList$design_formula)) {
      design_label <- paste(deparse(DAList$design_formula), collapse = "")
    } else if (!is.null(DAList$design) && inherits(DAList$design, "formula")) {
      design_label <- paste(deparse(DAList$design), collapse = "")
    } else {
      design_label <- "<design formula not available>"
    }
    # Expose to Rmd as `design`
    design <- design_label
    
    ## ---- Render interactive report -----------------------------------------
    
    rmarkdown::render(
      input            = template_file,
      knit_root_dir    = getwd(),
      intermediates_dir = tmp_subdir,
      output_file      = paste0(contrast, "_DA_report.html"),
      quiet            = TRUE
    )
    
    contrast_count <- contrast_count + 1
  }
  
  ## ---- Check all expected files exist --------------------------------------
  
  if (any(!file.exists(stringr::str_remove(
    expected_results,
    pattern = paste0("^", output_dir, "/")
  )))) {
    failed <- expected_results[!file.exists(stringr::str_remove(
      expected_results,
      pattern = paste0("^", output_dir, "/")
    ))]
    cli::cli_abort(c(
      "Failed to write the following {cli::qty(length(failed))} results file{?s}:",
      "!" = "{.path {failed}}"
    ))
  }
  
  ## ---- Restore optional slots & return -------------------------------------
  
  DAList[names(optional_preserved)] <- optional_preserved
  class(DAList) <- "DAList"
  invisible(DAList)
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
static_volcano_plot <- function(data,
                                lfc_thresh,
                                pval_thresh,
                                contrast,
                                pval_type,
                                control_proteins = NULL,
                                anno = NULL,
                                highlight_by = "uniprot_id") {
  
  data$uniprot_id <- rownames(data)
    
  if (!is.null(control_proteins)) {
    if (is.null(anno)) {
      stop("If control_proteins is provided, anno must also be provided.")
    }
    
    if (!(highlight_by %in% colnames(anno))) {
      warning(sprintf("The highlight_by column '%s' is not in the annotation data. Skipping control protein highlighting.", highlight_by))
      data$highlight <- FALSE
      data$highlight_label <- NA_character_
    } else if (!("uniprot_id" %in% colnames(anno))) {
      stop("Annotation data must contain 'uniprot_id' to align with DE results.")
    } else {
      anno <- anno[rownames(data), , drop = FALSE]  # align
      highlight_vals <- anno[[highlight_by]]
      data$highlight <- highlight_vals %in% control_proteins
      data$highlight_label <- highlight_vals  # for labeling later
    }
  } else {
    data$highlight <- FALSE
    data$highlight_label <- NA_character_
  }
  
  if (pval_type == "raw") {
    data$label_y <- -log10(data$P.Value)
  } else if (pval_type == "adjusted") {
    data$label_y <- -log10(data$adj.P.Val)
  } else {
    stop("invalid value for pval_type")
  }
  
  base <- ggplot(data = data) +
    geom_vline(xintercept = lfc_thresh * c(-1, 1),
               linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pval_thresh),
               linetype = "dashed", color = "grey50") +
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
      geom_point(aes(x = logFC,
                     y = label_y,
                     color = sig.pval.fct,
                     alpha = sig.pval.fct), na.rm = TRUE) +
      scale_y_continuous(breaks = seq(0, ceiling(max(data$label_y, na.rm = TRUE)) + 1, 1)) +
      ylab("-log10(P)")
  } else {
    final <- base +
      geom_point(aes(x = logFC,
                     y = label_y,
                     color = sig.FDR.fct,
                     alpha = sig.FDR.fct), na.rm = TRUE) +
      scale_y_continuous(breaks = seq(0, ceiling(max(data$label_y, na.rm = TRUE)) + 1, 1)) +
      ylab("-log10(adjP)")
  }
  
  if (any(data$highlight, na.rm = TRUE)) {
    highlight_df <- data[data$highlight, , drop = FALSE]
    final <- final +
      geom_point(data = highlight_df,
                 aes(x = logFC, y = label_y),
                 shape = 21, fill = "yellow", color = "black", size = 2.5, stroke = 0.6) +
      geom_text(data = highlight_df,
                aes(x = logFC, y = label_y, label = highlight_label),
                size = 2.5, vjust = -0.8)
  }
  
  return(final)
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
      geom_point(aes(x = .data$average_intensity,
                     y = .data$logFC,
                     color = .data$sig.pval.fct,
                     alpha = .data$sig.pval.fct), na.rm = T)
  } else if (pval_type == "adjusted") {
    final <- base +
      geom_point(aes(x = .data$average_intensity,
                     y = .data$logFC,
                     color = .data$sig.FDR.fct,
                     alpha = .data$sig.FDR.fct), na.rm = T)
  } else{
    stop("invalid value for pval_type") #nocov
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
  p1 <- ggplot(data) +
    geom_histogram(aes(x = .data$P.Value), binwidth = 0.025, color = "black", na.rm = TRUE) +
    xlim(c(-0.05, 1.05)) +
    theme_bw() +
    xlab("raw P value") +
    theme(panel.border = element_rect(fill = NA, color = "grey30"))
  
  p2 <- ggplot(data) +
    geom_histogram(aes(x = .data$adj.P.Val), binwidth = 0.025, color = "black", na.rm = TRUE) +
    xlim(c(-0.05, 1.05)) +
    theme_bw() +
    xlab("adjusted P value") +
    theme(panel.border = element_rect(fill = NA, color = "grey30"))
  
  # Combine with patchwork
 # combined <- p1 + p2 + patchwork::plot_annotation(
#    title = stringr::str_replace(contrast, "_vs_", " vs "))
    
    combined <- patchwork::wrap_plots(p1, p2) +
      patchwork::plot_annotation(title = stringr::str_replace(contrast, "_vs_", " vs "))
    
  
  return(combined)
}
