#' Make interactive reports on DE
#'
#' A wrapper function that creates and saves interactive HTML reports summarizing
#' differential expression analyses for each contrast. Creates one HTML report
#' for each contrast. Also creates and saves a set of static .pdf plots. Current
#' behavior is to overwrite any previous reports or other files with the same name.
#'
#' @param model_results A model results object returned by
#'   \code{\link{extract_limma_DE_results}}.
#' @param annotation A dataframe containing annotation data for this analysis.
#'   In our pipeline, usually the "annot" slot from the object returned by
#'   \code{\link{extract_data}}.
#' @param groups A vector describing the groups to which each sample belongs. Used
#'   for grouping expression data in the interactive expression plot.
#' @param output_dir The directory in which to create the reports and save the
#'   plot files. No defaults, must be specified.
#' @param tmp_subdir The subdirectory within the output directory in which to
#'   store temporary files. Deleted by default. Default is "tmp".
#' @param height The height of the interactive report objects, in pixels.
#'   Default is 1000.
#' @param width The width of the interactive report objects, in pixels.
#'   Default is 1000.
#'
#' @return If successful, invisible(TRUE).
#' @export
#'
#' @import ggplot2
#'
#' @examples
#' # No examples yet
make_limma_reports <- function(model_results = NULL,
                               annotation = NULL,
                               groups = NULL,
                               output_dir = NULL,
                               tmp_subdir = "tmp",
                               height = 1000,
                               width = 1000) {

  # TODO: add output filename validation once we decide on directory format
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }

  old_wd <- getwd()
  on.exit(expr = setwd(old_wd), add = T)

  cli::cli_rule()

  cli::cli_inform("Setting working directory to {.path {file.path(old_wd, output_dir)}}")
  setwd(output_dir)

  cli::cli_inform("Copying resources to output directory {.path {output_dir}}")
  if (!dir.exists("resources")) {
    dir.create("resources")
  }
  if (!dir.exists("static_plots")) {
    dir.create("static_plots")
  }

  file.copy(from = system.file("report_templates/glimma_xy_plot.Rmd",
                               package = "proteomicsDIA"),
            to = "plot_template.Rmd", overwrite = T)
  file.copy(from = system.file("report_templates/limma_report_per_contrast.Rmd",
                               package = "proteomicsDIA"),
            to = "report_template.Rmd", overwrite = T)
  file.copy(from = system.file("report_templates/logo_higherres.png",
                               package = "proteomicsDIA"),
            to = "resources/logo_higherres.png", overwrite = T)

  contrast_count <- 1
  num_contrasts <- length(names(model_results$stats_by_contrast))

  # TODO:
  # Add overwriting checks?

  # Loop over contrasts, making static plots and reports for each
  for (contrast in names(model_results$stats_by_contrast)) {

    # Prep data
    data <- prep_plot_model_data(model_results, contrast)
    counts <- model_results$data[rownames(data),]
    counts[which(is.na(counts))] <- -9 # reassign missing to -9, so we can filter out later when plotting in Vega
    anno <- annotation[rownames(data), ]
    cli::cli_inform("Writing report for contrast {contrast_count} of {num_contrasts}: {.val {contrast}}")

    # make and save static plots
    for (type in c("raw", "adjusted")) {
      volcano <- static_volcano_plot(data,
                                     lfc.thresh = model_results$lfc.thresh,
                                     pval.thresh = model_results$pval.thresh,
                                     contrast = contrast, pval.type = type)
      MD <- static_MD_plot(data,
                           lfc.thresh = model_results$lfc.thresh,
                           contrast = contrast, pval.type = type)
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
    # TODO: find better way.
    x <- knitr::rand_seed
    rm(x)

    # make and save report
    rmarkdown::render("report_template.Rmd",
                      knit_root_dir = getwd(),
                      intermediates_dir = tmp_subdir,
                      output_file = paste0(contrast, "_DE_report.html"),
                      quiet = T)
    contrast_count <- contrast_count + 1
  }

  # Clean up and exit

  # TODO: add checks that files exist?

  cli::cli_inform("Removing temporary files from {.path {output_dir}}")
  unlink(c("logo_higherres.png", "plot_template.Rmd", "report_template.Rmd", tmp_subdir), recursive = T)
  cli::cli_inform("Returning working directory to {.path {old_wd}}")
  setwd(old_wd)
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(TRUE)
}




#' Prepare per-contrast model data for plotting
#'
#' Internal function used to prepare a results dataframe for both static and interactive
#' plots in our reports.
#'
#' @param model_results A model results object from
#'   \code{\link{extract_limma_DE_results}}.
#' @param contrast The name of the contrast for which to prep the model data
#'
#' @return A data frame of model results for the given contrast
#'
#' @examples
#' # No examples yet
prep_plot_model_data <- function(model_results, contrast) {
  # Get just the contrast we're interested in,
  # rename cols,
  # convert missing values in sig cols to 0
  # and add factor columns for static plots
  data <- model_results$stats_by_contrast[[contrast]]
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
#' @param lfc.thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). LogFC are base 2.
#' @param pval.thresh The p-value threshold used to determine significance
#'   (significant when p < pval.thresh).
#' @param contrast The contrast being plotted. Used for generating the plot title.
#' @param pval.type The type of p-value to plot. Can be "raw" or "adjusted".
#'
#' @return A ggplot object.
#'
#' @examples
#' # No examples yet
static_volcano_plot <- function(data, lfc.thresh, pval.thresh, contrast, pval.type) {

  base <- ggplot(data = data) +
    geom_vline(xintercept = lfc.thresh*c(-1,1),
               linetype = "dashed",
               color = "grey50") +
    geom_hline(yintercept = -log10(pval.thresh),
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
    ggtitle(paste0(stringr::str_replace(contrast, "_vs_", " vs "), ", ", pval.type, " p-value"))

  if (pval.type == "raw") {
    final <- with(data, { # with is just to appease R CMD check for global variables
      base +
      geom_point(aes_string(x = logFC,
                     y = -log10(P.Value),
                     color = sig.pval.fct,
                     alpha = sig.pval.fct), na.rm = T) +
      scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(-log10(data$P.Value), na.rm = T)) + 1, by = 1)) +
      ylab("-log10(P)")
    })
  } else if (pval.type == "adjusted") {
    final <- with(data, {
      base +
      geom_point(aes(x = logFC,
                     y = -log10(adj.P.Val),
                     color = sig.FDR.fct,
                     alpha = sig.FDR.fct), na.rm = T) +
      scale_y_continuous(breaks = seq(from = 0, to = ceiling(max(-log10(data$adj.P.Val), na.rm = T)) + 1, by = 1)) +
      ylab("-log10(adjP)")
    })
  } else{
    stop("invalid value for pval.type")
  }
  final
}

#' Make a DE MD plot.
#'
#' Internal function for plotting static versions of MD plots.
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param lfc.thresh The logFC threshold used to determine significance
#'   (significant when |logFC| > lfc.tresh). LogFC are base 2.
#' @param contrast The contrast being plotted. Used for generating the plot title.
#' @param pval.type The type of p-value to plot. Can be "raw" or "adjusted".
#'
#' @return A ggplot object
#'
#' @examples
#' # No examples yet
static_MD_plot <- function(data, lfc.thresh, contrast, pval.type) {

  base <- ggplot(data = data) +
    geom_hline(yintercept = lfc.thresh*c(-1,1),
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
    ggtitle(paste0(stringr::str_replace(contrast, "_vs_", " vs "), ", ", pval.type, " p-value"))

  if (pval.type == "raw") {
    final <- with(data, {
      base +
      geom_point(aes(x = AveExpr,
                     y = logFC,
                     color = sig.pval.fct,
                     alpha = sig.pval.fct), na.rm = T)
    })
  } else if (pval.type == "adjusted") {
    final <- with(data, {base +
      geom_point(aes(x = AveExpr,
                     y = logFC,
                     color = sig.FDR.fct,
                     alpha = sig.FDR.fct), na.rm = T)
    })
  } else{
    stop("invalid value for pval.type")
  }
  final
}

#' Make a p-value histogram plot
#'
#' Internal function for plotting p-value histograms of raw- and adjusted p-values
#'
#' @param data Per-contrast DE results to be plotted, as prepared by
#'   \code{\link{prep_plot_model_data}}.
#' @param contrast The contrast being plotted. Used for generating the plot title.
#'
#' @return A ggplot object
#'
#' @examples
#' # No examples yet
static_pval_histogram <- function(data, contrast) {
  output <- with(data, {
  ggplot(data) +
      geom_histogram(aes(x = P.Value), binwidth = 0.025, color = "black", na.rm = T) +
      xlim(c(-0.05,1.05)) +
      theme_bw() +
      xlab("raw P value") +
      theme(panel.border = element_rect(fill = NA, color = "grey30")) +
      ggplot(data) +
      geom_histogram(aes(x = adj.P.Val), binwidth = 0.025, color = "black", na.rm = T) +
      xlim(c(-0.05,1.05)) +
      theme_bw() +
      xlab("adjusted P value") +
      theme(panel.border = element_rect(fill = NA, color = "grey30")) +
      patchwork::plot_annotation(
        title = stringr::str_replace(contrast, "_vs_", " vs ")
        )
  })
 output
}







