make_limma_reports <- function(model_results = NULL,
                               annotation = NULL,
                               groups = NULL,
                               output_dir = NULL,
                               tmp_dir = "tmp",
                               height = 1000,
                               width = 1000,
                               overwrite = F) {

  # TODO: add checks

  if (is.null(output_dir)) {
    output_dir <- "test_output"
  }

  # issue: can't do non-self-contained .html files.
  # This is partially because of issues that are baked into pandoc.
  # there is possibly a fix on the way: https://github.com/rstudio/rmarkdown/pull/2199
  # But, for now, we might be able to get around this differently, and possibly get around some
  # of the working directory issues.

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
  for (contrast in names(model_results$stats_by_contrast)) {
    data <- model_results$stats_by_contrast[[contrast]]
    data$`P value` <- data$P.Value
    data$`Adjusted P value` <- data$adj.P.Val
    data$negLog10rawP <- -log(data$P.Value, 10)
    data$negLog10adjP <- -log(data$adj.P.Val, 10)
    data$sig.PVal <- ifelse(is.na(data$sig.PVal), 0, data$sig.PVal)
    data$sig.FDR <- ifelse(is.na(data$sig.FDR), 0, data$sig.FDR)

    counts <- model_results$data[rownames(data),]
    anno <- annotation[rownames(data), ]

    cli::cli_inform("Writing report for contrast {contrast_count} of {num_contrasts}: {.val {contrast}}")
    rmarkdown::render("report_template.Rmd",
                      knit_root_dir = getwd(),
                      intermediates_dir = tmp_dir,
                      output_file = paste0(contrast, "_DE_report.html"),
                      quiet = T)
    contrast_count <- contrast_count + 1
  }

  # Clean up and exit

  # TODO: add checks that files exist?

  cli::cli_inform("Removing temporary files from {.path {output_dir}}")
  unlink(c("logo_higherres.png", "plot_template.Rmd", "report_template.Rmd", tmp_dir), recursive = T)
  cli::cli_inform("Returning working directory to {.path {old_wd}}")
  setwd(old_wd)
  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(NULL)
}

