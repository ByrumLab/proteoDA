make_limma_reports <- function(model_results = NULL,
                               annotation = NULL,
                               groups = NULL,
                               output_dir = NULL,
                               tmp_dir = "tmp",
                               height = 1000,
                               width = 1000,
                               overwrite = F) {

  if (is.null(output_dir)) {
    output_dir <- "test_output"
  }

  # issue: can't do non-self-contained .html files.
  # This is partially because of issues that are baked into pandoc.
  # there is possibly a fix on the way: https://github.com/rstudio/rmarkdown/pull/2199
  # But, for now, we might be able to get around this differently, and possibly get around some
  # of the working directory issues.

  # if (!dir.exists(output_dir)) {
  #   dir.create(output_dir)
  # }
  # file.copy(from = system.file("report_templates/glimma_xy_plot.Rmd",
  #                              package = "proteomicsDIA"),
  #           to = file.path(output_dir, ))



  # system.file("report_templates/limma_report_per_contrast.Rmd",
  #             package = "proteomicsDIA", mustWork = T)

  # Could try:
  # create output directory structure
  # copy templates and resources in.
  # set wd to output
  # do all knitting
  # remove intermediates and templates
  # return to old wd on exit




  resources <- "resources"
  contrast_count <- 1
  num_contrasts <- length(names(model_results$stats_by_contrast))
  cli::cli_rule()

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
    # Need to suppress warnings, otherwise warned about empty title (which we want, so the logo can be first)
    # and changing workding directories during kniting (which works fine, but I can't figure out how to stop
    # the warning)
    suppressWarnings(rmarkdown::render(system.file("report_templates/limma_report_per_contrast.Rmd",
                                  package = "proteomicsDIA", mustWork = T),
                      knit_root_dir = getwd(),
                      output_dir = output_dir,
                      intermediates_dir = tmp_dir,
                      output_file = paste0(contrast, "_DE_report.html"),
                      quiet = T))
    contrast_count <- contrast_count + 1
  }


  # remove temporary directory.

  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))

  invisible(NULL)
}

