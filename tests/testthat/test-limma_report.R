test_that("prep_plot_model_data outputs the proper data and format", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  # data frame with expected colnames,
  # same number of rows as data,
  # and no missing data in the sig cols

  out_control <- prep_plot_model_data(model_results = input$results,
                                      contrast = "control")

  out_control <- prep_plot_model_data(model_results = input$results,
                                      contrast = "treatment")

  expect_cols <- c("logFC", "CI.L", "CI.R", "AveExpr", "t", "B", "P.Value",
                   "adj.P.Val", "sig.PVal", "sig.FDR", "P value",
                   "Adjusted P value", "negLog10rawP", "negLog10adjP",
                   "sig.pval.fct", "sig.FDR.fct")



  for (contrast in c("control", "treatment")) {
    output <- prep_plot_model_data(model_results = input$results,
                                   contrast = contrast)

    expect_equal(colnames(output), expect_cols)
    expect_equal(nrow(output), nrow(input$data))
    expect_equal(sum(is.na(output$sig.PVal)), 0)
    expect_equal(sum(is.na(output$FDR)), 0)
  }
})

# don't test ggplot functions



# tests for write_limma_plots should be similar to
# limma tables tests

test_that("write_limma_plots gives error when there are no results", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))


  expect_error(write_limma_plots(input, overwrite = T), "does not have results")
})

test_that("write_limma_plots checks grouping column", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  # must specify grouping col
  expect_error(write_limma_plots(input), "cannot be empty")

  # can't specify more than one
  expect_error(write_limma_plots(input, grouping_column = c("group", "batch")),
               "Only specify one")

  # must be present in the data
  expect_error(write_limma_plots(input, grouping_column = "xxx"),
               "not found")

})

test_that("write_limma_plots checks table columns", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      table_columns = c("xxx", "yyy")),
    "not found in annotation"
  )
})

test_that("write_limma_plots checks title column", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  # Must be of length 1
  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      title_column = c("xxx", "yyy")),
    "does not equal 1"
  )


  # must be present
  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      title_column = c("xxx")),
    "not found in annotation"
  )

  # values must be unique
  input$annotation$uniprot_id[2] <- input$annotation$uniprot_id[1]
  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment"),
    "not unique"
  )
})

test_that("write_limma_plots check height and width args", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      height = "xxx"),
    "numeric value greater"
  )

  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      height = -1),
    "numeric value greater"
  )

  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      width = "xxx"),
    "numeric value greater"
  )

  expect_error(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      width = -1),
    "numeric value greater"
  )
})


test_that("write_limma_plots informs user of output dir", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("control_DA_report.html",
                   "treatment_DA_report.html",
                   "static_plots"), recursive = T), add = T)

  suppressMessages(
    expect_message(
      write_limma_plots(
        input,
        grouping_column = "treatment",
        overwrite = T),
      "current working directory:"
    )
  )
})


test_that("write_limma_plots won't overwrite when false", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("control_DA_report.html",
                   "treatment_DA_report.html",
                   "static_plots"), recursive = T), add = T)


  suppressMessages(write_limma_plots(input, grouping_column = "treatment"))
  suppressMessages(
    expect_error(
      write_limma_plots(input,
                        grouping_column = "treatment",
                        overwrite = F),
      "Change \\`output_dir\\`"
    )
  )
})

test_that("write_limma_plots messages user when overwriting", {
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("control_DA_report.html",
                   "treatment_DA_report.html",
                   "static_plots"), recursive = T), add = T)


  suppressMessages(write_limma_plots(input, grouping_column = "treatment"))
  suppressMessages(
    expect_message(
      write_limma_plots(input,
                        grouping_column = "treatment",
                        overwrite = T),
      "Overwriting"
    )
  )
})

test_that("write_limma_plots does not alter user working directory", {
  in_wd <- getwd()
  input <- readRDS(test_path("fixtures", "final_output_input.rds"))

  on.exit(unlink(c("control_DA_report.html",
                   "treatment_DA_report.html",
                   "static_plots"), recursive = T), add = T)


  suppressMessages(write_limma_plots(input, grouping_column = "treatment"))
  expect_equal(getwd(), in_wd)


  # When specifying custom wd
  on.exit(unlink(c("tempfortesting/control_DA_report.html",
                   "tempfortesting/treatment_DA_report.html",
                   "tempfortesting/static_plots",
                   "tempfortesting/"), recursive = T), add = T)

  suppressMessages(
    expect_message(
      write_limma_plots(input,
                        output_dir = "tempfortesting",
                        grouping_column = "treatment"),
      "Setting working directory to output directory"
    )
  )
  expect_equal(getwd(), in_wd)
})


test_that("write_limma_plots creates output files", {

  input <- readRDS(test_path("fixtures", "final_output_input.rds"))


  expected_files_default <- c("control_DA_report.html",
                              "treatment_DA_report.html",
                              "static_plots/control-MD-adjusted-pval.pdf",
                              "static_plots/control-MD-raw-pval.pdf",
                              "static_plots/control-pval-hist.pdf",
                              "static_plots/control-volcano-adjusted-pval.pdf",
                              "static_plots/control-volcano-raw-pval.pdf",
                              "static_plots/treatment-MD-adjusted-pval.pdf",
                              "static_plots/treatment-MD-raw-pval.pdf",
                              "static_plots/treatment-pval-hist.pdf",
                              "static_plots/treatment-volcano-adjusted-pval.pdf",
                              "static_plots/treatment-volcano-raw-pval.pdf")

  on.exit(unlink(c(expected_files_default,
                   "static_plots"), recursive = T), add = T)

  suppressMessages(write_limma_plots(input,
                                     grouping_column = "treatment",
                                     overwrite = F))
  expect_true(all(file.exists(expected_files_default)))

  # Non-standard directory and files
  expected_files_custom <- c("tempfortest/control_DA_report.html",
                             "tempfortest/treatment_DA_report.html",
                             "tempfortest/static_plots/control-MD-adjusted-pval.pdf",
                             "tempfortest/static_plots/control-MD-raw-pval.pdf",
                             "tempfortest/static_plots/control-pval-hist.pdf",
                             "tempfortest/static_plots/control-volcano-adjusted-pval.pdf",
                             "tempfortest/static_plots/control-volcano-raw-pval.pdf",
                             "tempfortest/static_plots/treatment-MD-adjusted-pval.pdf",
                             "tempfortest/static_plots/treatment-MD-raw-pval.pdf",
                             "tempfortest/static_plots/treatment-pval-hist.pdf",
                             "tempfortest/static_plots/treatment-volcano-adjusted-pval.pdf",
                             "tempfortest/static_plots/treatment-volcano-raw-pval.pdf")
  on.exit(unlink(c(expected_files_custom,
                 "tempfortest/static_plots",
                 "tempfortest"), recursive = T), add = T)
  suppressMessages(
    write_limma_plots(input,
                      grouping_column = "treatment",
                      output_dir = "tempfortest",
                      overwrite = F)
  )
  expect_true(all(file.exists(expected_files_custom)))
})

