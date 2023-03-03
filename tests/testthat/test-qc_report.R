# No tests for qc_plotting.R, which is all plotting functions that should
# just output plots

test_that("write_qc_report warns for unnormalized data", {

  on.exit(unlink("QC_Report.pdf"), add = T)

  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  suppressMessages(
    expect_warning(
      write_qc_report(input,
                      overwrite = T),
      "not normalized"
    )
  )

  input$tags$normalized <- F

  suppressMessages(
    expect_warning(
      write_qc_report(input,
                      overwrite = T),
      "not normalized"
    )
  )
})


test_that("write_qc_report checks color column", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  on.exit(unlink("QC_Report.pdf"), add = T)

  # message when color column empty
  suppressMessages(
    expect_message(
      write_qc_report(input,
                      overwrite = T),
      "as one group"
    )
  )

  # can't specify more than one
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      color_column = c("group", "batch"),
                      overwrite = T),
      "Only specify one"
    )
  )

  # must be present in the data
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      color_column = "xxx",
                      overwrite = T),
      "not found"
    )
  )
})

test_that("write_qc_report checks label column", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  on.exit(unlink("QC_Report.pdf"), add = T)

  # can't specify more than one
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      label_column = c("group", "batch"),
                      overwrite = T),
      "Only specify one"
    )
  )


  # must be present in the data
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      label_column = "xxx",
                      overwrite = T),
      "not found"
    )
  )
})

test_that("write_qc_report warns when truncating sample labels", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  on.exit(unlink("QC_Report.pdf"), add = T)

  input$metadata$id_long <- stringr::str_pad(input$metadata$sample_ID,
                                             width = 25,
                                             side = "right",
                                             pad = "x")

  suppressMessages(
    expect_message(
      write_qc_report(input,
                      label_column = "id_long",
                      overwrite = T),
      "Truncating sample labels to 20"
    )
  )
})

test_that("write_qc_report gives error if sample labels aren't unique", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  on.exit(unlink("QC_Report.pdf"), add = T)

  input$metadata$id_long <- stringr::str_pad(input$metadata$sample_ID,
                                             width = 25,
                                             side = "left",
                                             pad = "x")
  # With truncation
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      label_column = "id_long",
                      overwrite = T),
      "not unique after truncation"
    )
  )
  # without truncation
  input$metadata$sample_ID[1] <- input$metadata$sample_ID[2]
  suppressMessages(
    expect_error(
      write_qc_report(input,
                      label_column = "sample_ID",
                      overwrite = T),
      "column are not unique"
    )
  )

})


test_that("write_qc_report informs user of output dir and filename when not supplied", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")

  on.exit(unlink("QC_Report.pdf"), add = T)

  # Msg for dir
  suppressMessages(
    expect_message(
      write_qc_report(input,
                      color_column = "group",
                      overwrite = T),
      "current working directory:"
    )
  )

  # mesage for file
  suppressMessages(
    expect_message(
      write_qc_report(input,
                      color_column = "group",
                      overwrite = T),
      "QC_Report.pdf"
    )
  )
})


test_that("write_qc_report validates filename when supplied", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")

  suppressMessages(
    expect_error(
      write_qc_report(input,
                      color_column = "group",
                      filename = "bad.xxx",
                      overwrite = T),
      "file extension"
    )
  )
})


test_that("write_qc_report won't overwrite when false", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  sink <- file.create("temp.pdf")
  sink <- dir.create("tempdirfortesting")
  sink <- file.create("tempdirfortesting/temp.pdf")
  on.exit(unlink("temp.pdf"), add = T)
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)

  suppressMessages(
    expect_error(
      write_qc_report(input,
                      color_column = "group",
                      filename = "temp.pdf",
                      overwrite = F),
      "unique name"
    )
  )

  suppressMessages(
    expect_error(
      write_qc_report(input,
                      color_column = "group",
                      output_dir = "tempdirfortesting",
                      filename = "temp.pdf",
                      overwrite = F),
      "unique name"
    )
  )
})


test_that("write_qc_report messages user when overwriting", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")

  sink <- file.create("temp.pdf")
  sink <- dir.create("tempdirfortesting")
  sink <- file.create("tempdirfortesting/temp.pdf")
  on.exit(unlink("temp.pdf"), add = T)
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)


  suppressMessages(
    expect_message(
      write_qc_report(input,
                      color_column = "group",
                      filename = "temp.pdf",
                      overwrite = T),
      "Overwriting."
    )
  )

  suppressMessages(
    expect_message(
      write_qc_report(input,
                      color_column = "group",
                      output_dir = "tempdirfortesting",
                      filename = "temp.pdf",
                      overwrite = T),
      "Overwriting."
    )
  )
})


test_that("write_qc_report creates proper output file", {
  # default filenames
  input <- readRDS(test_path("fixtures", "norm_report_input.rds")) %>%
    normalize_data("log2")
  on.exit(unlink("QC_Report.pdf"), add = T)

  suppressMessages(
    write_qc_report(input,
                    overwrite = T)
  )

  expect_true(file.exists("QC_Report.pdf"))


  # custom filenames
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)

  suppressMessages(
    write_qc_report(input,
                    color_column = "group",
                    output_dir = "tempdirfortesting",
                    filename = "temp.pdf",
                   overwrite = T)
  )

  expect_true(file.exists("tempdirfortesting/temp.pdf"))
})

