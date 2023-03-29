test_that("write_norm_report gives error for already-normalized data", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  input$tags$normalized <- T

  expect_error(write_norm_report(input), "already normalized")
})


test_that("write_norm_report checks grouping column", {

  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))

  # must specify grouping col
  expect_error(write_norm_report(input), "cannot be empty")

  # can't specify more than one
  expect_error(write_norm_report(input, grouping_column = c("group", "batch")), "Only specify one")

  # must be present in the data
  expect_error(write_norm_report(input, grouping_column = "xxx"), "not found")


  # must be at least two groups in the grouping col
  input$metadata$onegroup <- rep("a", nrow(input$metadata))
  expect_error(write_norm_report(input, grouping_column = "onegroup"), "at least two different groups")
})

test_that("write_norm_report informs user of output dir and filename when not supplied", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))

  on.exit(unlink("normalization_report.pdf"), add = T)
  suppressMessages(expect_message(write_norm_report(input, grouping_column = "group", overwrite = T), "current working directory:"))
  suppressMessages(expect_message(write_norm_report(input, grouping_column = "group", overwrite = T), "normalization_report.pdf"))

})

test_that("write_norm_report validates filename when supplied", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))

  suppressMessages(expect_error(write_norm_report(input,
                                                  grouping_column = "group",
                                                  filename = "norm_report.xxx",
                                                  overwrite = T), "file extension"))
})


test_that("write_norm_report won't overwrite when false", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  sink <- file.create("temp.pdf")
  sink <- dir.create("tempdirfortesting")
  sink <- file.create("tempdirfortesting/temp.pdf")
  on.exit(unlink("temp.pdf"), add = T)
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)

  suppressMessages(expect_error(write_norm_report(input,
                                                  grouping_column = "group",
                                                  filename = "temp.pdf",
                                                  overwrite = F), "unique name"))
  suppressMessages(expect_error(write_norm_report(input,
                                                  grouping_column = "group",
                                                  output_dir = "tempdirfortesting",
                                                  filename = "temp.pdf",
                                                  overwrite = F), "unique name"))

})


test_that("write_norm_report messages user when overwriting", {
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  sink <- file.create("temp.pdf")
  sink <- dir.create("tempdirfortesting")
  sink <- file.create("tempdirfortesting/temp.pdf")
  on.exit(unlink("temp.pdf"), add = T)
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)

  suppressMessages(expect_message(write_norm_report(input,
                                                  grouping_column = "group",
                                                  filename = "temp.pdf",
                                                  overwrite = T), "Overwriting."))
  suppressMessages(expect_message(write_norm_report(input,
                                                  grouping_column = "group",
                                                  output_dir = "tempdirfortesting",
                                                  filename = "temp.pdf",
                                                  overwrite = T), "Overwriting."))
})

test_that("write_norm_report creates proper output file", {
  # default filenames
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  on.exit(unlink("normalization_report.pdf"), add = T)

  suppressMessages(write_norm_report(input, grouping_column = "group", overwrite = T))
  expect_true(file.exists("normalization_report.pdf"))


  # custom filenames
  on.exit(unlink(c("tempdirfortesting/temp.pdf", "tempdirfortesting"), recursive = T), add = T)

  suppressMessages(write_norm_report(input,
                                     grouping_column = "group",
                                     output_dir = "tempdirfortesting",
                                     filename = "temp.pdf",
                                     overwrite = T))
  expect_true(file.exists("tempdirfortesting/temp.pdf"))

})


test_that("write_norm_report creates proper output file when using ggrastr", {
  # default filenames
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  on.exit(unlink("normalization_report.pdf"), add = T)

  skip_if_not_installed("ggrastr")
  suppressMessages(write_norm_report(input, grouping_column = "group", use_ggrastr = T, overwrite = T))
  expect_true(file.exists("normalization_report.pdf"))
})


test_that("write_norm_report creates proper output file when using suppressing zoom legend", {
  # default filenames
  input <- readRDS(test_path("fixtures", "norm_report_input.rds"))
  on.exit(unlink("normalization_report.pdf"), add = T)

  skip_if_not_installed("ggrastr")
  suppressMessages(write_norm_report(input, grouping_column = "group", suppress_zoom_legend = T, overwrite = T))
  expect_true(file.exists("normalization_report.pdf"))
})
