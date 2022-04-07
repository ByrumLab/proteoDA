test_that("DIA data extract fxns check arguments", {
  expect_error(
    extract_data(
      file = "XXX.txt",
      pipe = "DIX",
      enrich = "protein"
    ),
    "DIA"
  )
  expect_error(
    extract_data(
      file = "XXX.txt",
      pipe = "DIA",
      enrich = "Protein"
    ),
    "protein"
  )
})

test_that("DIA data extract fxns check existence and type of input file", {
  expect_error(
    import_data(
      file = "XXX.xml",
      pipe = "DIA",
      enrich = "protein"
    ),
    "must end"
  )
  expect_error(
    import_data(
      file = "XXX.txt",
      pipe = "DIA",
      enrich = "protein"
    ),
    "does not exist"
  )
})

test_that("DIA data extract fxns checks reqd cols", {
  expect_error(
    import_data(
      file = "test_data/DIA_bad-cols.csv",
      pipe = "DIA",
      enrich = "protein"
    ),
    "not present in"
  )
})

test_that("DIA data extract fxns check for sample IDs when provided", {
  expect_error(
    suppressMessages(extract_data(
      file = "test_data/DIA_good.csv",
      sampleIDs = c("not-there"),
      pipe = "DIA",
      enrich = "protein"
    )),
    "not found in column names"
  )
  expect_error(
    suppressMessages(extract_data(
      file = "test_data/DIA_good.csv",
      sampleIDs = c("not-there", "Kinter_120720_DIA_Q15_3D.mzML"),
      pipe = "DIA",
      enrich = "protein"
    )),
    "not-there"
  )
  file.remove("DIA_good_Samples_Report_BQ.csv")
})
