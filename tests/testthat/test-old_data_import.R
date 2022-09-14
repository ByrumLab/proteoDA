test_that("DIA data extract fxns check existence and type of input file", {
  expect_error(
    read_maxquant_delim(
      input_file = "XXX.xml"
    ),
    "must end"
  )
  expect_error(
    read_maxquant_delim(
      input_file = "XXX.txt"
    ),
    "does not exist"
  )
})

test_that("DIA data extract fxns checks reqd cols", {
  expect_error(
    read_DIA_data(
      input_file = "test_data/DIA_bad-cols.csv"
    ),
    "not present in"
  )
})

test_that("DIA data extract fxns check for sample IDs when provided", {
  expect_error(
    suppressMessages(read_DIA_data(
      input_file = "test_data/DIA_good.csv",
      sample_IDs = c("not-there")
    )),
    "not found in column names"
  )
  expect_error(
    suppressMessages(read_DIA_data(
      input_file = "test_data/DIA_good.csv",
      sample_IDs = c("not-there", "Kinter_120720_DIA_Q15_3D.mzML")
    )),
    "not-there"
  )
})
