# One file to test util functions

# remove_commas()
test_that("remove_commas works as expected", {
  expect_equal(remove_commas("1,000,000"),
               1000000)
  expect_equal(remove_commas("1000000"),
               1000000)
})

# %notin% operator
test_that("%notin% works as expected", {
  expect_true("a" %notin% c("b", "c"))
  expect_false("a" %notin% c("a", "b"))
})


# make_factor()
test_that("make_factor works as expected", {
  expect_equal(make_factor(c("A", "B", "C")),
               factor(c("A", "B", "C")))
  expect_equal(make_factor(c(1, 2, 3)),
               factor(c("X1", "X2", "X3")))
  expect_equal(make_factor(c(1, 2, 3), prefix = "Z"),
               factor(c("Z1", "Z2", "Z3")))
})


# all_pw_diffs()
test_that("all_pw_diffs works as expected", {
  expect_equal(all_pw_diffs(c(1,1,1)),
               c(0,0,0))
  expect_equal(all_pw_diffs(c(1,2,3)),
               c(1,2,1))
  expect_equal(all_pw_diffs(c(1,2,3,4)),
               c(1,2,3,1,2,1))
})

# colorGroup()
# NEED TO ADD

# validate_filename()
test_that("validate_filename works as expected", {
  expect_true(validate_filename("xx.csv", c("csv", "tsv")))
  expect_error(validate_filename("xx .csv", c("csv", "tsv")),
               "spaces")
  expect_true(validate_filename("xx .csv", c("csv", "tsv"), check_space = F))
  expect_error(validate_filename("xxx.txt", c("csv", "tsv")),
               "Permitted")
})


# file_extension()
test_that("file_extension works as expected", {
  expect_equal(file_extension("test/file/path.txt"),
               "txt")
  expect_equal(file_extension("test/file/path.CSV"),
               "CSV")
  expect_equal(file_extension("test/file/path"),
               "")
})



# rowVars
test_that("rowVars works as expected", {
  # Make a test matrix
  test_matrix <- matrix(data = c(0,1,3,4,
                                 1,1,1,1,
                                 1,2,3,NA,
                                 4,3,2,1,
                                 4,3,2,1),
                        nrow = 5,
                        ncol = 4,
                        byrow = T)

  expect_equal(rowVars(test_matrix),
               c(3 + 1/3, 0, NA, 1 + 2/3, 1 + 2/3))
  expect_equal(rowVars(test_matrix, na.rm = T),
               c(3 + 1/3, 0, 1, 1 + 2/3, 1 + 2/3))
})


# rowSds
# should just be the sqrt of the row variances
test_that("rowSds works as expected", {
  # Make a test matrix
  test_matrix <- matrix(data = c(0,1,3,4,
                                 1,1,1,1,
                                 1,2,3,NA,
                                 4,3,2,1,
                                 4,3,2,1),
                        nrow = 5,
                        ncol = 4,
                        byrow = T)

  expect_equal(rowSds(test_matrix),
               sqrt(rowVars(test_matrix)))
  expect_equal(rowSds(test_matrix, na.rm = T),
               sqrt(rowVars(test_matrix, na.rm = T)))
})


# rowMedians
test_that("rowMedians works as expected", {
  # Make a test matrix
  test_matrix <- matrix(data = c(0,1,3,4,
                                 1,1,1,1,
                                 1,2,3,NA,
                                 4,3,2,1,
                                 4,3,2,1),
                        nrow = 5,
                        ncol = 4,
                        byrow = T)
  expect_equal(rowMedians(test_matrix),
               c(2,1,NA,2.5,2.5))
  expect_equal(rowMedians(test_matrix, na.rm = T),
               c(2,1,2,2.5,2.5))

  test_matrix2 <- matrix(data = c(0,1,3,4,5,
                                 1,1,1,1,5,
                                 1,2,3,NA,5),
                        nrow = 3,
                        ncol = 5,
                        byrow = T)

  expect_equal(rowMedians(test_matrix2),
               c(3,1,NA))
  expect_equal(rowMedians(test_matrix2, na.rm = T),
               c(3,1,2.5))
})


# rowMads
test_that("rowMads works as expected", {
  # Make a test matrix
  test_matrix <- matrix(data = c(0,1,3,4,
                                 1,1,1,1,
                                 1,2,3,NA,
                                 4,3,2,1,
                                 4,3,2,1),
                        nrow = 5,
                        ncol = 4,
                        byrow = T)
  # Here, with small n, results are a little odd
  # and dependent on the constant used in stats::mad()
  # to ensure consistency
  expect_equal(rowMads(test_matrix),
               c(1.5 * 1.4826, 0, NA, 1 * 1.4826, 1 * 1.4826))
  expect_equal(rowMads(test_matrix, na.rm = T),
               c(1.5 * 1.4826,  0, 1 * 1.4826, 1 * 1.4826, 1 * 1.4826))
})

# check_rows_in
test_that("check_rows_in works as expected", {
  df1 <- data.frame(1:4, row.names = c("a", "b", "c", "d"))
  df2 <- data.frame(1:4, row.names = c("a", "b", "c", "e"))
  ref <- c("a", "b", "c", "d")
  expect_true(check_rows_in(list(df1), ref))
  expect_false(check_rows_in(list(df2), ref))
})
