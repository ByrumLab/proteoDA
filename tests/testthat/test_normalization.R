# For normalize data function, test:
# sensible error when asking for invalid normalization option
test_that("normalize_data gives good error with incorrect method", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  expect_error(normalize_data(input, norm_method = "xxx"), "norm_method")
})

# error when data are already normalized
test_that("normalize_data gives error when input data are already normalized", {

  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  normal <- normalize_data(input, norm_method = "log2")

  expect_error(normalize_data(normal, norm_method = "cycloess"), "already")
})



# proper output structure
test_that("output of normalize_data is correct", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  normal <- normalize_data(input, norm_method = "log2")

  expect_true(normal$tags$normalized)
  expect_equal(normal$tags$norm_method, "log2")
  expect_equal(normal$data, log2Norm(input$data))


  normal2 <- normalize_data(input, norm_method = "mean")
  expect_true(normal2$tags$normalized)
  expect_equal(normal2$tags$norm_method, "mean")
  expect_equal(normal2$data, meanNorm(log2Norm(input$data)))

})


# Individual normalization methods ----------------------------------------
# log2Norm
test_that("log2Norm works", {
  raw <- as.matrix(data.frame(V1 = c(10, 20, 30),
                              V2 = c(10, 20, NA),
                              V3 = c(30, 15, 5),
                              V4 = c(0, 10, NA)))
  normed <- as.matrix(data.frame(V1 = log2(c(10, 20, 30)),
                               V2 = log2(c(10, 20, NA)),
                               V3 = log2(c(30, 15, 5)),
                               V4 = log2(c(NA, 10, NA))))

  expect_equal(log2Norm(raw), normed)
})

# medianNorm
test_that("medianNorm works", {
  raw <- as.matrix(data.frame(V1 = c(10, 20, 30),
                              V2 = c(10, 20, NA),
                              V3 = c(30, 15, 5),
                              V4 = c(0, 10, NA)))

  int <- as.matrix(data.frame(V1 = c(10, 20, 30)/20,
                              V2 = c(10, 20, NA)/15,
                              V3 = c(30, 15, 5)/15,
                              V4 = c(0, 10, NA)/5))
  normed <- int*mean(c(5,15,15,20))

  expect_equal(medianNorm(raw), normed)
})

# meanNorm
test_that("meanNorm works", {
  raw <- as.matrix(data.frame(V1 = c(10, 20, 30),
                              V2 = c(10, 20, NA),
                              V3 = c(30, 15, 5),
                              V4 = c(0, 10, NA)))

  int <- as.matrix(data.frame(V1 = c(10, 20, 30)/20,
                              V2 = c(10, 20, NA)/15,
                              V3 = c(30, 15, 5)/(50/3),
                              V4 = c(0, 10, NA)/5))
  normed <- int*mean(c(5,15,(50/3),20))

  expect_equal(meanNorm(raw), normed)
})

# A little tricky to test our functions which use extrenal normaliation functions
# Will just test that they don't change


# vsnNorm
test_that("vsnNorm hasn't changed", {
  raw <- readRDS(test_path("fixtures", "normalization_external_fxn_input.rds"))


  expect_snapshot_value(vsnNorm(raw), style = "json2")
})

# quantileNorm
test_that("quantileNorm hasn't changed", {
  raw <- readRDS(test_path("fixtures", "normalization_external_fxn_input.rds"))


  expect_snapshot_value(quantileNorm(raw), style = "json2")
})


# cycloessNorm
test_that("cycloessNorm hasn't changed", {
  raw <- readRDS(test_path("fixtures", "normalization_external_fxn_input.rds"))


  expect_snapshot_value(cycloessNorm(raw), style = "json2")
})


# rlrNorm
# Also a little hard to calculate, will make sure it hasn't changed
test_that("rlrNorm hasn't changed", {
  raw <- readRDS(test_path("fixtures", "normalization_external_fxn_input.rds"))

  expect_snapshot_value(rlrNorm(raw), style = "json2")
})


# giNorm
test_that("giNorm works", {
  raw <- as.matrix(data.frame(V1 = c(10, 20, 30),
                              V2 = c(10, 20, NA),
                              V3 = c(30, 15, 5),
                              V4 = c(0, 10, NA)))

  int <- as.matrix(data.frame(V1 = c(10, 20, 30)/60,
                              V2 = c(10, 20, NA)/30,
                              V3 = c(30, 15, 5)/50,
                              V4 = c(0, 10, NA)/10))
  normed <- log2(int*40)
  # Convert -Inf to NA
  normed[1,4] <- NA

  expect_equal(giNorm(raw), normed)
})
