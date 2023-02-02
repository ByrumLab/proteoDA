test_that("PCV, PMAD, and PEV return a named vector of group names", {
  data <- data.frame(sampleA = c(10000, 15000, 30000),
                     sampleB = c(12000, 16000, 31000),
                     sampleC = c(13000, 17000, 32000),
                     sampleD = c(14000, 18000, 33000))
  groups <- c("control", "control", "treatment", "treatment")

  expect_equal(names(PCV(data, groups)), c("control", "treatment"))
  expect_equal(names(PMAD(data, groups)), c("control", "treatment"))
  expect_equal(names(PEV(data, groups)), c("control", "treatment"))
})

test_that("PCV returns expected values", {
  # control should be 0 (no variation)
  # treatment should be mean of 0 and 0.707 = 0.3525534
  # treatment1 should be the mean of 0, 0.942809, and 0.707 =
  data <- data.frame(sampleA = c(10000, 15000, 30000),
                     sampleB = c(10000, 15000, 30000),
                     sampleC = c(10000, 17000, 1),
                     sampleD = c(10000, NA, 3),
                     sampleE = c(10000, 10, 1),
                     sampleF = c(10000, 20, 3))
  groups <- c("control", "control", "treatment", "treatment", "treatment2", "treatment2")

  expect_result <- c(
    0,
    mean(c(0, sd(c(1,3))/2)),
    mean(c(0, sd(c(10,20))/15, sd(c(1,3))/2))
    )

  expect_equal(unname(PCV(data, groups)), expect_result)
})

test_that("PMAD returns expected values", {
  # control should be 0.
  # treatment 1 (C and D) should be based on just the first and third proteins
  # Should be 2/2 = 1
  # Treatment 2 (E and F) should be (0 + 50 + 2)/3 = 17.33

  # Waiting for a response from Charity. PMAD code and documentation are in conflict,
  # not sure which is right.

  # data <- data.frame(sampleA = c(10000, 15000, 30000),
  #                    sampleB = c(10000, 15000, 30000),
  #                    sampleC = c(10000, 17000, 1),
  #                    sampleD = c(10000, NA, 3),
  #                    sampleE = c(10000, 10, 1),
  #                    sampleF = c(10000, 20, 3))
  # groups <- c("control", "control", "treatment", "treatment", "treatment2", "treatment2")
  #
  # expect_equal(unname(PMAD(data, groups)), c(0, 1, 52/3))

})

test_that("PEV returns expected values", {
  # control should be 0 variance
  # treatment 1 (C and D) should be based on just the first and third proteins
  # Should be 2/2 = 1
  # Treatment 2 (E and F) should be (0 + 50 + 2)/3 = 17.33
  data <- data.frame(sampleA = c(10000, 15000, 30000),
                     sampleB = c(10000, 15000, 30000),
                     sampleC = c(10000, 17000, 1),
                     sampleD = c(10000, NA, 3),
                     sampleE = c(10000, 10, 1),
                     sampleF = c(10000, 20, 3))
  groups <- c("control", "control", "treatment", "treatment", "treatment2", "treatment2")

  expect_result <- c(0, 1, 52/3)

  expect_equal(unname(PEV(data, groups)), expect_result)
})

test_that("COR returns results of proper length", {

  data <- data.frame(c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000))
  single_groups <- c("a", "b", "c", "d", "e", "f")
  pair_groups <- c("a", "a", "b", "b", "c", "c")
  threes_groups <- c("a", "a", "a", "b", "b", "b")
  fourtwo_groups <- c("a", "a", "a", "a", "b", "b")
  six_groups <- c("a", "a", "a", "a", "a", "a")

  expect_length(COR(data, single_groups), 6)
  expect_length(COR(data, pair_groups), 3)
  expect_length(COR(data, threes_groups), 6)
  expect_length(COR(data, fourtwo_groups), 7)
  expect_length(COR(data, six_groups), 15)
})

test_that("COR returns expected values", {
  # handle NAs, single groups, and no variance correctly
  data <- data.frame(c(1000, 2000, 3000),
                     c(1000, 2000, 3000),
                     c(1000, 2000, NA),
                     c(2000, 1000, 3000),
                     c(1000, 2000, 3000))
  groups <- c("a", "a", "b", "b", "c")
  expect_result <- c(1, -1, 1)
  expect_equal(unname(COR(data, groups)), expect_result)
})

test_that("log2ratio returns results of proper length", {

  data <- data.frame(c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000),
                     c(10000, 20000))
  single_groups <- c("a", "b", "c", "d", "e", "f")
  pair_groups <- c("a", "a", "b", "b", "c", "c")
  threes_groups <- c("a", "a", "a", "b", "b", "b")
  fourtwo_groups <- c("a", "a", "a", "a", "b", "b")
  six_groups <- c("a", "a", "a", "a", "a", "a")

  # with no protein ID
  # length should be number of pw comparisons times # of proteins
  expect_length(log2ratio(data, single_groups), 30)
  expect_length(log2ratio(data, pair_groups), 6)
  expect_length(log2ratio(data, threes_groups), 2)
  expect_length(log2ratio(data, fourtwo_groups), 2)
  expect_length(log2ratio(data, six_groups), 0)


  # With protein ID
  # dims should be number of proteins (2)
  # and number of pw comparisons
  expect_equal(dim(log2ratio(data, single_groups, keep_protein_ID = T)), c(2, 15))
  expect_equal(dim(log2ratio(data, pair_groups, keep_protein_ID = T)), c(2, 3))
  expect_equal(dim(log2ratio(data, threes_groups, keep_protein_ID = T)), c(2, 1))
  expect_equal(dim(log2ratio(data, fourtwo_groups, keep_protein_ID = T)), c(2, 1))
  expect_equal(dim(log2ratio(data, six_groups, keep_protein_ID = T)), c(2, 0))
})


test_that("log2ratio calculates correct values", {

  # handle NAs, single groups, and no variance correctly
  data <- data.frame(c(500, 1500, 500),
                     c(500, 1500, 500),
                     c(2000, 2000, NA),
                     c(2000, 2000, 2000),
                     c(1500, 1500, 1500))
  groups <- c("a", "a", "b", "b", "c")
  expect_result <- c(1500, 500, 1500, 1000, 0, 1000, -500, -500, -500)
  expect_result_df <- data.frame(V1 = c(1500, 500, 1500),
                                 V2 = c(1000, 0, 1000),
                                 V3 = c(-500, -500, -500),
                                 row.names = c("X1", "X2", "X3"))
  expect_equal(log2ratio(data, groups), expect_result)
  expect_equal(log2ratio(data, groups, keep_protein_ID = T), expect_result_df)
})
