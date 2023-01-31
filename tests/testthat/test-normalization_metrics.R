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
  expect_length(COR(data, fourtwo_groups), 7)
  expect_length(COR(data, six_groups), 15)
})

test_that("COR returns expected values", {
  # handle NAs, single groups, and no variance correctly

})
