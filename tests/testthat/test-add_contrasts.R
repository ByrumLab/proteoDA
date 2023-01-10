
# Tests for add_contrasts -------------------------------------------------
# Main function is tested here, with some subfunctions tested below



# Tests for validate_contrasts --------------------------------------------
test_that("validate contrasts checks for the right number of equals signs", {
  good_input <- c("Treatment1_vs_Control= Treatment1 - Control",
                  "Treatment2_vs_Control= Treatment2-Control",
                  "Treatment2 vs Treatment1=Treatment2-Treatment1")
  no_equals <- c("Treatment1_vs_Control Treatment1 - Control",
                 "Treatment2_vs_Control = Treatment2-Control",
                 "Treatment2 vs Treatment1=Treatment2-Treatment1")
  multiple_equals <- c("Treatment1_vs_Control == Treatment1 - Control",
                       "Treatment2_vs_Control == Treatment2-Control",
                       "Treatment2 vs Treatment1=Treatment2-Treatment1")

  expect_equal(validate_contrasts(good_input), good_input)
  expect_error(validate_contrasts(no_equals), "1")
  expect_error(validate_contrasts(multiple_equals), "1 and 2")

})

# Tests for extract_contrast_groups ---------------------------------------

test_that("extract contrast groups extracts group names as expected", {

  # Right now, allow group names to be numbers, but only allowed to add and subtract them
  input <- c("Treatment1_vs_Control= Treatment1 - Control",
             "Treatment2_vs_Control= Treatment2-Control",
             "Treatment2 vs Treatment1=Treatment2-Treatment1")
  output <- c("Treatment1", "Control", "Treatment2")

  input2 <- c("something = 999 - 100",
              "something = XXX - YYY",
              "something else = (XXX-YYY)/2",
              "another one = (a+b+c+d+e+f+g+h+i+j)/10")

  output2 <- c("999", "100", "XXX", "YYY", letters[1:10])

  expect_equal(unique(unlist(extract_contrast_groups(input))), output)
  expect_equal(unique(unlist(extract_contrast_groups(input2))), output2)

})


