

# Tests for add_design ----------------------------------------------------

# Most formula things get tested by validate_formula, testing below
# For add design itself test:

# error message when term is not present
test_that("add_design gives helpful error when terms are not present in metadata", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  expect_error(add_design(input, "~ missing"), "not found in metadata")
})


# warning when overwriting an existing design
test_that("add_design gives warning when overwriting an existing statistical design", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  input$design <- list("design")

  expect_message(add_design(input, "~ treatment"), "Overwriting.")

})

# Error when calling an invalid formula
# Mostly tested below
test_that("add_design errors via validate_formula with improper formula", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))

  expect_error(add_design(input, "~ treatment + (1 | sex + group)"), "Multiple effects")
})


# structure of design with non-random factor
test_that("add_design outputs proper design matrix and formulas for models", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  out_treat_intercept <- readRDS(test_path("fixtures", "add_design_output_treat_intercept.rds"))
  out_treat_nointercept <- readRDS(test_path("fixtures", "add_design_output_treat_nointercept.rds"))
  out_int_intercept <- readRDS(test_path("fixtures", "add_design_output_interact_intercept.rds"))
  out_int_nointercept <- readRDS(test_path("fixtures", "add_design_output_interact_nointercept.rds"))




  expect_equal(add_design(input, "~ treatment"), out_treat_intercept)
  expect_equal(add_design(input, "~ 0 + treatment"), out_treat_nointercept)
  expect_equal(add_design(input, "~ treatment * sex"), out_int_intercept)
  expect_equal(add_design(input, "~0 + treatment * sex"), out_int_nointercept)


})

# structure of design with random factor,
# including random term
test_that("add_design outputs proper design matrix and formulas for mixed models", {
  input <- readRDS(test_path("fixtures", "add_design_input.rds"))
  out_treat_mixed <- readRDS(test_path("fixtures", "add_design_output_treat_mixed.rds"))

  expect_equal(add_design(input, "~ treatment + (1 | group)"), out_treat_mixed)

})



# Tests for validate_formula ----------------------------------------------

test_that("validate_formula rejects eqns with LHS", {
  expect_error(
    validate_formula("y ~ 0 + group + xyz + (1 | hello)"),
    "right-hand side"
  )
  expect_error(
    validate_formula(y ~ 0 + group + xyz + (1 | hello)),
    "include a response variable"
  )
})

test_that("validate_formula rejects eqns without tilde", {
  expect_error(
    validate_formula(0 + group + xyz + (1 | hello)),
    "~"
  )
  expect_error(
    validate_formula("0 + group + xyz + (1 | hello)"),
    "~"
  )
})

test_that("validate_formula rejects malformed random factors", {
  expect_error(validate_formula("~ 0 + group + xyz + (1 | )"),
               "Could not")
  expect_error(validate_formula("~ 0 + group + xyz + (| )"),
               "Could not")
  expect_error(validate_formula("~ 0 + group + xyz + (1 || hello)"),
               "||")
  expect_error(validate_formula(~ 0 + group + xyz + (1 || hello)),
               "||")
  expect_error(validate_formula("~ 0 + group + xyz + 1 | hello"),
               "parentheses")
  expect_error(validate_formula(~ 0 + group + xyz + 1 | hello),
               "parentheses")
})

test_that("validate_formula rejects multiple random factors", {
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello) + (1 | hello)"),
               "multiple")
  expect_error(validate_formula("~ 0 + group + xyz + 1 | hello + 1 | hello"),
               "multiple")
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello) + 1 | hello"),
               "multiple")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello) + (1 | hello)),
               "multiple")
  expect_error(validate_formula(~ 0 + group + xyz + 1 | hello + 1 | hello),
               "multiple")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello ) + 1 | hello),
               "multiple")
})

test_that("validate_formula rejects complex or non-intercept random factors", {
  expect_error(validate_formula("~ 0 + group + xyz + (hello | hello)"),
               "only influence the intercept")
  expect_error(validate_formula(~ 0 + group + xyz + (hello | hello)),
               "only influence the intercept")
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello:group)"),
               "Multiple effects")
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello*group)"),
               "Multiple effects")
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello/group)"),
               "Multiple effects")
  expect_error(validate_formula("~ 0 + group + xyz + (1 | hello + group)"),
               "Multiple effects")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello:group)),
               "Multiple effects")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello*group)),
               "Multiple effects")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello/group)),
               "Multiple effects")
  expect_error(validate_formula(~ 0 + group + xyz + (1 | hello + group)),
               "Multiple effects")
})

test_that("validate_formula rejects random factors that are also fixed factors", {
  expect_error(validate_formula("~ group + (1 | group )"),
               "fixed")
  expect_error(validate_formula(~ group + (1 | group   )),
               "fixed")
})

test_that("validate_formula accepts valid eqns", {
  expect_equal(class(validate_formula("~ 0 + group + xyz + (1 | hello)")),
               "formula")
  expect_equal(class(validate_formula("~ group + xyz + (1 | hello)")),
               "formula")
  expect_equal(class(validate_formula("~ 0 + group*xyz + (1 | hello)")),
               "formula")
  expect_equal(class(validate_formula("~ group*xyz")),
               "formula")
  expect_equal(class(validate_formula("~ group*xyz*hello")),
               "formula")
  expect_equal(class(validate_formula("~ group*xyz*hello - xyz:hello")),
               "formula")
  expect_equal(class(validate_formula(~ 0 + group + xyz + (1 | hello))),
               "formula")
  expect_equal(class(validate_formula(~ group + xyz + (1 | hello))),
               "formula")
  expect_equal(class(validate_formula(~ 0 + group*xyz + (1 | hello))),
               "formula")
  expect_equal(class(validate_formula(~ group*xyz)),
               "formula")
  expect_equal(class(validate_formula(~ group*xyz*hello)),
               "formula")
  expect_equal(class(validate_formula(~ group*xyz*hello - xyz:hello)),
               "formula")
})
