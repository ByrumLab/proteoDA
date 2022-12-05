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
