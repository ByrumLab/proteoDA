test_that("message prints", {
  expect_message(useless_function("hello"), "function")
})
