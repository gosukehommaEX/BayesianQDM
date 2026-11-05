test_that("code_node() stores a call", {
  expect_identical(code_node("f", list(a = "1")),
                   list(fun = "f", args = list(a = "1")))
  expect_identical(code_node("g"), list(fun = "g", args = list()))
})
