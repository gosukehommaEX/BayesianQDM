test_that("code_flat() writes a call on one line", {
  expect_identical(code_flat("1"), "1")
  expect_identical(code_flat(code_node("f")), "f()")
  node <- code_node("f", list(a = "1", b = code_node("g", list("x",
                                                              y = "2"))))
  expect_identical(code_flat(node), "f(a = 1, b = g(x, y = 2))")
  expect_identical(code_flat(code_node("h", list("x", "y"))), "h(x, y)")
})
