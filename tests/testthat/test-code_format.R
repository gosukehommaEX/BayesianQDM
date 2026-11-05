test_that("code_format() writes a call that fits on one line", {
  node <- code_node("f", list(a = "1", b = code_node("g", list("x"))))
  expect_identical(code_format(node), "f(a = 1, b = g(x))")
  expect_identical(code_format("abc", indent = 100L), "abc")
})

test_that("code_format() places the arguments on several lines", {
  # Expected lines: Python_reference/port_code_format.out
  long <- code_node("fun", list(alpha = strrep("1", 30),
                                beta = strrep("2", 30),
                                gamma = strrep("3", 30)))
  expect_identical(code_format(long), c(
    paste0("fun(alpha = ", strrep("1", 30), ","),
    paste0("    beta = ", strrep("2", 30), ","),
    paste0("    gamma = ", strrep("3", 30), ")")
  ))
  nested <- code_node("outer", list(code_node("inner", list(
    a = strrep("1", 40), b = strrep("2", 40)
  ))))
  expect_identical(code_format(nested), c(
    paste0("outer(inner(a = ", strrep("1", 40), ","),
    paste0("            b = ", strrep("2", 40), "))")
  ))
  # The indentation counts in the width
  expect_length(code_format(code_node("f", list(a = strrep("1", 70))), 0L),
                1L)
  expect_identical(code_format(code_node("f", list(a = strrep("1", 70),
                                                   b = "2")), 10L),
                   c(paste0("f(a = ", strrep("1", 70), ","),
                     paste0(strrep(" ", 12), "b = 2)")))
})

test_that("code_format() indents the arguments of a hanging call", {
  long <- code_node("fun", list(alpha = strrep("1", 30),
                                beta = strrep("2", 30),
                                gamma = strrep("3", 30)))
  expect_identical(code_format(long, hanging = TRUE), c(
    "fun(",
    paste0("  alpha = ", strrep("1", 30), ", beta = ", strrep("2", 30), ","),
    paste0("  gamma = ", strrep("3", 30)),
    ")"
  ))
})
