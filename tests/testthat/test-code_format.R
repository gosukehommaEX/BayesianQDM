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
  # The indentation counts in the width: the aligned arguments of a call at
  # column 10 do not fit, so they start on the next line
  expect_length(code_format(code_node("f", list(a = strrep("1", 70))), 0L),
                1L)
  expect_identical(code_format(code_node("f", list(a = strrep("1", 70),
                                                   b = "2")), 10L),
                   c("f(", paste0("  a = ", strrep("1", 70), ","), "  b = 2",
                     ")"))
})

test_that("code_format() counts the characters written after the call", {
  # Expected lines: Python_reference/port_code_format.out
  one <- code_node("f", list(a = strrep("1", 70)))
  expect_identical(code_format(one, tail = 3L), code_flat(one))
  expect_identical(code_format(one, tail = 4L),
                   c("f(", paste0("  a = ", strrep("1", 70)), ")"))
  # Without the fallback, the aligned arguments are kept
  expect_identical(code_format(one, tail = 4L, fallback = FALSE),
                   code_flat(one))
})

test_that("code_format() starts nested arguments on the next line", {
  # Expected lines: Python_reference/port_code_format.out (B7). The
  # arguments start two spaces to the right of the line of the call (base)
  hyp <- code_node("hypothetical", list(data = code_node(
    "data_binary2", list(counts = "c(2, 1, 1, 0)")
  )))
  node <- code_node("qdm_control", list(hypothetical_c = hyp))
  expect_identical(code_format(node, 2L, width = 77L, base = 2L), c(
    "qdm_control(",
    "    hypothetical_c = hypothetical(",
    "      data = data_binary2(counts = c(2, 1, 1, 0))",
    "    )",
    "  )"
  ))
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

test_that("code_format() starts the argument after a long one on a new line", {
  # Expected lines: Python_reference/port_code_format.out
  inner <- code_node("g", list(x = strrep("1", 40), y = strrep("2", 40)))
  node <- code_node("f", list(a = inner, b = "3"))
  expect_identical(code_format(node), c(
    paste0("f(a = g(x = ", strrep("1", 40), ","),
    paste0("        y = ", strrep("2", 40), "),"),
    "  b = 3)"
  ))
})
