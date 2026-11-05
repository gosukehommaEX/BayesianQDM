test_that("code_num() writes numbers with the fewest digits", {
  expect_identical(code_num(0.1), "0.1")
  expect_identical(code_num(c(0.25, 12, -10)), "c(0.25, 12, -10)")
  expect_identical(code_num(c(a = 1)), "1")
  expect_identical(code_num(1e6), "1000000")
  expect_identical(code_num(1e-20), "1e-20")
  expect_identical(code_num(0.1 + 0.2), "0.30000000000000004")
  # Integers keep their type unless the argument is converted anyway
  expect_identical(code_num(12L), "12L")
  expect_identical(code_num(c(3L, 4L), keep_type = FALSE), "c(3, 4)")
})

test_that("code_num() gives back the numbers exactly", {
  x <- c(1 / 3, 0.1 + 0.2, 1e-20, 123456.789, 1e6, -2.5e-7, pi, 1e300)
  for (v in x) {
    expect_identical(as.numeric(code_num(v)), v, label = format(v))
  }
  expect_identical(eval(parse(text = code_num(x))), x)
  # A changed number does not give the original
  expect_false(identical(as.numeric(code_num(x[1L])), x[1L] + 1e-15))
})
