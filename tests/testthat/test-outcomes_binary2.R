test_that("outcomes_binary2() lists all count vectors once", {
  for (n in c(0, 1, 4, 10)) {
    x <- outcomes_binary2(n)
    expect_identical(nrow(x), as.integer(choose(n + 3, 3)))
    expect_true(all(rowSums(x) == n))
    expect_true(all(x >= 0))
    expect_identical(anyDuplicated(x), 0L)
  }
  expect_identical(colnames(outcomes_binary2(2)), c("00", "01", "10", "11"))
  expect_type(outcomes_binary2(2), "integer")
})
