test_that("dbetabinom() equals the direct formula for small m", {
  k <- 0:10
  direct <- choose(10, k) * beta(k + 2.5, 10 - k + 1.5) / beta(2.5, 1.5)
  expect_equal(dbetabinom(k, 10, 2.5, 1.5), direct, tolerance = 1e-12)
})

test_that("dbetabinom() sums to 1 and stays finite for large m", {
  p <- dbetabinom(0:5000, 5000, 50.5, 60.5)
  expect_true(all(is.finite(p)))
  expect_equal(sum(p), 1, tolerance = 1e-10)
})
