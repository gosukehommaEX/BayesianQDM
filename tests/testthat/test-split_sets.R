test_that("split_sets() separates the Go and NoGo probabilities", {
  p <- cbind(1:4, 5:8)
  expect_identical(split_sets(p, 1L), list(go = 1:4, nogo = 5:8))
  expect_identical(split_sets(p, 1L, dim = c(2L, 2L))$go, matrix(1:4, 2L))
  p4 <- cbind(1:2, 3:4, 5:6, 7:8)
  expect_identical(split_sets(p4, 2L),
                   list(go = list(1:2, 3:4), nogo = list(5:6, 7:8)))
})
