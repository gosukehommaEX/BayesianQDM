test_that("posterior_binary() adds the data to the prior", {
  pr <- prior_beta(0.5, 1.5)
  post <- posterior_binary(pr, 10, 0:10)
  expect_equal(post$a, 0.5 + 0:10)
  expect_equal(post$b, 1.5 + 10 - 0:10)
})

test_that("posterior_binary() weights external data by a0", {
  pr <- prior_beta(1, 1)
  ext <- borrow(data_binary(20, 8), 0.25)
  post <- posterior_binary(pr, 12, 5, external = ext)
  expect_equal(post$a, 1 + 0.25 * 8 + 5)
  expect_equal(post$b, 1 + 0.25 * 12 + 7)
  # a0 = 0 ignores the external data; a0 = 1 pools it
  expect_equal(posterior_binary(pr, 12, 5, borrow(data_binary(20, 8), 0)),
               posterior_binary(pr, 12, 5))
  expect_equal(posterior_binary(pr, 12, 5, borrow(data_binary(20, 8), 1)),
               posterior_binary(pr, 32, 13))
})

test_that("posterior_binary() adds hypothetical data", {
  pr <- prior_beta(2, 3)
  hyp <- hypothetical(data_binary(15, 6))
  post <- posterior_binary(pr, 0, 0, hypothetical = hyp)
  expect_equal(post, list(a = 8, b = 12))
  post2 <- posterior_binary(pr, 0, 0, borrow(data_binary(10, 4), 0.5), hyp)
  expect_equal(post2, list(a = 10, b = 15))
})
