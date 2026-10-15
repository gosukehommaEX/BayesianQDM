test_that("posterior_tte() adds events and observed times to the prior", {
  pr <- prior_gamma(0.1, 1)
  expect_identical(posterior_tte(pr, 5, 40), list(shape = 5.1, rate = 41))
  expect_identical(posterior_tte(pr, c(0, 3), c(0, 10)),
                   list(shape = c(0.1, 3.1), rate = c(1, 11)))
  # External data with a0 = 0.5, hypothetical data with full weight
  post <- posterior_tte(pr, 6, 150, external = borrow(data_tte(4, 30), 0.5),
                        hypothetical = hypothetical(data_tte(8, 60)))
  expect_equal(post, list(shape = 0.1 + 2 + 8 + 6, rate = 1 + 15 + 60 + 150))
})
