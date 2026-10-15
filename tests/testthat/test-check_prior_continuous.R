test_that("check_prior_continuous() accepts NIX and vague priors", {
  p <- prior_nix(0, 1, 3, 2)
  expect_identical(check_prior_continuous(p, "prior_t"), p)
  expect_invisible(check_prior_continuous(prior_vague(), "prior_t"))
})

test_that("check_prior_continuous() rejects other objects", {
  expect_error(check_prior_continuous(prior_beta(1, 1), "prior_c"),
               paste("'prior_c' must be created by prior_nix\\(\\) or",
                     "prior_vague\\(\\)."))
  expect_error(check_prior_continuous(NULL, "prior_c"),
               "'prior_c' must be created by prior_nix")
})
