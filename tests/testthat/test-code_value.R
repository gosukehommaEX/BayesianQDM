test_that("code_value() writes the objects of the package", {
  flat <- function(x) code_flat(code_value(x))
  expect_identical(flat(prior_beta(0.5, 0.5)), "prior_beta(a = 0.5, b = 0.5)")
  expect_identical(flat(prior_beta(1L, 2)), "prior_beta(a = 1L, b = 2)")
  expect_identical(flat(prior_vague()), "prior_vague()")
  expect_identical(flat(prior_gamma(0.1, 1)),
                   "prior_gamma(shape = 0.1, rate = 1)")
  expect_identical(flat(prior_nix(5, 2, 4, 20)),
                   "prior_nix(mu0 = 5, kappa0 = 2, nu0 = 4, sigma0 = 20)")
  expect_identical(flat(prior_dirichlet(rep(0.25, 4))),
                   "prior_dirichlet(alpha = c(0.25, 0.25, 0.25, 0.25))")
  expect_identical(flat(prior_niw(c(0, 1), 0.1, 4, lambda_cs4)),
                   paste("prior_niw(mu0 = c(0, 1), kappa0 = 0.1, nu0 = 4,",
                         "Lambda0 = matrix(c(25, 1, 1, 9), nrow = 2))"))
  expect_identical(flat(data_binary(12, 4)), "data_binary(n = 12, y = 4)")
  expect_identical(flat(data_binary2(c(2, 1, 1, 0))),
                   "data_binary2(counts = c(2, 1, 1, 0))")
  expect_identical(flat(data_continuous(6, 20, 30)),
                   "data_continuous(n = 6, mean = 20, sd = 30)")
  expect_identical(flat(data_continuous2(8, c(2, 1),
                                         matrix(c(25, 3, 3, 9), 2L))),
                   paste("data_continuous2(n = 8, mean = c(2, 1),",
                         "cov = matrix(c(25, 3, 3, 9), nrow = 2))"))
  expect_identical(flat(borrow(data_tte(4, 30), 0.5)),
                   "borrow(data = data_tte(events = 4, exposure = 30), a0 = 0.5)")
  expect_identical(flat(hypothetical(mean = c(-10, -1), r = 1.5)),
                   "hypothetical(mean = c(-10, -1), r = 1.5)")
  expect_identical(flat(hypothetical(data_binary(12, 4))),
                   "hypothetical(data = data_binary(n = 12, y = 4))")
})

test_that("code_value() writes vectors and rejects other values", {
  expect_identical(code_value(c("greater", "less")), "c(\"greater\", \"less\")")
  expect_identical(code_value("a\"b"), "\"a\\\"b\"")
  expect_identical(code_value(TRUE), "TRUE")
  expect_identical(code_value(c(TRUE, FALSE)), "c(TRUE, FALSE)")
  expect_identical(code_value(3L), "3L")
  expect_identical(code_value(3L, keep_type = FALSE), "3")
  expect_error(code_value(list()), "R code cannot be written for this value.",
               fixed = TRUE)
})

test_that("code_value() gives code that recreates the objects", {
  objects <- list(prior_beta(1L, 2), prior_niw(c(0, 1), 0.1, 4, lambda_cs2),
                  prior_dirichlet(c(1, 2, 3, 4)),
                  borrow(data_continuous2(20, c(-10, -1), cov2(70, 5, 0.3)),
                         0.5),
                  hypothetical(data_tte(15, 100)), conduct_example)
  for (obj in objects) {
    expect_identical(eval(parse(text = code_flat(code_value(obj)))), obj)
  }
})
