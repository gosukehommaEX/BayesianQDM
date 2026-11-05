test_that("spec_set() replaces the arguments that it is given", {
  spec <- spec_set(qdm_endpoint("binary"), list(n_t = 12, n_c = 10), "f")
  expect_identical(spec$args[["n_t"]], 12)
  expect_identical(spec$args[["n_c"]], 10)
  spec <- spec_set(spec, list(n_t = 14, prior_c = NULL), "f")
  expect_identical(spec$args[["n_t"]], 14)
  expect_identical(spec$args[["n_c"]], 10)
  expect_setequal(names(spec$args), c("n_t", "n_c"))
})

test_that("spec_set() removes an argument given as NULL", {
  spec <- spec_set(qdm_endpoint("binary"), list(n_t = 12, n_c = 10), "f")
  spec <- spec_set(spec, list(n_c = NULL), "f")
  expect_identical(names(spec$args), "n_t")
})

test_that("spec_set() rejects arguments that the design function does not have", {
  expect_error(spec_set(qdm_endpoint("binary"),
                        list(nMC = 100L, method = NULL), "f"),
               "f() cannot set 'nMC', which design_binary() does not have.",
               fixed = TRUE)
  expect_error(spec_set(qdm_endpoint("continuous"),
                        list(enumerate = TRUE, milestone = 3), "f"),
               paste("f() cannot set 'enumerate', 'milestone', which",
                     "design_continuous() does not have."),
               fixed = TRUE)
  expect_error(spec_set(list(), list(n_t = 1), "f"),
               "'x' must be a design specification", fixed = TRUE)
})

test_that("spec_set() starts from the specification of a design", {
  spec <- spec_set(ref_design("D1"), list(n_t = 20), "f")
  expect_s3_class(spec, "qdm_spec", exact = TRUE)
  expect_identical(spec$args$n_t, 20)
  expect_identical(spec$args$n_c, 12L)
  expect_identical(spec$args$criterion, crit_post())
})
