test_that("print.qdm_criterion() prints each line and returns invisibly", {
  cr <- criterion_posterior(0.2, 0.1, gamma_go = 0.5, gamma_nogo = 0.8)
  expect_output(res <- withVisible(print(cr)),
                "Posterior probability criterion\n  Go:")
  expect_false(res$visible)
  expect_identical(res$value, cr)
})
