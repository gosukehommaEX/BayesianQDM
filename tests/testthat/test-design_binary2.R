test_that("design_binary2() stores the specification", {
  des <- ref_design_bin2("B6")
  expect_s3_class(des, c("qdm_design_binary2", "qdm_design"), exact = TRUE)
  expect_identical(des$endpoint, "binary2")
  expect_identical(des$n_t, 8L)
  expect_identical(des$control, "concurrent")
  expect_identical(des$nMC, 20000L)
  expect_true(des$enumerate)
  expect_equal(des$external_c$a0, 0.5)
  expect_identical(ref_design_bin2("B7")$control, "hypothetical")
  expect_identical(ref_design_bin2("B8")$control, "external")
  expect_false(ref_design_bin2("B1", enumerate = FALSE)$enumerate)
})

test_that("design_binary2() validates its arguments", {
  q <- dir_q()
  cr <- criterion_posterior(c(0.15, 0.15), c(0.1, 0.1), 0.3, 0.4)
  expect_error(design_binary2(0, 5, q, q, cr), "'n_t' must be at least 1.")
  expect_error(design_binary2(5, 5, prior_beta(1, 1), q, cr),
               "'prior_t' must be created by prior_dirichlet\\(\\).")
  expect_error(design_binary2(5, 5, q, NULL, cr),
               "'prior_c' must be created by prior_dirichlet\\(\\).")
  expect_error(design_binary2(5, 5, q, q, crit_post()),
               "for two endpoints.")
  expect_error(design_binary2(5, 5, q, q, cr,
                              external_t = borrow(data_binary(5, 2), 0.5)),
               "'external_t' must be NULL or created by borrow\\(\\) with")
  expect_error(design_binary2(5, 0, q, q, cr,
                              hypothetical_c = hypothetical(mean = 0,
                                                            r = 1)),
               "'hypothetical_c' must be NULL or created by hypothetical")
  expect_error(design_binary2(5, 0, q, q, cr),
               "When 'n_c' is 0, specify 'external_c'")
  expect_error(design_binary2(5, 5, q, q, cr, nMC = 0),
               "'nMC' must be at least 1 and at most 2147483647.")
  expect_error(design_binary2(5, 5, q, q, cr, nMC = 3e9),
               "'nMC' must be at least 1 and at most 2147483647.")
  expect_error(design_binary2(5, 5, q, q, cr, enumerate = NA),
               "'enumerate' must be TRUE or FALSE.")
})

test_that("design_binary2() stores and checks a fixed control", {
  des <- ref_design_bin2("BF1")
  expect_identical(des$control, "fixed")
  expect_identical(des$fixed_c, c(0.25, 0.5))
  expect_null(des$prior_c)
  q <- dir_q()
  cr <- criterion_posterior(c(0.15, 0.15), c(0.1, 0.1), 0.3, 0.4)
  expect_error(design_binary2(5, 0, q, q, cr, fixed_c = c(0.2, 0.3)),
               "'prior_c' must be NULL for a fixed control")
  expect_error(design_binary2(5, 0, q, criterion = cr, fixed_c = 0.2),
               "'fixed_c' must have length 2.")
  expect_error(design_binary2(5, 0, q, criterion = cr, fixed_c = c(0.2, 2)),
               "'fixed_c' must be at least 0 and at most 1.")
})
