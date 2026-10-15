test_that("design_binary() stores a design with a concurrent control", {
  pr <- prior_beta(0.5, 0.5)
  cr <- crit_post()
  ext <- borrow(data_binary(20, 8), 0.5)
  des <- design_binary(12, 10, pr, prior_beta(1, 1), cr, external_c = ext)
  expect_s3_class(des, c("qdm_design_binary", "qdm_design"), exact = TRUE)
  expect_identical(des$endpoint, "binary")
  expect_identical(des$n_t, 12L)
  expect_identical(des$n_c, 10L)
  expect_identical(des$prior_t, pr)
  expect_identical(des$criterion, cr)
  expect_identical(des$external_c, ext)
  expect_null(des$external_t)
  expect_null(des$hypothetical_c)
  expect_identical(des$control, "concurrent")
})

test_that("design_binary() recognises the three forms of the control group", {
  expect_identical(ref_design("D1")$control, "concurrent")
  expect_identical(ref_design("D3")$control, "hypothetical")
  expect_identical(ref_design("D4")$control, "external")
  expect_identical(ref_design("D5")$criterion$type, "predictive")
})

test_that("design_binary() rejects invalid specifications", {
  pr <- prior_beta(1, 1)
  cr <- crit_post()
  expect_error(design_binary(0, 10, pr, pr, cr),
               "'n_t' must be at least 1.")
  expect_error(design_binary(10, 2.5, pr, pr, cr),
               "'n_c' must contain whole numbers only.")
  expect_error(design_binary(10, 10, prior_nix(0, 1, 3, 1), pr, cr),
               "'prior_t' must be created by prior_beta().")
  expect_error(design_binary(10, 10, pr, list(a = 1, b = 1), cr),
               "'prior_c' must be created by prior_beta().")
  expect_error(design_binary(10, 10, pr, pr, list()),
               "'criterion' must be created by criterion_posterior()")
  expect_error(design_binary(10, 10, pr, pr,
                             crit_post(tv = c(0.2, 0.2), mav = c(0, 0))),
               "for one endpoint.")
  expect_error(design_binary(10, 10, pr, pr, cr,
                             external_t = data_binary(10, 5)),
               "'external_t' must be NULL or created by borrow\\(\\)")
  expect_error(design_binary(10, 10, pr, pr, cr,
                             external_c = borrow(data_continuous(10, 1, 1),
                                                 0.5)),
               "'external_c' must be NULL or created by borrow\\(\\)")
  expect_error(design_binary(10, 0, pr, pr, cr,
                             hypothetical_c = hypothetical(mean = 0, r = 1)),
               "'hypothetical_c' must be NULL or created by hypothetical\\(\\)")
  expect_error(design_binary(10, 0, pr, pr, cr,
                             hypothetical_c = data_binary(10, 3)),
               "'hypothetical_c' must be NULL or created by hypothetical\\(\\)")
})

test_that("design_binary() checks the combination of control arguments", {
  pr <- prior_beta(1, 1)
  cr <- crit_post()
  hyp <- hypothetical(data_binary(12, 4))
  ext <- borrow(data_binary(20, 6), 0.5)
  expect_error(design_binary(10, 10, pr, pr, cr, hypothetical_c = hyp),
               "'hypothetical_c' can be used only when 'n_c' is 0.")
  expect_error(design_binary(10, 0, pr, pr, cr),
               "When 'n_c' is 0, specify 'external_c'")
  expect_error(design_binary(10, 0, pr, pr, cr, external_c = ext,
                             hypothetical_c = hyp),
               "specify only one of 'external_c', 'hypothetical_c' and")
})

test_that("design_binary() stores a fixed control", {
  des <- ref_design("D8")
  expect_identical(des$control, "fixed")
  expect_identical(des$fixed_c, 0.2)
  expect_null(des$prior_c)
  expect_null(ref_design("D1")$fixed_c)
})

test_that("design_binary() checks the fixed control", {
  pr <- prior_beta(1, 1)
  cr <- crit_post()
  expect_error(design_binary(10, 0, pr, pr, cr, fixed_c = 0.2),
               "'prior_c' must be NULL for a fixed control")
  expect_error(design_binary(10, 0, pr, criterion = cr),
               "When 'n_c' is 0, specify 'external_c'")
  expect_error(design_binary(10, 0, pr, criterion = cr, fixed_c = 1.5),
               "'fixed_c' must be at least 0 and at most 1.")
  expect_error(design_binary(10, 10, pr, pr, cr, fixed_c = 0.2),
               "'fixed_c' can be used only when 'n_c' is 0.")
  expect_error(design_binary(10, 0, pr, pr, cr,
                             external_c = borrow(data_binary(20, 6), 0.5),
                             fixed_c = 0.2),
               "specify only one of")
  expect_error(design_binary(10, 10, pr, criterion = cr),
               "'prior_c' must be created by prior_beta().")
})
