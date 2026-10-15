test_that("design_continuous2() stores the specification", {
  des <- ref_design_cont2("K3")
  expect_s3_class(des, c("qdm_design_continuous2", "qdm_design"),
                  exact = TRUE)
  expect_identical(des$endpoint, "continuous2")
  expect_identical(des$n_t, 12L)
  expect_identical(des$control, "concurrent")
  expect_identical(des$method, "MM")
  expect_identical(des$nMC, 10000L)
  expect_equal(des$external_t$a0, 0.5)
  expect_identical(ref_design_cont2("K4")$control, "hypothetical")
  expect_identical(ref_design_cont2("K5")$control, "external")
  expect_identical(ref_design_cont2("K1", "MC", 500)$nMC, 500L)
  expect_null(ref_design_cont2("K6")$prior_c)
})

test_that("design_continuous2() validates its arguments", {
  v <- prior_vague()
  cr <- crit_cont2()
  expect_error(design_continuous2(2, 12, v, v, cr),
               "'n_t' must be at least 3.")
  expect_error(design_continuous2(12, 2, v, v, cr),
               "'n_c' must be 0 or a whole number of at least 3.")
  expect_error(design_continuous2(12, 2 + 1e-9, v, v, cr),
               "'n_c' must be 0 or a whole number of at least 3.")
  expect_error(design_continuous2(12, 12, prior_nix(0, 1, 2, 3), v, cr),
               "'prior_t' must be created by prior_niw\\(\\) or prior_vague")
  expect_error(design_continuous2(12, 12, v, NULL, cr),
               "'prior_c' must be created by prior_niw\\(\\) or prior_vague")
  expect_error(design_continuous2(12, 12, v, v, crit_cont()),
               "for two endpoints.")
  expect_error(design_continuous2(12, 12, v, v, cr,
                                  external_t = borrow(data_continuous(5, 1,
                                                                      2),
                                                      0.5)),
               "'external_t' must be NULL or created by borrow\\(\\) with")
  expect_error(design_continuous2(12, 0, v, NULL, cr,
                                  hypothetical_c = hypothetical(mean = 0,
                                                                r = 1)),
               "with 'mean' of length 2 and 'r'.")
  expect_error(design_continuous2(12, 0, v, v, cr,
                                  hypothetical_c = hypothetical(
                                    mean = c(0, 0), r = 1
                                  )),
               "'prior_c' must be NULL for a hypothetical control")
  expect_error(design_continuous2(12, 0, v, v, cr),
               "When 'n_c' is 0, specify 'external_c'")
  expect_error(design_continuous2(12, 12, v, v, cr, method = "NI"),
               "'arg' should be one of")
  expect_error(design_continuous2(12, 12, v, v, cr, nMC = 0),
               "'nMC' must be at least 1 and at most 2147483647.")
  # External control only with the vague prior needs a0 * n_e > 2
  ext <- borrow(data_continuous2(4, c(0, 0), cov2(1, 1, 0)), 0.5)
  expect_error(design_continuous2(12, 0, v, v, cr, external_c = ext,
                                  method = "MC"),
               "must be greater than 2, so that the distribution")
  # Moment matching needs more than 4 degrees of freedom (n - 2 > 4 for the
  # vague prior)
  expect_error(design_continuous2(6, 12, v, v, cr),
               "\\(here 4 and 10\\). Use method = \"MC\".")
  expect_s3_class(design_continuous2(6, 12, v, v, cr, method = "MC"),
                  "qdm_design_continuous2")
})

test_that("design_continuous2() stores and checks a fixed control", {
  des <- ref_design_cont2("KF2")
  expect_identical(des$control, "fixed")
  expect_identical(des$fixed_c, c(1, 0.5))
  expect_null(des$prior_c)
  vague <- prior_vague()
  expect_error(design_continuous2(10, 0, vague, vague, crit_cont2(),
                                  fixed_c = c(0, 0)),
               "'prior_c' must be NULL for a fixed control")
  expect_error(design_continuous2(10, 0, vague, criterion = crit_cont2(),
                                  fixed_c = 0),
               "'fixed_c' must have length 2.")
})
