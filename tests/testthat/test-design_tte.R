test_that("design_tte() stores the specification", {
  des <- ref_design_tte("T2")
  expect_s3_class(des, c("qdm_design_tte", "qdm_design"), exact = TRUE)
  expect_identical(des$endpoint, "tte")
  expect_identical(des$n_t, 20L)
  expect_identical(des$n_c, 20L)
  expect_identical(des$milestone, 3)
  expect_identical(des$control, "concurrent")
  expect_identical(des$conduct, conduct_example)
  expect_null(des$conduct_future)
  expect_identical(des$nMC, 10000L)
  expect_identical(ref_design_tte("T1")$control, "fixed")
  expect_identical(ref_design_tte("T5")$control, "hypothetical")
  expect_identical(ref_design_tte("T6")$control, "external")
  # The future study of a predictive criterion is conducted like the study
  expect_identical(ref_design_tte("T4")$conduct_future, conduct_complete)
  # The conduct is needed only for the operating characteristics
  expect_null(design_tte(10, 0, prior_gamma(1, 1), criterion = crit_tte(),
                         milestone = 3, fixed_c = 0.5)$conduct)
})

test_that("design_tte() validates its arguments", {
  g <- prior_gamma(0.1, 1)
  cr <- crit_tte()
  expect_error(design_tte(0, 10, g, g, cr, 3), "'n_t' must be at least 1.")
  expect_error(design_tte(10, 10, prior_beta(1, 1), g, cr, 3),
               "'prior_t' must be created by prior_gamma\\(\\).")
  expect_error(design_tte(10, 10, g, NULL, cr, 3),
               "'prior_c' must be created by prior_gamma\\(\\).")
  expect_error(design_tte(10, 10, g, g, list(), 3),
               "'criterion' must be created by criterion_posterior\\(\\) or")
  expect_error(design_tte(10, 10, g, g, cr, 0),
               "'milestone' must be greater than 0.")
  expect_error(design_tte(10, 10, g, g, cr, 3,
                          external_t = borrow(data_binary(5, 2), 0.5)),
               "'external_t' must be NULL or created by borrow\\(\\) with")
  expect_error(design_tte(10, 0, g, g, cr, 3,
                          hypothetical_c = hypothetical(data_binary(5, 2))),
               "'hypothetical_c' must be NULL or created by hypothetical")
  expect_error(design_tte(10, 0, g, g, cr, 3, fixed_c = 0.5),
               "'prior_c' must be NULL for a fixed control")
  expect_error(design_tte(10, 0, g, criterion = cr, milestone = 3,
                          fixed_c = 1.5),
               "'fixed_c' must be at least 0 and at most 1.")
  expect_error(design_tte(10, 10, g, g, cr, 3, conduct = list()),
               "'conduct' must be created by conduct_tte\\(\\).")
  expect_error(design_tte(10, 10, g, g, cr, 3,
                          conduct = conduct_tte(c(0, 5), cut_events = 21)),
               "is at event 21, but the study has only 20 subjects.")
  expect_error(design_tte(10, 10, g, g, cr, 3, nMC = 0),
               "'nMC' must be at least 1")
})

test_that("design_tte() checks the conduct of the future study", {
  g <- prior_gamma(0.1, 1)
  cr <- criterion_predictive(0, m_t = 10, m_c = 8, gamma_go = 0.8,
                             gamma_nogo = 0.3)
  expect_error(design_tte(10, 10, g, g, cr, 3),
               "A predictive criterion needs 'conduct_future'")
  expect_error(design_tte(10, 10, g, g, cr, 3,
                          conduct_future = conduct_tte(c(0, 5),
                                                       cut_events = 19)),
               "'conduct_future' is at event 19, but the study has only 18")
  # Without a control group in the future study for a fixed control
  expect_error(design_tte(10, 0, g, criterion = cr, milestone = 3,
                          fixed_c = 0.5,
                          conduct_future = conduct_tte(c(0, 5),
                                                       cut_events = 11)),
               "has only 10 subjects.")
  co <- conduct_tte(c(0, 5), cut_events = 10)
  des <- design_tte(10, 0, g, criterion = cr, milestone = 3, fixed_c = 0.5,
                    conduct_future = co)
  expect_identical(des$conduct_future, co)
  expect_null(des$conduct)
})

test_that("design_tte() stores the effect measure", {
  expect_identical(ref_design_tte("T1")$effect, "milestone")
  des <- ref_design_hr("H1")
  expect_identical(des$effect, "hazard_ratio")
  expect_null(des$milestone)
  expect_identical(ref_design_hr("H2")$fixed_c, 6)
  expect_identical(ref_design_hr("H2")$control, "fixed")
})

test_that("design_tte() validates the arguments of the hazard ratio", {
  g <- prior_gamma(0.1, 1)
  cr <- crit_hr()
  expect_error(design_tte(10, 10, g, g, cr),
               "'milestone' must be given for effect = \"milestone\".",
               fixed = TRUE)
  expect_error(design_tte(10, 10, g, g, cr, milestone = 3,
                          effect = "hazard_ratio"),
               "'milestone' is not used for effect = \"hazard_ratio\".",
               fixed = TRUE)
  expect_error(design_tte(10, 10, g, g, crit_tte(), effect = "hazard_ratio"),
               paste("The thresholds of the criterion must be positive",
                     "hazard ratios for effect = \"hazard_ratio\"."),
               fixed = TRUE)
  expect_error(design_tte(10, 10, g, g,
                          criterion_predictive(0, m_t = 10, gamma_go = 0.8,
                                               gamma_nogo = 0.3),
                          effect = "hazard_ratio"),
               "must be positive hazard ratios")
  expect_error(design_tte(10, 0, g, criterion = cr, effect = "hazard_ratio",
                          fixed_c = 0),
               "'fixed_c' must be greater than 0.")
  expect_error(design_tte(10, 10, g, g, cr, effect = "ratio"),
               "'arg' should be one of")
  # A median survival time larger than 1 is accepted as a fixed control
  expect_identical(design_tte(10, 0, g, criterion = cr,
                              effect = "hazard_ratio",
                              fixed_c = 8.5)$fixed_c, 8.5)
})
