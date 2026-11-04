comp_text <- function(des) report_computation(des)[[2L]][["text"]]
computed_by <- function(how) {
  sprintf("The probabilities of the decision criteria are computed by %s.",
          how)
}

test_that("report_computation() names the method of one binary endpoint", {
  expect_identical(
    report_computation(ref_design("D1")),
    list(list(type = "heading", text = "Computation of the probabilities"),
         list(type = "par", text = computed_by("numerical integration")))
  )
  expect_identical(comp_text(ref_design("D5")), computed_by("exact summation"))
  expect_identical(comp_text(ref_design("D8")),
                   computed_by("exact computation (fixed control)"))
})

test_that("report_computation() names the method of continuous endpoints", {
  expect_identical(comp_text(ref_design_cont("C1")),
                   computed_by("numerical integration"))
  expect_identical(comp_text(ref_design_cont("C1", method = "MM")),
                   computed_by("moment-matching approximation"))
  expect_identical(comp_text(ref_design_cont("C1", method = "MC",
                                             nMC = 2000L)),
                   computed_by("Monte Carlo simulation (nMC = 2000)"))
  expect_identical(comp_text(ref_design_cont("C7")),
                   computed_by("exact computation (fixed control)"))
  expect_identical(comp_text(ref_design_cont2("K1")),
                   computed_by("moment-matching approximation"))
})

test_that("report_computation() names the method of two binary endpoints", {
  expect_identical(comp_text(ref_design_bin2("B6")),
                   computed_by("Monte Carlo simulation (nMC = 20000)"))
  expect_identical(comp_text(ref_design_bin2("B1")),
                   computed_by("exact computation"))
})

test_that("report_computation() describes the conduct of a time-to-event study", {
  expect_identical(
    report_computation(ref_design_tte("T2")),
    list(list(type = "heading", text = "Computation of the probabilities"),
         list(type = "par", text = computed_by("numerical integration")),
         list(type = "par",
              text = paste("The operating characteristics are computed by",
                           "simulating the conduct of the study:")),
         list(type = "list", items = format(conduct_example)))
  )
  res <- report_computation(ref_design_tte("T4", nMC = 500L))
  expect_length(res, 6L)
  expect_identical(res[[2L]][["text"]],
                   computed_by("simulation of future studies (nMC = 500)"))
  expect_identical(res[[5L]][["text"]],
                   "The future studies are simulated with:")
  expect_identical(res[[6L]][["items"]], format(conduct_complete))
  g <- prior_gamma(0.1, 1)
  expect_length(report_computation(design_tte(20, 20, g, g, crit_tte(),
                                              milestone = 3)), 2L)
  expect_identical(comp_text(ref_design_hr("H1")),
                   computed_by("exact computation (F distribution)"))
  expect_identical(comp_text(ref_design_tte("T1")),
                   computed_by("exact computation (fixed control)"))
})
