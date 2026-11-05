test_that("decision_subtitle() gives the probabilities of the criteria", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  p <- function(v) formatC(unname(v), format = "f", digits = 3L)
  expect_identical(decision_subtitle(dec), paste(
    sprintf("Go: P(theta > 0.2 | data) = %s (threshold 0.5, met)",
            p(dec$prob_go)),
    sprintf("NoGo: P(theta <= 0.1 | data) = %s (threshold 0.8, not met)",
            p(dec$prob_nogo)),
    sep = "\n"
  ))
  expect_identical(
    strsplit(decision_subtitle(dec, digits = 2L), "\n")[[1]][1L],
    sprintf("Go: P(theta > 0.2 | data) = %s (threshold 0.5, met)",
            formatC(unname(dec$prob_go), format = "f", digits = 2L))
  )
  # Two endpoints, joint rule
  two <- qdm_decide(ref_design_cont2("K1"), ref_data_cont2$K1$t,
                    ref_data_cont2$K1$c)
  expect_match(decision_subtitle(two),
               "^Go: P\\(regions \\{1\\} \\| data\\) = ")
})

test_that("decision_subtitle() gives each endpoint of a marginal rule", {
  dec <- qdm_decide(ref_design_cont2("K5"), ref_data_cont2$K5$t)
  p <- function(v) formatC(unname(v), format = "f", digits = 3L)
  met <- function(m) if (m) "met" else "not met"
  expect_identical(strsplit(decision_subtitle(dec), "\n")[[1]], c(
    sprintf(paste("Go (rule \"and\", %s): P(theta1 > 10 | data) = %s",
                  "(threshold 0.6); P(theta2 > 2 | data) = %s",
                  "(threshold 0.5)"),
            met(dec$met_go), p(dec$prob_go[1L]), p(dec$prob_go[2L])),
    sprintf(paste("NoGo (rule \"and\", %s): P(theta1 <= 5 | data) = %s",
                  "(threshold 0.5); P(theta2 <= 1 | data) = %s",
                  "(threshold 0.6)"),
            met(dec$met_nogo), p(dec$prob_nogo[1L]), p(dec$prob_nogo[2L]))
  ))
})
