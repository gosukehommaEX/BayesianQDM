test_that("oc_two_grid() places the scenarios on the effects", {
  oc <- qdm_oc(ref_design_bin2("B2"), scen_bin2)
  g <- oc_two_grid(oc, ref_design_bin2("B2"))
  expect_identical(names(g), c("x", "y", "labs", "decisions"))
  expect_identical(as.character(g$x),
                   format(round(scen_bin2$pi_t1 - scen_bin2$pi_c1, 10),
                          trim = TRUE))
  expect_identical(levels(g$y),
                   format(sort(unique(round(scen_bin2$pi_t2 -
                                              scen_bin2$pi_c2, 10))),
                          trim = TRUE))
  expect_identical(g$labs, c("Treatment effect on endpoint 1",
                             "Treatment effect on endpoint 2"))
})

test_that("oc_two_grid() uses the treatment parameters without a control group", {
  x <- data.frame(pi_t1 = c(0.5, 0.3), pi_t2 = c(0.2, 0.2),
                  Miss = c(0, 0.1))
  g <- oc_two_grid(x, ref_design_bin2("B7"))
  expect_identical(g$labs, c("Treatment group, endpoint 1",
                             "Treatment group, endpoint 2"))
  expect_identical(levels(g$x), c("0.3", "0.5"))
  expect_identical(as.integer(g$x), c(2L, 1L))
  # Miss decisions are shown only if they have positive probability
  expect_identical(g$decisions, c("Go", "NoGo", "Gray", "Miss"))
  x$Miss <- c(0, NA)
  expect_identical(oc_two_grid(x, ref_design_bin2("B7"))$decisions,
                   c("Go", "NoGo", "Gray"))
  # Two continuous endpoints
  y <- data.frame(mu_t1 = c(1, 2), mu_t2 = c(0, 0), Miss = 0)
  expect_identical(levels(oc_two_grid(y, ref_design_cont2("K4"))$x),
                   c("1", "2"))
})

test_that("oc_two_grid() keeps levels that differ beyond 7 digits", {
  x <- data.frame(pi_t1 = c(0.3, 0.30000001), pi_t2 = c(0.2, 0.2), Miss = 0)
  g <- oc_two_grid(x, ref_design_bin2("B7"))
  expect_identical(nlevels(g$x), 2L)
  expect_identical(levels(g$x), c("0.30000000", "0.30000001"))
})

test_that("oc_two_grid() needs distinct values on the axes", {
  x <- data.frame(pi_t1 = c(0.5, 0.5), pi_t2 = c(0.2, 0.2), Miss = 0)
  expect_error(oc_two_grid(x, ref_design_bin2("B7")),
               "with different correlations or control")
})
