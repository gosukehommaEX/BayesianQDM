test_that("report_decision() gives the data, the criteria and the decision", {
  dec <- qdm_decide(ref_design("D1"), data_binary(12, 9), data_binary(12, 4))
  prob <- function(p) formatC(p, format = "f", digits = 3)
  expect_true(dec$met_go)
  expect_false(dec$met_nogo)
  expect_gt(dec$prob_go, 0.5)
  expect_lt(dec$prob_nogo, 0.8)
  expect_identical(dec$decision, "Go")
  expect_identical(
    report_decision(dec, 3L),
    list(list(type = "heading", text = "Decision for the observed data"),
         list(type = "par",
              text = paste("The observed data are 9 responders among 12",
                           "subjects in the treatment group and 4",
                           "responders among 12 subjects in the control",
                           "group.")),
         list(type = "list",
              items = c(sprintf(paste("Go criterion (met): $\\Pr(\\theta >",
                                      "0.2 \\mid D) = %s \\ge 0.5$"),
                                prob(dec$prob_go)),
                        sprintf(paste("NoGo criterion (not met):",
                                      "$\\Pr(\\theta \\le 0.1 \\mid D) = %s <",
                                      "0.8$"), prob(dec$prob_nogo)))),
         list(type = "par", text = "The decision is Go."),
         list(type = "par",
              text = paste("The posterior distribution of $\\pi_t$ is",
                           "$\\textrm{Beta}(10, 4)$, and that of $\\pi_c$ is",
                           "$\\textrm{Beta}(5, 9)$.")))
  )
})

test_that("report_decision() gives the data of the treatment group only", {
  dec <- qdm_decide(ref_design("D8"), data_binary(20, 9))
  res <- report_decision(dec, 3L)
  expect_identical(res[[2L]][["text"]],
                   paste("The observed data of the treatment group are 9",
                         "responders among 20 subjects."))
  # The fixed control value is given in the section on the statistical model
  expect_identical(res[[length(res)]],
                   list(type = "par",
                        text = paste("The posterior distribution of $\\pi_t$",
                                     "is $\\textrm{Beta}(10, 12)$.")))
})

test_that("report_decision() gives the probabilities of the regions", {
  d <- ref_data_cont2$K1
  dec <- qdm_decide(ref_design_cont2("K1"), d$t, d$c)
  prob <- function(p) formatC(p, format = "f", digits = 2)
  res <- report_decision(dec, 2L)
  expect_identical(
    res[[3L]][["items"]][1L],
    sprintf(paste("Go criterion (%s): $\\Pr\\{(\\theta_1, \\theta_2) \\in",
                  "R_{1} \\mid D\\} = %s %s 0.5$"),
            if (dec$met_go) "met" else "not met", prob(dec$prob_go),
            if (dec$prob_go >= 0.5) "\\ge" else "<")
  )
  expect_identical(res[[5L]][["text"]], "The probabilities of the regions are:")
  expect_identical(
    res[[6L]],
    list(type = "table", header = c("Region", "Probability"),
         rows = cbind(sprintf("$R_{%d}$", 1:9),
                      prob(unname(dec$prob_regions))))
  )
})

test_that("report_decision() combines the criteria of each endpoint", {
  d <- ref_data_cont2$K1
  vague <- prior_vague()
  des <- design_continuous2(12, 12, vague, vague,
                            crit_cont2(c(0.6, 0.5), c(0.5, 0.6), rule = "or"))
  dec <- qdm_decide(des, d$t, d$c)
  prob <- function(p) formatC(p, format = "f", digits = 3)
  op <- ifelse(dec$prob_go >= c(0.6, 0.5), "\\ge", "<")
  expect_identical(
    report_decision(dec, 3L)[[3L]][["items"]][1L],
    sprintf(paste("Go criterion (%s): $\\Pr(\\theta_1 > 10 \\mid D) = %s %s",
                  "0.6$ or $\\Pr(\\theta_2 > 2 \\mid D) = %s %s 0.5$"),
            if (dec$met_go) "met" else "not met", prob(dec$prob_go[1L]),
            op[1L], prob(dec$prob_go[2L]), op[2L])
  )
})

test_that("report_decision() gives the posterior distributions", {
  d <- ref_data_cont2$K1
  dec <- qdm_decide(ref_design_cont2("K1"), d$t, d$c)
  res <- report_decision(dec, 3L)
  post_t <- report_posterior(dec$design, dec$posterior_t, "t")
  post_c <- report_posterior(dec$design, dec$posterior_c, "c")
  expect_identical(
    res[[length(res)]][["text"]],
    sprintf("The posterior distribution of %s is %s, and that of %s is %s.",
            post_t[1L], post_t[2L], post_c[1L], post_c[2L])
  )
  expect_match(post_c[2L], "^the bivariate t distribution with ")
})
