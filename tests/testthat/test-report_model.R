model_text <- function(des) report_model(des)[[2L]][["text"]]
prior_text <- function(des) report_model(des)[[3L]][["text"]]

test_that("report_model() describes a binary endpoint and the priors", {
  expect_identical(
    report_model(ref_design("D1")),
    list(list(type = "heading", text = "Statistical model"),
         list(type = "par",
              text = paste("The number of responders among the $n_k$",
                           "subjects of group $k$ ($k = t, c$ for the",
                           "treatment and control groups) follows a binomial",
                           "distribution with response rate $\\pi_k$, and",
                           "the treatment effect is the difference in",
                           "response rates, $\\theta = \\pi_t - \\pi_c$.")),
         list(type = "par",
              text = paste("The prior distribution of the treatment group",
                           "is $\\pi_t \\sim \\textrm{Beta}(1, 1)$. The",
                           "prior distribution of the control group is",
                           "$\\pi_c \\sim \\textrm{Beta}(1, 1)$. The",
                           "posterior distributions have closed forms.")))
  )
})

test_that("report_model() describes the external data of both groups", {
  expect_identical(
    prior_text(ref_design("D2")),
    paste("The prior distribution of the treatment group is",
          "$\\pi_t \\sim \\textrm{Beta}(1, 1)$. External data of the",
          "treatment group (n = 10, y = 6) are borrowed through a power",
          "prior with weight $a_0 = 0.5$. The prior distribution of the",
          "control group is $\\pi_c \\sim \\textrm{Beta}(2, 3)$. External",
          "data of the control group (n = 20, y = 8) are borrowed through a",
          "power prior with weight $a_0 = 0.5$. The posterior distributions",
          "have closed forms.")
  )
})

test_that("report_model() describes the controls without a concurrent group", {
  expect_match(
    prior_text(ref_design("D3")),
    paste("There is no concurrent control group: the control group is a",
          "hypothetical control, whose parameters have their posterior",
          "distribution given the hypothetical control data (n = 12, y = 4),",
          "with the prior distribution $\\pi_c \\sim \\textrm{Beta}(1, 1)$."),
    fixed = TRUE
  )
  expect_match(
    prior_text(ref_design("D4")),
    paste("There is no concurrent control group: the distribution of the",
          "parameters of the control group is their posterior distribution",
          "given the external control data, with the prior distribution",
          "$\\pi_c \\sim \\textrm{Beta}(1, 1)$. External data of the control",
          "group (n = 20, y = 6) are borrowed through a power prior with",
          "weight $a_0 = 0.5$."),
    fixed = TRUE
  )
  expect_identical(
    prior_text(ref_design("D8")),
    paste("The prior distribution of the treatment group is",
          "$\\pi_t \\sim \\textrm{Beta}(1, 1)$. There is no concurrent",
          "control group: the control response rate is fixed at",
          "$\\pi_c = 0.2$. The posterior distributions have closed forms.")
  )
})

test_that("report_model() describes a continuous endpoint", {
  expect_identical(
    model_text(ref_design_cont("C2")),
    paste("The outcomes of group $k$ ($k = t, c$ for the treatment and",
          "control groups) follow a normal distribution",
          "$N(\\mu_k, \\sigma_k^2)$, and the treatment effect is the",
          "difference in means, $\\theta = \\mu_t - \\mu_c$.")
  )
  expect_match(prior_text(ref_design_cont("C2")),
               "$\\mu_t \\mid \\sigma_t^2 \\sim N(5, \\sigma_t^2 / 2)$",
               fixed = TRUE)
  c3 <- prior_text(ref_design_cont("C3"))
  expect_match(c3, paste("The prior distribution of the treatment group is",
                         "the vague prior $p(\\mu_t, \\sigma_t^2) \\propto",
                         "\\sigma_t^{-2}$."), fixed = TRUE)
  expect_match(c3, paste("There is no concurrent control group: the control",
                         "group is a hypothetical control, whose mean has",
                         "the same t distribution as the posterior",
                         "distribution of the treatment mean, with location",
                         "$\\mu_{0c} = -10$ and squared scale multiplied by",
                         "$r = 1.5$."), fixed = TRUE)
  expect_match(prior_text(ref_design_cont("C7")),
               "the control mean is fixed at $\\mu_c = -10$.", fixed = TRUE)
})

test_that("report_model() describes two binary endpoints", {
  expect_match(model_text(ref_design_bin2("B6")),
               paste("the probabilities of the responses 00, 01, 10 and 11",
                     "(the first digit for endpoint 1 and the second for",
                     "endpoint 2) are $p_k = (p_{k00}, p_{k01}, p_{k10},",
                     "p_{k11})$"), fixed = TRUE)
  b6 <- prior_text(ref_design_bin2("B6"))
  expect_match(b6, "$p_t \\sim \\textrm{Dirichlet}(0.25, 0.25, 0.25, 0.25)$",
               fixed = TRUE)
  expect_match(b6, paste("External data of the control group (n = 20,",
                         "counts (00, 01, 10, 11) = (14, 2, 3, 1)) are",
                         "borrowed through a power prior with weight",
                         "$a_0 = 0.5$."), fixed = TRUE)
  expect_match(prior_text(ref_design_bin2("BF1")),
               paste("the control response rates are fixed at",
                     "$(\\pi_{c1}, \\pi_{c2}) = (0.25, 0.5)$."), fixed = TRUE)
})

test_that("report_model() describes two continuous endpoints", {
  expect_identical(
    model_text(ref_design_cont2("K2")),
    paste("The outcomes of group $k$ ($k = t, c$ for the treatment and",
          "control groups) on the two endpoints follow a bivariate normal",
          "distribution $N_2(\\mu_k, \\Sigma_k)$ with means",
          "$\\mu_k = (\\mu_{k1}, \\mu_{k2})$, and the treatment effects are",
          "the differences in means, $\\theta_1 = \\mu_{t1} - \\mu_{c1}$ and",
          "$\\theta_2 = \\mu_{t2} - \\mu_{c2}$.")
  )
  expect_match(prior_text(ref_design_cont2("K2")),
               "$\\Sigma_t \\sim \\textrm{Inv-Wishart}(4, \\Lambda_{0t})$",
               fixed = TRUE)
  expect_match(prior_text(ref_design_cont2("K4")),
               paste("whose means have the same bivariate t distribution as",
                     "the posterior distribution of the treatment means, with",
                     "location $\\mu_{0c} = (-10, -1)$ and scale matrix",
                     "multiplied by $r = 1.5$."), fixed = TRUE)
  expect_match(prior_text(ref_design_cont2("KF1")),
               paste("the control means are fixed at",
                     "$(\\mu_{c1}, \\mu_{c2}) = (-10, -1)$."), fixed = TRUE)
})

test_that("report_model() describes a time-to-event endpoint", {
  expect_identical(
    model_text(ref_design_tte("T2")),
    paste("The event times of group $k$ ($k = t, c$ for the treatment and",
          "control groups) follow an exponential distribution with hazard",
          "$\\lambda_k$, and the treatment effect is the difference in the",
          "survival probabilities at the milestone time $t_0 = 3$,",
          "$\\theta = S_t(t_0) - S_c(t_0)$ with",
          "$S_k(t_0) = \\exp(-\\lambda_k t_0)$.")
  )
  expect_match(prior_text(ref_design_tte("T1")),
               paste("the control survival probability at the milestone time",
                     "is fixed at $S_c(t_0) = 0.72$."), fixed = TRUE)
  expect_match(prior_text(ref_design_tte("T5")),
               paste("given the hypothetical control data (events = 8, total",
                     "observed time = 60)"), fixed = TRUE)
  expect_match(model_text(ref_design_hr("H1")),
               paste("the treatment effect is the hazard ratio,",
                     "$\\theta = \\lambda_t / \\lambda_c$."), fixed = TRUE)
  expect_match(prior_text(ref_design_hr("H2")),
               paste("the control median survival time is fixed at 6, so",
                     "that the control hazard is",
                     "$\\lambda_c = \\log(2) / 6$."), fixed = TRUE)
})
