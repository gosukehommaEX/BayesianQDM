#' Statistical model section of a report
#'
#' @param design A design.
#'
#' @return A list of blocks: the model and the treatment effect, then the
#'   prior distributions, the control group and the external data.
#' @noRd
report_model <- function(design) {
  ep <- design$endpoint
  num <- report_num
  hr <- identical(design$effect, "hazard_ratio")
  groups <- "($k = t, c$ for the treatment and control groups)"
  model <- switch(
    ep,
    binary = paste(
      "The number of responders among the $n_k$ subjects of group $k$",
      groups, "follows a binomial distribution with response rate",
      "$\\pi_k$, and the treatment effect is the difference in response",
      "rates, $\\theta = \\pi_t - \\pi_c$."
    ),
    continuous = paste(
      "The outcomes of group $k$", groups, "follow a normal distribution",
      "$N(\\mu_k, \\sigma_k^2)$, and the treatment effect is the difference",
      "in means, $\\theta = \\mu_t - \\mu_c$."
    ),
    binary2 = paste0(
      "Each subject has a response (1) or no response (0) on each of two ",
      "binary endpoints. In group $k$ ", groups, ", the probabilities of ",
      "the responses 00, 01, 10 and 11 (the first digit for endpoint 1 and ",
      "the second for endpoint 2) are ",
      "$p_k = (p_{k00}, p_{k01}, p_{k10}, p_{k11})$, the response rates are ",
      "$\\pi_{k1} = p_{k10} + p_{k11}$ and $\\pi_{k2} = p_{k01} + p_{k11}$, ",
      "and the treatment effects are the differences in response rates, ",
      "$\\theta_1 = \\pi_{t1} - \\pi_{c1}$ and ",
      "$\\theta_2 = \\pi_{t2} - \\pi_{c2}$."
    ),
    continuous2 = paste(
      "The outcomes of group $k$", groups, "on the two endpoints follow a",
      "bivariate normal distribution $N_2(\\mu_k, \\Sigma_k)$ with means",
      "$\\mu_k = (\\mu_{k1}, \\mu_{k2})$, and the treatment effects are the",
      "differences in means, $\\theta_1 = \\mu_{t1} - \\mu_{c1}$ and",
      "$\\theta_2 = \\mu_{t2} - \\mu_{c2}$."
    ),
    tte = paste(
      "The event times of group $k$", groups, "follow an exponential",
      "distribution with hazard $\\lambda_k$, and the treatment effect is",
      if (hr) {
        "the hazard ratio, $\\theta = \\lambda_t / \\lambda_c$."
      } else {
        sprintf(paste("the difference in the survival probabilities at the",
                      "milestone time $t_0 = %s$, $\\theta = S_t(t_0) -",
                      "S_c(t_0)$ with $S_k(t_0) = \\exp(-\\lambda_k t_0)$."),
                num(design$milestone))
      }
    )
  )
  none <- "There is no concurrent control group:"
  prior_c <- if (!is.null(design$prior_c)) {
    report_prior(design$prior_c, "c", ep)
  }
  hyp <- design$hypothetical_c
  control <- switch(
    design$control,
    concurrent = sprintf("The prior distribution of the control group is %s.",
                         prior_c),
    external = paste(none, sprintf(paste(
      "the distribution of the parameters of the control group is their",
      "posterior distribution given the external control data, with the",
      "prior distribution %s."
    ), prior_c)),
    hypothetical = if (hyp$type == "data") {
      paste(none, sprintf(paste(
        "the control group is a hypothetical control, whose parameters have",
        "their posterior distribution given the hypothetical control data",
        "(%s), with the prior distribution %s."
      ), report_data(hyp$data), prior_c))
    } else if (ep == "continuous") {
      paste(none, sprintf(paste(
        "the control group is a hypothetical control, whose mean has the",
        "same t distribution as the posterior distribution of the treatment",
        "mean, with location $\\mu_{0c} = %s$ and squared scale multiplied",
        "by $r = %s$."
      ), num(hyp$mean), num(hyp$r)))
    } else {
      paste(none, sprintf(paste(
        "the control group is a hypothetical control, whose means have the",
        "same bivariate t distribution as the posterior distribution of the",
        "treatment means, with location $\\mu_{0c} = (%s)$ and scale matrix",
        "multiplied by $r = %s$."
      ), num(hyp$mean), num(hyp$r)))
    },
    fixed = paste(none, sprintf("the control %s.", switch(
      ep,
      binary = sprintf("response rate is fixed at $\\pi_c = %s$",
                       num(design$fixed_c)),
      continuous = sprintf("mean is fixed at $\\mu_c = %s$",
                           num(design$fixed_c)),
      binary2 = sprintf(
        "response rates are fixed at $(\\pi_{c1}, \\pi_{c2}) = (%s)$",
        num(design$fixed_c)
      ),
      continuous2 = sprintf(
        "means are fixed at $(\\mu_{c1}, \\mu_{c2}) = (%s)$",
        num(design$fixed_c)
      ),
      tte = if (hr) {
        sprintf(paste("median survival time is fixed at %s, so that the",
                      "control hazard is $\\lambda_c = \\log(2) / %s$"),
                num(design$fixed_c), num(design$fixed_c))
      } else {
        sprintf(paste("survival probability at the milestone time is fixed",
                      "at $S_c(t_0) = %s$"), num(design$fixed_c))
      }
    )))
  )
  external <- function(e, group) {
    sprintf(paste("External data of the %s group (%s) are borrowed through",
                  "a power prior with weight $a_0 = %s$."), group,
            report_data(e$data), num(e$a0))
  }
  priors <- c(
    sprintf("The prior distribution of the treatment group is %s.",
            report_prior(design$prior_t, "t", ep)),
    if (!is.null(design$external_t)) external(design$external_t, "treatment"),
    control,
    if (!is.null(design$external_c)) external(design$external_c, "control"),
    "The posterior distributions have closed forms."
  )
  list(report_block("heading", text = "Statistical model"),
       report_block("par", text = model),
       report_block("par", text = paste(priors, collapse = " ")))
}
