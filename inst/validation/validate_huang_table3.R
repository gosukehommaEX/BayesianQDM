# Validation against Huang, Yamaguchi, Homma, Takeda and Maruo (submitted),
# "Conjugate representation of power priors for efficient Bayesian analysis
# of normal data", Table 3 (pediatric extrapolation example). Vague initial
# priors; adult (external) data of 60 subjects per group (treatment mean
# -0.7, SD 3.8; control mean 2.5, SD 4.0) with power prior weight a0 for
# both groups; pediatric (current) data of 30 subjects per group (treatment
# mean -0.5, SD 3.9; control mean 2.8, SD 4.3); a0 = 0.25, 0.5, 0.75.
#
# The posterior distributions of the group means are t distributions, whose
# parameters are returned by qdm_decide(). The article summarizes the
# posterior distribution of the treatment effect from Monte Carlo samples
# (Section 6.3), so the numbers may differ from the exact values by Monte
# Carlo error: differences up to 0.01 beyond rounding are EXPLAINED (reason
# "huang_mc"). The credible interval is computed here from ptdiff() with
# uniroot().

validate_huang_table3 <- function() {
  study <- "Huang et al. (submitted), Table 3"
  pub <- list(
    mean_c = c(2.70, 2.65, 2.62), sd_c = c(0.640, 0.543, 0.480),
    mean_t = c(-0.57, -0.60, -0.62), sd_t = c(0.589, 0.504, 0.448),
    mean_d = c(-3.27, -3.25, -3.24), sd_d = c(0.870, 0.741, 0.656),
    prob = c(0.928, 0.954, 0.970), lower = c(-4.98, -4.71, -4.53),
    upper = c(-1.56, -1.80, -1.95)
  )
  a0 <- c(0.25, 0.5, 0.75)
  rows <- list()
  for (i in seq_along(a0)) {
    des <- design_continuous(
      30, 30, prior_vague(), prior_vague(),
      criterion_posterior(theta_tv = 0, theta_mav = -2, gamma_go = 0.5,
                          gamma_nogo = 0.95),
      external_t = borrow(data_continuous(60, -0.7, 3.8), a0[i]),
      external_c = borrow(data_continuous(60, 2.5, 4.0), a0[i])
    )
    res <- qdm_decide(des, data_continuous(30, -0.5, 3.9),
                      data_continuous(30, 2.8, 4.3))
    pt <- res$posterior_t
    pc <- res$posterior_c
    scale_t <- pt[["sigma"]] / sqrt(pt[["kappa"]])
    scale_c <- pc[["sigma"]] / sqrt(pc[["kappa"]])
    var_t <- scale_t ^ 2 * pt[["nu"]] / (pt[["nu"]] - 2)
    var_c <- scale_c ^ 2 * pc[["nu"]] / (pc[["nu"]] - 2)
    cdf <- function(q) {
      ptdiff(q, pt[["mu"]], scale_t, pt[["nu"]], pc[["mu"]], scale_c,
             pc[["nu"]])
    }
    quant <- function(p) {
      stats::uniroot(function(q) cdf(q) - p, c(-20, 20),
                     tol = 1e-10)$root
    }
    setting <- sprintf("a0 = %s", format(a0[i]))
    row <- function(quantity, published, computed, decimals) {
      number_row(study, quantity, setting, "NI", published, computed,
                 decimals, explain_tol = 0.01,
                 reason_explained = "huang_mc")
    }
    rows <- c(rows, list(
      row("Control posterior mean", pub$mean_c[i], pc[["mu"]], 2),
      row("Control posterior SD", pub$sd_c[i], sqrt(var_c), 3),
      row("Treatment posterior mean", pub$mean_t[i], pt[["mu"]], 2),
      row("Treatment posterior SD", pub$sd_t[i], sqrt(var_t), 3),
      row("Effect posterior mean", pub$mean_d[i], pt[["mu"]] - pc[["mu"]],
          2),
      row("Effect posterior SD", pub$sd_d[i], sqrt(var_t + var_c), 3),
      row("P(delta < -2 | data)", pub$prob[i], res$prob_nogo, 3),
      row("95% credible interval, lower", pub$lower[i], quant(0.025), 2),
      row("95% credible interval, upper", pub$upper[i], quant(0.975), 2)
    ))
  }
  do.call(rbind, rows)
}
