# Designs shared by the tests of two continuous endpoints. They are the
# designs K1 to K6 of Python_reference/ref_continuous2_qdm.py (closed-form
# posteriors, moment matching with adaptive quadrature, Monte Carlo
# simulation with 2e7 draws and simulation of 1e5 studies per scenario).
cov2 <- function(s1, s2, r) {
  matrix(c(s1 ^ 2, r * s1 * s2, r * s1 * s2, s2 ^ 2), 2L, 2L)
}

lambda_cs2 <- matrix(c(80 ^ 2, 15, 15, 6 ^ 2), 2L, 2L)
lambda_cs4 <- matrix(c(25, 1, 1, 9), 2L, 2L)

crit_cont2 <- function(gamma_go = 0.5, gamma_nogo = 0.3, ...) {
  criterion_posterior(theta_tv = c(10, 2), theta_mav = c(5, 1),
                      gamma_go = gamma_go, gamma_nogo = gamma_nogo, ...)
}

ref_design_cont2 <- function(name, method = "MM", nMC = 10000L) {
  vague <- prior_vague()
  switch(
    name,
    # Concurrent control, vague priors
    K1 = design_continuous2(12, 12, vague, vague, crit_cont2(),
                            method = method, nMC = nMC),
    # Case study 2 of Yamaguchi et al. (submitted)
    K2 = design_continuous2(30, 30, prior_niw(c(0, 1), 0.1, 4, lambda_cs2),
                            prior_niw(c(-10, -1), 0.1, 4, lambda_cs2),
                            crit_cont2(0.296, 0.244), method = method,
                            nMC = nMC),
    # NIW priors and external data for both groups, predictive criterion
    K3 = design_continuous2(
      12, 12, prior_niw(c(1.984, 1.191), 0.1, 4, lambda_cs4),
      prior_niw(c(-0.5, -1.5), 0.1, 4, lambda_cs4),
      criterion_predictive(c(1.984, 1.191), m_t = 50, m_c = 50,
                           gamma_go = 0.6, gamma_nogo = 0.5),
      external_t = borrow(data_continuous2(8, c(2, 1),
                                           matrix(c(25, 3, 3, 9), 2L)), 0.5),
      external_c = borrow(data_continuous2(20, c(-0.5, -1.5),
                                           matrix(c(24, 2, 2, 10), 2L)), 0.4),
      method = method, nMC = nMC
    ),
    # Hypothetical control, vague treatment prior
    K4 = design_continuous2(15, 0, vague,
                            criterion = crit_cont2(0.4, 0.3, go_regions = 1:2),
                            hypothetical_c = hypothetical(mean = c(-10, -1),
                                                          r = 1.5),
                            method = method, nMC = nMC),
    # External control only, vague control prior, marginal rule "and"
    K5 = design_continuous2(
      12, 0, prior_niw(c(0, 1), 0.1, 4, lambda_cs2), vague,
      crit_cont2(c(0.6, 0.5), c(0.5, 0.6), rule = "and"),
      external_c = borrow(data_continuous2(20, c(-10, -1), cov2(70, 5, 0.3)),
                          0.5),
      method = method, nMC = nMC
    ),
    # Predictive criterion with a hypothetical control
    K6 = design_continuous2(15, 0, vague,
                            criterion = criterion_predictive(
                              c(5, 1), m_t = 40, m_c = 30, gamma_go = 0.7,
                              gamma_nogo = 0.6
                            ),
                            hypothetical_c = hypothetical(mean = c(0, 0),
                                                          r = 2),
                            method = method, nMC = nMC),
    # Fixed control (Python_reference/ref_fixed_tte_qdm.py), posterior
    # criterion
    KF1 = design_continuous2(10, 0, vague, criterion = crit_cont2(),
                             fixed_c = c(-10, -1), method = method,
                             nMC = nMC),
    # Fixed control, predictive criterion, 2 degrees of freedom
    KF2 = design_continuous2(4, 0, vague,
                             criterion = criterion_predictive(
                               c(0, 0), m_t = 20, gamma_go = 0.5,
                               gamma_nogo = 0.2
                             ),
                             fixed_c = c(1, 0.5), method = method, nMC = nMC)
  )
}

# Observed data of the reference cases
ref_data_cont2 <- list(
  K1 = list(t = data_continuous2(12, c(18, 2.5), cov2(60, 5, 0.4)),
            c = data_continuous2(12, c(-10, -1), cov2(70, 6, 0.2))),
  K2 = list(t = data_continuous2(30, c(5, 2), cov2(80, 6, 0.4)),
            c = data_continuous2(30, c(-10, -1), cov2(75, 6.5, 0.3))),
  K3 = list(t = data_continuous2(12, c(2.5, 1), cov2(5, 3, 0.3)),
            c = data_continuous2(12, c(-0.8, -1.2), cov2(4.5, 3.2, 0.1))),
  K4 = list(t = data_continuous2(15, c(4, 1.5), cov2(30, 3, 0.5)), c = NULL),
  K5 = list(t = data_continuous2(12, c(3, 1), cov2(75, 5.5, 0.4)), c = NULL),
  K6 = list(t = data_continuous2(15, c(7, 2), cov2(10, 2, -0.3)), c = NULL),
  KF1 = list(t = data_continuous2(10, c(5, 2),
                                  matrix(c(6400, 190, 190, 36), 2L)),
             c = NULL),
  KF2 = list(t = data_continuous2(4, c(2, 1), matrix(c(4, 1, 1, 2), 2L)),
             c = NULL)
)

# Scenarios of the operating characteristics (SCEN of the Python reference)
scen_cont2 <- function(name) {
  conc <- data.frame(mu_t1 = c(-10, 0, 10), mu_t2 = c(-1, 1, 2), sd_t1 = 70,
                     sd_t2 = 6, rho_t = c(0, 0.4, 0.8), mu_c1 = -10,
                     mu_c2 = -1, sd_c1 = 70, sd_c2 = 6, rho_c = c(0, 0.4, 0.8))
  switch(
    name,
    K1 = conc, K2 = conc,
    K3 = data.frame(mu_t1 = c(-0.5, 2, 4), mu_t2 = c(-1.5, 0.5, 1.5),
                    sd_t1 = 5, sd_t2 = 3, rho_t = 0.3, mu_c1 = -0.5,
                    mu_c2 = -1.5, sd_c1 = 5, sd_c2 = 3, rho_c = 0.3),
    K4 = data.frame(mu_t1 = c(-10, 0, 10), mu_t2 = c(-1, 1, 3), sd_t1 = 30,
                    sd_t2 = 3, rho_t = 0.5),
    K5 = data.frame(mu_t1 = c(-10, 5, 20), mu_t2 = c(-1, 1.5, 3),
                    sd_t1 = 75, sd_t2 = 5.5, rho_t = 0.4),
    K6 = data.frame(mu_t1 = c(0, 5, 10), mu_t2 = c(0, 1, 2), sd_t1 = 10,
                    sd_t2 = 2, rho_t = -0.3)
  )
}
