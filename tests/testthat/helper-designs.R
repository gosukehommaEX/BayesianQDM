# Designs shared by the tests of the binary endpoint. D1 to D7 are the
# designs of Python_reference/ref_binary_qdm.py, whose exact results
# (rational arithmetic) are used as references; D8 and D9 (fixed control) are
# those of Python_reference/ref_fixed_tte_qdm.py.
crit_post <- function(gamma_go = 0.5, gamma_nogo = 0.8, tv = 0.2, mav = 0.1) {
  criterion_posterior(theta_tv = tv, theta_mav = mav, gamma_go = gamma_go,
                      gamma_nogo = gamma_nogo)
}

ref_design <- function(name) {
  unif <- prior_beta(1, 1)
  switch(
    name,
    D1 = design_binary(12, 12, unif, unif, crit_post()),
    D2 = design_binary(10, 8, unif, prior_beta(2, 3), crit_post(),
                       external_t = borrow(data_binary(10, 6), 0.5),
                       external_c = borrow(data_binary(20, 8), 0.5)),
    D3 = design_binary(12, 0, unif, unif, crit_post(),
                       hypothetical_c = hypothetical(data_binary(12, 4))),
    D4 = design_binary(12, 0, unif, unif, crit_post(),
                       external_c = borrow(data_binary(20, 6), 0.5)),
    D5 = design_binary(10, 10, unif, unif,
                       criterion_predictive(0.05, m_t = 30, m_c = 20,
                                            gamma_go = 0.9,
                                            gamma_nogo = 0.7)),
    D6 = design_binary(8, 8, unif, unif, crit_post(0.2, 0.3)),
    D7 = design_binary(6, 6, unif, unif, crit_post(tv = 0, mav = -0.1)),
    D8 = design_binary(20, 0, unif,
                       criterion = crit_post(0.8, 0.8, tv = 0.2, mav = 0.05),
                       fixed_c = 0.2),
    D9 = design_binary(10, 0, unif,
                       criterion = criterion_predictive(0, m_t = 10,
                                                        gamma_go = 0.6,
                                                        gamma_nogo = 0.6),
                       fixed_c = 0.3)
  )
}

scen_c <- data.frame(pi_t = c(0.2, 0.4, 0.6, 0.8), pi_c = 0.4)
scen_u <- data.frame(pi_t = c(0.2, 0.4, 0.6, 0.8))
