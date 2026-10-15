# Designs shared by the tests of two binary endpoints. B1 to B4 are the
# designs of Python_reference/ref_binary2_qdm.py (exact references: rational
# arithmetic for predictive probabilities, Beta differences for criteria on
# one endpoint at a time).
dir_q <- function() prior_dirichlet(rep(0.25, 4))

ref_design_bin2 <- function(name, nMC = 20000L, enumerate = TRUE) {
  switch(
    name,
    B1 = design_binary2(3, 2, dir_q(), dir_q(),
                        criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3,
                                             gamma_go = 0.4,
                                             gamma_nogo = 0.3),
                        nMC = nMC, enumerate = enumerate),
    B2 = design_binary2(3, 2, dir_q(), dir_q(),
                        criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3,
                                             gamma_go = c(0.6, 0.5),
                                             gamma_nogo = c(0.5, 0.6),
                                             rule = "and"),
                        nMC = nMC, enumerate = enumerate),
    B3 = design_binary2(4, 4, dir_q(), dir_q(),
                        criterion_posterior(c(0.15, 0.15), c(0.1, 0.1),
                                            gamma_go = 0.7, gamma_nogo = 0.7,
                                            go_regions = 1:3,
                                            nogo_regions = 7:9),
                        nMC = nMC, enumerate = enumerate),
    B4 = design_binary2(4, 4, dir_q(), dir_q(),
                        criterion_posterior(c(0.15, 0.15), c(0.1, 0.1),
                                            gamma_go = c(0.5, 0.4),
                                            gamma_nogo = c(0.6, 0.5),
                                            rule = "and"),
                        nMC = nMC, enumerate = enumerate),
    # Case study 1 of Yamaguchi et al. (submitted), Table 3, rho = 0
    B5 = design_binary2(10, 10, dir_q(), dir_q(),
                        criterion_posterior(c(0.15, 0.15), c(0.10, 0.10),
                                            gamma_go = 0.263,
                                            gamma_nogo = 0.347),
                        nMC = nMC, enumerate = enumerate),
    # Case study 3 (Supplemental Section S1.1): external data, regions
    # G = {1, 2} and N = {9}
    B6 = design_binary2(8, 8, dir_q(), dir_q(),
                        criterion_posterior(c(0.15, 0.25), c(0.10, 0.15),
                                            gamma_go = 0.453,
                                            gamma_nogo = 0.229,
                                            go_regions = 1:2),
                        external_t = borrow(data_binary2(c(2, 1, 1, 2)), 0.7),
                        external_c = borrow(data_binary2(c(14, 2, 3, 1)),
                                            0.5),
                        nMC = nMC, enumerate = enumerate),
    # Hypothetical control and external control only
    B7 = design_binary2(3, 0, dir_q(), dir_q(),
                        criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3,
                                             gamma_go = 0.4,
                                             gamma_nogo = 0.3),
                        hypothetical_c = hypothetical(data_binary2(c(2, 1, 1,
                                                                     0))),
                        nMC = nMC, enumerate = enumerate),
    B8 = design_binary2(3, 0, dir_q(), dir_q(),
                        criterion_predictive(c(0.1, 0), m_t = 4, m_c = 3,
                                             gamma_go = 0.4,
                                             gamma_nogo = 0.3),
                        external_c = borrow(data_binary2(c(14, 2, 3, 1)),
                                            0.5),
                        nMC = nMC, enumerate = enumerate),
    # Fixed control (Python_reference/ref_fixed_tte_qdm.py)
    BF1 = design_binary2(3, 0, dir_q(),
                         criterion = criterion_predictive(c(0.1, 0), m_t = 4,
                                                          gamma_go = 0.4,
                                                          gamma_nogo = 0.3),
                         fixed_c = c(0.25, 0.5), nMC = nMC,
                         enumerate = enumerate),
    BF2 = design_binary2(10, 0, dir_q(),
                         criterion = criterion_posterior(c(0.25, 0.2),
                                                         c(0.05, 0),
                                                         gamma_go = 0.5,
                                                         gamma_nogo = 0.3),
                         fixed_c = c(0.35, 0.3), nMC = nMC,
                         enumerate = enumerate)
  )
}

scen_bin2 <- data.frame(pi_t1 = c(0.3, 0.6, 0.7), pi_t2 = c(0.4, 0.5, 0.8),
                        rho_t = c(0.2, 0, 0.5), pi_c1 = c(0.3, 0.3, 0.2),
                        pi_c2 = c(0.3, 0.3, 0.3), rho_c = c(0, 0, -0.2))
scen_bin2_d <- data.frame(pi_t1 = c(0.3, 0.6, 0.8), pi_t2 = c(0.3, 0.5, 0.7),
                          rho_t = c(0, 0.3, 0.2), pi_c1 = c(0.3, 0.3, 0.2),
                          pi_c2 = c(0.3, 0.3, 0.3), rho_c = c(0, 0.3, 0.1))
