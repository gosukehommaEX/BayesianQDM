# nMC:    A positive integer representing the number of Monte Carlo iterations for
#         simulation. Typical values range from 10,000 (for quick estimates) to 100,000+
#         (for high precision). Larger values yield more accurate results but require
#         more computational time.
# U1:     Upper threshold for endpoint1
# L1:     Lower threshold for endpoint1
# U2:     Upper threshold for endpoint2
# L2:     Lower threshold for endpoint2
# alpha1: Vector of shape parameters for group 1
# alpha2: Vector of shape parameters for group 2
#
# Example usage
# pMC2binaryendpoints(
#   nMC <- 1e+4,
#   U1 <- 0.2,
#   L1 <- 0.1,
#   U2 <- 0.2,
#   L2 <- 0.1,
#   alpha1 <- c(0.25, 2.25, 2.25, 1.25),
#   alpha2 <- c(1.25, 2.25, 2.25, 1.25)
# )
pMC2binaryendpoints <- function(nMC, U1, L1, U2, L2, alpha1, alpha2) {

  # Validate thredholds
  if(U1 < L1) stop("U1 must be larger than L1")
  if(U2 < L2) stop("U2 must be larger than L2")

  # Draw Dirichlet random numbers for treatment group
  p1 <- rdirichlet(nMC, alpha1)
  # π_t1 = p10 + p11  => columns 3 + 4
  pi11 <- rowSums(p1[, c(3, 4)])
  # π_t2 = p01 + p11  => columns 2 + 4
  pi12 <- rowSums(p1[, c(2, 4)])

  # Draw Dirichlet random numbers for control group
  p2 <- rdirichlet(nMC, alpha2)
  # π_c1 = p10 + p11  => columns 3 + 4
  pi21 <- rowSums(p2[, c(3, 4)])
  # π_c2 = p01 + p11  => columns 2 + 4
  pi22 <- rowSums(p2[, c(2, 4)])

  # Treatment effects
  theta1 <- pi11 - pi21
  theta2 <- pi12 - pi22

  # Region index 1,...,9
  R1 <- 1 * (theta1 > U1) + 2 * ((U1 >= theta1) & (theta1 > L1)) + 3 * (L1 >= theta1)
  R2 <- 1 * (theta2 > U2) + 2 * ((U2 >= theta2) & (theta2 > L2)) + 3 * (L2 >= theta2)
  R <- (R1 - 1) * 3 + R2

  # Region probabilities Pr(R_{l}|D))
  Pr_R <- tabulate(R, nbins = 9) / nMC
  names(Pr_R) <- paste0("R", 1:9)

  return(Pr_R)
}
