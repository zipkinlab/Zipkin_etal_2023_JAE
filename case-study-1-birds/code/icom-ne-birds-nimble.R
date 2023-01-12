# icom-ne-birds-nimble.R: NIMBLE code for running the ICOM with 
#                         eBird and BBS data for a community of interior
#                         forest obligates. 
# Author: Jeffrey W. Doser (doserjef@msu.edu)
# Citation: 

require(nimble)
icom.code <- nimbleCode ({
  # Priors ----------------------------------------------------------------
  # Process ---------------------------
  for (i in 1:p.occ) {
    beta.comm[i] ~ dnorm(0, var = 2.72)
    sigma.sq.beta[i] ~ dinvgamma(0.1, 0.1)
  } # i

  # BBS Detection ---------------------
  for (i in 1:p.det.bbs) {
    alpha.comm.bbs[i] ~ dnorm(0, var = 2.72)
    sigma.sq.bbs[i] ~ dinvgamma(0.1, 0.1)
  } # i
  # eBird Detection -------------------
  for (i in 1:p.det.eb) {
    alpha.comm.eb[i] ~ dnorm(0, var = 2.72)
    sigma.sq.eb[i] ~ dinvgamma(0.1, 0.1)
  }

  # Species-specific coefficients -----------------------------------------
  for (i in 1:N) {
    for (q in 1:p.occ) {
      beta[i, q] ~ dnorm(beta.comm[q], var = sigma.sq.beta[q])
    } # q
    for (q in 1:p.det.bbs) {
      alpha.bbs[i, q] ~ dnorm(alpha.comm.bbs[q], var = sigma.sq.bbs[q])
    } # q
    for (q in 1:p.det.eb) {
      alpha.eb[i, q] ~ dnorm(alpha.comm.eb[q], var = sigma.sq.eb[q])
    } # q
  } # i

  # Likelihood and Process Models ------------------------------------------
  # Process ---------------------------
  for (i in 1:N) {
    for (j in 1:J) {
      logit(psi[i, j]) <- inprod(beta[i, 1:p.occ], X[j, 1:p.occ])
      z[i, j] ~ dbern(psi[i, j])
    } # j
  } # i

  # BBS Likelihood --------------------
  for (i in 1:n.vals.bbs) {
    logit(p.bbs[i]) <- inprod(alpha.bbs[sp.indx.bbs[i], 1:p.det.bbs], X.bbs[i, 1:p.det.bbs])
    # 5 sections of each route
    y.bbs[i] ~ dbinom(p.bbs[i] * z[sp.indx.bbs[i], cell.bbs[i]], 5)
  } # i

  # eBird Likelihood ------------------
  for (i in 1:n.vals.eb) {
    logit(p.eb[i]) <- inprod(alpha.eb[sp.indx.eb[i], 1:p.det.eb], X.eb[i, 1:p.det.eb]) 
    y.eb[i] ~ dbern(p.eb[i] * z[sp.indx.eb[i], cell.eb[i]])
  } # i

  # Richness as a derived quantity ----------------------------------------
  for (j in 1:J) {
    z.sum[j] <- sum(z[1:N, j])
  } # j

})
