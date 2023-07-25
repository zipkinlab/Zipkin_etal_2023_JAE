# main-sim-spAbundance: this script fits the integrated community model for the 
#                       simulated butterfly community using the spAbundance R 
#                       package. The data were generated similarly to the structure
#                       of the true butterfly data, and this script allows the 
#                       full model structure to be fit without any of the data
#                       sharing restrictions from the butterfly case study
# Author: Jeffrey W. Doser
rm(list = ls())
library(spAbundance)

# Load data formatted for spAbundance -------------------------------------
# This loads an object called data.list
load('data/sim-spAbundance-data.rda')
# Check out its structure
str(data.list)
# data.list is comprised of two components: 
#    y: the abundance count data for ten simulated species. This is 
#       specified in the form of a two-dimensional matrix, where the first
#       dimension corresponds to species and the second dimension corresponds
#       to "site", where here "site" is defined as a unique combination of 
#       year, week, and spatial location. 
#    covs: a list of covariates for inclusion in the community model. This
#          includes the following covariates:
#          (1) data.ind.2: indicator variable taking value 1 if data point is from 
#                          the second simulated data set, and 0 if not. 
#          (2): data.ind.3: same as (1), but for a third data set. 
#          (3): data.ind.4: same as (1) but for a fourth data set.
#          (4): data.ind.5: same as (1) but for a fifth data set. 
#          (5): year.ind: indicator for year of sampling. 
#          (6): week.ind: indicator for week of sampling.
#          (7): site.ind: site id. 
#          (8): county.ind: county id.
#          (9): survey.cov: the total number of surveys performed for the 
#                           given combination of year, week, spatial location. Note
#                           this is actually equal to number of surveys - 1, such 
#                           that a value of 0 corresponds to surveys with one 
#                           survey (aka the intercept is expected count for one survey).

# Get chain number from command line run ----------------------------------
# This is used to save the chain number in the resulting file name.
# This is an alternative approach to running spAbundance models with multiple
# chains, which makes it possible to run chains in parallel.
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# For testing
# chain <- 1
if(length(chain) == 0) base::stop('Need to tell spAbundance the chain number')

# Specify inputs for spAbundance ------------------------------------------
# Priors
prior.list <- list(beta.comm.normal = list(mean = 0, var = 100),
                   kappa.unif = list(a = 0, b = 10), 
                   tau.sq.beta.ig = list(a = .1, b = .1)) 
# Starting values
inits.list <- list(beta.comm = 0, 
		   beta = 0, 
		   tau.sq.beta = 0.5)
# Tuning values. 
tuning.list <- list(beta = 0.05, beta.star = 0.1, kappa = 0.2)

# Specify number of iterations. Can change as desired for testing out. 

# Largest
n.batch <- 1600
batch.length <- 25
n.burn <- 20000
n.thin <- 20
n.chains <- 1

# Run the model in spAbundance --------------------------------------------
# Approximate run time for 40,000 iterations and 1 chain (chains can be run
# in parallel outside of spAbundance as is done in this script): 16.25 hours.
out <- msAbund(formula = ~ data.ind.2 + data.ind.3 + data.ind.4 + data.ind.5 + survey.cov + 
                           scale(year.ind) + scale(week.ind) + I(scale(week.ind)^2) + 
			   (scale(week.ind) | year.ind) + (I(scale(week.ind)^2) | year.ind) +
                           (1 | county.ind) + (1 | site.ind) + (1 | year.ind),
	       data = data.list,
	       n.batch = n.batch,
	       inits = inits.list,
	       priors = prior.list,
	       tuning = tuning.list,
	       family = 'NB',
	       batch.length = batch.length,
	       n.omp.threads = 1,
	       verbose = TRUE,
	       n.report = 1,
	       n.burn = n.burn,
	       n.thin = n.thin,
	       n.chains = n.chains)

# Save results ------------------------------------------------------------
# Full model results file. Note that this is a large object that is not provided 
# on GitHub since it contains proprietary data and it exceeds the maximum
# file size on GitHub. 
save(out, file = paste("results/sim-msAbund-", batch.length * n.batch, 
		       "-samples-", chain, "-chain-", Sys.Date(), ".R", sep = ''))
# Save specific posterior samples, which aren't too big for GitHub.
# Community-level coefficients --------
beta.comm.samples <- out$beta.comm.samples
save(beta.comm.samples, file = paste('results/sim-beta-comm-samples-', chain, '-chain.rda', 
				     sep = ''))
# Species-level coefficients ----------
beta.samples <- out$beta.samples
save(beta.samples, file = paste('results/sim-beta-samples-', chain, '-chain.rda', 
				     sep = ''))
# Quantiles of expected abundance -----
mu.quants <- apply(out$mu.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
save(mu.quants, file = paste('results/sim-mu-quants-', chain, '-chain.rda', sep = ''))
# Quantiles of predicted abundance ----
y.rep.quants <- apply(out$y.rep.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
save(y.rep.quants, file = paste('results/sim-y.rep-quants-', chain, '-chain.rda', sep = ''))
# Community-level variances -----------
tau.sq.beta.samples <- out$tau.sq.beta.samples
save(tau.sq.beta.samples, file = paste('results/sim-tau-sq-beta-samples-', chain, '-chain.rda', 
				     sep = ''))
# Variances of random effects ---------
sigma.sq.mu.samples <- out$sigma.sq.mu.samples
save(sigma.sq.mu.samples, file = paste('results/sim-sigma-sq-mu-samples-', chain, '-chain.rda', 
				     sep = ''))
# NB overdispersion parameters --------
kappa.samples <- out$kappa.samples
save(kappa.samples, file = paste('results/sim-kappa-samples-', chain, '-chain.rda', 
				     sep = ''))
# Species-level random effects --------
beta.star.samples <- out$beta.star.samples
save(beta.star.samples, file = paste('results/sim-beta-star-samples-', chain, '-chain.rda', 
				     sep = ''))
