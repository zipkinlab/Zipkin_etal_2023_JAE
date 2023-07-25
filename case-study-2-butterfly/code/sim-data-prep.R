# sim-data-prep.R: this script generates a simulated data set reminiscent of the 
#                  butterfly case study for testing the ICM framework in this 
#                  scenario. Two of the 4 butterfly data sets are proprietary, 
#                  and thus the full butterfly model fit in the case study cannot
#                  be done.
# Author: Jeffrey W. Doser
rm(list = ls())

# Simulted data set description -------------------------------------------
# Data will be simulated from four monitoring programs, analogous to the case
# study. The data sets will be comprised of 30, 20, 35, 25, and 15 locations, 
# respectively, for a total of 125 locations within 30 counties. Each site will be sampled during
# 4 weeks within each of 5 years. Abundance for each species will be generated
# according to the equations in Supplemental Information S2.

# Generate a data set of 10 species observed at 200 sites, sampled over 10 years,
# during 10 weeks in each year. 

# Generate covariates and locations ---------------------------------------
# Set seed to get the same data set
set.seed(18273)
# Number of species 
n.sp <- 10
# Number of sites
J <- 125
# Data locations in each data set
sites.1 <- 1:30
sites.2 <- 31:50
sites.3 <- 51:85
sites.4 <- 86:110
sites.5 <- 111:125
# Number of years
n.years <- 5
# Number of weeks
n.weeks <- 4
# Total number of data points for each species
n.data <- J * n.years * n.weeks
# Indicator variable for site
site.ind <- rep(1:J, each = n.weeks * n.years)
# Indicator variable for week
week.ind <- rep(1:n.weeks, times = J * n.years)
# Indicator variable for year
year.ind <- rep(rep(1:n.years, each = n.weeks), times = J)
# Indicator variables for each data type
data.ind.1 <- rep(0, n.data)
data.ind.1[which(site.ind %in% sites.1)] <- 1
data.ind.2 <- rep(0, n.data)
data.ind.2[which(site.ind %in% sites.2)] <- 1
data.ind.3 <- rep(0, n.data)
data.ind.3[which(site.ind %in% sites.3)] <- 1
data.ind.4 <- rep(0, n.data)
data.ind.4[which(site.ind %in% sites.4)] <- 1
data.ind.5 <- rep(0, n.data)
data.ind.5[which(site.ind %in% sites.5)] <- 1
# Generate covariate indicating the number of surveys done in each week
props.surveys <- c(0.7, 0.2, 0.05, 0.05)
survey.cov <- rep(0, n.data)
tmp.indx.1 <- sample(1:n.data, props.surveys[2] * n.data, replace = FALSE)
survey.cov[tmp.indx.1] <- 1
tmp.indx.2 <- sample(which(survey.cov != 1), props.surveys[3] * n.data, replace = FALSE)
survey.cov[tmp.indx.2] <- 2
tmp.indx.3 <- sample(which(survey.cov == 0), props.surveys[4] * n.data, replace = FALSE)
survey.cov[tmp.indx.3] <- 3
# Get correct values
survey.cov <- survey.cov + 1
# County covariate
n.counties <- 30
county.cov <- sample(1:n.counties, n.data, replace = TRUE)
# Set parameter values ----------------------------------------------------
# Community-level parameters ----------
# Intercept (first data set average)
beta.comm.0 <- 0
tau.sq.beta.0 <- 0.5
# Data set 2 difference from first
beta.comm.1 <- -1
tau.sq.beta.1 <- 0.5
# Data set 3 difference from first
beta.comm.2 <- 0.5
tau.sq.beta.2 <- 0.6
# Data set 4 difference from first
beta.comm.3 <- 0.6
tau.sq.beta.3 <- 0.3
# Data set 5 difference from first
beta.comm.4 <- -2
tau.sq.beta.4 <- 0.6
# Average effect of survey
beta.comm.5 <- 0.75
tau.sq.beta.5 <- 0.05
# Average trend
beta.comm.6 <- -0.3
tau.sq.beta.6 <- 0.5
# Average linear week effect
beta.comm.7 <- 0.25
tau.sq.beta.7 <- 0.5
# Average quadratic week effect
beta.comm.8 <- -0.4
tau.sq.beta.8 <- 0.3
# Variance for random linear slope of week by year
sigma.sq.week.year <- 0.2
# Variance for random quadratic effect of week by year
sigma.sq.week.2.year <- 0.3
# Variance for random county level effect
sigma.sq.county <- 0.6
# Variance for random site level effect
sigma.sq.site <- 0.4
# Variance for random year effect
sigma.sq.year <- 0.1
# Species-level parameters ------------
# Intercept
beta.0 <- rnorm(n.sp, beta.comm.0, sqrt(tau.sq.beta.0))
# Data set 2 diff from first
beta.1 <- rnorm(n.sp, beta.comm.1, sqrt(tau.sq.beta.1))
# Data set 3 diff from first
beta.2 <- rnorm(n.sp, beta.comm.2, sqrt(tau.sq.beta.2))
# Data set 4 diff from first
beta.3 <- rnorm(n.sp, beta.comm.3, sqrt(tau.sq.beta.3))
# Data set 5 diff from first
beta.4 <- rnorm(n.sp, beta.comm.4, sqrt(tau.sq.beta.4))
# Average effect of survey
beta.5 <- rnorm(n.sp, beta.comm.5, sqrt(tau.sq.beta.5))
# Average trend
beta.6 <- rnorm(n.sp, beta.comm.6, sqrt(tau.sq.beta.6))
# Average linear week effect
beta.7 <- rnorm(n.sp, beta.comm.7, sqrt(tau.sq.beta.7))
# Average quadratic week effect
beta.8 <- rnorm(n.sp, beta.comm.8, sqrt(tau.sq.beta.8))
# Species/year specific week linear effects
beta.star.week.linear <- matrix(rnorm(n.sp * n.years, 0, sqrt(sigma.sq.week.year)), 
				n.sp, n.years)
# Species/year specific week quadratic effects
beta.star.week.quad <- matrix(rnorm(n.sp * n.years, 0, sqrt(sigma.sq.week.2.year)), 
				n.sp, n.years)
# Species/county level random effect
beta.star.county <- matrix(rnorm(n.counties * n.sp, 0, sqrt(sigma.sq.county)), n.sp, n.counties)
# Species/site level random effect
beta.star.site <- matrix(rnorm(J * n.sp, 0, sqrt(sigma.sq.site)), n.sp, J)
# Species/year level random effect
beta.star.year <- matrix(rnorm(n.years * n.sp, 0, sqrt(sigma.sq.year)), n.sp, n.years)
# Species-specific dispersion parameters (lower values = more overdispersion)
kappa <- runif(n.sp, 0.01, 5)
# Generate the observed data ----------------------------------------------
y <- array(NA, dim = c(n.sp, n.data))
mu <- array(NA, dim = c(n.sp, n.data))
year.s <- c(scale(year.ind))
week.s <- c(scale(week.ind))
# Fill in the array
for (j in 1:n.data) {
  for (i in 1:n.sp) {
    mu[i, j] <- exp(beta.0[i] + beta.1[i] * data.ind.2[j] + beta.2[i] * data.ind.3[j] + 
		    beta.3[i] * data.ind.4[j] + beta.4[i] * data.ind.5[j] + 
		    beta.5[i] * log(survey.cov[j]) + beta.6[i] * year.s[j] + 
		    beta.7[i] * week.s[j] + beta.8[i] * week.s[j]^2 + 
		    beta.star.week.linear[i, year.ind[j]] * week.s[j] + 
		    beta.star.week.quad[i, year.ind[j]] * week.s[j]^2 + 
		    beta.star.site[i, site.ind[j]] + beta.star.county[i, county.cov[j]] + 
		    beta.star.year[i, year.ind[j]])
    y[i, j] <- rnbinom(1, size = kappa[i], mu = mu[i, j])
  }
}	

# Format data for spAbundance and save to hard drive ----------------------
data.list <- list(y = y, 
		  covs = list(data.ind.2 = data.ind.2, 
			      data.ind.3 = data.ind.3, 
			      data.ind.4 = data.ind.4, 
			      data.ind.5 = data.ind.5, 
			      survey.cov = log(survey.cov),
			      year.ind = year.ind, 
			      week.ind = week.ind,
			      site.ind = site.ind, 
			      county.ind = county.cov))
beta.true <- cbind(beta.0, beta.1, beta.2, beta.3, beta.4, beta.5, beta.6, beta.7, beta.8)
save(data.list, beta.true, file = 'data/sim-spAbundance-data.rda')
