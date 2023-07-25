# main-sim-stan.R: this script fits the integrated community model for the 
#                  simulated "butterfly" data set using stan, called through 
#                  the rstan R package.
# Author: Jeffrey W. Doser and Wendy Leuenberger
rm(list = ls())
library(rstan)
# Load data formatted for spAbundance -------------------------------------
# This loads an object called data.list
load('data/sim-spAbundance-data.rda')
# Check out its structure
str(data.list)
# data.list is comprised of two components: 
#    y: the abundance count data for the ten simulated species. This is 
#       specified in the form of a two-dimensional matrix, where the first
#       dimension corresponds to species and the second dimension corresponds
#       to "site", where here "site" is defined as a unique combination of 
#       year, week, and spatial location. 
#    covs: a list of covariates for inclusion in the community model. This
#          includes the following covariates:
#          (1) data.ind.2: indicator variable taking value 1 if data point is 
#                          from the second data set, and 0 if not
#          (2) data.ind.3: sames as (1), but for third data set.
#          (3) data.ind.4: same as (2), but for fourth data set.
#          (4) data.ind.5: sames as (3), but for fifth data set.
#          (5): year.ind: indicator for year of sampling. 
#          (6): week.ind: indicator for week of sampling.
#          (7): site.ind: site id. 
#          (8): county.ind: county id.
#          (9): survey.cov: the total number of surveys performed for the 
#                           given combination of year, week, spatial location. Note
#                           this is actually equal to number of surveys - 1, such 
#                           that a value of 0 corresponds to surveys with one 
#                           survey (aka the intercept is expected count for one survey).

# Format data in list for stan --------------------------------------------
sp.names <- dimnames(data.list$y)[[1]]
stan.list <- data.list$covs
# Using same names as those to fit the actual butterfly model.
names(stan.list) <- c('ia_ind', 'il_ind', 'mi_ind', 'oh_ind', 
		      'survey_cov', 'year_cov', 'week_cov', 'site_id', 
		      'county_id')

# Order: ordered by species, then site within species. 
stan.list$y <- c(t(data.list$y))
stan.list$n_species <- nrow(data.list$y)
stan.list$n_counties <- length(unique(data.list$covs$county.ind))
stan.list$n_counties_long <- stan.list$n_counties * stan.list$n_species
stan.list$n_years <- length(unique(data.list$covs$year.ind))
stan.list$n_years_long <- stan.list$n_years * stan.list$n_species
stan.list$n_sites <- length(unique(data.list$covs$site.ind))
stan.list$n_sites_long <- stan.list$n_sites * stan.list$n_species
stan.list$n_surveys <- length(stan.list$year_cov)
stan.list$n_surveys_long <- stan.list$n_surveys * stan.list$n_species
stan.list$year_id <- stan.list$year_cov
stan.list$year_cov <- c(scale(stan.list$year_cov))
stan.list$week_cov <- c(scale(stan.list$week_cov))
stan.list$week2_cov <- stan.list$week_cov^2 

# Take a look at the stan data object
str(stan.list)

# Run the model in stan ---------------------------------------------------
# MCMC parameters
n.samples <- 3000
n.burn <- 2000
n.thin <- 5
n.chains <- 3
n.cores <- 3

# Initial values
set.seed(1234)
inits <- lapply(1:n.chains, function(i)
  list(muBetaInt = rnorm(1),
       muBetaIA = rnorm(1),
       muBetaIL = rnorm(1),
       muBetaMI = rnorm(1),
       muBetaOH = rnorm(1),
       muBetaYear = rnorm(1),
       muBetaWeek = rnorm(1),
       muBetaWeek2 = rnorm(1),
       muBetaSurvey = rnorm(1),
       tauSqInt = runif(1, 0, 1),
       tauSqIA = runif(1, 0, 1),
       tauSqIL = runif(1, 0, 1),
       tauSqMI = runif(1, 0, 1),
       tauSqOH = runif(1, 0, 1),
       tauSqYear = runif(1, 0, 1),
       tauSqWeek = runif(1, 0, 1),
       tauSqWeek2 = runif(1, 0, 1),
       tauSqSurvey = runif(1, 0, 1),
       sigma2_county = runif(1, 0, 1),
       sigma2_week = runif(1, 0, 1),
       sigma2_week2 = runif(1, 0, 1),
       sigma2_site = runif(1, 0, 1),
       sigma2_year = runif(1, 0, 1),
       r_count=runif(stan.list$n_species,0,2)
  )
)

# Parameters to monitor
params <- c('betaInt', 'betaIA', 'betaIL', 'betaMI', 
	    'betaOH', 'betaYear', 'betaWeek', 'betaWeek2', 
	    'betaSurvey', 'sigma2_county', 'sigma2_week', 
	    'sigma2_week2', 'sigma2_site', 'sigma2_year', 
	    'betaStarCounty', 'betaStarWeek', 'betaStarWeek2', 
	    'betaStarSite', 'betaStarYear', 'muBetaInt', 
            'muBetaIA', 'muBetaIL', 'muBetaMI', 'muBetaOH', 
            'muBetaYear', 'muBetaWeek', 'muBetaWeek2', 
            'muBetaSurvey', 'tauSqInt', 'tauSqIA', 'tauSqIL', 
            'tauSqMI', 'tauSqOH', 'tauSqYear', 'tauSqWeek', 
            'tauSqWeek2', 'tauSqSurvey')

# Call stan from R --------------------------------------------------------
stan.model <- rstan::stan_model('code/stan-code/icm.stan')
# Note that 1000 iterations takes a little over 24 hours.
out <- stan('code/stan-code/icm.stan',
            control = list(adapt_delta = 0.8),
            data = stan.list, init = inits, pars = params,
            chains = n.chains, iter = n.samples, warmup = n.burn, thin = n.thin,
            seed = 1, cores = n.cores, open_progress = TRUE)
# Save to hard drive
save(out, file = 'results/stan-icm-results.rda')
