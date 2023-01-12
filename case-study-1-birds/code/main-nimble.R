# main-nimble.R: main script for calling NIMBLE to run the ICOM for a 
#                community of interior forest obligates using eBird and 
#                BBS data. Note that running this model with these data 
#                takes a substantial amount of RAM (about 13GB).
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(nimble)

# Read in data and NIMBLE code --------------------------------------------
# Bird data + forest cover + elevation
load("data/ne-data-bundle.rda")
# Bioclim variables
load("data/climate-data-full.rda")
# The NIMBLE code. 
source("code/icom-ne-birds-nimble.R")

# Data prep ---------------------------------------------------------------
# Total number of cells in the study area.
J <- nrow(grid.ne)
# Number of cells with eBird data
J.ebird <- n_distinct(ebird.df$cell)
# Number of cells with BBS data
J.bbs <- n_distinct(y.bbs$cell)
# Number of species
N <- n_distinct(ebird.df$sp)

# Filter the eBird data to use data from 3 weeks at the beginning of June
# when closure is reasonable. This also corresponds to the timing of the BBS
# surveys. 
# reasonable. 
ebird.df <- ebird.df %>%
  filter(week %in% c(23, 24, 25))

# Combine percent forest, elevation, and climate variables in a single data
# frame. 
occ.covs <- data.frame(occ.covs, climate.vars)

# Get chain number from command line run ----------------------------------
# This is used to save the chain number in the resulting file name.
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# For testing
chain <- 1
if(length(chain) == 0) base::stop('Need to tell NIMBLE the chain number')

# Only extract the 5x5km cells with eBird or BBS data for model fitting ---
# Will predict richness at non-sampled cells post-hoc to make things
# faster with NIMBLE, rather than doing it within the NIMBLE code itself. 
unique.cells <- sort(unique(c(unique(y.bbs$cell), unique(ebird.df$cell))))
# Total number of cells with data
J <- length(unique.cells)
# New occurrence covariates. 
occ.covs.small <- occ.covs[unique.cells, ]
# Adjust the cell number to make sure things match up with only the cells
# with data with them. 
for (i in 1:J) {
  y.bbs$cell[which(y.bbs$cell == unique.cells[i])] <- i
  ebird.df$cell[which(ebird.df$cell == unique.cells[i])] <- i
}

# Get Data Prepped for NIMBLE ---------------------------------------------
# Occurrence design matrix
X <- cbind(1, 
	   c(scale(occ.covs.small$elev)),
	   c(scale(occ.covs.small$elev)^2),
	   c(scale(occ.covs.small$pf)), 
           c(scale(occ.covs.small$bio1)), 
           c(scale(occ.covs.small$bio2)), 
           c(scale(occ.covs.small$bio8)),
           c(scale(occ.covs.small$bio12)),
           c(scale(occ.covs.small$bio18)))
# Set missing values to their means
X[which(is.na(X))] <- 0
# Number of occurrence parameters
p.occ <- ncol(X)
# Prep BBS Data------------------------
# Reorder by all species, then cell within species
bbs.df <- y.bbs %>%
  arrange(sp, cell)
# The actual data
y.bbs <- bbs.df$binom
# Unique species names
sp.names <- unique(bbs.df$sp)
# Species index for each data point
sp.indx.bbs <- as.numeric(factor(bbs.df$sp))
# Cell index for each data point. 
cell.bbs <- bbs.df$cell
# Observer index
obs.indx <- as.numeric(factor(bbs.df$obs))
n.obs.bbs <- n_distinct(obs.indx)
# Fixed effects detection design matrix
X.bbs <- cbind(1, 
	       c(scale(bbs.df$julian)),
	       c(scale(bbs.df$julian)^2))
# Number of BBS detection fixed effects
p.det.bbs <- ncol(X.bbs)
# Total number of BBS data points
n.vals.bbs <- length(y.bbs)
# Prep eBird Data ---------------------
# Reorder by all species, then cell within species
ebird.df <- ebird.df %>%
  arrange(sp, cell)
# The actual data
y.eb <- ebird.df$y
# Species index
sp.indx.eb <- as.numeric(factor(ebird.df$sp))
# Cell index
cell.eb <- ebird.df$cell
# Fixed effects detection design matrix
X.eb <- cbind(1, 
	      c(scale(ebird.df$day)),
	      c(scale(ebird.df$day)^2),
	      c(scale(ebird.df$time)),
	      c(scale(ebird.df$length)),
	      c(scale(ebird.df$dist)),
	      c(scale(ebird.df$obsv)))
# Number of eBird detection fixed effects
p.det.eb <- ncol(X.eb)
# Total number of eBird data points
n.vals.eb <- length(y.eb)

# Constants ---------------------------------------------------------------
icom.consts <- list(p.occ = p.occ, p.det.bbs = p.det.bbs, p.det.eb = p.det.eb,
		    N = N, J = J, 
		    sp.indx.bbs = sp.indx.bbs, obs.indx = obs.indx,
		    cell.bbs = cell.bbs, n.vals.bbs = n.vals.bbs,
		    sp.indx.eb = sp.indx.eb, 
		    cell.eb = cell.eb, n.vals.eb = n.vals.eb)
# Data --------------------------------------------------------------------
icom.data <- list(y.bbs = y.bbs, y.eb = y.eb, X.bbs = X.bbs, X.eb = X.eb, X = X)
# Initial values ----------------------------------------------------------
z.init <- matrix(1, N, J)
icom.inits <- list(z = z.init, beta.comm = rnorm(p.occ),
		   sigma.sq.beta = runif(p.occ, 0.5, 3),
		   alpha.comm.bbs = rnorm(p.det.bbs),
		   sigma.sq.bbs = runif(p.det.bbs, 0.5, 3),
		   alpha.comm.eb = rnorm(p.det.eb),
		   sigma.sq.eb = runif(p.det.eb, 0.5, 3))
# Create the model --------------------------------------------------------
icom.model <- nimbleModel(code = icom.code, name = 'icom', constants = icom.consts,
		          data = icom.data, inits = icom.inits)
# Configure MCMC ----------------------------------------------------------
icom.conf <- configureMCMC(icom.model, monitors = c('beta.comm', 'sigma.sq.beta',
						    'alpha.comm.bbs', 'sigma.sq.bbs',
						    'alpha.comm.eb', 'sigma.sq.eb',
						    'beta', 'alpha.bbs', 'alpha.eb'))
# Create an MCMC function -------------------------------------------------
icom.mcmc <- buildMCMC(icom.conf)
# Compile model
icom.c.model <- compileNimble(icom.model)
icom.c.mcmc <- compileNimble(icom.mcmc, project = icom.model)
# Number of iterations --------------------------------------------------
n.iter <- 70000
n.burn <- 40000
n.thin <- 15 
n.chain <- 1
samples <- runMCMC(icom.c.mcmc, niter = n.iter, nburnin = n.burn,
	           thin = n.thin, nchains = n.chain, samplesAsCodaMCMC = TRUE)

# Save results ------------------------------------------------------------
save(samples, file = paste("results/ne-icom-", chain, "-chain-",
		           Sys.Date(), ".R", sep = ''))
