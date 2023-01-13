# main-spOccupancy.R: this is the main script for fitting an integrated
#                     community occupancy model using the spOccupancy
#                     R package. Note that fitting this model in spOccupancy
#                     takes a decent amount of RAM space (~10GB or so).
# Author: Jeffrey W. Doser

rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(devtools)
library(spOccupancy)

# Get chain number from command line run ----------------------------------
# This is used to save the chain number in the resulting file name.
# This is an alternative approach to running spOccupancy models with multiple
# chains, which makes it possible to run chains in parallel. 
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# For testing
# chain <- 1
if(length(chain) == 0) base::stop('Need to tell spOccupancy the chain number')

# Read in data formatted for spOccupancy ----------------------------------
load('data/spOccupancy-data.rda')

# Run model in spOccupancy ------------------------------------------------
# Running multiple chains across cores instead of within spOccupancy
# itself using the "n.chains" argument, as this approach makes it possible
# to run the chains in parallel. 
n.samples <- 10000
n.report <- 10
n.burn <- 6000
n.thin <- 8
n.chains <- 1
# Using default prior and initial values
out <- intMsPGOcc(occ.formula = ~ scale(elev) + I(scale(elev)^2) + scale(pf) + 
                                  scale(bio1) + scale(bio2) + scale(bio8) + 
				  scale(bio12) + scale(bio18), 
		  det.formula = list(bbs = ~ scale(julian) + I(scale(julian)^2), 
				     ebird = ~ scale(day) + I(scale(day)^2) + 
					       scale(time) + scale(length) + 
					       scale(dist) + scale(obsv)), 
		  data = data.list,
                  n.samples = n.samples, 
                  n.report = 10, 
                  n.burn = n.burn,
                  n.thin = n.thin, 
                  n.chains = n.chains)
# Save results ------------------------------------------------------------
save(out, file = paste("results/spOccupancy-icom-", chain, "-chain-",
		           Sys.Date(), ".R", sep = ''))
