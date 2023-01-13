# calculate-bioclim.R: this script extracts the five bioclim variables
#                      from the downloaded climate data from PRISM for
#                      use as covariates in an integrated community 
#                      occupancy model. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(dismo)

# Read in climate data ----------------------------------------------------
# Load 2016 and 2017 climate data
# Need to extract June 2016 - May 2017 data for calculation of bioclim 
# variables following Rushing et al. (2019) and Clement et al. (2016)
load("data/ppt-2017.rda")
ppt.2017 <- ppt.vals
load("data/ppt-2016.rda")
ppt.2016 <- ppt.vals
# Final ppt values for time period of interest
ppt.vals <- cbind(ppt.2017[, 6:12], ppt.2016[, 1:5])
load("data/tmax-2017.rda")
tmax.2017 <- tmax.vals
load("data/tmax-2016.rda")
tmax.2016 <- tmax.vals
# Final tmax values for time period of interest
tmax.vals <- cbind(tmax.2017[, 6:12], tmax.2016[, 1:5])
load("data/tmin-2017.rda")
tmin.2017 <- tmin.vals
load("data/tmin-2016.rda")
tmin.2016 <- tmin.vals
# Final tmin values for time period of interest
tmin.vals <- cbind(tmin.2017[, 6:12], tmin.2016[, 1:5])
# Calculate bioclim variables
bioclim.vars <- biovars(ppt.vals, tmin.vals, tmax.vals)
# Extract 5 bioclim variables that have low correlation and are good for 
# modeling species ranges. 
# bio1: mean annual temperature
# bio2: mean diurnal temperature range
# bio8: mean tempearature of the wettest quarter
# bio12: annual precipitation
# bio18: precipitation of warmest quarter. 
climate.vars <- bioclim.vars[, c('bio1', 'bio2', 'bio8', 'bio12', 'bio18')]

save(climate.vars, file = 'data/climate-data-full.rda')
