# spOccupancy-data-prep.R: this script preps the data for fitting an 
#                          integrated community occupancy model using the 
#                          spOccupancy R package. 
# Author: Jeffrey W. Doser 
rm(list = ls())
library(tidyverse)
library(coda)
library(sf)
library(devtools)

# Read in data ------------------------------------------------------------
load("data/ne-data-bundle.rda")
load('data/climate-data-full.rda')

# Data prep ---------------------------------------------------------------
# Total number of cells in the study area.
J <- nrow(grid.ne)
# Number of cells with eBird data
J.ebird <- n_distinct(ebird.df$cell)
# Number of cells with BBS data
J.bbs <- n_distinct(y.bbs$cell)
# Number of species
N <- n_distinct(ebird.df$sp)

# Filter the eBird data to use data from 3 weeks in June when closure is
# reasonable, and the time period matches up well with the BBS data. 
ebird.df <- ebird.df %>%
  filter(week %in% c(23, 24, 25))
occ.covs <- data.frame(occ.covs, climate.vars)
# Set missing values to their means. 
for (i in 1:ncol(occ.covs)) {
  occ.covs[which(is.na(occ.covs[, i])), i] <- mean(occ.covs[, i], na.rm = TRUE)
}

# Only extract the 5x5km cells with eBird or BBS data for model fitting ---
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

# Format data for spOccupancy ---------------------------------------------
# Occurrence covariates ---------------
occ.covs <- occ.covs.small
# Detection-nondetection data ---------
# BBS -------------
# Reorder by cell (spatial location) then species, then only grab the first
# observer for the couple of cells with multiple observers. 
bbs.df <- y.bbs %>%
  group_by(sp, cell) %>%
  summarize(julian = first(julian), 
	    obs = first(obs), 
	    Count10 = first(Count10),
	    Count20 = first(Count20), 
	    Count30 = first(Count30),
	    Count40 = first(Count40),
	    Count50 = first(Count50)) %>%
  arrange(cell, sp) %>%
  ungroup() %>%
  select(starts_with('Count'), sp, cell, julian, obs)
# 5 spatial replicates for each BBS data point. 
K.bbs <- 5
y.bbs <- array(NA, dim = c(N, J.bbs, K.bbs))
for (i in 1:K.bbs) {
  y.bbs[, , i] <- matrix(as.matrix(bbs.df[, i]), N, J.bbs)
}
sp.names <- unique(bbs.df$sp)
rownames(y.bbs) <- sp.names
# eBird -----------
# Reorder by week, then cell, then species within cell. 
ebird.df <- ebird.df %>%
  arrange(week, cell, sp)
K.ebird <- n_distinct(ebird.df$week)
sites.ebird <- unique(ebird.df$cell)
J.ebird <- length(sites.ebird)
y.ebird <- array(NA, dim = c(N, J.ebird, K.ebird))
week.ebird <- unique(ebird.df$week)
for (i in 1:nrow(ebird.df)) {
  if (i %% 10000 == 0) {
    print(i)
  }
  sp.indx <- which(sp.names == ebird.df$sp[i])
  week.indx <- which(week.ebird == ebird.df$week[i])
  site.indx <- which(sites.ebird == ebird.df$cell[i])
  y.ebird[sp.indx, site.indx, week.indx] <- ebird.df$y[i]
}
rownames(y.ebird) <- sp.names
y <- list(y.bbs, y.ebird)

# Detection covariates ----------------
obs.bbs <- bbs.df %>%
  filter(sp == sp.names[1]) %>%
  pull(obs) 
julian.bbs <- bbs.df %>%
  filter(sp == sp.names[1]) %>%
  pull(julian)
X.p.bbs <- list(observer = obs.bbs, 
		julian = julian.bbs)
single.sp.ebird <- ebird.df %>%
  filter(sp == sp.names[1])
day.ebird <- matrix(NA, J.ebird, K.ebird)
time.ebird <- matrix(NA, J.ebird, K.ebird)
length.ebird <- matrix(NA, J.ebird, K.ebird)
dist.ebird <- matrix(NA, J.ebird, K.ebird)
obsv.ebird <- matrix(NA, J.ebird, K.ebird)
for (i in 1:nrow(single.sp.ebird)) {
  site.indx <- which(sites.ebird == single.sp.ebird$cell[i])
  week.indx <- which(week.ebird == single.sp.ebird$week[i])
  day.ebird[site.indx, week.indx] <- single.sp.ebird$day[i]
  time.ebird[site.indx, week.indx] <- single.sp.ebird$time[i]
  length.ebird[site.indx, week.indx] <- single.sp.ebird$length[i]
  dist.ebird[site.indx, week.indx] <- single.sp.ebird$dist[i]
  obsv.ebird[site.indx, week.indx] <- single.sp.ebird$obsv[i]
}
X.p.ebird <- list(day = day.ebird, time = time.ebird, 
		  length = length.ebird, dist = dist.ebird, 
		  obsv = obsv.ebird)
det.covs <- list(X.p.bbs, X.p.ebird)

# Site indices ------------------------
sites.bbs <- bbs.df %>%
  filter(sp == sp.names[1]) %>%
  pull(cell)
sites <- list(sites.bbs, sites.ebird)
# Species index -----------------------
species <- list(sp.names, sp.names)

# Put data into a list for spOccupancy
data.list <- list(y = y, occ.covs = occ.covs, 
		  det.covs = det.covs, sites = sites, 
		  species = species)
# Save resulting data -----------------------------------------------------
save(data.list, file = 'data/spOccupancy-data.rda')

