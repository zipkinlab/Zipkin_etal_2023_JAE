# summary-nimble.R: this file summarizes results from the integrated 
#                   community occupancy model bird case study fit with NIMBLE. 
# Author: Jeffrey W. Doser	
rm(list = ls())
library(coda)
library(nimble)
library(tidyverse)
# For plotting
library(sf)
library(stars)
library(viridis)
library(RColorBrewer)
library(ggpubr)

# Functions ---------------------------------------------------------------
logit <- function(theta, a = 0, b = 1) {log((theta-a)/(b-theta))}
logit.inv <- function(z, a = 0, b = 1) {b-(b-a)/(1+exp(z))}

# Load data and results ---------------------------------------------------
# Note these files are large, and take up a substantial amount of RAM. 
load("data/ne-data-bundle.rda")
load('data/climate-data-full.rda')
load("results/ne-icom-1-chain-2023-01-12.R")
samples.1 <- samples
load("results/ne-icom-2-chain-2023-01-12.R")
samples.2 <- samples
load("results/ne-icom-3-chain-2023-01-12.R")
samples.3 <- samples

samples.list <- mcmc.list(samples.1, samples.2, samples.3)

# Estimate species richness post-hoc --------------------------------------
ebird.df <- ebird.df %>%
  filter(week %in% c(23, 24, 25))
# Add climate variables to other occurrence variables
occ.covs <- data.frame(occ.covs, climate.vars)
# Number of species
N <- n_distinct(bbs.sf$sp)
param.names <- unlist(dimnames(samples)[2])
# Covariate data used to fit the model
unique.cells <- sort(unique(c(unique(y.bbs$cell), unique(ebird.df$cell))))
J <- length(unique.cells)
occ.covs.small <- occ.covs[unique.cells, ]
# Extract species-specific occurrence parameters
curr.indx <- which(substr(param.names, 1, 5) == 'beta[')
beta.samples <- samples.list[, curr.indx]
# Intercept
beta.1.samples <- do.call('rbind', beta.samples[, 1:N])
# Linear elevation
beta.2.samples <- do.call('rbind', beta.samples[, (N + 1):(N * 2)])
# Quadratic elevation
beta.3.samples <- do.call('rbind', beta.samples[, (N * 2 + 1):(N * 3)])
# Forest cover
beta.4.samples <- do.call('rbind', beta.samples[, (N * 3 + 1):(N * 4)])
# Bioclim 1
beta.5.samples <- do.call('rbind', beta.samples[, (N * 4 + 1):(N * 5)])
# Bioclim 2
beta.6.samples <- do.call('rbind', beta.samples[, (N * 5 + 1):(N * 6)])
# Bioclim 8
beta.7.samples <- do.call('rbind', beta.samples[, (N * 6 + 1):(N * 7)])
# Bioclim 12
beta.8.samples <- do.call('rbind', beta.samples[, (N * 7 + 1):(N * 8)])
# Bioclim 18
beta.9.samples <- do.call('rbind', beta.samples[, (N * 8 + 1):(N * 9)])
# Predict across all locations (including fitted)
occ.covs.pred <- occ.covs
# Get mean and sd for covariates used to fit the model
mean.elev.fit <- mean(occ.covs.small$elev, na.rm = TRUE)
sd.elev.fit <- sd(occ.covs.small$elev, na.rm = TRUE)
mean.pf.fit <- mean(occ.covs.small$pf, na.rm = TRUE)
sd.pf.fit <- sd(occ.covs.small$pf, na.rm = TRUE)
mean.bio1.fit <- mean(occ.covs.small$bio1, na.rm = TRUE)
sd.bio1.fit <- sd(occ.covs.small$bio1, na.rm = TRUE)
mean.bio2.fit <- mean(occ.covs.small$bio2, na.rm = TRUE)
sd.bio2.fit <- sd(occ.covs.small$bio2, na.rm = TRUE)
mean.bio8.fit <- mean(occ.covs.small$bio8, na.rm = TRUE)
sd.bio8.fit <- sd(occ.covs.small$bio8, na.rm = TRUE)
mean.bio12.fit <- mean(occ.covs.small$bio12, na.rm = TRUE)
sd.bio12.fit <- sd(occ.covs.small$bio12, na.rm = TRUE)
mean.bio18.fit <- mean(occ.covs.small$bio18, na.rm = TRUE)
sd.bio18.fit <- sd(occ.covs.small$bio18, na.rm = TRUE)
# Design matrix for model predictions
X.0 <- cbind(1, (occ.covs.pred$elev - mean.elev.fit) / sd.elev.fit, 
	     ((occ.covs.pred$elev - mean.elev.fit) / sd.elev.fit)^2, 
	     (occ.covs.pred$pf - mean.pf.fit) / sd.pf.fit, 
             (occ.covs.pred$bio1 - mean.bio1.fit) / sd.bio1.fit,
             (occ.covs.pred$bio2 - mean.bio2.fit) / sd.bio2.fit,
             (occ.covs.pred$bio8 - mean.bio8.fit) / sd.bio8.fit,
             (occ.covs.pred$bio12 - mean.bio12.fit) / sd.bio12.fit,
             (occ.covs.pred$bio18 - mean.bio18.fit) / sd.bio18.fit)
# Set missing values to 0
X.0[which(is.na(X.0))] <- 0
colnames(X.0) <- c('intercept', 'elev', 'elev.2', 'forest', 'bio1', 
                   'bio2', 'bio8', 'bio12', 'bio18')
# Number of prediction locations
J.0 <- nrow(X.0)
# Number of posterior samples
n.post <- nrow(beta.1.samples)
# Temporary occurrence probability for each species (not saved)
psi.tmp <- matrix(NA, n.post, J.0)
z.tmp <- matrix(NA, N, J.0)
psi.means <- array(NA, dim = c(N, J.0))
# Species richness samples
rich.samples <- matrix(NA, n.post, J.0)
# This takes a few minutes. 
# Get means of species-specific occurrence
for (i in 1:N) {
print(paste("Currently on species ", i, " out of ", N, sep = ''))
  for (a in 1:n.post) {
    psi.tmp[a, ] <- logit.inv(beta.1.samples[a, i] +
                              beta.2.samples[a, i] * X.0[, 2] +
                              beta.3.samples[a, i] * X.0[, 3] +
                              beta.4.samples[a, i] * X.0[, 4] + 
                              beta.5.samples[a, i] * X.0[, 5] + 
                              beta.6.samples[a, i] * X.0[, 6] + 
                              beta.7.samples[a, i] * X.0[, 7] + 
                              beta.8.samples[a, i] * X.0[, 8] + 
                              beta.9.samples[a, i] * X.0[, 9])
  } # a (samples)
  psi.means[i, ] <- apply(psi.tmp, 2, mean)
} # i (species)
# Predict Richness --------------------
for (a in 1:n.post) {
  print(paste("Currently on iteration ", a, " out of ", n.post, sep = ''))
  for (i in 1:N) {
    z.tmp[i,] <- rbinom(J.0, 1, logit.inv(beta.1.samples[a, i] + 
			                  beta.2.samples[a, i] * X.0[, 2] + 
			                  beta.3.samples[a, i] * X.0[, 3] +
			                  beta.4.samples[a, i] * X.0[, 4] + 
                                          beta.5.samples[a, i] * X.0[, 5] +
			                  beta.6.samples[a, i] * X.0[, 6] + 
			                  beta.7.samples[a, i] * X.0[, 7] + 
			                  beta.8.samples[a, i] * X.0[, 8] + 
			                  beta.9.samples[a, i] * X.0[, 9]))
  } # i (species)
  rich.samples[a, ] <- apply(z.tmp, 2, sum)
} # a (iteration)

# Extract species richness estimates --------------------------------------
rich.mean <- apply(rich.samples, 2, mean, na.rm = TRUE)
rich.sd <- apply(rich.samples, 2, sd, na.rm = TRUE)
rich.low <- apply(rich.samples, 2, quantile, 0.025, na.rm = TRUE)
rich.high <- apply(rich.samples, 2, quantile, 0.975, na.rm = TRUE)

# Put in richness values directly into grid
grid.ne$rich.mean <- rich.mean
grid.ne$rich.sd <- rich.sd
grid.ne$ci.width <- rich.high - rich.low

# Filter the eBird data to those used to fit the model.  
ebird.df <- ebird.df %>%
  filter(week %in% c(23, 24, 25))

# Calculate the number of eBird checklists in each cell. 
grid.ne$eBirdCheck <- 0
# Number of checklists in a cell can be 0, 1, 2, or 3
tmp.eb <- ebird.df %>%
  group_by(cell) %>%
  summarize(n.check = n_distinct(listID))
for (i in 1:nrow(grid.ne)) {
  if (i %in% tmp.eb$cell) {
    grid.ne$eBirdCheck[i] <- tmp.eb$n.check[which(tmp.eb$cell == i)]
  }
}
grid.ne$eBirdCheck <- factor(grid.ne$eBirdCheck, levels = 0:3)

# Get map of US -----------------------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
ne.states <- usa %>% 
  filter(ID %in% c('connecticut', 'delaware', 'maine', 'maryland', 
		     'massachusetts', 'new hampshire', 'new jersey', 
		     'new york', 'pennsylvania', 'rhode island', 
		     'vermont'))

# Map of the data used to fit the model -----------------------------------
# Load BBS locations
bbs.locs <- bbs.sf %>%
  filter(sp == 'HAWO')
grid.plot <- st_intersection(grid.ne, st_make_valid(ne.states))
grid.plot$eBirdCheck <- factor(grid.plot$eBirdCheck, levels = 0:3)
counts.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = eBirdCheck, col = eBirdCheck)) + 
  geom_sf(data = ne.states, col = 'gray', fill = NA) + 
  scale_fill_viridis_d() +
  scale_color_viridis_d() + 
  guides(col = "none") + 
  geom_sf(data = bbs.locs, pch = 21, fill = 'white', size = 2, col = 'black') +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "# eBird Lists") +
  theme(legend.position = c(0.83, 0.19), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.title = element_text(size = 14), 
        legend.text = element_text(size = 12))
# counts.plot

# Species Richness Map ----------------------------------------------------
rich.mean.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = rich.mean, col = rich.mean)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  # scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  # scale_color_gradientn(colors = viridis(10), na.value = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Richness\nMean") + 
  guides(col = 'none') +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
rich.mean.plot
# Map of CI width
rich.ci.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = ci.width, col = ci.width)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  # scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  # scale_color_gradientn(colors = viridis(10), na.value = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  guides(col = 'none') + 
  labs(x = "Longitude", y = "Latitude", fill = "CI Width") +
  theme(legend.position = c(0.84, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
rich.ci.plot

# Two example species -----------------------------------------------------
sp.codes <- unique(y.bbs$sp)
grid.ne$BTNW.occ <- psi.means[which(sp.codes == 'BTNW'), ]
grid.ne$CERW.occ <- psi.means[which(sp.codes == 'CERW'), ]
grid.plot <- st_intersection(grid.ne, st_make_valid(ne.states))
BTNW.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = BTNW.occ, col = BTNW.occ)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "BTNW\nOccurrence") + 
  guides(col = 'none') +
  theme(legend.position = c(0.82, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
BTNW.plot
CERW.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = CERW.occ, col = CERW.occ)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "CERW\nOccurrence") + 
  guides(col = 'none') +
  theme(legend.position = c(0.82, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
CERW.plot

# Plot of covariate effects -----------------------------------------------
param.names <- attr(samples, 'dimnames')[[2]]
curr.indx <- which(substr(param.names, 1, 4) == 'beta')
# Samples of species and community-level regression coefficients
beta.samples <- do.call("rbind", samples.list[, curr.indx])
# Order: elevation, elevation^2, forest cover, then all the bioclim variables.  
# Hardcoding because I'm lazy
# Includes species-specific effects of forest cover as well as community-level. 
beta.for.samples <- beta.samples[, (N * 3 + 1):(N * 4)]
beta.comm.for.samples <- beta.samples[, ncol(beta.samples)]
# Species codes
sp.code.factor <- factor(unique(y.bbs$sp), 
		         levels = unique(y.bbs$sp))
sp.codes <- as.numeric(sp.code.factor)
N <- length(sp.codes)
cov.plot.df <- data.frame(for.mean = apply(beta.for.samples, 2, mean), 
			  for.low = apply(beta.for.samples, 2, quantile, 0.25), 
			  for.lowest = apply(beta.for.samples, 2, quantile, 0.025), 
			  for.high = apply(beta.for.samples, 2, quantile, 0.75), 
			  for.highest = apply(beta.for.samples, 2, quantile, 0.975), 
			  sp = sp.codes)
# Rearrange and add things to get the plot to display effects in increasing
# order. 
cov.plot.df <- cov.plot.df %>%
  arrange(for.mean)
cov.plot.df$sp.factor <- as.character(sp.code.factor[cov.plot.df$sp])
cov.plot.df$sort.sp <- 1:N

# Add in the community level covariate
comm.plot.df <- data.frame(for.mean = mean(beta.comm.for.samples),
			   for.low = quantile(beta.comm.for.samples, 0.25), 
			   for.lowest = quantile(beta.comm.for.samples, 0.025), 
			   for.high = quantile(beta.comm.for.samples, 0.75), 
			   for.highest = quantile(beta.comm.for.samples, 0.975), 
			   sp = 28, 
			   sp.factor = 'COMM',
                           sort.sp = N + 1)
cov.plot.df <- rbind(cov.plot.df, comm.plot.df)
cov.plot.df$sp.factor <- factor(cov.plot.df$sp.factor, levels = c(unique(cov.plot.df$sp.factor)))

for.cov.plot <- ggplot(data = cov.plot.df, aes(x = sort.sp, fill = for.mean, group = sp.factor)) + 
  geom_hline(yintercept = 0, col = 'black', size = 0.75, lty = 2) + 
  geom_vline(xintercept = 27.5, col = 'black', size = 0.5, lty = 1) + 
  geom_boxplot(aes(ymin = for.lowest, lower = for.low, middle = for.mean, 
		   upper = for.high, ymax = for.highest), stat = 'identity', col = 'black') + 
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC', 
  	               na.value = NA) + 
  theme_bw(base_size = 18) + 
  guides(fill = "none") + 
  labs(x = "Species", y = "Effect of Forest Cover") + 
  scale_x_continuous(breaks = 1:(N+1), labels = cov.plot.df$sp.factor) +  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  coord_flip()
for.cov.plot

# Put everything together to produce Figure 2 in the main text (note Fig 2 in the 
# main text was created using the spOccupancy model, so there may be very slight
# differences).
ggarrange(rich.mean.plot, for.cov.plot, BTNW.plot, CERW.plot, nrow = 2, ncol = 2, 
          labels = c('(a)', '(b)', '(c)', '(d)'), 
          font.label = list(size = 18))
