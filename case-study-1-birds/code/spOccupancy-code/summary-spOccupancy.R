# summary-spOccupancy.R: this file summarizes results from the integrated
#                        community occupancy model bird case study fit with 
#                        spOccupancy. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(spOccupancy)
library(coda)
library(tidyverse)
library(MCMCvis)
# For plotting
library(sf)
library(stars)
library(viridis)
library(RColorBrewer)
library(ggpubr)

# Note that this script requires a substantial amount of RAM to generate
# some of the prediction objects (i.e., the full posterior predictions
# for each species across the 22,513 grid cells). 

# Load data and results ---------------------------------------------------
load("data/ne-data-bundle.rda")
load('data/climate-data-full.rda')
# Note these files are too large to store on GitHub (a little over 1GB each). 
# To create these files, run the "main-spOccupancy.R" script. 
load('results/spOccupancy-icom-1-chain-2023-01-11.R')
out.1 <- out
load('results/spOccupancy-icom-1-chain-2023-01-11.R')
out.2 <- out
load('results/spOccupancy-icom-1-chain-2023-01-11.R')
out.3 <- out

# Generate predictions ----------------------------------------------------
occ.covs <- data.frame(occ.covs, climate.vars)
# Set missing values to their means.
for (i in 1:ncol(occ.covs)) {
  occ.covs[which(is.na(occ.covs[, i])), i] <- mean(occ.covs[, i], na.rm = TRUE)
}
# Load data list used to fit the model
load("data/spOccupancy-data.rda")
# Standardize prediction covariates by the values used to fit the model. 
elev.0 <- (occ.covs$elev - mean(data.list$occ.covs$elev)) / sd(data.list$occ.covs$elev)
elev.0.sq <- elev.0^2
pf.0 <- (occ.covs$pf - mean(data.list$occ.covs$pf)) / sd(data.list$occ.covs$pf)
bio1.0 <- (occ.covs$bio1 - mean(data.list$occ.covs$bio1)) / sd(data.list$occ.covs$bio1)
bio2.0 <- (occ.covs$bio2 - mean(data.list$occ.covs$bio2)) / sd(data.list$occ.covs$bio2)
bio8.0 <- (occ.covs$bio8 - mean(data.list$occ.covs$bio8)) / sd(data.list$occ.covs$bio8)
bio12.0 <- (occ.covs$bio12 - mean(data.list$occ.covs$bio12)) / sd(data.list$occ.covs$bio12)
bio18.0 <- (occ.covs$bio18 - mean(data.list$occ.covs$bio18)) / sd(data.list$occ.covs$bio18)
# Create the design matrix. 
X.0 <- as.matrix(data.frame(intercept = 1, 
			    elev = elev.0,
			    elev.sq = elev.0.sq, 
			    pf = pf.0, 
			    bio1 = bio1.0, 
			    bio2 = bio2.0,
			    bio8 = bio8.0,
			    bio12 = bio12.0,
			    bio18 = bio18.0))
# Generate predictions using spOccupancy for each of the three chains
out.pred.1 <- predict(out.1, X.0)
out.pred.2 <- predict(out.2, X.0)
out.pred.3 <- predict(out.3, X.0)

# Calculate richness
rich.samples.1 <- apply(out.pred.1$z.0.samples, c(1, 3), sum)
rich.samples.2 <- apply(out.pred.2$z.0.samples, c(1, 3), sum)
rich.samples.3 <- apply(out.pred.3$z.0.samples, c(1, 3), sum)
rich.samples <- rbind(rich.samples.1, rich.samples.2, rich.samples.3)
rich.mean <- apply(rich.samples, 2, mean, na.rm = TRUE)
rich.sd <- apply(rich.samples, 2, sd, na.rm = TRUE)
rich.low <- apply(rich.samples, 2, quantile, 0.025, na.rm = TRUE)
rich.high <- apply(rich.samples, 2, quantile, 0.975, na.rm = TRUE)
# Put in richness values directly into grid
grid.ne$rich.mean <- rich.mean
grid.ne$rich.sd <- rich.sd
grid.ne$ci.width <- rich.high - rich.low

# Generate map of study area and data -------------------------------------
# Get sf object of states in this region. 
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
ne.states <- usa %>%
  filter(ID %in% c('connecticut', 'delaware', 'maine', 'maryland',
		     'massachusetts', 'new hampshire', 'new jersey',
		     'new york', 'pennsylvania', 'rhode island',
		     'vermont'))
# Filter the eBird data to the data we used to fit the model. 
ebird.df <- ebird.df %>%
  dplyr::filter(week %in% c(23, 24, 25))

# Number of eBird checklists in each cell
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
# Get BBS locations
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
# ggsave(device = 'png', filename = 'figures/ne-data.png', height = 6, width = 6,
#        units = 'in')
# ggsave(device = 'pdf', filename = 'figures/ne-data.pdf', height = 6, width = 6,
#        units = 'in')

# Species Richness Map ----------------------------------------------------
rich.mean.plot <- ggplot() +
  geom_sf(data = grid.plot, aes(fill = rich.mean, col = rich.mean)) +
  geom_sf(data = ne.states, col = 'black', fill = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "Richness\nMean") +
  guides(col = 'none') +
  theme(legend.position = c(0.84, 0.17),
        legend.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))
# rich.mean.plot
# ggsave(device = 'png', filename = 'figures/ne-richness.png', height = 7, width = 10,
#        units = 'in')
# ggsave(device = 'pdf', filename = 'figures/ne-richness.pdf', height = 7, width = 10,
#        units = 'in')
# Map of CI width
rich.ci.plot <- ggplot() +
  geom_sf(data = grid.plot, aes(fill = ci.width, col = ci.width)) +
  geom_sf(data = ne.states, col = 'black', fill = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  guides(col = 'none') +
  labs(x = "Longitude", y = "Latitude", fill = "CI Width") +
  theme(legend.position = c(0.84, 0.17),
        legend.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))
# rich.ci.plot
# ggsave(device = 'png', filename = 'figures/ne-richness-ci.png', height = 7, width = 10,
#        units = 'in')
# ggsave(device = 'pdf', filename = 'figures/ne-richness-ci.pdf', height = 7, width = 10,
#        units = 'in')

# Species-specific maps of occurrence probability -------------------------
sp.codes <- unique(y.bbs$sp)
curr.sp <- which(sp.codes == 'BTNW')
BTNW.psi.samples <- rbind(out.pred.1$psi.0.samples[, curr.sp, ], 
			  out.pred.2$psi.0.samples[, curr.sp, ],
			  out.pred.3$psi.0.samples[, curr.sp, ])
grid.ne$BTNW.occ <- apply(BTNW.psi.samples, 2, mean)
grid.ne$BTNW.ci.width <- apply(BTNW.psi.samples, 2, quantile, 0.975) - 
                         apply(BTNW.psi.samples, 2, quantile, 0.025)
curr.sp <- which(sp.codes == 'CERW')
CERW.psi.samples <- rbind(out.pred.1$psi.0.samples[, curr.sp, ], 
			  out.pred.2$psi.0.samples[, curr.sp, ],
			  out.pred.3$psi.0.samples[, curr.sp, ])
grid.ne$CERW.occ <- apply(CERW.psi.samples, 2, mean)
grid.ne$CERW.ci.width <- apply(CERW.psi.samples, 2, quantile, 0.975) - 
                         apply(CERW.psi.samples, 2, quantile, 0.025)
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
BTNW.ci.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = BTNW.ci.width, col = BTNW.ci.width)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "BTNW\nCI Width") + 
  guides(col = 'none') +
  theme(legend.position = c(0.82, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
# BTNW.plot
BTNW.ci.plot
ggsave(device = 'png', filename = 'figures/BTNW-psi-ci.png', height = 7, width = 10,
       units = 'in')
CERW.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = CERW.occ, col = CERW.occ)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  # scale_fill_gradientn(colors = viridis(10), na.value = NA) +
  # scale_color_gradientn(colors = viridis(10), na.value = NA) +
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "CERW\nOccurrence") + 
  guides(col = 'none') +
  theme(legend.position = c(0.82, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
CERW.ci.plot <- ggplot() + 
  geom_sf(data = grid.plot, aes(fill = CERW.ci.width, col = CERW.ci.width)) + 
  geom_sf(data = ne.states, col = 'black', fill = NA) + 
  scale_fill_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  scale_color_gradientn(colors = brewer.pal(9, 'YlOrRd'), na.value = NA) +
  theme_bw(base_size = 18) +
  labs(x = "Longitude", y = "Latitude", fill = "CERW\nCI Width") + 
  guides(col = 'none') +
  theme(legend.position = c(0.82, 0.17), 
        legend.background = element_rect(fill = NA), 
        axis.text.x = element_text(angle = 45, hjust = 1))
# CERW.plot
CERW.ci.plot
ggsave(device = 'png', filename = 'figures/CERW-psi-ci.png', height = 7, width = 10,
       units = 'in')

# Plot effect of forest cover ---------------------------------------------
beta.samples <- mcmc.list(list(out.1$beta.samples, out.2$beta.samples, out.3$beta.samples))
beta.names <- colnames(out.1$beta.samples)
beta.for.samples <- MCMCchains(beta.samples, param = 'pf', exact = FALSE)
beta.comm.samples <- mcmc.list(out.1$beta.comm.samples, 
			       out.2$beta.comm.samples, 
			       out.3$beta.comm.samples)
beta.for.comm.samples <- MCMCchains(beta.comm.samples, param = 'pf', exact = FALSE)
# Species codes, plus the community code (COMM)
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
comm.plot.df <- data.frame(for.mean = mean(beta.for.comm.samples),
			   for.low = quantile(beta.for.comm.samples, 0.25), 
			   for.lowest = quantile(beta.for.comm.samples, 0.025), 
			   for.high = quantile(beta.for.comm.samples, 0.75), 
			   for.highest = quantile(beta.for.comm.samples, 0.975), 
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
ggsave(device = 'png', filename = 'figures/ne-forest-cover.png', height = 6, width = 6, 
       units = 'in')
ggsave(device = 'pdf', filename = 'figures/ne-forest-cover.pdf', height = 6, width = 6, 
       units = 'in')

ggarrange(rich.mean.plot, for.cov.plot, BTNW.plot, CERW.plot, nrow = 2, ncol = 2, 
          labels = c('(a)', '(b)', '(c)', '(d)'), 
          font.label = list(size = 18))
ggsave(device = 'pdf', filename = 'figures/icom-main-figure.pdf', height = 14, width = 12, 
       units = 'in')
ggsave(device = 'png', filename = 'figures/icom-main-figure.png', height = 14, width = 12, 
       units = 'in')

# Calculate raw data summaries --------------------------------------------
tmp <- apply(out$y[[1]], c(1, 2), function(a) sum(a) > 0)
bbs.raw.occ <- apply(tmp, 1, mean)
mean(bbs.raw.occ < .3)

tmp <- apply(out$y[[2]], c(1, 2), function(a) sum(a, na.rm = TRUE) > 0)
eb.raw.occ <- apply(tmp, 1, mean)
mean(eb.raw.occ < .3)

# Total number of detections for each species in each data set. 
apply(out$y[[1]], 1, sum, na.rm = TRUE)
apply(out$y[[2]], 1, sum, na.rm = TRUE)
