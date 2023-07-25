# summary-sim-spAbundance.R: this script summarizes the results from the 
#                            simulated example that is reflective of the 
#                            butterfly case study using the spAbundance
#                            R package.
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(coda)
library(spAbundance)
library(MCMCvis)
# For the breaks_pretty() function. 
library(scales)

# Read in results ---------------------------------------------------------
# Note the full object that comes from fitting a model in spAbundance is
# not provided on Github as it is too large. Instead, we 
# have separated out by different parameters and chains, which makes the 
# file sizes smaller and easy to store.
load("results/sim-beta-samples-1-chain.rda")
beta.samples.1 <- beta.samples
load("results/sim-beta-samples-2-chain.rda")
beta.samples.2 <- beta.samples
load("results/sim-beta-samples-3-chain.rda")
beta.samples.3 <- beta.samples
# Put into one large MCMC object (not an MCMC list)
beta.samples <- mcmc(rbind(beta.samples.1, beta.samples.2, beta.samples.3))
load("results/sim-beta-comm-samples-1-chain.rda")
beta.comm.samples.1 <- beta.comm.samples
load("results/sim-beta-comm-samples-2-chain.rda")
beta.comm.samples.2 <- beta.comm.samples
load("results/sim-beta-comm-samples-3-chain.rda")
beta.comm.samples.3 <- beta.comm.samples
beta.comm.samples <- mcmc(rbind(beta.comm.samples.1, beta.comm.samples.2, 
				beta.comm.samples.3))
load("results/sim-beta-star-samples-1-chain.rda")
beta.star.samples.1 <- beta.star.samples
load("results/sim-beta-star-samples-2-chain.rda")
beta.star.samples.2 <- beta.star.samples
load("results/sim-beta-star-samples-3-chain.rda")
beta.star.samples.3 <- beta.star.samples
beta.star.samples <- mcmc(rbind(beta.star.samples.1, beta.star.samples.2, 
				beta.star.samples.3))

# Load in spAbundance formatted data --------------------------------------
# Loads in an object called data.list and beta.true (the true beta values)
load("data/sim-spAbundance-data.rda")
N <- nrow(data.list$y)
sp.names <- paste0('sp', 1:N)

# Plot true beta values vs. estimated beta values -------------------------
beta.means <- apply(beta.samples, 2, mean)
plot(c(beta.true), beta.means, pch = 19)
abline(0, 1)

# Plot expected abundance over time ---------------------------------------
N <- length(sp.names)
# Species-specific rend values on the log scale.
trend.log.samples <- MCMCchains(beta.samples, params = 'year.ind', exact = FALSE)
# Community-level trend values on the log scale.
trend.comm.log.samples <- MCMCchains(beta.comm.samples, params = 'year.ind', exact = FALSE)
# Probably community level trend < 0
mean(trend.comm.log.samples < 0)

n.years <- n_distinct(data.list$covs$year.ind)
int.samples <- MCMCchains(beta.samples, params = 'Intercept', exact = FALSE)
int.comm.samples <- MCMCchains(beta.comm.samples, params = 'Intercept', exact = FALSE)
# Extract year random effects
trend.re.samples <- MCMCchains(beta.star.samples, params = 'year.ind', exact = FALSE)
# Order: sample, years, species.
trend.re.samples <- array(MCMCchains(trend.re.samples, params = 'Intercept', exact = FALSE),
			  dim = c(nrow(trend.re.samples), n.years, N))
# Average for the community-level values.
trend.re.comm.samples <- apply(trend.re.samples, c(1, 2), mean)
pred.df <- data.frame(data.ind.2 = 0, data.ind.3 = 0, data.ind.4 = 0,
		      data.ind.5 = 0, survey.cov = 0, year.cov = 1:n.years, week.cov = 0)
pred.df$year.s <- (1:n.years - mean(data.list$covs$year.ind, na.rm = TRUE)) / sd(data.list$covs$year.ind,
									      na.rm = TRUE)
n.samples <- nrow(trend.re.samples)
# Generate full posterior samples of abundance in each year, when using just
# the linear trend and when using all random effects.
# Including random effects
mu.samples <- array(NA, dim = c(n.samples, n.years, N))
mu.comm.samples <- array(NA, dim = c(n.samples, n.years))
# Just the linear trend portion
mu.trend.samples <- array(NA, dim = c(n.samples, n.years, N))
mu.comm.trend.samples <- array(NA, dim = c(n.samples, n.years))
for (j in 1:n.samples) {
  for (i in 1:N) {
    mu.samples[j, , i] <- exp(int.samples[j, i] +
			      trend.log.samples[j, i] * pred.df$year.s +
			      trend.re.samples[j, , i])
    mu.comm.samples[j, ] <- exp(int.comm.samples[j, ] +
				trend.comm.log.samples[j, ] * pred.df$year.s +
				trend.re.comm.samples[j, ])
    mu.trend.samples[j, , i] <- exp(int.samples[j, i] +
			      trend.log.samples[j, i] * pred.df$year.s)
    mu.comm.trend.samples[j, ] <- exp(int.comm.samples[j, ] +
				      trend.comm.log.samples[j, ] * pred.df$year.s)
  }
}
prob.neg <- apply(trend.log.samples, 2, function(a) mean(a < 0))
prob.neg <- c(prob.neg, mean(trend.comm.log.samples < 0))
names(prob.neg) <- c(sp.names, 'COMM')
col.ind <- ifelse(prob.neg > 0.6, 1,
		  ifelse(prob.neg < 0.4, 2, 0))
# Quantiles of species-specific yearly abundance values.
mu.quants <- apply(mu.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
mu.trend.quants <- apply(mu.trend.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
mu.comm.quants <- apply(mu.comm.samples, 2, quantile, c(0.025, 0.5, 0.975))
mu.comm.trend.quants <- apply(mu.comm.trend.samples, 2, quantile, c(0.025, 0.5, 0.975))
my.cols <- c('0' = 'gray', '1' = '#FFF4B3', '2' = '#2166AC')
plot.df <- data.frame(med = c(c(mu.quants[2, , ]), c(mu.comm.quants[2, ])),
		      low = c(c(mu.quants[1, , ]), c(mu.comm.quants[1, ])),
		      high = c(c(mu.quants[3, , ]), c(mu.comm.quants[3, ])),
		      med.trend = c(c(mu.trend.quants[2, , ]), c(mu.comm.trend.quants[2, ])),
		      low.trend = c(c(mu.trend.quants[1, , ]), c(mu.comm.trend.quants[1, ])),
		      high.trend = c(c(mu.trend.quants[3, , ]), c(mu.comm.trend.quants[3, ])),
		      sp = factor(rep(c(sp.names, 'COMM'), each = n.years),
				  levels = c(sort(sp.names), 'COMM'), ordered = TRUE),
		      col.ind = as.character(rep(col.ind, each = n.years)),
		      year = rep(1:n.years, times = N + 1))
# Make the plot.
ggplot(data = plot.df, aes(x = year, y = med)) +
  geom_ribbon(aes(ymin = low.trend, ymax = high.trend, fill = col.ind), alpha = 0.5) +
  # geom_line(aes(y = med.trend, col = col.ind), size = 1.2, lineend = 'butt') +
  geom_line(aes(y = med.trend), size = 1, lineend = 'butt') +
  geom_point(size = 3.5) +
  geom_segment(aes(y = low,
		   yend = high, x = year, xend = year),
		 lineend = 'butt', size = 1) +
  theme_bw(base_size = 22) +
  facet_wrap(vars(sp), scales = 'free_y', ncol = 4) +
  # scale_color_manual(values = my.cols) +
  scale_fill_manual(values = my.cols) +
  scale_y_continuous(breaks = breaks_pretty(n = 4)) +
  guides(fill = 'none', color = 'none') +
  theme(axis.text.x = element_text(size = 16),
	axis.text.y = element_text(size = 16)) +
  labs(x = 'Year', y = 'Relative Abundance Index')

