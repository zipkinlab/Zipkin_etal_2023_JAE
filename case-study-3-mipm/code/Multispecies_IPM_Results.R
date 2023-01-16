#####################################################
## Integrated Community Models - Case Study Example #
## Multi-species Integrated Population Models       #
## Plotting results for Figure 4B-D                 #
## Created by: Courtney L. Davis                    #
## Last modified: January 16, 2023                  #
#####################################################

#####################################################
# Load packages
#####################################################
library(here)
library(IPMbook)
library(jagsUI)
library(ggplot2)
library(viridis)
library(patchwork)
library(reshape2)
library(dplyr)


#####################################################
# Load simulation results from Rdata file
#####################################################
load(here("results","Simulation_Results.Rdata"))

#nsim = 100 
#n.species = 10
#n.occasions = 10

#####################################################
# Compare simulation results to truth under which the data were simulated
#####################################################

### Calculating bias for community-level mean of juvenile survival
median_bias_mu.sj <- c(); low95_bias_mu.sj <- c(); high95_bias_mu.sj <- c(); low50_bias_mu.sj <- c(); high50_bias_mu.sj <- c()
for(s in 1:nsim){
  median_bias_mu.sj_rep <- median(((results_sims[[s]]$mu.sj - mu.sj)/mu.sj))
  median_bias_mu.sj <- c(median_bias_mu.sj, median_bias_mu.sj_rep)
  
  low95_bias_mu.sj_rep <- quantile((results_sims[[s]]$mu.sj - mu.sj)/mu.sj, probs = 0.025)
  low95_bias_mu.sj <- c(low95_bias_mu.sj, low95_bias_mu.sj_rep)
  
  high95_bias_mu.sj_rep <- quantile((results_sims[[s]]$mu.sj - mu.sj)/mu.sj, probs = 0.975)
  high95_bias_mu.sj <- c(high95_bias_mu.sj, high95_bias_mu.sj_rep)
  
  low50_bias_mu.sj_rep <- quantile((results_sims[[s]]$mu.sj - mu.sj)/mu.sj, probs = 0.25)
  low50_bias_mu.sj <- c(low50_bias_mu.sj, low50_bias_mu.sj_rep)
  
  high50_bias_mu.sj_rep <- quantile((results_sims[[s]]$mu.sj - mu.sj)/mu.sj, probs = 0.75)
  high50_bias_mu.sj <- c(high50_bias_mu.sj, high50_bias_mu.sj_rep)
}

# Summarize mean and quantiles
bias_mu.sj <- data.frame("param" = "Juvenile Survival", "median" = mean(median_bias_mu.sj), 
                         "lower95" = mean(low95_bias_mu.sj), 
                         "upper95" = mean(high95_bias_mu.sj), 
                         "lower50" = mean(low50_bias_mu.sj), 
                         "upper50" = mean(high50_bias_mu.sj), row.names = "")



### Calculating bias for community-level mean of adult survival
median_bias_mu.sa <- c(); low95_bias_mu.sa <- c(); high95_bias_mu.sa <- c(); low50_bias_mu.sa <- c(); high50_bias_mu.sa <- c()
for(s in 1:nsim){
  median_bias_mu.sa_rep <- median((results_sims[[s]]$mu.sa - mu.sa)/mu.sa)
  median_bias_mu.sa <- c(median_bias_mu.sa, median_bias_mu.sa_rep)
  
  low95_bias_mu.sa_rep <- quantile((results_sims[[s]]$mu.sa - mu.sa)/mu.sa, probs = 0.025)
  low95_bias_mu.sa <- c(low95_bias_mu.sa, low95_bias_mu.sa_rep)
  
  high95_bias_mu.sa_rep <- quantile((results_sims[[s]]$mu.sa - mu.sa)/mu.sa, probs = 0.975)
  high95_bias_mu.sa <- c(high95_bias_mu.sa, high95_bias_mu.sa_rep)
  
  low50_bias_mu.sa_rep <- quantile((results_sims[[s]]$mu.sa - mu.sa)/mu.sa, probs = 0.25)
  low50_bias_mu.sa <- c(low50_bias_mu.sa, low50_bias_mu.sa_rep)
  
  high50_bias_mu.sa_rep <- quantile((results_sims[[s]]$mu.sa - mu.sa)/mu.sa, probs = 0.75)
  high50_bias_mu.sa <- c(high50_bias_mu.sa, high50_bias_mu.sa_rep)
}

# Summarize mean and quantiles
bias_mu.sa <- data.frame("param" = "Adult Survival", "median" = mean(median_bias_mu.sa), 
                         "lower95" = mean(low95_bias_mu.sa), 
                         "upper95" = mean(high95_bias_mu.sa),
                         "lower50" = mean(low50_bias_mu.sa), 
                         "upper50" = mean(high50_bias_mu.sa), row.names = "")


### Calculating bias for community-level mean of adult fecundity
median_bias_mu.f <- c(); low95_bias_mu.f <- c(); high95_bias_mu.f <- c(); low50_bias_mu.f <- c(); high50_bias_mu.f <- c()
for(s in 1:nsim){
  median_bias_mu.f_rep <- median((results_sims[[s]]$mu.f - mu.f)/mu.f)
  median_bias_mu.f <- c(median_bias_mu.f, median_bias_mu.f_rep)
  
  low95_bias_mu.f_rep <- quantile((results_sims[[s]]$mu.f - mu.f)/mu.f, probs = 0.025)
  low95_bias_mu.f <- c(low95_bias_mu.f, low95_bias_mu.f_rep)
  
  high95_bias_mu.f_rep <- quantile((results_sims[[s]]$mu.f - mu.f)/mu.f, probs = 0.975)
  high95_bias_mu.f <- c(high95_bias_mu.f, high95_bias_mu.f_rep)
  
  low50_bias_mu.f_rep <- quantile((results_sims[[s]]$mu.f - mu.f)/mu.f, probs = 0.25)
  low50_bias_mu.f <- c(low50_bias_mu.f, low50_bias_mu.f_rep)
  
  high50_bias_mu.f_rep <- quantile((results_sims[[s]]$mu.f - mu.f)/mu.f, probs = 0.75)
  high50_bias_mu.f <- c(high50_bias_mu.f, high50_bias_mu.f_rep)
}

# Summarize mean and quantiles
bias_mu.f <- data.frame("param" = "Fecundity", "median" = mean(median_bias_mu.f), 
                        "lower95" = mean(low95_bias_mu.f), 
                        "upper95" = mean(high95_bias_mu.f),
                        "lower50" = mean(low50_bias_mu.f), 
                        "upper50" = mean(high50_bias_mu.f), row.names = "")



### Calculating bias for species-specific juvenile survival
median_bias_mean.sj <- c(); low95_bias_mean.sj <- c(); high95_bias_mean.sj <- c(); low50_bias_mean.sj <- c(); high50_bias_mean.sj <- c()
for(s in 1:nsim){
  median_bias_mean.sj_rep <- mean(apply(((results_sims[[s]]$mean.sj - mean.sj)/mean.sj),2,median))
  median_bias_mean.sj <- c(median_bias_mean.sj, median_bias_mean.sj_rep)
  
  low95_bias_mean.sj_rep <- mean(apply(((results_sims[[s]]$mean.sj - mean.sj)/mean.sj),2,quantile, probs = 0.025))
  low95_bias_mean.sj <- c(low95_bias_mean.sj, low95_bias_mean.sj_rep)
  
  high95_bias_mean.sj_rep <- mean(apply(((results_sims[[s]]$mean.sj - mean.sj)/mean.sj),2,quantile, probs = 0.975))
  high95_bias_mean.sj <- c(high95_bias_mean.sj, high95_bias_mean.sj_rep)
  
  low50_bias_mean.sj_rep <- mean(apply(((results_sims[[s]]$mean.sj - mean.sj)/mean.sj),2,quantile, probs = 0.25))
  low50_bias_mean.sj <- c(low50_bias_mean.sj, low50_bias_mean.sj_rep)
  
  high50_bias_mean.sj_rep <- mean(apply(((results_sims[[s]]$mean.sj - mean.sj)/mean.sj),2,quantile, probs = 0.75))
  high50_bias_mean.sj <- c(high50_bias_mean.sj, high50_bias_mean.sj_rep)
}

# Summarize mean and quantiles
bias_mean.sj <- data.frame("param" = "Juvenile Survival", "median" = mean(median_bias_mean.sj), 
                           "lower95" = mean(low95_bias_mean.sj), 
                           "upper95" = mean(high95_bias_mean.sj), 
                           "lower50" = mean(low50_bias_mean.sj), 
                           "upper50" = mean(high50_bias_mean.sj), row.names = "")


### Calculating bias for species-specific adult survival
median_bias_mean.sa <- c(); low95_bias_mean.sa <- c(); high95_bias_mean.sa <- c(); low50_bias_mean.sa <- c(); high50_bias_mean.sa <- c()
for(s in 1:nsim){
  median_bias_mean.sa_rep <- mean(apply(((results_sims[[s]]$mean.sa - mean.sa)/mean.sa),2,median))
  median_bias_mean.sa <- c(median_bias_mean.sa, median_bias_mean.sa_rep)
  
  low95_bias_mean.sa_rep <- mean(apply(((results_sims[[s]]$mean.sa - mean.sa)/mean.sa),2,quantile, probs = 0.025))
  low95_bias_mean.sa <- c(low95_bias_mean.sa, low95_bias_mean.sa_rep)
  
  high95_bias_mean.sa_rep <- mean(apply(((results_sims[[s]]$mean.sa - mean.sa)/mean.sa),2,quantile, probs = 0.975))
  high95_bias_mean.sa <- c(high95_bias_mean.sa, high95_bias_mean.sa_rep)
  
  low50_bias_mean.sa_rep <- mean(apply(((results_sims[[s]]$mean.sa - mean.sa)/mean.sa),2,quantile, probs = 0.25))
  low50_bias_mean.sa <- c(low50_bias_mean.sa, low50_bias_mean.sa_rep)
  
  high50_bias_mean.sa_rep <- mean(apply(((results_sims[[s]]$mean.sa - mean.sa)/mean.sa),2,quantile, probs = 0.75))
  high50_bias_mean.sa <- c(high50_bias_mean.sa, high50_bias_mean.sa_rep)
}

# Summarize mean and quantiles
bias_mean.sa <- data.frame("param" = "Adult Survival", "median" = mean(median_bias_mean.sa), 
                           "lower95" = mean(low95_bias_mean.sa), 
                           "upper95" = mean(high95_bias_mean.sa),
                           "lower50" = mean(low50_bias_mean.sa), 
                           "upper50" = mean(high50_bias_mean.sa), row.names = "")


### Calculating bias for species-specific fecundity
median_bias_mean.f <- c(); low95_bias_mean.f <- c(); high95_bias_mean.f <- c(); low50_bias_mean.f <- c(); high50_bias_mean.f <- c()
for(s in 1:nsim){
  median_bias_mean.f_rep <- mean(apply(((results_sims[[s]]$mean.f - mean.f)/mean.f),2,median))
  median_bias_mean.f <- c(median_bias_mean.f, median_bias_mean.f_rep)
  
  low95_bias_mean.f_rep <- mean(apply(((results_sims[[s]]$mean.f - mean.f)/mean.f),2,quantile, probs = 0.025))
  low95_bias_mean.f <- c(low95_bias_mean.f, low95_bias_mean.f_rep)
  
  high95_bias_mean.f_rep <- mean(apply(((results_sims[[s]]$mean.f - mean.f)/mean.f),2,quantile, probs = 0.975))
  high95_bias_mean.f <- c(high95_bias_mean.f, high95_bias_mean.f_rep)
  
  low50_bias_mean.f_rep <- mean(apply(((results_sims[[s]]$mean.f - mean.f)/mean.f),2,quantile, probs = 0.25))
  low50_bias_mean.f <- c(low50_bias_mean.f, low50_bias_mean.f_rep)
  
  high50_bias_mean.f_rep <- mean(apply(((results_sims[[s]]$mean.f - mean.f)/mean.f),2,quantile, probs = 0.75))
  high50_bias_mean.f <- c(high50_bias_mean.f, high50_bias_mean.f_rep)
}

# Summarize mean and quantiles
bias_mean.f <- data.frame("param" = "Fecundity", "median" = mean(median_bias_mean.f), 
                          "lower95" = mean(low95_bias_mean.f), 
                          "upper95" = mean(high95_bias_mean.f), 
                          "lower50" = mean(low50_bias_mean.f), 
                          "upper50" = mean(high50_bias_mean.f), row.names = "")



### Pull the estimates together into a single dataframe for plotting
full_bias_species_df <- cbind(rbind(bias_mean.f, bias_mean.sa, bias_mean.sj), data.frame("level" = "species"))
full_bias_community_mean_df <- cbind(rbind(bias_mu.f, bias_mu.sa, bias_mu.sj), data.frame("level" = "community"))

full_bias <- rbind(full_bias_species_df, full_bias_community_mean_df)


# Plotting Figure 4B
full_bias$param <- factor(full_bias$param, levels = c("Juvenile Survival", "Adult Survival", "Fecundity"))
full_bias$param <- as.numeric(full_bias$param)

ggplot(data = full_bias, aes(group = level)) +
  geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") + 
  geom_errorbar(aes(x = param, ymin = lower95, ymax = upper95), col = "black", width = 0, position = position_dodge(0.5)) +
  geom_rect(aes(xmin=param-0.2, xmax=param+0.2, ymin=lower50, ymax=upper50, fill=level), color="black", position= position_dodge(0.5)) +
  geom_rect(aes(xmin=param-0.2, xmax=param+0.2, ymin=median, ymax=median, color = level, fill = level), position= position_dodge(0.5)) +
  scale_color_manual(values = c("black", "black")) +
  scale_fill_manual(values = c("#A9B9DB", "white")) +
  theme_bw(base_size = 16) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title=element_text(margin=margin(t=40,b=-30), hjust=0.5, size = 12)) +
  labs(x = "Parameter", y = "Relative Bias") +
  scale_x_continuous(breaks = c(1,2,3), labels = c("Juvenile \nsurvival", "Adult \nsurvival", "Fecundity")) +
  scale_y_continuous(limits = c(-0.3, 0.3))



#####################################################
# Compare simulation results for annual population growth rates truth under which the data were simulated
#####################################################

### Pulling out annual population growth rates for each species (ann.growth.rate)
median_growth_rep <- low95_growth_rep <- high95_growth_rep  <- low50_growth_rep <- high50_growth_rep <- array(dim=c(n.species, n.occasions-1, nsim))
for(s in 1:nsim){
  for(i in 1:n.species){
    median_growth_rep[i,,s] <- apply(results_sims[[s]]$ann.growth.rate[,i,],2,median)
    
    low95_growth_rep[i,,s] <- apply(results_sims[[s]]$ann.growth.rate[,i,],2,quantile, probs = 0.025)
    high95_growth_rep[i,,s] <- apply(results_sims[[s]]$ann.growth.rate[,i,],2,quantile, probs = 0.975)
    
    low50_growth_rep[i,,s] <- apply(results_sims[[s]]$ann.growth.rate[,i,],2,quantile, probs = 0.25)
    high50_growth_rep[i,,s] <- apply(results_sims[[s]]$ann.growth.rate[,i,],2,quantile, probs = 0.75)
  }
}

# Summarize mean and quantiles of species population growth rates
Growth_median_df <- apply(median_growth_rep,c(1,2),mean)[,1]
Growth_lower95_df <- apply(low95_growth_rep,c(1,2),mean)[,1]
Growth_upper95_df <- apply(high95_growth_rep,c(1,2),mean)[,1]
Growth_lower50_df <- apply(low50_growth_rep,c(1,2),mean)[,1]
Growth_upper50_df <- apply(high50_growth_rep,c(1,2),mean)[,1]

Truth_df <- apply(lam, 1, mean) # Calculate the simulated mean growth rate for each species across all years

# Format the species-level population growth rate data for plotting
Growth_full_df <- data.frame("median" = Growth_median_df, "lower95" = Growth_lower95_df, "upper95" = Growth_upper95_df, 
                             "lower50" = Growth_lower50_df, "upper50" = Growth_upper50_df,
                             "species" = c("Species 1", "Species 2", "Species 3", "Species 4", "Species 5", 
                                           "Species 6", "Species 7", "Species 8", "Species 9", "Species 10"),
                             "truth" = Truth_df)

Growth_full_df$species <- factor(Growth_full_df$species, 
                                 levels = c("Species 1", "Species 2", "Species 3", "Species 4", "Species 5", 
                                            "Species 6", "Species 7", "Species 8", "Species 9", "Species 10"))
Growth_full_df$species <- as.numeric(Growth_full_df$species)


Growth_full_df$color <- ifelse(Growth_full_df$truth <= 1, "below", "above")


# Calculate the community-level population growth rate (geometric mean of species-level population growth rates)
Community_geomeans_df <- mean(apply(median_growth_rep,3,function(x) {exp(mean(log(x)))}))
Community_geolower95_df <- mean(apply(low95_growth_rep,3,function(x) {exp(mean(log(x)))}))
Community_geoupper95_df <- mean(apply(high95_growth_rep,3,function(x) {exp(mean(log(x)))}))
Community_geolower50_df <- mean(apply(low50_growth_rep,3,function(x) {exp(mean(log(x)))}))
Community_geoupper50_df <- mean(apply(high50_growth_rep,3,function(x) {exp(mean(log(x)))}))

Truth_community_geomeans_df <- mean(apply(lam,1,function(x) {exp(mean(log(x)))})) # Calculate the simulated geometric mean growth rate for all species across all years

Growth_community <- data.frame("median" = Community_geomeans_df,"lower95" = Community_geolower95_df, "upper95" = Community_geoupper95_df, 
                               "lower50" = Community_geolower50_df, "upper50" = Community_geoupper50_df, "species" = 11, "truth" = Truth_community_geomeans_df, color = "none")

Growth_full_df <- rbind(Growth_full_df, Growth_community)


# Plotting Figure 4C
ggplot(Growth_full_df, aes(group = species)) + 
  geom_hline(yintercept = 1, linetype = 2, col = "darkgrey") + 
  geom_vline(xintercept = 10.4, linetype = 1, col = "black") + 
  geom_errorbar(aes(x = species,  ymin = lower95, ymax = upper95), col = "black", width = 0) +
  geom_rect(aes(xmin=species-0.3, xmax=species+0.3, ymin=lower50, ymax=upper50, fill=color), color="black", position= position_dodge(0.5)) +
  geom_rect(aes(xmin=species-0.3, xmax=species+0.3, ymin=median, ymax=median), lwd = 1, col = "black", fill = "black", position= position_dodge(0.5)) +
  geom_point(aes(x = species-0.4, y = truth), col = "black", fill = "black", size = 3, shape = 21) + 
  scale_fill_manual(values = c("#5B7BB4", "#FFF4B3", "#A9B9DB")) +
  scale_color_manual(values = c("#5B7BB4", "#FFF4B3", "#A9B9DB")) +
  theme_bw(base_size = 16) + 
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = 1:11, labels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "COMM")) + 
  labs(y = "Mean population growth rate")




#####################################################
# Compare simulation results for species 4 and 7 as examples
#####################################################

# First calculate the mean estimates with 95% CI for survival and fecundity of all species
# Juvenile survival
median_mean.sj <- c(); low95_mean.sj <- c(); high95_mean.sj <- c(); low50_mean.sj <- c(); high50_mean.sj <- c()

for(s in 1:nsim){
  median_mean.sj_rep <- apply(results_sims[[s]]$mean.sj,2,median)
  median_mean.sj <- rbind(median_mu.sj, median_mu.sj_rep)

  low95_mean.sj_rep <- apply(results_sims[[s]]$mean.sj,2,quantile, probs = 0.025)
  low95_mean.sj <- rbind(low95_mean.sj, low95_mean.sj_rep)
  
  high95_mean.sj_rep <- apply(results_sims[[s]]$mean.sj,2,quantile, probs = 0.975)
  high95_mean.sj <- rbind(high95_mean.sj, high95_mean.sj_rep)
  
  low50_mean.sj_rep <- apply(results_sims[[s]]$mean.sj,2,quantile, probs = 0.25)
  low50_mean.sj <- rbind(low50_mean.sj, low50_mean.sj_rep)
  
  high50_mean.sj_rep <- apply(results_sims[[s]]$mean.sj,2,quantile, probs = 0.75)
  high50_mean.sj <- rbind(high50_mean.sj, high50_mean.sj_rep)
}

# Summarize mean and quantiles
estimated_mean.sj <- data.frame("param" = "Juvenile Survival", "species" = 1:nspecies,
                                "median" = apply(median_mean.sj,2,mean), 
                                "lower95" = apply(low95_mean.sj,2,mean), 
                                "upper95" = apply(high95_mean.sj,2,mean), 
                                "lower50" = apply(low50_mean.sj,2,mean), 
                                "upper50" = apply(high50_mean.sj,2,mean))

  
# Adult survival
median_mean.sa <- c(); low95_mean.sa <- c(); high95_mean.sa <- c(); low50_mean.sa <- c(); high50_mean.sa <- c()

for(s in 1:nsim){
  median_mean.sa_rep <- apply(results_sims[[s]]$mean.sa,2,median)
  median_mean.sa <- rbind(median_mean.sa, median_mean.sa_rep)
  
  low95_mean.sa_rep <- apply(results_sims[[s]]$mean.sa,2,quantile, probs = 0.025)
  low95_mean.sa <- rbind(low95_mean.sa, low95_mean.sa_rep)
  
  high95_mean.sa_rep <- apply(results_sims[[s]]$mean.sa,2,quantile, probs = 0.975)
  high95_mean.sa <- rbind(high95_mean.sa, high95_mean.sa_rep)
  
  low50_mean.sa_rep <- apply(results_sims[[s]]$mean.sa,2,quantile, probs = 0.25)
  low50_mean.sa <- rbind(low50_mean.sa, low50_mean.sa_rep)
  
  high50_mean.sa_rep <- apply(results_sims[[s]]$mean.sa,2,quantile, probs = 0.75)
  high50_mean.sa <- rbind(high50_mean.sa, high50_mean.sa_rep)
}

# Summarize mean and quantiles
estimated_mean.sa <- data.frame("param" = "Adult Survival", "species" = 1:nspecies,
                                "median" = apply(median_mean.sa,2,mean), 
                                "lower95" = apply(low95_mean.sa,2,mean), 
                                "upper95" = apply(high95_mean.sa,2,mean), 
                                "lower50" = apply(low50_mean.sa,2,mean), 
                                "upper50" = apply(high50_mean.sa,2,mean))


# Fecundity
median_mean.f <- c(); low95_mean.f <- c(); high95_mean.f <- c(); low50_mean.f <- c(); high50_mean.f <- c()

for(s in 1:nsim){
  median_mean.f_rep <- apply(results_sims[[s]]$mean.f,2,median)
  median_mean.f <- rbind(median_mean.f, median_mean.f_rep)
  
  low95_mean.f_rep <- apply(results_sims[[s]]$mean.f,2,quantile, probs = 0.025)
  low95_mean.f <- rbind(low95_mean.f, low95_mean.f_rep)
  
  high95_mean.f_rep <- apply(results_sims[[s]]$mean.f,2,quantile, probs = 0.975)
  high95_mean.f <- rbind(high95_mean.f, high95_mean.f_rep)
  
  low50_mean.f_rep <- apply(results_sims[[s]]$mean.f,2,quantile, probs = 0.25)
  low50_mean.f <- rbind(low50_mean.f, low50_mean.f_rep)
  
  high50_mean.f_rep <- apply(results_sims[[s]]$mean.f,2,quantile, probs = 0.75)
  high50_mean.f <- rbind(high50_mean.f, high50_mean.f_rep)
}

# Summarize mean and quantiles
estimated_mean.f <- data.frame("param" = "Fecundity", "species" = 1:n.species,
                                "median" = apply(median_mean.f,2,mean), 
                                "lower95" = apply(low95_mean.f,2,mean), 
                                "upper95" = apply(high95_mean.f,2,mean), 
                                "lower50" = apply(low50_mean.f,2,mean), 
                                "upper50" = apply(high50_mean.f,2,mean))

# Pull all of the information together into a single dataframe for plotting
estimated_species <- rbind(estimated_mean.sj, estimated_mean.sa, estimated_mean.f)

truth_species <- data.frame("param" = c(rep("Juvenile Survival",n.species),rep("Adult Survival",n.species),rep("Fecundity",n.species)),
                            "species" = 1:n.species,
                            "truth" = c(mean.sj, mean.sa, mean.f))

full_species_df <- left_join(estimated_species, truth_species, by = c("param", "species"))


# Pull out species 4 and species 7 as examples to highlight
subset_species_df <- full_species_df[full_species_df$species == 4 | full_species_df$species == 7,]

subset_species_df$param <- factor(subset_species_df$param, levels = c("Juvenile Survival", "Adult Survival", "Fecundity"))
subset_species_df$param <- as.numeric(subset_species_df$param)

# remove fecundity because it is the same for both species
subset_species_df <- subset_species_df[subset_species_df$param == 1 | subset_species_df$param == 2,]
subset_species_df$species <- factor(subset_species_df$species, labels = c("4" = "S4", "5" = "S7"))
  
# Plotting Figure 4D
ggplot(subset_species_df, aes(group = param)) + 
  geom_errorbar(aes(x = param,  ymin = lower95, ymax = upper95), col = "black", width = 0) +
  geom_rect(aes(xmin=param-0.2, xmax=param+0.2, ymin=lower50, ymax=upper50, fill = factor(species)), color="black", position= position_dodge(0.5)) +
  geom_rect(aes(xmin = param-0.2, xmax=param+0.2, ymin=median, ymax=median), lwd = 1, col = "black", fill = "black", position= position_dodge(0.5)) +
  geom_point(aes(x = param-0.3, y = truth), col = "black", fill = "black", size = 3, shape = 21) + 
  scale_fill_manual(values = c("#5B7BB4", "#FFF4B3", "#A9B9DB")) +
  scale_color_manual(values = c("#5B7BB4", "#FFF4B3", "#A9B9DB")) +
  facet_wrap(~species) +
  theme_bw(base_size = 16) + 
  theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        panel.spacing = unit(0, "lines")) +
  scale_x_continuous(breaks = 1:2, labels = c("Juvenile \nsurvival", "Adult \n survival"), limits = c(0.65,2.25)) + 
  labs(y = "Estimated survival")
  
