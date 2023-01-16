#####################################################
## Integrated Community Models - Case Study Example #
## Multi-species Integrated Population Models       #
## Created by: Courtney L. Davis                    #
## Last modified: January 3, 2023                   #
#####################################################

# Note: This code uses functions from the R package "IPMbook" to simulate count, productivity, and capture-recapture data
# These functions were originally created and are managed by Mike Meredith (https://github.com/mikemeredith/IPM_code) 
# in association with the "Integrated Population Models" book written by Michael Schaub and Marc Kery.


#####################################################
# Load packages
#####################################################
library(here)
library(IPMbook)
library(jagsUI)


#####################################################
# Definition of simulation parameters
#####################################################

# Number of simulations
nsim <- 100    # ~90 minutes
#nsim <- 5    # ~6 minutes for testing

# Number of years or occasions
n.occasions <- 10

# Number of species
n.species <- 10

# Observation error for the population survey - keep relatively constant for all species
sigma <- rnorm(n.species, 10, 2)

# Community-level mean parameters - probability scale
mu.sj <- 0.30
mu.sa <- 0.55
#mu.f.lower <- 2
#mu.f.upper <- 4
mu.f <- 3
mu.p <- 0.6

# Community-level variance parameters
sd.lsj <- 0.3   # logit scale
sd.lsa <- 0.3   # logit scale
sd.lp <- 0.3    # logit scale
sd.f <- 0.1     # log scale

# Species-specific effects - convert probability scale mean values to logit scale
set.seed(100)
lsj <- rnorm(n.species, qlogis(mu.sj), sd.lsj)
lsa <- rnorm(n.species, qlogis(mu.sa), sd.lsa)
lp <- rnorm(n.species, qlogis(mu.p), sd.lp)

# Convert logit scale species-specific effects to probability scale
mean.sj <- plogis(lsj)
mean.sa <- plogis(lsa)
mean.p <- plogis(lp)

mean.ef <- rnorm(n.species, log(mu.f), 0.1)
mean.f <- round(exp(mean.ef))

# Probability to find a brood whose reproductive output is recorded
pprod <- 0.9

# Set initial population size per species and age class
N <- array(NA, dim = c(n.species,2,1))
N[,1,1] <- round(rnorm(n.species, 30, 5))
N[,2,1] <- round(rnorm(n.species, 30, 10))


# Create empty matrices to store data and results
population_data_1 <- list(); population_data_2 <- list(); population_data_3 <- list()

counts <- matrix(NA, n.species, n.occasions)

productivity <- list()
sJ <- c(); nJ <- c()

captures <- list(); marr <- list(); 
marr.j <- array(NA, dim = c(n.species, n.occasions-1, n.occasions))
marr.a <- array(NA, dim = c(n.species, n.occasions-1, n.occasions))
rel.j <- matrix(NA, n.species, n.occasions-1)
rel.a <- matrix(NA, n.species, n.occasions-1)
lam <- matrix(NA, n.species, nsim)

results <- array(NA, dim=c(279, 11, nsim))
results_sims <- list()


#####################################################
# Settings for JAGS
#####################################################

# Initial values
inits.ipm <- function(){list()}

# Parameters monitored
parameters.ipm <- c("mean.sj", "mu.sj", "lsj", "sd.lsj", "mean.sa", "mu.sa", "lsa", "sd.lsa", "mean.p", "mu.p", "lsp", "sd.lp", "mean.f", "mu.f", 
                    "sd.f", "sigma", "ann.growth.rate", "Ntot")

# MCMC settings
ni <- 60000; nb <- 15000; nc <- 3; nt <- 5; na <- 5000


#####################################################
# Conduct the simulations
#####################################################

system.time(
for(s in 1:nsim){
  set.seed(s)
  print(paste('Currently on simulation ', s, ' out of ', nsim, sep = ''))
  
  for (i in 1:n.species){
    # Create three independent populations
    population_data_1[[i]] <- simPop(Ni = N[i,,], phi = c(mean.sj[i], mean.sa[i]), f = mean.f[i], nYears = n.occasions)
    population_data_2[[i]] <- simPop(Ni = N[i,,], phi = c(mean.sj[i], mean.sa[i]), f = mean.f[i], nYears = n.occasions)
    population_data_3[[i]] <- simPop(Ni = N[i,,], phi = c(mean.sj[i], mean.sa[i]), f = mean.f[i], nYears = n.occasions)
  
    # Simulate the count survey data  
    counts[i,] <- simCountNorm(N=population_data_1[[i]]$totB, sigma=sigma[i])$count
  
    lam[i,s] <- mean(population_data_1[[i]]$totB[-1] / population_data_1[[i]]$totB[-n.occasions]) # record the mean population growth rate
    
    # Simulate the productivity data
    productivity[[i]] <- simProd(reprod=population_data_2[[i]]$reprod, pInclude=pprod, verbose = FALSE)
    sJ[i] <- colSums(productivity[[i]]$prod.agg)[1]
    nJ[i] <- colSums(productivity[[i]]$prod.agg)[2]
  
    # Simulate the capture-recapture data
    captures[[i]] <- simCapHist(state=population_data_3[[i]]$state, cap=mean.p[i], recap=mean.p[i], maxAge=2, verbose = FALSE)
  
    # translate capture histories into m-arrays
    marr[[i]] <- marrayAge(captures[[i]]$ch, captures[[i]]$age)
  
    marr.j[i,,] <- marr[[i]][,,1]
    marr.a[i,,] <- marr[[i]][,,2]
  
    rel.j[i,] <- rowSums(marr[[i]][,,1])
    rel.a[i,] <- rowSums(marr[[i]][,,2])
  } # species


  # Package the data
  jags.data.ipm <- list(marr.j=marr.j, marr.a=marr.a, n.occasions = n.occasions,
        rel.j=rel.j, rel.a=rel.a, sJ = sJ, nJ = nJ, C=counts, n.species = n.species)

  # Fit the model
  m1 <- try(jags(jags.data.ipm, inits.ipm, parameters.ipm, here("code","Multispecies_Simulations.txt"),
               n.iter=ni, n.burnin=nb, n.chains=nc, n.thin=nt, n.adapt=na, parallel=TRUE))
  if(!inherits(m1, "try-error"))
    results[,,s] <- m1$summary
    results_sims[[s]] <- m1$sims.list
  
  print(s)
  } ) # simulations


# Save simulation results
rownames(results) <- rownames(m1$summary)
colnames(results) <- colnames(m1$summary)

save(results, results_sims, mu.sj, mu.sa, mu.f, mu.p, sd.lsj, sd.lsa, sd.f, sd.lp, lsj, lsa, lp, mean.sj, mean.sa, mean.p, mean.f, lam, n.species, n.occasions, n.sim, file = here("Results","Simulation_Results.Rdata"))
