# code

This directory contains all code to simulate the data, fit models, and produce the results shown in the manuscript for the multispecies integrated population model case study.

1. `Mutispecies_IPM_Simulations.R`: This R script contains all of the necessary code to simulate 100 independent annual population count (i.e., census), capture-recapture and productivity datasets and analyze these data using the multispecies integrated population model described in the main text.
2. `Multspecies_Simulations.txt`: This text file contains the mutlispecies integrated population model code for JAGS implementation in R using the jagsUI package.
3. `Multispecies_IPM_Results.R`: This R script contains all of the necessary code to analyze the simulation results (e.g., calculate relative bias of estimated parameters at both the species- and community-levels) and recreate Figure 4, panels b-d.
