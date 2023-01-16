# Case Study 3: Estimating species- and community-level demographic rates and population growth 

This repository contains all code, data, and documents related to case study 3 in which we developed a multi-species integrated population model to estimate demographic rates and annual population dynamics of 10 species over a 10-year time period using a data simulation approach. The directory is split into the following structure. Each sub-directory contains its own README file with details on each of the files in that specific sub-directory.

+ `code`: contains all scripts for simulating the data, running the multispecies integrated population model in `JAGS`, and processing of the results.
+ `example_data`: contains a set of example data files using values that were generated during the data simulation. These files can be viewed as example input files, but the code to format these files for analysis is not provided directly. The full set of data generated and analyzed in this case study can be recreated by keeping the set.seed functions as specified in `Multispecies_IPM_Simulations.R`.
