# code 

This directory contains all code for raw data analysis, fitting models, and producing the results shown in the manuscript for the bird case study. The scripts directly in the `code` directory are data prep scripts used to prepare the bird data and covariate data for fitting the model in both `spOccupancy` and `NIMBLE`. The scripts in `spOccupancy-code` and `nimble-code/` contain the specific scripts for fitting the model in `spOccupancy` and `NIMBLE`, respectively. 

1. `bird-data-prep.R`: preps the raw BBS and eBird data for analysis. The raw eBird and BBS data files are not currently on Github, but can be downloaded from the BBS and eBird websites. The output file from this script is `data/ne-data-bundle.rda`.  
2. `calculate-bioclim.R`: this script calculates the five bioclimatic variables using the `dismo` R package and the climate data obtained from PRISM. The output file from this script is `data/climate-data-full.rda`. 
3. `ppt-data-prep.R`: script that downloads and calculates precipitation data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/ppt-2016.rda` and `data/ppt-2017.rda`. 
4. `tmax-data-prep.R`: script that downloads and calculates maximum temperature data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/tmax-2016.rda` and `data/tmax-2017.rda`. 
5. `tmin-data-prep.R`: script that downloads and calculates minimum temperature data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/tmin-2016.rda` and `data/tmin-2017.rda`. 


## spOccupancy-code

1. `main-spOccupancy.R`: main script for running the integrated community occupancy model through `spOccupancy`. 
2. `spOccupancy-data-prep.R`: script that prepares the data for use in `spOccupancy`. The output file from this script is `data/spOccupancy-data.rda`. 
3. `summary-spOccupancy.R`: script that summarizes results from the model fit in `spOccupancy` and creates the figures shown in the main text and appendix for the bird case study.

## nimble-code

1. `icom-ne-birds-nimble.R`: NIMBLE code for running the integrated community occupancy model with the BBS and eBird data. 
2. `main-nimble.R`: main script for running the integrated community occupancy model through `NIMBLE`. 
3. `summary-nimble.R`: script that summarizes results from the model fit in NIMBLE and creates the figures shown in the main text and appendix for the bird case study. Note that the figures used in the main text were created using the model fit in `spOccupancy`. 

