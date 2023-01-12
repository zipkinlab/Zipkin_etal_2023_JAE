# code 

1. `bird-data-prep.R`: preps the raw BBS and eBird data for analysis. The raw eBird and BBS data files are not currently on Github, but can be downloaded from the BBS and eBird websites. The output file from this script is `data/ne-data-bundle.rda`.  
2. `calculate-bioclim.R`: this script calculates the five bioclimatic variables using the `dismo` R package and the climate data obtained from PRISM. The output file from this script is `data/climate-data-full.rda`. 
3. `icom-ne-birds-nimble.R`: NIMBLE code for running the integrated community occupancy model with the BBS and eBird data. 
4. `main-nimble.R`: main script for running the integrated community occupancy model through `NIMBLE`. 
5. `main-spOccupancy.R`: main script for running the integrated community occupancy model through `spOccupancy`. 
6. `ppt-data-prep.R`: script that downloads and calculates precipitation data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/ppt-2016.rda` and `data/ppt-2017.rda`. 
7. `spOccupancy-data-prep.R`: script that prepares the data for use in `spOccupancy`. The output file from this script is `data/spOccupancy-data.rda`. 
8. `summary-nimble.R`: script that summarizes results from the model fit in NIMBLE and creates the figures shown in the main text and appendix for the bird case study. Note that the figures used in the main text were created using the model fit in `spOccupancy`. 
9. `summary-spOccupancy.R`: script that summarizes results from the model fit in `spOccupancy` and creates the figures shown in the main text and appendix for the bird case study.
8. `tmax-data-prep.R`: script that downloads and calculates maximum temperature data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/tmax-2016.rda` and `data/tmax-2017.rda`. 
9. `tmin-data-prep.R`: script that downloads and calculates minimum temperature data from PRISM for use when calculating the five bioclimatic variables. Output files are `data/tmin-2016.rda` and `data/tmin-2017.rda`. 
