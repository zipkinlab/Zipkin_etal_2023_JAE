# data

1. `climate-data-full.rda`: five bioclimatic variables calculated for each 5 x 5 km grid cell across the Northeastern US for use as covariates in the integrated community occupancy model.
2. `guild-info.txt`: table that contains information on different structural, functional, and compositional bird guilds. This classification comes from [O'Connell et al. 2000](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/1051-0761%282000%29010%5B1706%3ABGAIOE%5D2.0.CO%3B2). 
3. `ne-data-bundle.rda`: eBird, BBS, forest cover, and elevation data across the Northeastern U.S. These data are in a format suitable for modeling with NIMBLE.
4. `ppt-2016.rda`, `ppt-2017.rda`, `tmax-2016.rda`, `tmax-2017.rda`, `tmin-2016.rda`, `tmin-2015.rda`: PRISM climate data from 2016 to 2017 used to generate bioclimatic data. The final five bioclimatic variables used to fit the model are in the `climate-data-full.rda` file. 
5. `spatial-covariates.rda`: forest cover and elevation data. 
6. `spOccupancy-data.rda`: a list of the BBS data, eBird data, occupancy covariates, detection covariates, and other information in the required format for fitting the integrated community occupancy model with the `spOccupancy` R package. 
