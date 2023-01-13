# ppt-data-prep.R: this script extracts precipitation data from PRISM
#                  at the locations of the 22,513 5x5km grid cells in the 
#                  Northeastern US for use in both fitting and predicting
#                  species-specific occupancy and richness using an 
#                  integrated community occupancy model. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(sf)
library(stars)
library(prism)

# Read in data to get prediction grid -------------------------------------
load("data/ne-data-bundle.rda")
coords.sf.albers <- st_centroid(grid.ne)
coords.albers <- st_coordinates(st_centroid(grid.ne))

# Set PRISM directory -----------------------------------------------------
prism_set_dl_dir("data/prism-climate")

# Extract data at BBS routes ----------------------------------------------
# 12 months in a year. 
ppt.vals <- matrix(NA, nrow = nrow(coords.albers), ncol = 12)
# NOTE: can change this to extract values for different years. For the case 
#       study in the associated manuscript, we extracted the values for 2016
#       and 2017, and then used data spanning from June 2016 - May 2017.
curr.year <- 2016
get_prism_monthlys(type = "ppt", year = curr.year, mon = 1:12, keepZip = FALSE)
# Get file name to data of interest
# Loop through all the months
for (i in 1:ncol(ppt.vals)) {
  ppt.curr.path <- prism_archive_subset("ppt", "monthly", years = curr.year, mon = i)
  # Get absolute file path
  ppt.curr.abs <- pd_to_file(ppt.curr.path)
  # Download the raster
  ppt.curr <- read_stars(ppt.curr.abs)
  # Get coordinates
  coords.proj <- coords.sf.albers %>%
    st_transform(crs = st_crs(ppt.curr))
  # Extract precipitation value at the center of each 5x5km grid cell. 
  tmp <- st_extract(ppt.curr, at = coords.proj)
  ppt.vals[, i] <- tmp[[1]]
}

# Save values to hard drive -----------------------------------------------
save(ppt.vals, file = paste("data/ppt-", curr.year, ".rda", sep = ''))

