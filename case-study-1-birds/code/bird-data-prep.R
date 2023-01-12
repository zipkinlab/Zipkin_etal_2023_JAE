# bird-data-prep.R: this script prepares the eBird and BBS data for analysis
#                   in an ICOM framework. The raw data files necessary for
#                   this script are not included on GitHub, but the resulting
#                   object from this script "data/ne-data-bundle.rda" is. 
# Author: Jeffrey W. Doser
rm(list = ls())
library(tidyverse)
library(lubridate)
library(sf)
# For reading in eBird data
library(auk)
# For US states data as sf objects
library(spData)
# For linking bird codes with species info. 
library(wildlifeR)
# For some exploratory plots
library(ggnewscale)
# For grabbing NED and NLCD data
library(FedData)
# For linking with NED and NLCD. 
library(sp)
library(rgdal)
library(raster)

# Get BBS Data ------------------------------------------------------------
# Read in BBS Data from NE states. Note these raw data are not included on GitHub,
# but can be downloaded from the BBS website. 
bbs.full.dat <- list.files(path = "data/BBS/ne-states", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows()
# Get associated route data
route.dat <- read.csv("data/BBS/routes.csv")
# Get associated weather data
weather.dat <- read.csv("data/BBS/weather.csv")
# Only grab weather data from NE states
weather.dat <- weather.dat %>%
  filter(StateNum %in% unique(bbs.full.dat$StateNum))
# Join BBS data with route data
bbs.dat <- left_join(bbs.full.dat, route.dat, by = c('Route', 'CountryNum', 'StateNum'))
# Only grab data from (2017)
bbs.2017 <- bbs.dat %>%
  filter(Year == 2017)
# Select columns of interest
bbs.dat <- bbs.dat %>%
  dplyr::select(Route, Year, AOU, starts_with('Stop'), Latitude, Longitude)
# Select columns of interest
bbs.2017 <- bbs.2017 %>%
  dplyr::select(RouteDataID, Latitude, Longitude, AOU, starts_with("Count")) %>%
  dplyr::select(-CountryNum)
# Fill in implicit zeros
bbs.2017 <- bbs.2017 %>%
  complete(AOU, nesting(RouteDataID, Latitude, Longitude))
# Replace NAs with 0s for all columns at once.
bbs.2017 <- bbs.2017 %>%
  replace(is.na(.), 0)
# Link AOU codes with species info. The AOU_species_codes object is 
# from wildlifeR. 
aou.info <- AOU_species_codes %>%
  mutate(AOU = spp.num)
bbs.2017 <- left_join(bbs.2017, aou.info, by = c("AOU"))

# Read in guild information -----------------------------------------------
guild.info <- read.table("data/guild-info.txt")
interior.for <- guild.info %>%
  dplyr::filter(Response_Guild == 'InteriorForestObligate') %>%
  pull(AOU_Code)

# Get detection covariates in proper format -------------------------------
weather.dat <- weather.dat %>% 
  unite('date', sep = '-', Year, Month, Day, remove = FALSE)
weather.dat$date <- as.Date(weather.dat$date, tz = "America/New_York")
# Get julian date of each survey
weather.dat$julian <- as.numeric(format(weather.dat$date, '%j'))
# Join detection covariates with observations
bbs.2017 <- left_join(bbs.2017, weather.dat, by = c('RouteDataID'))

# Filter out for only interior forest obligates ---------------------------
bbs.2017 <- bbs.2017 %>%
  separate(spp, c('orderName', 'sppName'), sep = ' ', remove = FALSE) %>%
  mutate(orderName = tolower(orderName))
bbs.2017 <- bbs.2017 %>%
  filter(alpha.code %in% interior.for)

# Number of routes in each state
bbs.2017 %>%
  group_by(StateNum) %>%
  summarize(unique.routes = n_distinct(RouteDataID))
# Get observation data in long format
y.bbs.long <- bbs.2017 %>%
  dplyr::select(alpha.code, starts_with('Count'), Route, Longitude, Latitude, 
		ObsN, julian) %>%
  mutate(sp = as.character(alpha.code)) %>%
  dplyr::select(-alpha.code, -CountryNum)
n.route <- n_distinct(y.bbs.long$Route)
n.sp <- n_distinct(y.bbs.long$sp)
n.visit <- 5
y.bbs.long <- y.bbs.long %>%
  mutate(across(starts_with('Count'), function(a) {ifelse(a > 0, 1, 0)}))
# Add summary binomial count. 
y.bbs.long$binom <- apply(y.bbs.long[, c('Count10', 'Count20', 'Count30', 
					 'Count40', 'Count50')], 1, sum)

# Get eBird data ----------------------------------------------------------
# States we're using: CT, DE, MA, MD, ME, NH, NJ, NY, PA, RI, VT
# Takes a few minutes to load these all in. 
# These data can be obtained from eBird, but are not included on the GitHub.
ebird.CT <- read_ebd('data/eBird-test/ebd_US-CT_201705_201707_relJan-2022.txt')
ebird.DE <- read_ebd('data/eBird-test/ebd_US-DE_201705_201707_relJan-2022.txt')
ebird.ME <- read_ebd('data/eBird-test/ebd_US-ME_201705_201707_relJan-2022.txt')
ebird.MA <- read_ebd('data/eBird-test/ebd_US-MA_201705_201707_relJan-2022.txt')
ebird.MD <- read_ebd('data/eBird-test/ebd_US-MD_201705_201707_relJan-2022.txt')
ebird.NH <- read_ebd('data/eBird-test/ebd_US-NH_201705_201707_relJan-2022.txt')
ebird.NJ <- read_ebd('data/eBird-test/ebd_US-NJ_201705_201707_relJan-2022.txt')
ebird.NY <- read_ebd('data/eBird-test/ebd_US-NY_201705_201707_relJan-2022.txt')
ebird.PA <- read_ebd('data/eBird-test/ebd_US-PA_201705_201707_relSep-2021.txt')
ebird.RI <- read_ebd('data/eBird-test/ebd_US-RI_201705_201707_relJan-2022.txt')
ebird.VT <- read_ebd('data/eBird-test/ebd_US-VT_201705_201707_relJan-2022.txt')
# Bind everything together
ebird.tmp <- rbind(ebird.CT, ebird.DE, ebird.ME, ebird.MA, ebird.MD, ebird.NH, 
		   ebird.NJ, ebird.NY, ebird.RI, ebird.VT)
# Bind everything together. 
ebird.ne <- full_join(ebird.PA, ebird.tmp)
# Filter for records in May and June
ebird.ne <- ebird.ne %>%
  mutate(month = month(observation_date)) %>%
  filter(month %in% c(5, 6))

# Filter eBird data set to only those that have auxiliary information
# and not extreme information
ebird.dat <- ebird.ne %>%
  filter(month %in% c(5, 6), # May and June "closure" is assumed
	 protocol_type %in% c('Stationary', 'Traveling'), # stationary or traveling protocols
	 all_species_reported == TRUE, # complete checklists
	 number_observers < 10, # restrict number of observers
	 duration_minutes < 300) %>% # remove super long surveys
  mutate(effort_distance_km = ifelse(protocol_type == 'Stationary', 0,
				     effort_distance_km)) %>%
  filter(effort_distance_km < 5) # restrict spatial scale. This should inform the size of grid cell.

# Select columns of interest, create new columns of interest
ebird.dat <- ebird.dat %>%
  dplyr::select(checklist_id, common_name, observation_count, latitude, longitude,
	 observation_date, time_observations_started, observer_id,
	 duration_minutes, effort_distance_km, month,
	 number_observers, all_species_reported) %>%
  mutate(year = year(observation_date),
	 week = week(observation_date))

# Grid up study area ------------------------------------------------------
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
ne.states <- usa %>% 
  filter(ID %in% c('connecticut', 'delaware', 'maine', 'maryland', 
		     'massachusetts', 'new hampshire', 'new jersey', 
		     'new york', 'pennsylvania', 'rhode island', 
		     'vermont'))
# Can change the grid cell size as needed. Grid cell size here is 5km x 5km.  
grid <- ne.states %>%
  st_make_grid(cellsize = c(5000, 5000)) # in meters

# Connnect eBird to the gridded area --------------------------------------
ebird.sf <- st_as_sf(ebird.dat, coords = c('longitude', 'latitude'), 
		     crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
# Convert ebird locations to albers equal area
ebird.sf <- ebird.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Convert grid object to sf data frame
grid.sf <- grid %>% st_as_sf()
grid.ne <- grid.sf %>%
  st_intersection(st_buffer(ne.states, 5000))
# Get cell each eBird observation falls in
cell.eBird <- st_contains(grid.ne, ebird.sf)
# Get number of eBird observations in each cell
cell.counts <- sapply(cell.eBird, length)
# Number of cells
n.cells <- length(cell.counts)
# Cell corresponding to each value in cells.eBird.vec
cell.by.val <- unlist(sapply(1:n.cells, function(a) {rep(a, cell.counts[a])}))
# Rows in eBird that are in the grid. 
cells.eBird.vec <- unlist(cell.eBird)
grid.ne$counts <- cell.counts
# Assign each eBird observation to a cell
ebird.dat$cell <- NA
ebird.dat$cell[cells.eBird.vec] <- cell.by.val
# Remove records not in the study area
ebird.dat <- ebird.dat %>%
  filter(!is.na(cell))

# Connect BBS data to the gridded area ------------------------------------
# Convert BBS data to sf object.
bbs.sf <- st_as_sf(y.bbs.long, coords = c('Longitude', 'Latitude'),
	           crs = "+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
# Convert to albers equal area
bbs.sf <- bbs.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")

# Get cell each BBS data point falls in
cell.bbs <- st_contains(grid.ne, bbs.sf)
bbs.sf$cell <- NA
for (i in 1:n.cells) {
  if (length(cell.bbs[[i]]) > 0) {
    bbs.sf$cell[cell.bbs[[i]]] <- i    
  }
}
# Total possible counts (denominator of binomial)
bbs.sf$denom <- 5

# BBS cells regardless of species
bbs.sf.all <- bbs.sf %>%
  group_by(Route) %>% 
  summarize(all.count = sum(binom))

# Get elevation and forest cover within each grid cell --------------------
# Split up grabbing the covariates by each state
coords.sp <- grid.ne %>% 
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame()
J <- nrow(coords.sp)
coordinates(coords.sp) <- ~X + Y
proj4string(coords.sp) <- '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs'
# Loop through all the sites
vals <- split(1:J, ceiling(seq_along(1:J)/30))
# Average elevation within a cell
elev.vals <- rep(0, J)
# Proportion of forest cover within the cell
for.vals <- rep(0, J)
# Get proportion of forest
props <- function(a, na.rm = TRUE) {
  my.sum <- sum(!is.na(a))	
  prop.for <- sum(a %in% c(41, 42, 43), na.rm = na.rm) / my.sum
  return(prop.for)
}
ned.dat <- list()
nlcd.dat <- list()
# This takes a day or so to calculate, so I have commented it out and saved
# the resulting object in the file "data/spatial-covariates.rda", which is read in 
# below. 
# for (i in 1:length(vals)) {
#   print(paste("Currently on iteration ", i, " out of ", length(vals), sep = ''))
#   ned.dat[[i]] <- get_ned(template = coords.sp[vals[[i]], ], label = paste('birds', i))
#   nlcd.dat[[i]] <- get_nlcd(template = coords.sp[vals[[i]], ], label = paste('birds', i), year = 2016)
#   # Slightly faster. 
#   elev.vals[vals[[i]]] <- raster::extract(ned.dat[[i]], coords.sp[vals[[i]], ])
#   for.vals[vals[[i]]] <- raster::extract(nlcd.dat[[i]], coords.sp[vals[[i]], ], buffer = 2500, fun = props)
# }

# Save covariates so you don't have to rerun all the time.
# save(elev.vals, for.vals, file = "data/spatial-covariates.rda")

# Load in the spatial covariates
load("data/spatial-covariates.rda")
# This brings in elev.vals and for.vals
# 
# # Get BBS data in analysis format -----------------------------------------
y.bbs <- bbs.sf %>% 
  as.data.frame() %>%
  dplyr::select(binom, Count10, Count20, Count30, Count40, 
		Count50, sp, cell, obs = ObsN, julian) %>%
  arrange(sp)
# Keep data in long format. This is how the data will be formatted for use in 
# NIMBLE. 
 
# Get eBird Data in analysis format ---------------------------------------
# Number of species
N <- n_distinct(y.bbs$sp)
# Julian Date. 
ebird.dat$julian <- as.numeric(format(ebird.dat$observation_date, '%j'))
# Time survey started in minutes since midnight.  
ebird.dat$time.started <- ebird.dat$time_observations_started %>%
  hms() %>%
  period_to_seconds() / 60

# Spatio-temporal subset the eBird data -----------------------------------
# Select one checklist per week/cell combination to reduce the effects
# of preferential sampling on the resulting estimates
ebird.dat <- ebird.dat %>%
  unite(col = "weekCell", c(week, cell), remove = FALSE)

n.combos <- n_distinct(ebird.dat$weekCell)
week.cell.combos <- unique(ebird.dat$weekCell)
curr.indx <- vector(mode = 'list', length = n.combos)
# Takes a few minutes
for (i in 1:n.combos) {
  print(i)
  tmp <- ebird.dat %>%
    filter(weekCell == week.cell.combos[i])
  check.tmp <- unique(tmp$checklist_id)
  curr.val <- sample(1:length(check.tmp), 1)
  curr.indx[[i]] <- which(ebird.dat$checklist_id == check.tmp[curr.val])
} # i

# Filter out eBird data
ebird.filt.dat <- ebird.dat[unlist(curr.indx), ]

checklist.data <- ebird.filt.dat %>%
  group_by(cell) %>%
  summarize(n.lists = n_distinct(checklist_id)) %>%
  ungroup()

# Get unique checklist id in each cell. 
checklist.by.cell <- ebird.filt.dat %>%
  group_by(cell) %>%
  summarize(listID = unique(checklist_id)) %>%
  ungroup()

checklist.by.cell <- checklist.by.cell %>%
  slice(rep(row_number(), N))

# Total number of data set rows (species x rep/site combos)
n.ebird <- nrow(checklist.by.cell)

# Species in alphabetical order of 4-letter code
sp <- unique(y.bbs$sp)
sp.names <- aou.info %>%
  filter(alpha.code %in% sp) %>%
  arrange(alpha.code) %>%
  pull(name) %>%
  as.character()

# Create data frame to hold data and auxiliary information
ebird.df <- checklist.by.cell
ebird.df$sp <- rep(sp, each = sum(checklist.data$n.lists))
ebird.df$day <- NA
ebird.df$time <- NA
ebird.df$length <- NA
ebird.df$dist <- NA
ebird.df$obsv <- NA
ebird.df$week <- NA
ebird.df$y <- NA

# ebird data are ordered by species, then cell within species
# Number of lists
n.lists <- sum(checklist.data$n.lists)

# Get all eBird data for analysis. This takes a few minutes to run. 
for (i in 1:(n.ebird/N)) {
  if(i %% 1000 == 0) print(i)
  tmp <- ebird.filt.dat %>%
    filter(cell == checklist.by.cell$cell[i],
	   checklist_id == checklist.by.cell$listID[i])
  indx <- seq(i, n.ebird, by = n.lists)
  ebird.df$day[indx] <- max(tmp$julian)
  ebird.df$time[indx] <- max(tmp$time.started)
  ebird.df$week[indx] <- max(tmp$week)
  ebird.df$length[indx] <- max(tmp$duration_minutes)
  ebird.df$dist[indx] <- max(tmp$effort_distance_km)
  ebird.df$obsv[indx] <- max(tmp$number_observers)
  tmp.2 <- tmp %>% filter(str_to_upper(common_name) %in% str_to_upper(sp.names))
  curr.indices <- which(str_to_upper(sp.names) %in% str_to_upper(tmp.2$common_name))
  curr.vals <- ifelse(1:N %in% curr.indices, 1, 0)
  ebird.df$y[indx] <- curr.vals
}

# Note that in this script, we have extracted the eBird data for all 
# dates between May and July. We further restrict this to the first three
# weeks in June during subsequent scripts. 

# Save results ------------------------------------------------------------
occ.covs <- data.frame(elev = elev.vals, 
		       pf = for.vals)
# Also saving bbs.sf and grid.ne for making plots. 
save(y.bbs, ebird.df, occ.covs, grid.ne, bbs.sf, file = "data/ne-data-bundle.rda")

