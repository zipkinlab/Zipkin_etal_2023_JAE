# butterfly-data-prep.R: this script prepares and cleans the raw butterfly data
#                        into a standardized format for use in the integrated
#                        community models. 
# Authors: Wendy Leuenberger (vast majority) and Jeffrey W. Doser
# NOTE: the data from the North American Butterfly Associate (NABA) are not publicly 
#       available, so this script will not run successfully using the resources
#       provided on the Github Repository. If you would like to use the NABA
#       data, please contact NABA: https://www.naba.org/.
# Load packages -----------------------------------------------------------
# Clear workspace
rm(list = ls())
library(magrittr)
library(tidyverse)
library(ggthemes)
library(tidybayes)
# For county/state classification and/or maps
library(sf)
library(sp)
library(maps)
library(maptools)
# Different color palettes
library(RColorBrewer)
library(wesanderson)
library(ggsci)


# Subroutines -------------------------------------------------------------
# Coordinates to county ---------------
# Function from: 
#     https://stackoverflow.com/questions/13316185/r-convert-zipcode-or-lat-long-to-county
# The single argument to this function, pointsDF, is a data.frame in which:
#   - column 1 contains the longitude in degrees (negative in the US)
#   - column 2 contains the latitude in degrees
latlong2county <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per county
  counties <- map('county', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(counties$names, ":"), function(x) x[1])
  counties_sp <- map2SpatialPolygons(counties, IDs=IDs,
                                     proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, counties_sp)
  
  # Return the county names of the Polygons object containing each point
  countyNames <- sapply(counties_sp@polygons, function(x) x@ID)
  countyNames[indices]
}

# Fill in the state if NA -------------
# Some States are NA. Use this function to fill in the state if
# the county is one of the ones recorded in the state
# Won't work if the survey with state = NA is the only one in that
# county or if the survey is in a different state
FillState <- function(Data, StateCode){
  StateCounties <- Data %>% 
    # Grab just rows recorded as in the particular state
    filter(State == StateCode) %>%
    # Grab all the counties in that state with recorded surveys
    use_series(state.county) %>% 
    # Take unique values so that it is simpler
    unique
  # Find the unique values and assign them to that state if the 
  # county has other recorded surveys in that county/state
  Data$State[is.na(Data$State) &
               Data$state.county %in% StateCounties] <- StateCode
  return(Data)
}

# Standard error ----------------------
se <- function(x, na.rm = FALSE){ 
  sqrt(var(x, na.rm = na.rm) / length(x))
}


# Load data ---------------------------------------------------------------
# Big csv file with information from IL, IA, OH, NFJ, and MI
# Takes a few seconds
Butterfly <- read_csv('data/masterButterflyData.csv')
Bfly <- Butterfly  

# Species names/information
Names <- read_csv('data/BflyNames.csv')
# Ordinal day and week relative to March 1
DayWeek <- read_csv('data/DayWeek.csv')

# Restrict analysis to the following programs and case study years/months
CSPrograms <- c('Illinois Butterfly Monitoring Network', 
                'Iowa Butterfly Survey Network',
                'Michigan Butterfly Network', 'NFJ', 'OhioLeps')
CSYears <- 2008:2017
CSMonths <- 6:7

  


# Set CSStates ------------------------
# Some programs have surveys in other states
# TRUE = Only states with the programs
# FALSE = Include surveys in other states from approved programs
StatesStrict <- FALSE 

if(StatesStrict == TRUE){
  CSStates <- c('OH', 'IL', 'IA', 'MI', NA)
} else {
  CSStates <- c('OH', 'IL', 'IN', 'IA', 'MI', 'WI', NA) 
}

# Select case study species -----------
# 10 resident, multivoltine, widespread, and common
CaseStudySpp <- Names %>% 
  filter(UMD_Code %in% c('ANCNUM', 'CELLAD', 'EVECOM', 'EPACLA',
                         'PAPGLA', 'PAPPOL', 'PHYTHA', 'PIERAP',
                         'POLIPEC', 'VANVIR'))
CaseStudySpp %<>% rename('Code' = 'UMD_Code')


# Summarize data by Subspecies or by Species?
SummaryLevel <- 'Species'
# SummaryLevel <- 'SubSpecies'

# Just a single species code to help with some data processing. 
OneSpp <- 'EVECOM'  # Eastern-tailed blue


# Take a look -------------------------
head(Bfly)
Bfly %>% glimpse
Names %>% glimpse
DayWeek %>% head
Bfly$Program %>% table

# Data Processing ---------------------------------------------------------
# Remove unnecessary columns ----------
Bfly[,c('Completed', '...1', 'Observer')] <- NULL

# Remove NABACode if they are all NAs
if(all(is.na(Bfly$NABACode))){
  Bfly$NABACode <- NULL
}

# Remove duplicated surveys -----------
Bfly %<>% 
  mutate(Check = str_extract(GUEventID, '[:digit:]+$'),
         Same = ifelse(Check == UniqueSurveyID, TRUE, FALSE))
Bfly %<>% 
  filter(Same == TRUE | is.na(Same))

# Add Week values relative to March ---
DayWeekMerge <- DayWeek %>% 
  select(`Day of month`, Month, OrdinalDay, Mar2FebWeeks) %>% 
  rename(Day = `Day of month`)
Bfly %<>% 
  left_join(DayWeekMerge) 

# Change Ohio to OH
Bfly %<>% 
  mutate(State = str_replace(State, 'Ohio', 'OH'))

# Add County names --------------------
Bfly %<>% 
  mutate(state.county = 
           latlong2county(data.frame(GULongitude, GULatitude)))

# There are some NA's. I used google maps to "ground-truth" these points
Bfly %>% 
  filter(is.na(state.county)) %>% 
  select(GULocationID, Program, state.county, GULatitude, GULongitude) %>% 
  unique %>% dim
# There are 1126 NA's as of now. Not all are in the case study though
# Ground truth just for case study/subset in use

# Make indicator variables
# May or may not be necessary
# Right now, state.county.ind and county.ind are the same. 
Bfly %<>% 
  mutate(state.county.ind = state.county %>% as.factor %>% as.numeric,
         county.ind = state.county %>% as.factor %>% as.numeric,
         site.ind = GULocationID %>% as.factor %>% as.numeric) 

# Check states for is.na(State)
Bfly %>% 
  filter(is.na(State)) %>% 
  select(Program, state.county) %>% 
  unique 

# Fill in the values assuming there are surveys in that county and
# the survey is within the expected state
Bfly %<>% 
  FillState(., StateCode = 'IL') %>% 
  FillState(., StateCode = 'MI') %>% 
  FillState(., StateCode = 'IA') %>% 
  FillState(., StateCode = 'OH') 
  
# Calculate if there are still NAs in state
StateNA <- is.na(Bfly$State) %>% sum == 0 

# Throw a stop if there are still NAs in states
if(StateNA != TRUE){
  stop('Some states are NAs')
}

# Filter for CaseStudy ----------------
CaseStudyData <- Bfly %>% 
  filter(Program %in% CSPrograms,
         Year %in% CSYears,
         Month %in% CSMonths,
         State %in% CSStates) 

# Check for NA's in state.county
StateCountyNA <- CaseStudyData %>% 
  filter(is.na(state.county)) %>% 
  select(GULocationID, Program, state.county, GULatitude, GULongitude) %>% 
  unique

# GULocations and the ground truthed values
GoogleMapGroundTruth <- 
  data.frame(GULocationID = c("1095", 
                              "OH-46", "OH-10", 
                              "IL-1322", "IL-6219", "IL-6265", 
                              "MI-10109", "MI-10377", 
                              "MI-12507", "MI-6586",
                              "MI-6623", "MI-7119",
                              "MI-9059"),
             StateCounty = c('illinois,cook',
                             'ohio,erie', 'ohio,erie',
                             'illinois,cook', 'illinois,cook', 'illinois,cook',
                             'michigan,st clair', 'michigan,emmet',
                             'michigan,wayne', 'michigan,keweenaw',
                             'michigan,wayne', 'michigan,keweenaw',
                             'michigan,keweenaw'))

# Add counties to CaseStudyData dataframe
CaseStudyData %<>% 
  left_join(GoogleMapGroundTruth)

# Replace state.county values of NA with StateCounty actual values 
CaseStudyData %<>% 
  mutate(state.county = ifelse(is.na(state.county), StateCounty, state.county))

# Remove StateCounty
CaseStudyData$StateCounty <- NULL

# Redo the state.county.ind and county.ind 
# These are used when creating dataNew
CaseStudyData %<>% 
  mutate(state.county.ind = state.county %>% as.factor %>% as.numeric,
         county.ind = state.county %>% as.factor %>% as.numeric)

# Put an error if there are any NAs in state.county.ind
SCNA <- is.na(CaseStudyData$state.county.ind) %>% sum == 0

if(SCNA != TRUE){
  CaseStudyData %>% 
    filter(is.na(state.county)) %>% 
    select(GULocationID, Program, state.county.ind, GULatitude, GULongitude) %>% 
    unique
  stop('Some state.county values are NA')
}

# Filter out survey and species case study data
CaseStudySurveys <- CaseStudyData %>% 
  select(Program, GUEventID, GULocationID, EventType, 
         GULatitude, GULongitude, Country, State,
         Day, Month, Mar2FebWeeks, Year, Duration, Temp, Wind,
         state.county, state.county.ind, county.ind, site.ind) %>% 
  unique

CaseStudyDataSpp <- CaseStudyData %>% 
  filter(Code %in% CaseStudySpp$Code)

# Look at state breakdown by program
CaseStudySurveys %>% 
  group_by(Program, State) %>% 
  summarize(N = n())

# Check for same event or location IDs within dataset ####
# Location IDs
LocIDUnique <- CaseStudySurveys %>% 
  select(Program, GULocationID) %>% 
  unique %>%
  ungroup %>% 
  group_by(GULocationID) %>% 
  filter(n()>1) 
# Event IDs
EventIDUnique <- CaseStudySurveys %>% 
  select(Program, GUEventID) %>% 
  unique %>%
  ungroup %>% 
  group_by(GUEventID) %>% 
  filter(n()>1) 
if(nrow(LocIDUnique) > 0 |
   nrow(EventIDUnique) > 0){
  stop('Duplicated location or event IDs')
}

# Some surveys have multiple recordings of a species ####
Spp <- unique(CaseStudyDataSpp$Code)
DF <- tibble(Program = NA, GUEventID = NA, ScientificName = NA,
             Code = NA, Count = NA)

for(ss in 1:length(Spp)){
  Events <- CaseStudyDataSpp %>% 
    filter(Code == Spp[ss]) %>% 
    group_by(GUEventID) %>% 
    filter(n()>1) %>% 
    select(Program, GUEventID, ScientificName, Code, Count)
  DF %<>% bind_rows(Events)
}
# Remove first row
DF %<>% filter(!is.na(Program))

# Number of affected surveys/counts
Affected <- DF %>% 
  group_by(Code) %>% 
  summarize(NumSurveys = n_distinct(GUEventID),
            NumRows = n(),
            NumCount = sum(Count),
            MeanCount = mean(Count))
Affected

# How many rows should we end up with after summing these
ExpectedRowsByCode <- nrow(CaseStudyDataSpp) - 
  sum(Affected$NumRows) + 
  sum(Affected$NumSurveys)

# Sum up these multiple entries by EventID 
NRowByCode <- CaseStudyDataSpp %>% 
  # Currently keeps subspecies separate (two CELLAD)
  group_by(Program, GUEventID, Code) %>% 
  summarize(Count = sum(Count)) %>% 
  nrow

# Expected and actual rows should line up (TRUE)
OK <- ExpectedRowsByCode == NRowByCode
if(OK != TRUE){
  stop("Expected rows don't line up with actual rows")
}

# Sum up multiple entries by SummaryLevel (set earlier)
if(SummaryLevel == 'Species'){
  CaseStudyDataSpp %<>%
    group_by(Program, GUEventID, GULocationID, 
             EventType, GULatitude, GULongitude,
             state.county, state.county.ind, county.ind, site.ind,
             Country, State, Day, Month, Year, 
             Code,
             Duration, Temp, Wind, OrdinalDay, Mar2FebWeeks) %>% 
    summarize(Count = sum(Count)) %>% 
    ungroup
} else {
  CaseStudyDataSpp %<>% 
    mutate(OriginalCode = Code,
           Code = ScientificName) %>% 
    group_by(Program, GUEventID, GULocationID, 
             EventType, GULatitude, GULongitude,
             state.county, state.county.ind, county.ind, site.ind,
             Country, State, Day, Month, Year, 
             Code, OriginalCode,
             Duration, Temp, Wind, OrdinalDay, Mar2FebWeeks) %>% 
    summarize(Count = sum(Count)) %>% 
    ungroup
}

# Fill in Zeros for species that weren't observed ####
CaseStudyDataSpp %>% 
  select(GUEventID, Code) %>%
  table() %>% 
  as.data.frame %>% 
  filter(Freq != 1) %>% 
  select(Freq) %>% 
  table
# Lots of zeros, none larger than that
# Throw an error if there are more than one row of species in an event
Problem <- CaseStudyDataSpp %>% 
  select(GUEventID, Code) %>%
  table() %>% 
  as.data.frame %>% 
  filter(Freq > 1) %>% 
  dim
if(Problem[1] > 0){
  stop('Multiple rows of a species during a survey')
}

# Grab the unique surveys
# Note that a survey without any of the ten species won't be here
Uniques <- CaseStudySurveys$GUEventID %>% unique
# Grab the unique species (not CELLAD twice)
EachSpp <- CaseStudySpp$Code %>% unique
# Create a dataframe to hold the survey/spp combos that aren't recorded
AddTo <- data.frame(Code = NA,
                    GUEventID = NA,
                    Count = NA)
# Takes a couple of minutes. 
# Loop through the unique surveys
for(uu in 1:length(Uniques)){
  print(uu)
  # Pull out just one survey
  US <- Uniques[uu]
  # Pull out the observations for that survey
  CSUS <- CaseStudyDataSpp %>% filter(GUEventID == US)
  # Loop through the species
  for(ss in 1:length(EachSpp)){
    # Check if the species is already present
    if(!(EachSpp[ss] %in% CSUS$Code)){
      # If the species isn't present, add the spp name, survey, and count of 0 
      # to the holder data frame
      AddTo %<>% add_row(Code = EachSpp[ss],
                         GUEventID = Uniques[uu],
                         Count = 0)
    }
  }
}

# Filter out the first row (NA placeholder) and rename column
# Original column name didn't seem to work in the loop, but it may
# have been due to something else. This was an easy solution anyways though
AddTo %<>% 
  filter(!is.na(GUEventID))

# Add the AddTo data to the original data frame and check if
# all species/surveys have data. All should be 1
CaseStudyDataSpp %>% 
  bind_rows(AddTo) %>% 
  select(GUEventID, Code) %>% table %>% table

# Add rows and save
# Add survey information 
AddTo %<>% left_join(CaseStudySurveys)

# Join observations with zeros ####
CaseStudyDataSpp %<>% bind_rows(AddTo)

# Should all have ones
Solution <- CaseStudyDataSpp %>% 
  select(GUEventID, Code) %>% table %>% table %>% 
  as.data.frame()
if(nrow(Solution) > 1 |
   Solution$Freq[1] != nrow(CaseStudySurveys) * 10){
  stop("Adding zero's didn't work right")
}


# Reformat data for use in different software -----------------------------
ColumnCompare <- tibble(
  new = c('usiteID', 'program', 'state.county', 
            'lat', 'long', 'yr', 'wk', 
            'monarch', 'duration',
            'site.ind', 'county.ind', 'perc.open'),
  old = c('GULocationID', 'Program', 'state.county.ind',
           'GULatitude', 'GULongitude', 'Year', 'Mar2FebWeeks',
           'Count', 'Not included yet', 
           'site.ind', 'county.ind', 'Not included yet'))
ColumnCompare

dataNew <- CaseStudyDataSpp %>% 
  select(GULocationID, Program, state.county.ind, 
         GULatitude, GULongitude,
         Year, Mar2FebWeeks, Count,
         site.ind, county.ind, Code, GUEventID) %>% 
  rename(usiteID = GULocationID, 
         program = Program,
         state.county = state.county.ind,
         lat = GULatitude,
         long = GULongitude,
         yr = Year,
         wk = Mar2FebWeeks,
         count = Count,
         species = Code)


# Pull out the observations for just one species if desired ---------------
# OneSppObs <- dataNew %>% 
#   filter(species == OneSpp)

# Final data prep for fitting the models ----------------------------------
# #Add indicators for state monitoring programs (using NABA as a reference level)
dataNew$ia.ind <- ifelse(dataNew$program=='Iowa Butterfly Survey Network',1,0)
dataNew$il.ind <- ifelse(dataNew$program=='Illinois Butterfly Monitoring Network',1,0)
dataNew$mi.ind <- ifelse(dataNew$program=='Michigan Butterfly Network',1,0)
dataNew$oh.ind <- ifelse(dataNew$program=='OhioLeps',1,0)

#------------------------------------------------------------#
# Combine ecological covariates into long form (nrows = n_counties * n_years * n_weeks)
#------------------------------------------------------------# 
uyears <- sort(unique(dataNew$yr))
ucounties <- sort(unique(dataNew$county.ind))
n_counties <- length(ucounties)
uweeks <- sort(unique(dataNew$wk))
uspecies <- sort(unique(dataNew$species))

cyw <- expand.grid(wk=uweeks,county.ind=ucounties,yr=uyears,species=uspecies)
cyw$yr.ind <- cyw$yr - min(cyw$yr) + 1
#Standarize week
wk.m <- mean(cyw$wk)
wk.sd <- sd(cyw$wk)
cyw$wk.st <- (cyw$wk-wk.m)/wk.sd
# Create an index for unique combinations of county, year, and week
cyw$cyw.ind <- 1:nrow(cyw)

# Join to the observation dataframe
dataNew %<>% left_join(cyw) 

# Change county indicator to start at one
CountyInd <- tibble(county.ind = cyw$county.ind %>% unique,
                    co.ind = 1:length(unique(cyw$county.ind)))
cyw %<>% left_join(CountyInd)

#Survey-level covariates in matrix
X_survey <- as.matrix(dataNew[,c('ia.ind','il.ind','mi.ind','oh.ind')])

# Format all of the Survey level data
SurveyList <- dataNew %>% 
  mutate(years = yr %>% factor,
         counties = county.ind %>% factor,
         surveys_s = paste(GUEventID, species, sep = '_'),
         cyws.ind = paste(cyw.ind, species) %>% as.factor,
	 wk = wk,
         cyws_id = as.numeric(cyws.ind)) %>% 
  rename(sites = usiteID,
         programs = program,
         y_count = count,
         surveys = GUEventID,
	 wk = wk,
         cyw_id = cyw.ind) %>% 
  select(years, programs, sites, counties, surveys, surveys_s, cyws_id,
         species, wk, y_count) %>% 
  compose_data()

# Format all site level data
SiteList <- cyw %>% 
  mutate(cyw = cyw.ind %>% factor,
         cyws = paste(cyw, species, sep = '_'),
         year_st = (yr.ind - mean(yr.ind))/sd(yr.ind),) %>% 
  rename(year_id = yr.ind,
         county_id = co.ind,
         week_st = wk.st,
         species_id = species) %>% 
  select(cyw, cyws, year_id, year_st, county_id, week_st, species_id) %>% 
  compose_data()
# SiteList %>% str

# Other data needed that aren't in one of the others
OtherList <- list(X_survey = X_survey,
                  n_cov_beta=ncol(X_survey))

# Combine them all
long.data <- c(SurveyList,
               SiteList,
               OtherList)
# The long.data object is formatted for potential use in Stan. 

# Format data for spAbundance ---------------------------------------------
# Species names
sp.names <- dataNew$species %>% unique %>% sort()
# Number of years
n.years <- long.data$n_years
# Number counties
n.counties <- long.data$n_counties
# Number of sites
n.sites <- long.data$n_sites
# Number of species
n.species <- long.data$n_species
# Number of weeks
n.weeks <- n_distinct(long.data$week_st)
# Total size of potential data points
n.total <- n.years * n.counties * n.species * n.weeks
# Get linking indexes
sp.indx.obs <- long.data$species
county.indx.obs <- long.data$counties
year.indx.obs <- long.data$years
site.indx.obs <- long.data$sites
# Convert week to a factor then back to numeric to get it as a linking index
week.indx.obs <- as.numeric(factor(long.data$wk))
# Total number of observations
n.obs.total <- length(long.data$y_count)
# Get all observational level data in a data frame
obs.data.df <- data.frame(years = year.indx.obs,
			  sites = site.indx.obs,
			  surveys = long.data$surveys_s,
			  counties = long.data$counties,
			  species = sp.indx.obs,
			  wk = week.indx.obs,
                          counts = long.data$y_count,
                          long.data$X_survey, 
                          lat = dataNew$lat, 
                          long = dataNew$long, 
                          program = dataNew$program)
# Total number of butterflies observed at a given site, in a given year, during
# a given week. Some of the state-wide surveys have multiple observations per week, 
# which is accounted for by using the number of surveys in a given week as a variable
# in the model. 
count.df <- obs.data.df %>%
  group_by(years, sites, counties, species, wk) %>%
  summarize(sum.counts = sum(counts), 
	    n.surveys = n_distinct(surveys), 
	    ia.ind = unique(ia.ind), 
	    il.ind = unique(il.ind), 
	    mi.ind = unique(mi.ind),
	    oh.ind = unique(oh.ind)) %>%
  ungroup()

# Create arrays of data for use in spAbundance
y <- array(NA, dim = c(n.species, n.sites, n.years, n.weeks))
p.abund <- long.data$n_cov_beta
program.covs <- array(NA, dim = c(n.sites, n.years, n.weeks, p.abund))
year.cov <- array(NA, dim = c(n.sites, n.years, n.weeks))
week.cov <- array(NA, dim = c(n.sites, n.years, n.weeks))
site.cov <- array(NA, dim = c(n.sites, n.years, n.weeks))
county.cov <- array(NA, dim = c(n.sites, n.years, n.weeks))
survey.cov <- array(NA, dim = c(n.sites, n.years, n.weeks))
for (i in 1:nrow(count.df)) {
  y[count.df$species[i], count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$sum.counts[i]
  if (count.df$species[i] == 1) {
    program.covs[count.df$sites[i], count.df$years[i], count.df$wk[i], 1] <- count.df$ia.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wk[i], 2] <- count.df$il.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wk[i], 3] <- count.df$mi.ind[i]
    program.covs[count.df$sites[i], count.df$years[i], count.df$wk[i], 4] <- count.df$oh.ind[i]
    year.cov[count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$years[i]
    week.cov[count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$wk[i]
    site.cov[count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$sites[i]
    county.cov[count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$counties[i]
    survey.cov[count.df$sites[i], count.df$years[i], count.df$wk[i]] <- count.df$n.surveys[i]
  }
  # y[sp.indx.obs[i], site.indx.obs[i], year.indx.obs[i], week.indx.obs[i]]
  # long.data$y_count[i]
}

# Collapse the data into a species x site x "replicate" three-dimensional array. 
# Order of the replicates is week, then year within week. 
y <- array(y, dim = c(n.species, n.sites, n.years * n.weeks))
ia.ind <- array(program.covs[, , , 1], dim = c(n.sites, n.years * n.weeks))
il.ind <- array(program.covs[, , , 2], dim = c(n.sites, n.years * n.weeks))
mi.ind <- array(program.covs[, , , 3], dim = c(n.sites, n.years * n.weeks))
oh.ind <- array(program.covs[, , , 4], dim = c(n.sites, n.years * n.weeks))
year.cov <- array(year.cov, dim = c(n.sites, n.years * n.weeks))
week.cov <- array(week.cov, dim = c(n.sites, n.years * n.weeks))
site.cov <- array(site.cov, dim = c(n.sites, n.years * n.weeks))
county.cov <- array(county.cov, dim = c(n.sites, n.years * n.weeks))
survey.cov <- array(survey.cov, dim = c(n.sites, n.years * n.weeks))

# Remove missing values to make size of stored objects smaller. In spAbundance, 
# this works for non-spatial models to reduce the size of model objects. 
y <- matrix(y, n.species, n.sites * n.years * n.weeks)
y.na <- is.na(y[1, ])
y <- y[, !y.na] 
ia.ind <- c(ia.ind)[!y.na]
il.ind <- c(il.ind)[!y.na]
mi.ind <- c(mi.ind)[!y.na]
oh.ind <- c(oh.ind)[!y.na]
year.cov <- c(year.cov)[!y.na]
week.cov <- c(week.cov)[!y.na]
site.cov <- c(site.cov)[!y.na]
county.cov <- c(county.cov)[!y.na]
# Subtract one from the survey variable such that the intercept is equal to 
# abundance for a single survey. 
survey.cov <- c(survey.cov)[!y.na] - 1

# Get data in format for msAbund
abund.covs <- list(ia.ind = ia.ind, 
		   il.ind = il.ind, 
		   mi.ind = mi.ind, 
		   oh.ind = oh.ind,
		   year.cov = year.cov,
		   week.cov = week.cov,
		   site.cov = site.cov,
		   county.cov = county.cov, 
                   survey.cov = survey.cov)
# Add species names
rownames(y) <- sp.names
data.list <- list(y = y, 
		  covs = abund.covs)
save(data.list, file = 'data/spAbundance-data.rda')

# Generate plot of study locations ----------------------------------------
counts.by.lat.long <- CaseStudyData %>%
  filter(Code %in% CaseStudySpp$Code) %>%
  group_by(Program, GULocationID, GULatitude, GULongitude) %>%
  summarize(Count = sum(Count)) %>%
  ungroup() %>%
  mutate(Program = ifelse(Program == 'NFJ', 
                          'NABA', 
                          ifelse(Program == 'OhioLeps',
                                 'OH',
                                 ifelse(Program == 'Illinois Butterfly Monitoring Network', 
			  'IL', ifelse(Program == 'Iowa Butterfly Survey Network', 
			  'IA', ifelse(Program == 'Michigan Butterfly Network', 'MI', NA)))))) %>%
  mutate(Program = factor(Program, levels = c('IL', 'IA', 'MI', 'OH', 'NABA'), ordered = TRUE))

coords.sf <- st_as_sf(counts.by.lat.long, 
		      coords = c('GULongitude', 'GULatitude'), 
                      crs = 4326)
coords.sf <- coords.sf %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa <- usa %>%
  st_transform(crs = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
# Full data
mw.states <- usa %>%
  filter(ID %in% c('iowa', 'wisconsin', 'michigan', 'indiana', 
		   'illinois', 'ohio'))
my.cols <- c('#440154FF', pal_tron()(5))
my.cols <- my.cols[1:5]
names(my.cols) <- c('NABA', 'IA', 'IL', 'MI', 'OH')
ggplot() + 
  geom_sf(data = mw.states, col = 'black', fill = 'white', alpha = 1) + 
  geom_sf(data = coords.sf, aes(col = Program, size = Count)) + 
  theme_bw(base_size = 18) + 
  scale_color_manual(values = my.cols) +
  guides(color = guide_legend(override.aes = list(size = 2), 
			      order = 1)) +
  scale_size(breaks = c(1, 1000, 5000, 10000), labels = c('1', '1,000', '5,000', '10,000')) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))
# Save plot to hard drive.
ggsave(device = 'png', filename = 'figures/mw-data.png', height = 6, width = 8,
        units = 'in')
