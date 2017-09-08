# Copyright (C) 2017 Jacqueline May.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariable analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.
# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.

# Contributions & Acknowledgements #
# Matt Orton (https://github.com/m-orton/R-Scripts) for contributions to the 
# latitude trait section (lines 52-63).
# Adapted lines 199-200 and 319-323 from code shared in Stack Overflow 
# discussion:
# Author: https://stackoverflow.com/users/403310/matt-dowle.
# https://stackoverflow.com/questions/13273833/merging-multiple-data-table.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# There is a copy of the GNU General Public License
# along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

################################################################################

##### SECTION 2: TRAIT ASSIGNMENT #####
# This section is designed to assign trait data from BOLD and FishBase and match
# it to BOLD sequence data. The version here is a shortened version and only
# considers a subset of traits. Many more traits are available on FishBase
# and can be easily accessed using the package rfishbase.

### PACKAGES REQUIRED ###
# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages("plyr")
library(plyr)
# For FishBase data:
#install.packages("rfishbase")
library(rfishbase)
# Load the function(s) designed for this script:
source("GetTraitSpecificDataBIN.R")
source("GetTraitSpecificData.R")

################################################################################

### TRAIT: MEDIAN LATITUDE ###
# Currently, median latitude is the only trait whose information is taken from 
# BOLD. The rest of the data will be obtained from FishBase.
# Filtering for presence of a latitude value.
dfLatitudeSpecies <- dfFiltered[grep("[0-9]", lat)]
# Convert the latitude (lat) column to number instead of character type
dfLatitudeSpecies[, lat_num := as.numeric(lat)]
# Conversion to absolute values before median latitude values are calculated.
dfLatitudeSpecies[, abs_lat_num := abs(lat_num)]
# Determine a median latitude for each BIN using absolute values.
dfLatitudeSpecies[, median_lat := median(abs_lat_num), keyby = bin_uri]

# While considering traits for eventual multivariable analyses, it is necessary
# for them to have an adequate sample size (i.e. over x # rows of data, 
# depending on your purposes). In addition, they should exhibit some amount of 
# variation across the observations. For example, categorical traits shouldn't 
# have all of their data in one category.

# TRAIT: Latitude.
# Use the GetTraitSpecificDataBIN function to obtain latitude data on a per 
# BIN basis.
# TEST 1: Does trait have data for at least 300 species?
# Note: I am also converting to datatable format in process.
dfLatitude <- setDT(GetTraitSpecificDataBIN(dfLatitudeSpecies, 17))
nrow(dfLatitude)
# TEST 2: Does the trait have enough data variation?
hist(dfLatitude$median_lat)
range(dfLatitude$median_lat)

# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name = order_label, 
                             family_name = family_label, 
                             genus_name = genus_label, 
                             species_name = species_label, nucleotides, 
                             initial_bin_size, filtered_bin_size)]

################################################################################

### FISHBASE TRAITS ###
# In this section, traits from FishBase are extracted using the rfishbase 
# package and matched against the information obtained from BOLD.
# Extract all of the species names that are available on FishBase.
allFish <- fishbase
# Paste Genus and Species columns together to get the species name.
allFish$fullName <- paste(allFish$Genus, allFish$Species)
fishBaseSpecies <- allFish$fullName  # 33104 species names.
# Match the species labels from BOLD with the species names from FishBase.
# Make fishBaseSpecies into a dataframe first so it can be merged with 
# dfFiltered.
dfFishBaseSpecies <- data.table(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
# These are the names that match.
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by = "species_name")
# Extract species' name as a vector if trying to access trait information
# for first time (aka if you haven't saved trait info in your current working
# directory (CWD) yet).
speciesNames <- dfBoldBase[, unique(species_name)]

### TRAIT ASSIGNMENT AND RECODING SECTION ###
# As there are multiple entries per species for some traits, I want to take the 
# median or mode value of some traits. This will depend on the type of trait 
# i.e. continuous vs. categorical).
# Note: Only some traits are shown here as examples.

### SPECIES TRAITS ###
dfSpecies <- data.frame(species(speciesNames))  # Could take a while.
# Storing this as a file in CWD.
write.csv(dfSpecies, file = "species_data.csv")
# Read it back in as a datatable.
dfSpecies <- fread("species_data.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline. Here, I am only extracting the 
# trait data I am interested in from this datatable: longevity, body length,
# and salinity.
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, 
                                 longevity = LongevityWild, 
                                 max_length = Length, 
                                 type_length = LTypeMaxM,
                                 freshwater = Fresh,
                                 saltwater = Saltwater,
                                 brackish = Brack)]

# Use the GetTraitSpecificData function to obtain trait data on a per species
# basis. This function will also remove species that have no data available for
# the trait of interest.

# TRAIT: Longevity.
# First, make sure it is a numeric variable so I can perform statistical tests 
# on the data down the line.
dfSpeciesTraits[, longevity := as.double(longevity)]
# TEST 1: Does the trait have an adequate sample size?
dfLongWild <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
nrow(dfLongWild)
# TEST 2: Does the trait have enough data variation?
hist(dfLongWild$longevity)
range(dfLongWild$longevity)

# TRAIT: Maximum length.
dfSpeciesTraits[, max_length := as.double(max_length)]
# We only want total length measurements.
dfMaxLength <- dfSpeciesTraits[which(type_length == "TL")]
# TEST 1: Does the trait have an adequate sample size?
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 3))
nrow(dfMaxLength)
# TEST 2: Does the trait have enough data variation?
hist(dfMaxLength$max_length)
range(dfMaxLength$max_length)

# TRAIT: Salinity.
# I have to do some recoding for this trait. I only want purely saltwater,
# freshwater and brackish species. This will simplify comparison.
# Subset the salinity data.
dfSalinity <- dfSpeciesTraits[, .(species_name, saltwater, freshwater, brackish)]
# Which columns are integers?
integerVars <- dfSalinity[, sapply(.SD, is.integer)]
integerVars <- which(integerVars == "TRUE")
# Recode the variables (from "-1" to "1"...more intuitive).
dfSalinity[, (integerVars) := lapply(.SD, function(x) as.numeric(revalue(as.character(x), c("-1" = "1")))), 
           .SDcols = integerVars]
# I want to only look at purely saltwater, brackish, or freshwater species. So,
# remove species that inhabit more than 1 water type.
# Count number of water types by summing the rows.
dfSalinity[, num_salinity_types := rowSums(.SD), .SDcols = integerVars]
# Only keep species with one water type.
dfSalinity <- dfSalinity[which(num_salinity_types == 1)]
# Now let's combine them into one "Salinity" variable.
dfSalinity <- as.data.frame(dfSalinity)
dfSalinity$salinity <- colnames(dfSalinity[2:4])[max.col(dfSalinity[2:4])]
# TEST 1: Does the trait have an adequate sample size?
dfSalinity <- setDT(GetTraitSpecificData(dfSalinity, 6))
nrow(dfSalinity)
# TEST 2: Does the trait have enough data variation?
# Convert to factor type.
dfSalinity[, salinity := as.factor(salinity)]
table(dfSalinity$salinity)
# To remove a rare variable (e.g. brackish), you can use the following lines.
#dfSalinity <- dfSalinity[!salinity == "brackish"]
# Make sure to also drop the level from the factor.
#dfSalinity[, salinity := droplevels(salinity)]

# Finally, prepare the dfSpeciesGenMV datatable by merging all single variable
# datatables. This datatable will be used for the eventual multivariable 
# analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(..., all = T), 
                         list(dfLongWild, dfMaxLength, dfSalinity))

### LIFE HISTORY RELATED ###
# Maturity.
dfMaturity <- data.frame(maturity(speciesNames))
write.csv(dfMaturity, file = "maturity_data.csv") 
dfMaturity <- fread("maturity_data.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfMaturityTraits <- dfMaturity[, .(species_name = sciname, 
                                   age_at_maturity = tm)]

# Median trait(s).
# TRAIT: Age at maturity.
# The column must first be converted to double (numeric) type.
dfAgeMaturity <- dfMaturityTraits[, age_at_maturity := 
                                    as.double(age_at_maturity)]
# The median value is then determined for each species as there can be multiple
# entries per species in this dataframe.
dfAgeMaturity[, age_at_maturity := median(age_at_maturity, na.rm = TRUE), 
              keyby = species_name]
# TEST 1: Does the trait have an adequate sample size?
dfAgeMaturity <- setDT(GetTraitSpecificData(dfAgeMaturity, 2))
nrow(dfAgeMaturity)
# TEST 2: Does the trait have enough data variation?
hist(dfAgeMaturity$age_at_maturity)
range(dfAgeMaturity$age_at_maturity) 

### ECOLOGY TRAITS ###
dfEcology <- data.frame(ecology(speciesNames))
write.csv(dfEcology, file = "ecology_data.csv") 
dfEcology <- fread("ecology_data.csv")
# Get rid of columns I do not need for the regression analysis. We are
# just doing to experiment with a few traits right now.
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, .(species_name = sciname, lakes = Lakes, 
                                 streams = Stream)]
# As some of the variables are coded as integers right now, we need to recode 
# them to the types needed for regression analysis (i.e. factor type).
# Which columns are integers?
integerVars <- dfEcologyTraits[, sapply(.SD, is.integer)]
integerVars <- which(integerVars == "TRUE")
# Which columns are characters?
characterVars <- dfEcologyTraits[, sapply(.SD, is.character)]
characterVars <- which(characterVars == "TRUE")
# Except for species_name.
characterVars <- tail(characterVars, -1)
# Combine integerVars and characterVars.
changeVars <- c(integerVars, characterVars)
# Change all of the character and integer columns to factor type.
dfEcologyTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]
# Recode the variables (from "-1" to "1").
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) 
  revalue(x, c("-1" = "1"))), .SDcols = integerVars]

### Categorical traits. ###
# Binary trait(s).
# TRAIT: Lakes.
# TEST 1: Does the trait have an adequate sample size?
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
nrow(dfLakes)
# TEST 2: Does the trait have enough data variation?
table(dfLakes$lakes)

# TRAIT: Streams.
# TEST 1: Does the trait have an adequate sample size?
dfStreams <- setDT(GetTraitSpecificData(dfEcologyTraits, 3))
nrow(dfStreams)
# TEST 2: Does the trait have enough data variation?
table(dfStreams$streams)

# Finally, prepare the dfEcologyMV datatable by merging all single variable
# datatables.
dfEcologyMV <- merge(dfLakes, dfStreams, all = T)

# Reproduction.
dfReproduction <- data.frame(reproduction(speciesNames))
write.csv(dfReproduction, file = "reproduction_data.csv") 
# Read in the reproduction information.
dfReproduction <- fread("reproduction_data.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfReproTraits <- dfReproduction[, .(species_name = sciname, parental_care = ParentalCare)]
# Unlikely to vary within species.
dfReproTraits <- dfReproTraits[!duplicated(species_name)]

# Multi-level trait.
# TRAIT: Parental Care.
# TEST 1: Does trait have data for at least 500 species?
dfParentalCare <- setDT(GetTraitSpecificData(dfReproTraits, 2))
# TEST 2: Does the trait have enough data variation?
table(dfParentalCare$parental_care)
# If there are rare categories, remove them! For example:
#dfParentalCare <- dfParentalCare[!parental_care == "biparental"]
# Also dropping these levels from the factor.
dfParentalCare[, parental_care := as.factor(parental_care)]
dfParentalCare[, parental_care := droplevels(parental_care)]

# Construction of dfTraits datatable.
# This table will contain all of the potential traits for multivariable 
# analysis.NA/NULL/blank for those species that don't have info for that 
# particular trait.
# Let's merge the trait information back to dfFiltered.
# Note: I only want a single row per BIN from dfFiltered for this merging 
# process (i.e. I just want information for the BIN URI, species name, and size
# of the bin). bin_uri
dfFilteredSingle <- dfFiltered[!duplicated(bin_uri)]
# Let's take the columns we need to construct the dfTraits datatable.
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_name, 
                                         initial_bin_size, filtered_bin_size)]
# Now merge to all of the single variable trait datatables.
dfTraits <- merge(dfFilteredSingle, dfLatitude, all = T, by = "bin_uri")
# Set the keys for datatable merging.
setkey(dfTraits, species_name)
setkey(dfSpeciesGenMV, species_name)
setkey(dfEcologyMV, species_name)
setkey(dfAgeMaturity, species_name)
setkey(dfParentalCare, species_name)
dfTraits <- Reduce(function(...) merge(..., all = T), list(dfTraits, 
                                                           dfSpeciesGenMV,
                                                           dfEcologyMV,
                                                           dfAgeMaturity,
                                                           dfParentalCare))
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:4, 7:14)]
colnames(dfTraits)[4] <- "filtered_bin_size"
# Remove species that don't have ANY trait information available.
dfMissing <- dfTraits[apply(dfTraits[, c(5:12)], 1, function(x) all(is.na(x)))]
dfTraits <- dfTraits[!dfTraits$bin_uri %in% dfMissing$bin_uri]

# Merge back to dfFiltered to obtain all of the sequence information for 
# each BIN. This is for creation of the master phylogeny. 
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[6] <- "species_name"
colnames(dfPreCentroid)[8] <- "initial_bin_size"
colnames(dfPreCentroid)[9] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:9)]