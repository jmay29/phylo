# Copyright (C) 2017 Jacqueline May.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariable analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.
# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.

# Contributions & Acknowledgements #
# Matt Orton (https://github.com/m-orton/R-Scripts) for contributions to the 
# latitude trait section (lines 53-63).
# Adapted lines 177-178 and 235-238 from code shared in Stack Overflow 
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
source("GetTraitInfo.R")

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
# Note: I am also converting to datatable format in process.
dfLatitude <- setDT(GetTraitSpecificDataBIN(dfLatitudeSpecies, 16))
# Get information for the trait.
GetTraitInfo(dfLatitude, "median_lat", type = "continuous")

# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, filtered_bin_size, recordID, order_name = order_label,
                             family_name = family_label, genus_name = genus_label,
                             species_name = species_label, nucleotides)]

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
dfFishBaseSpecies <- data.frame(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
# These are the names that match.
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by = "species_name")
# Extract species' name as a vector if trying to access trait information
# for first time (aka if you haven't saved trait info in your current working
# directory (CWD) yet).
speciesNames <- unique(dfBoldBase$species_name)
speciesNames <- as.character(speciesNames)

### TRAIT ASSIGNMENT AND RECODING SECTION ###
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
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, longevity = LongevityWild, 
                                 max_length = Length, type_length = LTypeMaxM,
                                 freshwater = Fresh, saltwater = Saltwater,
                                 brackish = Brack)]

# Use the GetTraitSpecificData function to obtain trait data on a per species
# basis. This function will also remove species that have no data available for
# the trait of interest.

# TRAIT: Longevity.
# Get the trait specific datatable.
dfLongWild <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
# Get the trait information.
GetTraitInfo(dfLongWild, "longevity", type = "continuous")

# TRAIT: Maximum length.
# We only want total length measurements for this trait.
dfMaxLength <- dfSpeciesTraits[type_length == "TL"]
# TEST 1: Does the trait have an adequate sample size?
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 3))
# Get the trait information.
GetTraitInfo(dfMaxLength, "max_length", type = "continuous")

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
dfSalinity <- dfSalinity[num_salinity_types == 1]
# Now let's combine them into one "Salinity" variable.
dfSalinity <- as.data.frame(dfSalinity)
dfSalinity$salinity <- colnames(dfSalinity[2:4])[max.col(dfSalinity[2:4])]
# Get trait specific datatable.
dfSalinity <- setDT(GetTraitSpecificData(dfSalinity, 6))
# Convert to factor type.
dfSalinity[, salinity := as.factor(salinity)]
# Get the trait information.
GetTraitInfo(dfSalinity, "salinity", type = "discrete")
# To remove a rare variable (e.g. brackish), you can use the following lines.
#dfSalinity <- dfSalinity[!salinity == "brackish"]
# Make sure to also drop the level from the factor.
#dfSalinity[, salinity := droplevels(salinity)]

# Finally, prepare the dfSpeciesGenMV datatable by merging all single variable
# datatables. This datatable will be used for the eventual multivariable analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(..., all = T), 
                         list(dfLongWild, dfMaxLength, dfSalinity))

### ECOLOGY TRAITS ###
dfEcology <- data.frame(ecology(speciesNames))
write.csv(dfEcology, file = "ecology_data.csv") 
dfEcology <- fread("ecology_data.csv")
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, .(species_name = sciname, lakes = Lakes)]
# Recode as factor.
dfEcologyTraits[, lakes := as.factor(lakes)]
# Recode the variables (from "-1" to "1").
dfEcologyTraits[, lakes := revalue(lakes, c("-1" = "1"))]

# Binary trait(s).
# TRAIT: Lakes.
# Get the trait specific datatable.
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
# Get trait information.
GetTraitInfo(dfLakes, "lakes", type = "discrete")

### REPRODUCTION TRAIT(S) ###
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
# Get the trait specific datatable.
dfParentalCare <- setDT(GetTraitSpecificData(dfReproTraits, 2))
# Convert to factor.
dfParentalCare[, parental_care := as.factor(parental_care)]
# Get trait information.
GetTraitInfo(dfParentalCare, "parental_care", type = "discrete")

# Construction of dfTraits datatable.
# This table will contain all of the potential traits for multivariable 
# analysis. NA/NULL/blank for those species that don't have info for that 
# particular trait.
# Let's merge the trait information back to dfFiltered.
# Note: I only want a single row per BIN from dfFiltered for this merging 
# process (i.e. I just want information for the BIN URI, species name, and size
# of the bin).
dfFilteredSingle <- dfFiltered[!duplicated(bin_uri)]
# Let's take the columns we need to construct the dfTraits datatable.
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_name, filtered_bin_size)]
# Now merge to all of the single variable trait datatables.
dfTraits <- merge(dfFilteredSingle, dfLatitude, all = T, by = "bin_uri")
# Set the keys for datatable merging.
setkey(dfTraits, species_name)
setkey(dfLakes, species_name)
setkey(dfParentalCare, species_name)
dfTraits <- Reduce(function(...) merge(..., all = T), list(dfTraits, 
                                                           dfSpeciesGenMV,
                                                           dfLakes,
                                                           dfParentalCare))
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:3, 6:11)]
setnames(dfTraits, "filtered_bin_size.x", "filtered_bin_size")
# Remove species that don't have ANY trait information available.
dfMissing <- dfTraits[apply(dfTraits[, c(4:9)], 1, function(x) all(is.na(x)))]
dfTraits <- dfTraits[!dfTraits$bin_uri %in% dfMissing$bin_uri]
# Merge back to dfFiltered to obtain all of the sequence information for 
# each BIN. This is for creation of the master phylogeny. 
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")
# Dataframe reorganization and renaming.
setnames(dfPreCentroid, c("species_name.x", "filtered_bin_size.x"), 
         c("species_name", "filtered_bin_size"))
dfPreCentroid <- dfPreCentroid[, c(1:8)]
