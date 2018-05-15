# Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Matt Orton (https://github.com/m-orton/R-Scripts) for contributions to the latitude trait section (lines 53-63).
# Adapted lines 177-178 and 235-238 from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/403310/matt-dowle.
# https://stackoverflow.com/questions/13273833/merging-multiple-data-table.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

##################################################################################################################

##### SECTION 2: TRAIT ASSIGNMENT #####
# This section is designed to assign trait data from BOLD and FishBase and match it to BOLD sequence data. The version here is a 
# shortened version and only considers a subset of traits. Many more traits are available on FishBase and can be easily accessed using 
# the package rfishbase.

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

##################################################################################################################

### TRAIT: MEDIAN LATITUDE ###
# Currently, median latitude is the only trait whose information is taken from BOLD. The rest of the data will be obtained from FishBase.
# Filtering for presence of a latitude value.
dfLatitudeSpecies <- dfFiltered[grep("[0-9]", lat)]
# Convert the latitude (lat) column to number instead of character type
dfLatitudeSpecies[, lat_num := as.numeric(lat)]
# Conversion to absolute values before median latitude values are calculated.
dfLatitudeSpecies[, abs_lat_num := abs(lat_num)]
# Determine a median latitude for each BIN using absolute values.
dfLatitudeSpecies[, median_lat := median(abs_lat_num), keyby = bin_uri]

# While considering traits for eventual multivariate analyses, it is necessary for them to have an adequate sample size 
# (i.e. over x # rows of data, depending on your purposes). In addition, they should exhibit some amount of variation across the observations.

# Use the GetTraitSpecificDataBIN function to obtain a subset of data for those species that have latitude data available.
# Get the trait specific datatable.
dfLatitude <- setDT(GetTraitSpecificDataBIN(dfLatitudeSpecies, 16))
# Get information for the trait.
GetTraitInfo(dfLatitude, "median_lat", type = "continuous")

# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, filtered_bin_size, recordID, order_name = order_label, family_name = family_label, genus_name = genus_label,
                             species_name = species_label, nucleotides)]

##################################################################################################################

### FISHBASE TRAITS ###
# In this section, traits from FishBase are extracted using the rfishbase package and matched against the information obtained from BOLD.
# Extract all of the species names that are available on FishBase.
allFish <- fishbase
# Paste Genus and Species columns together to get the species name.
allFish$fullName <- paste(allFish$Genus, allFish$Species)
fishBaseSpecies <- allFish$fullName  # 33104 species names.
# Match the species labels from BOLD with the species names from FishBase.
# Make fishBaseSpecies into a dataframe first so it can be merged with dfFiltered.
dfFishBaseSpecies <- data.frame(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
# These are the names that match.
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by = "species_name")
# Extract species' name as a vector if trying to access trait information for first time (aka if you haven't saved trait info in your current working
# directory (CWD) yet).
speciesNames <- unique(dfBoldBase$species_name)
speciesNames <- as.character(speciesNames)

### TRAIT ASSIGNMENT AND RECODING SECTION ###
# As there are multiple entries per species for some traits, I want to take the median or mode value of some traits. This will depend on the type of trait 
# i.e. continuous vs. categorical).
# Note: I am just using a subset of traits from FishBase here as example traits.

### SPECIES TRAITS ###
dfSpecies <- data.frame(rfishbase::species(speciesNames))  # Could take a while.
# Storing this as a file in CWD.
write.csv(dfSpecies, file = "species_info.csv")
# Read it back in as a datatable.
dfSpecies <- fread("species_info.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable syntax consistent throughout the pipeline.
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, body_shape = BodyShapeI, male_max_length = Length, male_type_length = LTypeMaxM, 
                                 female_max_length = LengthFemale, female_type_length = LTypeMaxF, freshwater = Fresh, saltwater = Saltwater, 
                                 brackish = Brack)]

# TRAIT: Body Shape.
# Get the trait specific datatable.
dfBodyShape <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
# Get the trait information.
GetTraitInfo(dfBodyShape, "body_shape", type = "discrete")
# To remove rare categories, uncomment the following lines (i.e. "other" could be a rare category that you want to remove).
#rareVars <- which(dfBodyShape$body_shape == "other")
#dfBodyShape <- dfBodyShape[-rareVars, ]
# Make it a factor variable.
dfBodyShape[, body_shape := as.factor(body_shape)]
# If you removed a trait, make sure to drop the level from the factor with the following line.
#dfBodyShape$body_shape <- droplevels(dfBodyShape$body_shape)

# TRAIT: Male/unsexed maximum length.
# We only want total length measurements for this trait.
dfMaleMaxLength <- dfSpeciesTraits[male_type_length == "TL"]
dfMaleMaxLength <- setDT(GetTraitSpecificData(dfMaleMaxLength, 3))
GetTraitInfo(dfMaleMaxLength, "male_max_length", type = "continuous")

# TRAIT: Female maximum length.
# We only want total length measurements for this trait.
dfFemaleMaxLength <- dfSpeciesTraits[female_type_length == "TL"]
dfFemaleMaxLength <- setDT(GetTraitSpecificData(dfFemaleMaxLength, 5))
GetTraitInfo(dfFemaleMaxLength, "female_max_length", type = "continuous")

# Let's take the female max length when it is available, otherwise take male/unsexed max length.
dfMaxLength <- merge(dfMaleMaxLength, dfFemaleMaxLength, all = T, by = "species_name")
dfMaxLength$original_fem <- dfMaxLength$female_max_length
# Replace values that are missing in female_max_length with those from male_max_length.
dfMaxLength$female_max_length[is.na(dfMaxLength$female_max_length)] <- dfMaxLength$male_max_length[is.na(dfMaxLength$female_max_length)]
setnames(dfMaxLength, "female_max_length", "max_length")
# Get the trait specific datatable.
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 3))
# Get the trait information.
GetTraitInfo(dfMaxLength, "max_length", type = "continuous")

# TRAIT: Salinity.
# I have to do some recoding for this trait. I only want purely saltwater, freshwater and brackish species. This will simplify comparison.
# Subset the salinity data.
dfSalinity <- dfSpeciesTraits[, .(species_name, saltwater, freshwater, brackish)]
# Which columns are integers?
integerVars <- dfSalinity[, sapply(.SD, is.integer)]
integerVars <- which(integerVars == "TRUE")
# Recode the variables (from "-1" to "1"...more intuitive).
dfSalinity[, (integerVars) := lapply(.SD, function(x) as.numeric(revalue(as.character(x), c("-1" = "1")))), .SDcols = integerVars]
# I want to only look at purely saltwater, brackish, or freshwater species. So, remove species that inhabit more than 1 water type.
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

# Finally, prepare the dfSpeciesGenMV datatable by merging all univariate 
# datatables. This datatable will be used for the eventual multivariate analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(..., all = T), list(dfBodyShape, dfMaxLength, dfSalinity))

### ECOLOGY TRAITS ###
dfEcology <- data.frame(ecology(speciesNames))
write.csv(dfEcology, file = "ecology_info.csv") 
dfEcology <- fread("ecology_info.csv")
# Get rid of columns I do not need for the regression analysis. We are
# just doing to experiment with a few traits right now.
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, .(species_name = sciname, lakes = Lakes, oceanic = Oceanic, benthic = Benthic, diet_troph = DietTroph)]
# As some of the variables are coded as integers right now, we need to recode them to the types needed for regression analysis (i.e. factor type).
# Which columns are integers?
integerVars <- dfEcologyTraits[, lapply(.SD, is.integer)]
integerVars <- which(integerVars == "TRUE")
# Which columns are characters?
characterVars <- dfEcologyTraits[, lapply(.SD, is.character)]
characterVars <- which(characterVars == "TRUE")
# Except for species_name.
characterVars <- tail(characterVars, -1)
# Combine integerVars and characterVars.
changeVars <- c(integerVars, characterVars)
# Change all of the character and integer columns to factor type.
dfEcologyTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]
# Recode the variables (from "-1" to "1").
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), .SDcols = integerVars]

### Categorical traits. ###
# TRAIT: Lakes (binary).
# Get the trait specific datatable.
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
# Get trait information.
GetTraitInfo(dfLakes, "lakes", type = "discrete")

# TRAIT: Oceanic (binary).
# Get the trait specific datatable.
dfOceanic <- setDT(GetTraitSpecificData(dfEcologyTraits, 3))
# Get trait information.
GetTraitInfo(dfOceanic, "oceanic", type = "discrete")

# TRAIT: Benthic (binary).
# Get the trait specific datatable.
dfBenthic <- setDT(GetTraitSpecificData(dfEcologyTraits, 4))
# Get trait information.
GetTraitInfo(dfBenthic, "benthic", type = "discrete")

# TRAIT: DietTroph (continuous).
# Get the trait specific datatable.
dfDietTroph <- setDT(GetTraitSpecificData(dfEcologyTraits, 5))
# Get trait information.
GetTraitInfo(dfDietTroph, "diet_troph", type = "continuous")

# Finally, prepare the dfEcologyMV datatable by merging all univariate datatables.
# Merging multiple datatables at once.
dfEcologyMV <- Reduce(function(...) merge(..., all = T), list(dfLakes, dfOceanic, dfBenthic, dfDietTroph))

# Reproduction.
dfReproduction <- data.frame(reproduction(speciesNames))
write.csv(dfReproduction, file = "reproduction_info.csv") 
# Read in the reproduction information.
dfReproduction <- fread("reproduction_info.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable syntax consistent throughout the pipeline.
dfReproTraits <- dfReproduction[, .(species_name = sciname, repro_mode = ReproMode, parental_care = ParentalCare)]
# Unlikely to vary within species.
dfReproTraits <- dfReproTraits[!duplicated(dfReproTraits$species_name), ]
# First, change all of the traits to factor type.
dfReproTraits[, 2:3 := lapply(.SD, as.factor), .SDcols = 2:3]

# TRAIT: Repro mode.
# Get the trait specific datatable.
dfReproMode <- setDT(GetTraitSpecificData(dfReproTraits, 2))
# Get trait information.
GetTraitInfo(dfReproMode, "repro_mode", type = "discrete")
# Remove rare variables. Example:
#dfReproMode <- dfReproMode[!(repro_mode == "true hermaphroditism")]
# Drop levels.
#dfReproMode$repro_mode <- droplevels(dfReproMode$repro_mode)

# TRAIT: Parental care.
# Get the trait specific datatable.
dfParCare <- setDT(GetTraitSpecificData(dfReproTraits, 3))
# Get trait information.
GetTraitInfo(dfParCare, "parental_care", type = "discrete")

# Merging the datatables.
dfReproMV <- merge(dfReproMode, dfParCare, by = "species_name")

# Construction of dfTraits datatable.
# This table contains all the potential traits for multivariate analysis.
# Let's first merge the trait information back to dfFiltered.
# NA/NULL/blank for those species that don't have info for that particular trait.
# Note: I only want a single row per species from dfFiltered for this merging process (i.e. I just want the bin name, species name, and size of the bin). 
dfFilteredSingle <- dfFiltered[!duplicated(dfFiltered$species_name), ][, .(bin_uri, species_name, filtered_bin_size)]
# Have to merge back with latitude using bin_uri and we obtained this data from the BOLD records.
dfTraits <- merge(dfFilteredSingle, dfLatitude, all = T, by = "bin_uri")
# Set the keys for datatable merging.
dfTraits <- Reduce(function(...) merge(..., all = T, by = "species_name"), list(dfTraits, dfSpeciesGenMV, dfEcologyMV, dfReproMV))
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:3, 6:15)]
# Remove species that don't have ANY trait information available.
dfMissing <- dfTraits[apply(dfTraits[, c(4:13)], 1, function(x) all(is.na(x))), ]
dfTraits <- dfTraits[!dfTraits$bin_uri %in% dfMissing$bin_uri, ]
# Merge back to dfFiltered to obtain all of the sequence information for each BIN. This is for creation of the master phylogeny.
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")[, 1:8]
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[7] <- "species_name"
