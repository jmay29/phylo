# Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Matt Orton (https://github.com/m-orton/R-Scripts) for contributions to the latitude trait section.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

#############################################################################################################################

##### SECTION 2: TRAIT ASSIGNMENT #####
# This section is designed to assign trait data from BOLD and FishBase and match it to BOLD sequence data. The version here 
# is a shortened version and only considers a subset of traits. Many more traits are available on FishBase and can be easily 
# accessed using the package rfishbase.

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
source("GetTraitInfo.R")

#############################################################################################################################

### TRAIT: MEDIAN LATITUDE ###
# Currently, latitude is the only trait whose information is taken from BOLD. The rest of the data will be obtained from FishBase.
# Filtering for presence of a latitude value.
dfLatitude <- dfFiltered[lat %like% "[0-9]"]
# Conversion to absolute values before median latitude values are calculated.
dfLatitude[, abs_lat := abs(lat)]
# Determine a median latitude for each BIN using absolute values.
dfLatitude[, median_lat := median(abs_lat), keyby = bin_uri]
# The GetTraitSpecificDataBIN can be used to clean up the datatable. It removes duplicate values per BIN and only keeps the BIN URIs,
# species names, and the column specified (median latitude is column 15 here).
dfLatitude <- as.data.table(GetTraitSpecificDataBIN(dfLatitude, 15))
# Dataframe reorganization and renaming.
setnames(dfLatitude, "species_label", "species_name")

# In this pipeline, we will eventually be constructing a multivariable datatable once we have also downloaded the FishBase data.
# While considering traits for eventual multivariable analyses, it is necessary for them to have a decent sample size 
# (i.e. > nrows of data, depending on your purposes). In addition, they should exhibit some amount of variation 
# across the observations (all of the observations shouldn't be in one category)!

# The GetTraitInfo function can be used to obtain some information about the trait. 
# Are there a decent number of observations for this trait? Is there variation in the trait?
GetTraitInfo(dfLatitude$median_lat)

# Datatable reorganization for dfFiltered. We are only keeping the labels we assigned and not the original taxonomic classifications.
dfFiltered <- dfFiltered[, .(bin_uri, filtered_bin_size, recordID, order_name = order_label, family_name = family_label, 
                             genus_name = genus_label, species_name = species_label, nucleotides)]

### FISHBASE TRAITS ###
# In this section, traits from FishBase are extracted using the rfishbase package and matched against the information obtained from BOLD.
# First, extract all of the species names that are available on FishBase.
dfFishBase <- as.data.table(fishbase)
# Paste Genus and Species columns together to get the full species name.
dfFishBase[, species_name := paste(Genus, Species)]
# We only want the species name column.
dfFishBase <- dfFishBase[, .(species_name)]
# Match the species labels from BOLD with the species names from FishBase. These are the species with records available on both BOLD and FishBase.
dfBoldBase <- merge(dfFiltered, dfFishBase, by = "species_name")
# Extract species' name as a vector if trying to access trait information for first time (aka if you haven't saved trait info in 
# your current working directory (CWD) yet).
speciesNames <- unique(dfBoldBase$species_name)

### TRAIT ASSIGNMENT AND RECODING SECTION ###
# Note: I am just using a subset of 10 traits in this pipeline as example traits.

### SPECIES TRAITS ###
# Downloading trait information from FishBase's Species data table using the species function. 
dfSpecies <- data.frame(rfishbase::species(speciesNames))  ## Could take a while.
# Storing this as a csv file in CWD.
write.csv(dfSpecies, file = "species_info.csv")
# Read it back in as a datatable using the fread function.
dfSpecies <- fread("species_info.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable syntax consistent throughout the pipeline.
# We are looking at the traits body shape, maximum length (male_max_length and female_max_length) and salinity traits (freshwater, saltwater
# and brackish) from this table. We also need to take the type of length that was measured (male_type_length and female_max_length).
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, body_shape = BodyShapeI, male_max_length = Length, male_type_length = LTypeMaxM, 
                                 female_max_length = LengthFemale, female_type_length = LTypeMaxF, freshwater = Fresh, 
                                 saltwater = Saltwater, brackish = Brack)]

# TRAIT: Body shape.
# Convert body_shape to factor type (the type that categorical traits should be coded for regression analysis).
dfSpeciesTraits[, body_shape := as.factor(body_shape)]
# Get information about the trait body_shape.
GetTraitInfo(dfSpeciesTraits$body_shape)
# To remove rare categories, uncomment the following line (e.g. "eel-like" could be a rare category that you want to remove).
#dfSpeciesTraits[body_shape == "eel-like", body_shape := NA][, body_shape := droplevels(body_shape)]

# TRAIT: Maximum length.
# For the length data, we have both male and female data (male_max_length and female_max_length). 
# First, let's remove observations that are not "TL" or total length under male_type_length and female_type_length, to keep things consistent.
dfSpeciesTraits[male_type_length != "TL", male_max_length := NA][female_type_length != "TL", female_max_length := NA]
# Since we are dealing with mitochondrial DNA, let's take the female data when it is available, otherwise take the male observation for the species.
dfSpeciesTraits[, max_length := ifelse(is.na(female_max_length), male_max_length, female_max_length)]
# Remove the extra length columns.
dfSpeciesTraits[, c("male_type_length","female_type_length", "male_max_length", "female_max_length") := NULL]
# Get information about the trait.
GetTraitInfo(dfSpeciesTraits$max_length)

# TRAIT: Salinity.
# We have to do some recoding for the salinity trait
# Recode the salinity variables (from "-1" to "1"...more intuitive).
salinities <- c("saltwater", "freshwater", "brackish")
dfSpeciesTraits[, (salinities) := lapply(.SD, function(x) as.numeric(revalue(as.character(x), c("-1" = "1")))), .SDcols = salinities]
# I want to only look at purely saltwater, brackish, or freshwater species. So, remove species that inhabit more than 1 water type.
# Count number of water types by summing the rows.
dfSpeciesTraits[, num_salinity_types := rowSums(.SD), .SDcols = salinities]
# For those species that inhabit only one water type, select the column name that was assigned a "1" to represent that species.
dfSpeciesTraits[num_salinity_types == 1, salinity := as.factor(colnames(.SD)[max.col(.SD)]), .SDcols = salinities]
# Remove the extra salinity columns.
dfSpeciesTraits[, c("freshwater","saltwater", "brackish", "num_salinity_types") := NULL]
# Get information about the trait.
GetTraitInfo(dfSpeciesTraits$salinity)

### ECOLOGY TRAITS ###
# Now, we are downloading data from the FishBase Ecology data table using the ecology function.
dfEcology <- data.frame(ecology(speciesNames))
write.csv(dfEcology, file = "ecology_info.csv") 
dfEcology <- fread("ecology_info.csv")
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, .(species_name = sciname, lakes = Lakes, oceanic = Oceanic, benthic = Benthic, diet_troph = DietTroph)]
# As only some of the variables are coded as integers right now, we need to recode them factor types.
integerVars <- dfEcologyTraits[, lapply(.SD, is.integer)]
integerVars <- which(integerVars == "TRUE")
# Change all of the integer columns to factor type.
dfEcologyTraits[, (integerVars) := lapply(.SD, as.factor), .SDcols = integerVars]
# Recode the variables (from "-1" to "1").
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), .SDcols = integerVars]

### Categorical traits. ###
# TRAIT: Lakes (binary).
# Get trait information.
GetTraitInfo(dfEcologyTraits$lakes)

# TRAIT: Oceanic (binary).
GetTraitInfo(dfEcologyTraits$oceanic)

# TRAIT: Benthic (binary).
GetTraitInfo(dfEcologyTraits$benthic)

# TRAIT: DietTroph (continuous).
GetTraitInfo(dfEcologyTraits$diet_troph)

### REPRODUCTION TRAITS ###
# Now, we are downloading data from the FishBase Reproduction data table using the reproduction function.
dfReproduction <- data.frame(reproduction(speciesNames))
write.csv(dfReproduction, file = "reproduction_info.csv") 
dfReproduction <- fread("reproduction_info.csv")
# Datatable reorganization and renaming.
dfReproTraits <- dfReproduction[, .(species_name = sciname, repro_mode = ReproMode, parental_care = ParentalCare)]
# Unlikely to vary within species.
dfReproTraits <- dfReproTraits[!duplicated(species_name)]
# First, change the traits to factor type.
dfReproTraits[, 2:3 := lapply(.SD, as.factor), .SDcols = 2:3]

# TRAIT: Repro mode.
GetTraitInfo(dfReproTraits$repro_mode)

# TRAIT: Parental care.
GetTraitInfo(dfReproTraits$parental_care)

# Construction of dfTraits datatable.
# This table contains all the potential traits for multivariable analysis.
# Let's first merge the trait information back to dfFiltered.
# Note: I only want a single row per species from dfFiltered for this merging process (i.e. I just want the bin name, species name, and size of the bin). 
dfFilteredSingle <- dfFiltered[!duplicated(species_name)][, .(bin_uri, species_name, filtered_bin_size)]
# Merge dfTraits with all of the FishBase data.
dfTraits <- Reduce(function(...) merge(..., all = T, by = "species_name"), list(dfFilteredSingle, dfLatitude, dfSpeciesTraits, dfEcologyTraits, dfReproTraits))
# Remove species that don't have ANY trait information available on either BOLD or FishBase (these BINs passed the filters but
# do not have any trait data).
missing <- dfTraits[, apply(.SD, 1, function(x) all(is.na(x))), .SDcols = 4:13]
missing <- which(missing == TRUE)
dfTraits <- dfTraits[!missing]
# Merge back to dfFiltered to obtain all of the sequence information for each BIN. This is for creation of the master phylogeny.
# dfPreCentroid contains only information needed for selection of the centroid sequences. dfTraits has all of the trait information
# needed for downstream PGLS analyses.
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "species_name")[, 1:8]
# Dataframe reorganization and renaming.
setnames(dfPreCentroid, old = c('bin_uri.x', 'filtered_bin_size.x'), new = c('bin_uri', 'filtered_bin_size'))

# Remove objects that are not needed for Section 3.
rm(speciesNames, salinities, integerVars, missing); rm(dfBoldBase, dfFishBase, dfFilteredSingle)
rm(dfLatitude, dfSpecies, dfSpeciesTraits, dfEcology, dfEcologyTraits, dfReproduction, dfReproTraits)
