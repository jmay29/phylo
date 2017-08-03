################################################################################

# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
#                  determination/reference sequence trimming (lines TBD).

################################################################################

##### SECTION 2: TRAIT ASSIGNMENT #####

### PACKAGES REQUIRED ###

# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages("plyr")
library(plyr)

# For FishBase data:
#install.packages("rfishbase")
library(rfishbase)

# Load the function(s) that I designed for this script:
source("GetTraitSpecificDataBIN.R")
source("GetTraitSpecificData.R")

################################################################################

### TRAIT: MEDIAN LATITUDE ###
# Currently, median latitude and latitudinal range are the only traits whose 
# information is taken from BOLD. The rest of the data will be obtained from 
# FishBase.
# Filtering for presence of a latitude value.
containLat <- dfFiltered[, grep("[0-9]", lat)]
dfLatitudeSpecies <- dfFiltered[containLat, ]
# Convert the latitude (lat) column to number instead of character type. This is
# necessary for median and range calculations.
dfLatitudeSpecies[, lat_num := as.numeric(lat)]
# Conversion to absolute values before median latitude values are calculated.
dfLatitudeSpecies[, abs_lat_num := abs(lat_num)]
# Determine a median latitude for each BIN using absolute values.
dfLatitudeSpecies[, median_lat := median(abs_lat_num), keyby = bin_uri]

# While considering traits for eventual multivariate analyses, it is necessary
# for them to have an adequate sample size (i.e. over x # rows of data, depending
# on your purposes).
# In addition, they should exhibit some amount of variation across the 
# observations.
# TRAIT: Latitude.
# Use the GetTraitSpecificDataBIN function to obtain a subset of data for those 
# species that have latitude data available.
# TEST 1: Does trait have data for at least 300 species?
dfLatitude <- setDT(GetTraitSpecificDataBIN(dfLatitudeSpecies, 17))
nrow(dfLatitude)
# TEST 2: Does the trait have enough data variation?
hist(dfLatitude$median_lat)

# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name = order_label, 
                             family_name = family_label, 
                             genus_name = genus_label, 
                             species_name = species_label, nucleotides, 
                             initial_bin_size, filtered_bin_size)]

################################################################################

### FISHBASE TRAITS ###
# In this section, traits from FishBase are extracted using the rfishbase package
# and matched against the information obtained from BOLD.
# Extract all of the species names that are available on FishBase.
allFish <- fishbase
allFish$fullName <- paste(allFish$Genus, allFish$Species)
fishBaseSpecies <- allFish$fullName  # 33104 species names.
# Match the species labels from BOLD with the species names from FishBase.
# Make into a dataframe first so they can be merged.
dfFishBaseSpecies <- data.table(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by = "species_name")
# Extract species' name as a vector if trying to access trait information
# for first time (aka haven't saved trait info in CWD yet).
speciesNames <- dfBoldBase[, unique(species_name)]

### TRAIT ASSIGNMENT AND RECODING SECTION ###
# As there are multiple entries per species for some traits, I want to take the 
# median or modevalue of some traits. This will depend on the type of trait 
# i.e. continuous vs. categorical).
# Note: Only some traits are shown here as examples.

### SPECIES TRAITS ###
dfSpecies <- data.frame(species(speciesNames))
# Storing this as a file.
write.csv(dfSpecies, file = "species_datatable.csv")
dfSpecies <- fread("species_datatable.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, 
                                 body_shape = BodyShapeI, 
                                 longevity = LongevityWild, 
                                 max_length = Length, 
                                 type_length = LTypeMaxM,
                                 max_depth = DepthRangeDeep)]
# TRAIT: Body Shape.
# Filtering for presence of body shape data. This is for the univariate analyses 
# section.
# TEST 1: Does the trait have an adequate sample size?
# Use the GetTraitSpecificData function to obtain a subset of species that
# have data available for the target trait.
dfBodyShape <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
nrow(dfBodyShape)
# TEST 2: Does the trait have enough data variation?
table(dfBodyShape$body_shape)
# "other" is a rare category.
rareVars <- which(dfBodyShape$body_shape == "eel-like")
dfBodyShape <- dfBodyShape[-rareVars, ]
# Make it a factor variable.
dfBodyShape[, body_shape := as.factor(body_shape)]
# Also dropping this level from the factor.
dfBodyShape$body_shape <- droplevels(dfBodyShape$body_shape)

# TRAIT: Longevity.
# TEST 1: Does the trait have an adequate sample size?
dfLongWild <- setDT(GetTraitSpecificData(dfSpeciesTraits, 3))
nrow(dfLongWild)
# TEST 2: Does the trait have enough data variation?
hist(dfLongWild$longevity)
range(dfLongWild$longevity)
# Make sure it is a numeric variable.
dfLongWild[, longevity := as.double(longevity)]

# TRAIT: Maximum length.
# The column must first be converted to double (numeric) type.
dfMaxLength <- dfSpeciesTraits[, max_length := as.double(max_length)]
# We only want total length measurements.
keep <- dfMaxLength[, which(type_length == "TL")]
dfMaxLength <- dfMaxLength[keep, ]
# TEST 1: Does the trait have an adequate sample size?
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 4))
nrow(dfMaxLength)
# TEST 2: Does the trait have enough data variation?
hist(dfMaxLength$max_length)
range(dfMaxLength$max_length)

# TRAIT: Maximum depth.
# The column must first be converted to double (numeric) type.
dfMaxDepth <- dfSpeciesTraits[, max_depth := as.double(max_depth)]
# TEST 1: Does the trait have an adequate sample size?
dfMaxDepth <- setDT(GetTraitSpecificData(dfMaxDepth, 6))
nrow(dfMaxDepth)
# TEST 2: Does the trait have enough data variation?
hist(dfMaxDepth$max_depth)
range(dfMaxDepth$max_depth)

# Finally, prepare the dfSpeciesGenMV datatable by merging all univariate 
# datatables. This datatable will be used for the eventual multivariate analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(..., all = T), 
                         list(dfBodyShape, dfLongWild, dfMaxLength,
                              dfMaxDepth))

### LIFE HISTORY RELATED ###
# Maturity.
dfMaturity <- data.frame(maturity(speciesNames))
write.csv(dfMaturity, file = "maturity_datatable.csv") 
# Read in the maturity information.
dfMaturity <- fread("maturity_datatable.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfMaturityTraits <- dfMaturity[, .(species_name = sciname, age_at_maturity = tm)]

# Median trait(s).
# TRAIT: Age at maturity.
# The column must first be converted to double (numeric) type.
dfAgeMaturity <- dfMaturityTraits[, age_at_maturity := as.double(age_at_maturity)]
# The median value is then determined for each species.
dfAgeMaturity[, age_at_maturity := median(age_at_maturity, na.rm = TRUE), 
              keyby = species_name]
# Filtering for presence of average depth data. This is for the univariate 
# analyses section.
# TEST 1: Does the trait have an adequate sample size?
dfAgeMaturity <- setDT(GetTraitSpecificData(dfAgeMaturity, 2))
# TEST 2: Does the trait have enough data variation?
hist(dfAgeMaturity$age_at_maturity)
range(dfAgeMaturity$age_at_maturity) 

### ECOLOGY TRAITS ###
dfEcology <- data.frame(ecology(speciesNames))
# Storing this as a file.
write.csv(dfEcology, file = "ecology_datatable.csv") 
# Read in the ecology information.
dfEcology <- fread("ecology_datatable.csv")
# Get rid of columns I do not need for the regression analysis. We are
# just doing to experiment with a few traits right now.
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, .(species_name = sciname, streams = Stream, 
                                 lakes = Lakes, feeding_type = FeedingType)]
# As some of the variables are coded as integers right now, we need to recode them to the 
# types needed for regression analysis (i.e. factor type).
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
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), 
                .SDcols = integerVars]

### Categorical traits. ###
# Binary traits.

# TRAIT: Streams.
# TEST 1: Does the trait have an adequate sample size?
dfStreams <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
nrow(dfStreams)
# TEST 2: Does the trait have enough data variation?
table(dfStreams$streams)

# TRAIT: Lakes.
# TEST 1: Does the trait have an adequate sample size?
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 3))
nrow(dfLakes)
# TEST 2: Does the trait have enough data variation?
table(dfLakes$lakes)

# Non-binary traits.
# TRAIT: Feeding Type.
# TEST 1: Does the trait have an adequate sample size?
dfFeedingType <- setDT(GetTraitSpecificData(dfEcologyTraits, 4))
# TEST 2: Does the trait have enough data variation?
# Three categories do not reach the 1% threshold and are removed.
table(dfFeedingType$feeding_type)
rareVars <- which(dfFeedingType$feeding_type == "filtering plankton" |
                    dfFeedingType$feeding_type == "other" | 
                    dfFeedingType$feeding_type == "sucking food-containing material")
dfFeedingType <- dfFeedingType[-rareVars, ]
# Also dropping these levels from the factor.
dfFeedingType$feeding_type <- droplevels(dfFeedingType$feeding_type)

# Finally, prepare the dfEcologyMV datatable by merging all univariate 
# datatables.
# Merging multiple datatables at once.
dfEcologyMV <- Reduce(function(...) merge(..., all = T), 
                      list(dfStreams, dfLakes, dfFeedingType))

# Construction of dfTraits datatable.
# This table contains all the potential traits for multivariate analysis.
# Let's first merge the trait information back to dfFiltered.
# NA/NULL/blank for those species that don't have info for that particular trait.
# Note: I only want a single row per BIN from dfFiltered for this merging 
# process (i.e. I just want the bin name, species name, and size of the bin). 
dfFilteredSingle <- dfFiltered[!duplicated(dfFiltered$bin_uri), ]
# Let's take the columns we need to construct the dfTraits datatable.
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_name, 
                                         initial_bin_size, filtered_bin_size)]
# Now merge to all of the trait MV datatables.
dfTraits <- merge(dfFilteredSingle, dfLatitude, all = T, by = "bin_uri")
# Set the keys for datatable merging.
setkey(dfTraits, species_name)
setkey(dfSpeciesGenMV, species_name)
setkey(dfEcologyMV, species_name)
dfTraits <- Reduce(function(...) merge(..., all = T), list(dfTraits, 
                                                           dfSpeciesGenMV,
                                                           dfEcologyMV,
                                                           dfAgeMaturity))
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:4, 7:15)]
colnames(dfTraits)[4] <- "filtered_bin_size"

# Merge back to dfFiltered to obtain all of the sequence information for 
# each BIN. This is for creation of the master phylogeny.
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[6] <- "species_name"
colnames(dfPreCentroid)[8] <- "initial_bin_size"
colnames(dfPreCentroid)[9] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:9)]