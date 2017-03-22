###################
# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
# correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
# determination/reference sequence trimming (lines TBD).
# Last version saved: Mar 8th 2017 (R_Pipeline_PHYLO_V7_Feb25th.R)

### PACKAGES REQUIRED ###

# For looping purposes.
#install.packages("foreach")
library(foreach)

# For efficiency purposes.
#install.packages("data.table")
library(data.table)

# For multiple sequence alignments:
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("DECIPHER")
#biocLite("muscle")
library(Biostrings)
library(DECIPHER)
library(muscle)
#install.packages("seqinr")
#install.packages("ade4")
library(seqinr)

# For tree building/PGLS.
#install.packages("ape")
library(ape)
#install.packages("phangorn")
library(phangorn)
#install.packages("caper")
library(caper)
#install.packages("geiger")
library(geiger)
#install.packages("nlme")
library(nlme)
#install.packages("expm")
#install.packages("phytools")
library(phytools)

# For BOLD and FishBase data:
#install.packages("bold")
library(bold)
#install.packages("rfishbase")
library(rfishbase)

# For stats/graphs.
#install.packages("plyr")
library(plyr)
library(lattice)
library(gtools)
#install.packages("caret")
library(caret)
#install.packages("Boruta")
library(Boruta)
library(randomForest)

# For testing:
#install.packages("rbenchmark")
#library(rbenchmark)
#install.packages("pryr")
#library(pryr)
#install.packages("devtools")
#library(devtools)
#library(compiler)

# For trait data purposes:
#install.packages("adephylo")
#install.packages("jsonlite")
#install.packages("spdep")
library(adephylo)
library(phylobase)
library(car)

# Missing data:
library(dplyr)
#install.packages("mice")
#library(mice)
#install.packages("VIM")
#library(VIM)
#install.packages("picante")
#library(picante)
#install.packages("Rphylopars")
#library(Rphylopars)

update.packages()

# Start the clock!
#ptm <- proc.time()


##### SECTION 1: DATA PROCESSING #####
# Download sequences from BOLD using the function bold_seqspec() for sequence
# and specimen data. In addition, I am only selecting those columns needed for
# downstream analysis.
dfInitial <- bold_seqspec(taxon = "Actinopterygii", 
                          geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                           "species_name", "lat", "nucleotides", 
                                           "markercode")] 

# Convert the dataframe to a datatable for time and memory saving purposes.
dfSequences <- setDT(dfInitial)
rm(dfInitial)

# Download outgroup species data from BOLD.
# Note: I am currently using the great white shark BIN AAA9092 as an outgroup. 
# This choice of outgroup species may be altered in the near future.
dfOutgroup <- bold_seqspec(bin = "BOLD:AAA4391", 
                           geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                            "species_name", "lat", "nucleotides", 
                                            "markercode")] 
dfOutgroup <- setDT(dfOutgroup)  # Again, convert to datatable.

# Combine dfOutgroup and dfSequences datatables so that they are in one useable 
# datatable.
l <- list(dfSequences, dfOutgroup)
dfFiltered <- rbindlist(l)

# To remove unneeded objects (mostly for memory saving purposes).
rm(dfSequences)
rm(dfOutgroup)
rm(l)

# Loaded TSV files from here, always. #

### FILTER 1 ###
# Filtering for presence of a BIN URI. Assignment of a BIN URI is an indicator 
# of sequence quality. BIN URI is also necessary for species assignment to the 
# BINs later on in the pipeline.
containBIN <- dfFiltered[, grep("[:]", bin_uri)]
dfFiltered <- dfFiltered[containBIN, ]
rm(containBIN)
# Modifying BIN column slightly to remove "BOLD:".
dfFiltered[, bin_uri := substr(bin_uri, 6, 13)]

### FILTER 2 ###
# Filtering for presence of a sequence.
containSeq <- dfFiltered[, grep("[CTG]", nucleotides)]
dfFiltered <- dfFiltered[containSeq, ]
rm(containSeq)

################################################################################
### TRAIT: INITIAL BIN SIZE ###
# Determine how many sequences are in a BIN in total prior to any sequence 
# filtering.
dfFiltered <- dfFiltered[, initial_bin_size := length(recordID), keyby = bin_uri]
################################################################################

# Datatable reorganization. Removing redundant columns and columns that are not
# used in the downstream analysis.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, species_name,
                             markercode, initial_bin_size, nucleotides, lat)]

### FILTER 3 ###
# Filtering for COI-5P as these are the only markers we are looking at.
setkey(dfFiltered, markercode)  # Setting a key facilitates filtering of datatables.
dfFiltered <- dfFiltered["COI-5P"]  # Now I can just specify "COI-5P" without indicating which column.

### FILTER 4 ###
# N and gap content will interfere with the multiple sequence alignment and the 
# alignment will give warning messages, so we need to trim sequences with high N
# and gap content at their terminal ends.
# First, convert nucleotides to chr type. This is necessary for regular 
# expression (regex) searches.
dfFiltered[, nucleotides := as.character(nucleotides)]
# Let's first trim large portions of Ns and gaps at the start of a sequence. 
# First, find sequences that begin (^) with an N or a gap.
startGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("^[-N]", dfFiltered$nucleotides)), length)
startGapN <- foreach(i = 1:nrow(dfFiltered)) %do%  # Loop through all of the sequences.
# If at least one sequence is found that begins with a gap or N...
  if (startGapN[[i]] > 0) { 
    # Split the sequence up using strsplit!
    # Using a regex to find gaps and/or Ns that may or may not have trailing gaps and/or Ns.
    split <- strsplit(dfFiltered$nucleotides[i], "^[-N]+")
    # Take only the second half of the element (the sequence without the gaps/Ns at the start!).
    dfFiltered$nucleotides[i] <- split[[1]][2]
  }
rm(startGapN)
# Now, let's trim large portions of Ns and gaps at the end of a sequence. 
# First, find sequences that end ($) with an N or a gap.
endGapN <- sapply(regmatches(dfFiltered$nucleotides, 
                             gregexpr("[-N]$", dfFiltered$nucleotides)), length)
endGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
  if (endGapN[[i]] > 0) {
    # Using a regex to find gaps and/or Ns that may or may not have trailing 
    # gaps and/or Ns and that *end* with one of those characters.
    split <- strsplit(dfFiltered$nucleotides[i], "[-N]+$")
    # Take only the first half of the element (the sequence without the gaps/Ns at the end!).
    dfFiltered$nucleotides[i] <- split[[1]][1]
  }
rm(endGapN)
rm(split)


### FILTER 5 ###
# Remove sequences with N/gap content above a certain threshold. In this case,
# we will be removing sequences with greater than 1% gap/N content.
# First, let's determine the number of positions where an *internal* N or gap is
# found for each sequence.
internalGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("[-N]", dfFiltered$nucleotides)), length)
# We then iterate over each sequence and see if the number of Ns or gaps is 
# greater than 1% (0.01) of the total length of the sequence.
internalGapN <- foreach(i = 1:nrow(dfFiltered)) %do% 
  which((internalGapN[[i]]/nchar(dfFiltered$nucleotides[i]) > 0.01))
# Here, we are basically "flagging" the sequences with high N/gap content. Those
# sequence will have values greater than 0.
checkGapN <- sapply(internalGapN, function (x) length(x))
checkGapN <- which(checkGapN > 0)  # Identify the sequences with high N/gap content.
# Remove these higher gap/N content sequences.
dfFiltered <- dfFiltered[-checkGapN, ]
rm(checkGapN)
rm(i)
rm(internalGapN)

### FILTER 6 ###
# Filter out sequences that are less than 640 bp and greater than 1000 bp. 
# This is because extreme long or short sequence lengths can interfere with the alignment.
# First, determine the lengths of the sequences without gaps.
seqLengths <- dfFiltered[, nchar(gsub("-", "", nucleotides))] 
# Which sequences are greater than 1000 bp and less than 640 bp in length?
seqLengthCheck <- which(seqLengths > 1000 | seqLengths < 640)
dfFiltered <- dfFiltered[-seqLengthCheck, ]
rm(seqLengths)
rm(seqLengthCheck)

# BIN Species Information. #
# Here, we are obtaining information on a per BIN basis to facilitate trait 
# matching later on in the pipeline.
### FILTER 7 ###
# Remove rows with no species information (removes BINs with no species data). 
# BINs without species data would not match with any trait information down the line.
containSpecies <- dfFiltered[, grep("[A-Z]", species_name)]
# Create a new datatable containing only sequences baring species-level identification. 
# This is necessary because NA values are considered when counting the number of species.
dfSpecies <- dfFiltered[containSpecies, ]
rm(containSpecies)
# Determine the number of unique species-level identifications assigned to each BIN. 
dfNumberOfSpecies <- dfSpecies[, .(number_of_species = length(unique(species_name))), keyby = bin_uri]
# Determine the most common species in each BIN (this can be used later down the
# road to deal with BINs with sequences from more than one species). Now, it's 
# just used to assign a species_label to sequences that have no species 
# information (giving them a "name"!).
# Note: The following line is an example of chaining. First the number of
# observations grouped by bin_uri and then by species_name is determined (.N),
# then, the number of observations are ordered and the species name with the 
# highest number of observations in each is selected.
dfSpeciesLabel <- dfSpecies[, .N, by = .(bin_uri, species_name)][order(-N), .(species_label = species_name[1L]), keyby = bin_uri]

# MERGING DATATABLES.
# Merge datatables containing BIN species information.
# Note: bin_uri is the key for each of these datatables due to the use of "keyby" 
# instead of just "by". This ultimately facilitates datatable merging.
dfSpeciesInfo <- merge(dfNumberOfSpecies, dfSpeciesLabel)
rm(dfNumberOfSpecies)
rm(dfSpeciesLabel)
rm(dfSpecies)
# Merge dfSpeciesInfo with dfFiltered to obtain all of the relevant taxonomic information.
setkey(dfFiltered, bin_uri)
dfFiltered <- merge(dfFiltered, dfSpeciesInfo)
rm(dfSpeciesInfo)


### FILTER 8 ###
# Filtering for BINs that contain only ONE species name. This removes BINs with
# species level conflicts.
# EDIT: I could reorder this to before species_label so it is more straightforward.
dfFiltered <- dfFiltered[number_of_species == 1]

################################################################################
### TRAIT: POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfFiltered[, filtered_bin_size := length(recordID), by = bin_uri]
################################################################################

################################################################################
### TRAIT: MEDIAN LATITUDE/LATITUDINAL RANGE ###
# Currently, median latitude and latitudinal range are the only traits whose 
# information is being taken from BOLD. The rest of the data will be obtained 
# from FishBase.

# Filtering for presence of a latitude value.
containLat <- dfFiltered[, grep("[0-9]", lat)]
dfLatitudeSpecies <- dfFiltered[containLat, ]
rm(containLat)
# Convert the latitude (lat) column to number instead of character type. This is
# necessary for median and range calculations.
dfLatitudeSpecies[, lat_num := as.numeric(lat)]
# Conversion to absolute values before median latitude values are calculated.
dfLatitudeSpecies[, abs_lat_num := abs(lat_num)]
# Determine a median latitude for each BIN using absolute values.
dfLatitudeSpecies[, median_lat := median(abs_lat_num), keyby = bin_uri]
# Find the minimum latitude.
dfLatitudeSpecies[, lat_min := min(lat_num), keyby = bin_uri]
# Find the maximum latitude.
dfLatitudeSpecies[, lat_max := max(lat_num), keyby = bin_uri]
# Determine a latitude range for each BIN.
dfLatitudeSpecies[, lat_range := abs(lat_max - lat_min), keyby = bin_uri]
################################################################################

# Datatable reorganzation for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, species_name = species_label, 
                             nucleotides, initial_bin_size, filtered_bin_size)]
# Building a datatable containing nucleotide data required for a "master" phylogeny.
#dfSpeciesForTree <- dfLatitudeSpecies[, c(1:4, 10, 7, 6, 11)]


################################################################################
### SECTION 2: FISHBASE TRAITS ###
# In this section, traits from FishBase are extracted using the rfishbase package
# and matched against the information obtained from BOLD.

# Extract all of the species names that are available on FishBase.
allFish <- fishbase
allFish$fullName <- paste(allFish$Genus, allFish$Species)
fishBaseSpecies <- allFish$fullName  # 33104 species names.
rm(allFish)
# Match the species labels from BOLD with the species names from FishBase.
# Make into a dataframe first so they can be merged.
dfFishBaseSpecies <- data.table(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
rm(fishBaseSpecies)
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by = "species_name")
rm(dfFishBaseSpecies)
# Extract species' name as a vector if trying to access trait information
# for first time (aka haven't saved trait info in CWD yet).
speciesNames <- dfBoldBase[, unique(species_name)]
rm(dfBoldBase)


### TRAIT ASSIGNMENT AND RECODING SECTION ###
# As there are multiple entries per species, I want to take the median or mode 
# value of some traits. This will depend on the type of trait (i.e. continuous 
# vs. categorical). If necessary, the trait data is recoded with data that is 
# useable in a downstream regression model. For instance, categorical data is 
# replaced with ranked or binary data. For instance, categorical traits that are 
# assumed to beget higher rates of molecular evolution are assigned a higher 
# rank. 

### ECOSYSTEM TRAITS ###
# The entries in this datatable are from different references/sources.
#dfEcosystem <- data.frame(ecosystem(speciesNames))
# Storing this as a file.
#write.csv(dfEcosystem, file = "ecosystem_information.csv")
# Read in the ecosystem information.
dfEcosystem <- fread("ecosystem_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfEcosystemTraits <- dfEcosystem[, .(species_name = sciname, status = Status, 
                                     ecosystem_type = EcosystemType, 
                                     salinity = Salinity, 
                                     average_depth = AverageDepth,
                                     max_depth = MaxDepth,
                                     temp_surface = TempSurface,
                                     temp_depth = TempDepth,
                                     climate = Climate)]

# Mode trait(s).
# TRAIT: Ecosystem Type.
# A count column is created to count the number of ecosystem_type observations per species.
dfEcosystemTraits[, count := .N, by = .(species_name, ecosystem_type)]
# The ecosystem_type with the highest number of observations per species is then selected.
# This data replaces the original data.
dfEcosystemTraits[order(-count), ecosystem_type := ecosystem_type[1L], by = species_name]
# Recode ecosystem_type for regression analyses.
# Change the data to factor type.
dfEcosystemTraits[, ecosystem_type := as.factor(ecosystem_type)]
# Five-level factor type trait.
# EDIT: Need error/flag statements when a value isn't present
dfEcosystemTraits[, ecosystem_type := revalue(ecosystem_type, 
                                              c("Sea/Bay/Gulf" = "1",
                                                "River (basin)" = "2",
                                                "Lake" = "3",
                                                "Lagoon" = "4",
                                                "Seamount" = "5",
                                                "Zoogeographic realm" = NA))]

# Univariate section.
GetTraitSpecificData <- function(x, y) {
  # Filters a dataframe for only those data related to a specified trait.
  # x = dataframe of species and trait information.
  # y = trait
  
  x <- as.data.frame(x)  # Still want to figure this out for datatable.
  # Find rows without data for column.
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  # Construct the univariate trait datatable. This datatable will be used in the eventually 
  # univariate analysis.
  # If there are rows without data for column...
  if (length(noY) > 0) {
    # Remove the species without data for column.
    z <- x[-noY, ]
    # Reorganize datatable.
    z <- z[c(1, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$species_name), ] 
    # If all rows have data for column...
  } else {
    # If no entries need to be removed, just rename and reorganize x.
    z <- x[c(1, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$species_name), ]  
  }
  rm(noY)
  return(z)
}
# # Filtering for presence of ecosystem_type data. This is for the univariate analyses section.
# First, does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfEcosystemType <- setDT(GetTraitSpecificData(dfEcosystemTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: Lagoon and Seamount categories too rare so they will be removed.
table(dfEcosystemType$ecosystem_type)
rareVars <- which(dfEcosystemType$ecosystem_type == 4 | dfEcosystemType$ecosystem_type == 5)
dfEcosystemType <- dfEcosystemType[-rareVars, ]
# Also dropping these levels from the factor.
dfEcosystemType$ecosystem_type <- droplevels(dfEcosystemType$ecosystem_type)
# Create a master species list for eventually master phylogeny.
dfMasterSpecies <- dfEcosystemType[, 1]


# TRAIT: Environmental salinity level.
dfEcosystemTraits[, count := .N, by = .(species_name, salinity)]
dfEcosystemTraits[order(-count), env_salinity_level := salinity[1L], by = species_name]
dfEcosystemTraits$count <- NULL
# Recode env_salinity_level for regression analyses.
# Change the data to factor type.
dfEcosystemTraits[, env_salinity_level := as.factor(env_salinity_level)]
# Three-level factor type trait.
# EDIT: Need error/flag statements when a value isn't present
dfEcosystemTraits[, env_salinity_level := revalue(env_salinity_level, 
                                                  c("freshwater" = "0", 
                                                    "saltwater" = "1"))]
# Filtering for presence of env_salinity_level data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfEnvironSalLevel <- setDT(GetTraitSpecificData(dfEcosystemTraits, 10))
# TEST 2: Does the trait have enough data variation?
# Answer: Brackish water was rare and was removed.
table(dfEnvironSalLevel$env_salinity_level)
rareVars <- which(dfEnvironSalLevel$env_salinity_level == "brackish")
dfEnvironSalLevel <- dfEnvironSalLevel[-rareVars, ]
# Also dropping this levels from the factor.
dfEnvironSalLevel$env_salinity_level <- droplevels(dfEnvironSalLevel$env_salinity_level)
# Also remove the "brackish" level.
# I will be adding unique species to the master phylogeny a lot.
AppendToMasterSpeciesDf <- function(trait_data, df_to_append_to) {
  # Collect unique species.
  dfAddSpecies <- trait_data[, 1]
  # Merge back to master species list.
  df_to_append_to <- rbind(df_to_append_to, dfAddSpecies)
  # Keep only unique species.
  df_to_append_to <- df_to_append_to[!duplicated(df_to_append_to$species_name), ]
  # Return the updated master species df.
  return(df_to_append_to)
}
# Collect any new species from dfEnvironSalLevel.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfEnvironSalLevel, dfMasterSpecies)


# TRAIT: Climate.
# First, revalue the category names so the syntax is consistent (i.e. "Polar" should be "polar").
dfEcosystemTraits[, climate := revalue(climate, c("Boreal" = "boreal", "Polar" = "polar",
                                                  "Subtropical" = "subtropical", 
                                                  "Temperate" = "temperate",
                                                  "Tropical" = "tropical"))]
dfEcosystemTraits[, count := .N, by = .(species_name, climate)]
dfEcosystemTraits[order(-count), climate := climate[1L], by = species_name]
dfEcosystemTraits$count <- NULL
# Recode env_salinity_level for regression analyses.
# Change the data to factor type.
dfEcosystemTraits[, climate := as.factor(climate)]
# Filtering for presence of climate data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: 
dfClimate <- setDT(GetTraitSpecificData(dfEcosystemTraits, 9))
# TEST 2: Does the trait have enough data variation?
# Answer: Boreal is too rare and is removed.
table(dfClimate$climate)
rareVars <- which(dfClimate$climate == "boreal")
dfClimate <- dfClimate[-rareVars, ]
#Also dropping this levels from the factor.
dfClimate$climate <- droplevels(dfClimate$climate)
# Collect any new species from dfClimate.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfClimate, dfMasterSpecies)



# Median trait(s).
# TRAIT: Average depth.
# Note: These values don't seem to match up with those found on fish base.
# The column average_depth must first be converted to double (numeric) type.
dfEcosystemTraits[, average_depth := as.double(average_depth)]
# The median average_depth is then determined for each species.
dfEcosystemTraits[, average_depth := median(average_depth, na.rm = TRUE), keyby = species_name]
# Filtering for presence of average depth data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfAvgDepth <- setDT(GetTraitSpecificData(dfEcosystemTraits, 5))
# TEST 2: Does the trait have enough data variation?
# Answer: Continuous trait. Limited variation but this can be noted in discussion.
hist(dfAvgDepth$average_depth)
# Collect any new species from dfAvgDepth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfAvgDepth, dfMasterSpecies)


# TRAIT: Maximum depth.
# The column max_depth must first be converted to double (numeric) type.
dfEcosystemTraits[, max_depth := as.double(max_depth)]
# The median max_depth is then determined for each species.
dfEcosystemTraits[, max_depth := median(max_depth, na.rm = TRUE), keyby = species_name]
# Filtering for presence of max depth data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfMaxDepth <- setDT(GetTraitSpecificData(dfEcosystemTraits, 6))
# TEST 2: Does the trait have enough data variation?
# Answer: Continuous trait. Limited variation but this can be noted in discussion.
hist(dfMaxDepth$max_depth)
# Collect any new species from dfMaxDepth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfMaxDepth, dfMasterSpecies)


# TRAIT: Surface temperature.
# The column temp_surface must first be converted to double (numeric) type.
dfEcosystemTraits[, temp_surface := as.double(temp_surface)]
# The median temp_surface is then determined for each species.
dfEcosystemTraits[, temp_surface := median(temp_surface, na.rm = TRUE), keyby = species_name]
# Filtering for presence of max depth data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfTempSurface <- setDT(GetTraitSpecificData(dfEcosystemTraits, 7))
# TEST 2: Does the trait have enough data variation?
# Answer: Continuous trait. Limited variation but this can be noted in discussion.
hist(dfTempSurface$temp_surface)
# Collect any new species from dfTempSurface.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfTempSurface, dfMasterSpecies)

# TRAIT: Depth temperature.
# The column temp_depth must first be converted to double (numeric) type.
dfEcosystemTraits[, temp_depth := as.double(temp_depth)]
# The median temp_depth is then determined for each species.
dfEcosystemTraits[, temp_depth := median(temp_depth, na.rm = TRUE), keyby = species_name]
# Filtering for presence of max depth data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfTempDepth <- setDT(GetTraitSpecificData(dfEcosystemTraits, 8))
# TEST 2: Does the trait have enough data variation?
# Answer: Continuous trait. Limited variation but this can be noted in discussion.
hist(dfTempDepth$temp_depth)
# Collect any new species from dfTempDepth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfTempDepth, dfMasterSpecies)



# Special cases.
# These are traits for which additional steps are required to pick the representative value.
# INTRODUCED VS. NOT INTRODUCED.
#DetermineIfSpeciesIntroduced <- function(x) {
 # answer <- NULL
  #containsIntroduced <- length(grep("introduced", x))
  # If the species has ever been introduced in an area...
  #if (containsIntroduced > 0) {
    # It will be assigned "Yes".
   # answer <- "Yes" 
  #} else {
    # Else, it will be assigned "No".
   # answer <- "No"
  #} 
  #return(answer)
#}

# Apply the DetermineIfSpeciesIntroduced function to the appropriate column in the datatable.
#dfEcosystemTraits[, introduced := DetermineIfSpeciesIntroduced(status), by = species_name]

# ENDEMICITY.
DetermineIfSpeciesEndemic <- function(x) {
  answer <- NULL
  containsEndemic <- length(grep("endemic", x))
  # If the species has ever been introduced in an area..
  if (containsEndemic > 0) {
    # It will be assigned "Yes".
    answer <- "Yes" 
  } else {
    # Else, it will be assigned "No".
    answer <- "No"
  } 
  return(answer)
}

# Apply the DetermineIfSpeciesEndemic function to the appropriate column in the datatable.
dfEcosystemTraits[, endemic := DetermineIfSpeciesEndemic(status), by = species_name]
# Filtering for presence of endemicity data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfEndemic <- setDT(GetTraitSpecificData(dfEcosystemTraits, 11))
# TEST 2: Does the trait have enough data variation?
# Answer: PASS.
table(dfEndemic$endemic)
# Collect any new species from dfTempDepth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfEndemic, dfMasterSpecies)
# Don't need this column anymore.
dfEcosystemTraits$status <- NULL


# SALINITY TOLERANCE
DetermineSalinityTolerance <- function(x) {
  answer <- NULL
  typesOfWater <- unique(x)
  number_of_water_types <- length(unique(x))
  # If species is ever present in more than 1 type of salinity level...
  if (number_of_water_types > 1) {
    # It will be assigned "euryhaline".
    answer <- "euryhaline"
    # Else, it will be assigned "stenohaline".
  } else if (number_of_water_types == 1) {
    answer <- "stenohaline"
  } else {
    answer <- NA
  }
  return(answer)
}

# Apply the DetermineSalinityTolerance function to the appropriate column in the datatable.
dfEcosystemTraits[, salinity_tolerance := DetermineSalinityTolerance(salinity), by = species_name]
# Remove salinity column now as it is no longer needed.
dfEcosystemTraits$salinity <- NULL
# Recode env_salinity_level for regression analyses.
# Change the data to factor type.
dfEcosystemTraits[, salinity_tolerance := as.factor(salinity_tolerance)]
# Two-level factor type trait.
# EDIT: Need error/flag statements when a value isn't present
dfEcosystemTraits[, salinity_tolerance := revalue(salinity_tolerance, 
                                                  c("stenohaline" = "0", 
                                                    "euryhaline" = "1"))]
# Filtering for presence of salinity tolerance data. This is for the univariate analyses section.
# 500 SPECIES TEST: PASS 
dfSalTolerance <- setDT(GetTraitSpecificData(dfEcosystemTraits, 10))
# TEST 2: Does the trait have enough data variation?
# Answer: PASS.
table(dfSalTolerance$salinity_tolerance)
# Collect any new species from dfEnvironSalLevel.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfEnvironSalLevel, dfMasterSpecies)

# Finally, prepare the dfEcosystemTraits datatable for the multivariate analysis
# by taking only a single row per species.
dfEcosystemTraits <- dfEcosystemTraits[!duplicated(dfEcosystemTraits$species_name), ] 



### MORPHOMETRICS TRAITS ###
#dfMorphometrics <- data.frame(morphometrics(speciesNames))
# Storing this as a file.
#write.csv(dfMorphometrics, file = "morphometrics_information.csv") 
# Read in the ecology information.
dfMorphometrics <- fread("morphometrics_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfMorphometricTraits <- dfMorphometrics[, .(species_name = sciname, total_length = TL, standard_length = SL, 
                                            fork_length = FL, head_length = HL, body_depth = BD, 
                                            aspect_ratio = AspectRatio)]
# As the morphometric information is entered on an individual basis, it is 
# necessary to find the median trait for each species. The following are 
# morphometric traits that I am including in the regression analysis.
# First, change all of the traits to numeric type.
changeVars <- c(2:6)
dfMorphometricTraits[, (changeVars) := lapply(.SD, as.double), .SDcols = changeVars]

# Median traits.
# Since I want to apply the function median() to multiple columns I can do it in one statement.
# TRAIT: Total length.
medianVars <- c(2:7) # TL, SL, FL, HL, BD, AR columns.
dfMorphometricTraits[, (medianVars) := lapply(.SD, function(x) median(x, na.rm = TRUE)),
                     by = species_name, .SDcols = medianVars]
# I only want one row per species.
dfMorphometricTraits <- dfMorphometricTraits[!duplicated(dfMorphometricTraits$species_name), ] 
# TRAIT: Total length.
# Filtering for presence of total length data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfTotalLength <- setDT(GetTraitSpecificData(dfMorphometricTraits, 2))
# TEST 2: Does the trait have enough data variation?
hist(dfTotalLength$total_length)
# Collect any new species from dfTotalLength.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfTotalLength, dfMasterSpecies)

# TRAIT: Standard length.
# Filtering for presence of standard length data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfStandardLength <- setDT(GetTraitSpecificData(dfMorphometricTraits, 3))
# TEST 2: Does the trait have enough data variation?
hist(dfStandardLength$standard_length)
# Collect any new species from dfStandardLength.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfStandardLength, dfMasterSpecies)

# TRAIT: Fork length.
# Filtering for presence of fork length data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfForkLength <- setDT(GetTraitSpecificData(dfMorphometricTraits, 4))
# TEST 2: Does the trait have enough data variation?
hist(dfForkLength$fork_length)
# Collect any new species from dfForkLength.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfForkLength, dfMasterSpecies)

# TRAIT: Head length.
# Filtering for presence of head length data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfHeadLength <- setDT(GetTraitSpecificData(dfMorphometricTraits, 5))
# TEST 2: Does the trait have enough data variation?
hist(dfHeadLength$head_length)
# Collect any new species from dfHeadLength.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfHeadLength, dfMasterSpecies)

# TRAIT: Body depth.
# Filtering for presence of body depth data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfBodyDepth <- setDT(GetTraitSpecificData(dfMorphometricTraits, 6))
# TEST 2: Does the trait have enough data variation?
hist(dfBodyDepth$body_depth)
# Collect any new species from dfBodyDepth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfBodyDepth, dfMasterSpecies)

# TRAIT: Aspect Ratio.
# Filtering for presence of aspect ratio data. This is for the univariate analyses section.
# Does the trait pass the tests for inclusion in analyses?
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes! 
dfAspectRatio <- setDT(GetTraitSpecificData(dfMorphometricTraits, 7))
# TEST 2: Does the trait have enough data variation?
hist(dfAspectRatio$aspect_ratio)
# Collect any new species from dfAspectRatio.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfAspectRatio, dfMasterSpecies)

# IDEA: Apply these functions across dataframe. Esepcially important for eclogy
# traits since there are so many!


### MORPHOLOGY TRAITS ###
#dfMorphology <- data.frame(morphology(speciesNames))
# Storing this as a file.
#write.csv(dfMorphology, file = "morphology_information.csv") 
# Read in the ecology information.
dfMorphology <- fread("morphology_information.csv")
# Get rid of columns I do not need for the regression analysis.
dfMorphologyTraits <- dfMorphology[, .(species_name = sciname, body_shape_I = BodyShapeI, 
                                       body_shape_II = BodyShapeII, operculum_present = OperculumPresent, 
                                       pos_of_mouth = PosofMouth, type_of_scales = TypeofScales)]
# These traits are unlikely to differ between different stocks of fish so I can 
# take a single row per species.
dfMorphologyTraits <- dfMorphologyTraits[!duplicated(dfMorphologyTraits$species_name), ] 

# Recoding dfMorphologyTraits.
# Converting to factor type as these are discrete traits.
changeVars <- c(2:6)
dfMorphologyTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]
# Also, recode operculum_present.
dfMorphologyTraits[, operculum_present := revalue(operculum_present, c("-1" = "1"))]

# TRAIT: Body Shape I.
# Filtering for presence of body shape I data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBodyShapeI <- setDT(GetTraitSpecificData(dfMorphologyTraits, 2))
# TEST 2: Does the trait have enough data variation?
# Answer: "other" is too rare and will be removed.
table(dfBodyShapeI$body_shape_I)
rareVars <- which(dfBodyShapeI$body_shape_I == "other")
dfBodyShapeI <- dfBodyShapeI[-rareVars, ]
#Also dropping this levels from the factor.
dfBodyShapeI$body_shape_I <- droplevels(dfBodyShapeI$body_shape_I)
# Collect any new species from dfBodyShapeI.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfBodyShapeI, dfMasterSpecies)

# TRAIT: Body Shape II.
# Filtering for presence of body shape II data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBodyShapeII <- setDT(GetTraitSpecificData(dfMorphologyTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: "other" and "angular" are too rare and will be removed.
table(dfBodyShapeII$body_shape_II)
rareVars <- which(dfBodyShapeII$body_shape_II == "other (see Diagnosi" | dfBodyShapeII$body_shape_II == "angular")
dfBodyShapeII <- dfBodyShapeII[-rareVars, ]
#Also dropping this levels from the factor.
dfBodyShapeII$body_shape_II <- droplevels(dfBodyShapeII$body_shape_II)
# Collect any new species from dfBodyShapeII.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfBodyShapeII, dfMasterSpecies)

# TRAIT: Operculum Present.
# Filtering for presence of operculum present data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfOperPresent <- setDT(GetTraitSpecificData(dfMorphologyTraits, 4))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfOperPresent$operculum_present)
# Collect any new species from dfOperPresent.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfOperPresent, dfMasterSpecies)

# TRAIT: Position of Mouth.
# Filtering for presence of position of mouth data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfPosMouth <- setDT(GetTraitSpecificData(dfMorphologyTraits, 5))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfPosMouth$pos_of_mouth)
# Collect any new species from dfPosMouth.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfPosMouth, dfMasterSpecies)

# TRAIT: Type of Scales.
# First, revalue trait so "cycloid" is equivalent to "cycloid scales".
dfMorphologyTraits[, type_of_scales := revalue(type_of_scales, c("cycloid" = "cycloid scales"))]
# Filtering for presence of type of scales data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfScaleType <- setDT(GetTraitSpecificData(dfMorphologyTraits, 6))
# TEST 2: Does the trait have enough data variation?
# Answer: "other" and "rhombic scales" are too rare and are removed.
table(dfScaleType$type_of_scales)
rareVars <- which(dfScaleType$type_of_scales == "other (see remark)" | dfScaleType$type_of_scales == "rhombic scales")
dfScaleType <- dfScaleType[-rareVars, ]
# Also dropping this levels from the factor.
dfScaleType$type_of_scales <- droplevels(dfScaleType$type_of_scales)
# Collect any new species from dfScaleType.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfScaleType, dfMasterSpecies)



### ECOLOGY TRAITS ###
dfEcology <- data.frame(ecology(speciesNames))
# Storing this as a file.
write.csv(dfEcology, file = "ecology_information.csv") 
# Read in the ecology information.
dfEcology <- fread("ecology_information.csv")
colnames(dfEcology)[3] <- "species_name"
# Get rid of columns I do not need for the regression analysis. 
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, c(3, 7:27, 29, 31, 33, 39, 52:57, 66:78,
                                 80:88, 90:110, 117:118, 121)]
# Recode ecology traits to the types needed for regression analysis.
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
# Recode the integer variables.
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), .SDcols = integerVars]
rm(integerVars)
rm(characterVars)
rm(changeVars)

# Removal of near zero variance traits. These traits are not informative enough for
# regression analyses.
dfEcologyTraits <- as.data.frame(dfEcologyTraits)
dfEcologyTraits <- dfEcologyTraits[, -nearZeroVar(dfEcologyTraits)]

# TRAIT: Stream.
# Filtering for presence of type of stream data. This is for the univariate analyses section.
dfStream <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
# Collect any new species from dfScaleType.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfScaleType, dfMasterSpecies)

# TRAIT: Lakes.
# Filtering for presence of type of lake data. This is for the univariate analyses section.
dfLake <- setDT(GetTraitSpecificData(dfEcologyTraits, 3))
# Collect any new species from dfLake.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfLake, dfMasterSpecies)

# TRAIT: Herbivory.
# Filtering for presence of type of herbivory data. This is for the univariate analyses section.
dfHerbivory <- setDT(GetTraitSpecificData(dfEcologyTraits, 4))
# Collect any new species from dfHerbivory.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfHerbivory, dfMasterSpecies)

# TRAIT: Feeding Type.
# Filtering for presence of type of feeding type data. This is for the univariate analyses section.
dfFeedingType <- setDT(GetTraitSpecificData(dfEcologyTraits, 5))
# Collect any new species from dfFeedingType.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfFeedingType, dfMasterSpecies)

# TRAIT: Diet Troph.
# Filtering for presence of type of diet troph data. This is for the univariate analyses section.
dfDietTroph <- setDT(GetTraitSpecificData(dfEcologyTraits, 6))
# Collect any new species from dfDietTroph.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfDietTroph, dfMasterSpecies)

# TRAIT: Food Troph.
# Filtering for presence of type of diet troph data. This is for the univariate analyses section.
dfFoodTroph <- setDT(GetTraitSpecificData(dfEcologyTraits, 7))
# Collect any new species from dfFoodTroph.
dfMasterSpecies <- AppendToMasterSpeciesDf(dfFoodTroph, dfMasterSpecies)


# POTENTIAL NEW TRAITS.

### STATUS RELATED ###
# Introductions.
dfIntroductions <- data.frame(introductions(speciesNames))
write.csv(dfIntroductions, file = "introductions_information.csv") 
# Read in the introductions information.
dfIntroductions <- fread("introductions_information.csv")
colnames(dfIntroductions)[5] <- "species_name"

# Reason, Estabwild, EstabType, EstabAqua, ReproMode, EcolEff, EcolEffType, Impacts, Invasive




### DIET RELATED ###
# Diet data.
dfDiet <- data.frame(diet(speciesNames))
write.csv(dfDiet, file = "diet_information.csv") 
# Read in the diet information.
dfDiet <- fread("diet_information.csv")
colnames(dfDiet)[4] <- "species_name"

# Troph, FoodI


### LIFE HISTORY RELATED ###
# Fecundity.
dfFecundity <- data.frame(fecundity(speciesNames))
write.csv(dfFecundity, file = "fecundity_information.csv") 
# Read in the fecundity information.
dfFecundity <- fread("fecundity_information.csv")
colnames(dfFecundity)[4] <- "species_name"

# FecundityMin, WeightMin, LengthFecunMin, FecundityMax, WeightMax, LengthFecunMax, LengthTypeFecMac, SpawningCycles


# Maturity.
dfMaturity <- data.frame(maturity(speciesNames))
write.csv(dfMaturity, file = "maturity_information.csv") 
# Read in the maturity information.
dfMaturity <- fread("maturity_information.csv")
colnames(dfMaturity)[3] <- "species_name"

# Sex, AgeMatMin (lower range), AgeMatMin2 (upper range), tm (age at maturity)



# Reproduction.
dfReproduction <- data.frame(reproduction(speciesNames))
write.csv(dfReproduction, file = "reproduction_information.csv") 
# Read in the reproduction information.
dfReproduction <- fread("reproduction_information.csv")
colnames(dfReproduction)[4] <- "species_name"

# ReproMode, Fertilization, MatingSystem, Spawning, BatchSpawner, RepGuild1, RepGuild2,
# ParentalCare, ParentalCareQ


# Spawning.
dfSpawning <- data.frame(spawning(speciesNames))
write.csv(dfSpawning, file = "spawning_information.csv") 
# Read in the spawning information.
dfSpawning <- fread("spawning_information.csv")
colnames(dfSpawning)[3] <- "species_name"

# SpawningGround, TempLow, TempHigh, WeightMin, WeightMax




### MORPHOMETRIC RELATED ###
# Popchar.
dfPopChar <- data.frame(popchar(speciesNames))
write.csv(dfPopChar, file = "popchar_information.csv") 
# Read in the popchar information.
dfPopChar <- fread("popchar_information.csv")
colnames(dfPopChar)[3] <- "species_name"

# Wmax, TypeWeight, Lmax, Type, tmax, 


### PHYSIOLOGY RELATED ###
# Oxygen.
dfOxygen <- data.frame(oxygen(speciesNames))
write.csv(dfOxygen, file = "oxygen_information.csv") 
# Read in the oxygen information.
dfOxygen <- fread("oxygen_information.csv")
colnames(dfOxygen)[4] <- "species_name"

# Weight, Temp, Oxygen, Oxygenmgl, OxygenCons, MetabolicLevel


# Speed.
dfSpeed <- data.frame(speed(speciesNames))
write.csv(dfSpeed, file = "speed_information.csv") 
# Read in the speed information.
dfSpeed <- fread("speed_information.csv")
colnames(dfSpeed)[4] <- "species_name"

# Speedms, Mode


# Swimming.
dfSwimming <- data.frame(swimming(speciesNames))
write.csv(dfSwimming, file = "swimming_information.csv") 
# Read in the swimming information.
dfSwimming <- fread("swimming_information.csv")
colnames(dfSwimming)[2] <- "species_name"

# AdultType, AdultMode





# Construction of dfTraits datatable.
# This table contains all potential traits for multivariate analysis.
# Let's merge the trait information back to dfFiltered.
# NA/NULL/blank for those species that don't have info for that particular trait.
# I only want a single row per BIN for this merging process.
dfFilteredSingle <- dfFiltered[!duplicated(dfFiltered$bin_uri), ]
# Let's take the columns we need to construct the dfTraits datatable.
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_label, median_lat, lat_range, 
                                         initial_bin_size, filtered_bin_size)]
# Merge dfFiltered with dfEcosystemTraits to create dfTraits.
dfTraits <- merge(dfFilteredSingle, dfEcosystemTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
rm(dfFilteredSingle)
# Merge dfTraits with dfMorphometricTraits.
dfTraits <- merge(dfTraits, dfMorphometricTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
# Merge dfTraits with dfMorphologyTraits.
dfTraits <- merge(dfTraits, dfMorphologyTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
# Merge dfTraits with dfEcologyTraits.
dfTraits <- merge(dfTraits, dfEcologyTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)

### COMPLETE CASES ###
# Let's try to maximize the amount of complete cases in our dataset (species that have no missing data). #
# First, order the columns in dfTraits by the amount of missing data (NA values).
dfTraitsNA <- sort(dfTraits[, lapply(.SD, function(x) sum(is.na(x)))])
# Reorder the original dfTraits. The columns with the least amount of NA values are now first.
setcolorder(dfTraits, names(dfTraitsNA))
# Now, count the number of NAs in each row.
dfTraits[, count := rowSums(is.na(dfTraits))]
# Sort the rows by this order.
dfTraits <- dfTraits[order(count), ]
# Remove this count column as it is no longer needed.
dfTraits[, count := NULL]
# Removing the first three columns and creating a new datatable called dfTemp. 
# These columns do not need to be considered in the complete cases as I am only looking at traits.
dfTemp <- dfTraits[, 3:123]

# Write the trait data in order to completeness to file
#write.csv(dfTraits, file = "TraitDataCompleteness.csv")

# The datatable is now sorted: columns with the least amount of NAs to the left, 
# rows with with the least amount of NAs on top.
# Now I want to loop through the traits, removing one column (trait) at a time
# and count the number of complete cases. This will provide us some information
# as to which traits would provide an adequate sample size for downstream analysis.
# First, take the number of columns in dfTemp.
len <- ncol(dfTemp)
# Create a temporary variable to hold this number. This variable will hold the 
# number of subsets to check at each iteration.
tempLen <- len
# Create a vector to hold the results of the loop.
all.cc <- NULL
# Start the loop:
for (i in 1:len) {
  # Works best if you set dfTemp back to a dataframe.
  x <- as.data.frame(dfTemp)
  # x is the placeholder dataframe in the loop.
  x <- x[, 1:tempLen]
  # Determine which rows are "complete" using the "tempLen" subset of traits.
  x <- complete.cases(x)
  # Complete rows of data will be "TRUE".
  x <- which(x == "TRUE")
  # Find the number of complete cases.
  x <- length(x)
  # Add it to the all.cc variable that's holding all of the results of the loop.
  all.cc[i] <- x
  # Minus 1 from tempLen so we can check the next subset of traits.
  tempLen <- tempLen - 1
}

# Decide where to cut the datatable. (i.e. Pick an adequate subset of traits that
# maximize sample size).
# First, name it according to the trait columns.
names(all.cc) <- rev(colnames(dfTemp))
# Look at the results.
all.cc
# fork_length seems like a good cutoff point.
which(colnames(dfTraits) == "fork_length")
dfTraits <- dfTraits[, 1:109] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfCompleteCases <- dfTraits %>% filter(complete.cases(.))

# Removal of near zero variance traits. These are traits that are invariant and 
# that wouldn't contribute anything to the multivariate analysis.
dfCompleteCases <- dfCompleteCases[, -nearZeroVar(dfCompleteCases)]

# Merge back to dfFiltered to obtain all of the sequence information for each BIN.
dfPreCentroid <- merge(dfFiltered, dfCompleteCases, by = "bin_uri")
rm(dfFiltered)
rm(dfCompleteCases)
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[6] <- "initial_bin_size"
colnames(dfPreCentroid)[8] <- "median_lat"
colnames(dfPreCentroid)[9] <- "lat_range"
colnames(dfPreCentroid)[11] <- "species_name"
colnames(dfPreCentroid)[12] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:3, 11, 6, 12, 7:9, 18:46)]


### SECTION 3: CENTROID SEQUENCE DETERMINATION ###
# Centroid Sequence: BIN sequence with minimum average pairwise distance to all 
# other sequences in a given BIN.

# Subset dataframe to find BINs with more than one sequence.These are the BINs 
# that will require centroid sequences.
largeBins <- which(dfPreCentroid$filtered_bin_size > 1)
  
# If there is at least one BIN with more than one member, then a dataframe 
# dfCentroidSeqs will be created with those BINs.
if (length(largeBins) > 0) {
  
  # Subset out the larger BINs.
  dfCentroidSeqs <- dfPreCentroid[largeBins, ]
  
  # Find the number of unique BINs in dfCentroidSeqs. 
  binNumberCentroid <- unique(dfCentroidSeqs$bin_uri)
  binNumberCentroid <- length(binNumberCentroid)

  # We also have to create another separate dataframe with BINs that only have 
  # one member called dfSingletons.
  dfSingletons <- dfPreCentroid[-largeBins, ]
  rm(largeBins)
    
  # We then take the dfCentroidSeqs sequences and break it down into a list with 
  # each element being a unique BIN. 
  largeBinList <- lapply(unique(dfCentroidSeqs$bin_uri), 
                        function(x) dfCentroidSeqs[dfCentroidSeqs$bin_uri == x,])
    
  # Extract record id from each BIN.
  largeBinRecordId <- lapply(largeBinList, function(x) (x$recordID))
  
  # Convert all of the sequences in the largeBinList to dnaStringSet format for 
  # the alignment step.
  DNAStringSet1 <- lapply(largeBinList, function(x) DNAStringSet(x$nucleotides))
  rm(largeBinList)
  
  # Name DNAStringSet with the record ids.
  for (i in seq(from = 1, to = binNumberCentroid, by = 1)) {
    names(DNAStringSet1[[i]]) <- largeBinRecordId[[i]]
  }
  rm(largeBinRecordId)
  
  # Alignment using MUSCLE.
  alignment1 <- lapply(DNAStringSet1, function(x) 
    muscle::muscle(x, diags = TRUE, gapopen = -3000))

  # We can then convert each alignment to DNAbin format.
  dnaBINCentroid <- lapply(alignment1, function(x) as.DNAbin(x))
  rm(alignment1)
    
  # Then, we perform genetic distance determination with the TN93 model on each 
  # DNAbin list.
  geneticDistanceCentroid <- lapply(dnaBINCentroid, function(x) 
    dist.dna(x, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE))
  rm(dnaBINCentroid)
    
  # The centroid sequence can be determined from the distance matrix alone.
  # It is the sequence in a BIN with minimum average pairwise distance to all 
  # other sequences in its BIN.
  centroidSeq <- lapply(geneticDistanceCentroid, function(x) 
                        which.min(rowSums(x)))
  rm(binNumberCentroid)
  rm(geneticDistanceCentroid)
  centroidSeq <- unlist(centroidSeq)
  centroidSeq <- names(centroidSeq)
  centroidSeq <- as.numeric(centroidSeq)
    
  # Subset dfCentroidSeqs by the record ids on this list.
  dfCentroidSeqs <- subset(dfCentroidSeqs, recordID %in% centroidSeq)
  rm(centroidSeq)
    
  # Append our singletons to our centroid sequences. We will make a new 
  # dataframe for this - dfAllSeqs. 
  # Now we have all of the sequences we need for the next alignment of all 
  # sequences, with one sequence representing each BIN.
  dfAllSeqs <- rbind(dfCentroidSeqs, dfSingletons)
  rm(dfCentroidSeqs)
  rm(dfSingletons)
    
} else {
    
  dfAllSeqs <- dfPreCentroid
    
}

# Make sure there is only a single row per BIN.
dfAllSeqs <- dfAllSeqs[!duplicated(dfAllSeqs$bin_uri), ] 

# For testing purposes of refSeqTrim.
data <- dfAllSeqs

################################################################################
# REFERENCE SEQUENCE TRIMMING #
# In this section we establish reference sequences that will be used for 
# appropriate trimming of the alignment to a standard sequence length of 620 bp.
# That length was chosen as many sequences are slightly shorter than the 
# standard barcode region. Moreover, base calling errors are more frequent 
# towards the ends of the sequences. These reference sequences are manually 
# curated from the literature, thereby ensuring high sequence quality. They 
# are also specific to the taxon being analyzed.
# Currently, a standard length (658 bp) COI-5P sequence from Perca flavescens 
# (yellow perch) is being used to trim Actinopterygii barcode sequences.

# refSeqTrim function.
refSeqTrim <- function(data) {
  
  dfRefSeqs <- data.frame(taxa = c("Actinopterygii"),
                          nucleotides = c("CACCCTTTATCTAGTATTTGGTGCTTGAGCCGGAATAGTGGGCACTGCCCTAAGCCTGCTTATCCGAGCAGAACTAAGCCAGCCCGGCGCTCTCCTAGGAGACGACCAGATTTATAACGTAATTGTTACAGCACATGCCTTCGTAATAATTTTCTTTATAGTAATACCAATTATGATTGGGGGCTTTGGAAACTGACTAATTCCACTTATGATCGGTGCCCCTGACATAGCTTTCCCTCGAATAAATAATATGAGCTTTTGGCTCCTGCCTCCTTCTTTCCTTCTCCTCCTTGCTTCCTCAGGAGTTGAAGCCGGAGCTGGTACCGGATGAACTGTTTATCCCCCTCTTGCTGGGAACTTAGCACATGCTGGAGCATCTGTTGATTTAACCATTTTCTCTTTACACTTAGCAGGGGTTTCCTCAATTCTAGGTGCTATTAATTTTATTACAACCATCATTAATATAAAACCCCCTGCCATTTCCCAATATCAAACTCCCTTGTTCGTATGGGCTGTATTAATTACCGCCGTTCTTCTCCTTCTTTCACTACCTGTTCTTGCCGCTGGCATTACAATGCTTCTTACAGACCGAAATTTGAACACCACTTTCTTCGATCCTGCAGGAGGGGGTGATCCCATCCTTTACCAACACTTATTC"))
  colnames(dfRefSeqs)[2] <- "nucleotides"
  dfRefSeqs <- setDT(dfRefSeqs)
  dfRefSeqs[, nucleotides := as.character(nucleotides)]

  # Symmetrical trimming of the references to a standard 620 bp from 658 bp.
  # A different trimming length could be used, depending upon the distribution of 
  # sequence lengths in a particular taxon.
  dfRefSeqs[, nucleotides := substr(nucleotides, 20, nchar(nucleotides) - 19)]

  # Check sequence length.
  dfRefSeqs[, seq_length := nchar(nucleotides)]
  
  # Remove gaps prior to aligning (when using DECIPHER)
  data$nucleotides <- data[, gsub("-", "", nucleotides)]

  # We must ensure that the sequences are of the chr type when all of the sequences 
  # PLUS the reference sequence(s) are combined into a vector. The reference 
  # sequence is added as the first sequence.
  alignmentSeqs <- as.character(data$nucleotides)

  # Name our sequences according to BIN URI.
  names(alignmentSeqs) <- data$bin_uri
  alignmentRef <- as.character(dfRefSeqs$nucleotides[1])

  # Name our reference sequence "REFERENCE".
  names(alignmentRef) <- "REFERENCE"

  # Append our sequences together.
  alignmentSeqsPlusRef <- append(alignmentRef, alignmentSeqs)
  rm(alignmentRef)
  rm(alignmentSeqs)

  # Convert all sequences in alignmentSeqsPlusRef to DNAStringSet format. 
  # This is the format required for the alignment.
  DNAStringSet2 <- DNAStringSet(alignmentSeqsPlusRef)
  rm(alignmentSeqsPlusRef)

  # Run a multiple sequence alignment of all sequences including the reference 
  # using MUSCLE. This could take several minutes depending on the number of 
  # sequences and computer speed.
  gc()
  alignment2 <- muscle::muscle(DNAStringSet2, diags = TRUE, gapopen = -3000)
  
  # If you want to save the alignment as a FASTA file to your current working
  # directory, uncomment the following lines. The file will be named according to 
  # the class of organisms whose barcode sequences you are you currently analyzing.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentUntrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  # Convert to DNAStringSet format.
  alignmentUntrimmed <- DNAStringSet(alignment2)
  writeXStringSet(alignmentUntrimmed, file = classFileNames[[1]], 
                  format = "fasta")

  # For trimming of the sequences, we have to determine where in the alignment the 
  # reference sequence is and determine its start and stop positions relative to 
  # the other sequences. We can then use these positions to trim the rest of the 
  # sequences in the alignment.
  refSeqPos <- which(alignment2@unmasked@ranges@NAMES == "REFERENCE")
  refSeqPos <- alignment2@unmasked[refSeqPos]

  # Find the start position by searching for the first nucleotide position of the 
  # reference sequence.
  refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
  refSeqPosStart <- as.numeric(refSeqPosStart)

  # Find last nucleotide position of the reference sequence.
  refSeqPosEnd <- nchar(dfRefSeqs$nucleotides[1]) + refSeqPosStart
  refSeqPosEnd <- as.numeric(refSeqPosEnd)

  # Then we can take a substring of the alignment by using 
  # these positions to effectively "trim" the alignment.
  alignment2Trimmed <- substr(alignment2, refSeqPosStart, refSeqPosEnd)

  # Checking alignment.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentTrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  # Convert to DNAStringSet format.
  alignmentTrimmed <- DNAStringSet(alignment2)
  writeXStringSet(alignmentTrimmed, file = classFileNames[[1]], 
                  format = "fasta", width = 1500)

  # Again, convert to DNAStringSet format.
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)

  # Remove the reference sequence from this, as we dont want it to be included in 
  # further analysis.
  refSeqRm <- which(DNAStringSet3@ranges@NAMES == "REFERENCE")
  DNAStringSet3 <- subset(DNAStringSet3[-refSeqRm])

  # Reorder dfAllSeqs according to the order of species produced by the alignment, 
  # which are now contained in the DNA_StringSet object.
  # Make a variable with the ordering.
  alignmentOrder <- DNAStringSet3@ranges@NAMES

  # Order dfAllSeqs according to this.
  data <- data[match(alignmentOrder, data$bin_uri), ]

  # Repopulate dfAllSeqs with the newly trimmed sequences instead of the raw 
  # sequences.
  trimmedSeqs <- as.character(DNAStringSet3)
  data$nucleotides <- trimmedSeqs
  # Make sure species and sequences are correctly matched up!!!
  
  return(data)
}

dfAllSeqs <- refSeqTrim(dfAllSeqs)


### ALIGNMENT QUALITY CHECKING ###
# Here, extremely gappy/ungappy sequences are removed. These sequences are assumed
# to contribute to misalignment of the sequences or may even be pseudogenes.
# Manually checking of the alignment is recommended.

# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGaps <- sapply(regmatches(dfAllSeqs$nucleotides, gregexpr("[-+]", dfAllSeqs$nucleotides)), length)
  
# Mean gap length and range.
meanGap <- mean(internalGaps)
extremeHighGap <- meanGap + 7  # Upper range.
extremeLowGap <- meanGap - 7  # Lower range.

# We then loop through each sequence to see if the number of gaps 
# deviates greatly from the mean.
# Which sequences exceed the range of meanGap +/- 7?
extremeSeqs <- foreach(i = 1:nrow(dfAllSeqs)) %do% 
  which(internalGaps[[i]] > extremeHighGap | internalGaps[[i]] < extremeLowGap)
# The "deviant" sequences will be flagged with a 1.
extremeBins <- which(extremeSeqs > 0)

# Subset out these sequences to look at them if desired.
dfExtreme <- dfAllSeqs[extremeBins, ]
# If you decide to remove all from your data:
dfAllSeqs <- dfAllSeqs[-extremeBins, ]
# Align and trim again without the extreme BINs/sequences.
dfAllSeqs <- refSeqTrim(dfAllSeqs)  # Check over sequences/alignment.



### SECTION 4: WORKING PHYLOGENY ###
# This is the section where a working phylogeny is constructed based on the
# centroid sequence data. This is for use in phylogenetic regression in the 
# downstream analysis. The tree built here will also be compared to one built
# using external programs.
# Adapting code provided by Paradis (2012) and Schliep (2016).

# First, assign the sequences to a vector.
alignmentSeqs <- as.character(dfAllSeqs$nucleotides)
# Name our sequences according to the BIN URIs.
names(alignmentSeqs) <- dfAllSeqs$bin_uri
# Convert sequences to DNAStringSet format.
DNAStringSet4 <- DNAStringSet(alignmentSeqs)
# Convert to a DNAbin object so that a distance matrix may be calculated.
dnaBIN <- as.DNAbin(DNAStringSet4)
# Create a phyDat object for use with pml() functions (maximum likelihood methods).
alignmentPhyDat <- as.phyDat(dnaBIN)

# Create a distance matrix to build a starting neighbour joining tree.
# Cite: Schliep (2016).
distK80 <- dist.dna(dnaBIN, pairwise.deletion = TRUE)
hist(distK80) # Distribution of pairwise distances.
# Explore use of different substitution models.
#distJC69 <- dist.dna(dnaBIN, model = "JC69", pairwise.deletion = TRUE)
#distF84 <- dist.dna(dnaBIN, model = "F84", pairwise.deletion = TRUE)
#distTN93 <- dist.dna(dnaBIN, model = "TN93", pairwise.deletion = TRUE)
#distGG95 <- dist.dna(dnaBIN, model = "GG95", pairwise.deletion = TRUE)
#distRaw <- dist.dna(dnaBIN, model = "raw", pairwise.deletion = TRUE)

# Compare these matrices by looking at their correlations.
# Which are different?
#round(cor(cbind(distK80, distJC69, distF84, distTN93, distGG95, distRaw)), 3)

# Compare the distances estimated from the different methods.
# You must take a subset of the data as they are unable to plot otherwise.
#layout(matrix(1:2, 1))
#plot(distJC69, distRaw); abline(b = 1, a = 0)  # Effects of multiple substituion.
#plot(distK80, distJC69); abline(b = 1, a = 0)  # Ts/Tv does not seem to greatly effect.
#plot(distK80, distF84); abline(b = 1, a = 0)  # Base frequencies no effect.
#plot(distK80, distTN93); abline(b = 1, a = 0)  # Diff rates for Ts/Tv/base freqs different no effect.
#plot(distK80, distGG95); abline(b = 1, a = 0)  # Mess!

# Create some neighbour joining trees.
njK80 <- nj(distK80)
njGG95 <- nj(distGG95)

# Check the topological distance between the trees.
dist.topo(njK80, njGG95)

# Root the tree.
f <- function(xx) root(nj(dist.dna(xx, p = TRUE)), "AAA9092")
tr <- f(dnaBIN)

# Call boot.phylo to perform bootstrap. This may take a while.
#nj.boot <- boot.phylo(tr, as.matrix(dnaBIN), f, 100, rooted = TRUE)

# Plot the NJ tree with bootstrap values.
#plot(njK80, no.margin = TRUE)
#nodelabels(round(nj.boot / 100, 2), bg = "white")
#add.scale.bar(length = 0.01)


### MAXIMUM LIKELIHOOD TREE ###

# First, let's perform a modelTest to estimate the substitution model that best
# fits our data.
mt <- modelTest(alignmentPhyDat)  # Takes a while.
print(mt)  # Print the results.
mt[order(mt$BIC), ]  # Order by BIC value.

# Build a ML Tree using phangorn (pml() functions).
# Note: These may take a while.
# Using the NJ tree as a starting tree.
treeML <- pml(njK80, alignmentPhyDat)  
# Optimize the tree using the JC model.
treeJC <- optim.pml(treeML, optNni = TRUE)
# Optimize the tree using the HKY model.
treeHKY <- optim.pml(treeML, model = "HKY", optNni = TRUE)
# Optimize the tree using the HKY+G+I model.
treeHKY.G.I <- optim.pml(update(treeML, k = 4), model = "HKY", optNni = TRUE, 
                         optGamma = TRUE, optInv = TRUE)
# Optimize the tree using the GTR model.
treeGTR <- optim.pml(treeML, model = "GTR", optNni = TRUE)
# Optimize the tree using the GTR+G+I model.
treeGTR.G.I <- optim.pml(update(treeML, k = 4), model = "GTR", optNni = TRUE, 
                         optGamma = TRUE, optInv = TRUE)

# Comparing the models.
anova(treeJC, treeHKY)
anova(treeHKY, treeHKY.G.I)
anova(treeHKY.G.I, treeGTR)
anova(treeHKY.G.I, treeGTR.G.I)  # Interesting.
AIC(treeHKY.G.I, treeGTR.G.I) 
BIC(treeHKY.G.I)
BIC(treeGTR.G.I)

# You are able to optimize other parameters using the pml() function.
# Let's optimize base frequency and rate matrix.
treeGTR.G.I.2 <- optim.pml(update(treeML, k = 4), model = "GTR", optNni = TRUE, 
                           optGamma = TRUE, optInv = TRUE, optBf = TRUE, optQ = TRUE)
anova(treeGTR.G.I, treeGTR.G.I.2)  # No significant difference compared to the first model.

# SH test.
SH.test(treeGTR.G.I, treeHKY.G.I, treeJC)
SH.test(treeGTR.G.I, treeHKY.G.I)

# Partitioned models.
# Estimating phylogenies based on the different codon positions.
alignmentPhyDat
# Create weight table based on codons.
w <- table(attr(alignmentPhyDat, "index"), rep(1:3, length.out = 620))
dim(w)
# Optimize the rate locally (for each codon position).
part1 <- pmlPart(~ rate, treeGTR.G.I, weight = w, control = ctr)  # loglikelihood: -305226.2  AIC:  618698.4  BIC:  636962.1  
# Optimize rate and shape locally.
part2 <- pmlPart(~ rate + shape, treeGTR.G.I, weight = w, control = ctr)  # loglikelihood: -302897.2  AIC:  614044.3  BIC:  632316.9  
# Optimize rate, shape and nni locally. Takes a long time - haven't finished it yet!
part3 <- pmlPart(~ rate + shape + nni, treeGTR.G.I, weight = w, control = ctr) 
# Extract trees from the pmlPart object.
trees <- pmlPart2multiPhylo(part3)

# Mixture models.
# Higher AIC and BIC than partioned models.
mix1 <- pmlMix(~ rate, treeJC, m = 4, control = ctr)  
mix2 <- logLik(optim.pml(update(treeJC, k = 4), optGamma = TRUE, control = ctr)) 

# Read the ML tree back in so I don't have to make a new tree everytime!
treeGTR.G.I$tree <- read.tree(file = "MLtree.tre")

################################################################################
# TRAIT: NUMBER OF NODES
# Let's first determine the number of nodes for each species.
phy4ML <- as(treeML, "phylo4")
#plot(phy4ML, show.node = TRUE)
root <- rootNode(phy4ML)
nodeList <- lapply(1:nTips(phy4ML), function(i) .tipToRoot(phy4ML, i, root))  # A bit slow.
numberOfNodes <- sapply(1:nTips(phy4ML), function(i) length(nodeList[[i]]))  # Check if these are accurate.
names(numberOfNodes) <- tipLabels(phy4ML)
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(numberOfNodes)
dfNumberOfNodes$bin_uri <- row.names(dfNumberOfNodes)
################################################################################

# Merge with the node_number df.
dfAllSeqs <- merge(dfAllSeqs, dfNumberOfNodes, by = "bin_uri")

# Let's calculate the branch lengths of the tree now.
# Still need to plot these on the tree and determine accuracy.
branchLengths <- diag(vcv.phylo(treeML))  # Root to tip distances.

# Add the branch lengths to the dfRecodedPreCent dataframe.
dfBranchLengths <- data.frame(branchLengths)
colnames(dfBranchLengths)[1] <- "branchLength"
dfBranchLengths$bin_uri <- names(branchLengths)

hist(dfBranchLengths$branchLength)

# Merge dfBranchLengths and dfRecodedPreCent.
dfRegression <- merge(dfBranchLengths, dfAllSeqs, by = "bin_uri")
hist(dfRegression$branchLength)
rm(branchLengths)
rm(dfBranchLengths)

# Dataframe reorganization
dfRegression <- dfRegression[, -c(3, 5)]
# Write the data to external file.
write.csv(dfRegression, file = "RegressionData.csv")

# Read the data back in so I don't have to run through the whole pipeline 
# everytime.
dfRegression <- read.csv(file="RegressionData.csv", header=TRUE, sep=",")



### SECTION 5: STATISTICS ###
# Order the data according to the tree.
dfRegression <- dfRegression[match(treeML$tip.label, dfRegression$bin_uri), ]

# Continuous variables.
# Median latitude.
hist(dfRegression$median_lat)
hist(scale(dfRegression$median_lat))
dfRegression$median_lat <- scale(dfRegression$median_lat)

# Lat range.
hist(dfRegression$lat_range) # Don't think I can fix this skew?
dfRegression$lat_range <- NULL  # Removing for now.

# Average depth.
hist(dfRegression$avg_depth) # This one is....questionable.
dfRegression$avg_depth <- NULL  # Removing for now.

# Total length.
hist(dfRegression$total_length) # Leave it as is?

# Standard length.
hist(dfRegression$standard_length) # Leave it as is?

# Fork length.
hist(dfRegression$fork_length) # Leave it as is?

# Head length.
hist(dfRegression$head_length) # Leave it as is?

# Body depth.
hist(dfRegression$body_depth) # Leave it as is?

# Aspect ratio.
hist(dfRegression$aspect_ratio)  # Skewed!
hist(log(dfRegression$aspect_ratio))
dfRegression$aspect_ratio <- scale(dfRegression$aspect_ratio)

# Food troph.
#plot(dfRegression$FoodTroph) 

# Food troph SE.
#hist(dfRegression$FoodSeTroph)  # Leave it as is?

# Categorical traits.
# Inspect the frequency distribution of factor type predictors.

# Endemicity.
table(dfRegression$endemic)  # 1 is rare compared to 0.
dfRegression$endemic <- NULL  # Removing for now.

# Introduced
table(dfRegression$introduced)  # 1 is rare compared to 0.

# Ecosystem type.
#table(dfRegression$ecosystem_type)  # 1 is very dominant.
#dfRegression$ecosystemType <- NULL  # Removing for now.

# Salinity tolerance.
table(dfRegression$salinity_tolerance)

# Salinity tolerance.
table(dfRegression$env_salinity_level)  # 1 too rare & NO 2s!
dfRegression$env_salinity_level <- NULL  # Removing for now.

# Body Shape I.
table(dfRegression$BodyShapeI)  # Should remove eel-like species as they were be influential.
# Omit them from the data:
dfRegression <- subset(x = dfRegression, BodyShapeI!="eel-like")
nrow(dfRegression)
dfRegression <- droplevels(x = dfRegression)
treeML <- drop.tip(phy = treeML, 
                   tip = treeML$tip.label[!treeML$tip.label%in%as.character(dfRegression$bin_uri)], rooted = T)
str(treeML)

# Operculum present.
table(dfRegression$OperculumPresent) # Good!
dfRegression$OperculumPresent <- as.factor(dfRegression$OperculumPresent)  # Make sure it is a factor.

# Neritic.
table(dfRegression$Neritic) # Is this OK?
dfRegression$Neritic <- as.factor(dfRegression$Neritic)  # Make sure it is a factor.

# Intertidal.
table(dfRegression$Intertidal) # Is this OK?
dfRegression$Intertidal <- as.factor(dfRegression$Intertidal)  # Make sure it is a factor.

# Oceanic.
table(dfRegression$Oceanic) # Is this OK?
dfRegression$Oceanic <- as.factor(dfRegression$Oceanic)  # Make sure it is a factor.

# Mesopelagic
table(dfRegression$Mesopelagic) # Is this OK?
dfRegression$Mesopelagic <- as.factor(dfRegression$Mesopelagic)  # Make sure it is a factor.

# Estuaries.
table(dfRegression$Estuaries) # Is this OK?
dfRegression$Estuaries <- as.factor(dfRegression$Estuaries)  # Make sure it is a factor.

# Mangroves.
table(dfRegression$Mangroves) # Is this OK?
dfRegression$Mangroves <- as.factor(dfRegression$Mangroves)  # Make sure it is a factor.

# Stream.
table(dfRegression$Stream) # Is this OK?
dfRegression$Stream <- as.factor(dfRegression$Stream)  # Make sure it is a factor.

# Lakes.
table(dfRegression$Lakes) # Is this OK?
dfRegression$Lakes <- NULL  # Removing for now.

# Herbivory2.
#table(dfRegression$Herbivory2)

# Benthic.
table(dfRegression$Benthic) # Is this OK?
dfRegression$Benthic <- as.factor(dfRegression$Benthic)  # Make sure it is a factor.

# Soft bottom.
table(dfRegression$SoftBottom) # Is this OK?
dfRegression$SoftBottom <- as.factor(dfRegression$SoftBottom)  # Make sure it is a factor.

# Sand.
table(dfRegression$Sand) # Is this OK?
dfRegression$Sand <- as.factor(dfRegression$Sand)  # Make sure it is a factor.

# Mud.
table(dfRegression$Mud) # Is this OK?
dfRegression$Mud <- as.factor(dfRegression$Mud)  # Make sure it is a factor.

# Hard bottom.
table(dfRegression$HardBottom)
dfRegression$HardBottom <- as.factor(dfRegression$HardBottom)  # Make sure it is a factor.

# Rocky.
table(dfRegression$Rocky)
dfRegression$Rocky <- as.factor(dfRegression$Rocky)  # Make sure it is a factor.

# Seagrass beds.
table(dfRegression$SeaGrassBeds)
dfRegression$SeaGrassBeds <- as.factor(dfRegression$SeaGrassBeds)  # Make sure it is a factor.

# Coral reefs.
table(dfRegression$CoralReefs)
dfRegression$CoralReefs <- as.factor(dfRegression$CoralReefs)  # Make sure it is a factor.

# Response variable (branch length).
hist(dfRegression$branchLength) # Looks good.
# Some summary statistics.
mean(dfRegression$branchLength)
summary(dfRegression$branchLength)

# Checking for collinearity.
lm.res <- lm(branchLength ~ initial_bin_size + filtered_bin_size + median_lat +
             introduced + salinity_tolerance + OperculumPresent + BodyShapeI +
             head_length + standard_length + total_length + body_depth + 
             aspect_ratio + Neritic + Intertidal + Oceanic + Mesopelagic +
             Estuaries + Mangroves + Stream + Benthic + SoftBottom + Sand + 
             Mud + HardBottom + Rocky + SeaGrassBeds + CoralReefs +
             fork_length + numberOfNodes, data = dfRegression)
vif(mod = lm.res)

# Removing variables with high ( > 5) VIF values. 
lm.res <- lm(branchLength ~ initial_bin_size + filtered_bin_size + median_lat +
             introduced + salinity_tolerance + OperculumPresent + BodyShapeI +
             head_length + total_length + body_depth + aspect_ratio + Neritic + 
             Intertidal + Oceanic + Mesopelagic + Estuaries + Mangroves + 
             Stream + Benthic + SoftBottom + Sand + Mud + HardBottom + Rocky + 
             SeaGrassBeds + CoralReefs + numberOfNodes, data = dfRegression)
vif(mod = lm.res)
# Remove high VIF variables from dfRegression.
dfRegression$standard_length <- NULL
dfRegression$fork_length <- NULL

# PGLS using caper.
cDat <- comparative.data(phy = treeML, data = dfRegression, names.col = bin_uri, 
                         vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
full <- pgls(branchLength ~ initial_bin_size + filtered_bin_size + median_lat +
             introduced + salinity_tolerance + OperculumPresent + BodyShapeI +
             head_length + total_length + body_depth + aspect_ratio + Neritic + 
             Intertidal + Oceanic + Mesopelagic + Estuaries + Mangroves + 
             Stream + Benthic + SoftBottom + Sand + Mud + HardBottom + Rocky + 
             SeaGrassBeds + CoralReefs + numberOfNodes, data = cDat, 
             lambda = "ML")

# Inspect phylogenetic residuals for normality:
# Inspect a histogram:
par(mar = rep(2, 4))
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)

# Inspect phylogenetic residuals for homogeneity.
plot(x=fitted(full), y=full$phyres, pch=19)

# Some more plots.
source("pgls_check_fncs.r")
diagnostics.plot(full)

# Removal of outliers with studentized residuals >?3.
# Adapted from: http://www.anthrotree.info/wiki/pages/w5j9m8/7.5.html.
res <- residuals(full, phylo = TRUE) # Extracts phylogenetic residuals from the pgls model.
res <- res/sqrt(var(res))[1] # Standardizes residuals by sqrt of their variance.
rownames(res) <- rownames(full$residuals) # Matches the residuals up with the row name.
rownames(res)[(abs(res) > 3)] # Gives the names of the outliers.

# Remove if any outliers.
cDat <- cDat[-which(abs(res) > 3), ]
dfRegression <- dfRegression[-which(abs(res) > 3), ]
treeML <- drop.tip(phy = treeML, 
                   tip = treeML$tip.label[!treeML$tip.label%in%as.character(dfRegression$bin_uri)], rooted = T)
str(treeML)
full2 <- pgls(branchLength ~ initial_bin_size + filtered_bin_size + median_lat +
                introduced + salinity_tolerance + OperculumPresent + BodyShapeI +
                head_length + total_length + body_depth + aspect_ratio + Neritic + 
                Intertidal + Oceanic + Mesopelagic + Estuaries + Mangroves + 
                Stream + Benthic + SoftBottom + Sand + Mud + HardBottom + Rocky + 
                SeaGrassBeds + CoralReefs + numberOfNodes, data = cDat, 
              lambda="ML")

# Look at the plots again.
par(mar = rep(2, 4))
hist(full2$phyres)
qqnorm(full2$phyres)
qqline(full2$phyres)
# Inspect phylogenetic residuals for homogeneity.
plot(x = fitted(full2), y = full2$phyres, pch = 19)
diagnostics.plot(full2)

# FEATURE SELECTION!
# Feature selection using Boruta.
# Adapted from: https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/.
boruta <- Boruta(branchLength ~ initial_bin_size + filtered_bin_size + median_lat +
                   introduced + salinity_tolerance + OperculumPresent + BodyShapeI +
                   head_length + total_length + body_depth + aspect_ratio + Neritic + 
                   Intertidal + Oceanic + Mesopelagic + Estuaries + Mangroves + 
                   Stream + Benthic + SoftBottom + Sand + Mud + HardBottom + Rocky + 
                   SeaGrassBeds + CoralReefs + numberOfNodes, data = dfRegression, 
                 doTrace = 2)
plot(boruta, xlab = "", xaxt = "n")
lz <- lapply(1:ncol(boruta$ImpHistory), function(i) 
  boruta$ImpHistory[is.finite(boruta$ImpHistory[, i]), i])
names(lz) <- colnames(boruta$ImpHistory)
Labels <- sort(sapply(lz, median))
axis(side = 1, las=2, labels = names(Labels), at = 1:ncol(boruta$ImpHistory), 
     cex.axis = 0.7)

# Decide on tentative attributes.
final.boruta <- TentativeRoughFix(boruta)
print(final.boruta)

# Confirmed attributes.
getSelectedAttributes(final.boruta, withTentative = F)
boruta.df <- attStats(final.boruta)
print(boruta.df)


# Recursive feature selection.
# VARIABLE SELECTION!
# Recursive feature selection.
dfRegression$X <- NULL
predictors <- dfRegression[, c(7:33)]
outcomes <- dfRegression$branchLength
control <- rfeControl(functions = rfFuncs, method = "cv", number = 10)
rfe.train <- rfe(predictors, outcomes, sizes = 1:26, rfeControl = control)
rfe.train
plot(rfe.train, type = c("g", "o"), cex = 1.0)
predictors(rfe.train)
rfe.train$optsize


# Model selection/averaging using the most important predictors.
# Adapted from: http://www.mpcm-evolution.org/practice/online-practical-material-chapter-12/chapter-12-1-evolutionary-models-different-combinations-predictors

# Create a variable containing all possible predictors in order of importance.
# Note: numberOfNodes is a control variable.
pred.vars <- c("numberOfNodes", "aspectRatio", "bodyDepth", "headLength", 
               "FoodTroph", "forkLength", "median_lat", "BodyShapeI", 
               "totalLength", "OperculumPresent", "CoralReefs", "Rocky",
               "Oceanic", "Sand", "HardBottom", "Mangroves", "Benthic", 
               "Herbivory2", "Neritic", "SoftBottom", "SeaGrassBeds", 
               "Intertidal", "Stream", "Estuaries")
len <- length(pred.vars)
tempLen <- len
models <- NULL

# All possible models.
for (i in 1:len) {
  tempModels <- NULL
  for (n in 1:tempLen) {
    if (n == 1) {
      tempModels <- paste(tempModels, pred.vars[n], sep = "1 + ")
    } else {
      tempModels <- paste(tempModels, pred.vars[n], sep = " + ")
    }
  }
  models <- append(models, tempModels)
  tempLen <- tempLen - 1
}
# Clean up.
rm(len)
rm(tempLen)
# Complete the models.
models <- paste("branchLength", models, sep = " ~ ")

# AIC of models.
all.aic <- rep(NA, length(models))
# Estimated lambdas
all.lambda <- rep(NA, length(models))

# Run PGLS for each model.
for (k in 1:length(models)) {
  res <- try(pgls(as.formula(models[k]), data = cDat, lambda = "ML"))
  if (class(res) != "try-error") {
    all.aic[k] <- AIC(res)
    all.lambda[k] <- summary(res)$param[2]
    #xx <- coefficients(res)
    #coefs[res, match(names(xx), colnames(dfRegression))] <- xx
  }
}

# Which model to pick?
min(all.aic)
models[which(all.aic == min(all.aic))]
coefs[which(all.aic == min(all.aic)), ]
all.lambda[which(all.aic == min(all.aic))]

# PGLS using ape and nlme packages.
# Cite: Modern Phylogenetic Comparative Methods and Their Application in 
# Evolutionary Biology. Chapters 5 & 6.
# I want to make sure dfRegression matches the order of treeNJ for downstream
# pgls analysis.
dfRegression <- dfRegression$bin_uri[treeML$tip.label, ]

# PGLS model assuming Brownian Motion.
fitBM <- gls(branchLength ~ median_lat + numberOfNodes, 
             corBrownian(phy=treeML), data = dfRegression)

# Check the residuals.
plot(fitBM, resid(., type = "n") ~ fitted(.), 
     main = "Normalized Residuals vs. Fitted Values", abline = c(0, 0))
res <- resid(fitBM, type="n")
qqnorm(res)
qqline(res)

# Check for outliers and remove them.
res[which.max(res)]
#dat3 <- dat[-which(rownames(dat) == "Pg"),]
#tree3 <- drop.tip(tree, "Pg")
#fit3 <- gls(matur.L ~ age.mat, correlation=corBrownian(phy=tree3), data=dat3)

# Summary information.
summary(fitBM)

# Experiment with incorporating Pagel's lamdbda.
fitPagel <- gls(branchLength ~ median_lat + numberOfNodes, 
                correlation = corPagel(value = 0.8, phy = treeML), 
                data = dfRegression)
intervals(fitPagel, which = "var-cov")

# Testing hypotheses of lambda = 0 vs. lambda = 0 using a likelihood ratio test.
# Independence - lambda = 0.
fitPagel0 <- gls(branchLength ~ median_lat + numberOfNodes, 
                 correlation = corPagel(value = 0, phy = treeML, fixed = TRUE), 
                 data = dfRegression)
# Brownian motion - lamdba = 1.
fitPagel1 <- gls(branchLength ~ median_lat + numberOfNodes, 
                 correlation = corPagel(value = 1, phy = treeML, fixed = TRUE), 
                 data = dfRegression) 
anova(fitPagel, fitPagel0)
anova(fitPagel, fitPagel1)

# Estimating lambda.
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda, function(lambda) 
  logLik(gls(branchLength ~ median_lat + numberOfNodes, 
             correlation = corPagel(value = lambda, phy = treeML, fixed = TRUE),
             data = dfRegression)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
abline(v = fitPagel$modelStruct, col = "red")

# Plot branch length vs. median lat.
with(dfRegression, plot(branchLength ~ median_lat, xlab = "Median Latitude",
                        ylab = "Branch Length", 
                        main = "Branch length vs. median latitude"))
abline(fitPagel)

# Estimate lambda for each trait.
# Set up matrix to hold the results.
mat <- matrix(NA, ncol = 3, nrow = (dim(dat3)[2]), 
              dimnames = list(variable = names(dfRegression[blahtoblah]), 
                              c("Lower" , "Estimate", "Upper")))

# Loop through each variable.
for (i in 1:(dim(dat)[dfRegression])) {  
  form <- formula(paste(names(dat)[i], " ~ 1")) 
  this.fit <- gls(form, data = dfRegression, 
                  correlation = corPagel(0.1, phy = treeML),
                  control = glsControl(opt = "optim"))  # Fit the model.
  ints <- try(intervals(this.fit)$corStruct)  # Get lambda and CI.
# If there is an error, just get the point estimate and set the confidence 
  # limits to NA.
  if (inherits(ints, "try-error")) ints <- c(NA, coef(this.fit$modelStruct), NA) 
  mat[i, ] <- ints # save the results in the matrix
}   # End loop.


# Look at plots of variables.
#pairs(branchLength~latitude + avg_depth + endemicity + introducedNon + 
#ecosystemType + salinity, data = dfRegression,  
#main = "Simple Scatterplot Matrix")
plot(branchLength ~ bodyDepth, data = dfRegression)
reg <- gls(branchLength ~ bodyDepth, data = dfRegression)
abline(reg, col = "red")

# Phylogenetic generalized least squares analysis.
# Utilizing the pgls() function in caper.
# "Rows of data need to be matched carefully to the tips to ensure that the 
# phylogenetic structure in the variables is represented correctly"
# - comparative.data function.
# vcv = TRUE to include covariance matrix.

### Using different branch length transformations of the covariance matrix. ###
# Using lambda.
# comparative.data() make take some time.
dfRegressionOpt <- dfRegression[, c(1:2, 5, 12:13, 16:18, 20:21, 23, 25, 29, 31, 
                                35:36, 38)]
cDat <- comparative.data(data = dfRegressionOpt, phy = treeML, 
                         names.col = "bin_uri", vcv = TRUE)

# Stepwise regression.
fit <- pgls(branchLength ~ ., cDat, lambda = "ML")
summary(fit)
step <- step(fit, direction="backward")
step$anova # display results
# Set lambda = "ML" for maximum likelihood estimation.
# Latitude.
modL <- pgls(branchLength ~ numberOfNodes, cDat, lambda = "ML") 
summary(modL)
anova(modL)
# Lambda profile.
#L <- pgls.profile(modL, which = c("lambda"), N = 50, param.CI = NULL) 
#pgls.confint(modL)$opt  # Optimal lambda value = 0.9794417

# Lat range.
modL <- pgls(branchLength ~ bodyDepth + numberOfNodes + aspect Ratio + forkLength 
             + FoodTroph + totalLength + median_lat + HardBottom + Food, cDat, lambda = "ML") 
summary(modL)

# Std length.
modL <- pgls(branchLength ~ standardLength + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Std length.
modL <- pgls(branchLength ~ totalLength + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Body depth.
modL <- pgls(branchLength ~ bodyDepth + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# FoodTroph
modL <- pgls(branchLength ~ FoodTroph + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Oceanic
modL <- pgls(branchLength ~ Oceanic + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Lakes
modL <- pgls(branchLength ~ Lakes + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Benthic
modL <- pgls(branchLength ~ Benthic + numberOfNodes, cDat, lambda = "ML") 
summary(modL)

# Mangroves
modL <- pgls(branchLength ~ Mangroves + numberOfNodes, cDat, lambda = "ML") 
summary(modL)




# Using kappa.
# Note: vcv.dim must be = 3 for kappa branch length transformation.
#cDat2 <- comparative.data(data = dfRegression, phy = treeNJ, 
                          #names.col = "bin_uri", vcv = TRUE, vcv.dim = 3)
#modK <- pgls(branchLength ~ latitude, cDat2, kappa = "ML")
#summary(modK)
# pgls.profile() = likelihood profiles for branch length transformations.
# Kappa profile.
#K <- pgls.profile(modK, which = c("kappa"), N = 50, param.CI = NULL) 
#pgls.confint(modK, which = c('kappa'))$opt  # 0.4876088 

# Using delta.
#cDat3 <- comparative.data(data = dfRegression, phy = treeNJ, 
                          #names.col = "bin_uri", vcv = TRUE)
#modD <- pgls(log(branchLength) ~ latitude, cDat3, delta = "ML")
#summary(modD)
# Delta profile.
#D <- pgls.profile(modD, which = c("delta"), N = 50, param.CI = NULL)
#pgls.confint(modD, which = c('delta'))$opt  # 0.3491073

# A potential means for model selection?
#AIC(modL, modK, modD)


### Models using a different number of predictors. ###

# This is a function that allows for use of discrete predictors with more than 
# 2 levels in pgls.
anova.pgls.fixed <- function (object) 
{
  data <- object$data
  tlabels <- attr(terms(object$formula), "term.labels")
  k <- object$k
  n <- object$n
  NR <- length(tlabels) + 1
  rss <- resdf <- rep(NA, NR)
  rss[1] <- object$NSSQ
  resdf[1] <- n - 1
  lm <- object$param["lambda"]
  dl <- object$param["delta"]
  kp <- object$param["kappa"]
  for (i in 1:length(tlabels)) {
    fmla <- as.formula(paste(object$namey, " ~ ", paste(tlabels[1:i], 
                                                        collapse = "+")))
    plm <- pgls(fmla, data, lambda = lm, delta = dl, kappa = kp)
    rss[i + 1] <- plm$RSSQ
    resdf[i + 1] <- (n - 1) - plm$k + 1
  }
  ss <- c(abs(diff(rss)), object$RSSQ)
  df <- c(abs(diff(resdf)), n - k)
  ms <- ss/df
  fval <- ms/ms[NR]
  P <- pf(fval, df, df[NR], lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f = fval, P)
  table[length(P), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Residuals"), c("Df", "Sum Sq", "Mean Sq", 
                                                     "F value", "Pr( > F)"))
  structure(table, heading = c("Analysis of Variance Table", sprintf("Sequential SS for pgls: lambda = %0.2f, delta = %0.2f, kappa = %0.2f\n", lm, dl, kp), paste("Response:", deparse(formula(object)[[2L]]))), class = c("anova", "data.frame"))
}


# 1 predictor.
mod1 <- pgls(log(branchLength) ~ sampleSize, cDat, lambda = 'ML')
summary(mod1)

# 2 predictors.
#mod2 <- pgls(log(branchLength) ~ aspectRatio + numberOfNodes, cDat, lambda = 'ML') 
#summary(mod2)  # envSalinityLevel has 3 levels.

# 3 predictors.
#mod3 <- pgls(log(branchLength) ~ latitude + envSalinityLevel + avg_depth, cDat, lambda = 'ML') 
#summary(mod3)

# 4 predictors.
#mod4 <- pgls(log(branchLength) ~ latitude + envSalinityLevel + lat_range + numberOfNodes, cDat,
 #            lambda = 'ML') 
#anova.pgls.fixed(mod2)

# 5 predictors.
#mod5 <- pgls(log(branchLength) ~ latitude + envSalinityLevel + endemicity + introducedAnywhere + 
 #            lat_range, cDat, lambda = 'ML') 
#anova.pgls.fixed(mod5)

# Compare different models (model selection).
#anova(mod1, mod2, mod3, mod4, mod5)
#AIC(mod1, mod2, mod4, mod5)

# Stop the clock!
print(proc.time() - ptm)
print("DONE")
