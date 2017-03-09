###################
# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
# correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
# determination/reference sequence trimming (lines TBD).
# source() and library() statements***
# Function definitions***
# Executed statements, if applicable (e.g., print, plot)***
# Last version saved: Mar 8th 2017 (R_Pipeline_PHYLO_V7_Feb25th.R)

### PACKAGES REQUIRED ###
# For variable selection purposes.
#install.packages("caret")
library(caret)

# For looping purposes.
library(foreach)

# For efficiency purposes.
#install.packages("data.table")
library(data.table)

# For multiple sequence alignments:
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
#biocLite("DECIPHER")
library(Biostrings)
library(DECIPHER)
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
#biocLite("muscle")
library(muscle)

# For BOLD and FishBase data:
#install.packages("bold")
library(bold)
#install.packages("rfishbase")
library(rfishbase)

# For stats
#install.packages("plyr")
library(plyr)

# For testing
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

# Missing data imputation:
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
# In this section,  the Initial TSV of the desired taxon is parsed from the BOLD 
# database API. First, we download the TSV and convert it into a datatable. 
# The URL below is what is modified by the user and will determine the taxon, 
# geographic region, etc. 
# Example: taxon=Aves$geo=all

# Note: read_tsv has been modified to select only certain columns to save on 
# downloading time.
dfInitial <- bold_seqspec(taxon = "Actinopterygii", 
                          geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                           "species_name", "lat", "nucleotides", 
                                           "markercode")] 

# If you want to run pre downloaded BOLD TSVs to avoid downloading of the same 
# tsv multiple times, this will let you choose a path to that TSV and parse:
#tsvParseDoc <- file.choose()
#dfSequences <- read_tsv(tsvParseDoc)[ ,
 #                         c("recordID", "bin_uri", "order_name","species_name", 
  #                           "lat", "nucleotides", "markercode")]
#rm(tsvParseDoc)

# Convert to data.table for efficiency purposes.
dfSequences <- setDT(dfInitial)
rm(dfInitial)

# Download outgroup species data from BOLD.
# I am currently using the great white shark BIN AAA9092 as an outgroup. This 
# choice of outgroup species may be altered in the near future.
dfOutgroup <- bold_seqspec(bin = "BOLD:AAA9092", 
                           geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                            "species_name", "lat", "nucleotides", 
                                            "markercode")] 
dfOutgroup <- setDT(dfOutgroup) # Again, convert to data.table.

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

################################################################################
### TRAIT: MEDIAN LATITUDE/LATITUDINAL RANGE ###
# Currently, latitude is the only trait whose information is being taken 
# from BOLD. The rest of the data will be obtained from FishBase.
# Gather the median latitude, minimum latitude, maximum latitude and 
# latitudinal range for each BIN.

### FILTER 3 ###
# Filtering for presence of a latitude value.
containLat <- dfFiltered[, grep("[0-9]", lat)]
dfFiltered <- dfFiltered[containLat, ]
rm(containLat)
# Convert the latitude (lat) column to number instead of character type. 
dfFiltered[, lat_num := as.numeric(lat)]
# Conversion to absolute values before median latitude values are calculated.
dfFiltered[, abs_lat_num := abs(lat_num)]
# Determine a median latitude for each BIN using absolute values.
dfFiltered[, median_lat := median(abs_lat_num), keyby = bin_uri]
# Find the minimum latitude.
dfFiltered[, lat_min := min(lat_num), keyby = bin_uri]
# Find the maximum latitude.
dfFiltered[, lat_max := max(lat_num), keyby = bin_uri]
# Determine a latitude range for each BIN.
dfFiltered[, lat_range := abs(lat_max - lat_min), keyby = bin_uri]
################################################################################

# Datatable reorganization. Removing redundant columns and columns that are not
# used in the downstream analysis.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, species_name,
                             markercode, initial_bin_size, nucleotides, 
                             median_lat, lat_range)]

### FILTER 4 ###
# Filtering for COI-5P as these are the only markers we are looking at.
setkey(dfFiltered, markercode)
dfFiltered <- dfFiltered["COI-5P"]

### FILTER 5 ###
# N and gap content will affect the multiple sequence alignment and the 
# alignment will give warning messages, so we need to filter out sequences with 
# high N and gap content.
# Convert nucleotides to chr type.
dfFiltered[, nucleotides := as.character(nucleotides)]
# Cut off large portions of Ns and gaps at the start of a sequence. 
# First, find sequences that begin with an N or a gap.
startGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("^[-N]", dfFiltered$nucleotides)), length)
startGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
# If at least one sequence is found that begins with gap or N...
  if (startGapN[[i]] > 0) { 
    # Split the sequence up!
    split <- strsplit(dfFiltered$nucleotides[i], "^[-N]+") 
    dfFiltered$nucleotides[i] <- split[[1]][2] 
  }
rm(startGapN)
# Cut off large portions of Ns and gaps at the end of a sequence.
# Same thing for sequences that end with N or a gap.
endGapN <- sapply(regmatches(dfFiltered$nucleotides, 
                             gregexpr("[-N]$", dfFiltered$nucleotides)), length)
endGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
  if (endGapN[[i]] > 0) {
    split <- strsplit(dfFiltered$nucleotides[i], "[-N]+$")
    dfFiltered$nucleotides[i] <- split[[1]][1]
  }
rm(endGapN)
rm(split)


### FILTER 6 ###
# Remove sequences with N/gap content above a certain threshold. 
# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("[-N]", dfFiltered$nucleotides)), length)
# We then iterate over each sequence and see if the number of Ns or gaps is 
# greater than 1% (0.01) of the total sequence length.
internalGapN <- foreach(i = 1:nrow(dfFiltered)) %do% 
  which((internalGapN[[i]]/nchar(dfFiltered$nucleotides[i]) > 0.01))
# Basically "flagging" the sequences with high N/gap content.
checkGapN <- sapply(internalGapN, function (x) length(x))
checkGapN <- which(checkGapN > 0)  # Identify them.
# Remove these higher gap and N content sequences.
dfFiltered <- dfFiltered[-checkGapN, ]
rm(checkGapN)
rm(i)
rm(internalGapN)

### FILTER 7 ###
# Filter out sequences that are less than or equal to 640 bp and greater than
# or equal to 1000 bp. This is because widely differing sequence lengths can 
# interfere with the alignment.
# Get the length of the sequences without gaps.
seqLengths <- dfFiltered[, nchar(gsub("-", "", nucleotides))] 
# Which sequences are greater than 1000 bp and less than 640 bp in length?
seqLengthCheck <- which(seqLengths > 1000 | seqLengths < 640)
dfFiltered <- dfFiltered[-seqLengthCheck, ]
rm(seqLengths)
rm(seqLengthCheck)

# BIN Species Information. #
### FILTER 8 ###
# Remove rows with no species information (removes BINs with no species data).
containSpecies <- dfFiltered[, grep("[A-Z]", species_name)]
# Create a new datatable containing only sequences baring species-level 
# identification. Doing this because NA values are considered in number of 
# species.
dfSpecies <- dfFiltered[containSpecies, ]
rm(containSpecies)
# Determine the number of unique species-level identifications assigned to 
# each BIN. 
dfNumberOfSpecies <- dfSpecies[, .(number_of_species = length(unique(species_name))), keyby = bin_uri]
# Determine the most common species in each BIN (this can be used later down the
# road to deal with BINs with sequences from more than one species). Now, it's 
# just used to assign species_label to sequences that have no species 
# information (giving them a "name"!).
# Note: Example of chaining.
dfSpeciesLabel <- dfSpecies[, .N, by = .(bin_uri, species_name)][order(-N), .(species_label = species_name[1L]), keyby = bin_uri]
dfSpeciesInfo <- merge(dfNumberOfSpecies, dfSpeciesLabel)
rm(dfNumberOfSpecies)
rm(dfSpeciesLabel)
rm(dfSpecies)

# MERGING DATAFRAMES #
# Note: When using merge(), "x" represents the first dataframe specified, 
# and "y" the second dataframe. "by" specifies which column you are using to 
# merge the two dataframes.
# Merge species info with the dataframe dfFiltered to obtain all of the relevant
# taxonomic information.
setkey(dfFiltered, bin_uri)
dfFiltered <- merge(dfFiltered, dfSpeciesInfo)
rm(dfSpeciesInfo)

# "CLEAN" DATA FILTERS #
# Filtering for BINs that contain more than one sequence (excluding singletons).
#atleastTwoSeqs <- which(dfFiltered$binSize > 1)
#dfFiltered <- dfFiltered[atleastTwoSeqs, ]

### FILTER 9 ###
# Filtering for BINs that contain only ONE species name.
dfFiltered <- dfFiltered[number_of_species == 1]

################################################################################
### TRAIT: POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfFiltered[, filtered_bin_size := length(recordID), by = bin_uri]
################################################################################


################################################################################
### SECTION 2: TRAITS ###

# FISH BASE TRAITS! #
# I need to extract all of the species names that are available on FishBase.
allFish <- fishbase
allFish$fullName <- paste(allFish$Genus, allFish$Species)
fishBaseSpecies <- allFish$fullName  # 33104 species names.
rm(allFish)

# Then, I need to match the species obtained from BOLD with the species obtained 
# from FishBase.
# Make into a dataframe first so I can merge them.
dfFishBaseSpecies <- data.table(fishBaseSpecies)
colnames(dfFishBaseSpecies)[1] <- "species_name"
rm(fishBaseSpecies)

# Let's merge the BOLD names with the FishBase names.
dfBoldBase <- merge(dfFiltered, dfFishBaseSpecies, by.x = "species_label",
                    by.y = "species_name")
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
# replaced with ranked or binary data. The rankings for each trait are assigned 
# based on biological reasoning. For instance, categorical traits that are 
# assumed to beget higher rates of molecular evolution are assigned a higher 
# rank. Continuous/numerical traits do not need to be transformed (except 
# perhaps a log or ln transformation). As each categorical trait is unique with 
# regards to its ranking, the code at this point is long and repetitive. Working
# on ways to make it more efficient!

### ECOSYSTEM TRAITS ###
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
                                     average_depth = AverageDepth)]

# Mode trait(s).
# Ecosystem Type.
dfEcosystemTraits[, count := .N, by = .(species_name, ecosystem_type)]
dfEcosystemTraits[order(-count), ecosystem_type := ecosystem_type[1L], by = species_name]

# Environmental salinity level.
dfEcosystemTraits[, count := .N, by = .(species_name, salinity)]
dfEcosystemTraits[order(-count), env_salinity_level := salinity[1L], by = species_name]
dfEcosystemTraits$count <- NULL

# Median trait(s).
# Average depth.
dfEcosystemTraits[, average_depth := as.double(average_depth)]
dfEcosystemTraits[, average_depth := median(average_depth, na.rm = TRUE), keyby = species_name]


# Special cases
# INTRODUCED VS. NOT INTRODUCED
# I need a better name for this function.
DetermineIfSpeciesIntroduced <- function(x) {
  answer <- NULL
  containsIntroduced <- length(grep("introduced", x))
  if (containsIntroduced > 0) {
    answer <- "Yes" 
  } else {
    answer <- "No"
  } 
  return(answer)
}

# Apply to datatable.
dfEcosystemTraits[, introduced := DetermineIfSpeciesIntroduced(status), by = species_name]

# ENDEMICITY
DetermineIfSpeciesEndemic <- function(x) {
  answer <- NULL
  containsEndemic <- length(grep("endemic", x))
  if (containsEndemic > 0) {
    answer <- "Yes" 
  } else {
    answer <- "No"
  } 
  return(answer)
}

# Apply to datatable.
dfEcosystemTraits[, endemic := DetermineIfSpeciesEndemic(status), by = species_name]

# SALINITY TOLERANCE
DetermineSalinityTolerance <- function(x) {
  answer <- NULL
  typesOfWater <- unique(x)
  number_of_water_types <- length(unique(x))
  # If species ever present in more than 1 type of salinity level,
  if (number_of_water_types > 1) {
    # then flag it as "euryhaline".
    answer <- "euryhaline"
    # Else, flag it as "stenohaline".
  } else if (number_of_water_types == 1) {
    answer <- "stenohaline"
  } else {
    answer <- NA
  }
  return(answer)
}

# Apply to datatable.
dfEcosystemTraits[, salinity_tolerance := DetermineSalinityTolerance(salinity), 
                  by = species_name]
# Remove salinity column now as it is no longer needed.
dfEcosystemTraits$salinity <- NULL

# Recoding dfEcosystemTraits.
# Make sure all of the traits are of the factor class.
changeVars <- c(2:3, 5:8)
dfEcosystemTraits[, (changeVars) := lapply(.SD, as.factor), 
                  .SDcols = changeVars]

# TRAITs: introduced & endemic.
# Recoding them in the same way.
changeVars <- c("introduced", "endemic")
dfEcosystemTraits[, (changeVars) := lapply(.SD, function(x) revalue(x, c("No" = "0", "Yes" = "1"))), 
                  .SDcols = changeVars]  # Use revalue function to recode. Double check that this works!

# TRAIT: env_salinity_level.
# 3 levels
# Fix brackish vs. euryhaline
dfEcosystemTraits[, env_salinity_level := revalue(env_salinity_level, 
                                                  c("freshwater" = "1",
                                                    "brackish" = "2",
                                                    "saltwater" = "3"))]

# TRAIT: salinity_tolerance.
dfEcosystemTraits[, salinity_tolerance := revalue(salinity_tolerance, 
                                                  c("stenohaline" = "0",
                                                    "euryhaline" = "1"))]

# TRAIT: ecosystem_type.
# 4 levels
# EDIT: Need error/flag statements when a value isn't present
dfEcosystemTraits[, ecosystem_type := revalue(ecosystem_type, 
                                              c("Sea/Bay/Gulf" = "1",
                                                "River (basin)" = "2",
                                                "Lake" = "3",
                                                "Lagoon" = "4",
                                                "Zoogeographic realm" = NA))]

# I only want one row per species.
dfEcosystemTraits <- dfEcosystemTraits[!duplicated(dfEcosystemTraits$species_name), ] 


### MORPHOMETRICS TRAITS ###
#dfMorphometrics <- data.frame(morphometrics(speciesNames))
# Storing this as a file.
#write.csv(dfMorphometrics, file = "morphometrics_information.csv") 
# Read in the ecology information.
dfMorphometrics <- fread("morphometrics_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfMorphometricTraits <- dfMorphometrics[, .(species_name = sciname, 
                                            total_length = TL, 
                                            standard_length = SL, 
                                            fork_length = FL,
                                            head_length = HL,
                                            body_depth = BD,
                                            aspect_ratio = AspectRatio)]

# As the morphometric information is entered on an individual basis, it is 
# necessary to find the median trait for each species. The following are 
# morphometrics traits that I am including in the regression analysis.

# First, change all of the traits to numeric type.
changeVars <- c(2:6)
dfMorphometricTraits[, (changeVars) := lapply(.SD, as.double), 
                     .SDcols = changeVars]

# Median traits.
# Since I want to apply the function median() to multiple columns I can do it 
# in one statement.
# Total length.
medianVars <- c(2:7)
dfMorphometricTraits[, (medianVars) := lapply(.SD, function(x) median(x, na.rm = TRUE)),
                     by = species_name, .SDcols = medianVars]

# I only want one row per species.
dfMorphometricTraits <- dfMorphometricTraits[!duplicated(dfMorphometricTraits$species_name), ] 



### MORPHOLOGY TRAITS ###
#dfMorphology <- data.frame(morphology(speciesNames))
# Storing this as a file.
#write.csv(dfMorphology, file = "morphology_information.csv") 
# Read in the ecology information.
dfMorphology <- fread("morphology_information.csv")
# Get rid of columns I do not need for the regression analysis.
dfMorphologyTraits <- dfMorphology[, .(species_name = sciname, 
                                       body_shape_I = BodyShapeI, 
                                       body_shape_II = BodyShapeII, 
                                       operculum_present = OperculumPresent, 
                                       pos_of_mouth = PosofMouth, 
                                       type_of_scales = TypeofScales)]
# These traits are unlikely to differ between different stocks of fish so I can 
# take only a single row per species.
dfMorphologyTraits <- dfMorphologyTraits[!duplicated(dfMorphologyTraits$species_name), ] 

# Recoding dfMorphologyTraits.
# Converting to factor type as these are discrete traits.
changeVars <- c(2:6)
dfMorphologyTraits[, (changeVars) := lapply(.SD, as.factor), 
                   .SDcols = changeVars]
# Also recode operculum_present.
dfMorphologyTraits[, operculum_present := revalue(operculum_present, c("-1" = "1"))]


### ECOLOGY TRAITS ###
#dfEcology <- data.frame(ecology(speciesNames))
# Storing this as a file.
#write.csv(dfEcology, file = "ecology_information.csv") 
# Read in the ecology information.
dfEcology <- fread("ecology_information.csv")
colnames(dfEcology)[3] <- "species_name"
# Get rid of columns I do not need for the regression analysis. 
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, c(3, 7:29, 31, 33:36, 39:40, 45:57, 61, 63, 66:88,
                                 90:111, 113, 115:123)]
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
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), 
                .SDcols = integerVars]
rm(integerVars)
rm(characterVars)
rm(changeVars)


# Construction of dfTraits datatable.
# Let's merge the trait information back to dfFiltered.
# NA/NULL/blank for those species who don't have info for that particular trait.
# I only want a single row per BIN.
dfFilteredSingle <- dfFiltered[!duplicated(dfFiltered$bin_uri), ]
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_label, median_lat, 
                                         lat_range, initial_bin_size, 
                                         filtered_bin_size)]
dfTraits <- merge(dfFilteredSingle, dfEcosystemTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
rm(dfFilteredSingle)
dfTraits <- merge(dfTraits, dfMorphometricTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
dfTraits <- merge(dfTraits, dfMorphologyTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)
dfTraits <- merge(dfTraits, dfEcologyTraits, by.x = "species_label", 
                  by.y = "species_name", all.x = TRUE)


# I only want complete cases (no missing data).
library(dplyr) # 0.4
# "Most" complete dataset.
y <- sort(dfTraits[, lapply(.SD, function(x) sum(is.na(x)))]) # nr of NA in columns, increasing
setcolorder(dfTraits, names(y)) # now the columns are ordered - less NA first
dfTraits[, idx := rowSums(is.na(dfTraits))] # count nr of NA in rows
dfTraits <- dfTraits[order(idx), ] # sort by nr of NA in rows
dfTraits[, idx := NULL] # idx not needed anymore
dfTemp <- dfTraits[, 3:123]


# Write the trait data in order to completeness to file
write.csv(dfTraits, file = "TraitDataCompleteness.csv")
# now your data.table is sorted: columns with least NA to the left,  
# rows with with least NA on top
len <- ncol(dfTemp)
tempLen <- len
all.cc <- NULL
# All possible complete cases.
for (i in 1:len) {
  x <- as.data.frame(dfTemp)
  x <- x[, 1:tempLen]
  x <- complete.cases(x)
  x <- which(x == "TRUE")
  x <- length(x)
  all.cc[i] <- x
  tempLen <- tempLen - 1
}

# Decide where to cut the datatable.
names(all.cc) <- rev(colnames(dfTemp))
all.cc
which(colnames(dfTraits) == "fork_length")
dfTraits <- dfTraits[, 1:109] 
# Filter the datatable so only complete cases are kept.
dfCompleteCases <- dfTraits %>% filter(complete.cases(.))

# Removal of near zero variance traits.
dfCompleteCases <- dfCompleteCases[, -nearZeroVar(dfCompleteCases)]

# Merge back to dfFiltered to get all of the sequence information for each BIN.
dfPreCentroid <- merge(dfFiltered, dfCompleteCases, by = "bin_uri")
rm(dfFiltered)
rm(dfCompleteCases)
# Dataframe reorganization.
colnames(dfPreCentroid)[6] <- "initial_bin_size"
colnames(dfPreCentroid)[8] <- "median_lat"
colnames(dfPreCentroid)[9] <- "lat_range"
colnames(dfPreCentroid)[11] <- "species_name"
colnames(dfPreCentroid)[12] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:3, 11, 6, 12, 7:9, 18:46)]


### SECTION 3: CENTROID SEQUENCE DETERMINATION ###
# A single Barcode Index Number (BIN) can contain many record ids and sequences.
# In order to simplify the analyses, we select one sequence per BIN.
# To do this we are determining the centroid sequence for each BIN.
# Centroid Sequence: BIN sequence with minimum average pairwise distance to all 
# other sequences in a given BIN.

# Subset dataframe to find BINs with more than one sequence.
 # These are the BINs that will require centroid sequences.
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

# TRAIT: Latitudinal range.
# Since singletons are included, their range is automically 0. Thus, I need to
# account for this discrepency by assigning singleton BINs an NA value.
# First, I need a recursive vector to deal with.
#recodedList <- lapply(unique(dfPreCentroid$bin_uri), 
 #                     function(x) dfPreCentroid[dfPreCentroid$bin_uri == x, ])
#RecodeLatRange <- function(listOfTraitInfo) {
  #answer <- listOfTraitInfo$lat_range
  #if (listOfTraitInfo$sample_size == 1) {
   # answer <- NA 
  #} 
  #return(answer)
#}
#recodedLatRanges <- sapply(recodedList, function(x) RecodeLatRange(x))
#rm(recodedList)
#dfPreCentroid$lat_range <- recodedLatRanges
data <- dfAllSeqs

################################################################################
# FUNCTION IDEA
# Function for reference sequence trimming.
refSeqTrim <- function(data) {
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
  
  SGC1 <- getGeneticCode("SGC1")  # Vertebrate Mitochondrial code
  Biostrings::translate(DNAStringSet1[[1]], genetic.code=SGC1, if.fuzzy.codon = "solve")
  
  # Get the right genetic code.
  verty <- getGeneticCode("SGC1")
  
  
  # Testing 
  alignment2 <- AlignTranslation(DNAStringSet2, sense = "+", direction = "5' to 3'",
                                 readingFrame = NA,
                                 type = "DNAStringSet",
                                 geneticCode = SGC1, gapOpening = -3000, iterations = 2)
  
  
  
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

# Remove extremely gappy/ungappy sequences.
# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGaps <- sapply(regmatches(dfAllSeqs$nucleotides, gregexpr("[-+]", dfAllSeqs$nucleotides)), length)
  
# Mean gap length and range.
meanGap <- mean(internalGaps)
extremeHighGap <- meanGap + 7
extremeLowGap <- meanGap - 7

# We then loop through each sequence to see if the number of gaps 
# deviates greatly from the mean.
extremeSeqs <- foreach(i = 1:nrow(dfAllSeqs)) %do% 
  which(internalGaps[[i]] > extremeHighGap | internalGaps[[i]] < extremeLowGap)
extremeBins <- which(extremeSeqs > 0)

# Subset out these higher gap and N content sequences.
# Make sure to check the alignment for these BINs to decide whether to remove
# them or not.
dfExtreme <- dfAllSeqs[extremeBins, ]
# If you decide to remove all...
dfAllSeqs <- dfAllSeqs[-extremeBins, ]
# Align and trim again without the extreme BINs/sequences.
dfAllSeqs <- refSeqTrim(dfAllSeqs) # Check over sequences/alignment.



### SECTION 4: WORKING PHYLOGENY ###
# This is the section where a working phylogeny is constructed based on the
# centroid sequence data.
# Adapting code provided by Paradis (2012) and Schliep (2016).

# Convert sequences into a vector.
alignmentSeqs <- as.character(dfAllSeqs$nucleotides)
# Name our sequences according to BIN URIs.
names(alignmentSeqs) <- dfAllSeqs$bin_uri
# Convert sequences to DNAStringSet format.
DNAStringSet4 <- DNAStringSet(alignmentSeqs)
# Convert to DNAbin object.
dnaBIN <- as.DNAbin(DNAStringSet4)
# Create a phyDat object for use with pml() functions.
alignmentPhyDat <- as.phyDat(dnaBIN)
.
# Create a distance matrix to build a starting neighbour joining tree.
# Cite: Schliep (2016).
distK80 <- dist.dna(dnaBIN, pairwise.deletion = TRUE)
hist(distK80) # Distribution of pairwise distances.
#distJC69 <- dist.dna(dnaBINsmall, model = "JC69", pairwise.deletion = TRUE)
#distF84 <- dist.dna(dnaBINsmall, model = "F84", pairwise.deletion = TRUE)
#distTN93 <- dist.dna(dnaBINsmall, model = "TN93", pairwise.deletion = TRUE)
#distGG95 <- dist.dna(dnaBINsmall, model = "GG95", pairwise.deletion = TRUE)
#distRaw <- dist.dna(dnaBINsmall, model = "raw", pairwise.deletion = TRUE)

# Compare these matrices by looking at their correlations.
#round(cor(cbind(distK80, distJC69, distF84, distTN93, distGG95, distRaw)), 3)

# Compare the distances estimated from the different methods.
# Taking a subset as they are unable to plot otherwise.
#layout(matrix(1:2, 1))
#plot(distJC69, distRaw); abline(b = 1, a = 0)  # Effects of multiple substituion.
#plot(distK80, distJC69); abline(b = 1, a = 0)  # Ts/Tv does not seem to greatly effect.
#plot(distK80, distF84); abline(b = 1, a = 0)  # Base frequencies no effect.
#plot(distK80, distTN93); abline(b = 1, a = 0)  # Diff rates for Ts/Tv/base freqs different no effect.
#plot(distK80, distGG95); abline(b = 1, a = 0)  # Mess!

align <- read.dna("TestActinopterygii.fas", format = "fasta")
align <- as.matrix(align)

# Histogram of pairwise distances for each codon position.
# Checking for heterogeneity along the sequences.
library(lattice)
layout(matrix(1:3, 1))
y <- numeric()
for (i in 1:3) {
  s <- logical(3); s[i] <- TRUE
  y <- c(y, dist.dna(align[, s], pairwise.deletion = TRUE))
}
g <- gl(3, length(y) / 3)
histogram(~ y | g, breaks = 30)

# NJ Trees.
njK80 <- nj(distK80)
njGG95 <- nj(distGG95)

# Check the topological distance between the trees.
dist.topo(njK80, njGG95)

# Root the tree.
f <- function(xx) root(nj(dist.dna(xx, p = TRUE)), "AAA5492") # Experimenting.
tr <- f(dnaBIN)

# Call boot.phylo.
nj.boot <- boot.phylo(tr, as.matrix(dnaBIN), f, 100, rooted = TRUE)

# Plot the NJ tree with bootstrap values.
plot(njK80, no.margin = TRUE)
nodelabels(round(nj.boot / 100, 2), bg = "white")
add.scale.bar(length = 0.01)

# Model selection.
mt <- modelTest(alignmentPhyDat) # Takes a while.
print(mt) # Print the results.
mt[order(mt$BIC), ] # Order by BIC value.

# ML Tree using phangorn.
treeML <- pml(njK80, alignmentPhyDat)
treeJC <- optim.pml(treeML, optNni = TRUE)
treeHKY <- optim.pml(treeML, model = "HKY", optNni = TRUE)
treeHKY.G.I <- optim.pml(update(treeML, k = 4), model = "HKY", optNni = TRUE, 
                         optGamma = TRUE, optInv = TRUE)
treeGTR <- optim.pml(treeML, model = "GTR", optNni = TRUE)
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

# Optimize base frequency and rate matrix.
treeGTR.G.I.2 <- optim.pml(update(treeML, k = 4), model = "GTR", optNni = TRUE, 
                           optGamma = TRUE, optInv = TRUE, optBf = TRUE, optQ = TRUE)
anova(treeGTR.G.I, treeGTR.G.I.2)  # No significant difference.

# SH test.
SH.test(treeGTR.G.I, treeHKY.G.I, treeJC)
SH.test(treeGTR.G.I, treeHKY.G.I)

# Partioned models.
# Estimating phylogenies based on the different codon positions.
alignmentPhyDat
# Create weight table based on codons.
w <- table(attr(alignmentPhyDat, "index"), rep(1:3, length.out = 620))
dim(w)
# Optimize rate locally.
part1 <- pmlPart(~ rate, treeGTR.G.I, weight = w, control = ctr)  # loglikelihood: -305226.2  AIC:  618698.4  BIC:  636962.1  
# Optimize rate and shape locally.
part2 <- pmlPart(~ rate + shape, treeGTR.G.I, weight = w, control = ctr)  # loglikelihood: -302897.2  AIC:  614044.3  BIC:  632316.9  
# Optimize rate, shape and nni locally.
part3 <- pmlPart(~ rate + shape + nni, treeGTR.G.I, weight = w, control = ctr) 
# Optimize rate, shape and nni locally.
part3 <- pmlPart(~ rate + shape + nni, treeGTR.G.I, weight = w, control = ctr) 
# Extract trees from the pmlPart object.
trees <- pmlPart2multiPhylo(part3)
tree1 <- trees[[1]]
tree2 <- trees[[2]]
tree3 <- trees[[3]]

# Let's calculate the branch lengths of the codon trees.
# Still need to plot these on the tree and determine accuracy.
branchLengths1 <- diag(vcv.phylo(tree1))  # Root to tip distances.
branchLengths2 <- diag(vcv.phylo(tree2))  # Root to tip distances.
branchLengths3 <- diag(vcv.phylo(tree3))  # Root to tip distances.
all.equal(branchLengths1, branchLengths3) # Hmm they are equal.

# Mixture models.
# Higher AIC and BIC than partioned models.
mix1 <- pmlMix(~ rate, treeJC, m = 4, control = ctr)  # loglikelihood: -6164.47 AIC:  12428.94  BIC:  12604.21 
mix2 <- logLik(optim.pml(update(treeJC, k = 4), optGamma = TRUE, control = ctr))  # 'log Lik.' -6168.579 (df=48)


# Read the ML tree back in so I don't have to make a new tree everytime!
treeML <- read.tree(file="MLtree.tre")

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
library(car)
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
#install.packages("Boruta")
library(Boruta)
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
library(randomForest)
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
library(gtools)

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

plot(treeML, main = "Maximum Likelihood Phylogeny", )

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
