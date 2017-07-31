###################
# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
# correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
# determination/reference sequence trimming (lines TBD).
# Last version saved: July 19th 2017 (PHYLO_Final_July5th.R)

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
require(reshape)

# For trait data purposes:
#install.packages("adephylo")
library(adephylo)
library(phylobase)
library(car)
#install.packages("phylometrics")
library(phylometrics)

# Missing data:
library(dplyr)


##### SECTION 1: DATA PROCESSING #####
# Download sequences from BOLD using the function bold_seqspec() for sequence
# and specimen data. In addition, I am only selecting those columns needed for
# downstream analysis.
# Note" Using Perciformes as an example taxon.
dfInitial <- bold_seqspec(taxon = "Perciformes", 
                          geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                           "family_name", "genus_name", 
                                           "species_name", "lat", "nucleotides", 
                                           "markercode")]
# Convert to datatable. Datatables have useful features for data manipulation.
dfInitial <- setDT(dfInitial)

# Download outgroup species data from BOLD. These sequences may be used to root 
# phylogenetic trees (depending if the taxa are an appropriate outgroup for
# the organismal group under study).
# Note: These are just example taxa.
outgroups <- c("Acanthuriformes", "Gadiformes")
dfOutgroup <- bold_seqspec(taxon = outgroups, 
                           geo = "all")[, c("recordID", "bin_uri", "order_name", 
                                            "family_name", "genus_name", 
                                            "species_name", "lat", 
                                            "nucleotides", "markercode")]
dfOutgroup <- setDT(dfOutgroup)

# Combine dfOutgroup and dfInitial datatables so that they are in one useable 
# datatable.
l <- list(dfInitial, dfOutgroup)
dfFiltered <- rbindlist(l)

### FILTER 1 ###
# Filters are used for quality control purposes.
# Filtering for presence of a BIN URI. Assignment of a BIN URI is an indicator 
# of sequence quality. A BIN URI is also necessary for species assignment to the 
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
# filtering (i.e. # of original sequences).
dfFiltered <- dfFiltered[, initial_bin_size := length(recordID), keyby = bin_uri]
################################################################################

### FILTER 3 ###
# Filtering for COI-5P as these are the only markers we are looking at.
# Setting a key facilitates filtering of datatables.
setkey(dfFiltered, markercode)
# Now I can just specify "COI-5P" without indicating which column.
dfFiltered <- dfFiltered["COI-5P"]

### FILTER 4 ###
# N and gap content will interfere with the multiple sequence alignment and the 
# alignment will give warning messages, so we need to trim sequences with high N
# and gap content at their terminal ends.
# First, make sure nucleotides are "chr" type. This is necessary for regular 
# expression (regex) searches.
dfFiltered[, nucleotides := as.character(nucleotides)]
# Let's first trim large portions of Ns and gaps at the start of a sequence. 
# First, find sequences that begin (^) with an N or a gap.
startGapN <- sapply(regmatches(dfFiltered$nucleotides, 
                               gregexpr("^[-N]", 
                                        dfFiltered$nucleotides)), length)
# Loop through all of the sequences.
startGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
# If at least one sequence is found that begins with a gap or N (is flagged as a
# 1)...
  if (startGapN[[i]] > 0) { 
    # Split the sequence up using strsplit!
    # Using a regex to find sequences that begin with gaps and/or Ns.
    split <- strsplit(dfFiltered$nucleotides[i], "^[-N]+")
    # Take only the second half of the element (the sequence without the gaps/Ns
    # at the start!).
    dfFiltered$nucleotides[i] <- split[[1]][2]
  }
rm(startGapN)
# Now, let's trim large portions of Ns and gaps at the end of a sequence. 
# First, find sequences that end ($) with an N or a gap.
endGapN <- sapply(regmatches(dfFiltered$nucleotides, 
                             gregexpr("[-N]$", dfFiltered$nucleotides)), length)
endGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
  if (endGapN[[i]] > 0) {
    # Using a regex to find sequences that end with gaps and/or Ns.
    split <- strsplit(dfFiltered$nucleotides[i], "[-N]+$")
    # Take only the first half of the element (the sequence without the gaps/Ns 
    # at the end!).
    dfFiltered$nucleotides[i] <- split[[1]][1]
  }
rm(endGapN)
rm(split)

### FILTER 5 ###
# Remove sequences with N/gap content above a certain threshold. In this case,
# we will be removing sequences with greater than 1% gap/N content.
# First, let's determine the number of positions where an *internal* N or gap is
# found for each sequence.
internalGapN <- sapply(regmatches(dfFiltered$nucleotides, 
                                  gregexpr("[-N]", dfFiltered$nucleotides)), length)
# We then iterate over each sequence and see if the number of Ns or gaps is 
# greater than 1% (0.01) of the total length of the sequence.
internalGapN <- foreach(i = 1:nrow(dfFiltered)) %do% 
  which((internalGapN[[i]]/nchar(dfFiltered$nucleotides[i]) > 0.01))
# Here, we are basically "flagging" the sequences with high N/gap content. Those
# sequences will have values greater than 0.
checkGapN <- sapply(internalGapN, function (x) length(x))
# Identify the sequences with high N/gap content.
checkGapN <- which(checkGapN > 0) 
# Remove these higher gap/N content sequences.
dfFiltered <- dfFiltered[-checkGapN, ]
rm(checkGapN)
rm(i)
rm(internalGapN)

### FILTER 6 ###
# Filter out sequences that are less than 640 bp and greater than 1000 bp. 
# This is because extremely long or short sequence lengths can interfere with 
# multiple sequence alignments.
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
# Remove rows with no species information. This will thereby remove BINs without
# any species information (no rows with BINs without species information will 
# remain). BINs without species data would not match with any trait information 
# down the line.
containSpecies <- dfFiltered[, grep("[A-Z]", species_name)]
# Create a new datatable containing only sequences baring species-level 
# identification. This is necessary so we can extract the BIN URIs that contain 
# species-level identification and remove those without! 
dfSpecies <- dfFiltered[containSpecies, ]
rm(containSpecies)
# Now we have the BIN URIs of BINs that contain sequences with species 
# information.
speciesBins <- unique(dfSpecies$bin_uri)
# Subset out these BINs in dfFiltered.
dfResolve <- subset(dfFiltered, bin_uri %in% speciesBins)

# RESOLVING TAXONOMIC CONFLICTS.
# These steps are performed to improve BIN reliability and ensure we are 
# matching the appropriate sequence information to the appropriate trait data 
# down the line.
# Now, I want to resolve BINs with more than 1 order and/or family.
# First, I need to replace all blanks with NA values in the taxonomy columns.
# This is to ensure that blanks are not counted as their own taxa.
dfResolve[, order_name := revalue(order_name, c(" " = NA))]
dfResolve[, family_name := revalue(family_name, c(" " = NA))]
dfResolve[, genus_name := revalue(genus_name, c(" " = NA))]
dfResolve[, species_name := revalue(species_name, c(" " = NA))]

# Now find the number of orders/families/genera/species in each BIN.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), 
          keyby = bin_uri]

# FUNCTION: For removing BINs with taxonomic conflicts.
ResolveBIN <- function(x, y, method = c("bin_uri","recordID")){
  # x = Either the BIN or recordID to be removed.
  # y = Dataframe to apply function to.
  # method = Specify whether the item to be removed is a BIN or a recordID.
  method <- match.arg(method)
  if(method == "bin_uri") {
     rmThisBIN <- which(y$bin_uri == x)
     resolved <- y[-rmThisBIN, ]
  } else if(method == "recordID") {
     rmThisRecord <- which(y$recordID == x)
     resolved <- y[-rmThisRecord, ]
  }
  return(resolved)
}

# First, let's deal with order level conflicts.
# Checking these manually so I can double check the BOLD record online.
orderConflicts <- dfResolve[, which(number_of_orders > 1), by = bin_uri]
orderConflicts <- unique(orderConflicts$bin_uri)
length(orderConflicts)
# There doesn't seem to be any orderConflicts for this group.
# Update number_of_orders column and make sure there are no more order conflicts.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), 
          keyby = bin_uri]
# Making sure there aren't any order level conflicts.
orderConflicts <- dfResolve[, which(number_of_orders > 1), by = bin_uri]
orderConflicts <- unique(orderConflicts$bin_uri)
rm(orderConflicts)

# Family level resolving.
# Checking these manually.
# Make sure to update number_of_families column first.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), 
          keyby = bin_uri]
# Find number of BINs with family level conflicts.
familyConflicts <- dfResolve[, which(number_of_families > 1), by = bin_uri]
familyConflicts <- unique(familyConflicts$bin_uri)
length(familyConflicts)
# Removing BOLD:AAB2488. No clear consensus.
dfResolve <- ResolveBIN("AAB2488", dfResolve, method = "bin_uri")
# BOLD:AAC3188 removing 1 deviant seq from different fam.
dfResolve <- ResolveBIN(338753, dfResolve, method = "recordID")
# BOLD:AAD2909 removing 4 deviant seqs from different fam.
dfResolve <- ResolveBIN(210644, dfResolve, method = "recordID")
dfResolve <- ResolveBIN(210643, dfResolve, method = "recordID")
dfResolve <- ResolveBIN(210645, dfResolve, method = "recordID")
dfResolve <- ResolveBIN(210653, dfResolve, method = "recordID")
# Update number_of_orders column and make sure there are no more order conflicts.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), 
          keyby = bin_uri]
familyConflicts <- dfResolve[, which(number_of_families > 1), by = bin_uri]
familyConflicts <- unique(familyConflicts$bin_uri)
rm(familyConflicts)
rm(speciesBins)

# Genus level resolving.
# Update the column.
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), 
          keyby = bin_uri]
# At least 10+ records and 80% consistency.
genusConflicts <- dfResolve[, which(number_of_genera > 1), by = bin_uri]
genusConflictBins <- unique(genusConflicts$bin_uri)
length(genusConflictBins)
# Create a new datatable for BINs with genus level conflicts.
dfGenusConflicts <- dfResolve[genusConflictBins, ]
# Now we must determine the most common genus and if it has at least 80% 
# consistency in sequences that DO have genus level information.
# Only looking at sequences with genus classifications.
containGenus <- dfGenusConflicts[, grep("[A-Z]", genus_name)]
# Subset only those sequences baring genus-level identification. 
dfGenusConflicts <- dfGenusConflicts[containGenus, ]
# Create a new column for the number of sequences with genus level information 
# per BIN.
dfGenusConflicts[, number_of_seqs := length(recordID), by = bin_uri]
# Which bins have more than 10 sequences? These are probably more reliable.
dfGenusConflicts <- dfGenusConflicts[which(number_of_seqs >= 10)]
# A count column is created to count the number of rows per genus per BIN.
# This is necessary to calculate the percentage of sequences from each genus per
# BIN.
dfGenusConflicts[, count := .N, by = .(bin_uri, genus_name)]
# Calculate the percentage of sequences from each genus per BIN.
dfGenusConflicts[, genus_percentage := .(count / number_of_seqs)]
dfGenusConflicts[order(-count), majority_genus := genus_name[1L], by = bin_uri]
# Reorder to take a closer look.
dfGenusConflicts <- dfGenusConflicts[, .(bin_uri, genus_name, majority_genus, 
                                         number_of_seqs, count, 
                                         genus_percentage)]
# Make a column for majority species percentage to test if it is over 80%.
# This is the percentage for the genus with the majority of entries.
dfGenusConflicts[order(-genus_percentage), 
                 majority_genus_percentage := genus_percentage[1L], 
                 by = bin_uri]
# Subset out those BINs that have a majority genera over 80%.
dfAcceptedGenus <- dfGenusConflicts[which(majority_genus_percentage > 0.80)]
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# bad = BINs in genusConflicts which were not accepted.
bad <- anti_join(genusConflicts, dfAcceptedGenus, by = "bin_uri")
badBins <- unique(bad$bin_uri)
dfResolve <- dfResolve[!dfResolve$bin_uri %in% badBins, ]

# Species level resolving.
# Update the column.
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), 
          keyby = bin_uri]
# At least 10+ records and 80% consistency.
speciesConflicts <- dfResolve[, which(number_of_species > 1), by = bin_uri]
speciesConflictBins <- unique(speciesConflicts$bin_uri)
length(speciesConflictBins)
# Create a new datatable for BINs with genus level conflicts.
dfSpeciesConflicts <- dfResolve[speciesConflictBins, ]
# Now we must determine the most common species and if it has at least 80% 
# consistency in sequences that DO have species level information.
# Only looking at sequences with species classifications.
containSpecies <- dfSpeciesConflicts[, grep("[A-Z]", species_name)]
# Subset only those sequences baring genus-level identification.  
dfSpeciesConflicts <- dfSpeciesConflicts[containSpecies, ]
# Create a new column for the number of sequences with species level information
# per BIN.
dfSpeciesConflicts[, number_of_seqs := length(recordID), by = bin_uri]
# Which bins have more than 10 sequences? These are probably more reliable.
dfSpeciesConflicts <- dfSpeciesConflicts[which(number_of_seqs >= 10)]
# A count column is created to count the number of rows per species per BIN.
# This is necessary to calculate the percentage of sequences from each species 
# per BIN.
dfSpeciesConflicts[, count := .N, by = .(bin_uri, species_name)]
# Calculate the percentage of sequences from each species per BIN.
dfSpeciesConflicts[, species_percentage := .(count / number_of_seqs)]
dfSpeciesConflicts[order(-count), majority_species := species_name[1L], by = bin_uri]
# Reorder to take a closer look.
dfSpeciesConflicts <- dfSpeciesConflicts[, .(bin_uri, species_name, 
                                             majority_species, number_of_seqs, 
                                             count, species_percentage)]
# Make a column for majority species percentage to test if it is over 80.
# # This is the percentage for the genus with the majority of entries.
dfSpeciesConflicts[order(-species_percentage),
                   majority_species_percentage := species_percentage[1L], 
                   by = bin_uri]
# Subset out those BINs that have a majority species over 80%.
dfAcceptedSpecies <- dfSpeciesConflicts[which(majority_species_percentage > 0.80)]
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# bad = BINs in speciesConflicts which are not accepted.
bad <- anti_join(speciesConflicts, dfAcceptedSpecies, by = "bin_uri")
badBins <- unique(bad$bin_uri)
dfResolve <- dfResolve[!dfResolve$bin_uri %in% badBins, ]

################################################################################
### TRAIT: POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfResolve[, filtered_bin_size := length(recordID), by = bin_uri]
################################################################################

### TAXONOMIC LABELS ###
# Determine the most common taxonomic classifications in each BIN. 
# Create a new datatable containing only sequences baring taxonomic 
# identification at the corresponding level.
# This is necessary because NA values are considered when counting the number 
# of species.
containSpecies <- dfResolve[, grep("[A-Z]", species_name)]
dfSpecies <- dfResolve[containSpecies, ]
rm(containSpecies)
dfSpeciesLabel <- dfSpecies[, .N, by = .(bin_uri, species_name)][order(-N), .(species_label = species_name[1L]), keyby = bin_uri]
# Genus label.
containGenus <- dfResolve[, grep("[A-Z]", genus_name)]
dfGenus <- dfResolve[containGenus, ]
rm(containGenus)
dfGenusLabel <- dfGenus[, .N, by = .(bin_uri, genus_name)][order(-N), .(genus_label = genus_name[1L]), keyby = bin_uri]
rm(dfGenus)
# Family label. Technically, there should be only one family/order name per BIN
# after manually resolving these cases.
# But doing this step just in case!!!
containFamily <- dfResolve[, grep("[A-Z]", family_name)]
dfFamily <- dfResolve[containFamily, ]
rm(containFamily)
dfFamilyLabel <- dfFamily[, .N, by = .(bin_uri, family_name)][order(-N), .(family_label = family_name[1L]), keyby = bin_uri]
rm(dfFamily)
# Order label.
containOrder <- dfResolve[, grep("[A-Z]", order_name)]
dfOrder <- dfResolve[containOrder, ]
rm(containOrder)
dfOrderLabel <- dfOrder[, .N, by = .(bin_uri, order_name)][order(-N), .(order_label = order_name[1L]), keyby = bin_uri]
rm(dfOrder)

# MERGING DATATABLES.
# Merge datatables containing BIN species information.
# Set the key for dfResolve.
setkey(dfResolve, bin_uri)
# Note: bin_uri is the key for each of these datatables already due to the use 
# of "keyby" instead of just "by". This ultimately facilitates datatable merging.
# Merging multiple datatables at once.
dfFiltered <- Reduce(function(...) merge(...), list(dfResolve, dfSpeciesLabel, 
                                                    dfGenusLabel, dfFamilyLabel, 
                                                    dfOrderLabel))
# Datatable reorganization. Double check this datatable.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, order_label, 
                             family_name, family_label, genus_name, genus_label, 
                             species_name, species_label, nucleotides, 
                             initial_bin_size, filtered_bin_size, lat)]


################################################################################
### TRAITS: MEDIAN LATITUDE/LATITUDINAL RANGE ###
# Currently, median latitude and latitudinal range are the only traits whose 
# information is taken from BOLD. The rest of the data will be obtained from 
# FishBase.
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

# FUNCTION: Filtering for presence of trait data. This function is for the 
# univariate analyses section..
GetTraitSpecificData <- function(x, y) {
  # Filters a dataframe for only those data related to a specified trait.
  # x = dataframe of species and trait information.
  # y = trait.
  
  # Make sure x is in dataframe format.
  x <- as.data.frame(x)
  # Find rows without data for column.
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  # Construct the univariate trait datatable. This datatable will be used in the 
  # eventual univariate analysis.
  # If there are rows without data for column...
  if (length(noY) > 0) {
    # Remove the species without data for column.
    z <- x[-noY, ]
    # Reorganize datatable.
    # Column 1 = bin_uri, column 10 = species_label, 
    # column 13 = filtered_bin_size
    z <- z[c(1, 10, 13, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$bin_uri), ] 
    # If all rows have data for column...
  } else {
    # If no entries need to be removed, just rename and reorganize x.
    z <- x[c(1, 10, 13, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$bin_uri), ]  
  }
  rm(noY)
  return(z)
}

# While considering traits for eventual multivariate analyses, it is necessary
# for them to have an adequate sample size (i.e. over x # rows of data, depending
# on your purposes).
# In addition, they should exhibit some amount of variation across the 
# observations.
# TRAIT: Latitude.
# TEST 1: Does trait have data for at least 300 species?
dfLatitude <- setDT(GetTraitSpecificData(dfLatitudeSpecies, 17))
nrow(dfLatitude)
# TEST 2: Does the trait have enough data variation?
hist(dfLatitude$median_lat)

################################################################################
# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name = order_label, 
                             family_name = family_label, 
                             genus_name = genus_label, 
                             species_name = species_label, nucleotides, 
                             initial_bin_size, filtered_bin_size)]

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
# As there are multiple entries per species for some traits, I want to take the 
# median or modevalue of some traits. This will depend on the type of trait 
# i.e. continuous vs. categorical).
# Note: Only some traits are shown here as examples.

# Species version of GetTraitSpecificData (vs. BIN version).
GetTraitSpecificData <- function(x, y) {
  # Filters a dataframe for only those data related to a specified trait.
  # x = dataframe of species and trait information.
  # y = trait
  
  x <- as.data.frame(x)  # Still want to figure this out for datatable.
  # Find rows without data for column.
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  # Construct the univariate trait datatable. This datatable will be used in the eventual 
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

### SPECIES TRAITS ###
dfSpecies <- data.frame(species(speciesNames))
# Storing this as a file.
write.csv(dfSpecies, file = "species_datatable.csv")
dfSpecies <- fread("species_datatable.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, 
                                 BodyShapeI, LongevityWild, Length, LTypeMaxM,
                                 DepthRangeDeep)]
# TRAIT: Body Shape I.
# Filtering for presence of body shape data. This is for the univariate analyses 
# section.
# TEST 1: Does the trait have an adequate sample size?
dfBodyShapeI <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
nrow(dfBodyShapeI)
# TEST 2: Does the trait have enough data variation?
table(dfBodyShapeI$BodyShapeI)
# "other" is a rare category.
rareVars <- which(dfBodyShapeI$BodyShapeI == "eel-like")
dfBodyShapeI <- dfBodyShapeI[-rareVars, ]
# Make it a factor variable.
dfBodyShapeI[, BodyShapeI := as.factor(BodyShapeI)]
# Also dropping this level from the factor.
dfBodyShapeI$BodyShapeI <- droplevels(dfBodyShapeI$BodyShapeI)

# TRAIT: LongevityWild.
# TEST 1: Does the trait have an adequate sample size?
dfLongWild <- setDT(GetTraitSpecificData(dfSpeciesTraits, 3))
nrow(dfLongWild)
# TEST 2: Does the trait have enough data variation?
hist(dfLongWild$LongevityWild)
range(dfLongWild$LongevityWild)
# Make sure it is a numeric variable.
dfLongWild[, LongevityWild := as.double(LongevityWild)]

# TRAIT: Maximum length.
# The column must first be converted to double (numeric) type.
dfMaxLength <- dfSpeciesTraits[, Length := as.double(Length)]
# We only want total length measurements.
keep <- dfMaxLength[, which(LTypeMaxM == "TL")]
dfMaxLength <- dfMaxLength[keep, ]
# TEST 1: Does the trait have an adequate sample size?
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 4))
nrow(dfMaxLength)
# TEST 2: Does the trait have enough data variation?
hist(dfMaxLength$Length)
range(dfMaxLength$Length)

# TRAIT: Maximum depth.
# The column must first be converted to double (numeric) type.
dfMaxDepth <- dfSpeciesTraits[, DepthRangeDeep := as.double(DepthRangeDeep)]
# TEST 1: Does the trait have an adequate sample size?
dfMaxDepth <- setDT(GetTraitSpecificData(dfMaxDepth, 6))
nrow(dfMaxDepth)
# TEST 2: Does the trait have enough data variation?
hist(dfMaxDepth$DepthRangeDeep)
range(dfMaxDepth$DepthRangeDeep)

# Finally, prepare the dfSpeciesGenMV datatable by merging all univariate 
# datatables. This datatable will be used for the eventual multivariate analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(..., all = T), 
                         list(dfBodyShapeI, dfLongWild, dfMaxLength,
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
rm(dfMaturityTraits)
# The median value is then determined for each species.
dfAgeMaturity[, age_at_maturity := median(age_at_maturity, na.rm = TRUE), keyby = species_name]
# Filtering for presence of average depth data. This is for the univariate analyses section.
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
colnames(dfEcology)[3] <- "species_name"
# Get rid of columns I do not need for the regression analysis. 
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, c(3, 7, 15, 25:26, 31, 66)]

# As the columns are coded as integers right now, we need to recode them to the 
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
# Recode the integer variables (from "-1" to "1").
dfEcologyTraits[, (integerVars) := lapply(.SD, function(x) revalue(x, c("-1" = "1"))), 
                .SDcols = integerVars]
rm(integerVars)
rm(characterVars)
rm(changeVars)

### Categorical traits. ###
# Binary traits.
# TRAIT: Neritic.
# TEST 1: Does the trait have an adequate sample size?
dfNeritic <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))
nrow(dfNeritic)
# TEST 2: Does the trait have enough data variation?
table(dfNeritic$Neritic)

# TRAIT: Oceanic.
# TEST 1: Does the trait have an adequate sample size?
dfOceanic <- setDT(GetTraitSpecificData(dfEcologyTraits, 3))
nrow(dfOceanic)
# TEST 2: Does the trait have enough data variation?
table(dfOceanic$Oceanic)

# TRAIT: Stream.
# TEST 1: Does the trait have an adequate sample size?
dfStreams <- setDT(GetTraitSpecificData(dfEcologyTraits, 4))
nrow(dfStreams)
# TEST 2: Does the trait have enough data variation?
table(dfStreams$Stream)

# TRAIT: Lakes.
# TEST 1: Does the trait have an adequate sample size?
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 5))
nrow(dfLakes)
# TEST 2: Does the trait have enough data variation?
table(dfLakes$Lakes)

# TRAIT: Benthic.
# TEST 1: Does the trait have an adequate sample size?
dfBenthic <- setDT(GetTraitSpecificData(dfEcologyTraits, 7))
nrow(dfBenthic)
# TEST 2: Does the trait have enough data variation?
table(dfBenthic$Benthic)

# Non-binary traits.
# TRAIT: Feeding Type.
# TEST 1: Does the trait have an adequate sample size?
dfFeedingType <- setDT(GetTraitSpecificData(dfEcologyTraits, 6))
# TEST 2: Does the trait have enough data variation?
# Many categories do not reach the 1% threshold and are removed.
table(dfFeedingType$FeedingType)
rareVars <- which(dfFeedingType$FeedingType == "filtering plankton" |
                  dfFeedingType$FeedingType == "other" | 
                  dfFeedingType$FeedingType == "sucking food-containing material")
dfFeedingType <- dfFeedingType[-rareVars, ]
# Also dropping these levels from the factor.
dfFeedingType$FeedingType <- droplevels(dfFeedingType$FeedingType)

# Finally, prepare the dfEcologyMV datatable by merging all univariate 
# datatables.
# Merging multiple datatables at once.
dfEcologyMV <- Reduce(function(...) merge(..., all = T), 
                      list(dfNeritic, dfOceanic, dfStreams, 
                           dfLakes, dfFeedingType, dfBenthic))

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
dfTraits <- dfTraits[, c(1:4, 7:18)]
colnames(dfTraits)[4] <- "filtered_bin_size"

# Merge back to dfFiltered to obtain all of the sequence information for 
# each BIN. This is for creation of the master phylogeny.
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[6] <- "species_name"
colnames(dfPreCentroid)[8] <- "initial_bin_size"
colnames(dfPreCentroid)[9] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:9)]

################################################################################
### SECTION 3: CENTROID SEQUENCE DETERMINATION ###
# Centroid Sequence: BIN sequence with minimum average pairwise distance to all 
# other sequences in a given BIN.

# Subset dataframe to find BINs with more than one sequence.These are the BINs 
# that will require centroid sequences.
largeBins <- which(dfPreCentroid$filtered_bin_size > 1)
  
# If there is at least one BIN with more than one member, then a dataframe 
# dfCentroidSeqs will be created with those BINs.
if (length(largeBins) > 0) {
  
  # Remove gaps from the sequences.
  dfPreCentroid[, nucleotides := gsub("-", "", nucleotides)] 
  
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
  
  # Alignment using muscle.
  alignment1 <- lapply(DNAStringSet1, function(x) 
    muscle::muscle(x, diags = TRUE, gapopen = -3000))
  
  # We can then convert each alignment to DNAbin format.
  dnaBINCentroid <- lapply(alignment1, function(x) as.DNAbin(x))
  #rm(alignment1)
    
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

# Make sure there is only a single row per record/BIN.
dfAllSeqs <- dfAllSeqs[!duplicated(dfAllSeqs$recordID), ] 
dfAllSeqs <- dfAllSeqs[!duplicated(dfAllSeqs$bin_uri), ]

# Take BINs with highest number of sequences with species info in order to avoid
# using BINs that were assigned the same species label.
dfAllSeqs <- merge(aggregate(filtered_bin_size ~ species_name, data = dfAllSeqs, 
                             max), dfAllSeqs, all.x = T, sort = TRUE)
dup_majority_species <- which(duplicated(dfAllSeqs$species_name)) # Dealing with ties.
dfAllSeqs <- dfAllSeqs[-dup_majority_species,]
rm(dup_majority_species)

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
  #data$nucleotides <- data[, gsub("-", "", nucleotides)]

  # We must ensure that the sequences are of the chr type when all of the 
  # sequences PLUS the reference sequence(s) are combined into a vector. The
  # reference sequence is added as the first sequence.
  alignmentSeqs <- as.character(data$nucleotides)

  # Name our sequences according to species names.
  names(alignmentSeqs) <- data$species_name
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
  #alignment2 <- AlignSeqs(DNAStringSet2, gapOpening = -3000)
  
  # If you want to save the alignment as a FASTA file to your current working
  # directory, uncomment the following lines. The file will be named according to 
  # the class of organisms whose barcode sequences you are you currently 
  # analyzing.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentUntrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  # Convert to DNAStringSet format.
  alignmentUntrimmed <- DNAStringSet(alignment2)
  writeXStringSet(alignmentUntrimmed, file = classFileNames[[1]], 
                  format = "fasta")

  # For trimming of the sequences, we have to determine where in the alignment 
  # the reference sequence is and determine its start and stop positions 
  # relative to the other sequences. We can then use these positions to trim 
  # the rest of the sequences in the alignment.
  refSeqPos <- which(alignment2@unmasked@ranges@NAMES == "REFERENCE")
  refSeqPos <- alignment2@unmasked[refSeqPos]

  # Find the start position by searching for the first nucleotide position of 
  # the reference sequence.
  refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
  refSeqPosStart <- as.numeric(refSeqPosStart)

  # Find last nucleotide position of the reference sequence.
  refSeqPosEnd <- nchar(dfRefSeqs$nucleotides[1]) + refSeqPosStart
  refSeqPosEnd <- as.numeric(refSeqPosEnd)

  # Then we can take a substring of the alignment by using 
  # these positions to effectively "trim" the alignment.
  alignment2Trimmed <- substr(alignment2, refSeqPosStart, refSeqPosEnd)

  # Again, convert to DNAStringSet format.
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)
  
  # Checking alignment.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentTrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  writeXStringSet(DNAStringSet3, file = classFileNames[[1]], 
                  format = "fasta", width = 1500)

  # Remove the reference sequence from this, as we dont want it to be included 
  # in further analysis.
  refSeqRm <- which(DNAStringSet3@ranges@NAMES == "REFERENCE")
  DNAStringSet3 <- subset(DNAStringSet3[-refSeqRm])

  # Reorder dfAllSeqs according to the order of species produced by the 
  # alignment, which are now contained in the DNA_StringSet object.
  # Make a variable with the ordering.
  alignmentOrder <- DNAStringSet3@ranges@NAMES

  # Order dfAllSeqs according to this.
  data <- data[match(alignmentOrder, data$species_name), ]

  # Repopulate dfAllSeqs with the newly trimmed sequences instead of the raw 
  # sequences.
  trimmedSeqs <- as.character(DNAStringSet3)
  data$nucleotides <- trimmedSeqs
  # Make sure species and sequences are correctly matched up!!!
  
  return(data)
}

dfCheckAllSeqs <- refSeqTrim(dfAllSeqs)

################################################################################
### SECTION 4: ALIGNMENT QUALITY CHECKING ###
# Here, extremely gappy/ungappy sequences are removed. These sequences are 
# assumed to contribute to misalignment of the sequences or may even be 
# pseudogenes. Manually checking of the alignment is recommended.
# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGaps <- sapply(regmatches(dfCheckAllSeqs$nucleotides, 
                                  gregexpr("[-+]", 
                                           dfCheckAllSeqs$nucleotides)), length)
# Mean gap length and range.
meanGap <- mean(internalGaps)
extremeHighGap <- meanGap + 7  # Upper range.
extremeLowGap <- meanGap - 7  # Lower range.
rm(meanGap)
# We then loop through each sequence to see if the number of gaps deviates 
# greatly from the mean.
# Which sequences exceed the range of meanGap +/- 7?
extremeSeqs <- foreach(i = 1:nrow(dfCheckAllSeqs)) %do% 
  which(internalGaps[[i]] > extremeHighGap | internalGaps[[i]] < extremeLowGap)
rm(extremeHighGap)
rm(extremeLowGap)
rm(internalGaps)
# The "deviant" sequences will be flagged with a 1.
extremeBins <- which(extremeSeqs > 0)
# Subset out these sequences to look at them if desired.
dfExtreme <- dfCheckAllSeqs[extremeBins, ] # Need BIN names not index numbers.
# Make sure outgroups are not removed.
goodBins <- which(dfExtreme$order_name == "Acanthuriformes")
dfExtreme <- dfExtreme[-goodBins, ]
extremeBins <- dfExtreme$bin_uri
# Remove the gappy sequences from dfAllSeqs as we will be retrimming these 
# sequences again once troublesome cases are removed.
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% extremeBins, ]
rm(extremeBins)

### OUTLIER CHECK ###
# Remove centroid sequences whose genetic distances to all other sequences fall 
# outside the typical range of genetic divergence for this group of organisms.
# Convert each alignment to DNAbin format.
dnaBinNN <- DNAStringSet(dfCheckAllSeqs$nucleotides)
names(dnaBinNN) <- dfCheckAllSeqs$bin_uri
dnaBinNN <- as.DNAbin(dnaBinNN)
# Then, we perform genetic distance determination with the TN93 model.
geneticDistanceCentroid <- dist.dna(dnaBinNN, model = "TN93", as.matrix = TRUE, 
                                    pairwise.deletion = TRUE)
# Using the IQR to detect outliers.
lowerq <- quantile(geneticDistanceCentroid)[2]
upperq <- quantile(geneticDistanceCentroid)[4]
iqr <- upperq - lowerq
mild.threshold.upper <- (iqr * 1.5) + upperq
# Remove 0 values (when a BIN is compared to itself - the diagonal values).
geneticDistanceCentroid[geneticDistanceCentroid == 0] <- NA
# Convert to dataframe format.
dfOutlierCheck <- as.data.frame(geneticDistanceCentroid)
# Identify BINs with no relatives within "typical" range of genetic divergence
# (i.e. all of their genetic distances fall outside 1.5 x IQR upper threshold.)
outliers <- apply(dfOutlierCheck, MARGIN = 1, 
                  function(x) all(x > mild.threshold.upper))
typical <- which(outliers == "FALSE")
# Take a closer look.
dfOutliers <- dfOutlierCheck[-typical, ]
outliers <- row.names(dfOutliers)
outliers <- which(dfAllSeqs$bin_uri == outliers)
# Remove the outliers from dfAllSeqs as we will be retrimming these sequences 
# again once troublesome cases are removed.
dfAllSeqs <- dfAllSeqs[-outliers, ]

### NEAREST NEIGHBOUR CHECK ###
# Remove centroid sequences whose nearest neighbours are in a different order or 
# family. The nearest neighbour can be determined from the distance matrix alone.
# It is the sequence with minimum pairwise distance to the sequence in question.

# FUNCTION: Identifies the nearest neighbour of each BIN.
# Applying this function to geneticDistanceCentroid.
NearestNeighbour <- t(sapply(seq(nrow(geneticDistanceCentroid)), function(i) {
  # Identify the sequence with the minimum distance in the row (the nearest
  # neighbour of i).
  j <- which.min(geneticDistanceCentroid[i, ])
  # Assign the name of the BIN and the name of its nearest neighbour to a 
  # separate column.
  c(paste(rownames(geneticDistanceCentroid)[i], 
          colnames(geneticDistanceCentroid)[j], sep='/'), 
    geneticDistanceCentroid[i, j])
}))
# Convert to dataframe.
dfNearestNeighbour <- as.data.frame(NearestNeighbour, stringsAsFactors = FALSE)
# Split the columns using the transform() function.
dfNearestNeighbour <- transform(dfNearestNeighbour, 
                                V1 = colsplit(V1, split = "/", 
                                              names = c('bin_uri', 
                                                        'nearest_neighbour')))
# Make sure the genetic distances are numeric.
dfNearestNeighbour$V2 <- as.double(dfNearestNeighbour$V2)
# Now we need to get the family and order names of the BINs and their nearest 
# neighbours.
# Convert into another dataframe first (we need to do this to fully separate
# the bin_uri and nearest_neighbour columns).
df2NearestNeighbour <- dfNearestNeighbour$V1
df2NearestNeighbour$pairwise_distance <- dfNearestNeighbour$V2
# Now extract the bin_uri and nearest_neighbour columns as vectors.
bins <- as.character(df2NearestNeighbour$bin_uri)
nn <- as.character(df2NearestNeighbour$nearest_neighbour)
# Make a dataframe called dfBinOrd and match the order of the bin_uris in this 
# dataframe to the order of the bin_uris in dfCheckAllSeqs.
dfBinOrd <- dfCheckAllSeqs[match(df2NearestNeighbour$bin_uri, 
                                dfCheckAllSeqs$bin_uri), ]
# Extract the order and family names.
bin_orders <- dfBinOrd$order_name
bin_families <- dfBinOrd$family_name
# Do the same for the nearest neighbour bins.
dfNeighbourOrd <- dfCheckAllSeqs[match(df2NearestNeighbour$nearest_neighbour, 
                                       dfCheckAllSeqs$bin_uri), ] 
nn_orders <- dfNeighbourOrd$order_name
nn_families <- dfNeighbourOrd$family_name
# Add these columns to df2NearestNeighbour.
df2NearestNeighbour$bin_order <- bin_orders
df2NearestNeighbour$bin_family <- bin_families
df2NearestNeighbour$nn_order <- nn_orders
df2NearestNeighbour$nn_family <- nn_families
df2NearestNeighbour <- setDT(df2NearestNeighbour)
# Non-matching orders.
nonmatchOrd <- which(df2NearestNeighbour$bin_order != df2NearestNeighbour$nn_order)
# Look closely at the cases where the orders do not match.
dfNonmatchOrd <- df2NearestNeighbour[nonmatchOrd, ]
# Which are very close neighbours (i.e. there may be something weird going on)?
dfNonmatchOrd <- dfNonmatchOrd[which(pairwise_distance < 0.05)]
# There doesn't seem to be any weird cases with the current taxa. But if there 
# were, lines 1137-1138 may be used to remove weird BINs (you should double 
# check these cases on the BOLD website, too) from dfAllSeqs.
#nonMatchOrdBins <- dfNonmatchOrd$bin_uri
#dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% nonMatchOrdBins, ]
# Now, let's take a look at the BINs where the nearest neighbour is in a 
# different family.
nonmatchFam <- which(df2NearestNeighbour$bin_family != df2NearestNeighbour$nn_family)
# Which are very close neighbours?
dfNonmatchFam <- df2NearestNeighbour[nonmatchFam, ]
dfNonmatchFam <- dfNonmatchFam[which(pairwise_distance < 0.05)]
# Again, there doesn't seem to be any weird cases with the current taxa. But if 
# there were, lines 1149-1150 may be used to remove weird BINs (you should 
# doublecheck these cases on the BOLD website, too) from dfAllSeqs.
#nonMatchFamBins <- dfNonmatchFam$bin_uri
#dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% nonMatchFamBins, ]

# Now, we will be performing a more comprehensive check, just in case. Do the
# bins have ANY close neighbours within a genetic distance of 0.05 that are in a
# different order or family? Again, this may be indicative of something weird
# going on in either the sequence data or taxonomic assignment.
# First, subset out all close neighbour pairings that share a genetic distance
# under 0.05 to any other sequence.
dfGeneticDistance <- as.data.frame(geneticDistanceCentroid)
closeNeighbours <- which(dfGeneticDistance < 0.05, arr.ind = TRUE)
# Extract the row names (the names of the BINs).
bins <- row.names(dfGeneticDistance[closeNeighbours[, 1],] )
# Extract the column names (the names of the neighbour BINs).
nn <- colnames(dfGeneticDistance[closeNeighbours[, 2]])
# Convert to a dataframe.
dfCloseNeighbours <- as.data.frame(bins)
dfCloseNeighbours$bins <- as.character(dfCloseNeighbours$bins)
dfCloseNeighbours$neighbour <- nn
# Following the same protocol as before, extract the order and family names of 
# the BINs and their close neighbours.
dfBinOrd <- dfCheckAllSeqs[match(dfCloseNeighbours$bins, 
                                 dfCheckAllSeqs$bin_uri), ] 
bin_orders <- dfBinOrd$order_name
bin_families <- dfBinOrd$family_name
dfNeighbourOrd <- dfCheckAllSeqs[match(dfCloseNeighbours$neighbour, 
                                       dfCheckAllSeqs$bin_uri), ] 
nn_orders <- dfNeighbourOrd$order_name
nn_families <- dfNeighbourOrd$family_name
# Add these columns to dfCloseNeighbours.
dfCloseNeighbours$bin_order <- bin_orders
dfCloseNeighbours$bin_family <- bin_families
dfCloseNeighbours$nn_order <- nn_orders
dfCloseNeighbours$nn_family <- nn_families
# Convert to datatable.
dfCloseNeighbours <- setDT(dfCloseNeighbours)
# Some BINs have multiple close neighbours and were assigned extra numbers in 
# the dataframe (e.g. ACFGGF, ACFGGF.1). To resolve this, use the function 
# na.locf() from the "zoo" package to fill in missing data based on the data in 
# the row above.
# Make sure they are in the correct order in the "bins" column.
dfCloseNeighbours <- dfCloseNeighbours[order(dfCloseNeighbours$bins), ]
#install.packages("zoo")
library(zoo)
# Fill in the bin_order column.
dfCloseNeighbours$bin_order <- na.locf(dfCloseNeighbours$bin_order)
# Fill in the bin_family column.
dfCloseNeighbours$bin_family <- na.locf(dfCloseNeighbours$bin_family)
# Do the same for the neighbour BINs.
dfCloseNeighbours <- dfCloseNeighbours[order(dfCloseNeighbours$neighbour), ]
dfCloseNeighbours$nn_order <- na.locf(dfCloseNeighbours$nn_order)
dfCloseNeighbours$nn_family <- na.locf(dfCloseNeighbours$nn_family)
# Now, which orders do not match between bin_order and nn_order?
nonmatchOrd <- which(dfCloseNeighbours$bin_order != dfCloseNeighbours$nn_order)
# Which families do not match between bin_family and nn_family?
nonmatchFam <- which(dfCloseNeighbours$bin_family != dfCloseNeighbours$nn_family)
# Looking good, no mismatches!

### OUTGROUP CHECK ###
# Which outgroups made the cut? Remove them from dfAllSeqs so I can build tree
# just using the ingroup (so that inclusion of outgroups in the tree building
# process doesn't affect the branch length estimates of the in-group).
outgroupSpecies <- unique(dfOutgroup$species_name)
dfGoodOutgroups <- dfAllSeqs[dfAllSeqs$species_name %in% outgroupSpecies, ]
# Remove the outgroups from dfAllSeqs.
dfAllSeqsNO <- dfAllSeqs[!dfAllSeqs$species_name %in% outgroupSpecies, ]
# Now, re-trim and align the sequences without the outgroups.
dfAllSeqsNO <- refSeqTrim(dfAllSeqsNO)
# Once finished, make sure to check over sequences/alignment, and make sure they
# are in the correct reading frame. Make sure to save the resulting alignments
# under a different name, or save in a new directory so they are not replaced.

# Now re-run the alignment including outgroups (pick outgroup species that are
# well represented and that serve as appropriate outgroup to your taxa).
goodOG <- which(dfAllSeqs$species_name == "Boreogadus saida" |
                dfAllSeqs$species_name == "Merluccius merluccius")
mergeOG <- dfAllSeqs[goodOG, ]
# Add them back.
dfAllSeqsWithOG <- rbind(dfAllSeqsNO, mergeOG)
# Run the alignment.
dfAllSeqsWithOG <- refSeqTrim(dfAllSeqsWithOG)

################################################################################
### SECTION 5: STATISTICAL ANALYSES ###
# This is the section where the univariate analyses are performed as a form
# of model selection.

# First, read in your phylogenetic tree.
tree <- read.tree(file = "RAxML_labelledTree.EPAPercTree3")
tree$tip.label <- gsub("_", " ", tree$tip.label)

# Root the tree.
outgroups <- c("QUERY   Boreogadus saida")
tree <- root(tree, outgroup = outgroups, resolve.root = TRUE)
# Re-name the outgroup.
tree$tip.label <- gsub("QUERY   Boreogadus saida", "Boreogadus saida", tree$tip.label)

### TRAIT: NUMBER OF NODES.
# Let's first determine the number of nodes for each species. This will be used
# as a control variable in the multiple regression analysis (to account for the
# node density effect).
number_of_nodes <- distRoot(tree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(number_of_nodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)
# Merge with dfAllSeqsNode.
dfAllSeqsNode <- merge(dfAllSeqsWithOG, dfNumberOfNodes, by = "species_name")

### TRAIT: BRANCH LENGTHS.
# Let's calculate the sum of branch lengths now (from root to tip). These values
# will serve as out measurement of molecular evolution rate.
branch_lengths <- distRoot(tree, tips = "all", method = "patristic")
# Make into a dataframe.
dfBranchLengths <- data.frame(branch_lengths)
colnames(dfBranchLengths)[1] <- "branch_length"
dfBranchLengths$species_name <- names(branch_lengths)
# Check the distribution of the branch lengths.
hist(dfBranchLengths$branch_length, main = "", xlab = "Branch Length")
# Some more stats.
median(dfBranchLengths$branch_length)
mean(dfBranchLengths$branch_length)
range(dfBranchLengths$branch_length)
# Merge dfBranchLengths and dfAllSeqsNode.
dfRegression <- merge(dfBranchLengths, dfAllSeqsNode, by = "species_name")
# Merge back to the traits dataframe.
dfRegression <- merge(dfRegression, dfTraits, all.x = T, by = "bin_uri")
# Some reordering.
dfRegression <- dfRegression[, c(1, 6:8, 2:4, 11, 15:26)]
# Renaming for consistency purposes.
colnames(dfRegression)[5] <- "species_name"
colnames(dfRegression)[6] <- "branch_length"
colnames(dfRegression)[7] <- "bin_size"
colnames(dfRegression)[8] <- "number_of_nodes"
colnames(dfRegression)[10] <- "body_shape"
colnames(dfRegression)[11] <- "longevity"
colnames(dfRegression)[12] <- "max_length"
colnames(dfRegression)[13] <- "max_depth"
colnames(dfRegression)[14] <- "neritic"
colnames(dfRegression)[15] <- "oceanic"
colnames(dfRegression)[16] <- "streams"
colnames(dfRegression)[17] <- "lakes"
colnames(dfRegression)[18] <- "feeding_type"
colnames(dfRegression)[19] <- "benthic"


### Univariate analyses ###
# Running a univariate PGLS regression analysis for each trait to determine
# whether significance can be detected. If so, they will be included in the 
# multiple regression model selection process.

# First, make sure the trait data and phylo tree are in the same order.
# Make sure the order of the data matches the tree.
tree <- drop.tip(phy = tree, 
                     tip = tree$tip.label[!tree$tip.label%in%dfRegression$species_name], 
                     rooted = T)
dfRegression <- dfRegression[match(tree$tip.label, 
                                   dfRegression$species_name), ]

## 1. Branch length.
# As branch length is our response variable, we will only be estimating Pagel's 
# lambda, which is a measure of phylogenetic signal.
branch_length <- dfRegression$branch_length
names(branch_length) <- tree$tip.label
sigBL <- phylosig(tree, branch_length, method = "lambda", test = TRUE)

## 2. Number of nodes.
# Estimation of Pagel's lambda.
number_of_nodes <- dfRegression$number_of_nodes
names(number_of_nodes) <- tree$tip.label
sigNodes <- phylosig(tree, number_of_nodes, method = "lambda", test = TRUE)
# Make a univariate dataframe for number_of_nodes.
dfNodes <- dfRegression[, c(5:6, 8)]
# Perform a PGLS analysis using caper.
# First, we must make a comparative.data object to ensure the tree tips and data
# match.
c_data <- comparative.data(tree, dfNodes, "species_name")
caperNodes <- pgls(branch_length ~ number_of_nodes, c_data, lambda = "ML")

## 3. Median latitude.
# Pagel's lambda. A test for phylogenetic signal.
median_lat <- dfRegression$median_lat
names(median_lat) <- tree$tip.label
sigLat <- phylosig(tree, median_lat, method = "lambda", test = TRUE)
# Merge univariate dataframe to dfRegression to get branch length data.
dfLatitude <- merge(dfLatitude, dfRegression, by.x = "species_label", 
                    by.y = "species_name")
dfLatitude <- dfLatitude[, c(1, 5, 10:11)]
# As latitude was determined by BIN and not by species, take the BIN with
# the highest number of sequences with species info in order to avoid
# using BINs that were assigned the same species label.
dfLatitude <- merge(aggregate(filtered_bin_size ~ species_label, data = dfLatitude, max), 
                    dfLatitude, all.x = T, sort = TRUE)
# Dealing with ties.
dup_majority_species <- which(duplicated(dfLatitude$species_label))
dfLatitude <- dfLatitude[-dup_majority_species,]
# Test using caper.
c_data <- comparative.data(tree, dfLatitude, "species_label")
caperLatitude <- pgls(branch_length ~ median_lat.x + number_of_nodes, c_data, 
                      lambda = "ML")

## 4. Maximum length.
# Pagel's lambda. A test for phylogenetic signal.
maximum_length <- dfRegression$max_length
names(maximum_length) <- tree$tip.label
sigLength <- phylosig(tree, maximum_length, method = "lambda", test = TRUE)
# Merge univariate dataframe to dfRegression to get branch length data.
dfMaxLength <- merge(dfMaxLength, dfRegression, by = "species_name")
dfMaxLength <- dfMaxLength[, c(1, 7, 9, 13)]
# Make sure it is of the dataframe format.
dfMaxLength <- as.data.frame(dfMaxLength)
# Test using caper.
c_data <- comparative.data(tree, dfMaxLength, "species_name")
caperLength <- pgls(branch_length ~ number_of_nodes + max_length, c_data, 
                    lambda = "ML")

## 5. Longevity.
# Pagel's lambda. A test for phylogenetic signal.
longevity <- dfRegression$longevity
names(longevity) <- tree$tip.label
sigLong <- phylosig(tree, longevity, method = "lambda", test = TRUE)
# Merge univariate dataframe to dfRegression to get branch length data.
dfLongWild <- merge(dfLongWild, dfRegression, by = "species_name")
dfLongWild <- dfLongWild[, c(1, 7, 9, 12)]
# Make sure it is of the dataframe format.
dfLongWild <- as.data.frame(dfLongWild)
# Test using caper.
c_data <- comparative.data(tree, dfLongWild, "species_name")
caperLong <- pgls(branch_length ~ number_of_nodes + longevity, c_data, 
                  lambda = "ML")

## 6. Age at maturity.
# Pagel's lambda. A test for phylogenetic signal.
age_at_mat <- dfRegression$age_at_maturity
names(age_at_mat) <- tree$tip.label
sigAge <- phylosig(tree, age_at_mat, method = "lambda", test = TRUE)
# Merge univariate dataframe to dfRegression to get branch length data.
dfAgeMaturity <- merge(dfAgeMaturity, dfRegression, by = "species_name")
dfAgeMaturity <- dfAgeMaturity[, c(1, 7, 9, 21)]
# Make sure it is of the dataframe format.
dfAgeMaturity <- as.data.frame(dfAgeMaturity)
# Test using caper.
c_data <- comparative.data(tree, dfAgeMaturity, "species_name")
caperAge <- pgls(branch_length ~ number_of_nodes + age_at_maturity.y, c_data, 
                 lambda = "ML")

# LEFT OFF HERE!

# Discrete traits.
## 7. Neritic.
# Merge so can get branch length info for UV.
dfNeritic <- merge(dfNeritic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
neriticTree <- drop.tip(phy = tree, 
                           tip = tree$tip.label[!tree$tip.label%in%dfNeritic$species_name], 
                        rooted = T)
dfNeritic <- dfNeritic[match(neriticTree$tip.label, dfNeritic$species_name), ]
dfNeritic <- as.data.frame(dfNeritic)
row.names(dfNeritic) <- dfNeritic$species_name
dfNeritic <- dfNeritic[, c(1, 15, 7, 9)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfNeritic$neritic <- relevel(dfNeritic$neritic, ref = "0")
dfRegression$neritic <- relevel(dfRegression$neritic, ref = "0")
# D metric for phylogenetic signal.
sigNer <- phylo.d(dfNeritic, neriticTree, names.col = species_name, 
                  binvar = neritic)
# Try using caper.
c_data <- comparative.data(tree, dfNeritic, "species_name")
caperNeritic <- pgls(branch_length ~ neritic + number_of_nodes, c_data, 
                     lambda = "ML")

## 8. Oceanic.
#oceanic <- dfRegression$Oceanic
#names(oceanic) <- tree$tip.label
#fit <- fitDiscrete(tree, oceanic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, 
                                   dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfOceanic <- merge(dfOceanic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
oceanicTree <- drop.tip(phy = tree, 
                        tip = tree$tip.label[!tree$tip.label%in%dfOceanic$species_name], 
                        rooted = T)
dfOceanic <- dfOceanic[match(oceanicTree$tip.label, dfOceanic$species_name), ]
dfOceanic <- as.data.frame(dfOceanic)
row.names(dfOceanic) <- dfOceanic$species_name
dfOceanic <- dfOceanic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfOceanic$Oceanic.x <- relevel(dfOceanic$Oceanic.x, ref = "0")
dfRegression$Oceanic <- relevel(dfRegression$Oceanic, ref = "0")
# Try using caper.
c_data <- comparative.data(tree, dfOceanic, "species_name")
caperOceanic <- pgls(branch_length ~ Oceanic.x + number_of_nodes, c_data, 
                     lambda = "ML")
# D metric for phylogenetic signal.
sigOceanic <- phylo.d(dfOceanic, oceanicTree, names.col = species_name, 
                      binvar = Oceanic.x)

## 9. Benthic.
#benthic <- dfRegression$Benthic
#names(benthic) <- tree$tip.label
#fit <- fitDiscrete(tree, benthic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, 
                                   dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBenthic <- merge(dfBenthic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
benthicTree <- drop.tip(phy = tree, 
                        tip = tree$tip.label[!tree$tip.label%in%dfBenthic$species_name], 
                        rooted = T)
dfBenthic <- dfBenthic[match(benthicTree$tip.label, dfBenthic$species_name), ]
dfBenthic <- as.data.frame(dfBenthic)
row.names(dfBenthic) <- dfBenthic$species_name
dfBenthic <- dfBenthic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBenthic$Benthic.x <- relevel(dfBenthic$Benthic.x, ref = "0")
dfRegression$Benthic <- relevel(dfRegression$Benthic, ref = "0")
# Try using caper.
c_data <- comparative.data(tree, dfBenthic, "species_name")
caperBenthic <- pgls(branch_length ~ Benthic.x + number_of_nodes, c_data, 
                     lambda = "ML")
# D metric for phylogenetic signal.
sigBenthic <- phylo.d(dfBenthic, benthicTree, names.col = species_name, 
                      binvar = Benthic.x)

## 10. Lakes.
#lakes <- dfRegression$Lakes
#names(lakes) <- tree$tip.label
#fit <- fitDiscrete(tree, lakes, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfLakes <- merge(dfLakes, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
lakeTree <- drop.tip(phy = tree, 
                      tip = tree$tip.label[!tree$tip.label%in%dfLakes$species_name], rooted = T)
dfLakes <- dfLakes[match(lakeTree$tip.label, dfLakes$species_name), ]
dfLakes <- as.data.frame(dfLakes)
row.names(dfLakes) <- dfLakes$species_name
dfLakes <- dfLakes[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfLakes$Lakes.x <- relevel(dfLakes$Lakes.x, ref = "0")
dfRegression$Lakes <- relevel(dfRegression$Lakes, ref = "0")
# Try using caper.
c_data <- comparative.data(tree, dfLakes, "species_name")
caperLakes <- pgls(branch_length ~ Lakes.x + number_of_nodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigLakes <- phylo.d(dfLakes, lakeTree, names.col = species_name, binvar = Lakes.x)

## 11. Streams.
#streams <- dfRegression$Stream
#names(streams) <- tree$tip.label
#fit <- fitDiscrete(tree, streams, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfStreams <- merge(dfStreams, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
streamTree <- drop.tip(phy = tree, 
                      tip = tree$tip.label[!tree$tip.label%in%dfStreams$species_name], rooted = T)
dfStreams <- dfStreams[match(streamTree$tip.label, dfStreams$species_name), ]
dfStreams <- as.data.frame(dfStreams)
row.names(dfStreams) <- dfStreams$species_name
dfStreams <- dfStreams[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfStreams$Stream.x <- relevel(dfStreams$Stream.x, ref = "0")
dfRegression$Stream <- relevel(dfRegression$Stream, ref = "0")
# Try using caper.
c_data <- comparative.data(tree, dfStreams, "species_name")
caperStreams <- pgls(branch_length ~ Stream.x + number_of_nodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigStreams <- phylo.d(dfStreams, streamTree, names.col = species_name, binvar = Stream.x)

## 12. Body Shape.
#body_shape_I <- dfRegression$BodyShapeI
#names(body_shape_I) <- tree$tip.label
#fit <- fitDiscrete(tree, body_shape_I, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBodyShapeI <- merge(dfBodyShapeI, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
bodyTree <- drop.tip(phy = tree, 
                       tip = tree$tip.label[!tree$tip.label%in%dfBodyShapeI$species_name], rooted = T)
dfBodyShapeI <- dfBodyShapeI[match(bodyTree$tip.label, dfBodyShapeI$species_name), ]
dfBodyShapeI <- as.data.frame(dfBodyShapeI)
row.names(dfBodyShapeI) <- dfBodyShapeI$species_name
dfBodyShapeI <- dfBodyShapeI[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBodyShapeI$BodyShapeI.x <- relevel(dfBodyShapeI$BodyShapeI.x, ref = "fusiform / normal")
dfRegression$BodyShapeI <- relevel(dfRegression$BodyShapeI, ref = "fusiform / normal")
# Try using caper.
c_data <- comparative.data(tree, dfBodyShapeI, "species_name")
caperBodyShapeI <- pgls(branch_length ~ BodyShapeI.x + number_of_nodes, c_data, lambda = "ML")
anovaBSI <- anova.pgls.fixed(caperBodyShapeI)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyBs <- model.matrix(~ BodyShapeI.x - 1, data=dfBodyShapeI)
row.names(dfDummyBs) <- dfBodyShapeI$species_name
dfDummyBs <- as.data.frame(dfDummyBs)
dfDummyBs$species_name <- row.names(dfDummyBs)
colnames(dfDummyBs)[1] <- "Fusiform"
colnames(dfDummyBs)[2] <- "EelLike"
colnames(dfDummyBs)[3] <- "Elongated"
colnames(dfDummyBs)[4] <- "Short"
sigBsFusiform <- phylo.d(dfDummyBs, bodyTree, names.col = species_name, binvar = Fusiform)
sigBsEel <- phylo.d(dfDummyBs, bodyTree, names.col = species_name, binvar = EelLike)
sigBsElong <- phylo.d(dfDummyBs, bodyTree, names.col = species_name, binvar = Elongated)
sigBsShort <- phylo.d(dfDummyBs, bodyTree, names.col = species_name, binvar = Short)

## 13. Feeding Type.
#feeding_type <- dfRegression$FeedingType
#names(feeding_type) <- tree$tip.label
#fit <- fitDiscrete(tree, feeding_type, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(tree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfFeedingType <- merge(dfFeedingType, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
feedingTypeTree <- drop.tip(phy = tree, 
                        tip = tree$tip.label[!tree$tip.label%in%dfFeedingType$species_name], rooted = T)
dfFeedingType <- dfFeedingType[match(feedingTypeTree$tip.label, dfFeedingType$species_name), ]
dfFeedingType <- as.data.frame(dfFeedingType)
row.names(dfFeedingType) <- dfFeedingType$species_name
dfFeedingType <- dfFeedingType[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfFeedingType$FeedingType.x <- relevel(dfFeedingType$FeedingType.x, ref = "hunting macrofauna (predator)")
dfRegression$FeedingType <- relevel(dfRegression$FeedingType, ref = "hunting macrofauna (predator)")
# Try using caper.
c_data <- comparative.data(tree, dfFeedingType, "species_name")
caperFeedingType <- pgls(branch_length ~ FeedingType.x + number_of_nodes, c_data, lambda = "ML")
anovaFeedingType <- anova.pgls.fixed(caperFeedingType)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyFT <- model.matrix(~ FeedingType.x - 1, data=dfFeedingType)
row.names(dfDummyFT) <- dfFeedingType$species_name
dfDummyFT <- as.data.frame(dfDummyFT)
dfDummyFT$species_name <- row.names(dfDummyFT)
colnames(dfDummyFT)[1] <- "Predator"
colnames(dfDummyFT)[2] <- "BrowsingSubstrate"
colnames(dfDummyFT)[3] <- "FilteringPlankton"
colnames(dfDummyFT)[4] <- "GrazingPlants"
colnames(dfDummyFT)[5] <- "SelPlankton"
colnames(dfDummyFT)[6] <- "Variable"
sigPred <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = Predator)
sigBrowse <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = BrowsingSubstrate)
sigFilt <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = FilteringPlankton)
sigGraze <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = GrazingPlants)
sigSel <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = SelPlankton)
sigVar <- phylo.d(dfDummyFT, feedingTypeTree, names.col = species_name, binvar = Variable)

### Multivariate (multiple regression) analyses ###
# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfRegression[, c(1:2, 4, 46, 13:14, 5, 23, 7, 10, 16, 12, 32, 43)]
colnames(dfMultivariate)[2] <- "branch_length"
colnames(dfMultivariate)[3] <- "number_of_nodes"
dfMultivariate <- setDT(dfMultivariate)
# First, order the columns by the amount of missing data (NA values).
dfTraitsNA <- sort(dfMultivariate[, lapply(.SD, function(x) sum(is.na(x)))])
# Reorder the original dfTraits. The columns with the least amount of NA values are now first.
setcolorder(dfMultivariate, names(dfTraitsNA))
rm(dfTraitsNA)
# Now, count the number of NAs in each row.
dfMultivariate[, count := rowSums(is.na(dfMultivariate))]
# Sort the rows by this order.
dfMultivariate <- dfMultivariate[order(count), ]
# Remove this count column as it is no longer needed.
dfMultivariate[, count := NULL]
# The datatable is now sorted: columns with the least amount of NAs to the left, 
# rows with with the least amount of NAs on top.
# Now I want to loop through the traits, removing one column (trait) at a time
# and count the number of complete cases. This will provide us some information
# as to which traits would provide an adequate sample size for downstream analysis.
# First, take the number of columns in dfTemp.
len <- ncol(dfMultivariate)
# Create a temporary variable to hold this number. This variable will hold the 
# number of subsets to check at each iteration.
tempLen <- len
# Create a vector to hold the results of the loop.
all.cc <- NULL
# Start the loop:
for (i in 1:len) {
  # Works best if you set dfTemp back to a dataframe.
  x <- as.data.frame(dfMultivariate)
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
names(all.cc) <- rev(colnames(dfMultivariate))
# Look at the results.
all.cc
# rep_guild_2 seems like a good cutoff point.
which(colnames(dfMultivariate) == "age_at_maturity")
dfMultivariateCut <- dfMultivariate[, 1:8] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariateCut <- dfMultivariateCut %>% filter(complete.cases(.))

# Now for forward selection PGLS.
# Prune constrained tree so only those tips that are needed are present.
firstTree <- drop.tip(phy = tree, 
                          tip = tree$tip.label[!tree$tip.label%in%dfMultivariateCut$species_name], rooted = T)
# Make sure the order of the data matches the constrained tree.
dfMultivariateCut <- dfMultivariateCut[match(firstTree$tip.label, dfMultivariateCut$species_name), ]

# Now check for data variability in this subset.
hist(dfMultivariateCut$branch_length)
range(dfMultivariateCut$branch_length)

hist(dfMultivariateCut$number_of_nodes)
range(dfMultivariateCut$number_of_nodes)

hist(dfMultivariateCut$Length)
range(dfMultivariateCut$Length)

hist(dfMultivariateCut$DepthRangeDeep)
range(dfMultivariateCut$DepthRangeDeep)

hist(dfMultivariateCut$median_lat)
range(dfMultivariateCut$median_lat)

table(dfMultivariateCut$Saltwater)

table(dfMultivariateCut$BodyShapeI)

table(dfMultivariateCut$Lakes)

hist(dfMultivariateCut$age_at_maturity)
range(dfMultivariateCut$age_at_maturity)

lm.res=lm(branch_length ~ Length + DepthRangeDeep + median_lat
          + age_at_maturity, data = dfMultivariateCut)
library(car)
vif(mod=lm.res)  # Longevity removed.

# Backward selection using BIC.
c_data <- comparative.data(tree, dfMultivariateCut, names.col = "species_name", vcv=TRUE)

# Full model.
full <- pgls(branch_length ~ number_of_nodes + DepthRangeDeep
             + Length + median_lat + age_at_maturity, data=c_data, lambda="ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x=fitted(full), y=full$phyres, pch=5)

# Test for effect of bs.
bs <- pgls(branch_length ~ number_of_nodes + DepthRangeDeep
             + Length + median_lat, 
              data=c_data, lambda="ML")
anova(full, bs)
BIC(full, bs)
# Confounding? No.
summary(full)
summary(bs)

# Test for effect of body shape.
bs <- pgls(branch_length ~ number_of_nodes + DepthRangeDeep
           + Length + median_lat + Saltwater, 
               data=c_data, lambda="ML")
BIC(lakes, bs)
# Confounding? 
summary(lakes)$coefficients
summary(bs)$coefficients

null <- pgls(branch_length ~ number_of_nodes, data=c_data, lambda="ML")
anova(null, full, test="F")  # As a whole predictors are significant.
summary(full)$coefficients

# PLOTS!
# Plot branch length vs. median lat.
plot(branch_length ~ log(Length.x), dfMaxLength, pch = 19, 
     xlab = "Log(Maximum length (cm))", ylab = "Branch length")
lm(dfDepthDeep$branch_length ~ dfDepthDeep$DepthRangeDeep.x)
abline(a = 1.108e+00, b = -0.0000018, col = "red")
with(dfLakes, plot(branch_length ~ Lakes.x, xlab = "Presence/absence in lakes",
                        ylab = "Branch length", 
                        main = " "))
axis(2, cex.axis = 1)

ggplot(dfDepthDeep, aes(x = log(DepthRangeDeep.x), y=branch_length)) + geom_point(size=2, shape=19) + 
  geom_abline(intercept = 1.108, slope = 1.037e-05, size = 0.75, color = "red") + 
  labs(x="Log(Maximum depth (m))", y=expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

ggplot(dfRegression, aes(branch_length)) + 
  geom_histogram(col="black",fill="deepskyblue1", alpha = .2) + 
  labs(x="Branch length", y="Frequency") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

qplot(dfRegression$branch_length,
      geom="histogram", 
      binwidth = 0.01)

ggplot(dfCoralReefs, aes(x = CoralReefs.x, y = branch_length.y)) +
  geom_boxplot() + labs(x="Coral Reef Presence/Absence", y="Branch length") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

# Save file in current directory.
png(filename="depthWhole.png", 
    units="in", 
    width=7, 
    height=5, 
    pointsize=12, 
    res=300)
my_sc_plot(data)
dev.off()