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
                                 BodyShapeI, LongevityWild, Length, LTypeMaxM)]
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

# Finally, prepare the dfSpeciesGenMV datatable by merging all univariate 
# datatables. This datatable will be used for the eventual multivariate analysis.
# Merging multiple datatables at once.
dfSpeciesGenMV <- Reduce(function(...) merge(...), 
                         list(dfBodyShapeI, dfLongWild, dfMaxLength))

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
dfEcologyMV <- Reduce(function(...) merge(...), list(dfNeritic, dfOceanic,
                                                     dfStreams, dfLakes,
                                                     dfFeedingType, dfBenthic))

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
dfTraits <- merge(dfFilteredSingle, dfLatitude, all = TRUE, by = "bin_uri")
# Set the keys for datatable merging.
setkey(dfTraits, species_name)
setkey(dfSpeciesGenMV, species_name)
setkey(dfEcologyMV, species_name)
dfTraits <- Reduce(function(...) merge(..., all = T), list(dfTraits, 
                                                           dfSpeciesGenMV,
                                                           dfEcologyMV))
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:4, 7:16)]
colnames(dfTraits)[4] <- "filtered_bin_size"

# Merge back to dfFiltered to obtain all of the sequence information for 
# each BIN. This is for creation of the master phylogeny.
dfPreCentroid <- merge(dfFiltered, dfTraits, by = "bin_uri")
# Dataframe reorganization and renaming.
colnames(dfPreCentroid)[6] <- "species_name"
colnames(dfPreCentroid)[8] <- "initial_bin_size"
colnames(dfPreCentroid)[9] <- "filtered_bin_size"
dfPreCentroid <- dfPreCentroid[, c(1:9)]

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

### ALIGNMENT QUALITY CHECKING ###
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
# If you decide to remove all from your data:
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
# Remove 0 values.
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
# Remove them from dfAllSeqs, if desired.
dfAllSeqs <- dfAllSeqs[-outliers, ]

# LEFT OFF HERE!

### NEAREST NEIGHBOUR CHECK ###
# Remove centroid sequences whose nearest neighbours are in a different order or 
# family. The nearest neighbour can be determined from the distance matrix alone.
# It is the sequence with minimum pairwise distance to the sequence in question.

# FUNCTION: Identifies the nearest neighbour of each BIN.
# Applying this function to geneticDistanceCentroid.
NearestNeighbour <- t(sapply(seq(nrow(geneticDistanceCentroid)), function(i) {
  # Identify the sequence with the minimum distance in the row (the nearest
  # neighbour of i)
  j <- which.min(geneticDistanceCentroid[i, ])
  # Assign the name of the BIN and the name of its nearest neighbour to a 
  # separate column.
  c(paste(rownames(geneticDistanceCentroid)[i], 
          colnames(geneticDistanceCentroid)[j], sep='/'), 
    geneticDistanceCentroid[i, j])
}))
# Convert to dataframe.
NearestNeighbour <- as.data.frame(NearestNeighbour, stringsAsFactors = FALSE)
# Split the columns using the transform() function.
NearestNeighbour <- transform(NearestNeighbour,  V1 = colsplit(V1, split = "/", 
                                           names = c('bin_uri', 
                                                     'nearest_neighbour')))
# Re-add the genetic distances to the dataframe.
NearestNeighbour$V2 <- as.numeric(NearestNeighbour$V2)
colnames(NearestNeighbour)[2] <- "pairwise_distance"
# Get the family and order names of BINs and nearest neighbours.
dfNN <- NearestNeighbour$V1
dfNN$pairwise_distance <- NearestNeighbour$pairwise_distance
bins <- as.character(dfNN$bin_uri)
nn <- as.character(dfNN$nearest_neighbour)
# Now extract the orders and families from dfAllSeqs.
bin_ord <- dfAllSeqs[match(dfNN$bin_uri, dfCheckAllSeqs$bin_uri), ] 
bin_orders <- bin_ord$order_name
bin_families <- bin_ord$family_name
nn_ord <- dfAllSeqs[match(dfNN$nearest_neighbour, dfCheckAllSeqs$bin_uri), ] 
nn_orders <- nn_ord$order_name
nn_families <- nn_ord$family_name
# Add these columns to dfNN.
dfNN$bin_order <- bin_orders
dfNN$bin_family <- bin_families
dfNN$nn_order <- nn_orders
dfNN$nn_family <- nn_families
dfNN <- setDT(dfNN)
# Non-matching orders.
nonmatchOrd <- which(dfNN$bin_order != dfNN$nn_order)
# Look closely.
dfNonmatchOrd <- dfNN[nonmatchOrd, ]
dfNonmatchOrd <- dfNonmatchOrd[which(pairwise_distance < 0.05)]
# Remove these BINs.
nonMatchOrdBins <- dfNonmatchOrd$bin_uri
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% nonMatchOrdBins, ]
# Non-matching fams.
nonmatchFam <- which(dfNN$bin_family != dfNN$nn_family)
# LOook closely.
dfNonmatchFam <- dfNN[nonmatchFam, ]
dfNonmatchFam <- dfNonmatchFam[which(pairwise_distance < 0.05)]
# Remove these BINs.
nonMatchFamBins <- dfNonmatchFam$bin_uri
# If you decide to remove all from your data:
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% nonMatchFamBins, ]

# Now check for any neighbours under 0.05 in different ord or fam.
# Subset out all close neighbour pairings.
dfGeneticDistance <- as.data.frame(geneticDistanceCentroid)
closeNeighbours <- which(dfGeneticDistance < 0.05, arr.ind = TRUE)
bins <- row.names(dfGeneticDistance[closeNeighbours[, 1],] )
nn <- colnames(dfGeneticDistance[closeNeighbours[, 2]])
closeNeighbours <- as.data.frame(bins)
closeNeighbours$bins <- as.character(closeNeighbours$bins)
closeNeighbours$neighbour <- nn
# Get family and order names of BINs and nearest neighbours.
# Make dfs following these orders we can extract order and family names.
bin_ord <- dfAllSeqs[match(closeNeighbours$bins, dfAllSeqs$bin_uri), ] 
bin_orders <- bin_ord$order_name
bin_families <- bin_ord$family_name
nn_ord <- dfAllSeqs[match(closeNeighbours$neighbour, dfAllSeqs$bin_uri), ] 
nn_orders <- nn_ord$order_name
nn_families <- nn_ord$family_name
# Add these columns to the dfNN.
closeNeighbours$bin_order <- bin_orders
closeNeighbours$bin_family <- bin_families
closeNeighbours$nn_order <- nn_orders
closeNeighbours$nn_family <- nn_families
closeNeighbours <- setDT(closeNeighbours)
# Fill in missing data in bin columnn.
# BINs columns.
closeNeighbours <- closeNeighbours[order(closeNeighbours$bins), ]
#install.packages("zoo")
library(zoo)
closeNeighbours$bin_order <- na.locf(closeNeighbours$bin_order)
closeNeighbours$bin_family <- na.locf(closeNeighbours$bin_family)
# NN columns.
closeNeighbours <- closeNeighbours[order(closeNeighbours$neighbour), ]
closeNeighbours$nn_order <- na.locf(closeNeighbours$nn_order)
closeNeighbours$nn_family <- na.locf(closeNeighbours$nn_family)
# Non matching orders.
nonmatchOrd <- which(closeNeighbours$bin_order != closeNeighbours$nn_order)
# Look closely.
dfNonmatchOrd <- dfNN[nonmatchOrd, ]
dfNonmatchOrd <- dfNonmatchOrd[which(pairwise_distance < 0.05)]
# Non matching fams.
nonmatchFam <- which(closeNeighbours$bin_family != closeNeighbours$nn_family)
# LOook closely.
dfNonmatchFam <- closeNeighbours[nonmatchFam, ]  # Some ones I already removed!
# If you decide to remove all from your data:
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% dfNonmatchFam$bins, ]

# Align and trim dfAllSeqs again without the extreme BINs and conflicted BINs.
dfCheck2AllSeqs <- refSeqTrim(dfAllSeqs)

# If dfCheck2AllSeqs alignment acceptable, proceed!

# Which outgroups made the cut? Remove them from MSA so I can build tree
# just using the ingroup.
outgroupSpecies <- unique(dfOutgroup$species_name)
dfGoodOutgroups <- dfAllSeqs[dfAllSeqs$species_name %in% outgroupSpecies, ]
outgroupBins <- unique(dfGoodOutgroups$species_name)
# Remove the outgroups from dfCheck2AllSeqs.
dfCheck3AllSeqs <- dfCheck2AllSeqs[!dfCheck2AllSeqs$bin_uri %in% outgroupBins, ]
# Run the alignment without the outgroups.
# Check over sequences/alignment, make sure it is in correct reading frame.
dfCheck4AllSeqs <- refSeqTrim(dfCheck3AllSeqs)

# Now run with the outgroups.
goodOG <- which(dfAllSeqs$species_name == "Raja montagui" |
                dfAllSeqs$species_name == "Raja polystigma" |
                dfAllSeqs$species_name == "Neoceratodus forsteri" |
                dfAllSeqs$species_name == "Castor fiber" |
                dfAllSeqs$species_name == "Cratogeomys perotensis")
mergeOG <- dfAllSeqs[goodOG, ]
# Add them back.
dfCheck5AllSeqs <- rbind(dfCheck4AllSeqs, mergeOG)
# Run the alignment.
dfFinalAllSeqs <- refSeqTrim(dfCheck5AllSeqs)

################################################################################
# CHECK THAT ALIGNMENT AND dfCheck2AllSeqs are the same.
x <- read.fasta(file="withNOoutgroupsMAY23rd.fas")
fish <- names(x)
fish <- gsub("_", " ", fish)

dfNamesNoOG <- dfCheckAllSeqs[dfCheckAllSeqs$species_name %in% fish, ]

# TREE SECTION! #
# BINARY CONSTRAINT TREE PRUNING #
# First, let's download the tree we are using as a topological constraint.
fishTree <- read.tree(file = "fishTree.tre")
# Change the tip labels to just species names.
fishTree$tip.label <- gsub("_", " ", fishTree$tip.label)
fishTree$tip.label <- substr(fishTree$tip.label, 7, nchar(fishTree$tip.label))
fishTree$node.label <- NULL

# Prune the constraint tree so only those tips that are match with names in 
# dfNoOutgroups remain.
prunedFishTree <- drop.tip(phy = fishTree, 
                           tip = fishTree$tip.label[!fishTree$tip.label%in%dfCheck4AllSeqs$species_name], 
                           rooted = T)
prunedFishTree$edge.length <- NULL  # Don't need branch lengths for binary constraint file (just relationships.)
prunedFishTree$node.label <- NULL
write.tree(prunedFishTree, file = "prunedFishTreeFINAL.txt")
write.tree(prunedFishTree, file = "prunedFishTreeFINAL.tre")


### SECTION 4: UNIVARIATE ANALYSES ###
# This is the section where the univariate analyses are performed as a form
# of model selection.

### WHOLE ALIGNMENT: CONSTRAINED ###
# Read in the tree.
rootedWholeTree <- read.tree(file = "RAxML_labelledTree.FinalTryEPA")
rootedWholeTree$tip.label <- gsub("_", " ", rootedWholeTree$tip.label)

# Root the tree.
outgroups <- c("QUERY   Raja montagui", "QUERY   Raja polystigma", 
               "Echinorhinus cookei", "Echinorhinus brucus")
rootedWholeTree <- root(rootedWholeTree, outgroup = outgroups, 
                        resolve.root = TRUE)

################################################################################
# TRAIT: NUMBER OF NODES
# Let's first determine the number of nodes for each species.
#phy4ML <- as(unconstrainedTree, "phylo4")
#plot(phy4ML, show.node = TRUE)
#root <- rootNode(phy4ML)
# A bit slow.
#nodeList <- lapply(1:nTips(phy4ML), function(i) .tipToRoot(phy4ML, i, root))
#numberOfNodes2 <- sapply(1:nTips(phy4ML), function(i) length(nodeList[[i]]))
#names(numberOfNodes2) <- tipLabels(phy4ML)
numberOfNodes <- distRoot(rootedWholeTree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(numberOfNodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)

################################################################################
# Merge with the node_number df.
# Match the species names in data to tip labels.
dfAllSeqsNode <- merge(dfAllSeqs, dfNumberOfNodes, by = "species_name", 
                       all = TRUE)

# Let's calculate the sum of branch lengths now (from root to tip).
#is.rooted(unconstrainedTree)
#branchLengths <- distRoot(unconstrainedTree, tips = "all", method = "patristic")
branchLengths <- diag(vcv.phylo(rootedWholeTree)) 
# Make into a dataframe.
dfBranchLengths <- data.frame(branchLengths)
colnames(dfBranchLengths)[1] <- "branchLength"
dfBranchLengths$species_name <- names(branchLengths)
hist(dfBranchLengths$branchLength, main = "", xlab = "Branch Length")
median(dfBranchLengths$branchLength)
mean(dfBranchLengths$branchLength)
range(dfBranchLengths$branchLength)
# Merge dfBranchLengths and dfRecodedPreCent.
dfRegression <- merge(dfBranchLengths, dfAllSeqsNode, by = "species_name")
# Merge back to the trait dataframe.
dfRegression <- merge(dfRegression, dfTraits, all.x = TRUE, by = "bin_uri")
# Some reordering.
dfRegression <- dfRegression[, c(2:4, 11, 15:56)]
colnames(dfRegression)[1] <- "species_name"
colnames(dfRegression)[3] <- "bin_size"


### SECTION 5: STATISTICS ###
# Check the traits for inclusion in MV analysis.
# First, make sure the trait data and phylo tree are in the same order.
# Make sure the order of the data matches the tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, 
                                   dfRegression$species_name), ]

# TRAIT: Branch length.
# Pagel's lambda. A test for phylogenetic signal.
branch_length <- dfRegression$branchLength
names(branch_length) <- rootedWholeTree$tip.label
sigBL <- phylosig(rootedWholeTree, branch_length, method = "lambda", 
                  test = TRUE)

# TRAIT: Number of nodes.
# Pagel's lambda. A test for phylogenetic signal.
number_of_nodes <- dfRegression$numberOfNodes
names(number_of_nodes) <- rootedWholeTree$tip.label
sigNodes <- phylosig(rootedWholeTree, number_of_nodes, method = "lambda", 
                     test = TRUE)
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, 
                                   dfRegression$species_name), ]
dfRegression <- as.data.frame(dfRegression)
dfNodes <- dfRegression[, c(1:2, 4)]
dfNodes <- dfNodes[sample(nrow(dfNodes), 4000), ] # 5929
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfNodes, "species_name")
caperNodes <- pgls(branchLength ~ numberOfNodes, c_data, lambda = "ML")

# TRAIT: Median latitude.
# Pagel's lambda. A test for phylogenetic signal.
median_lat <- dfRegression$median_lat
names(median_lat) <- rootedWholeTree$tip.label
sigLat <- phylosig(rootedWholeTree, median_lat, method = "lambda", 
                   test = TRUE)
# Pruning tree for latitude.
# Merge so can get branch length info for UV.
dfLatitude <- merge(dfLatitude, dfRegression, by.x = "species_label", 
                    by.y = "species_name")
# Take BINs with highest number of sequences with species info in order to avoid
# using BINs that were assigned the same species label.
dfLatitude <- merge(aggregate(bin_size ~ species_label, data = dfLatitude, max), 
                    dfLatitude, all.x = T, sort = TRUE)
# Dealing with ties.
dup_majority_species <- which(duplicated(dfLatitude$species_label))
dfLatitude <- dfLatitude[-dup_majority_species,]
rm(dup_majority_species)
# Prune constrained tree so only those tips for lat are present.
latTree <- drop.tip(phy = rootedWholeTree, 
                    tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfLatitude$species_name], 
                    rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfLatitude <- dfLatitude[match(latTree$tip.label, dfLatitude$species_name), ]
dfLatitude <- as.data.frame(dfLatitude)
row.names(dfLatitude) <- dfLatitude$species_name
dfLatitude <- dfLatitude[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLatitude, "species_name")
caperLatitude <- pgls(branchLength ~ median_lat.x + numberOfNodes, c_data, 
                      lambda = "ML")

# TRAIT: Maximum length.
# Pagel's lambda. A test for phylogenetic signal.
maximum_length <- dfRegression$Length
names(maximum_length) <- rootedWholeTree$tip.label
sigLength <- phylosig(rootedWholeTree, maximum_length, method = "lambda", 
                      test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfMaxLength <- merge(dfMaxLength, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
lengthTree <- drop.tip(phy = rootedWholeTree, 
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMaxLength$species_name], 
                       rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfMaxLength <- dfMaxLength[match(lengthTree$tip.label, 
                                 dfMaxLength$species_name), ]
dfMaxLength <- as.data.frame(dfMaxLength)
row.names(dfMaxLength) <- dfMaxLength$species_name
dfMaxLength <- dfMaxLength[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfMaxLength, "species_name")
caperLength <- pgls(branchLength ~ Length.x + numberOfNodes, c_data, 
                    lambda = "ML")

# TRAIT: Longevity.
# Pagel's lambda. A test for phylogenetic signal.
longevity <- dfRegression$LongevityWild
names(longevity) <- rootedWholeTree$tip.label
sigLong <- phylosig(rootedWholeTree, longevity, method = "lambda", test = TRUE)
# Pruning tree for latitude.
# Merge so can get branch length info for UV.
dfLongWild <- merge(dfLongWild, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
longevityTree <- drop.tip(phy = rootedWholeTree, 
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfLongWild$species_name], 
                       rooted = T)
tiph <- diag(vcv.phylo(longevityTree))
dfLongWild <- dfLongWild[match(longevityTree$tip.label, dfLongWild$species_name), ]
dfLongWild <- as.data.frame(dfLongWild)
row.names(dfLongWild) <- dfLongWild$species_name
dfLongWild <- dfLongWild[, c(1:3, 5)]
dfLongWild <- na.omit(dfLongWild)
fitLong <- gls(branchLength ~ LongevityWild.x + numberOfNodes, 
               correlation=corPagel(value = 0, phy = longevityTree),
           weights=varFixed(~tiph), method = "ML", data = dfLongWild)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLongWild, "species_name")
caperLong <- pgls(branchLength ~ LongevityWild.x + numberOfNodes, c_data, 
                  lambda = "ML")

# TRAIT: Age at maturity
# Pagel's lambda. A test for phylogenetic signal.
age_at_mat <- dfRegression$age_at_maturity
names(age_at_mat) <- rootedWholeTree$tip.label
sigAge <- phylosig(rootedWholeTree, age_at_mat, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfAgeMaturity <- merge(dfAgeMaturity, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
ageTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfAgeMaturity$species_name], 
                    rooted = T)
tiph <- diag(vcv.phylo(ageTree))
dfAgeMaturity <- dfAgeMaturity[match(ageTree$tip.label, dfAgeMaturity$species_name), ]
dfAgeMaturity <- as.data.frame(dfAgeMaturity)
row.names(dfAgeMaturity) <- dfAgeMaturity$species_name
dfAgeMaturity <- dfAgeMaturity[, c(1:3, 5)]
fitAge <- gls(branchLength ~ age_at_maturity.x + numberOfNodes, 
              correlation=corPagel(value = 0, phy = ageTree), 
              weights=varFixed(~tiph), method = "ML", data = dfAgeMaturity)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfAgeMaturity, "species_name")
caperAge <- pgls(branchLength ~ age_at_maturity.x + numberOfNodes, c_data, 
                 lambda = "ML")

# Discrete traits.
# TRAIT: NERITIC
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, 
                                   dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfNeritic <- merge(dfNeritic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
neriticTree <- drop.tip(phy = rootedWholeTree, 
                           tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfNeritic$species_name], 
                        rooted = T)
dfNeritic <- dfNeritic[match(neriticTree$tip.label, dfNeritic$species_name), ]
dfNeritic <- as.data.frame(dfNeritic)
row.names(dfNeritic) <- dfNeritic$species_name
dfNeritic <- dfNeritic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfNeritic$Neritic.x <- relevel(dfNeritic$Neritic.x, ref = "0")
dfRegression$Neritic <- relevel(dfRegression$Neritic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfNeritic, "species_name")
caperNeritic <- pgls(branchLength ~ Neritic.x + numberOfNodes, c_data, 
                     lambda = "ML")
# D metric for phylogenetic signal.
sigNer <- phylo.d(dfNeritic, neriticTree, names.col = species_name, 
                  binvar = Neritic.x)
#fit$opt$lambda

# TRAIT: OCEANIC
#oceanic <- dfRegression$Oceanic
#names(oceanic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, oceanic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, 
                                   dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfOceanic <- merge(dfOceanic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
oceanicTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfOceanic$species_name], 
                        rooted = T)
dfOceanic <- dfOceanic[match(oceanicTree$tip.label, dfOceanic$species_name), ]
dfOceanic <- as.data.frame(dfOceanic)
row.names(dfOceanic) <- dfOceanic$species_name
dfOceanic <- dfOceanic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfOceanic$Oceanic.x <- relevel(dfOceanic$Oceanic.x, ref = "0")
dfRegression$Oceanic <- relevel(dfRegression$Oceanic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfOceanic, "species_name")
caperOceanic <- pgls(branchLength ~ Oceanic.x + numberOfNodes, c_data, 
                     lambda = "ML")
# D metric for phylogenetic signal.
sigOceanic <- phylo.d(dfOceanic, oceanicTree, names.col = species_name, 
                      binvar = Oceanic.x)

# TRAIT: BENTHIC
#benthic <- dfRegression$Benthic
#names(benthic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, benthic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, 
                                   dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBenthic <- merge(dfBenthic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
benthicTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfBenthic$species_name], 
                        rooted = T)
dfBenthic <- dfBenthic[match(benthicTree$tip.label, dfBenthic$species_name), ]
dfBenthic <- as.data.frame(dfBenthic)
row.names(dfBenthic) <- dfBenthic$species_name
dfBenthic <- dfBenthic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBenthic$Benthic.x <- relevel(dfBenthic$Benthic.x, ref = "0")
dfRegression$Benthic <- relevel(dfRegression$Benthic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfBenthic, "species_name")
caperBenthic <- pgls(branchLength ~ Benthic.x + numberOfNodes, c_data, 
                     lambda = "ML")
# D metric for phylogenetic signal.
sigBenthic <- phylo.d(dfBenthic, benthicTree, names.col = species_name, 
                      binvar = Benthic.x)

# TRAIT: LAKES
#lakes <- dfRegression$Lakes
#names(lakes) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, lakes, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfLakes <- merge(dfLakes, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
lakeTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfLakes$species_name], rooted = T)
dfLakes <- dfLakes[match(lakeTree$tip.label, dfLakes$species_name), ]
dfLakes <- as.data.frame(dfLakes)
row.names(dfLakes) <- dfLakes$species_name
dfLakes <- dfLakes[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfLakes$Lakes.x <- relevel(dfLakes$Lakes.x, ref = "0")
dfRegression$Lakes <- relevel(dfRegression$Lakes, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfLakes, "species_name")
caperLakes <- pgls(branchLength ~ Lakes.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigLakes <- phylo.d(dfLakes, lakeTree, names.col = species_name, binvar = Lakes.x)

# TRAIT: STREAMS
#streams <- dfRegression$Stream
#names(streams) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, streams, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfStreams <- merge(dfStreams, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
streamTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfStreams$species_name], rooted = T)
dfStreams <- dfStreams[match(streamTree$tip.label, dfStreams$species_name), ]
dfStreams <- as.data.frame(dfStreams)
row.names(dfStreams) <- dfStreams$species_name
dfStreams <- dfStreams[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfStreams$Stream.x <- relevel(dfStreams$Stream.x, ref = "0")
dfRegression$Stream <- relevel(dfRegression$Stream, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfStreams, "species_name")
caperStreams <- pgls(branchLength ~ Stream.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigStreams <- phylo.d(dfStreams, streamTree, names.col = species_name, binvar = Stream.x)

# TRAIT: BODY SHAPE I
#body_shape_I <- dfRegression$BodyShapeI
#names(body_shape_I) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, body_shape_I, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBodyShapeI <- merge(dfBodyShapeI, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
bodyTree <- drop.tip(phy = rootedWholeTree, 
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfBodyShapeI$species_name], rooted = T)
dfBodyShapeI <- dfBodyShapeI[match(bodyTree$tip.label, dfBodyShapeI$species_name), ]
dfBodyShapeI <- as.data.frame(dfBodyShapeI)
row.names(dfBodyShapeI) <- dfBodyShapeI$species_name
dfBodyShapeI <- dfBodyShapeI[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBodyShapeI$BodyShapeI.x <- relevel(dfBodyShapeI$BodyShapeI.x, ref = "fusiform / normal")
dfRegression$BodyShapeI <- relevel(dfRegression$BodyShapeI, ref = "fusiform / normal")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfBodyShapeI, "species_name")
caperBodyShapeI <- pgls(branchLength ~ BodyShapeI.x + numberOfNodes, c_data, lambda = "ML")
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

# FeedingType
#feeding_type <- dfRegression$FeedingType
#names(feeding_type) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, feeding_type, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfFeedingType <- merge(dfFeedingType, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
feedingTypeTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfFeedingType$species_name], rooted = T)
dfFeedingType <- dfFeedingType[match(feedingTypeTree$tip.label, dfFeedingType$species_name), ]
dfFeedingType <- as.data.frame(dfFeedingType)
row.names(dfFeedingType) <- dfFeedingType$species_name
dfFeedingType <- dfFeedingType[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfFeedingType$FeedingType.x <- relevel(dfFeedingType$FeedingType.x, ref = "hunting macrofauna (predator)")
dfRegression$FeedingType <- relevel(dfRegression$FeedingType, ref = "hunting macrofauna (predator)")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfFeedingType, "species_name")
caperFeedingType <- pgls(branchLength ~ FeedingType.x + numberOfNodes, c_data, lambda = "ML")
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

# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfRegression[, c(1:2, 4, 46, 13:14, 5, 23, 7, 10, 16, 12, 32, 43)]
colnames(dfMultivariate)[2] <- "branchLength"
colnames(dfMultivariate)[3] <- "numberOfNodes"
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
firstTree <- drop.tip(phy = rootedWholeTree, 
                          tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMultivariateCut$species_name], rooted = T)
# Make sure the order of the data matches the constrained tree.
dfMultivariateCut <- dfMultivariateCut[match(firstTree$tip.label, dfMultivariateCut$species_name), ]

# Now check for data variability in this subset.
hist(dfMultivariateCut$branchLength)
range(dfMultivariateCut$branchLength)

hist(dfMultivariateCut$numberOfNodes)
range(dfMultivariateCut$numberOfNodes)

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

lm.res=lm(branchLength ~ Length + DepthRangeDeep + median_lat
          + age_at_maturity, data = dfMultivariateCut)
library(car)
vif(mod=lm.res)  # Longevity removed.

# Backward selection using BIC.
c_data <- comparative.data(rootedWholeTree, dfMultivariateCut, names.col = "species_name", vcv=TRUE)

# Full model.
full <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep
             + Length + median_lat + age_at_maturity, data=c_data, lambda="ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x=fitted(full), y=full$phyres, pch=5)

# Test for effect of bs.
bs <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep
             + Length + median_lat, 
              data=c_data, lambda="ML")
anova(full, bs)
BIC(full, bs)
# Confounding? No.
summary(full)
summary(bs)

# Test for effect of body shape.
bs <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep
           + Length + median_lat + Saltwater, 
               data=c_data, lambda="ML")
BIC(lakes, bs)
# Confounding? 
summary(lakes)$coefficients
summary(bs)$coefficients

null <- pgls(branchLength ~ numberOfNodes, data=c_data, lambda="ML")
anova(null, full, test="F")  # As a whole predictors are significant.
summary(full)$coefficients

### THIRD CODON POSITION: CONSTRAINED ###
# Read in the tree.
rootedThirdTree <- read.tree(file = "RAxML_bestTree.Third.PARTITION.1")
rootedThirdTree$tip.label <- gsub("_", " ", rootedThirdTree$tip.label)
# Root the tree.
outgroups <- c("Echinorhinus cookei", "Echinorhinus brucus")
rootedThirdTree <- root(rootedThirdTree, outgroup = outgroups, resolve.root = TRUE)

################################################################################
# TRAIT: NUMBER OF NODES
# Let's first determine the number of nodes for each species.
numberOfNodes <- distRoot(rootedThirdTree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(numberOfNodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)

################################################################################
# Merge with the node_number df.
dfAllSeqsNode <- merge(dfAllSeqs, dfNumberOfNodes, by = "species_name", all = TRUE)

# Let's calculate the sum of branch lengths now (from root to tip).
#is.rooted(rootedThirdTree)
#branchLengths <- distRoot(rootedThirdTree, tips = "all", method = "patristic")  # Takes a while.
branchLengths <- diag(vcv.phylo(rootedThirdTree))  # Much faster.
# Make into a dataframe.
dfBranchLengths <- data.frame(branchLengths)
colnames(dfBranchLengths)[1] <- "branchLength"
dfBranchLengths$species_name <- names(branchLengths)
hist(dfBranchLengths$branchLength, main = "", xlab = "Branch Length")
median(dfBranchLengths$branchLength)
mean(dfBranchLengths$branchLength)
range(dfBranchLengths$branchLength)
# Merge dfBranchLengths and dfRecodedPreCent.
dfThirdRegression <- merge(dfBranchLengths, dfAllSeqsNode, by = "species_name")
# Merge back to the trait dataframe.
dfThirdRegression <- merge(dfThirdRegression, dfTraits, all.x = TRUE, by = "bin_uri")
# Some reordering.
dfThirdRegression <- dfThirdRegression[, c(2:4, 11, 15:56)]
colnames(dfThirdRegression)[1] <- "species_name"
colnames(dfThirdRegression)[3] <- "bin_size"


### SECTION 5: STATISTICS ###
# Check the traits for inclusion in MV analysis.
# First, make sure the trait data and phylo tree are in the same order.
# Make sure the order of the data matches the tree.
dfThirdRegression <- dfThirdRegression[match(rootedWholeTree$tip.label, dfThirdRegression$species_name), ]

# TRAIT: Branch length.
# Pagel's lambda. A test for phylogenetic signal.
branch_length <- dfThirdRegression$branchLength
names(branch_length) <- rootedThirdTree$tip.label
sigBL3 <- phylosig(rootedThirdTree, branch_length, method = "lambda", test = TRUE)

# TRAIT: Number of nodes.
dfThirdRegression <- as.data.frame(dfThirdRegression)
dfNodes <- dfThirdRegression[, c(1:2, 4)]
dfNodes <- as.data.frame(dfNodes)
dfNodes <- dfNodes[!duplicated(row.names(dfNodes)), ]
# Test using caper.
row.names(dfNodes) <- dfNodes$species_name
c_data <- comparative.data(rootedWholeTree, dfNodes, "species_name")
caperNodes <- pgls(branchLength ~ numberOfNodes, c_data, lambda = "ML")

# TRAIT: Median latitude.
dfLatitude <- merge(dfLatitude, dfThirdRegression, by = "species_name")
dfLatitude <- dfLatitude[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLatitude, "species_name")
caperLatitude <- pgls(branchLength.y ~ median_lat.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: Maximum length.
dfMaxLength <- merge(dfMaxLength, dfThirdRegression, by = "species_name")
dfMaxLength <- dfMaxLength[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfMaxLength, "species_name")
caperLength <- pgls(branchLength.y ~ Length.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: Longevity.
dfLongWild <- merge(dfLongWild, dfThirdRegression, by = "species_name")
dfLongWild <- dfLongWild[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLongWild, "species_name")
caperLong <- pgls(branchLength.y ~ LongevityWild.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: Maximum depth
# Merge so can get branch length info for third codon.
dfDepthDeep <- merge(dfDepthDeep, dfThirdRegression, by = "species_name")
dfDepthDeep <- dfDepthDeep[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDepthDeep, "species_name")
caperDeep <- pgls(branchLength.y ~ DepthRangeDeep.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: Age at maturity
# Merge so can get branch length info for third codon.
dfAgeMaturity <- merge(dfAgeMaturity, dfThirdRegression, by = "species_name")
dfAgeMaturity <- dfAgeMaturity[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfAgeMaturity, "species_name")
caperAge <- pgls(branchLength.y ~ age_at_maturity.x + numberOfNodes.y, c_data, lambda = "ML")

# Discrete traits.
# TRAIT: NERITIC
# Merge so can get branch length info for third codon.
dfNeritic <- merge(dfNeritic, dfThirdRegression, by = "species_name")
dfNeritic <- dfNeritic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfNeritic, "species_name")
caperNer <- pgls(branchLength.y ~ Neritic.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: OCEANIC
# Merge so can get branch length info for third codon.
dfOceanic <- merge(dfOceanic, dfThirdRegression, by = "species_name")
dfOceanic <- dfOceanic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfOceanic, "species_name")
caperOceanic <- pgls(branchLength.y ~ Oceanic.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: BENTHIC
# Merge so can get branch length info for third codon.
dfBenthic <- merge(dfBenthic, dfThirdRegression, by = "species_name")
dfBenthic <- dfBenthic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfBenthic, "species_name")
caperBenthic <- pgls(branchLength.y ~ Benthic.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: LAKES
# Merge so can get branch length info for third codon.
dfLakes <- merge(dfLakes, dfThirdRegression, by = "species_name")
dfLakes <- dfLakes[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLakes, "species_name")
caperLakes <- pgls(branchLength ~ Lakes.x + numberOfNodes, c_data, lambda = "ML")

# TRAIT: STREAMS
# Merge so can get branch length info for third codon.
dfStreams <- merge(dfStreams, dfThirdRegression, by = "species_name")
dfStreams <- dfStreams[, c(1:2, 5, 7)]
dfStreams <- as.data.frame(dfStreams)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfStreams, "species_name")
caperStreams <- pgls(branchLength.y ~ Stream.x + numberOfNodes.y, c_data, lambda = "ML")

# TRAIT: BODY SHAPE I
# Merge so can get branch length info for third codon.
dfBodyShapeI <- merge(dfBodyShapeI, dfThirdRegression, by = "species_name")
dfBodyShapeI <- dfBodyShapeI[, c(1:3, 7)]
dfBodyShapeI <- as.data.frame(dfBodyShapeI)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfBodyShapeI, "species_name")
caperBodyShapeI <- pgls(branchLength.y ~ BodyShapeI.x + numberOfNodes.y, c_data, lambda = "ML")
anovaBS <- anova.pgls.fixed(caperBodyShapeI)

# FeedingType
# Merge so can get branch length info for third codon.
dfFeedingType <- merge(dfFeedingType, dfThirdRegression, by = "species_name")
dfFeedingType <- dfFeedingType[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfFeedingType, "species_name")
caperFeedingType <- pgls(branchLength.y ~ FeedingType.x + numberOfNodes.y, c_data, lambda = "ML")
anovaFT <- anova.pgls.fixed(caperFeedingType)

# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfThird[, c(1:2, 46, 4:5, 14, 23)]
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
dfMultivariateCut <- dfMultivariate[, 1:7] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariateCut <- dfMultivariateCut %>% filter(complete.cases(.))

# Now for forward selection PGLS.
# Prune constrained tree so only those tips that are needed are present.
firstTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMultivariateCut$species_name], rooted = T)
# Make sure the order of the data matches the constrained tree.
dfMultivariateCut <- dfMultivariateCut[match(firstTree$tip.label, dfMultivariateCut$species_name), ]

# Now check for data variability in this subset.
hist(dfMultivariateCut$branchLength)
range(dfMultivariateCut$branchLength)

hist(dfMultivariateCut$numberOfNodes)
range(dfMultivariateCut$numberOfNodes)

hist(dfMultivariateCut$Length)
range(dfMultivariateCut$Length)

hist(dfMultivariateCut$DepthRangeDeep)
range(dfMultivariateCut$DepthRangeDeep)

table(dfMultivariateCut$Lakes)

table(dfMultivariateCut$BodyShapeI)

hist(dfMultivariateCut$age_at_maturity)
range(dfMultivariateCut$age_at_maturity)

hist(dfMultivariateCut$LongevityWild)
range(dfMultivariateCut$LongevityWild)

lm.res=lm(branchLength ~ median_lat + Length +
            CoralReefs + DepthRangeDeep +
            age_at_maturity + LongevityWild, data = dfMultivariateCut)
library(car)
vif(mod=lm.res)  # Longevity removed.

# Relevel traits.
dfMultivariateCut$CoralReefs <- relevel(dfMultivariateCut$CoralReefs, ref = "0")
dfMultivariateCut$Mangroves <- relevel(dfMultivariateCut$Mangroves, ref = "0")

# Backward selection using BIC.
c_data <- comparative.data(rootedWholeTree, dfMultivariateCut, names.col = "species_name", vcv=TRUE)

# Full model.
full <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep + median_lat + 
               Length + age_at_maturity, data=c_data, lambda = "ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x=fitted(full), y=full$phyres, pch=5)

# Test for effect of mangroves.
mg <- pgls(branchLength ~ numberOfNodes + median_lat + Length +
               CoralReefs + DepthRangeDeep + BodyShapeI, 
           data=c_data, lambda="ML")
BIC(full, mg)
summary(full)
summary(mg)
# No effect so remove this trait.

# Test for effect of body shape.
bs <- pgls(branchLength ~ numberOfNodes + median_lat + Length +
               CoralReefs + DepthRangeDeep, 
             data=c_data, lambda="ML")
BIC(mg, bs)
summary(mg)
summary(bg)
# No effect so remove this trait.

null <- pgls(branchLength ~ numberOfNodes, data=c_data, lambda="ML")
anova(null, full, test="F")  # As a whole predictors are significant.

# PLOTS!
# Plot branch length vs. median lat.
plot(branchLength ~ log(Length.x), dfMaxLength, pch = 19, 
     xlab = "Log(Maximum length (cm))", ylab = "Branch length")
lm(dfDepthDeep$branchLength ~ dfDepthDeep$DepthRangeDeep.x)
abline(a = 1.108e+00, b = -0.0000018, col = "red")
with(dfLakes, plot(branchLength ~ Lakes.x, xlab = "Presence/absence in lakes",
                        ylab = "Branch length", 
                        main = " "))
axis(2, cex.axis = 1)

ggplot(dfDepthDeep, aes(x = log(DepthRangeDeep.x), y=branchLength)) + geom_point(size=2, shape=19) + 
  geom_abline(intercept = 1.108, slope = 1.037e-05, size = 0.75, color = "red") + 
  labs(x="Log(Maximum depth (m))", y=expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

ggplot(dfRegression, aes(branchLength)) + 
  geom_histogram(col="black",fill="deepskyblue1", alpha = .2) + 
  labs(x="Branch length", y="Frequency") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

qplot(dfRegression$branchLength,
      geom="histogram", 
      binwidth = 0.01)

ggplot(dfCoralReefs, aes(x = CoralReefs.x, y = branchLength.y)) +
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


# FUNCTION for factors with multiple levels.
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