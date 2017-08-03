################################################################################

# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
#                  determination/reference sequence trimming (lines TBD).

################################################################################

##### SECTION 1: DATA PROCESSING #####

### PACKAGES REQUIRED ###

# For BOLD data:
#install.packages("bold")
library(bold)

# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages("dplyr")
library(dplyr)
#install.packages("foreach")
library(foreach)

# Load the function(s) that I designed for this script:
source("ResolveBIN.R")

################################################################################

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
# Modifying BIN column slightly to remove "BOLD:".
dfFiltered[, bin_uri := substr(bin_uri, 6, 13)]

### FILTER 2 ###
# Filtering for presence of a sequence.
containSeq <- dfFiltered[, grep("[CTG]", nucleotides)]
dfFiltered <- dfFiltered[containSeq, ]

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

### FILTER 6 ###
# Filter out sequences that are less than 640 bp and greater than 1000 bp. 
# This is because extremely long or short sequence lengths can interfere with 
# multiple sequence alignments.
# First, determine the lengths of the sequences without gaps.
seqLengths <- dfFiltered[, nchar(gsub("-", "", nucleotides))] 
# Which sequences are greater than 1000 bp and less than 640 bp in length?
seqLengthCheck <- which(seqLengths > 1000 | seqLengths < 640)
dfFiltered <- dfFiltered[-seqLengthCheck, ]

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
# This is to ensure that empty cells are not counted as their own taxa.
dfResolve[dfResolve == ""] <- NA

# Now find the number of orders/families/genera/species in each BIN.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), 
          keyby = bin_uri]
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), 
          keyby = bin_uri]

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
length(orderConflicts)

# Family level resolving.
# Checking these manually.
# Make sure to update number_of_families column first.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), 
          keyby = bin_uri]
# Find number of BINs with family level conflicts.
familyConflicts <- dfResolve[, which(number_of_families > 1), by = bin_uri]
familyConflicts <- unique(familyConflicts$bin_uri)
length(familyConflicts)
# Now use the ResolveBIN function to resolve these bin conflicts.
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
length(familyConflicts)

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
dfSpeciesLabel <- dfSpecies[, .N, by = .(bin_uri, species_name)][order(-N), .(species_label = species_name[1L]), keyby = bin_uri]
# Genus label.
containGenus <- dfResolve[, grep("[A-Z]", genus_name)]
dfGenus <- dfResolve[containGenus, ]
dfGenusLabel <- dfGenus[, .N, by = .(bin_uri, genus_name)][order(-N), .(genus_label = genus_name[1L]), keyby = bin_uri]
# Family label. Technically, there should be only one family/order name per BIN
# after manually resolving these cases.
# But doing this step just in case!!!
containFamily <- dfResolve[, grep("[A-Z]", family_name)]
dfFamily <- dfResolve[containFamily, ]
dfFamilyLabel <- dfFamily[, .N, by = .(bin_uri, family_name)][order(-N), .(family_label = family_name[1L]), keyby = bin_uri]
# Order label.
containOrder <- dfResolve[, grep("[A-Z]", order_name)]
dfOrder <- dfResolve[containOrder, ]
dfOrderLabel <- dfOrder[, .N, by = .(bin_uri, order_name)][order(-N), .(order_label = order_name[1L]), keyby = bin_uri]

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