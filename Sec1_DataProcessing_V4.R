# Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Matthew Orton (https://github.com/m-orton/R-Scripts) for design/testing/contributions to the sequence filters (lines 79-88, 132-153).
# Dr. Robert Hanner for recommendations about how to deal with BIN data.
# Adapted lines 295-297 from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/403310/matt-dowle.
# https://stackoverflow.com/questions/13273833/merging-multiple-data-table.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

##################################################################################################################

##### SECTION 1: DATA PROCESSING #####
# This section is primarily for quality control purposes of DNA barcode data obtained from the BOLD API. 
# Filter arguments may be altered to meet the user's needs.

### PACKAGES REQUIRED ###
# For BOLD data:
#install.packages("bold")
library(bold)
# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages("foreach")
library(foreach)
# Load the function(s) designed for this script:
source("ResolveBIN.R")
source("CountConflicts.R")
source("AssignLabel.R")

##################################################################################################################

# Download sequences from BOLD using the function bold_seqspec() for sequence and specimen data. 
# In addition, I am only selecting those columns needed for downstream analysis.
# Enter your taxon between the "".
dfRawSeqs <- bold_seqspec(taxon = "", geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", "genus_name", "species_name", "lat", "nucleotides", "markercode")]

# Download outgroup species data from BOLD. These sequences may be used to root phylogenetic trees (depending if the taxa are an appropriate outgroup for
# the organismal group under study).
# Enter your outgroup name between the "".
outgroups <- c("") 
dfOutgroup <- bold_seqspec(taxon = outgroups, geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", "genus_name", "species_name", "lat", "nucleotides", "markercode")]
# Combine dfOutgroup and dfRawSeqs so that they are in one useable dataframe. Also, convert to datatable as datatables have useful features for data manipulation.
dfFiltered <- as.data.table(rbind(dfRawSeqs, dfOutgroup))

### FILTER 1 ###
# Filters are used for quality control purposes.
# Filtering for presence of a BIN URI, a form of BIN identification.
dfFiltered <- dfFiltered[grep("[:]", bin_uri)]
# Remove "BOLD:" from the BIN URIs.
dfFiltered[, bin_uri := substr(bin_uri, 6, 13)]

### FILTER 2 ###
# Filtering for presence of a sequence.
dfFiltered <- dfFiltered[grep("[CTG]", nucleotides)]

### INITIAL BIN SIZE ###
# Determine how many sequences are in a BIN in total prior to any sequence filtering (i.e. # of original raw sequences).
dfFiltered[, initial_bin_size := length(recordID), keyby = bin_uri]

### FILTER 3 ###
# Filtering for COI-5P as these are the only markers we are looking at.
dfFiltered <- dfFiltered[markercode == "COI-5P"]

### FILTER 4 ###
# Trim sequences with high N and gap content at their terminal ends.
# First, make sure nucleotides are "chr" type. 
dfFiltered[, nucleotides := as.character(nucleotides)]
# Trim large portions of Ns and gaps at the start of a sequence. 
# Find sequences that begin (^) with an N or a gap.
startGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("^[-N]", dfFiltered$nucleotides)), length)
# Loop through all of the sequences.
startGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
  # If at least one sequence is found that begins with a gap or N (is flagged as a 1)...
  if (startGapN[[i]] > 0) { 
    # Split the sequence up using strsplit!
    # Using a regex to find sequences that begin with gaps and/or Ns.
    split <- strsplit(dfFiltered$nucleotides[i], "^[-N]+")
    # Take only the second half of the element (the sequence without the gaps/Ns
    # at the start!).
    dfFiltered$nucleotides[i] <- split[[1]][2]
  }
# Trim large portions of Ns and gaps at the end of a sequence. 
# Find sequences that end ($) with an N or a gap.
endGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("[-N]$", dfFiltered$nucleotides)), length)
endGapN <- foreach(i = 1:nrow(dfFiltered)) %do%
  if (endGapN[[i]] > 0) {
    # Using a regex to find sequences that end with gaps and/or Ns.
    split <- strsplit(dfFiltered$nucleotides[i], "[-N]+$")
    # Take only the first half of the element (the sequence without the gaps/Ns 
    # at the end!).
    dfFiltered$nucleotides[i] <- split[[1]][1]
  }

### FILTER 5 ###
# Remove sequences with N/gap content above a certain threshold (1%).
# Determine the number of positions where an *internal* N or gap is found for each sequence.
internalGapN <- sapply(regmatches(dfFiltered$nucleotides, gregexpr("[-N]", dfFiltered$nucleotides)), length)
# Iterate over each sequence and see if the number of Ns or gaps is greater than 1% (0.01) of the total length of the sequence.
internalGapN <- foreach(i = 1:nrow(dfFiltered)) %do% which((internalGapN[[i]]/nchar(dfFiltered$nucleotides[i]) > 0.01))
# Identify the sequences with high N/gap content.
checkGapN <- which(internalGapN > 0) 
# Remove these sequences.
checkGapN <- dfFiltered[-checkGapN, ]

### FILTER 6 ###
# Filter out sequences that are less than 640 bp and greater than 1000 bp. 
# Determine the lengths of the sequences without gaps.
seqLengths <- dfFiltered[, nchar(gsub("-", "", nucleotides))] 
# Which sequences are greater than 1000 bp and less than 640 bp in length?
seqLengthCheck <- which(seqLengths > 1000 | seqLengths < 640)
# Remove these sequences.
dfFiltered <- dfFiltered[-seqLengthCheck, ]

# BIN Species Information. #
# Here, we are obtaining information on a per BIN basis to facilitate trait matching later on.

### FILTER 7 ###
# Remove rows with no species information. This will remove BINs without any species information. BINs without species data would 
# not match with any trait information down the line.
# Create a new datatable containing only sequences baring species-level identification. This is necessary so we can extract the 
# BIN URIs that contain species-level identification and remove those without. 
dfSpecies <- dfFiltered[grep("[A-Z]", species_name)]
# Now we have the URIs of BINs that contain sequences with species information.
speciesBins <- unique(dfSpecies$bin_uri)
# Subset out these BINs in dfFiltered.
dfResolve <- dfFiltered[dfFiltered$bin_uri %in% speciesBins]

# RESOLVING TAXONOMIC CONFLICTS (MERGES).
# These steps are performed to improve BIN reliability and ensure we are matching the appropriate sequence information to the 
# appropriate trait data seqLengthCheckdown the line.

# First, I need to replace all blanks with NA values in the taxonomy columns. This is to ensure that empty cells are not counted as taxa.
dfResolve[dfResolve == ""] <- NA

# Order level conflicts.
# Find the number of orders in each BIN.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), keyby = bin_uri]
# Which BINs have more than 1 order assigned to them?
orderConflicts <- CountConflicts(dfResolve, "number_of_orders")
orderConflicts
# If there are more than 0 conflicts, you may apply the ResolveBIN function to remove a deviant record from a BIN. 
# You may want to do this if there are a few weird records in your BIN, but the BIN seems reliable otherwise.

# EXAMPLE: the following line would remove a specific record from the dataframe.
# dfResolve <- ResolveBIN(dfResolve, 8464523, method = "recordID")

# EXAMPLE: Or you can remove an entire BIN from the dataset. You may want to do this if
# there is no clear consensus in the BIN for an order level assignment.
# The following lines would remove entire BINs from the dataframe.:
# dfResolve <- ResolveBIN(dfResolve, "AAD3116", method = "bin_uri")

# Family level conflicts.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), keyby = bin_uri]
familyConflicts <- CountConflicts(dfResolve, "number_of_families")
familyConflicts

# Genus level conflicts.
# There are probably going to be a lot more genus and species level conflicts, so we will not be able to check them manually. 
# Instead, we will keep only those BINs that have AT LEAST 10 records and that have 80% consistency for genus or species level assignment 
# (i.e. at least 8 out of 10 records share the same genus or species level assignment). You can change these thresholds if
# you want to make them more or less strict.
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), keyby = bin_uri]
genusConflicts <- CountConflicts(dfResolve, "number_of_genera")
# Create a new datatable for BINs with genus level conflicts.
dfGenusConflicts <- dfResolve[bin_uri %in% genusConflicts]
# Now we must determine the most common genus and if it has at least 80% consistency in sequences that DO have genus level information.
# Only looking at sequences with genus classifications.
dfGenusConflicts <- dfGenusConflicts[grep("[A-Z]", genus_name)]
# Create a new column for the number of sequences with genus level information per BIN.
dfGenusConflicts[, number_of_seqs := length(recordID), by = bin_uri]
# Which bins have more than 10 sequences? These are probably more reliable.
dfGenusConflicts <- dfGenusConflicts[number_of_seqs >= 10]
# A count column is created to count the number of rows per genus per BIN.
# This is necessary to calculate the percentage of sequences from each genus per BIN.
dfGenusConflicts[, count := .N, by = .(bin_uri, genus_name)]
# Calculate the percentage of sequences from each genus per BIN.
dfGenusConflicts[, genus_percentage := .(count / number_of_seqs)]
dfGenusConflicts[order(-count), majority_genus := genus_name[1L], by = bin_uri]
# Reorder to take a closer look.
dfGenusConflicts <- dfGenusConflicts[, .(bin_uri, genus_name, majority_genus, number_of_seqs, count, genus_percentage)]
# Make a column for majority species percentage to test if it is over 80%. This is the percentage for the genus with the majority of entries.
dfGenusConflicts[order(-genus_percentage), majority_genus_percentage := genus_percentage[1L], by = bin_uri]
# Subset out those BINs that have a majority genera over 80%.
dfAcceptedGenus <- dfGenusConflicts[majority_genus_percentage > 0.80]
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# unacceptedBins = BINs in genusConflicts which were not accepted.
unacceptedBins <- setdiff(genusConflicts, unique(dfAcceptedGenus$bin_uri))
dfResolve <- dfResolve[!dfResolve$bin_uri %in% unacceptedBins]

# # Species level conflicts.
# Repeat the same process for species as we did for genera.
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), keyby = bin_uri]
speciesConflicts <- CountConflicts(dfResolve, "number_of_species")
dfSpeciesConflicts <- dfResolve[bin_uri %in% speciesConflicts]
dfSpeciesConflicts <- dfSpeciesConflicts[grep("[A-Z]", species_name)]
dfSpeciesConflicts[, number_of_seqs := length(recordID), by = bin_uri]
dfSpeciesConflicts <- dfSpeciesConflicts[number_of_seqs >= 10]
dfSpeciesConflicts[, count := .N, by = .(bin_uri, species_name)]
dfSpeciesConflicts[, species_percentage := .(count / number_of_seqs)]
dfSpeciesConflicts[order(-count), majority_species := species_name[1L], by = bin_uri]
dfSpeciesConflicts <- dfSpeciesConflicts[, .(bin_uri, species_name, majority_species, number_of_seqs, count, species_percentage)]
dfSpeciesConflicts[order(-species_percentage), majority_species_percentage := species_percentage[1L], by = bin_uri]
dfAcceptedSpecies <- dfSpeciesConflicts[majority_species_percentage > 0.80]
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# unacceptedBins = BINs in genusConflicts which were not accepted.
unacceptedBins <- setdiff(speciesConflicts, unique(dfAcceptedSpecies$bin_uri))
dfResolve <- dfResolve[!dfResolve$bin_uri %in% unacceptedBins]

##################################################################################################################
### TRAIT: POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfResolve[, filtered_bin_size := length(recordID), by = bin_uri]
##################################################################################################################

### SPECIES LABEL ###
# Now, we want to assign every sequence in a BIN a taxonomic label at the species level. This will ensure that even those sequences 
# with discordant taxonomic classifications will share a common name with the "accepted" taxonomic assignment for their BIN.
# First, create a new datatable containing only sequences bearing taxonomic identification at the species level. This is necessary because 
# NA values are considered when counting the number of species.
# Species label.
dfSpeciesLabel <- AssignLabel(dfResolve, "species_name", "species_label")
# Genus label.
dfGenusLabel <- AssignLabel(dfResolve, "genus_name", "genus_label")
# Family label.
dfFamilyLabel <- AssignLabel(dfResolve, "family_name", "family_label")
# Order label.
dfOrderLabel <- AssignLabel(dfResolve, "order_name", "order_label")

# MERGING DATATABLES.
# Merge datatables containing BIN species information.
dfFiltered <- Reduce(function(...) merge(..., all = T), list(dfResolve, dfOrderLabel, dfFamilyLabel, dfGenusLabel, dfSpeciesLabel))
# Datatable organization. Check this datatable to make sure it is accurate!
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, order_label, family_name, family_label, genus_name, genus_label,
                             species_name, species_label, nucleotides, filtered_bin_size, lat)]

# RESOLVING TAXONOMIC CONFLICTS (SPLITS).
# Now we will select only a single BIN per species name to resolve "splits".
# We only want the BIN names so we will filter dfFiltered to only one sequence for now.
dfSplits <- dfFiltered[!duplicated(bin_uri)]
dfSplits <- dfSplits[, .SD[which.max(filtered_bin_size)], keyby = species_label]
# Now subset dfFiltered and remove the smaller SPLIT BINs.
dfFiltered <- dfFiltered[dfFiltered$bin_uri %in% dfSplits$bin_uri]
