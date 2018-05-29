# Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Matthew Orton (https://github.com/m-orton/R-Scripts) for testing/contributions to the sequence filters.
# Dr. Robert Hanner for recommendations about how to deal with BIN data.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

#############################################################################################################################

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
#install.packages(stringr)
library(stringr)
# Load the function(s) designed for this script:
source("ResolveBIN.R")
source("CountConflicts.R")
source("AssignLabel.R")

#############################################################################################################################

# Download sequences from BOLD using the function bold_seqspec() for sequence and specimen data. 
# In addition, I am only selecting those columns needed for downstream analysis.
# Enter your taxon between the "".
dfRawSeqs <- bold_seqspec(taxon = "Scorpaeniformes", geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", 
                                                                      "genus_name", "species_name", "lat", "nucleotides", 
                                                                      "markercode")]

# Download outgroup species data from BOLD. These sequences may be used to root phylogenetic trees (depending if the taxa 
# are an appropriate outgroup for the organismal group under study).
# Enter your outgroup name between the "".
outgroups <- c("Neoceratodus forsteri") 
dfOutgroup <- bold_seqspec(taxon = outgroups, geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", 
                                                               "genus_name", "species_name", "lat", "nucleotides", 
                                                               "markercode")]
# Combine dfOutgroup and dfRawSeqs so that they are in one useable dataframe. Also, convert to datatable as datatables 
# have useful features for data manipulation.
dfFiltered <- as.data.table(rbind(dfRawSeqs, dfOutgroup))

### FILTER 1 ###
# Filters are used for quality control purposes.
# Filtering for presence of a BIN URI, a form of BIN identification.
dfFiltered <- dfFiltered[bin_uri %like% "[:]"]
# Remove "BOLD:" from the BIN URIs.
dfFiltered[, bin_uri := substr(bin_uri, 6, 13)]

### FILTER 2 ###
# Filtering for presence of a sequence.
dfFiltered <- dfFiltered[nucleotides %like% "[CTG]"]

### INITIAL BIN SIZE ###
# Determine how many sequences are in a BIN in total prior to any sequence filtering (i.e. # of original raw sequences).
dfFiltered[, initial_bin_size := length(recordID), keyby = bin_uri]

### FILTER 3 ###
# Filtering for COI-5P as these are the only markers we are looking at.
dfFiltered <- dfFiltered[markercode == "COI-5P"]

### FILTER 4 ###
# Trim sequences with high N and gap content at their terminal ends.
# Find sequences that begin (^) with an N or a gap and split them using tstrsplit. Take only the second half of the split.
# strsplit gives a blank first element here when splitting at the beginning of the string, so that is why the chained 
# command is neccessary.
dfFiltered[, c("split1", "split2") := tstrsplit(nucleotides, "^[-N]+")][split1 == "", nucleotides := split2]
# Remove extra columns.
dfFiltered[, c("split1","split2") := NULL]
# Trim large portions of Ns and gaps at the end of a sequence. 
# Find sequences that end ($) with an N or a gap and split them using tstrsplit.
dfFiltered[, nucleotides := tstrsplit(nucleotides, "[-N]+$")]

### FILTER 5 ###
# Remove sequences with N/gap content above a certain threshold (1%).
# Determine the number of positions where an *internal* N or gap is found for each sequence.
dfFiltered[, internal_gapN := str_count(nucleotides, c("[N-]"))]
# Remove sequence if the number of Ns or gaps is greater than 1% (0.01) of the total length of the sequence.
dfFiltered <- dfFiltered[, percentage_gapN := internal_gapN/nchar(nucleotides)][!percentage_gapN > 0.01]
# Remove these columns as they are no longer needed.
dfFiltered[, c("internal_gapN", "percentage_gapN") := NULL]

### FILTER 6 ###
# Filter out sequences (when gaps are not included) that are less than 640 bp and greater than 1000 bp. 
dfFiltered <- dfFiltered[nchar(gsub("-", "", nucleotides)) %between% c(640, 1000)]

# BIN Species Information. #
# Here, we are obtaining information on a per BIN basis to facilitate trait matching later on.

### FILTER 7 ###
# Remove rows with no species information. This will remove BINs without any species information. BINs without species data 
# would not match with any trait information down the line.
# Create a new datatable containing only sequences baring species-level identification. This is so we can extract the 
# BIN URIs that contain species-level identification and remove BIN URIs without species information. 
dfSpecies <- dfFiltered[species_name %like% "[A-Z]"]
# Subset out the BINs containing species-level information in dfFiltered.
dfResolve <- dfFiltered[bin_uri %in% dfSpecies$bin_uri]

# RESOLVING TAXONOMIC CONFLICTS (MERGES).
# These steps are performed to improve BIN reliability and ensure we are matching the appropriate sequence information to the 
# appropriate trait data down the line.

# First, I need to replace all blanks with NA values in the taxonomy columns. This is to ensure that empty cells are not 
# counted as taxa.
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
# Instead, we will keep only those BINs that have AT LEAST 10 records and that have 80% consistency for genus or species 
# level assignment (i.e. at least 8 out of 10 records share the same genus or species level assignment). You can change these 
# thresholds if you want to make them more or less strict.
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), keyby = bin_uri]
genusConflicts <- CountConflicts(dfResolve, "number_of_genera")
# Create a new datatable for BINs with genus level conflicts.
dfGenusConflicts <- dfResolve[bin_uri %in% genusConflicts]
# Now we must determine the most common genus and if it has at least 80% consistency in sequences that DO have genus level information.
# Only looking at records with genus classifications.
dfGenusConflicts <- dfGenusConflicts[genus_name %like% "[A-Z]"]
# Create a new column for the number of sequences with genus level information per BIN and only take BINs with more than 
# 10 sequences (they are probably more reliable).
dfGenusConflicts <- dfGenusConflicts[, number_of_seqs := length(recordID), by = bin_uri][number_of_seqs >= 10]
# A count column is created to count the number of sequences per genus per BIN. This is necessary to calculate the 
# percentage of sequences from each genus per BIN.
dfGenusConflicts[, count := .N, by = .(bin_uri, genus_name)]
# Calculate the percentage of sequences from each genus per BIN.
dfGenusConflicts[, genus_percentage := .(count / number_of_seqs)]
# The majority genus is the one with the largest percentage (in descending order i.e. -count).
dfGenusConflicts[order(-count), majority_genus := genus_name[1L], by = bin_uri]
# Make a column for majority species percentage to test if it is over 80%. This is the percentage for the genus with the 
# majority of entries.
dfGenusConflicts[order(-genus_percentage), majority_genus_percentage := genus_percentage[1L], by = bin_uri]
# Reorganize to double check.
dfGenusConflicts <- dfGenusConflicts[, .(bin_uri, genus_name, majority_genus, number_of_seqs, count, 
                                         genus_percentage, majority_genus_percentage)]
# Subset out those BINs that have a majority genera over 80%.
dfAcceptedGenus <- dfGenusConflicts[majority_genus_percentage > 0.80]
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# unacceptedBins = BINs in genusConflicts which were not accepted.
unacceptedBins <- setdiff(genusConflicts, unique(dfAcceptedGenus$bin_uri))
dfResolve <- dfResolve[!bin_uri %in% unacceptedBins]

# # Species level conflicts.
# Repeat the same process for species as we did for genera.
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), keyby = bin_uri]
speciesConflicts <- CountConflicts(dfResolve, "number_of_species")
dfSpeciesConflicts <- dfResolve[bin_uri %in% speciesConflicts]
dfSpeciesConflicts <- dfSpeciesConflicts[species_name %like% "[A-Z]"]
dfSpeciesConflicts <- dfSpeciesConflicts[, number_of_seqs := length(recordID), by = bin_uri][number_of_seqs >= 10]
dfSpeciesConflicts[, count := .N, by = .(bin_uri, species_name)]
dfSpeciesConflicts[, species_percentage := .(count / number_of_seqs)]
dfSpeciesConflicts[order(-count), majority_species := species_name[1L], by = bin_uri]
dfSpeciesConflicts[order(-species_percentage), majority_species_percentage := species_percentage[1L], by = bin_uri]
dfSpeciesConflicts <- dfSpeciesConflicts[, .(bin_uri, species_name, majority_species, number_of_seqs, count, 
                                             species_percentage, majority_species_percentage)]
dfAcceptedSpecies <- dfSpeciesConflicts[majority_species_percentage > 0.80]
unacceptedBins <- setdiff(speciesConflicts, unique(dfAcceptedSpecies$bin_uri))
dfResolve <- dfResolve[!bin_uri %in% unacceptedBins]

### POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfResolve[, filtered_bin_size := length(recordID), by = bin_uri]

### SPECIES LABEL ###
# Now, we want to assign every sequence in a BIN taxonomic labels. This will ensure that even those sequences 
# with discordant taxonomic classifications will share a common name with the "accepted" taxonomic assignment for their BIN. 
# This is necessary for matching trait information using species names.
# First, create a new datatable containing only sequences bearing taxonomic identification at the species level. 
# This is necessary because NA values are considered when counting the number of species.
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
# We only want the BIN URIs so we will filter dfFiltered to only one sequence for now.
dfSplits <- dfFiltered[!duplicated(bin_uri)]
dfSplits <- dfSplits[, .SD[which.max(filtered_bin_size)], keyby = species_label]
# Now subset dfFiltered and remove the smaller SPLIT BINs.
dfFiltered <- dfFiltered[bin_uri %in% dfSplits$bin_uri]