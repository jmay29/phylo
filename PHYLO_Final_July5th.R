###################
# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
# correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
# determination/reference sequence trimming (lines TBD).
# Last version saved: Apr 2nd 2017 (R_Pipeline_PHYLO_V7_Mar24th.R)

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
#library(gtools)
#install.packages("caret")
library(caret)
#install.packages("Boruta")
library(Boruta)
#install.packages("randomForest")
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
library(adephylo)
library(phylobase)
library(car)
#install.packages("phylometrics")
library(phylometrics)

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

# Start the clock!
#ptm <- proc.time()


##### SECTION 1: DATA PROCESSING #####
# Download sequences from BOLD using the function bold_seqspec() for sequence
# and specimen data. In addition, I am only selecting those columns needed for
# downstream analysis.
dfInitial <- bold_seqspec(taxon = "Actinopterygii", geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", 
                                                                     "genus_name", "species_name", "lat", "nucleotides", 
                                                                     "markercode")] 

# Download outgroup species data from BOLD.
outgroups <- c("Elasmobranchii", "Sarcopterygii", 
               "	Rodentia", "Anura")
dfOutgroup <- bold_seqspec(taxon = outgroups, 
                           geo = "all")[, c("recordID", "bin_uri", "order_name", "family_name", 
                                            "genus_name", "species_name", "lat", "nucleotides", 
                                            "markercode")] 
dfOutgroup <- setDT(dfOutgroup)  # Again, convert to datatable.

# Combine dfOutgroup and dfSequences datatables so that they are in one useable 
# datatable.
l <- list(dfInitial, dfOutgroup)
dfFiltered <- rbindlist(l)


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
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, family_name, genus_name, species_name,
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

# GONNA WANT TO DO LATITUDE INFO HERE INSTEAD...


# BIN Species Information. #
# Here, we are obtaining information on a per BIN basis to facilitate trait 
# matching later on in the pipeline.
### FILTER 7 ###
# Remove rows with no species information (removes BINs with no species data). 
# BINs without species data would not match with any trait information down the line.
containSpecies <- dfFiltered[, grep("[A-Z]", species_name)]
# Create a new datatable containing only sequences baring species-level identification. 
# This is necessary so we can extract BINs that contain species-level identification 
# and remove those without! Can maybe make this easier by doing it the datatable way...
dfSpecies <- dfFiltered[containSpecies, ]
rm(containSpecies)
# Now we have the BINs that contain species.
speciesBins <- unique(dfSpecies$bin_uri)
# Subset out these BINs in dfFiltered.
dfResolve <- subset(dfFiltered, bin_uri %in% speciesBins)


# RESOLVING TAXONOMIC CONFLICTS.
# Now, I want to resolve BINs with more than 1 order and/or family.
# First, I need to replace all blanks with NA values in the taxonomy columns.
dfResolve[, order_name := revalue(order_name, c(" " = NA))]
dfResolve[, family_name := revalue(family_name, c(" " = NA))]
dfResolve[, genus_name := revalue(genus_name, c(" " = NA))]
dfResolve[, species_name := revalue(species_name, c(" " = NA))]

# Now find the number of orders/families/genera/species in each BIN.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), keyby = bin_uri]
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), keyby = bin_uri]
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), keyby = bin_uri]
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), keyby = bin_uri]

# First, let's deal with order level conflicts.
# Checking these manually.
# Upating number_of_orders column each time (and orderConflicts).
orderConflicts <- dfResolve[, which(number_of_orders > 1), by = bin_uri]
orderConflicts <- unique(orderConflicts$bin_uri)
# BOLD:AAA6100 removing 1 deviant seq from different order.
x <- which(dfResolve$recordID == 212145)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAB5446. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAB5446")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAB5867. Note on BOLD.
rmThese <- which(dfResolve$bin_uri == "AAB5867")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAD1382. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAD1382")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAE7526. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAE7526")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAE8407. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAE8407")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAF0799. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAF0799")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAV2783. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAV2783")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:ADE3647. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "ADE3647")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAB4180 removing 3 deviant seqs from different order.
x <- which(dfResolve$recordID == 244928 |
           dfResolve$recordID == 244926 |
           dfResolve$recordID == 244929 )
dfResolve <- dfResolve[-x, ]
# Update number_of_orders column and make sure there are no more order conflicts.
dfResolve[, number_of_orders := length(unique(order_name[!is.na(order_name)])), keyby = bin_uri]
orderConflicts <- dfResolve[, which(number_of_orders > 1), by = bin_uri]
orderConflicts <- unique(orderConflicts$bin_uri)
rm(orderConflicts)

# Family level resolving.
# Checking these manually.
# Make sure to update number_of_families column first.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), keyby = bin_uri]
# Find number of BINs with family level conflicts.
familyConflicts <- dfResolve[, which(number_of_families > 1), by = bin_uri]
familyConflicts <- unique(familyConflicts$bin_uri)
# BOLD:AAB2488 removing 1 deviant seq from different fam.
x <- which(dfResolve$recordID == 552843)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAB4315. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAB4315")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAB4318. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAB4318")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAB5416 removing 1 deviant seq from different fam.
x <- which(dfResolve$recordID == 718811)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAB5825. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAB5825")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAB5827. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAB5827")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAB5842. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAB5842")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAB5416 removing 1 deviant seq from different fam.
x <- which(dfResolve$recordID == 553287)
dfResolve <- dfResolve[-x, ]
# BOLD:AAC1457 removing 1 deviant seq from different fam.
x <- which(dfResolve$recordID == 718927)
dfResolve <- dfResolve[-x, ]
# BOLD:AAC3188 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 338753)
dfResolve <- dfResolve[-x, ]
# BOLD:AAC4093 removing 2 deviant seqs from different fam.
x <- which(dfResolve$recordID == 803396 | dfResolve$recordID == 803397)
dfResolve <- dfResolve[-x, ]
# BOLD:AAC9496 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 338648)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAD1426. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAD1426")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAD2909 removing 4 deviant seqs from different fam.
x <- which(dfResolve$recordID == 210644 | dfResolve$recordID == 210643 |
           dfResolve$recordID == 210645 | dfResolve$recordID == 210653)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAD4188. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAD4188")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAD4397. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAD4397")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAD8781. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAD8781")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAD9212 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 3020192)
dfResolve <- dfResolve[-x, ]
# BOLD:AAE0699 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 2479339)
dfResolve <- dfResolve[-x, ]
# BOLD:AAE2823 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 338687)
dfResolve <- dfResolve[-x, ]
# BOLD:AAE4499 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 338935)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:AAE6275. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "AAE6275")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:AAF3989. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "AAF3989")
dfResolve <- dfResolve[-rmThese, ]
# BOLD:AAF8738 removing 1 deviant seq from different fams.
x <- which(dfResolve$recordID == 1131517)
dfResolve <- dfResolve[-x, ]
# BOLD:AAO9420 removing 2 deviant seqs from different fam.
x <- which(dfResolve$recordID == 3020853 | dfResolve$recordID == 1508430)
dfResolve <- dfResolve[-x, ]
# Removing BOLD:ABW2356. Not enough sequences.
rmThese <- which(dfResolve$bin_uri == "ABW2356")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:ABX1827. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "ABX1827")
dfResolve <- dfResolve[-rmThese, ]
# Removing BOLD:ACD2391. No clear consensus.
rmThese <- which(dfResolve$bin_uri == "ACD2391")
dfResolve <- dfResolve[-rmThese, ]
# Update number_of_orders column and make sure there are no more order conflicts.
dfResolve[, number_of_families := length(unique(family_name[!is.na(family_name)])), keyby = bin_uri]
familyConflicts <- dfResolve[, which(number_of_families > 1), by = bin_uri]
familyConflicts <- unique(familyConflicts$bin_uri)
rm(familyConflicts)
rm(rmThese)
rm(speciesBins)
rm(x)

# Genus level resolving.
# Update the column.
dfResolve[, number_of_genera := length(unique(genus_name[!is.na(genus_name)])), keyby = bin_uri]
# At least 10+ records and 80% consistency.
genusConflicts <- dfResolve[, which(number_of_genera > 1), by = bin_uri]
genusConflictBins <- unique(genusConflicts$bin_uri)
# Create a new datatable for BINs with genus level conflicts.
dtGenus <- dfResolve[genusConflictBins, ]
#rm(genusConflictBins)
# Now we must determine the most common genus and if it has at least 80% consistency
# in sequences that DO have genus level information.
# Only looking at sequences with genus classifications.
containGenus <- dtGenus[, grep("[A-Z]", genus_name)]
# Create a new datatable containing only sequences baring genus-level identification. 
dtGenus10 <- dtGenus[containGenus, ]
# Create a new column for the number of sequences with genus level information per BIN.
dtGenus10[, number_of_seqs := length(recordID), by = bin_uri]
# Which bins have more than 10 sequences?
dtGenus10 <- dtGenus10[which(number_of_seqs >= 10)]
#rm(dtGenus)
# A count column is created to count the number of rows per genus per BIN.
dtGenus10[, count := .N, by = .(bin_uri, genus_name)]
# Majority genus percentage.
dtGenus10[, genus_percentage := .(count / number_of_seqs)]
dtGenus10[order(-count), majority_genus := genus_name[1L], by = bin_uri]
# Reorder to take a closer look.
dtGenus10 <- dtGenus10[, .(bin_uri, genus_name, majority_genus, number_of_seqs, count, genus_percentage)]
# Make a column for majority species percentage to test if it is over 80.
# The genus with the majority of entries.
dtGenus10[order(-genus_percentage), majority_genus_percentage := genus_percentage[1L], by = bin_uri]
dtAcceptedGenus <- dtGenus10[which(majority_genus_percentage > 0.80)]
#rm(dtGenus10)
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# bad = BINs in genusConflicts which are not accepted.
bad <- anti_join(genusConflicts, dtAcceptedGenus, by = "bin_uri")
#rm(genusConflicts)
#rm(dtAcceptedGenus)
badBins <- unique(bad$bin_uri)
dfResolve <- dfResolve[!dfResolve$bin_uri %in% badBins, ]

# Species level resolving.
# Update the column.
dfResolve[, number_of_species := length(unique(species_name[!is.na(species_name)])), keyby = bin_uri]
# At least 10+ records and 80% consistency.
speciesConflicts <- dfResolve[, which(number_of_species > 1), by = bin_uri]
speciesConflictBins <- unique(speciesConflicts$bin_uri)
# Create a new datatable for BINs with genus level conflicts.
dtSpecies <- dfResolve[speciesConflictBins, ]
#rm(speciesConflictBins)
# Now we must determine the most common species and if it has at least 80% consistency
# in sequences that DO have species level information.
# Only looking at sequences with genus classifications.
containSpecies <- dtSpecies[, grep("[A-Z]", species_name)]
# Create a new datatable containing only sequences baring species-level identification. 
dtSpecies10 <- dtSpecies[containSpecies, ]
# Create a new column for the number of sequences with species level information per BIN.
dtSpecies10[, number_of_seqs := length(recordID), by = bin_uri]
# Which bins have more than 10 sequences?
dtSpecies10 <- dtSpecies10[which(number_of_seqs >= 10)]
#rm(dtSpecies)
# A count column is created to count the number of rows per species per BIN.
dtSpecies10[, count := .N, by = .(bin_uri, species_name)]
# Majority species percentage.
dtSpecies10[, species_percentage := .(count / number_of_seqs)]
dtSpecies10[order(-count), majority_species := species_name[1L], by = bin_uri]
# Reorder to take a closer look.
dtSpecies10 <- dtSpecies10[, .(bin_uri, species_name, majority_species, number_of_seqs, count, species_percentage)]
# Make a column for majority species percentage to test if it is over 80.
# The genus with the majority of entries.
dtSpecies10[order(-species_percentage), majority_species_percentage := species_percentage[1L], by = bin_uri]
dtAcceptedSpecies <- dtSpecies10[which(majority_species_percentage > 0.80)]
#rm(dtGenus10)
# Find the UNACCEPTED conflicted bins and remove them from dfResolve.
# bad = BINs in genusConflicts which are not accepted.
bad <- anti_join(speciesConflicts, dtAcceptedSpecies, by = "bin_uri")
#rm(genusConflicts)
#rm(dtAcceptedGenus)
badBins <- unique(bad$bin_uri)
dfResolve <- dfResolve[!dfResolve$bin_uri %in% badBins, ]


################################################################################
### TRAIT: POST FILTER BIN SIZE ###
# Determine how many sequences are in a BIN in total after sequence filtering.
dfResolve[, filtered_bin_size := length(recordID), by = bin_uri]
################################################################################

### TAXONOMIC LABELS ###
# Determine the most common taxonomic classifications in each BIN. 
# Create a new datatable containing only sequences baring taxonomic identification at
# the corresponding level.
# This is necessary because NA values are considered when counting the number of species.
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
# Family label.
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
# Note: bin_uri is the key for each of these datatables due to the use of "keyby" 
# instead of just "by". This ultimately facilitates datatable merging.
setkey(dfResolve, bin_uri)
dfFiltered <- merge(dfResolve, dfSpeciesLabel)
rm(dfSpeciesLabel)
dfFiltered <- merge(dfFiltered, dfGenusLabel)
rm(dfGenusLabel)
dfFiltered <- merge(dfFiltered, dfFamilyLabel)
rm(dfFamilyLabel)
dfFiltered <- merge(dfFiltered, dfOrderLabel)
rm(dfOrderLabel)
# Datatable reorganization.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name, order_label, family_name, 
                             family_label, genus_name, genus_label, species_name, 
                             species_label, nucleotides, initial_bin_size,
                             filtered_bin_size, lat)]
rm(dfSpecies)

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

# Creating datatables for univariate analyses of latitude and latitudinal range traits.
# TRAIT: Latitude.
# Filtering for presence of trait data. This is for the univariate analyses section.
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
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfLatitude <- setDT(GetTraitSpecificData(dfLatitudeSpecies, 17))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes.
hist(dfLatitude$median_lat)

# Trait: Latitudinal range.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfLatRange <- setDT(GetTraitSpecificData(dfLatitudeSpecies, 20))
rm(dfLatitudeSpecies)
# TEST 2: Does the trait have enough data variation?
# First, remove all singleton BINs.
singletons <- which(dfLatRange$filtered_bin_size == 1)
dfLatRange <- dfLatRange[-singletons, ]
hist(dfLatRange$lat_range)

# Finally, prepare the dfLatitudeMV datatable by merging all univariate datatables.
# This datatable will be used for the eventual multivariate analysis.
dfLatitudeMV <- merge(dfLatitude, dfLatRange, by = "bin_uri", all = TRUE)
dfLatitudeMV <- dfLatitudeMV[, .(bin_uri, species_name = species_label.x, 
                                 filtered_bin_size = filtered_bin_size.x, 
                                 median_lat, lat_range)]


################################################################################

# Datatable reorganization for dfFiltered.
dfFiltered <- dfFiltered[, .(bin_uri, recordID, order_name = order_label, family_name = family_label, 
                             genus_name = genus_label, species_name = species_label, nucleotides, initial_bin_size, 
                             filtered_bin_size)]

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
#dfSpecies <- data.frame(species(speciesNames))
# Storing this as a file.
#write.csv(dfSpecies, file = "species_information.csv")
dfSpecies <- fread("species_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfSpeciesTraits <- dfSpecies[, .(species_name = sciname, 
                                 BodyShapeI, Fresh, Brack, Saltwater,
                                 DemersPelag, AnaCat, DepthRangeShallow,
                                 DepthRangeDeep, DepthRangeComShallow,
                                 DepthRangeComDeep, LongevityWild, Length, LTypeMaxM,
                                 CommonLength, LTypeComM, Weight)]
# TRAIT: Body Shape I.
# Filtering for presence of body shape I data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBodyShapeI <- setDT(GetTraitSpecificData(dfSpeciesTraits, 2))
# TEST 2: Does the trait have enough data variation?
# Answer: "other" is too rare and will be removed.
table(dfBodyShapeI$BodyShapeI)
rareVars <- which(dfBodyShapeI$BodyShapeI == "other")
dfBodyShapeI <- dfBodyShapeI[-rareVars, ]
# Make it a factor variable.
dfBodyShapeI[, BodyShapeI := as.factor(BodyShapeI)]
# Also dropping this levels from the factor.
dfBodyShapeI$BodyShapeI <- droplevels(dfBodyShapeI$BodyShapeI)

# TRAIT: Fresh.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfFreshwater <- setDT(GetTraitSpecificData(dfSpeciesTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfFreshwater$Fresh)
# Make it a factor variable.
dfFreshwater[, Fresh := as.factor(Fresh)]

# TRAIT: Fresh.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBrackish <- setDT(GetTraitSpecificData(dfSpeciesTraits, 4))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfBrackish$Brack)
# Make it a factor variable.
dfBrackish[, Brack := as.factor(Brack)]

# TRAIT: Saltwater.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfSaltwater <- setDT(GetTraitSpecificData(dfSpeciesTraits, 5))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfSaltwater$Saltwater)
# Make it a factor variable.
dfSaltwater[, Saltwater := as.factor(Saltwater)]

# TRAIT: Demers/Pelag.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDemersPelag <- setDT(GetTraitSpecificData(dfSpeciesTraits, 6))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfDemersPelag$DemersPelag)
# Make it a factor variable.
dfDemersPelag[, DemersPelag := as.factor(DemersPelag)]

# TRAIT: Ana/Cat.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfAnaCat <- setDT(GetTraitSpecificData(dfSpeciesTraits, 7))
# TEST 2: Does the trait have enough data variation?
# Answer: Rare categories removed.
table(dfAnaCat$AnaCat)
# Answer: "other" is too rare and will be removed.
rareVars <- which(dfAnaCat$AnaCat == "" |
                  dfAnaCat$AnaCat == "amphidromous?" |
                  dfAnaCat$AnaCat == "diadromous")
dfAnaCat <- dfAnaCat[-rareVars, ]
# Make it a factor variable.
dfAnaCat[, AnaCat := as.factor(AnaCat)]
# Also dropping this levels from the factor.
dfAnaCat$AnaCat <- droplevels(dfAnaCat$AnaCat)

# TRAIT: DepthRangeShallow.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDepthShallow <- setDT(GetTraitSpecificData(dfSpeciesTraits, 8))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes.
hist(dfDepthShallow$DepthRangeShallow)
# Make it a numeric variable.
dfDepthShallow[, DepthRangeShallow := as.double(DepthRangeShallow)]

# TRAIT: DepthRangeDeep.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDepthDeep <- setDT(GetTraitSpecificData(dfSpeciesTraits, 9))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfDepthDeep$DepthRangeDeep)
# Make it a numeric variable.
dfDepthDeep[, DepthRangeDeep := as.double(DepthRangeDeep)]

# TRAIT: DepthRangeComShallow.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDepthComShallow <- setDT(GetTraitSpecificData(dfSpeciesTraits, 10))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfDepthComShallow$DepthRangeComShallow)
# Make it a numeric variable.
dfDepthComShallow[, DepthRangeComShallow := as.double(DepthRangeComShallow)]

# TRAIT: DepthRangeComDeep.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDepthComDeep <- setDT(GetTraitSpecificData(dfSpeciesTraits, 11))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfDepthComDeep$DepthRangeComDeep)
# Make it a numeric variable.
dfDepthComDeep[, DepthRangeComDeep := as.double(DepthRangeComDeep)]

# TRAIT: LongevityWild.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfLongWild <- setDT(GetTraitSpecificData(dfSpeciesTraits, 12))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfLongWild$LongevityWild)
# Make it a numeric variable.
dfLongWild[, LongevityWild := as.double(LongevityWild)]
range(dfLongWild$LongevityWild)

# TRAIT: Maximum length.
# The column must first be converted to double (numeric) type.
dfMaxLength <- dfSpeciesTraits[, Length := as.double(Length)]
# Only want total length measurements.
keep <- dfMaxLength[, which(LTypeMaxM == "TL")]
dfMaxLength <- dfMaxLength[keep, ]
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfMaxLength <- setDT(GetTraitSpecificData(dfMaxLength, 13))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfMaxLength$Length)

# TRAIT: Common length.
# The column must first be converted to double (numeric) type.
dfComLength <- dfSpeciesTraits[, CommonLength := as.double(CommonLength)]
# Only want total length measurements.
keep <- dfComLength[, which(LTypeComM == "TL")]
dfComLength <- dfComLength[keep, ]
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfComLength <- setDT(GetTraitSpecificData(dfComLength, 15))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfComLength$CommonLength)

# TRAIT: Weight.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfWeight <- setDT(GetTraitSpecificData(dfSpeciesTraits, 17))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfWeight$Weight)
# Make it a numeric variable.
dfWeight[, Weight := as.double(Weight)]


### STOCKS TRAITS ###
#dfStocks <- data.frame(stocks(speciesNames))
# Storing this as a file.
#write.csv(dfStocks, file = "stocks_information.csv")
dfStocks <- fread("stocks_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfStocksTraits <- dfStocks[, .(species_name = sciname, 
                                 Level, TempMin, TempMax, EnvTemp)]
# Only want species-level information.
dfStocksTraits[, Level := revalue(Level, c("Species in general" = "species in general"))]
keep <- dfStocksTraits[, which(Level == "species in general")]
dfStocksTraits <- dfStocksTraits[keep, ]

# TRAIT: TempMin.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfTempMin <- setDT(GetTraitSpecificData(dfStocksTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfTempMin$TempMin)
# Make it a numeric variable.
dfTempMin[, TempMin := as.double(TempMin)]

# TRAIT: TempMax.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfTempMax <- setDT(GetTraitSpecificData(dfStocksTraits, 4))
# TEST 2: Does the trait have enough data variation?
# Answer: 
hist(dfTempMax$TempMax)
# Make it a numeric variable.
dfTempMax[, TempMax := as.double(TempMax)]

# TRAIT: EnvTemp.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfEnvTemp <- setDT(GetTraitSpecificData(dfStocksTraits, 5))
# TEST 2: Does the trait have enough data variation?
# Answer: 
table(dfEnvTemp$EnvTemp)
rareVars <- which(dfEnvTemp$EnvTemp == "boreal" |
                    dfEnvTemp$EnvTemp == "high altitude")
dfEnvTemp <- dfEnvTemp[-rareVars, ]
# Make it a factor variable.
dfEnvTemp[, EnvTemp := as.factor(EnvTemp)]
# Also dropping this levels from the factor.
dfEnvTemp$EnvTemp <- droplevels(dfEnvTemp$EnvTemp)

# Finally, prepare the dfSpeciesGenMV datatable by merging all univariate datatables.
# This datatable will be used for the eventual multivariate analysis.
dfSpeciesGenMV <- merge(dfBodyShapeI, dfFreshwater, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfBrackish, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfSaltwater, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfAnaCat, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfEnvTemp, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfLongWild, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfMaxLength, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfWeight, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfTempMin, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfTempMax, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfComLength, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfDemersPelag, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfDepthComShallow, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfDepthComDeep, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfDepthShallow, all = TRUE, by = "species_name")
dfSpeciesGenMV <- merge(dfSpeciesGenMV, dfDepthDeep, all = TRUE, by = "species_name")


### MORPHOLOGY TRAITS ###
#dfMorphology <- data.frame(morphology(speciesNames))
# Storing this as a file.
#write.csv(dfMorphology, file = "morphology_information.csv") 
# Read in the ecology information.
dfMorphology <- fread("morphology_information.csv")
# Get rid of columns I do not need for the regression analysis.
dfMorphologyTraits <- dfMorphology[, .(species_name = sciname, body_shape_II = BodyShapeII, 
                                       pos_of_mouth = PosofMouth, type_of_scales = TypeofScales)]
# These traits are unlikely to differ between different stocks of fish so I can 
# take a single row per species.
dfMorphologyTraits <- dfMorphologyTraits[!duplicated(dfMorphologyTraits$species_name), ] 

# Recoding dfMorphologyTraits.
# Converting to factor type as these are discrete traits.
changeVars <- c(2:4)
dfMorphologyTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]

# TRAIT: Body Shape II.
# Filtering for presence of body shape II data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBodyShapeII <- setDT(GetTraitSpecificData(dfMorphologyTraits, 2))
# TEST 2: Does the trait have enough data variation?
# Answer: "other" and "angular" are too rare and will be removed.
table(dfBodyShapeII$body_shape_II)
rareVars <- which(dfBodyShapeII$body_shape_II == "other (see Diagnosi" | dfBodyShapeII$body_shape_II == "angular")
dfBodyShapeII <- dfBodyShapeII[-rareVars, ]
# Also dropping this levels from the factor.
dfBodyShapeII$body_shape_II <- droplevels(dfBodyShapeII$body_shape_II)

# TRAIT: Position of Mouth.
# Filtering for presence of position of mouth data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfPosMouth <- setDT(GetTraitSpecificData(dfMorphologyTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: Yes!
table(dfPosMouth$pos_of_mouth)


# TRAIT: Type of Scales.
# First, revalue trait so "cycloid" is equivalent to "cycloid scales".
dfMorphologyTraits[, type_of_scales := revalue(type_of_scales, c("cycloid" = "cycloid scales"))]
# Filtering for presence of type of scales data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfScaleType <- setDT(GetTraitSpecificData(dfMorphologyTraits, 4))
rm(dfMorphologyTraits)
# TEST 2: Does the trait have enough data variation?
# Answer: "other" and "rhombic scales" are too rare and are removed.
table(dfScaleType$type_of_scales)
rareVars <- which(dfScaleType$type_of_scales == "other (see remark)" | dfScaleType$type_of_scales == "rhombic scales")
dfScaleType <- dfScaleType[-rareVars, ]
# Also dropping this levels from the factor.
dfScaleType$type_of_scales <- droplevels(dfScaleType$type_of_scales)


# Finally, prepare the dfMorphologyMV datatable by merging all univariate datatables.
# This datatable will be used for the eventual multivariate analysis.
dfMorphologyMV <- merge(dfBodyShapeII, dfOperPresent, all = TRUE, by = "species_name")
dfMorphologyMV <- merge(dfMorphologyMV, dfPosMouth, all = TRUE, by = "species_name")
dfMorphologyMV <- merge(dfMorphologyMV, dfScaleType, all = TRUE, by = "species_name")


### ECOLOGY TRAITS ###
#dfEcology <- data.frame(ecology(speciesNames))
# Storing this as a file.
#write.csv(dfEcology, file = "ecology_information.csv") 
# Read in the ecology information.
dfEcology <- fread("ecology_information.csv")
colnames(dfEcology)[3] <- "species_name"
# Get rid of columns I do not need for the regression analysis. 
# Note: There is only one row per species in dfEcology.
dfEcologyTraits <- dfEcology[, c(3, 7:27, 29, 31, 33, 39, 52, 66:78,
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

# There are a lot of traits here so I am going to use nearZeroVar to cut a lot of them at once.
# AKA Cutting traits that have little variation (rarer categories less than 1%).
# The remaining traits automically pass the variation test but I will still look more
# closely at non-binary traits.
dfEcologyTraits <- as.data.frame(dfEcologyTraits)
dfEcologyTraits <- dfEcologyTraits[, -nearZeroVar(dfEcologyTraits, freqCut = 99/1)]


# Categorical traits.
# Binary traits.

# TRAIT: Neritic
dfNeritic <- setDT(GetTraitSpecificData(dfEcologyTraits, 2))

# TRAIT: Oceanic
dfOceanic <- setDT(GetTraitSpecificData(dfEcologyTraits, 4))

# TRAIT: Epipelagic
dfEpipelagic <- setDT(GetTraitSpecificData(dfEcologyTraits, 5))

# TRAIT: Mesopelagic
dfMesopelagic <- setDT(GetTraitSpecificData(dfEcologyTraits, 6))

# TRAIT: Bathypelagic
dfBathypelagic <- setDT(GetTraitSpecificData(dfEcologyTraits, 7))

# TRAIT: Estuaries
dfEstuaries <- setDT(GetTraitSpecificData(dfEcologyTraits, 8))

# TRAIT: Mangroves
dfMangroves <- setDT(GetTraitSpecificData(dfEcologyTraits, 9))

# TRAIT: Stream.
# Filtering for presence of type of stream data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfStreams <- setDT(GetTraitSpecificData(dfEcologyTraits, 10))

# TRAIT: Lakes.
# Filtering for presence of type of lake data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfLakes <- setDT(GetTraitSpecificData(dfEcologyTraits, 11))

# TRAIT: Schooling.
# Filtering for presence of type of lake data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfSchooling <- setDT(GetTraitSpecificData(dfEcologyTraits, 16))

# TRAIT: Benthic.
# Filtering for presence of type of lake data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfBenthic <- setDT(GetTraitSpecificData(dfEcologyTraits, 17))

# TRAIT: Coral reefs.
# Filtering for presence of type of lake data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfCoralReefs <- setDT(GetTraitSpecificData(dfEcologyTraits, 27))


# Non-binary traits.
# TRAIT: Feeding Type.
# Filtering for presence of type of feeding type data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfFeedingType <- setDT(GetTraitSpecificData(dfEcologyTraits, 13))
# TEST 2: Does the trait have enough data variation?
# Answer: Many categories do not reach the 1% threshold and are removed.
table(dfFeedingType$FeedingType)
rareVars <- which(dfFeedingType$FeedingType == "feeding on a host (parasite)" |
                  dfFeedingType$FeedingType == "feeding on dead animals (scavenger)" |
                  dfFeedingType$FeedingType == "feeding on the prey of a host (commensal)" |
                  dfFeedingType$FeedingType == "other" | 
                  dfFeedingType$FeedingType == "picking parasites off a host (cleaner)" |
                  dfFeedingType$FeedingType == "plants/detritus+animals (troph. 2.2-2.79)" |
                  dfFeedingType$FeedingType == "sucking food-containing material")
dfFeedingType <- dfFeedingType[-rareVars, ]
# Also dropping these levels from the factor.
dfFeedingType$FeedingType <- droplevels(dfFeedingType$FeedingType)


# Continuous traits.
# TRAIT: Diet Troph.
# Filtering for presence of type of diet troph data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfDietTroph <- setDT(GetTraitSpecificData(dfEcologyTraits, 14))
# TEST 2: Does the trait have enough data variation?
hist(dfDietTroph$DietTroph)

# TRAIT: Food Troph.
# Filtering for presence of type of diet troph data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfFoodTroph <- setDT(GetTraitSpecificData(dfEcologyTraits, 15))
# TEST 2: Does the trait have enough data variation?
hist(dfFoodTroph$FoodTroph)
rm(dfEcologyTraits)

# Finally, prepare the dfEcologyMV datatable by merging all univariate datatables.
dfEcologyMV <- merge(dfNeritic, dfOceanic, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfEpipelagic, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfMesopelagic, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfBathypelagic, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfEstuaries, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfMangroves, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfStreams, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfLakes, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfHerbivory, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfFeedingType, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfDietTroph, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfFoodTroph, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfSchooling, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfBenthic, all = TRUE, by = "species_name")
dfEcologyMV <- merge(dfEcologyMV, dfCoralReefs, all = TRUE, by = "species_name")



### LIFE HISTORY RELATED ###
# Maturity.
#dfMaturity <- data.frame(maturity(speciesNames))
#write.csv(dfMaturity, file = "maturity_information.csv") 
# Read in the maturity information.
dfMaturity <- fread("maturity_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfMaturityTraits <- dfMaturity[, .(species_name = sciname, age_at_maturity = tm)]
# Should sex be a control variable here?

# Median trait(s).
# TRAIT: Age at maturity.
# The column must first be converted to double (numeric) type.
dfAgeMaturity <- dfMaturityTraits[, age_at_maturity := as.double(age_at_maturity)]
rm(dfMaturityTraits)
# The median value is then determined for each species.
dfAgeMaturity[, age_at_maturity := median(age_at_maturity, na.rm = TRUE), keyby = species_name]
# Filtering for presence of average depth data. This is for the univariate analyses section.
# 500 SPECIES TEST: 374 but keeping because interesting trait.
dfAgeMaturity <- setDT(GetTraitSpecificData(dfAgeMaturity, 2))
# TEST 2: Does the trait have enough data variation?
# Answer: Skewed.
hist(dfAgeMaturity$age_at_maturity)
# Only one maturity trait currently unless we use sex as control variable.
dfMaturityMV <- dfAgeMaturity 


# Reproduction.
#dfReproduction <- data.frame(reproduction(speciesNames))
#write.csv(dfReproduction, file = "reproduction_information.csv") 
# Read in the reproduction information.
dfReproduction <- fread("reproduction_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfReproTraits <- dfReproduction[, .(species_name = sciname, repro_mode = ReproMode, 
                                    fertilization = Fertilization, rep_guild_1 = RepGuild1, 
                                    rep_guild_2 = RepGuild2, parental_care = ParentalCare)]
# Unlikely to vary within species.
dfReproTraits <- dfReproTraits[!duplicated(dfReproTraits$species_name), ]

x <- dfReproTraits[duplicated(dfReproTraits$species_name), ]
# First, change all of the traits to factor type.
changeVars <- c(2:6)
dfReproTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]


# TRAIT: Repro mode.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfReproMode <- setDT(GetTraitSpecificData(dfReproTraits, 2))
# TEST 2: Does the trait have enough data variation?
# Answer:
table(dfReproMode$repro_mode)
rareVars <- which(dfReproMode$repro_mode == "true hermaphroditism")
dfReproMode <- dfReproMode[-rareVars, ]
# Also dropping these levels from the factor.
dfReproMode$repro_mode <- droplevels(dfReproMode$repro_mode)

# TRAIT: Fertilization.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfFert <- setDT(GetTraitSpecificData(dfReproTraits, 3))
# TEST 2: Does the trait have enough data variation?
# Answer: A couple rare categories will be removed.
table(dfFert$fertilization)
rareVars <- which(dfFert$fertilization == "in mouth" | dfFert$fertilization == "other")
dfFert <- dfFert[-rareVars, ]
# Also dropping these levels from the factor.
dfFert$fertilization <- droplevels(dfFert$fertilization)

# TRAIT: Rep Guild 1.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfRepGuild1 <- setDT(GetTraitSpecificData(dfReproTraits, 4))
# TEST 2: Does the trait have enough data variation?
# Answer: Looks good.
table(dfRepGuild1$rep_guild_1)

# TRAIT: Rep Guild 2.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfRepGuild2 <- setDT(GetTraitSpecificData(dfReproTraits, 5))
# Revalue the category names so the syntax is consistent (i.e. "Polar" should be "polar").
dfRepGuild2[, rep_guild_2 := revalue(rep_guild_2, c("Brood hiders" = "brood hiders", 
                                                    "Clutch tenders" = "clutch tenders",
                                                    "External brooders" = "external brooders", 
                                                    "Nesters" = "nesters",
                                                    "Open water/substratum egg scatterers" = "open water/substratum egg scatterers"))]
# TEST 2: Does the trait have enough data variation?
# Answer: Yes.
table(dfRepGuild2$rep_guild_2)

# TRAIT: Parental care.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfParentalCare <- setDT(GetTraitSpecificData(dfReproTraits, 6))
rm(dfReproTraits)
# TEST 2: Does the trait have enough data variation?
# Answer: Yes.
table(dfParentalCare$parental_care)

# Finally, prepare the dfReproductionMV datatable by merging all univariate datatables.
dfReproductionMV <- merge(dfReproMode, dfFert, all = TRUE, by = "species_name")
dfReproductionMV <- merge(dfReproductionMV, dfRepGuild1, all = TRUE, by = "species_name")
dfReproductionMV <- merge(dfReproductionMV, dfRepGuild2, all = TRUE, by = "species_name")
dfReproductionMV <- merge(dfReproductionMV, dfParentalCare, all = TRUE, by = "species_name")



### PHYSIOLOGY RELATED ###
# Swimming.
#dfSwimming <- data.frame(swimming(speciesNames))
#write.csv(dfSwimming, file = "swimming_information.csv") 
# Read in the swimming information.
dfSwimming <- fread("swimming_information.csv")
# Datatable reorganization and renaming. Renaming column names to keep variable
# syntax consistent throughout the pipeline.
dfSwimmingTraits <- dfSwimming[, .(species_name = sciname, adult_type = AdultType, 
                                   adult_mode = AdultMode)]
# First, change all of the traits to factor type.
changeVars <- c(2:3)
dfSwimmingTraits[, (changeVars) := lapply(.SD, as.factor), .SDcols = changeVars]
rm(changeVars)
rm(speciesNames)

# TRAIT: Adult Type.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfSwimmingType <- setDT(GetTraitSpecificData(dfSwimmingTraits, 2))
# TEST 2: Does the trait have enough data variation?
# Answer: One rare category.
table(dfSwimmingType$adult_type)
# Revalue the category names so the syntax is consistent (i.e. "Polar" should be "polar").
dfSwimmingType[, adult_type := revalue(adult_type, 
                                       c("Undulation of median or pectoral fins" = "undulation of median or pectoral fins"))]

# TRAIT: Adult Mode.
# Filtering for presence of trait data. This is for the univariate analyses section.
# TEST 1: Does trait have data for at least 500 species?
# Answer: Yes!
dfSwimmingMode <- setDT(GetTraitSpecificData(dfSwimmingTraits, 3))
rm(dfSwimmingTraits)
# Revalue the category names so the syntax is consistent (i.e. "Polar" should be "polar").
dfSwimmingMode[, adult_mode := revalue(adult_mode, c("Balistiform" = "balistiform", 
                                                     "Subcarangiform" = "subcarangiform"))]
# TEST 2: Does the trait have enough data variation?
# Answer: Several rare categories.
table(dfSwimmingMode$adult_mode)
rareVars <- which(dfSwimmingMode$adult_mode == "Bathypteroiform" | 
                  dfSwimmingMode$adult_mode == "tetraodontiform")
dfSwimmingMode <- dfSwimmingMode[-rareVars, ]
rm(rareVars)
# Also dropping these levels from the factor.
dfSwimmingMode$adult_mode <- droplevels(dfSwimmingMode$adult_mode)

# Finally, prepare the dfSwimmingMV datatable by merging all univariate datatables.
dfSwimmingMV <- merge(dfSwimmingType, dfSwimmingMode, all = TRUE, by = "species_name")


# Construction of dfTraits datatable.
# This table contains all potential traits for multivariate analysis.
# Let's merge the trait information back to dfFiltered.
# NA/NULL/blank for those species that don't have info for that particular trait.
# I only want a single row per BIN for this merging process.
dfFilteredSingle <- dfFiltered[!duplicated(dfFiltered$bin_uri), ]
# Let's take the columns we need to construct the dfTraits datatable.
dfFilteredSingle <- dfFilteredSingle[, .(bin_uri, species_name, initial_bin_size, filtered_bin_size)]
# Now merge to all of the trait MV datatables.
dfTraits <- merge(dfFilteredSingle, dfLatitudeMV, all = TRUE, by = "bin_uri")
rm(dfFilteredSingle)
rm(dfLatitudeMV)
# Dataframe reorganization.
dfTraits <- dfTraits[, c(1:4, 7:8)]
colnames(dfTraits)[2] <- "species_name"
colnames(dfTraits)[4] <- "filtered_bin_size"
# Merge the FishBase traits.
dfTraits <- merge(dfTraits, dfSpeciesGenMV, all = TRUE, by = "species_name")
dfTraits <- merge(dfTraits, dfEcologyMV, all = TRUE, by = "species_name")
dfTraits <- merge(dfTraits, dfReproductionMV, all = TRUE, by = "species_name")
dfTraits <- merge(dfTraits, dfSwimmingMV, all = TRUE, by = "species_name")
dfTraits <- merge(dfTraits, dfMaturityMV, all = TRUE, by = "species_name")
dfTraits <- merge(dfTraits, dfMorphologyMV, all = TRUE, by = "species_name")
rm(dfEcologyMV)
rm(dfReproductionMV)
rm(dfSwimmingMV)
rm(dfMaturityMV)
rm(dfMorphologyMV)

# Collect the species names needed for the master phylogeny.
dfMasterSpecies <- as.data.frame(unique(dfTraits$species_name))
colnames(dfMasterSpecies)[1] <- "species_name"

# Merge back to dfFiltered to obtain all of the sequence information for each BIN.
# This is for creation of the master phylogeny.
# For univariate tree, dfTraits here instead of dfCompleteCases.
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
  #alignment1 <- lapply(DNAStringSet1, function(x) 
    #muscle::muscle(x, diags = TRUE, gapopen = -3000))
  
  # Alignment using DECIPHER.
  alignment1 <- lapply(DNAStringSet1, function(x) 
    AlignSeqs(x, iterations = 0, refinements = 0, gapOpening = -3000))

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
# This will be changed in the future but for tree building purposes I just need one
# species name for each BIN.
dfAllSeqs <- merge(aggregate(filtered_bin_size ~ species_name, data = dfAllSeqs, max), dfAllSeqs, all.x = T, sort = TRUE)
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

  # We must ensure that the sequences are of the chr type when all of the sequences 
  # PLUS the reference sequence(s) are combined into a vector. The reference 
  # sequence is added as the first sequence.
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
  
  # Writing this to file so I can perform alignment in MEGA.
  #writeXStringSet(DNAStringSet2, file = "centroidsToBeAligned.fas", 
   #               format = "fasta")

  # Run a multiple sequence alignment of all sequences including the reference 
  # using MUSCLE. This could take several minutes depending on the number of 
  # sequences and computer speed.
  gc()
  alignment2 <- muscle::muscle(DNAStringSet2, diags = TRUE, gapopen = -3000)
  #alignment2 <- AlignSeqs(DNAStringSet2, gapOpening = -3000)
  
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

  # Again, convert to DNAStringSet format.
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)
  
  # Checking alignment.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentTrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  writeXStringSet(DNAStringSet3, file = classFileNames[[1]], 
                  format = "fasta", width = 1500)

  # Remove the reference sequence from this, as we dont want it to be included in 
  # further analysis.
  refSeqRm <- which(DNAStringSet3@ranges@NAMES == "REFERENCE")
  DNAStringSet3 <- subset(DNAStringSet3[-refSeqRm])

  # Reorder dfAllSeqs according to the order of species produced by the alignment, 
  # which are now contained in the DNA_StringSet object.
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
# Here, extremely gappy/ungappy sequences are removed. These sequences are assumed
# to contribute to misalignment of the sequences or may even be pseudogenes.
# Manually checking of the alignment is recommended.

# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGaps <- sapply(regmatches(dfCheckAllSeqs$nucleotides, gregexpr("[-+]", dfCheckAllSeqs$nucleotides)), length)
  
# Mean gap length and range.
meanGap <- mean(internalGaps)
extremeHighGap <- meanGap + 7  # Upper range.
extremeLowGap <- meanGap - 7  # Lower range.
rm(meanGap)

# We then loop through each sequence to see if the number of gaps 
# deviates greatly from the mean.
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
# Make sure any outgroups are not removed.
goodBins <- which(dfExtreme$order_name == "Anura" | 
                  dfExtreme$order_name == "Rajiformes" |
                  dfExtreme$order_name == "Rodentia" |
                  dfExtreme$order_name == "Carcharhiniformes")
dfExtreme <- dfExtreme[-goodBins, ]
extremeBins <- dfExtreme$bin_uri
# If you decide to remove all from your data:
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% extremeBins, ]
rm(extremeBins)


### NEAREST NEIGHBOUR CHECK ###
# Remove centroid sequences whose nearest neighbours are in a different order or family.
# Convert each alignment to DNAbin format.
dnaBinNN <- DNAStringSet(dfCheckAllSeqs$nucleotides)
names(dnaBinNN) <- dfCheckAllSeqs$bin_uri
dnaBinNN <- as.DNAbin(dnaBinNN)
# Then, we perform genetic distance determination with the TN93 model.
geneticDistanceCentroid <- dist.dna(dnaBinNN, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

mean(geneticDistanceCentroid)
median(geneticDistanceCentroid)
range(geneticDistanceCentroid)
sd(geneticDistanceCentroid)
# IQR
lowerq = quantile(geneticDistanceCentroid)[2]
upperq = quantile(geneticDistanceCentroid)[4]
iqr = upperq - lowerq
mild.threshold.upper = (iqr * 1.5) + upperq
mild.threshold.lower = lowerq - (iqr * 1.5)

# Remove 0 values.
geneticDistanceCentroid[geneticDistanceCentroid == 0] <- NA

test <- as.data.frame(geneticDistanceCentroid)
# Identify BINs with no relatives within range of divergence.
x <- apply(test, MARGIN = 1, function(x) all(x > 0.338))
x <- which(x == "FALSE")
test2 <- test[-x, ]

# The nearest neighbour can be determined from the distance matrix alone.
# It is the sequence with minimum pairwise distance to the sequence in question.
result <- t(sapply(seq(nrow(geneticDistanceCentroid)), function(i) {
  j <- which.min(geneticDistanceCentroid[i,])
  c(paste(rownames(geneticDistanceCentroid)[i], colnames(geneticDistanceCentroid)[j], sep='/'), geneticDistanceCentroid[i, j])
}))

require(reshape)
result <- as.data.frame(result, stringsAsFactors = FALSE)
result <- transform(result, V1 = colsplit(V1, split = "/", names = c('bin_uri', 'nearest_neighbour')))
result$V2 <- as.numeric(result$V2)

# Get family and order names of BINs and nearest neighbours.
dfNN <- result$V1
dfNN$distance <- result$V2
bins <- as.character(dfNN$bin_uri)
nn <- as.character(dfNN$nearest_neighbour)
# Make dfs following these orders we can extract order and family names.
bin_ord <- dfCheckAllSeqs[match(dfNN$bin_uri, dfCheckAllSeqs$bin_uri), ] 
bin_orders <- bin_ord$order_name
bin_families <- bin_ord$family_name
nn_ord <- dfCheckAllSeqs[match(dfNN$nearest_neighbour, dfCheckAllSeqs$bin_uri), ] 
nn_orders <- nn_ord$order_name
nn_families <- nn_ord$family_name
# Add these columns to the dfNN.
dfNN$bin_order <- bin_orders
dfNN$bin_family <- bin_families
dfNN$nn_order <- nn_orders
dfNN$nn_family <- nn_families
dfNN <- setDT(dfNN)
# Non matching orders.
nonmatchOrd <- which(dfNN$bin_order != dfNN$nn_order)
# LOook closely.
dfNonmatchOrd <- dfNN[nonmatchOrd, ]
dfNonmatchOrd <- dfNonmatchOrd[which(distance < 0.05)]
# Non matching fams.
nonmatchFam <- which(dfNN$bin_family != dfNN$nn_family)
# LOook closely.
dfNonmatchFam <- dfNN[nonmatchFam, ]
dfNonmatchFam <- dfNonmatchFam[which(distance < 0.05)]
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
bin_ord <- dfCheckAllSeqs[match(closeNeighbours$bins, dfCheckAllSeqs$bin_uri), ] 
bin_orders <- bin_ord$order_name
bin_families <- bin_ord$family_name
nn_ord <- dfCheckAllSeqs[match(closeNeighbours$neighbour, dfCheckAllSeqs$bin_uri), ] 
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
# LOook closely.
#dfNonmatchOrd <- dfNN[nonmatchOrd, ]
#dfNonmatchOrd <- dfNonmatchOrd[which(distance < 0.05)]
# Non matching fams.
nonmatchFam <- which(closeNeighbours$bin_family != closeNeighbours$nn_family)
# LOook closely.
dfNonmatchFam <- closeNeighbours[nonmatchFam, ]  # Some ones I already removed!
# If you decide to remove all from your data:
dfAllSeqs <- dfAllSeqs[!dfAllSeqs$bin_uri %in% dfNonmatchFam$bins, ]


# Align and trim dfAllSeqs again without the extreme BINs and conflicted BINs.
dfCheck2AllSeqs <- refSeqTrim(dfAllSeqs)  # Check over sequences/alignment, make sure it is in correct reading frame.

# If dfCheck2AllSeqs alignment acceptable, proceed with tree building!

# Which outgroups made the cut? Remove them from MSA so I can build tree
# just using the ingroup.
outgroupSpecies <- unique(dfOutgroup$species_name)
dfGoodOutgroups <- dfAllSeqs[dfAllSeqs$species_name %in% outgroupSpecies, ]
outgroupBins <- unique(dfGoodOutgroups$species_name)
# Remove the outgroups from dfCheck2AllSeqs.
dfCheck3AllSeqs <- dfCheck2AllSeqs[!dfCheck2AllSeqs$bin_uri %in% outgroupBins, ]

# Run the alignment without the outgroups.
dfCheck4AllSeqs <- refSeqTrim(dfCheck3AllSeqs)  # Check over sequences/alignment, make sure it is in correct reading frame.

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

###########################################################################################################


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

# Prune the constraint tree so only those tips that are match with names in dfNoOutgroups remain.
prunedFishTree <- drop.tip(phy = fishTree, 
                           tip = fishTree$tip.label[!fishTree$tip.label%in%dfCheck4AllSeqs$species_name], rooted = T)
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
outgroups <- c("QUERY   Raja montagui", "QUERY   Raja polystigma", "Echinorhinus cookei", "Echinorhinus brucus")
rootedWholeTree <- root(rootedWholeTree, outgroup = outgroups, resolve.root = TRUE)



################################################################################
# TRAIT: NUMBER OF NODES
# Let's first determine the number of nodes for each species.
#phy4ML <- as(unconstrainedTree, "phylo4")
#plot(phy4ML, show.node = TRUE)
#root <- rootNode(phy4ML)
#nodeList <- lapply(1:nTips(phy4ML), function(i) .tipToRoot(phy4ML, i, root))  # A bit slow.
#numberOfNodes2 <- sapply(1:nTips(phy4ML), function(i) length(nodeList[[i]]))  # Check if these are accurate.
#names(numberOfNodes2) <- tipLabels(phy4ML)
numberOfNodes <- distRoot(rootedWholeTree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(numberOfNodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)

################################################################################
# Merge with the node_number df.
# Match the species names in data to tip labels.
dfAllSeqsNode <- merge(dfAllSeqs, dfNumberOfNodes, by = "species_name", all = TRUE)

# Let's calculate the sum of branch lengths now (from root to tip).
#is.rooted(unconstrainedTree)
#branchLengths <- distRoot(unconstrainedTree, tips = "all", method = "patristic")  # Takes a while.
branchLengths <- diag(vcv.phylo(rootedWholeTree))  # Much faster.
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
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]

# TRAIT: Branch length.
# Pagel's lambda. A test for phylogenetic signal.
branch_length <- dfRegression$branchLength
names(branch_length) <- rootedWholeTree$tip.label
sigBL <- phylosig(rootedWholeTree, branch_length, method = "lambda", test = TRUE)


# TRAIT: Number of nodes.
# Pagel's lambda. A test for phylogenetic signal.
number_of_nodes <- dfRegression$numberOfNodes
names(number_of_nodes) <- rootedWholeTree$tip.label
sigNodes <- phylosig(rootedWholeTree, number_of_nodes, method = "lambda", test = TRUE)
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
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
sigLat <- phylosig(rootedWholeTree, median_lat, method = "lambda", test = TRUE)
# Pruning tree for latitude.
# Merge so can get branch length info for UV.
dfLatitude <- merge(dfLatitude, dfRegression, by.x = "species_label", by.y = "species_name")
# Take BINs with highest number of sequences with species info in order to avoid
# using BINs that were assigned the same species label.
dfLatitude <- merge(aggregate(bin_size ~ species_label, data = dfLatitude, max), dfLatitude, all.x = T, sort = TRUE)
dup_majority_species <- which(duplicated(dfLatitude$species_label)) # Dealing with ties.
dfLatitude <- dfLatitude[-dup_majority_species,]
rm(dup_majority_species)
# Prune constrained tree so only those tips for lat are present.
latTree <- drop.tip(phy = rootedWholeTree, 
                    tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfLatitude$species_name], rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfLatitude <- dfLatitude[match(latTree$tip.label, dfLatitude$species_name), ]
dfLatitude <- as.data.frame(dfLatitude)
row.names(dfLatitude) <- dfLatitude$species_name
dfLatitude <- dfLatitude[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLatitude, "species_name")
caperLatitude <- pgls(branchLength ~ median_lat.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Maximum length.
# Pagel's lambda. A test for phylogenetic signal.
maximum_length <- dfRegression$Length
names(maximum_length) <- rootedWholeTree$tip.label
sigLength <- phylosig(rootedWholeTree, maximum_length, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfMaxLength <- merge(dfMaxLength, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
lengthTree <- drop.tip(phy = rootedWholeTree, 
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMaxLength$species_name], rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfMaxLength <- dfMaxLength[match(lengthTree$tip.label, dfMaxLength$species_name), ]
dfMaxLength <- as.data.frame(dfMaxLength)
row.names(dfMaxLength) <- dfMaxLength$species_name
dfMaxLength <- dfMaxLength[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfMaxLength, "species_name")
caperLength <- pgls(branchLength ~ Length.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Maximum weight.
# Pagel's lambda. A test for phylogenetic signal.
maximum_weight <- dfRegression$Weight
names(maximum_weight) <- rootedWholeTree$tip.label
sigWeight <- phylosig(rootedWholeTree, maximum_weight, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfWeight <- merge(dfWeight, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
weightTree <- drop.tip(phy = rootedWholeTree, 
                     tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfWeight$species_name], rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfWeight <- dfWeight[match(weightTree$tip.label, dfWeight$species_name), ]
dfWeight <- as.data.frame(dfWeight)
row.names(dfWeight) <- dfWeight$species_name
dfWeight <- dfWeight[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfWeight, "species_name")
caperWeight <- pgls(branchLength ~ Weight.x + numberOfNodes, c_data, lambda = "ML")


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
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfLongWild$species_name], rooted = T)
tiph <- diag(vcv.phylo(longevityTree))
dfLongWild <- dfLongWild[match(longevityTree$tip.label, dfLongWild$species_name), ]
dfLongWild <- as.data.frame(dfLongWild)
row.names(dfLongWild) <- dfLongWild$species_name
dfLongWild <- dfLongWild[, c(1:3, 5)]
dfLongWild <- na.omit(dfLongWild)
fitLong <- gls(branchLength ~ LongevityWild.x + numberOfNodes, correlation=corPagel(value = 0, phy = longevityTree),
           weights=varFixed(~tiph), method = "ML", data = dfLongWild)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfLongWild, "species_name")
caperLong <- pgls(branchLength ~ LongevityWild.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Temp min.
# Pagel's lambda. A test for phylogenetic signal.
temp_min <- dfRegression$TempMin
names(temp_min) <- rootedWholeTree$tip.label
sigTempMin <- phylosig(rootedWholeTree, temp_min, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfTempMin <- merge(dfTempMin, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
tempMinTree <- drop.tip(phy = rootedWholeTree, 
                          tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfTempMin$species_name], rooted = T)
tiph <- diag(vcv.phylo(tempMinTree))
dfTempMin <- dfTempMin[match(tempMinTree$tip.label, dfTempMin$species_name), ]
dfTempMin <- as.data.frame(dfTempMin)
row.names(dfTempMin) <- dfTempMin$species_name
dfTempMin <- dfTempMin[, c(1:3, 5)]
fitTempMin <- gls(branchLength.x ~ TempMin.x + numberOfNodes.x, correlation=corPagel(value = 0, phy = tempMinTree),
               weights=varFixed(~tiph), method = "ML", data = dfTempMin)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfTempMin, "species_name")
caperTempMin <- pgls(branchLength.x ~ TempMin.x + numberOfNodes.x, c_data, lambda = "ML")


# TRAIT: Temp max.
# Pagel's lambda. A test for phylogenetic signal.
temp_max <- dfRegression$TempMax
names(temp_max) <- rootedWholeTree$tip.label
sigTempMax <- phylosig(rootedWholeTree, temp_max, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfTempMax <- merge(dfTempMax, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
tempMaxTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfTempMax$species_name], rooted = T)
tiph <- diag(vcv.phylo(tempMaxTree))
dfTempMax <- dfTempMax[match(tempMaxTree$tip.label, dfTempMax$species_name), ]
dfTempMax <- as.data.frame(dfTempMax)
row.names(dfTempMax) <- dfTempMax$species_name
dfTempMax <- dfTempMax[, c(1:3, 5)]
fitTempMax <- gls(branchLength.x ~ TempMax.x + numberOfNodes.x, correlation=corPagel(value = 0, phy = tempMaxTree),
                  weights=varFixed(~tiph), method = "ML", data = dfTempMax)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfTempMax, "species_name")
caperTempMax <- pgls(branchLength.x ~ TempMax.x + numberOfNodes.x, c_data, lambda = "ML")



# TRAIT: Minimum depth
# Pagel's lambda. A test for phylogenetic signal.
shallow_depth <- dfRegression$DepthRangeShallow
names(shallow_depth) <- rootedWholeTree$tip.label
sigShallowDepth <- phylosig(rootedWholeTree, shallow_depth, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfDepthShallow <- merge(dfDepthShallow, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
shallowTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfDepthShallow$species_name], rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfDepthShallow <- dfDepthShallow[match(shallowTree$tip.label, dfDepthShallow$species_name), ]
dfDepthShallow <- as.data.frame(dfDepthShallow)
row.names(dfDepthShallow) <- dfDepthShallow$species_name
dfDepthShallow <- dfDepthShallow[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDepthShallow, "species_name")
caperShallow <- pgls(branchLength ~ DepthRangeShallow.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Maximum depth
# Pagel's lambda. A test for phylogenetic signal.
deep_depth <- dfRegression$DepthRangeDeep
names(deep_depth) <- rootedWholeTree$tip.label
sigDeepDepth <- phylosig(rootedWholeTree, shallow_depth, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfDepthDeep <- merge(dfDepthDeep, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
deepTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfDepthDeep$species_name], rooted = T)
# Prune constrained tree so only those tips for lat are present.
dfDepthDeep <- dfDepthDeep[match(deepTree$tip.label, dfDepthDeep$species_name), ]
dfDepthDeep <- as.data.frame(dfDepthDeep)
row.names(dfDepthDeep) <- dfDepthDeep$species_name
dfDepthDeep <- dfDepthDeep[, c(1:3, 5)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDepthDeep, "species_name")
caperDeep <- pgls(branchLength ~ DepthRangeDeep.x + numberOfNodes, c_data, lambda = "ML")

# TRAIT: Diet Troph
# Pagel's lambda. A test for phylogenetic signal.
diet_troph <- dfRegression$DietTroph
names(diet_troph) <- rootedWholeTree$tip.label
sigDiet <- phylosig(rootedWholeTree, diet_troph, method = "lambda", test = TRUE)
# Pruning tree.
# Merge so can get branch length info for UV.
dfDietTroph <- merge(dfDietTroph, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
dietTree <- drop.tip(phy = rootedWholeTree, 
                    tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfDietTroph$species_name], rooted = T)
tiph <- diag(vcv.phylo(dietTree))
dfDietTroph <- dfDietTroph[match(dietTree$tip.label, dfDietTroph$species_name), ]
dfDietTroph <- as.data.frame(dfDietTroph)
row.names(dfDietTroph) <- dfDietTroph$species_name
dfDietTroph <- dfDietTroph[, c(1:3, 5)]
fitDiet <- gls(branchLength ~ DietTroph.x + numberOfNodes, correlation=corPagel(value = 0, phy = dietTree),
              weights=varFixed(~tiph), method = "ML", data = dfDietTroph)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDietTroph, "species_name")
caperDiet <- pgls(branchLength ~ DietTroph.x + numberOfNodes, c_data, lambda = "ML")


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
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfAgeMaturity$species_name], rooted = T)
tiph <- diag(vcv.phylo(ageTree))
dfAgeMaturity <- dfAgeMaturity[match(ageTree$tip.label, dfAgeMaturity$species_name), ]
dfAgeMaturity <- as.data.frame(dfAgeMaturity)
row.names(dfAgeMaturity) <- dfAgeMaturity$species_name
dfAgeMaturity <- dfAgeMaturity[, c(1:3, 5)]
fitAge <- gls(branchLength ~ age_at_maturity.x + numberOfNodes, correlation=corPagel(value = 0, phy = ageTree),
                  weights=varFixed(~tiph), method = "ML", data = dfAgeMaturity)
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfAgeMaturity, "species_name")
caperAge <- pgls(branchLength ~ age_at_maturity.x + numberOfNodes, c_data, lambda = "ML")



# Discrete traits.

# TRAIT: NERITIC
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfNeritic <- merge(dfNeritic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
neriticTree <- drop.tip(phy = rootedWholeTree, 
                           tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfNeritic$species_name], rooted = T)
dfNeritic <- dfNeritic[match(neriticTree$tip.label, dfNeritic$species_name), ]
dfNeritic <- as.data.frame(dfNeritic)
row.names(dfNeritic) <- dfNeritic$species_name
dfNeritic <- dfNeritic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfNeritic$Neritic.x <- relevel(dfNeritic$Neritic.x, ref = "0")
dfRegression$Neritic <- relevel(dfRegression$Neritic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfNeritic, "species_name")
caperNeritic <- pgls(branchLength ~ Neritic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigNer <- phylo.d(dfNeritic, neriticTree, names.col = species_name, binvar = Neritic.x)
#fit$opt$lambda


# TRAIT: OCEANIC
#oceanic <- dfRegression$Oceanic
#names(oceanic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, oceanic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfOceanic <- merge(dfOceanic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
oceanicTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfOceanic$species_name], rooted = T)
dfOceanic <- dfOceanic[match(oceanicTree$tip.label, dfOceanic$species_name), ]
dfOceanic <- as.data.frame(dfOceanic)
row.names(dfOceanic) <- dfOceanic$species_name
dfOceanic <- dfOceanic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfOceanic$Oceanic.x <- relevel(dfOceanic$Oceanic.x, ref = "0")
dfRegression$Oceanic <- relevel(dfRegression$Oceanic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfOceanic, "species_name")
caperOceanic <- pgls(branchLength ~ Oceanic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigOceanic <- phylo.d(dfOceanic, oceanicTree, names.col = species_name, binvar = Oceanic.x)


# TRAIT: BENTHIC
#benthic <- dfRegression$Benthic
#names(benthic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, benthic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBenthic <- merge(dfBenthic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
benthicTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfBenthic$species_name], rooted = T)
dfBenthic <- dfBenthic[match(benthicTree$tip.label, dfBenthic$species_name), ]
dfBenthic <- as.data.frame(dfBenthic)
row.names(dfBenthic) <- dfBenthic$species_name
dfBenthic <- dfBenthic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBenthic$Benthic.x <- relevel(dfBenthic$Benthic.x, ref = "0")
dfRegression$Benthic <- relevel(dfRegression$Benthic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfBenthic, "species_name")
caperBenthic <- pgls(branchLength ~ Benthic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigBenthic <- phylo.d(dfBenthic, benthicTree, names.col = species_name, binvar = Benthic.x)


# TRAIT: CORAL REEFS
#coral_reefs <- dfRegression$CoralReefs
#names(coral_reefs) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, coral_reefs, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfCoralReefs <- merge(dfCoralReefs, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
coralReefsTree <- drop.tip(phy = rootedWholeTree, 
                        tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfCoralReefs$species_name], rooted = T)
dfCoralReefs <- dfCoralReefs[match(coralReefsTree$tip.label, dfCoralReefs$species_name), ]
dfCoralReefs <- as.data.frame(dfCoralReefs)
row.names(dfCoralReefs) <- dfCoralReefs$species_name
dfCoralReefs <- dfCoralReefs[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfCoralReefs$CoralReefs.x <- relevel(dfCoralReefs$CoralReefs.x, ref = "0")
dfRegression$CoralReefs <- relevel(dfRegression$CoralReefs, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfCoralReefs, "species_name")
caperCoralReefs <- pgls(branchLength ~ CoralReefs.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigCoral <- phylo.d(dfCoralReefs, coralReefsTree, names.col = species_name, binvar = CoralReefs.x)


# TRAIT: EPIPELAGIC
#epipelagic <- dfRegression$Epipelagic
#names(epipelagic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, epipelagic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfEpipelagic <- merge(dfEpipelagic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
epiTree <- drop.tip(phy = rootedWholeTree, 
                           tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfEpipelagic$species_name], rooted = T)
dfEpipelagic <- dfEpipelagic[match(epiTree$tip.label, dfEpipelagic$species_name), ]
dfEpipelagic <- as.data.frame(dfEpipelagic)
row.names(dfEpipelagic) <- dfEpipelagic$species_name
dfEpipelagic <- dfEpipelagic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfEpipelagic$Epipelagic.x <- relevel(dfEpipelagic$Epipelagic.x, ref = "0")
dfRegression$Epipelagic <- relevel(dfRegression$Epipelagic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfEpipelagic, "species_name")
caperEpi <- pgls(branchLength ~ Epipelagic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigEpipelagic <- phylo.d(dfEpipelagic, epiTree, names.col = species_name, binvar = Epipelagic.x)


# TRAIT: MESOPELAGIC
#mesopelagic <- dfRegression$Mesopelagic
#names(epipelagic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, epipelagic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfMesopelagic <- merge(dfMesopelagic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
mesoTree <- drop.tip(phy = rootedWholeTree, 
                    tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMesopelagic$species_name], rooted = T)
dfMesopelagic <- dfMesopelagic[match(mesoTree$tip.label, dfMesopelagic$species_name), ]
dfMesopelagic <- as.data.frame(dfMesopelagic)
row.names(dfMesopelagic) <- dfMesopelagic$species_name
dfMesopelagic <- dfMesopelagic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfMesopelagic$Mesopelagic.x <- relevel(dfMesopelagic$Mesopelagic.x, ref = "0")
dfRegression$Mesopelagic <- relevel(dfRegression$Mesopelagic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfMesopelagic, "species_name")
caperMeso <- pgls(branchLength ~ Mesopelagic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigMeso <- phylo.d(dfMesopelagic, mesoTree, names.col = species_name, binvar = Mesopelagic.x)



# TRAIT: BATHYPELAGIC
#bathypelagic <- dfRegression$Bathypelagic
#names(bathypelagic) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, bathypelagic, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfBathypelagic <- merge(dfBathypelagic, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
bathyTree <- drop.tip(phy = rootedWholeTree, 
                     tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfBathypelagic$species_name], rooted = T)
dfBathypelagic <- dfBathypelagic[match(bathyTree$tip.label, dfBathypelagic$species_name), ]
dfBathypelagic <- as.data.frame(dfBathypelagic)
row.names(dfBathypelagic) <- dfBathypelagic$species_name
dfBathypelagic <- dfBathypelagic[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfBathypelagic$Bathypelagic.x <- relevel(dfBathypelagic$Bathypelagic.x, ref = "0")
dfRegression$Bathypelagic <- relevel(dfRegression$Bathypelagic, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfBathypelagic, "species_name")
caperBathy <- pgls(branchLength ~ Bathypelagic.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigBathy <- phylo.d(dfBathypelagic, bathyTree, names.col = species_name, binvar = Bathypelagic.x)

# TRAIT: MANGROVES
#mangroves <- dfRegression$Mangroves
#names(mangroves) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, mangroves, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfMangroves <- merge(dfMangroves, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
mangrovesTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfMangroves$species_name], rooted = T)
dfMangroves <- dfMangroves[match(mangrovesTree$tip.label, dfMangroves$species_name), ]
dfMangroves <- as.data.frame(dfMangroves)
row.names(dfMangroves) <- dfMangroves$species_name
dfMangroves <- dfMangroves[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfMangroves$Mangroves.x <- relevel(dfMangroves$Mangroves.x, ref = "0")
dfRegression$Mangroves <- relevel(dfRegression$Mangroves, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfMangroves, "species_name")
caperMangroves <- pgls(branchLength ~ Mangroves.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigMangroves <- phylo.d(dfMangroves, mangrovesTree, names.col = species_name, binvar = Mangroves.x)



# TRAIT: ESTUARIES
#estuaries <- dfRegression$Estuaries
#names(estuaries) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, estuaries, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfEstuaries <- merge(dfEstuaries, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
estuariesTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfEstuaries$species_name], rooted = T)
dfEstuaries <- dfEstuaries[match(estuariesTree$tip.label, dfEstuaries$species_name), ]
dfEstuaries <- as.data.frame(dfEstuaries)
row.names(dfEstuaries) <- dfEstuaries$species_name
dfEstuaries <- dfEstuaries[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfEstuaries$Estuaries.x <- relevel(dfEstuaries$Estuaries.x, ref = "0")
dfRegression$Estuaries <- relevel(dfRegression$Estuaries, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfEstuaries, "species_name")
caperEstuaries <- pgls(branchLength ~ Estuaries.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigEstuaries <- phylo.d(dfEstuaries, estuariesTree, names.col = species_name, binvar = Estuaries.x)


# TRAIT: SCHOOLING
#schooling <- dfRegression$Schooling
#names(schooling) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, schooling, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfSchooling <- merge(dfSchooling, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
schoolingTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfSchooling$species_name], rooted = T)
dfSchooling <- dfSchooling[match(schoolingTree$tip.label, dfSchooling$species_name), ]
dfSchooling <- as.data.frame(dfSchooling)
row.names(dfSchooling) <- dfSchooling$species_name
dfSchooling <- dfSchooling[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfSchooling$Schooling.x <- relevel(dfSchooling$Schooling.x, ref = "0")
dfRegression$Schooling <- relevel(dfRegression$Schooling, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfSchooling, "species_name")
caperSchooling <- pgls(branchLength ~ Schooling.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigSchool <- phylo.d(dfSchooling, schoolingTree, names.col = species_name, binvar = Schooling.x)



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







# TRAIT: Freshwater
#freshwater <- dfRegression$Fresh
#names(freshwater) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, freshwater, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfFreshwater <- merge(dfFreshwater, dfRegression, by = "species_name")
dfFreshwater <- dfFreshwater[sample(nrow(dfFreshwater), 4393), ]
# Prune constrained tree so only those tips for lat are present.
freshTree <- drop.tip(phy = rootedWholeTree, 
                       tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfFreshwater$species_name], rooted = T)
dfFreshwater <- dfFreshwater[match(freshTree$tip.label, dfFreshwater$species_name), ]
dfFreshwater <- as.data.frame(dfFreshwater)
row.names(dfFreshwater) <- dfFreshwater$species_name
dfFreshwater <- dfFreshwater[, c(1:2, 5, 7)]
# Take random sample.
# Relevel body shape as I want fusiform / normal to be the reference.
dfFreshwater$Fresh.x <- relevel(dfFreshwater$Fresh.x, ref = "0")
dfRegression$Fresh <- relevel(dfRegression$Fresh, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfFreshwater, "species_name")
caperFresh <- pgls(branchLength ~ Fresh.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigFresh <- phylo.d(dfFreshwater, freshTree, names.col = species_name, binvar = Fresh.x)


# TRAIT: Saltwater
#saltwater <- dfRegression$Saltwater
#names(saltwater) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, saltwater, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfSaltwater <- merge(dfSaltwater, dfRegression, by = "species_name")
dfSaltwater <- dfSaltwater[sample(nrow(dfSaltwater), 4393), ]
# Prune constrained tree so only those tips for lat are present.
saltTree <- drop.tip(phy = rootedWholeTree, 
                      tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfSaltwater$species_name], rooted = T)
dfSaltwater <- dfSaltwater[match(saltTree$tip.label, dfSaltwater$species_name), ]
dfSaltwater <- as.data.frame(dfSaltwater)
row.names(dfSaltwater) <- dfSaltwater$species_name
dfSaltwater <- dfSaltwater[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfSaltwater$Saltwater.x <- relevel(dfSaltwater$Saltwater.x, ref = "0")
dfRegression$Saltwater <- relevel(dfRegression$Saltwater, ref = "0")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfSaltwater, "species_name")
caperSalt <- pgls(branchLength ~ Saltwater.x + numberOfNodes, c_data, lambda = "ML")
# D metric for phylogenetic signal.
sigSalt <- phylo.d(dfSaltwater, saltTree, names.col = species_name, binvar = Saltwater.x)



# TRAIT: Env temp
# get rid of deep-water, too confusing
deep_water <- which(dfEnvTemp$EnvTemp == "deep-water")
dfEnvTemp <- dfEnvTemp[-deep_water,]
dfEnvTemp$EnvTemp <- droplevels(dfEnvTemp$EnvTemp)
#env_temp <- dfRegression$EnvTemp
#names(env_temp) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, env_temp, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfEnvTemp <- merge(dfEnvTemp, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
envTempTree <- drop.tip(phy = rootedWholeTree, 
                     tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfEnvTemp$species_name], rooted = T)
dfEnvTemp <- dfEnvTemp[match(envTempTree$tip.label, dfEnvTemp$species_name), ]
dfEnvTemp <- as.data.frame(dfEnvTemp)
row.names(dfEnvTemp) <- dfEnvTemp$species_name
dfEnvTemp <- dfEnvTemp[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfEnvTemp$EnvTemp.x <- relevel(dfEnvTemp$EnvTemp.x, ref = "polar")
dfRegression$EnvTemp <- relevel(dfRegression$EnvTemp, ref = "polar")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfEnvTemp, "species_name")
caperEnvTemp <- pgls(branchLength ~ EnvTemp.x + numberOfNodes, c_data, lambda = "ML")
anovaEnvTemp <- anova.pgls.fixed(caperEnvTemp)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyEnvTemp <- model.matrix(~ EnvTemp.x - 1, data=dfEnvTemp)
row.names(dfDummyEnvTemp) <- dfEnvTemp$species_name
dfDummyEnvTemp <- as.data.frame(dfDummyEnvTemp)
dfDummyEnvTemp$species_name <- row.names(dfDummyEnvTemp)
colnames(dfDummyEnvTemp)[1] <- "Polar"
colnames(dfDummyEnvTemp)[2] <- "Subtropical"
colnames(dfDummyEnvTemp)[3] <- "Temperate"
colnames(dfDummyEnvTemp)[4] <- "Tropical"
sigPolar <- phylo.d(dfDummyEnvTemp, envTempTree, names.col = species_name, binvar = Polar)
sigTemperate <- phylo.d(dfDummyEnvTemp, envTempTree, names.col = species_name, binvar = Temperate)
sigSubTrop <- phylo.d(dfDummyEnvTemp, envTempTree, names.col = species_name, binvar = Subtropical)
sigTrop <- phylo.d(dfDummyEnvTemp, envTempTree, names.col = species_name, binvar = Tropical)





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



# Repro mode
#repro_mode <- dfRegression$repro_mode
#names(repro_mode) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, repro_mode, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfReproMode <- merge(dfReproMode, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
reproModeTree <- drop.tip(phy = rootedWholeTree, 
                            tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfReproMode$species_name], rooted = T)
dfReproMode <- dfReproMode[match(reproModeTree$tip.label, dfReproMode$species_name), ]
dfReproMode <- as.data.frame(dfReproMode)
row.names(dfReproMode) <- dfReproMode$species_name
dfReproMode <- dfReproMode[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfReproMode$repro_mode.x <- relevel(dfReproMode$repro_mode.x, ref = "dioecism")
dfRegression$repro_mode <- relevel(dfRegression$repro_mode, ref = "dioecism")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfReproMode, "species_name")
caperReproMode <- pgls(branchLength ~ repro_mode.x + numberOfNodes, c_data, lambda = "ML")
anovaReproMode <- anova.pgls.fixed(caperReproMode)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyRM <- model.matrix(~ repro_mode.x - 1, data=dfReproMode)
row.names(dfDummyRM) <- dfReproMode$species_name
dfDummyRM <- as.data.frame(dfDummyRM)
dfDummyRM$species_name <- row.names(dfDummyRM)
colnames(dfDummyRM)[1] <- "Dioecism"
colnames(dfDummyRM)[2] <- "Protandry"
colnames(dfDummyRM)[3] <- "Protogyny"
sigDio <- phylo.d(dfDummyRM, reproModeTree, names.col = species_name, binvar = Dioecism)
sigProta <- phylo.d(dfDummyRM, reproModeTree, names.col = species_name, binvar = Protandry)
sigProto <- phylo.d(dfDummyRM, reproModeTree, names.col = species_name, binvar = Protogyny)






# Fertilization
#fert <- dfRegression$fertilization
#names(fert) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, fert, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfFert <- merge(dfFert, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
fertTree <- drop.tip(phy = rootedWholeTree, 
                          tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfFert$species_name], rooted = T)
dfFert <- dfFert[match(fertTree$tip.label, dfFert$species_name), ]
dfFert <- as.data.frame(dfFert)
row.names(dfFert) <- dfFert$species_name
dfFert <- dfFert[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfFert$fertilization.x <- relevel(dfFert$fertilization.x, ref = "external")
dfRegression$fertilization <- relevel(dfRegression$fertilization, ref = "external")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfFert, "species_name")
caperFert <- pgls(branchLength ~ fertilization.x + numberOfNodes, c_data, lambda = "ML")
anovaFert <- anova.pgls.fixed(caperFert)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyFert <- model.matrix(~ fertilization.x - 1, data=dfFert)
row.names(dfDummyFert) <- dfFert$species_name
dfDummyFert <- as.data.frame(dfDummyFert)
dfDummyFert$species_name <- row.names(dfDummyFert)
colnames(dfDummyFert)[1] <- "External"
colnames(dfDummyFert)[2] <- "BroodPouch"
colnames(dfDummyFert)[3] <- "Internal"
sigExt <- phylo.d(dfDummyFert, fertTree, names.col = species_name, binvar = External)
sigBrood <- phylo.d(dfDummyFert, fertTree, names.col = species_name, binvar = BroodPouch)
sigInt <- phylo.d(dfDummyFert, fertTree, names.col = species_name, binvar = Internal)






# Reproductive guild I
#rep_guild_1 <- dfRegression$rep_guild_1
#names(rep_guild_1) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, rep_guild_1, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfRepGuild1 <- merge(dfRepGuild1, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
repGuild1Tree <- drop.tip(phy = rootedWholeTree, 
                     tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfRepGuild1$species_name], rooted = T)
dfRepGuild1 <- dfRepGuild1[match(repGuild1Tree$tip.label, dfRepGuild1$species_name), ]
dfRepGuild1 <- as.data.frame(dfRepGuild1)
row.names(dfRepGuild1) <- dfRepGuild1$species_name
dfRepGuild1 <- dfRepGuild1[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfRepGuild1$rep_guild_1.x <- relevel(dfRepGuild1$rep_guild_1.x, ref = "nonguarders")
dfRegression$rep_guild_1 <- relevel(dfRegression$rep_guild_1, ref = "nonguarders")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfRepGuild1, "species_name")
caperRepGuild1 <- pgls(branchLength ~ rep_guild_1.x + numberOfNodes, c_data, lambda = "ML")
anovaRepGuild1 <- anova.pgls.fixed(caperRepGuild1)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyRP1 <- model.matrix(~ rep_guild_1.x - 1, data=dfRepGuild1)
row.names(dfDummyRP1) <- dfRepGuild1$species_name
dfDummyRP1 <- as.data.frame(dfDummyRP1)
dfDummyRP1$species_name <- row.names(dfDummyRP1)
colnames(dfDummyRP1)[1] <- "Nonguarders"
colnames(dfDummyRP1)[2] <- "Bearers"
colnames(dfDummyRP1)[3] <- "Guarders"
sigNG <- phylo.d(dfDummyRP1, repGuild1Tree, names.col = species_name, binvar = Nonguarders)
sigBear <- phylo.d(dfDummyRP1, repGuild1Tree, names.col = species_name, binvar = Bearers)
sigG <- phylo.d(dfDummyRP1, repGuild1Tree, names.col = species_name, binvar = Guarders)



# Reproductive guild II
#rep_guild_2 <- dfRegression$rep_guild_2
#names(rep_guild_2) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, rep_guild_2, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfRepGuild2 <- merge(dfRepGuild2, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
repGuild2Tree <- drop.tip(phy = rootedWholeTree, 
                          tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfRepGuild2$species_name], rooted = T)
dfRepGuild2 <- dfRepGuild2[match(repGuild2Tree$tip.label, dfRepGuild2$species_name), ]
dfRepGuild2 <- as.data.frame(dfRepGuild2)
row.names(dfRepGuild2) <- dfRepGuild2$species_name
dfRepGuild2 <- dfRepGuild2[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfRepGuild2$rep_guild_2.x <- relevel(dfRepGuild2$rep_guild_2.x, ref = "open water/substratum egg scatterers")
dfRegression$rep_guild_2 <- relevel(dfRegression$rep_guild_2, ref = "open water/substratum egg scatterers")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfRepGuild2, "species_name")
caperRepGuild2 <- pgls(branchLength ~ rep_guild_2.x + numberOfNodes, c_data, lambda = "ML")
anovaRepGuild2 <- anova.pgls.fixed(caperRepGuild2)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyRP2 <- model.matrix(~ rep_guild_2.x - 1, data=dfRepGuild2)
row.names(dfDummyRP2) <- dfRepGuild2$species_name
dfDummyRP2 <- as.data.frame(dfDummyRP2)
dfDummyRP2$species_name <- row.names(dfDummyRP2)
colnames(dfDummyRP2)[1] <- "OpenWater"
colnames(dfDummyRP2)[2] <- "BroodHiders"
colnames(dfDummyRP2)[3] <- "ClutchTenders"
colnames(dfDummyRP2)[4] <- "ExtBrooders"
colnames(dfDummyRP2)[5] <- "IntLive"
colnames(dfDummyRP2)[6] <- "Nesters"
sigOW <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = OpenWater)
sigBH <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = BroodHiders)
sigCT <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = ClutchTenders)
sigEB <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = ExtBrooders)
sigIL <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = IntLive)
sigN <- phylo.d(dfDummyRP2, repGuild2Tree, names.col = species_name, binvar = Nesters)



# Parental care
#parental_care <- dfRegression$parental_care
#names(parental_care) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, parental_care, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfParentalCare <- merge(dfParentalCare, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
parentalCareTree <- drop.tip(phy = rootedWholeTree, 
                          tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfParentalCare$species_name], rooted = T)
dfParentalCare <- dfParentalCare[match(parentalCareTree$tip.label, dfParentalCare$species_name), ]
dfParentalCare <- as.data.frame(dfParentalCare)
row.names(dfParentalCare) <- dfParentalCare$species_name
dfParentalCare <- dfParentalCare[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfParentalCare$parental_care.x <- relevel(dfParentalCare$parental_care.x, ref = "none")
dfRegression$parental_care <- relevel(dfRegression$parental_care, ref = "none")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfParentalCare, "species_name")
caperParentalCare <- pgls(branchLength ~ numberOfNodes + parental_care.x, c_data, lambda = "ML")
anovaParentalCare <- anova.pgls.fixed(caperParentalCare)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyPC <- model.matrix(~ parental_care.x - 1, data=dfParentalCare)
row.names(dfDummyPC) <- dfParentalCare$species_name
dfDummyPC <- as.data.frame(dfDummyPC)
dfDummyPC$species_name <- row.names(dfDummyPC)
colnames(dfDummyPC)[1] <- "None"
colnames(dfDummyPC)[2] <- "Biparental"
colnames(dfDummyPC)[3] <- "Maternal"
colnames(dfDummyPC)[4] <- "Paternal"
sigNone <- phylo.d(dfDummyPC, parentalCareTree, names.col = species_name, binvar = None)
sigBP <- phylo.d(dfDummyPC, parentalCareTree, names.col = species_name, binvar = Biparental)
sigM <- phylo.d(dfDummyPC, parentalCareTree, names.col = species_name, binvar = Maternal)
sigP <- phylo.d(dfDummyPC, parentalCareTree, names.col = species_name, binvar = Paternal)



# Adult type swimming
#adult_type <- dfRegression$adult_type
#names(adult_type) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, adult_type, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfSwimmingType <- merge(dfSwimmingType, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
swimmingTypeTree <- drop.tip(phy = rootedWholeTree, 
                             tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfSwimmingType$species_name], rooted = T)
dfSwimmingType <- dfSwimmingType[match(swimmingTypeTree$tip.label, dfSwimmingType$species_name), ]
dfSwimmingType <- as.data.frame(dfSwimmingType)
row.names(dfSwimmingType) <- dfSwimmingType$species_name
dfSwimmingType <- dfSwimmingType[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfSwimmingType$adult_type.x <- relevel(dfSwimmingType$adult_type.x, ref = "movements of body and/or caudal fin")
dfRegression$adult_type <- relevel(dfRegression$adult_type, ref = "movements of body and/or caudal fin")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfSwimmingType, "species_name")
caperSwimmingType <- pgls(branchLength ~ adult_type.x + numberOfNodes, c_data, lambda = "ML")
anovaSwimmingType <- anova.pgls.fixed(caperSwimmingType)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummySwim <- model.matrix(~ adult_type.x - 1, data = dfSwimmingType)
row.names(dfDummySwim) <- dfSwimmingType$species_name
dfDummySwim <- as.data.frame(dfDummySwim)
dfDummySwim$species_name <- row.names(dfDummySwim)
colnames(dfDummySwim)[1] <- "Movements"
colnames(dfDummySwim)[2] <- "Oscillation"
colnames(dfDummySwim)[3] <- "Undulation"
sigMove <- phylo.d(dfDummySwim, swimmingTypeTree, names.col = species_name, binvar = Movements)
sigOsc <- phylo.d(dfDummySwim, swimmingTypeTree, names.col = species_name, binvar = Oscillation)
sigUnd <- phylo.d(dfDummySwim, swimmingTypeTree, names.col = species_name, binvar = Undulation)






# Position of mouth
#pos_of_mouth <- dfRegression$pos_of_mouth
#names(pos_of_mouth) <- rootedWholeTree$tip.label
#fit <- fitDiscrete(rootedWholeTree, pos_of_mouth, transform = "lambda")
#fit$opt$lambda
# Make sure the order of the data matches the constrained tree.
dfRegression <- dfRegression[match(rootedWholeTree$tip.label, dfRegression$species_name), ]
# Merge so can get branch length info for UV.
dfPosMouth <- merge(dfPosMouth, dfRegression, by = "species_name")
# Prune constrained tree so only those tips for lat are present.
posMouthTree <- drop.tip(phy = rootedWholeTree, 
                             tip = rootedWholeTree$tip.label[!rootedWholeTree$tip.label%in%dfPosMouth$species_name], rooted = T)
dfPosMouth <- dfPosMouth[match(posMouthTree$tip.label, dfPosMouth$species_name), ]
dfPosMouth <- as.data.frame(dfPosMouth)
row.names(dfPosMouth) <- dfPosMouth$species_name
dfPosMouth <- dfPosMouth[, c(1:3, 5)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfPosMouth$pos_of_mouth <- relevel(dfPosMouth$pos_of_mouth, ref = "terminal")
dfRegression$pos_of_mouth <- relevel(dfRegression$pos_of_mouth, ref = "terminal")
# Try using caper.
c_data <- comparative.data(rootedWholeTree, dfPosMouth, "species_name")
caperPosMouth <- pgls(branchLength ~ pos_of_mouth + numberOfNodes, c_data, lambda = "ML")
anovaPosMouth <- anova.pgls.fixed(caperPosMouth)
# D metric for phylogenetic signal.
# Need dummy variables first.
dfDummyMouth <- model.matrix(~ pos_of_mouth - 1, data=dfPosMouth)
row.names(dfDummyMouth) <- dfPosMouth$species_name
dfDummyMouth <- as.data.frame(dfDummyMouth)
dfDummyMouth$species_name <- row.names(dfDummyMouth)
colnames(dfDummyMouth)[1] <- "Inferior"
colnames(dfDummyMouth)[2] <- "Superior"
colnames(dfDummyMouth)[3] <- "Terminal"
sigInferior <- phylo.d(dfDummyMouth, posMouthTree, names.col = species_name, binvar = Inferior)
sigTerminal <- phylo.d(dfDummyMouth, posMouthTree, names.col = species_name, binvar = Superior)
sigSuperior <- phylo.d(dfDummyMouth, posMouthTree, names.col = species_name, binvar = Terminal)






# Plot significant traits.
dfPlotLength <- dfMaxLength[sample(nrow(dfMaxLength), 100), ]
dfPlotLength$Scaled <- scale(dfPlotLength$Length.x, scale = FALSE)
plot(branchLength ~ numberOfNodes, data = dfPlotLength)
# Make sure the order of the data matches the constrained tree.
c_data <- comparative.data(rootedWholeTree, dfPlotLength, "species_name")
plotlength <- lm(branchLength ~ Length.x, data=dfPlotLength)
####
abline(plotlength)



# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfRegression[, c(1:2, 4, 13:14, 46, 5, 23, 46)]
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

hist(dfMultivariateCut$age_at_maturity)
range(dfMultivariateCut$age_at_maturity)

hist(dfMultivariateCut$LongevityWild)
range(dfMultivariateCut$LongevityWild)


lm.res=lm(branchLength ~ Length + DepthRangeDeep + median_lat
          + age_at_maturity, data = dfMultivariateCut)
library(car)
vif(mod=lm.res)  # Longevity removed.

# Backward selection using BIC.
c_data <- comparative.data(rootedWholeTree, dfMultivariateCut, names.col = "species_name", vcv=TRUE)

# Full model.
full <- pgls(branchLength ~ numberOfNodes + age_at_maturity + LongevityWild, data=c_data, lambda="ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x=fitted(full), y=full$phyres, pch=5)

# Test for effect of lat.
x <- pgls(branchLength ~ numberOfNodes + median_lat + DepthRangeDeep + Length
          + age_at_maturity, 
              data=c_data, lambda="ML")
BIC(full, lat)
# Confounding?
summary(full)$coefficients
summary(lat)$coefficients

# Test for effect of longevity.
long <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep + Length
            + age_at_maturity, 
             data=c_data, lambda="ML")
BIC(lat, long)
# Confounding? Seems to be, on length.
summary(lat)$coefficients
summary(long)$coefficients

# Test for effect of body shape.
bs <- pgls(branchLength ~ numberOfNodes + Length + DepthRangeDeep
              + age_at_maturity + Saltwater, 
               data=c_data, lambda="ML")
BIC(lat, bs)

# Test for effect of length.
length <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep
            + age_at_maturity + Saltwater, 
              data=c_data, lambda="ML")
BIC(bs, length)
# Confounding? Seems to be.
summary(bs)$coefficients
summary(length)$coefficients

# Test for effect of saltwater.
sw <- pgls(branchLength ~ numberOfNodes + DepthRangeDeep
               + age_at_maturity, 
               data=c_data, lambda="ML")
BIC(length, sw)
# Confounding? Seems to be.
summary(length)$coefficients
summary(sw)$coefficients


# Test for effect of depth.
depth <- pgls(branchLength ~ numberOfNodes
           + age_at_maturity, 
           data=c_data, lambda="ML")
BIC(sw, depth)
# Confounding? Seems to be.
summary(sw)$coefficients
summary(depth)$coefficients



null <- pgls(branchLength ~ numberOfNodes, data=c_data, lambda="ML")
anova(null, full, test="F")  # As a whole predictors are significant.
summary(full)$coefficients



tiph <- diag(vcv.phylo(firstTree))
# Run the PGLS.
fit <- gls(branchLength ~ numberOfNodes + BodyShapeI + Length + DepthRangeDeep + median_lat, 
           correlation=corPagel(value = 0, phy = firstTree), weights=varFixed(~tiph), 
           method = "ML", data = dfMultivariateCut)


# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + fertilization + FeedingType, 
           correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)



anova.pgls.fixed(fit2)


# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Length + Fresh + Saltwater + EnvTemp 
           + BodyShapeI + DepthRangeDeep + repro_mode + FeedingType +
             rep_guild_1 + pos_of_mouth + parental_care, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)




# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization + FeedingType, 
           correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization + FeedingType + rep_guild_1, 
           correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)




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

# TRAIT: Temp min.
dfTempMin <- merge(dfTempMin, dfThirdRegression, by = "species_name")
dfTempMin <- dfTempMin[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfTempMin, "species_name")
caperTempMin <- pgls(branchLength ~ TempMin.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Temp max.
dfTempMax <- merge(dfTempMax, dfThirdRegression, by = "species_name")
dfTempMax <- dfTempMax[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfTempMax, "species_name")
caperTempMax <- pgls(branchLength ~ TempMax.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Minimum depth
dfDepthShallow <- merge(dfDepthShallow, dfThirdRegression, by = "species_name")
dfDepthShallow <- dfDepthShallow[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDepthShallow, "species_name")
caperShallow <- pgls(branchLength.y ~ DepthRangeShallow.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: Maximum depth
# Merge so can get branch length info for third codon.
dfDepthDeep <- merge(dfDepthDeep, dfThirdRegression, by = "species_name")
dfDepthDeep <- dfDepthDeep[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDepthDeep, "species_name")
caperDeep <- pgls(branchLength.y ~ DepthRangeDeep.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: Diet Troph
# Merge so can get branch length info for third codon.
dfDietTroph <- merge(dfDietTroph, dfThirdRegression, by = "species_name")
dfDietTroph <- dfDietTroph[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfDietTroph, "species_name")
caperDiet <- pgls(branchLength.y ~ DietTroph.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: Age at maturity
# Merge so can get branch length info for third codon.
dfAgeMaturity <- merge(dfAgeMaturity, dfThirdRegression, by = "species_name")
dfAgeMaturity <- dfAgeMaturity[, c(1:2, 5, 7)]
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


# TRAIT: CORAL REEFS
# Merge so can get branch length info for third codon.
dfCoralReefs <- merge(dfCoralReefs, dfThirdRegression, by = "species_name")
dfCoralReefs <- dfCoralReefs[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfCoralReefs, "species_name")
caperCoral <- pgls(branchLength.y ~ CoralReefs.x + numberOfNodes.y, c_data, lambda = "ML")



# TRAIT: EPIPELAGIC
# Merge so can get branch length info for third codon.
dfEpipelagic <- merge(dfEpipelagic, dfThirdRegression, by = "species_name")
dfEpipelagic <- dfEpipelagic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfEpipelagic, "species_name")
caperEpi <- pgls(branchLength.y ~ Epipelagic.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: MESOPELAGIC
# Merge so can get branch length info for third codon.
dfMesopelagic <- merge(dfMesopelagic, dfThirdRegression, by = "species_name")
dfMesopelagic <- dfMesopelagic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfMesopelagic, "species_name")
caperMeso <- pgls(branchLength.y ~ Mesopelagic.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: BATHYPELAGIC
# Merge so can get branch length info for third codon.
dfBathypelagic <- merge(dfBathypelagic, dfThirdRegression, by = "species_name")
dfBathypelagic <- dfBathypelagic[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfBathypelagic, "species_name")
caperBathy <- pgls(branchLength.y ~ Bathypelagic.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: MANGROVES
# Merge so can get branch length info for third codon.
dfMangroves <- merge(dfMangroves, dfThirdRegression, by = "species_name")
dfMangroves <- dfMangroves[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfMangroves, "species_name")
caperMangroves <- pgls(branchLength.y ~ Mangroves.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: ESTUARIES
# Merge so can get branch length info for third codon.
dfEstuaries <- merge(dfEstuaries, dfThirdRegression, by = "species_name")
dfEstuaries <- dfEstuaries[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfEstuaries, "species_name")
caperEstuaries <- pgls(branchLength.y ~ Estuaries.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: SCHOOLING
# Merge so can get branch length info for third codon.
dfSchooling <- merge(dfSchooling, dfThirdRegression, by = "species_name")
dfSchooling <- dfSchooling[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfSchooling, "species_name")
caperSchooling <- pgls(branchLength.y ~ Schooling.x + numberOfNodes.y, c_data, lambda = "ML")


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


# TRAIT: Freshwater
# Merge so can get branch length info for third codon.
dfFreshwater <- merge(dfFreshwater, dfThirdRegression, by = "species_name")
dfFreshwater <- dfFreshwater[, c(1:2, 5, 7)]
dfFreshwater <- dfFreshwater[sample(nrow(dfFreshwater), 4393), ]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfFreshwater, "species_name")
caperFresh <- pgls(branchLength ~ Fresh.x + numberOfNodes, c_data, lambda = "ML")


# TRAIT: Saltwater
# Merge so can get branch length info for third codon.
dfSaltwater <- merge(dfSaltwater, dfThirdRegression, by = "species_name")
dfSaltwater <- dfSaltwater[, c(1:4)]
dfSaltwater <- dfSaltwater[sample(nrow(dfSaltwater), 4393), ]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfSaltwater, "species_name")
caperSaltwater <- pgls(branchLength.y ~ Saltwater.x + numberOfNodes.y, c_data, lambda = "ML")


# TRAIT: Env temp
# Merge so can get branch length info for third codon.
dfEnvTemp <- merge(dfEnvTemp, dfThirdRegression, by = "species_name")
dfEnvTemp <- dfEnvTemp[, c(1:4)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfEnvTemp, "species_name")
caperEnvTemp <- pgls(branchLength.y ~ EnvTemp.x + numberOfNodes.y, c_data, lambda = "ML")
anovaET <- anova.pgls.fixed(caperEnvTemp)


# FeedingType
# Merge so can get branch length info for third codon.
dfFeedingType <- merge(dfFeedingType, dfThirdRegression, by = "species_name")
dfFeedingType <- dfFeedingType[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfFeedingType, "species_name")
caperFeedingType <- pgls(branchLength.y ~ FeedingType.x + numberOfNodes.y, c_data, lambda = "ML")
anovaFT <- anova.pgls.fixed(caperFeedingType)


# Repro mode
# Merge so can get branch length info for third codon.
dfReproMode <- merge(dfReproMode, dfThirdRegression, by = "species_name")
dfReproMode <- dfReproMode[, c(1:3, 5)]
dfReproMode <- as.data.frame(dfReproMode)
dfReproMode <- dfReproMode[sample(nrow(dfReproMode), 1692), ]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfReproMode, "species_name")
caperReproMode <- pgls(branchLength ~ repro_mode.x + numberOfNodes, c_data, lambda = "ML")
anovaRM <- anova.pgls.fixed(caperReproMode)


# Fertilization
# Merge so can get branch length info for third codon.
dfFert <- merge(dfFert, dfThirdRegression, by = "species_name")
dfFert <- dfFert[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfFert, "species_name")
caperFert <- pgls(branchLength ~ fertilization.x + numberOfNodes, c_data, lambda = "ML")
anovaFert <- anova.pgls.fixed(caperFert)



# Reproductive guild I
# Merge so can get branch length info for third codon.
dfRepGuild1 <- merge(dfRepGuild1, dfThirdRegression, by = "species_name")
dfRepGuild1 <- dfRepGuild1[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfRepGuild1, "species_name")
caperRepGuild1 <- pgls(branchLength.y ~ rep_guild_1.x + numberOfNodes.y, c_data, lambda = "ML")
anovaRG1 <- anova.pgls.fixed(caperRepGuild1)


# Reproductive guild II
dfRepGuild2 <- merge(dfRepGuild2, dfThirdRegression, by = "species_name")
dfRepGuild2 <- dfRepGuild2[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfRepGuild2, "species_name")
caperRepGuild2 <- pgls(branchLength.y ~ rep_guild_2.x + numberOfNodes.y, c_data, lambda = "ML")
anovaRG2 <- anova.pgls.fixed(caperRepGuild2)


# Parental care
dfParentalCare <- merge(dfParentalCare, dfThirdRegression, by = "species_name")
dfParentalCare <- dfParentalCare[, c(1:2, 5, 7)]
dfParentalCare <- as.data.frame(dfParentalCare)
# Relevel body shape as I want fusiform / normal to be the reference.
dfParentalCare$parental_care.x <- relevel(dfParentalCare$parental_care.x, ref = "none")
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfParentalCare, "species_name")
caperParentalCare <- pgls(branchLength.y ~ parental_care.x + numberOfNodes.y, c_data, lambda = "ML")
anovaPC <- anova.pgls.fixed(caperParentalCare)


# Adult type swimming
dfSwimmingType <- merge(dfSwimmingType, dfThirdRegression, by = "species_name")
dfSwimmingType <- dfSwimmingType[, c(1:2, 5, 7)]
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfSwimmingType, "species_name")
caperSwimmingType <- pgls(branchLength.y ~ adult_type.x + numberOfNodes.y, c_data, lambda = "ML")
anovaST <- anova.pgls.fixed(caperSwimmingType)


# Position of mouth
dfPosMouth <- merge(dfPosMouth, dfThirdRegression, by = "species_name")
dfPosMouth <- dfPosMouth[, c(1:2, 5, 7)]
# Relevel body shape as I want fusiform / normal to be the reference.
dfPosMouth$pos_of_mouth <- relevel(dfPosMouth$pos_of_mouth, ref = "terminal")
# Test using caper.
c_data <- comparative.data(rootedWholeTree, dfPosMouth, "species_name")
caperPosMouth <- pgls(branchLength.y ~ pos_of_mouth + numberOfNodes.y, c_data, lambda = "ML")
anovaPM <- anova.pgls.fixed(caperPosMouth)




# Plot significant traits.
plot(log(DepthRangeDeep)~ log(age_at_maturity), data = dfMultivariateCut)
abline(model.pgls)
profile_lambda=pgls.profile(model.pgls, which="lambda") # vary lambda
plot(profile_lambda)




# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfThird[, c(1:2, 4:5, 10, 14, 23, 32, 38, 46)]
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
dfMultivariateCut <- dfMultivariate[, 1:10] 
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

table(dfMultivariateCut$Saltwater)

table(dfMultivariateCut$Lakes)

table(dfMultivariateCut$CoralReefs)

hist(dfMultivariateCut$age_at_maturity)
range(dfMultivariateCut$age_at_maturity)

lm.res=lm(branchLength ~ Length + DepthRangeDeep
          + Saltwater + Lakes + CoralReefs +
            age_at_maturity, data = dfMultivariateCut)
library(car)
vif(mod=lm.res)  # Longevity removed.

# Relevel traits.
dfMultivariateCut$Saltwater <- relevel(dfMultivariateCut$Saltwater, ref = "0")
dfMultivariateCut$Lakes <- relevel(dfMultivariateCut$Lakes, ref = "0")
dfMultivariateCut$CoralReefs <- relevel(dfMultivariateCut$CoralReefs, ref = "0")

# Backward selection using BIC.
c_data <- comparative.data(rootedWholeTree, dfMultivariateCut, names.col = "species_name", vcv=TRUE)

# Full model.
full <- pgls(branchLength ~ numberOfNodes + Length + DepthRangeDeep + median_lat +
             + Saltwater + Lakes + CoralReefs + age_at_maturity, data=c_data, lambda = "ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x=fitted(full), y=full$phyres, pch=5)

# Test for effect of mesopelagic.
meso <- pgls(branchLength ~ numberOfNodes + Length, 
           data=c_data, lambda="ML")
BIC(full, meso)
summary(full)$coefficients
summary(meso)$coefficients
# No effect so remove this trait.



null <- pgls(branchLength ~ numberOfNodes, data=c_data, lambda="ML")
anova(null, full, test="F")  # As a whole predictors are significant.



# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization, correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization + FeedingType, 
           correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)
# Add a variable.
fit <- gls(branchLength ~ numberOfNodes + Fresh + Saltwater + BodyShapeI + Length + fertilization + FeedingType + rep_guild_1, 
           correlation=corPagel(value = 0, phy = firstTree),
           weights=varFixed(~tiph), method = "ML", data = dfMultivariateCut)


# Checking for collinearity.
lm.res <- lm(branchLength ~ numberOfNodes
             + salinity_tolerance + env_salinity_level 
             + median_lat + lat_range
             + operculum_present + body_shape_I + Neritic + Intertidal
             + Oceanic + Epipelagic + Estuaries + Mangroves
             + Stream + Lakes + Schooling + Benthic + SoftBottom +
               HardBottom + Macrophyte + SeaGrassBeds
             + CoralReefs + temp_surface + Herbivory2 + repro_mode
             + FoodTroph + fertilization + rep_guild_1
             + length_max, data = dfRegression)
vif(mod = lm.res)

# Removing variables with high ( > 5) VIF values. 
lm.res <- lm(branchLength ~ numberOfNodes
             + endemic + salinity_tolerance + env_salinity_level 
             + median_lat + operculum_present
             + body_shape_I + max_depth + Neritic + Intertidal
             + Oceanic + Epipelagic + Mesopelagic + Estuaries + Mangroves
             + Stream + Lakes + Schooling + Benthic + SoftBottom + Sand + Silt
             + Mud + HardBottom + Rocky + Rubble + Macrophyte + SeaGrassBeds
             + CoralReefs + temp_surface + repro_mode
             + FoodTroph + rep_guild_1, 
             data = dfRegression)
vif(mod = lm.res)

# FEATURE SELECTION!
# Feature selection using Boruta.
# Adapted from: https://www.analyticsvidhya.com/blog/2016/03/select-important-variables-boruta-package/.
boruta <- Boruta(branchLength ~ numberOfNodes
                 + salinity_tolerance + env_salinity_level 
                 + median_lat + lat_range
                 + operculum_present + body_shape_I + Neritic 
                 + Oceanic + Epipelagic + Estuaries + Mangroves
                 + Stream + Lakes + Benthic + SoftBottom +
                   HardBottom
                 + CoralReefs + temp_surface + Herbivory2 + repro_mode
                 + FoodTroph + fertilization + rep_guild_1
                 + length_max, data = dfRegression, 
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
# Winner Boruta:
branchLength ~ numberOfNodes + endemic + salinity_tolerance + env_salinity_level 
+ median_lat + operculum_present + body_shape_I + max_depth + Neritic + 
  Oceanic + Estuaries + Mangroves + Stream + Lakes + Benthic + Macrophyte + SeaGrassBeds
+ CoralReefs + temp_surface + repro_mode
+ FoodTroph + rep_guild_1



## Using stepwise selection.
# Order the data according to the tree.
dfRegression <- dfRegression[match(overallMLtree, dfRegression$species_name), ]
fit <- gls(branchLength ~ numberOfNodes
           + salinity_tolerance + env_salinity_level 
           + median_lat + lat_range
           + operculum_present + body_shape_I + Neritic 
           + Oceanic + Epipelagic + Estuaries + Mangroves
           + Stream + Lakes + Benthic + SoftBottom + HardBottom
           + CoralReefs + temp_surface + Herbivory2 + repro_mode
           + FoodTroph + fertilization + rep_guild_1
           + length_max, data = dfRegression, 
             correlation = corBrownian(phy = treeGTR.G.I$tree), method = "ML")


# PGLS using features that were selected.
overallMLtree <- drop.tip(phy = overallMLtree, 
                                     tip = overallMLtree$tip.label[!overallMLtree$tip.label%in%as.character(dfRegression$species_name)], rooted = T)
tiph <- diag( vcv.phylo(overallMLtree) )
dfRegression <- dfRegression[match(overallMLtree$tip.label, dfRegression$species_name.x),]
fit <- gls(branchLength ~ numberOfNodes
           + salinity_tolerance + env_salinity_level 
           + median_lat + body_shape_I + Intertidal
           + Oceanic + Epipelagic + Estuaries + Mangroves
           + Stream + Lakes + Schooling + Macrophyte 
           + CoralReefs + temp_surface + Herbivory2 + repro_mode + 
             fertilization + rep_guild_1
           + length_max, correlation=corBrownian(phy=overallMLtree),
           weights=varFixed(~tiph), method = "ML", data=dfRegression)


           
full <- pgls(branchLength ~ numberOfNodes
             + salinity_tolerance + env_salinity_level 
             + median_lat + lat_range
             + operculum_present + body_shape_I + Neritic 
             + Oceanic + Epipelagic + Estuaries + Mangroves
             + Stream + Lakes + Benthic + SoftBottom + HardBottom
             + CoralReefs + temp_surface + Herbivory2 + repro_mode
             + FoodTroph + fertilization + rep_guild_1
             + length_max, data = cDat, 
             lambda = "ML")
summary(full)
# Inspect phylogenetic residuals for normality:
# Inspect a histogram:
par(mar = rep(2, 4))
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)

# Inspect phylogenetic residuals for homogeneity.
plot(x=fitted(full), y=full$phyres, pch=19)

# Some more plots.
#source("pgls_check_fncs.r")
#diagnostics.plot(full)

# Removal of outliers with studentized residuals >?3.
# Adapted from: http://www.anthrotree.info/wiki/pages/w5j9m8/7.5.html.
res <- residuals(full, phylo = TRUE) # Extracts phylogenetic residuals from the pgls model.
res <- res/sqrt(var(res))[1] # Standardizes residuals by sqrt of their variance.
rownames(res) <- rownames(full$residuals) # Matches the residuals up with the row name.
rownames(res)[(abs(res) > 3)] # Gives the names of the outliers.

# Remove if any outliers.
cDat <- cDat[-which(abs(res) > 3), ]
dfRegression <- dfRegression[-which(abs(res) > 3), ]
WholeConstrainedMLtreeMV <- drop.tip(phy = WholeConstrainedMLtreeMV, 
                   tip = WholeConstrainedMLtreeMV$tip.label[!WholeConstrainedMLtreeMV$tip.label%in%as.character(dfRegression$species_name)], rooted = T)
str(WholeConstrainedMLtreeMV)
full2 <- pgls(branchLength ~ numberOfNodes + endemic + salinity_tolerance + env_salinity_level 
              + median_lat + operculum_present + body_shape_I + max_depth + Neritic + 
                Oceanic + Estuaries + Mangroves + Stream + Lakes + Benthic + Macrophyte + SeaGrassBeds
              + CoralReefs + temp_surface + repro_mode
              + FoodTroph + rep_guild_1, data = cDat, 
              lambda = "ML")

# Look at the plots again.
par(mar = rep(2, 4))
hist(full2$phyres)
qqnorm(full2$phyres)
qqline(full2$phyres)
# Inspect phylogenetic residuals for homogeneity.
plot(x = fitted(full2), y = full2$phyres, pch = 19)
#diagnostics.plot(full2)

summary(full2)


# Endemic.
full2 <- pgls(branchLength ~ filtered_bin_size + number_of_nodes + rep_guild_1, data = cDat, 
              lambda = "ML")


anova.pgls.fixed(full2)

# Again using only univariate sig variables.
full2 <- pgls(branchLength ~ filtered_bin_size + number_of_nodes + 
                median_lat + operculum_present + body_shape_I + 
                max_depth + Neritic + Intertidal + Oceanic + Estuaries + Mangroves + 
                Stream + Lakes + Benthic + Macrophyte + CoralReefs + 
                temp_surface + repro_mode + FoodTroph + rep_guild_1, data = cDat, 
              lambda = "ML")




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
dfRegression <- dfRegression$species_name[WholeConstrainedMLtreeMV$tip.label, ]

# PGLS model assuming Brownian Motion.
fitBM <- gls(branchLength ~ median_lat + numberOfNodes, 
             corBrownian(phy=WholeConstrainedMLtreeMV), data = dfRegression)

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
                correlation = corPagel(value = 0.8, phy = WholeConstrainedMLtreeMV), 
                data = dfRegression)
intervals(fitPagel, which = "var-cov")

# Testing hypotheses of lambda = 0 vs. lambda = 0 using a likelihood ratio test.
# Independence - lambda = 0.
fitPagel0 <- gls(branchLength ~ median_lat + numberOfNodes, 
                 correlation = corPagel(value = 0, phy = WholeConstrainedMLtreeMV, fixed = TRUE), 
                 data = dfRegression)
# Brownian motion - lamdba = 1.
fitPagel1 <- gls(branchLength ~ median_lat + numberOfNodes, 
                 correlation = corPagel(value = 1, phy = WholeConstrainedMLtreeMV, fixed = TRUE), 
                 data = dfRegression) 
anova(fitPagel, fitPagel0)
anova(fitPagel, fitPagel1)

# Estimating lambda.
lambda <- seq(0, 1, length.out = 500)
lik <- sapply(lambda, function(lambda) 
  logLik(gls(branchLength ~ median_lat + numberOfNodes, 
             correlation = corPagel(value = lambda, phy = WholeConstrainedMLtreeMV, fixed = TRUE),
             data = dfRegression)))
plot(lik ~ lambda, type = "l", main = expression(paste("Likelihood Plot for ",
lambda)), ylab = "Log Likelihood", xlab = expression(lambda))
abline(v = fitPagel$modelStruct, col = "red")

# Plot branch length vs. median lat.
plot(branchLength ~ log(Length.x), dfMaxLength, pch = 19, 
     xlab = "Log(Maximum length (cm))", ylab = "Branch length")


lm(dfLongWild$branchLength ~ log(dfLongWild$LongevityWild.x))
abline(a = 1.108e+00, b = -0.0000018, col = "red")

with(dfLakes, plot(branchLength ~ Lakes.x, xlab = "Presence/absence in lakes",
                        ylab = "Branch length", 
                        main = " "))
axis(2, cex.axis = 1)


ggplot(dfParentalCare, aes(y=branchLength)) + geom_point(size=2, shape=19) + 
  geom_abline(intercept = 1.18785, slope = -0.03137, size = 0.75, color = "red") + 
  labs(x="Log(Longevity (years))", y="Branch length") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

ggplot(dfRegression, aes(branchLength)) + 
  geom_histogram(col="black",fill="deepskyblue1", alpha = .2) + 
  labs(x="Branch length", y="Frequency") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

qplot(dfRegression$branchLength,
      geom="histogram", 
      binwidth = 0.01)

ggplot(dfBodyShapeI, aes(x = BodyShapeI.x, y = branchLength)) +
  geom_boxplot() + labs(x="Body shape", y="Branch length") + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))


png(filename="bs.png", 
    units="in", 
    width=7, 
    height=5, 
    pointsize=12, 
    res=300)
my_sc_plot(data)
dev.off()



# Estimate lambda for each trait.
# Set up matrix to hold the results.
mat <- matrix(NA, ncol = 3, nrow = (dim(dat3)[2]), 
              dimnames = list(variable = names(dfRegression[blahtoblah]), 
                              c("Lower" , "Estimate", "Upper")))

# Loop through each variable.
for (i in 1:(dim(dat)[dfRegression])) {  
  form <- formula(paste(names(dat)[i], " ~ 1")) 
  this.fit <- gls(form, data = dfRegression, 
                  correlation = corPagel(0.1, phy = WholeConstrainedMLtreeMV),
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
cDat <- comparative.data(data = dfRegressionOpt, phy = WholeConstrainedMLtreeMV, 
                         names.col = "species_name", vcv = TRUE)

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
                          #names.col = "species_name", vcv = TRUE, vcv.dim = 3)
#modK <- pgls(branchLength ~ latitude, cDat2, kappa = "ML")
#summary(modK)
# pgls.profile() = likelihood profiles for branch length transformations.
# Kappa profile.
#K <- pgls.profile(modK, which = c("kappa"), N = 50, param.CI = NULL) 
#pgls.confint(modK, which = c('kappa'))$opt  # 0.4876088 

# Using delta.
#cDat3 <- comparative.data(data = dfRegression, phy = treeNJ, 
                          #names.col = "species_name", vcv = TRUE)
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
