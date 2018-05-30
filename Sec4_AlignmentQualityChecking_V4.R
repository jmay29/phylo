# # Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Adapted lines from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/1312519/by0.
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
# Author: https://stackoverflow.com/users/2474755/j-r.
# https://stackoverflow.com/questions/27892100/distance-matrix-to-pairwise-distance-list-in-r.

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

#############################################################################################################################

##### SECTION 4: ALIGNMENT QUALITY CHECKING #####
# This section performs alignment quality control checking by removing extremely gappy sequences, outliers, and BINs that have neighbours in a 
# different taxonomic group (i.e. they may be contaminated or may have been misidentified).

### PACKAGES REQUIRED ###
# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages(stringr)
library(stringr)
#install.packages("foreach")
library(foreach)
#install.packages("reshape")
library(reshape)
# For multiple sequence alignments:
#install.packages("ape")
library(ape)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
#biocLite("muscle")
library(muscle)
# For missing data:
#install.packages("zoo")
library(zoo)
# Load the function(s) designed for this script:
source("RefSeqTrim.R")
source("RemoveSequences.R")

#############################################################################################################################

### QUALITY CHECK 1: GAPPY SEQUENCES ###
# Here, extremely gappy/ungappy sequences are removed. These sequences are assumed to contribute to misalignment of the 
# sequences or may even be pseudogenes. Manual checking of the alignment is recommended.
# Determine the number of positions where an *internal* N or gap is found for each sequence.
dfCheckCentroidSeqs[, internal_gapN := str_count(nucleotides, c("[-+]"))]
# Which sequences are NOT within the range of mean number of gaps in the centroid sequences +/- 7? These represent extremely gappy sequences.
dfGappySeqs <- dfCheckCentroidSeqs[!internal_gapN %between% c(mean(internal_gapN) - 7, mean(internal_gapN) + 7)]
# Make sure outgroups are not included in the dfGappySeqs!
dfGappySeqs <- dfGappySeqs[!species_name %in% outgroups]
# Remove the gappy sequences from the original dfCentroidSeqs as we will be realigning these sequences again once troublesome cases are removed.
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, dfGappySeqs$species_name)

### QUALITY CHECK 2: OUTLIERS ###
# Remove centroid sequences whose genetic distances to all other sequences fall outside the typical range of genetic divergence for this group of organisms.
# First, convert the sequences to DNAbin format so we can build a distance matrix.
DNABinNN <- DNAStringSet(dfCheckCentroidSeqs$nucleotides)
names(DNABinNN) <- dfCheckCentroidSeqs$species_name
DNABinNN <- as.DNAbin(DNABinNN)
# Then, we construct a pairwise distance matrix using the TN93 model.
distanceMatrix <- dist.dna(DNABinNN, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)
# Use the upper threshold of the IQR to detect outliers.
lowerQuantile <- quantile(distanceMatrix)[2]
upperQuantile <- quantile(distanceMatrix)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 1.5) + upperQuantile
# Remove 0 values so that these are not considered (when a species is compared to itself - the diagonal values).
distanceMatrix[distanceMatrix == 0] <- NA
# Convert to datatable.
dfOutliers <- as.data.table(distanceMatrix, keep.rownames = T)
# Rename the "rn" column (row names).
setnames(dfOutliers, "rn", "species_name")
# Identify BINs with no relatives within "typical" range of genetic divergence (i.e. all of their genetic distances are greater than 1.5 x IQR upper threshold.)
dfOutliers <- dfOutliers[, outlier := apply(.SD, 1, function(x) all(x > upperThreshold, na.rm = T))][outlier == TRUE]
# Make sure outgroups are not included in the outliers!
dfOutliers <- dfOutliers[!species_name %in% outgroups]
# If desired, remove the outliers from dfCentroidSeqs.
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, dfOutliers$species_name)

### QUALITY CHECK 3: CLOSE NEIGHBOUR TAXONOMY ###
# Remove centroid sequences whose close neighbours are in a different order or family. Close neighbours can be determined 
# from the distance matrix. They are sequences that are within a genetic distance of 0.05. If these neighbours are in a different
# order or family, this may be indicative of something weird going on in either the sequence data or taxonomic assignment. 
dfGeneticDistance <- as.data.table(distanceMatrix)
# Convert the distance matrix to a datatable with the names of the species pairs and their distances.
dfGeneticDistance <- data.table(t(combn(names(dfGeneticDistance), 2)), distance = t(dfGeneticDistance)[lower.tri(dfGeneticDistance)])
setnames(dfGeneticDistance, old = c("V1", "V2"), new = c("species_1", "species_2"))
# Subset out all close neighbour pairings that share a genetic distance under 0.05 to any other sequence.
dfGeneticDistance <- dfGeneticDistance[distance < 0.05]
# Get the order and families names of the species from dfCentroidSeqs.
dfGeneticDistance <- merge(dfGeneticDistance, dfCentroidSeqs[, c(1, 5:6)], by.x = "species_1", by.y = "species_name")
dfGeneticDistance <- merge(dfGeneticDistance, dfCentroidSeqs[, c(1, 5:6)], by.x = "species_2", by.y = "species_name")
setnames(dfGeneticDistance, old = c("order_name.x", "family_name.x", "order_name.y", "family_name.y"), 
         new = c("order_1", "family_1", "order_2", "family_2"))
# Now, which orders do not match between order_1 and order_2?
dfMismatchOrders <- dfGeneticDistance[order_1 != order_2]
# Remove these species from dfCentroidSeqs if desired:
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, c(unique(dfMismatchOrders$species_1), unique(dfMismatchOrders$species_2)))
# Now, which families do not match between family_1 and family_2?
dfMismatchFamilies <- dfGeneticDistance[family_1 != family_2]
# Remove these species from dfCentroidSeqs if desired:
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, c(unique(dfMismatchFamilies$species_1), unique(dfMismatchFamilies$species_2)))

### OUTGROUP CHECK ###
# Which outgroups made the cut? Remove them from dfCentroidSeqs so I can build a tree just using the ingroup 
# (so that inclusion of outgroups in the tree building process doesn't affect the branch length estimates of the in-group).
outgroupSpecies <- unique(dfOutgroup$species_name)
dfGoodOutgroups <- dfCentroidSeqs[dfCentroidSeqs$species_name %in% outgroupSpecies]
# Remove the outgroups from dfCentroidSeqs.
dfCentroidSeqsNO <- dfCentroidSeqs[!dfCentroidSeqs$species_name %in% outgroupSpecies]
# Now, re-trim and align the sequences without the outgroups.
dfCentroidSeqsNO <- RefSeqTrim(dfCentroidSeqsNO)
# Once finished, make sure to check over sequences/alignment, and make sure they are in the correct reading frame. 
# Make sure to save the resulting alignments under a different name, or save in a new directory so they are not replaced.
# Now re-run the alignment including outgroups (pick outgroup species that are well represented and that serve as an appropriate 
# outgroup to your taxa).
# Enter your chosen outgroup species between the "".
goodOG <- which(dfCentroidSeqs$species_name == "Neoceratodus forsteri")
mergeOG <- dfCentroidSeqs[goodOG, ]
# Add them back.
dfCentroidSeqsWithOG <- rbind(dfCentroidSeqsNO, mergeOG)
# Run the alignment with outgroups included.
dfCentroidSeqsWithOG <- RefSeqTrim(dfCentroidSeqsWithOG)
