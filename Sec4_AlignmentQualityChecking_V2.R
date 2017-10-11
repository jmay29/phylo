# # Copyright (C) 2017 Jacqueline May.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariable analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.
# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.

# Contributions & Acknowledgements #
# Adapted lines 102-106, 126-132, and 135-137 from code shared in Stack Overflow 
# discussion:
# Author: https://stackoverflow.com/users/1312519/by0.
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
# Author: https://stackoverflow.com/users/271678/fotnelton.
# https://stackoverflow.com/questions/20210787/r-getting-the-minimum-value-for-each-row-in-a-matrix-and-returning-the-row-and.
# Author: https://stackoverflow.com/users/235349/ramnath.
# https://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame.

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# There is a copy of the GNU General Public License
# along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

################################################################################

##### SECTION 4: ALIGNMENT QUALITY CHECKING #####
# This section performs alignment quality control checking by removing 
# extremely gappy sequences, outliers, and BINs that have neighbours in a 
# different taxonomic group (i.e. they may be contaminated or may have been
# misidentified).

### PACKAGES REQUIRED ###
# For data manipulation:
#install.packages("data.table")
library(data.table)
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

################################################################################

# Here, extremely gappy/ungappy sequences are removed. These sequences are 
# assumed to contribute to misalignment of the sequences or may even be 
# pseudogenes. Manual checking of the alignment is recommended.
# This will give the number of positions where an *internal* N or gap is found 
# for each sequence.
internalGaps <- sapply(regmatches(dfCheckCentroidSeqs$nucleotides, gregexpr("[-+]", 
                                  dfCheckCentroidSeqs$nucleotides)), length)
# Mean gap length and range.
meanGap <- mean(internalGaps)
extremeHighGap <- meanGap + 7  # Upper range.
extremeLowGap <- meanGap - 7  # Lower range.
# We then loop through each sequence to see if the number of gaps deviates 
# greatly from the mean.
# Which sequences exceed the range of meanGap +/- 7?
extremeSeqs <- sapply(internalGaps, function(x) which(x > extremeHighGap | x < extremeLowGap))
# The "deviant" sequences will be flagged with a 1.
extremeBins <- which(extremeSeqs > 0)
# Subset out these sequences to look at them if desired.
dfExtreme <- dfCheckCentroidSeqs[extremeBins, ]
# Make sure outgroups are not removed.
goodBins <- which(dfExtreme$family_name == "")
if (length(goodBins) > 1) {
  dfExtreme <- dfExtreme[-goodBins, ]
}
# Remove the gappy sequences from dfCentroidSeqs as we will be retrimming these 
# sequences again once troublesome cases are removed.
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, dfExtreme$bin_uri)

### OUTLIER CHECK ###
# Remove centroid sequences whose genetic distances to all other sequences fall 
# outside the typical range of genetic divergence for this group of organisms.
# Convert the sequences to DNAbin format so we can build a distance matrix.
dnaBinNN <- DNAStringSet(dfCheckCentroidSeqs$nucleotides)
names(dnaBinNN) <- dfCheckCentroidSeqs$bin_uri
dnaBinNN <- as.DNAbin(dnaBinNN)
# Then, we perform genetic distance determination with the TN93 model.
geneticDistanceCentroid <- dist.dna(dnaBinNN, model = "TN93", as.matrix = TRUE, 
                                    pairwise.deletion = TRUE)
# Using the IQR to detect outliers.
lowerQuantile <- quantile(geneticDistanceCentroid)[2]
upperQuantile <- quantile(geneticDistanceCentroid)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 1.5) + upperQuantile
# Remove 0 values (when a BIN is compared to itself - the diagonal values).
geneticDistanceCentroid[geneticDistanceCentroid == 0] <- NA
# Convert to dataframe format.
dfOutlierCheck <- as.data.frame(geneticDistanceCentroid)
# Identify BINs with no relatives within "typical" range of genetic divergence
# (i.e. all of their genetic distances fall outside 1.5 x IQR upper threshold.)
dfOutlierCheck <- dfOutlierCheck[apply(dfOutlierCheck, 1, function(x) all(x > upperThreshold)), ]
# Take a closer look.
outliers <- row.names(dfOutlierCheck)
# If desired, remove the outliers from dfCentroidSeqs as we will be retrimming 
# these sequences again once troublesome cases are removed.
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, outliers)

### NEAREST NEIGHBOUR CHECK ###
# Remove centroid sequences whose nearest neighbours are in a different order or 
# family. The nearest neighbour can be determined from the distance matrix.
# It is the sequence with minimum pairwise distance to the sequence in question.

# Identify the nearest neighbour of each BIN using geneticDistanceCentroid.
NearestNeighbour <- t(sapply(seq(nrow(geneticDistanceCentroid)), function(i) {
  # Identify the sequence with the minimum distance in the row (the nearest neighbour of i).
  j <- which.min(geneticDistanceCentroid[i, ])
  # Assign the name of the BIN and the name of its nearest neighbour to a separate column.
  c(paste(rownames(geneticDistanceCentroid)[i], colnames(geneticDistanceCentroid)[j], sep = '/'), 
    geneticDistanceCentroid[i, j])
}))
# Convert to dataframe.
dfNearestNeighbour <- as.data.frame(NearestNeighbour, stringsAsFactors = FALSE)
# Split the columns using the colsplit() function.
dfNearestNeighbour <- transform(dfNearestNeighbour, V1 = colsplit(V1, split = "/", 
                                names = c('bin_uri', 'nearest_neighbour')))
# Make sure the genetic distances are numeric.
dfNearestNeighbour$V2 <- as.numeric(dfNearestNeighbour$V2)
# Now, we need to get the order and family names of the BINs and their nearest neighbours.
# Convert into another dataframe first (we need to do this to fully separate
# the bin_uri and nearest_neighbour columns).
df2NearestNeighbour <- dfNearestNeighbour$V1
df2NearestNeighbour$pairwise_distance <- dfNearestNeighbour$V2
# Convert to datatable since it is faster.
df2NearestNeighbour <- setDT(df2NearestNeighbour)
# Make a datatable called dfBinOrd and match the order of the bin_uris in this 
# dataframe to the order of the bin_uris in dfCheckCentroidSeqs so that we can
# extract the order and family names of the BINs.
dfBinOrd <- dfCheckCentroidSeqs[match(df2NearestNeighbour$bin_uri, dfCheckCentroidSeqs$bin_uri), ]
# Do the same for the nearest neighbour bins.
dfNeighbourOrd <- dfCheckCentroidSeqs[match(df2NearestNeighbour$nearest_neighbour, dfCheckCentroidSeqs$bin_uri), ] 
# Add these columns to df2NearestNeighbour.
df2NearestNeighbour[, bin_order := dfBinOrd$order_name][, bin_family := dfBinOrd$family_name]
df2NearestNeighbour[, nn_order := dfNeighbourOrd$order_name][, nn_family := dfNeighbourOrd$family_name]
# Non-matching orders.
nonmatchOrd <- which(df2NearestNeighbour$bin_order != df2NearestNeighbour$nn_order)
# Look closely at the cases where the orders do not match.
dfNonmatchOrd <- df2NearestNeighbour[nonmatchOrd, ]
# Which are very close neighbours (i.e. there may be something weird going on)?
dfNonmatchOrd <- dfNonmatchOrd[pairwise_distance < 0.05]
# The next few lines may be used to remove the sequences, if any are found! 
# (you should double check these cases on the BOLD website, too).
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, dfNonmatchOrd$bin_uri)
# Now, let's take a look at the BINs where the nearest neighbour is in a 
# different family.
nonmatchFam <- which(df2NearestNeighbour$bin_family != df2NearestNeighbour$nn_family)
# Which are very close neighbours?
dfNonmatchFam <- df2NearestNeighbour[nonmatchFam, ]
dfNonmatchFam <- dfNonmatchFam[pairwise_distance < 0.05]
# Again, remove non-matching BINs if desired:
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, dfNonmatchFam$bin_uri)

# Now, you can perform a more comprehensive check if desired. Do the
# BINs have ANY close neighbours within a genetic distance of 0.05 that are in 
# a different order or family? Again, this may be indicative of something weird
# going on in either the sequence data or taxonomic assignment.
# First, subset out all close neighbour pairings that share a genetic distance
# under 0.05 to any other sequence.
dfGeneticDistance <- as.data.frame(geneticDistanceCentroid)
closeNeighbours <- which(dfGeneticDistance < 0.05, arr.ind = TRUE)
# Take the names of the BINs and nearest neighbours and convert to dataframe.
dfCloseNeighbours <- as.data.frame(row.names(dfGeneticDistance[closeNeighbours[, 1],] ))
colnames(dfCloseNeighbours)[1] <- "bins"
dfCloseNeighbours$neighbour <- colnames(dfGeneticDistance[closeNeighbours[, 2]])
# Following the same protocol as before, extract the order and family names of 
# the BINs and their close neighbours.
dfBinOrd <- dfCheckCentroidSeqs[match(dfCloseNeighbours$bins, dfCheckCentroidSeqs$bin_uri), ] 
dfNeighbourOrd <- dfCheckCentroidSeqs[match(dfCloseNeighbours$neighbour, dfCheckCentroidSeqs$bin_uri), ] 
# Assign to dfCloseNeighbours.
dfCloseNeighbours <- setDT(dfCloseNeighbours)
dfCloseNeighbours[, bin_order := dfBinOrd$order_name][, bin_family := dfBinOrd$family_name]
dfCloseNeighbours[, nn_order := dfNeighbourOrd$order_name][, nn_family := dfNeighbourOrd$family_name]
# Some BINs have multiple close neighbours and were assigned extra numbers in 
# the dataframe (e.g. ACFGGF, ACFGGF.1). To resolve this, use the function 
# na.locf() from the "zoo" package to fill in missing data based on the data in 
# the row above.
# This section seems to agree better with dataframes at the moment.
dfCloseNeighbours <- as.data.frame(dfCloseNeighbours)
# Make sure they are in the correct order in the "bins" column.
dfCloseNeighbours <- dfCloseNeighbours[order(dfCloseNeighbours$bins), ]
# Fill in the bin_order column.
dfCloseNeighbours$bin_order <- na.locf(dfCloseNeighbours$bin_order)
# Fill in the bin_family column.
dfCloseNeighbours$bin_family <- na.locf(dfCloseNeighbours$bin_family)
# Do the same for the neighbour BINs.
dfCloseNeighbours <- dfCloseNeighbours[order(dfCloseNeighbours$neighbour), ]
dfCloseNeighbours$nn_order <- na.locf(dfCloseNeighbours$nn_order)
dfCloseNeighbours$nn_family <- na.locf(dfCloseNeighbours$nn_family)
# Convert back to datatable.
dfCloseNeighbours <- setDT(dfCloseNeighbours)
# Trim the BIN URIs to remove the added numbers.
dfCloseNeighbours[, bins := substr(bins, 1, 7)]
dfCloseNeighbours[, neighbour := substr(neighbour, 1, 7)]
# Now, which orders do not match between bin_order and nn_order?
nonmatchOrd <- which(dfCloseNeighbours$bin_order != dfCloseNeighbours$nn_order)
dfNonmatchOrd <- dfCloseNeighbours[nonmatchOrd, ]
# Again, remove non-matching BINs if desired:
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, unique(dfNonmatchOrd$bin_uri))
# Which families do not match between bin_family and nn_family?
dfNonmatchFam <- dfCloseNeighbours[which(dfCloseNeighbours$bin_family != dfCloseNeighbours$nn_family)]
# To remove:
dfCentroidSeqs <- RemoveSequences(dfCentroidSeqs, unique(dfNonmatchFam$bin_uri))

### OUTGROUP CHECK ###
# Which outgroups made the cut? Remove them from dfCentroidSeqs so I can build 
# a tree just using the ingroup (so that inclusion of outgroups in the tree 
# building process doesn't affect the branch length estimates of the in-group).
outgroupSpecies <- unique(dfOutgroup$species_name)
dfGoodOutgroups <- dfCentroidSeqs[dfCentroidSeqs$species_name %in% outgroupSpecies]
# Remove the outgroups from dfCentroidSeqs.
dfCentroidSeqsNO <- dfCentroidSeqs[!dfCentroidSeqs$species_name %in% outgroupSpecies]
# Now, re-trim and align the sequences without the outgroups.
dfCentroidSeqsNO <- RefSeqTrim(dfCentroidSeqsNO)
# Once finished, make sure to check over sequences/alignment, and make sure 
# they are in the correct reading frame. Make sure to save the resulting 
# alignments under a different name, or save in a new directory so they are not 
# replaced.
# Now re-run the alignment including outgroups (pick outgroup species that are
# well represented and that serve as an appropriate outgroup to your taxa).
# Enter your chosen outgroup species between the "".
goodOG <- which(dfCentroidSeqs$species_name == "")
mergeOG <- dfCentroidSeqs[goodOG, ]
# Add them back.
dfCentroidSeqsWithOG <- rbind(dfCentroidSeqsNO, mergeOG)
# Run the alignment with outgroups included.
dfCentroidSeqsWithOG <- RefSeqTrim(dfCentroidSeqsWithOG)
