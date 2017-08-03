################################################################################

# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements: Matt Orton for filtering steps/centroid sequence 
#                  determination/reference sequence trimming (lines TBD).

################################################################################

##### SECTION 4: ALIGNMENT QUALITY CHECKING #####

### PACKAGES REQUIRED ###

# For data manipulation:
#install.packages("data.table")
library(data.table)
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

# Load the function(s) that I designed for this script:
source("refSeqTrim.R")

################################################################################

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
# We then loop through each sequence to see if the number of gaps deviates 
# greatly from the mean.
# Which sequences exceed the range of meanGap +/- 7?
extremeSeqs <- foreach(i = 1:nrow(dfCheckAllSeqs)) %do% 
  which(internalGaps[[i]] > extremeHighGap | internalGaps[[i]] < extremeLowGap)
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

## FUNCTION: NearestNeighbour ##
# Purpose: Identifies the nearest neighbour of each BIN.
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
                  dfAllSeqs$species_name == "Gadus morhua" |
                  dfAllSeqs$species_name == "Larimichthys polyactis" |
                  dfAllSeqs$species_name == "Acanthurus triostegus")
mergeOG <- dfAllSeqs[goodOG, ]
# Add them back.
dfAllSeqsWithOG <- rbind(dfAllSeqsNO, mergeOG)
# Run the alignment with outgroups included.
dfAllSeqsWithOG <- refSeqTrim(dfAllSeqsWithOG)
