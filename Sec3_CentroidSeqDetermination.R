################################################################################

# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.
# Acknowlegements:  Centroid sequence selection (lines 38-111) designed by Matt Orton.

################################################################################

##### SECTION 3: CENTROID SEQUENCE DETERMINATION #####
# Centroid Sequence: BIN sequence with minimum average pairwise distance to all 
#                    other sequences in a given BIN.

### PACKAGES REQUIRED ###

# For data manipulation:
#install.packages("data.table")
library(data.table)

# For multiple sequence alignments:
#install.packages("ape")
library(ape)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Biostrings")
library(Biostrings)
#biocLite("muscle")
library(muscle)

# Load the function(s) that I designed for this script:
source("refSeqTrim.R")

################################################################################

# Subset dataframe to find BINs with more than one sequence.
largeBins <- which(dfPreCentroid$filtered_bin_size > 1)

# If there is at least one BIN with more than one member...
if (length(largeBins) > 0) {
  
  # Remove gaps from the sequences.
  dfPreCentroid[, nucleotides := gsub("-", "", nucleotides)] 
  
  # Subset out the BINs with more than 1 sequence.
  dfCentroidSeqs <- dfPreCentroid[largeBins, ]
  
  # How many unique BINs are in dfCentroidSeqs? 
  binNumberCentroid <- unique(dfCentroidSeqs$bin_uri)
  binNumberCentroid <- length(binNumberCentroid)
  
  # We also have to create another separate dataframe with BINs that only have 
  # one member, called dfSingletons.
  dfSingletons <- dfPreCentroid[-largeBins, ]
  
  # We then take the dfCentroidSeqs sequences and group them by BIN.
  largeBinList <- lapply(unique(dfCentroidSeqs$bin_uri), 
                         function(x) dfCentroidSeqs[dfCentroidSeqs$bin_uri == x,])
  
  # Extract the recordID from each BIN.
  largeBinRecordId <- lapply(largeBinList, function(x) (x$recordID))
  
  # Convert all the sequences in largeBinList to DNAStringSet format for 
  # multiple sequence alignment.
  DNAStringSet1 <- lapply(largeBinList, function(x) DNAStringSet(x$nucleotides))
  
  # Name DNAStringSet1 using the recordIDs.
  for (i in seq(from = 1, to = binNumberCentroid, by = 1)) {
    names(DNAStringSet1[[i]]) <- largeBinRecordId[[i]]
  }
  
  # Align the sequences in each BIN using MUSCLE.
  alignment1 <- lapply(DNAStringSet1, function(x) 
    muscle::muscle(x, diags = TRUE, gapopen = -3000))
  
  # Convert each BIN alignment to DNAbin format.
  dnaBINCentroid <- lapply(alignment1, function(x) as.DNAbin(x))
  
  # Estimates the genetic distance between sequences in each BIN with the TN93 model.
  geneticDistanceCentroid <- lapply(dnaBINCentroid, function(x) 
    dist.dna(x, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE))
  
  # Find the centroid sequence using the genetic distance matrix.
  # It is the sequence in a BIN with minimum average pairwise distance to all 
  # other sequences.
  centroidSeq <- lapply(geneticDistanceCentroid, function(x) 
    which.min(rowSums(x)))
  centroidSeq <- unlist(centroidSeq)
  centroidSeq <- names(centroidSeq)
  centroidSeq <- as.numeric(centroidSeq)
  
  # Subset dfCentroidSeqs by the recordIDs of the centroid sequences.
  dfCentroidSeqs <- subset(dfCentroidSeqs, recordID %in% centroidSeq)
  
  # Combine the singletons and centroid sequences into a new dataframe:
  # dfAllSeqs. Now each BIN has a representative sequence.
  dfAllSeqs <- rbind(dfCentroidSeqs, dfSingletons)
  
} else {
  
  # Centroid sequence selection not required if all BINs are singletons.
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

################################################################################
# REFERENCE SEQUENCE TRIMMING #
# Trim the centroid sequences according to a standardized reference sequence.
# Currently, a standard length (658 bp) COI-5P sequence from Perca flavescens 
# (yellow perch) is being used to trim Actinopterygii barcode sequences.

# Use the refSeqTrim function to trims nucleotide sequences in a dataframe 
# according to a given reference sequence.
dfCheckAllSeqs <- refSeqTrim(dfAllSeqs)
