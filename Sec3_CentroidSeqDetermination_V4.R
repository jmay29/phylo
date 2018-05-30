# # Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Centroid sequence selection designed by Matt Orton (https://github.com/m-orton/R-Scripts).

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

#############################################################################################################################

##### SECTION 3: CENTROID SEQUENCE DETERMINATION #####
# This section is designed to select a centroid sequence for each BIN. A centroid sequence is the sequence in a BIN with minimum sum of pairwise
# distance to all other sequences in the BIN. It will serve as a representative sequence for the BIN/species.

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
#install.packages("foreach")
library(foreach)
# Load the function(s) designed for this script:
source("RefSeqTrim.R")

#############################################################################################################################

# Subset dataframe to find BINs with more than one sequence.
largeBins <- dfPreCentroid[filtered_bin_size > 1]
# If there is at least one BIN with more than one sequence...
if (nrow(largeBins) > 0) {
  # Remove gaps from the sequences.
  dfPreCentroid[, nucleotides := gsub("-", "", nucleotides)] 
  # Subset out the BINs with more than 1 sequence.
  dfCentroidSeqs <- dfPreCentroid[bin_uri %in% largeBins$bin_uri]
  # We also have to create another separate dataframe with BINs that only have one sequence, called dfSingletons.
  dfSingletons <- dfPreCentroid[!bin_uri %in% largeBins$bin_uri]
  # We then take the dfCentroidSeqs sequences and group them by BIN.
  largeBinList <- split(dfCentroidSeqs, by = "bin_uri")
  # Convert all the sequences in largeBinList to DNAStringSet format for 
  # the multiple sequence alignment.
  DNAStringSetList <- lapply(largeBinList, function(x) DNAStringSet(x$nucleotides))
  # Name DNAStringSetList using the recordIDs.
  for (i in seq(from = 1, to = length(unique(dfCentroidSeqs$bin_uri)), by = 1)) {
    names(DNAStringSetList[[i]]) <- largeBinList[[i]]$recordID
  }
  # Align the sequences in each BIN using MUSCLE.
  alignmentList <- lapply(DNAStringSetList, function(x) muscle::muscle(x, diags = TRUE, gapopen = -3000))
  # Convert each BIN alignment to DNAbin format.
  alignmentList <- lapply(alignmentList, function(x) as.DNAbin(x))
  # Estimates the genetic distance between sequences in each BIN with the TN93 model.
  distanceMatrixList <- lapply(alignmentList, function(x) dist.dna(x, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE))
  # Find the centroid sequence using the genetic distance matrix. It is the sequence in a BIN with minimum average pairwise distance to all other sequences in the BIN.
  centroidSeqs <- sapply(distanceMatrixList, function(x) which.min(rowSums(x)))
  centroidSeqs <- names(centroidSeqs)
  centroidSeqs <- gsub("^.*\\.", "", centroidSeqs)
  centroidSeqs <- as.numeric(centroidSeqs)
  # Subset dfCentroidSeqs by the recordIDs of the centroid sequences.
  dfCentroidSeqs <- dfCentroidSeqs[dfCentroidSeqs$recordID %in% centroidSeqs]
  # Combine the singletons and centroid sequences into a new dataframe. Now each BIN has a representative sequence.
  dfCentroidSeqs <- rbind(dfCentroidSeqs, dfSingletons)
} else {
  # Centroid sequence selection not required if all BINs are singletons.
  dfCentroidSeqs <- dfPreCentroid
}

# REFERENCE SEQUENCE TRIMMING #
# Trim the centroid sequences according to a standardized reference sequence. Currently, a standard length (658 bp) COI-5P sequence from 
# Perca flavescens (yellow perch) is being used to trim Actinopterygii barcode sequences.

# Use the RefSeqTrim function to trim nucleotide sequences in a dataframe according to a given reference sequence.
dfCheckCentroidSeqs <- RefSeqTrim(dfCentroidSeqs)

# Remove objects that are not required for Section 4.
rm(alignmentList, centroidSeqs, first_time, i); rm(dfPreCentroid, dfLargeBins, dfSingletons); rm(largeBinList, distanceMatrixList, DNAStringSetList)

