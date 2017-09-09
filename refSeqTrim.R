# Function: RefSeqTrim ##
# Purpose: Trims nucleotide sequences in a dataframe according to a given
#          reference sequence.
# Acknowlegements: This function adapted from code provided by Matt Orton (https://github.com/m-orton/R-Scripts).
RefSeqTrim <- function(data) {
  # data = Dataframe containing sequence information.
  
  # Insert your reference sequence here! Between the "".
  dfRefSeqs <- data.frame(taxa = c(""), nucleotides = c(""))
  colnames(dfRefSeqs)[2] <- "nucleotides"
  # Convert to datatable for some data manipulation.
  dfRefSeqs <- setDT(dfRefSeqs)
  dfRefSeqs[, nucleotides := as.character(nucleotides)]
  # Symmetrical trimming of the references to a standard 620 bp.
  # The user can alter this if desired.
  dfRefSeqs[, nucleotides := substr(nucleotides, 20, nchar(nucleotides) - 19)]
  # Check sequence length.
  dfRefSeqs[, seq_length := nchar(nucleotides)]
  # We must ensure that the sequences are of the chr type.
  alignmentSeqs <- as.character(data$nucleotides)
  # Name our sequences according to species names.
  names(alignmentSeqs) <- data$species_name
  alignmentRef <- as.character(dfRefSeqs$nucleotides[1])
  # Name our reference sequence "REFERENCE".
  names(alignmentRef) <- "REFERENCE"
  # Append our sequences together.
  alignmentSeqsPlusRef <- append(alignmentRef, alignmentSeqs)
  # Convert all sequences in alignmentSeqsPlusRef to DNAStringSet format. 
  # This is the format required for the alignment.
  DNAStringSet2 <- DNAStringSet(alignmentSeqsPlusRef)
  rm(alignmentSeqsPlusRef)
  # Run a multiple sequence alignment using MUSCLE.
  alignment2 <- muscle::muscle(DNAStringSet2, diags = TRUE, gapopen = -3000)
  # If you want to save the alignment as a FASTA file:
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
  # Find the start position of the reference sequence.
  refSeqPosStart <- regexpr("[ACTG]", refSeqPos)
  refSeqPosStart <- as.numeric(refSeqPosStart)
  # Find last nucleotide position of the reference sequence.
  refSeqPosEnd <- nchar(dfRefSeqs$nucleotides[1]) + refSeqPosStart
  refSeqPosEnd <- as.numeric(refSeqPosEnd)
  # Take a substring of the alignment by using these positions to "trim" the 
  # alignment.
  alignment2Trimmed <- substr(alignment2, refSeqPosStart, refSeqPosEnd)
  # Again, convert to DNAStringSet format.
  DNAStringSet3 <- DNAStringSet(alignment2Trimmed)
  # Checking alignment.
  classFileNames <- foreach(i = 1:nrow(dfRefSeqs)) %do% 
    paste("alignmentTrimmed", dfRefSeqs$taxa[i], ".fas", sep = "")
  writeXStringSet(DNAStringSet3, file = classFileNames[[1]], 
                  format = "fasta", width = 1500)
  # Remove the reference sequence.
  refSeqRm <- which(DNAStringSet3@ranges@NAMES == "REFERENCE")
  DNAStringSet3 <- subset(DNAStringSet3[-refSeqRm])
  # Reorder dfAllSeqs according to the order of the alignment.
  alignmentOrder <- DNAStringSet3@ranges@NAMES
  # Order dfAllSeqs according to this.
  data <- data[match(alignmentOrder, data$species_name), ]
  # Replace the old sequences with the new sequences.
  trimmedSeqs <- as.character(DNAStringSet3)
  data$nucleotides <- trimmedSeqs
  # Return a dataframe with the newly trimmed sequences.
  return(data)
}