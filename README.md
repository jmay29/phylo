# phylo

Pipeline for detecting associations between traits and COI molecular evolutionary rates in ray-finned fishes.

Information about the different sections:

Section 1: This section is primarily for quality control purposes of DNA barcode data obtained from the BOLD API. Filter arguments may be altered to meet the user's needs.

Section 2: This section is designed to assign trait data from BOLD and FishBase and match it to BOLD sequence data. The version here is a shortened version and only considers a subset of traits. Many more traits are available on FishBase and can be easily accessed using the package rfishbase.

Section 3: This section is designed to select a centroid sequence for each BIN. A centroid sequence is the sequence in a BIN with minimum sum of pairwise distance to all other sequences in the BIN. It will serve as a representative sequence for the BIN/species.

Section 4: This section performs alignment quality control checking by removing extremely gappy sequences, outliers, and BINs that have neighbours in a different taxonomic group (i.e. they may be contaminated or may have been misidentified).

Section 5: This section is designed to perform single variable and multivariable analyses regression analyses while controlling for phylogeny. The main objective is to identify those variables that contribute most to variation in molecular evolution rate.

All code and required functions are available in the /Code folder. All associated data are available in the /Data folder.
