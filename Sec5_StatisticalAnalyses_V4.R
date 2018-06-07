# # Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Adapted lines from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/1312519/by0.
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
# Author: https://stackoverflow.com/users/580110/henk.
# https://stackoverflow.com/questions/24016612/subset-of-data-frame-columns-to-maximize-complete-observations.
# Used as general guide for PGLS analyses (lines 172-255 and 319-358):
# Mundry, R. (2014). Statistical Issues and Assumptions of Phylogenetic Generalized Least Squares. In L.Z. Garamszegi (Ed.), 
# Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Online Practice Materials.
# URL: http://www.mpcm-evolution.org/practice/online-practical-material-chapter-6/chapter-6-1exercises-testing-assumptions-statistical-issues-framework-phylogenetic-generalized-least-squares

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

#############################################################################################################################

##### SECTION 5: STATISTICAL ANALYSES #####
# This section is designed to perform single variable and multivariable analyses regression analyses while controlling for phylogeny. The main objective is to
# identify those variables that contribute most to variation in branch length (molecular evolution rate) in fish. Tests to detect the amount of 
# phylogenetic signal in a trait are also performed.

### PACKAGES REQUIRED ###
# For data manipulation:
#install.packages("data.table")
library(data.table)
# For phylogenetic tree manipulation and analysis:
#install.packages("adephylo")
library(adephylo)
#install.packages("ape")
library(ape)
#install.packages("caper")
library(caper)
#install.packages("phytools")
library(phytools)
# For statistical analysis/graphs:
#install.packages("car")
library(car)
#install.packages("plotly")
library(plotly)
#install.packages("Rmisc")
library(Rmisc)
# Load the function(s) designed for this script:
source("GetTraitInfo.R")
source("TestPhyloSig.R")
source("PGLS.R")
source("MergeAndPGLS.R")

#############################################################################################################################

# A phylogenetic tree containing branch length data for your species is required for this section.
# Read in your phylogenetic tree.
mainTree <- read.tree(file = "final_wholeFishTree.tree")
# Fixing the tip labels (if the species names contain underscores when they are read in).
mainTree$tip.label <- gsub("_", " ", mainTree$tip.label)
# Root the tree using your chosen outgroup species.
mainTree <- root(mainTree, outgroup = outgroups, resolve.root = T)
# Match mainTree with data subset. This will ensure the tree has only the tips we need for data analysis.
mainTree <- drop.tip(phy = mainTree, tip = mainTree$tip.label[!mainTree$tip.label %in% dfCentroidSeqsWO$species_name])

### TRAIT: NUMBER OF NODES.
# Order dfTraits by the mainTree tips.
dfTraits <- dfTraits[match(mainTree$tip.label, dfTraits$species_name), ]
# Let's first determine the number of nodes for each species. This will be used as a control variable in the multivariable regression analysis (to account 
# for the node density effect).
dfTraits[, number_of_nodes := distRoot(mainTree, method = "nNodes")]

### TRAIT: BRANCH LENGTHS.
# Let's calculate the sum of branch lengths now (from root to tip). These values will serve as our measurement of molecular evolution rate.
dfTraits[, branch_length := distRoot(mainTree, method = "patristic")]
# Get info about the branch lengths.
GetTraitInfo(dfTraits$branch_length)
# Range within which 95% of the values fall.
quantile(dfTraits$branch_length, probs = c(.025, .975))

# Take a closer look at branch length outliers. Some contaminated sequences might have STILL gotten through, so it is best to check!
# Using the IQR to detect statistical outliers.
lowerQuantile <- quantile(dfTraits$branch_length)[2]
upperQuantile <- quantile(dfTraits$branch_length)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 3) + upperQuantile
lowerThreshold <-  lowerQuantile - (iqr * 3)
# Extreme short branches.
dfShort <- dfTraits[branch_length < lowerThreshold][, c(1, 14:15)]
# Get the sequence information in case you want to BLAST the sequence (also, we aren't interested in outgroup species here,
# that's why we are using dfCentroidSeqsNO).
dfShort <- merge(dfShort, dfCentroidSeqsNO, by = "species_name")
# Do the same for the extreme long branches.
dfLong <- dfTraits[branch_length > upperThreshold][, c(1, 14:15)]
dfLong <- merge(dfLong, dfCentroidSeqsNO, by = "species_name")
# Remove from dataset, if desired.
dfTraits <- RemoveSequences(dfTraits, c(dfShort$species_name, dfLong$species_name))

### SINGLE VARIABLE REGRESSION ANALYSIS ###
# Running a single variable PGLS regression analysis for each trait to determine whether significance can be detected. If so, they will be included 
# in the multivariable regression model selection process.

# First, make sure the trait data and phylo tree match (in case species were removed).
mainTree <- drop.tip(phy = mainTree, tip = mainTree$tip.label[!mainTree$tip.label %in% dfTraits$species_name])
dfTraits <- dfTraits[match(mainTree$tip.label, dfTraits$species_name), ]

### SINGLE-VARIABLE PGLS ANALYSES ###
# Use the PGLS function to perform single-variable (with number of nodes as a control variable) for all of the traits. 
# e.g. branch_length ~ trait_of_interest + number_of_nodes
# We will do this by looping through all of the columns containing the trait data using lapply.
traits <- as.list(colnames(dfTraits[, 4:13]))
# Set to dataframe.
dfTraits <- as.data.frame(dfTraits)
# Start the loop.
singleVarResults <- lapply(traits, function(x) {
  # We only want the columns containing species name and dependent and independent variables.
  data <- dfTraits[, c("species_name", x, "branch_length", "number_of_nodes")]
  # Remove NA values.
  data <- data[complete.cases(data), ]
  # Perform PGLS. The trait of interest in this case will always be the 2nd column.
  caper <- PGLS(data, mainTree, branch_length ~ data[, 2] + number_of_nodes)
  # Take the summary of the results.
  caperSum <- summary(caper)
})
# Assign names to the list of results based on the trait of interest.
names(singleVarResults) <- traits

# Which traits have p-values 0.15 or below?
# For now, this is only taking the first p-value of the trait (I still need to change it to deal with multi-level factors).
sigVars <- lapply(singleVarCapers, function(x) (x$coefficients[2,4]))
names(sigVars) <- names(singleVarResults)
# Which are below 0.15?
keepVars <- names(which(sigVars <= 0.15))
keepVars <- c("species_name", "branch_length", "number_of_nodes", keepVars)

### MULTIVARIABLE REGRESSION ANALYSES ###
# First, we must find the most complete dataset that we can use for our model selection process. So we will remove traits 
# with smaller sample sizes that don't overlap with data from other traits.
# First, get a dataframe of only those traits I am considering.
dfMultivariable <- as.data.table(dfTraits[, keepVars])
# We want the most complete dataset possible given our traits.
# First, order the columns by the amount of missing data (NA values).
dfTraitsNA <- sort(dfMultivariable[, lapply(.SD, function(x) sum(is.na(x)))])
# Reorder the original dfTraits. The columns with the least amount of NA values will now be first.
setcolorder(dfMultivariable, names(dfTraitsNA))
# Now I want to loop through the traits, removing one column (trait) at a time and count the number of complete cases. 
# This will provide us some information as to which traits would provide an adequate sample size for downstream analysis.
# First, take the number of columns in dfMultivariable.
len <- ncol(dfMultivariable)
# Create a numeric vector to hold the results of the loop.
all.cc <- NULL
# Start the loop:
for (i in 1:len) {
  # Works best if you set dfMultivariable back to a dataframe.
  x <- as.data.frame(dfMultivariable)
  # x is the placeholder dataframe in the loop.
  x <- x[, 1:len]
  # Determine which rows are "complete" using the "len" subset of traits.
  x <- complete.cases(x)
  # Complete rows of data will be "TRUE".
  x <- which(x == "TRUE")
  # Find the number of complete cases.
  x <- length(x)
  # Add it to the all.cc variable that's holding all of the results of the loop.
  all.cc[i] <- x
  # Minus 1 from tempLen so we can check the next subset of traits (we started at the last column because the columns were 
  # ordered by number of NA values).
  len <- len - 1
}
# Now, decide where to cut the datatable. (i.e. pick an adequate subset of 
# traits that maximize sample size).
# First, name it according to the trait columns.
names(all.cc) <- rev(colnames(dfMultivariable))
# Look at the results.
all.cc
# What seems like a good cut off point?
len <- which(colnames(dfMultivariable) == "max_length")
dfMultivariableCut <- dfMultivariable[, 1:len] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariableCut <- dfMultivariableCut[complete.cases(dfMultivariableCut)]

# Now check for data variability in this subset so our tests will actually work.
# Examples:
GetTraitInfo(dfMultivariableCut$median_lat)
GetTraitInfo(dfMultivariableCut$salinity)
GetTraitInfo(dfMultivariableCut$max_length)

# Check for multicollinearity between variables using the variance inflation factor (vif), if desired. 
# Multicollinearity can lead to errors in the estimations of our coefficients.
fit <- lm(branch_length ~ max_length + median_lat + salinity, data = dfMultivariableCut)
vif(mod = fit)

# MODEL SELECTION #
# In this section we will perform manual stepwise model selection. We will remove one trait at a time (usually the one with the highest p-value
# and examine the R squared and BIC values for each model. The model with the lowest BIC value and/or highest R squared will serve as
# our best-fit model.

# Now, let's perform a PGLS regression analysis using all of the variables.
# This is our "global" model.
global <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + median_lat + salinity + max_length)
# Check that the phylogenetic residuals are normal.
hist(global$phyres)
qqnorm(global$phyres)
qqline(global$phyres)
plot(x = fitted(global), y = global$phyres, pch = 5)

# Remove variables one at a time see if the fit of the model improves. 
fit1 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + median_lat + salinity)
# Compare the model to the global model using BIC.
BIC(global, fit1)

# Remove another variable.
fit2 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + median_lat + max_length)
# Compare the model to the global model using BIC.
BIC(global, fit1, fit2)

# Remove another variable.
fit3 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + salinity + max_length)
# Compare the model to the global model using BIC.
BIC(global, fit1, fit2, fit3)

# Continue process until you find the model with the lowest BIC value and/or highest R squared value.

# GRAPHS #
# For example, plot branch_length against median_lat.
fit <- lm(branch_length ~ median_lat, data = dfMultivariableCut)
dfMultivariableCut %>%  plot_ly(x = ~median_lat) %>% 
  add_markers(y = ~branch_length) %>% 
  add_lines(x = ~median_lat, y = fitted(fit))
