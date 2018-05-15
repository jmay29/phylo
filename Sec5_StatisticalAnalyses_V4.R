# # Copyright (C) 2018 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting fish molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.
# Adapted lines 112-117 and 268-272 from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/1312519/by0.
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
# Author: https://stackoverflow.com/users/580110/henk.
# https://stackoverflow.com/questions/24016612/subset-of-data-frame-columns-to-maximize-complete-observations.
# Used as general guide for PGLS analyses (lines 172-255 and 319-358):
# Mundry, R. (2014). Statistical Issues and Assumptions of Phylogenetic Generalized Least Squares. In L.Z. Garamszegi (Ed.), 
# Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Online Practice Materials.
# URL: http://www.mpcm-evolution.org/practice/online-practical-material-chapter-6/chapter-6-1exercises-testing-assumptions-statistical-issues-framework-phylogenetic-generalized-least-squares

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. 
# Or view it directly here at http://www.gnu.org/licenses/

################################################################################

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
#install.packages("ggplot2")
library(ggplot2)
#install.packages("gtools")
library(gtools)
#install.packages("Rmisc")
library(Rmisc)
# Load the function(s) designed for this script:
source("TestPhyloSig.R")
source("PGLS.R")
source("MergeAndPGLS.R")

################################################################################

# A phylogenetic tree containing branch length data for your species is required for this section.
# Read in your phylogenetic tree.
mainTree <- read.tree(file = "final_wholeFishTree.tree")
# Fixing the tip labels.
mainTree$tip.label <- gsub("_", " ", mainTree$tip.label)
# Root the tree using your chosen outgroup species.
outgroups <- c("Neoceratodus forsteri")
mainTree <- root(mainTree, outgroup = outgroups, resolve.root = TRUE)
# Match mainTree with data subset. This will ensure the tree has only the tips we need for data analysis.
mainTree <- drop.tip(phy = mainTree, tip = mainTree$tip.label[!mainTree$tip.label %in% dfCentroidSeqs$species_name])

### TRAIT: NUMBER OF NODES.
# Let's first determine the number of nodes for each species. This will be used as a control variable in the multivariable regression analysis (to account 
# for the node density effect).
number_of_nodes <- distRoot(mainTree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(number_of_nodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)
# Merge with dfCentroidSeqsWithOG.
dfSeqsNodes <- merge(dfCentroidSeqsWithOG, dfNumberOfNodes, by = "species_name")

### TRAIT: BRANCH LENGTHS.
# Let's calculate the sum of branch lengths now (from root to tip). These values will serve as our measurement of molecular evolution rate.
branch_length <- distRoot(mainTree, tips = "all", method = "patristic")
# Convert to dataframe format.
dfBranchLengths <- data.frame(branch_length)
dfBranchLengths$species_name <- names(branch_length)
# Check the distribution of the branch lengths.
hist(dfBranchLengths$branch_length, main = "", xlab = "Branch Length")
# Some more stats.
median(dfBranchLengths$branch_length)
mean(dfBranchLengths$branch_length)
range(dfBranchLengths$branch_length)
# Range within which 95% of the values fall.
quantile(dfBranchLengths$branch_length, probs = c(.05, .975))

# Take a closer look at branch length outliers. Some contaminated sequences might have STILL gotten through, so it is best to check!
# Using the IQR to detect statistical outliers.
lowerQuantile <- quantile(dfBranchLengths$branch_length)[2]
upperQuantile <- quantile(dfBranchLengths$branch_length)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 3) + upperQuantile
lowerThreshold <-  lowerQuantile - (iqr * 3)
# Extreme short branches.
exShortBranchLengths <- which(dfBranchLengths$branch_length < lowerThreshold)
# Take a closer look.
dfCheckShort <- dfBranchLengths[exShortBranchLengths, ]
dfCheckShort <- merge(dfCheckShort, dfCentroidSeqsNO, by = "species_name")
# Do the same for the extreme long branches.
exLongBranchLengths <- which(dfBranchLengths$branch_length > upperThreshold)
dfCheckLong <- dfBranchLengths[exLongBranchLengths, ]
dfCheckLong <- merge(dfCheckLong, dfCentroidSeqsNO, by = "species_name")
# BLAST these and/or remove from dataset, if desired.
# Short branch outliers.
if (nrow(dfCheckShort) > 0) {
  dfBranchLengths <- dfBranchLengths[!dfBranchLengths$species_name %in% dfCheckShort$species_name, ]
}
# Long branch outliers.
if (nrow(dfCheckLong) > 0) {
  dfBranchLengths <- dfBranchLengths[!dfBranchLengths$species_name %in% dfCheckLong$species_name, ]
}

# Merge dfBranchLengths and dfSeqsNodes.
dfRegression <- merge(dfBranchLengths, dfSeqsNodes, by = "species_name")
# Merge back to the traits dataframe.
dfRegression <- merge(dfRegression, dfTraits, by = "bin_uri")[, c(1:3, 10, 13:22)]
# Renaming for consistency purposes.
colnames(dfRegression)[2] <- "species_name"

### SINGLE VARIABLE REGRESSION ANALYSIS ###
# Running a single variable PGLS regression analysis for each trait to determine whether significance can be detected. If so, they will be included 
# in the multivariable regression model selection process.

# First, make sure the trait data and phylo tree are in the same order.
# Make sure the order of the data matches the tree.
mainTree <- drop.tip(phy = mainTree, tip = mainTree$tip.label[!mainTree$tip.label %in% dfRegression$species_name])
dfRegression <- dfRegression[match(mainTree$tip.label, dfRegression$species_name), ]

## Response variable: Branch length.
# As branch length is our response variable, we will only be estimating Pagel's lambda, which is a measure of phylogenetic signal.
# Use the TestPhyloSig function to estimate phylogenetic signal.
sigBL1 <- TestPhyloSig(dfRegression, "branch_length", mainTree, type = "continuous")

## Control variable: Number of nodes.
# Estimate phylogenetic signal.
sigNodes1 <- TestPhyloSig(dfRegression, "number_of_nodes", mainTree, type = "continuous")
# Now, make a single variable dataframe for number_of_nodes.
dfNodes <- dfRegression[, c("species_name", "branch_length", "number_of_nodes")]
# PGLS. Does the trait have a significant effect on branch length?
# Note: Number of nodes will be included as a control variable in subsequent analyses due to its significant effect on branch length 
# (indicative of the node density effect).
caperNodes1 <- PGLS(dfNodes, mainTree, branch_length ~ number_of_nodes)

############## TRAITS ##############

## 1. MEDIAN LATITUDE.
# Estimate phylogenetic signal.
sigLat1 <- TestPhyloSig(dfRegression, "median_lat", mainTree, type = "continuous")
# Merge single variable dataframe to dfRegression to get branch length data.
dfLatitude <- merge(dfLatitude, dfRegression, by.x = "species_label", by.y = "species_name")[, .(species_name = species_label, branch_length, number_of_nodes, median_lat = median_lat.y)]
# Dealing with ties of species_label from the merging of dataframes since median_lat was determined by bin_uri and not by species name.
dfLatitude <- dfLatitude[!duplicated(species_name)]
# Perform a PGLS analysis.
caperLat1 <- PGLS(dfLatitude, mainTree, branch_length ~ number_of_nodes + median_lat)
# Perform a trait-only analysis as well to see how much variation is explained just by our trait.
caperLat2 <- PGLS(dfLatitude, mainTree, branch_length ~ median_lat)

## 2. MAXIMUM LENGTH.
sigMaxLength <- TestPhyloSig(dfRegression, "max_length", mainTree, type = "continuous")
caperMaxLength1 <- MergeAndPGLS(dfMaxLength, "max_length.x", mainTree, branch_length ~ number_of_nodes + max_length.x)
caperMaxLength2 <- MergeAndPGLS(dfMaxLength, "max_length.x", mainTree, branch_length ~ max_length.x)

## 14. DIET.
sigDiet <- TestPhyloSig(dfRegression, "diet_troph", mainTree, type = "continuous").
caperDiet1 <- MergeAndPGLS(dfDietTroph, "diet_troph.x", mainTree, branch_length ~ number_of_nodes + diet_troph.x)
caperDiet2 <- MergeAndPGLS(dfDietTroph, "diet_troph.x", mainTree, branch_length ~ diet_troph.x)

## 16. SALINITY.
dfSalinity <- as.data.frame(dfSalinity)
dfSalinity$salinity <- relevel(dfSalinity$salinity, ref = "freshwater")
dfRegression$salinity <- relevel(dfRegression$salinity, ref = "freshwater")
sigSalinity <- TestPhyloSig(dfSalinity, "salinity", mainTree, type = "discrete")
caperSalinity1 <- MergeAndPGLS(dfSalinity, "salinity.x", mainTree, branch_length ~ number_of_nodes + salinity.x)
caperSalinity2 <- MergeAndPGLS(dfSalinity, "salinity.x", mainTree, branch_length ~ salinity.x)

## 18. LAKES.
dfLakes <- as.data.frame(dfLakes)
dfLakes$lakes <- relevel(dfLakes$lakes, ref = "0")
dfRegression$lakes <- relevel(dfRegression$lakes, ref = "0")
sigLakes <- TestPhyloSig(dfLakes, "lakes", mainTree, type = "discrete")
caperLakes1 <- MergeAndPGLS(dfLakes, "lakes.x", mainTree, branch_length ~ number_of_nodes + lakes.x)
caperLakes2 <- MergeAndPGLS(dfLakes, "lakes.x", mainTree, branch_length ~ lakes.x)

## 20. OCEANIC.
dfOceanic <- as.data.frame(dfOceanic)
dfOceanic$oceanic <- relevel(dfOceanic$oceanic, ref = "0")
dfRegression$oceanic <- relevel(dfRegression$oceanic, ref = "0")
sigOceanic <- TestPhyloSig(dfOceanic, "oceanic", mainTree, type = "discrete")
caperOceanic1 <- MergeAndPGLS(dfOceanic, "oceanic.x", mainTree, branch_length ~ number_of_nodes + oceanic.x)
caperOceanic2 <- MergeAndPGLS(dfOceanic, "oceanic.x", mainTree, branch_length ~ oceanic.x)

# 21. BENTHIC.
dfBenthic <- as.data.frame(dfBenthic)
dfBenthic$benthic <- relevel(dfBenthic$benthic, ref = "0")
dfRegression$benthic <- relevel(dfRegression$benthic, ref = "0")
sigBenthic <- TestPhyloSig(dfBenthic, "benthic", mainTree, type = "discrete")
caperBenthic1 <- MergeAndPGLS(dfBenthic, "benthic.x", mainTree, branch_length ~ number_of_nodes + benthic.x)
caperBenthic2 <- MergeAndPGLS(dfBenthic, "benthic.x", mainTree, branch_length ~ benthic.x)

### Discrete + biological traits. ###

## 25. BODY SHAPE.

## 30. REPRODUCTIVE MODE.

## 35. PARENTAL CARE.
# Merge single variable dataframe to dfRegression to get branch length data.
dfParentalCare <- merge(dfParentalCare, dfRegression, by = "species_name")
dfParentalCare <- dfParentalCare[, .(species_name, branch_length, number_of_nodes, parental_care = parental_care.y)]
# Relevel the trait according to the category I want to be the reference.
dfParentalCare$parental_care <- relevel(dfParentalCare$parental_care, ref = "none")
# Also relevel the trait in dfRegression.
dfRegression$parental_care <- relevel(dfRegression$parental_care, ref = "none")
# Estimate phylogenetic signal.
# As there are several categories for this trait, I need to create a dataframe 
# of dummy variables first in order to estimate the D metric (it only considers
# binary traits).
dfDummyPc <- model.matrix(~ parental_care - 1, data = dfParentalCare)
dfDummyPc <- as.data.frame(dfDummyPc)
dfDummyPc$species_name <- dfParentalCare$species_name
# Change these according to the categories in your data:
colnames(dfDummyPc)[1] <- "none"
colnames(dfDummyPc)[2] <- "biparental"
colnames(dfDummyPc)[3] <- "maternal"
colnames(dfDummyPc)[4] <- "paternal"
# Now estimate phylogenetic signal separately for each category.
sigNone <- TestPhyloSig(dfDummyPc, "none", mainTree, type = "discrete")
sigBipar <- TestPhyloSig(dfDummyPc, "biparental", mainTree, type = "discrete")
sigMat <- TestPhyloSig(dfDummyPc, "maternal", mainTree, type = "discrete")
sigPat <- TestPhyloSig(dfDummyPc, "paternal", mainTree, type = "discrete")
# Perform a PGLS analysis.
caperParentalCare <- PGLS(dfParentalCare, mainTree, branch_length ~ number_of_nodes + parental_care)
# Since it is a multi-level factor, let's see what the p-value is for the entire trait by performing an ANOVA.
caperControl <- PGLS(dfParentalCare, mainTree, branch_length ~ number_of_nodes)
anova(caperControl, caperParentalCare)

### MULTIVARIABLE REGRESSION ANALYSIS ###
# First, we must find the most complete dataset that we can use for our model
# selection process. So we will remove traits with smaller sample sizes that 
# don't overlap with data from other traits.
# First, get a dataframe of only those traits I am considering.
dfMultivariable <- dfRegression[, c("species_name", "branch_length", "number_of_nodes",  
                                    "median_lat", "longevity", "max_length", 
                                    "salinity", "lakes")]
# Set to datatable.
dfMultivariable <- as.data.table(dfMultivariable)
# We want the most complete dataset possible given our traits.
# First, order the columns by the amount of missing data (NA values).
dfTraitsNA <- sort(dfMultivariable[, lapply(.SD, function(x) sum(is.na(x)))])
# Reorder the original dfTraits. The columns with the least amount of NA values 
# will now be first.
setcolorder(dfMultivariable, names(dfTraitsNA))
# Now I want to loop through the traits, removing one column (trait) at a time
# and count the number of complete cases. This will provide us some information
# as to which traits would provide an adequate sample size for downstream 
# analysis.
# First, take the number of columns in dfMultivariable.
len <- ncol(dfMultivariable)
# Create a temporary variable to hold this number. This variable will hold the 
# number of subsets to check at each iteration.
tempLen <- len
# Create a vector to hold the results of the loop.
all.cc <- NULL
# Start the loop:
for (i in 1:len) {
  # Works best if you set dfMultivariable back to a dataframe.
  x <- as.data.frame(dfMultivariable)
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
# Now, decide where to cut the datatable. (i.e. pick an adequate subset of 
# traits that maximize sample size).
# First, name it according to the trait columns.
names(all.cc) <- rev(colnames(dfMultivariable))
# Look at the results.
all.cc
# What seems like a good cut off point?
which(colnames(dfMultivariable) == "longevity")
dfMultivariableCut <- dfMultivariable[, 1:8] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariableCut <- dfMultivariableCut[complete.cases(dfMultivariableCut)]

# Now check for data variability in this subset so our tests will actually work.
# Examples:
GetTraitInfo(dfMultivariableCut, "branch_length", type = "continuous")
GetTraitInfo(dfMultivariableCut, "median_lat", type = "continuous")
GetTraitInfo(dfMultivariableCut, "salinity", type = "discrete")

# Check for multicollinearity between variables using the variance inflation
# factor (vif), if desired. Multicollinearity can lead to errors in the
# estimations of our coefficients.
fit <- lm(branch_length ~ max_length + median_lat + longevity + lakes + salinity, data = dfMultivariableCut)
vif(mod = fit)

# MODEL SELECTION #
# In this section we will perform manual stepwise model selection. 
# We will remove one trait at a time (usually the one with the highest p-value
# and examine the R squared and BIC values for each model.
# The model with the lowest BIC value and/or highest R squared will serve as
# our best-fit model.

# Now, let's perform a PGLS regression analysis using all of the variables.
# This is our "global" model.
global <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes +  
               median_lat + longevity + salinity + lakes + max_length)
# Check that the phylogenetic residuals are normal.
hist(global$phyres)
qqnorm(global$phyres)
qqline(global$phyres)
plot(x = fitted(global), y = global$phyres, pch = 5)

# Remove variables one at a time see if the fit of the model improves. 
fit1 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + 
               longevity + salinity + lakes + max_length)
# Compare the model to the global model using BIC.
BIC(global, fit1)

# Remove another variable.
fit2 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + 
               longevity + salinity + max_length)
# Compare the model to the global model using BIC.
BIC(fit1, fit2)

# Remove another variable.
fit3 <- PGLS(dfMultivariableCut, mainTree, branch_length ~ number_of_nodes + 
             longevity + max_length)
# Compare the model to the global model using BIC.
BIC(fit2, fit3)

# Continue process until you find the model with the lowest BIC value and/or
# highest R squared value.

# GRAPHS #
install.packages("plotly")
library(plotly)

plotLM <- plot_ly(
  df, x = ~dfColumn1) %>%
  add_markers(y = ~dfColumn2, color = ~dfColumn3, colors = c("#009e73", "#D55E00", "#0072B2", "#CC79A7","#7F7F7F")) %>% 
  add_lines(x = ~dfColumn1, y = fitted(lm(dfColumn2, dfColumn1, data = df)))

# Now you can plot some of your data using ggplot.
# Plot the branch length data.
ggplot(data = dfBranchLengths, aes(log(branch_length))) + 
  labs(x = expression(log(BL[WHOLE])), y = "Frequency") +
  geom_histogram(col="black", 
                 fill="skyblue2") +
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

# Now let's plot some of the trait data against branch length using ggplot.

# Maximum length.
# First calculate a line of best fit.
# Log transforming the length data to improve visualization of the data.
lm(dfRegression$branch_length ~ log(dfRegression$max_length))
# Scatterplot.
ggplot(dfRegression, aes(x = log(max_length), y = branch_length)) + 
  geom_point(size = 2, shape = 19) + geom_abline(intercept = 0.416729, 
                                                 slope = -0.004051, 
                                                 size = 0.75, color = "red") + 
  labs(x = "Log(Maximum length) (cm)", y = expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

# Median latitude.
# First calculate a line of best fit.
lm(dfRegression$branch_length ~ dfRegression$median_lat)
# Scatterplot.
ggplot(dfRegression, aes(x = median_lat, y = branch_length)) + 
  geom_point(size = 2, shape = 19) + geom_abline(intercept = 0.4049763, 
                                                 slope = -0.0001124, 
                                                 size = 0.75, color = "red") + 
  labs(x = "Median latitude", y = expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))