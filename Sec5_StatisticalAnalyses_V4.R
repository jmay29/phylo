# Copyright (C) 2020 Jacqueline May.
# Program Description: Multivariable analysis of environmental and biological correlates affecting molecular evolution rates.

# Contributions & Acknowledgements #
# Dr. Sarah J. Adamowicz and Dr. Zeny Feng for help with designing and structuring the pipeline.

# Adapted lines 87-96 for using CV to select lambda value:
# Paradis, E. (2012). Analysis of Phylogenetics and Evolution with R. Second Edition. New York, New York: Springer. Pgs 184-185, 198-199.

# Adapted lines 106-112 for extracting terminal rate data:
# http://blog.phytools.org/2016/02/extracting-terminal-edge-lengths-for.html
# Paradis, E. (2012). Analysis of Phylogenetics and Evolution with R. Second Edition. New York, New York: Springer. Pgs 198-199.

# Adapted lines 87-91 and 147-149 from code shared in Stack Overflow discussion:
# Author: https://stackoverflow.com/users/1312519/by0.
# https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r.
# Author: https://stackoverflow.com/users/580110/henk.
# https://stackoverflow.com/questions/24016612/subset-of-data-frame-columns-to-maximize-complete-observations.

# Used as general guide for PGLS analyses (lines 238-254):
# Mundry, R. (2014). Statistical Issues and Assumptions of Phylogenetic Generalized Least Squares. In L.Z. Garamszegi (Ed.), 
# Modern Phylogenetic Comparative Methods and Their Application in Evolutionary Biology. Online Practice Materials.
# URL: http://www.mpcm-evolution.org/practice/online-practical-material-chapter-6/chapter-6-1exercises-testing-assumptions-statistical-issues-framework-phylogenetic-generalized-least-squares

# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# There is a copy of the GNU General Public License along with this program in the repository where it is located. Or view it directly here at http://www.gnu.org/licenses/

#################################################################################################################

#### SECTION 5: STATISTICAL ANALYSES ####
# This section is designed to perform single variable and multivariable analyses regression analyses while controlling for phylogeny. The main objective is to identify those variables that contribute most to variation in molecular evolution rate.

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

################################################################################################

#### PREPARING RATE DATA ####

# A phylogenetic tree containing branch length data for your species is required for this section.
# Read in your phylogenetic tree.
mainTree <- read.tree(file = "fish_tree.tre")
# Fixing the tip labels (if the species names contain underscores when they are read in).
mainTree$tip.label <- gsub("_", " ", mainTree$tip.label)
# Match mainTree with data subset. This will ensure the tree has only the tips we need for data analysis.
mainTree <- drop.tip(phy = mainTree, tip = mainTree$tip.label[!mainTree$tip.label %in% dfCentroidSeqsWO$species_name])
# Root the tree using your chosen outgroup species.
mainTree <- root(mainTree, outgroup = outgroups, resolve.root = T)

# Calibrating tree to obtain absolute rates of substitution.
# Get the node number of the most recent common ancestor (or you could just specify the root in the chronos function).
nodeNum <- getMRCA(mainTree, tip = c(outgroups, "Trimma hoesei"))
# Using known fossil information (e.g. from https://fishtreeoflife.org/), set the estimated minimum age.
minAge <- c(49)
# Set the estimated maximum age.
maxAge <- c(93.62)

# Cross validation for selecting lambda.
# Create empty list to store CV results
cvResults <- vector("list", 9)
# Iterate through the different lambda values.
for (i in -4:4)
  cvResults[[i + 5]] <- chronopl(mainTree, 10^i, age.min = minAge, age.max = maxAge, node = nodeNum, CV = TRUE)
# Plot the sum of the CV values against lambda value.
lambda <- 10^(-4:4)
CV <- sapply(cvResults, function(x) sum(attr(x, "D2")))
# What is the lambda value that minimizes the cross-validation criterion?
plot(lambda, CV / 1e5, log = "x")

# Using the chosen lambda value to build the tree using chronos function (new version of chronopl).
# Making the calibration object.
mycalibration <- makeChronosCalib(mainTree, node = nodeNum, age.min = minAge, age.max = maxAge)
# Calibrate the tree.
fish.chrono <- chronos(mainTree, calibration = mycalibration, control = chronos.control())
# Extract the rates.
fishRates <- attr(fish.chrono, "rates")
summary(fishRates)
# Extract the tip labels.
tips <- fish.chrono$tip.label
# Get the node numbers of the tips.
nodes <- sapply(tips, function(x, y) which(y == x), y = fish.chrono$tip.label)
# Then get the rates for those nodes.
edge.lengths <- setNames(fish.chrono$edge.length[sapply(nodes, function(x, y) 
  which(y == x), y = fish.chrono$edge[, 2])], names(nodes))
# Extract the rates that correspond to the terminal edge lengths.
terminal_rate <- setNames(fishRates[sapply(nodes, function(x, y) which(y == x), y = fish.chrono$edge[, 2])], names(nodes))
# Plot the rates
histogram(terminal_rate)
# Put these values into a dataframe.
dfRates <- as.data.frame(terminal_rate)
dfRates$species_name <- row.names(dfRates)
# Get info about the rates.
GetTraitInfo(dfRates$terminal_rate)
# Range within which 95% of the values fall.
quantile(dfRates$terminal_rate, probs = c(.025, .975))
# Save chronogram to file.
write.tree(fish.chrono, "chronoTree.tre")

# Take a closer look at outliers. Some contaminated sequences might have STILL gotten through, so it is best to check!
# Using the IQR to detect statistical outliers.
lowerQuantile <- quantile(dfRates$terminal_rate)[2]
upperQuantile <- quantile(dfRates$terminal_rate)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 3) + upperQuantile
lowerThreshold <-  lowerQuantile - (iqr * 3)
# Extreme slow rates.
dfSlow <- dfRates[terminal_rate < lowerThreshold]
# Get the sequence information in case you want to BLAST the sequence (also, we aren't interested in outgroup species here, that's why we are using dfCentroidSeqsNO).
#dfSlow <- merge(dfSlow, dfCentroidSeqsNO, by = "species_name")
# Do the same for the extremely fast rates.
dfFast <- dfRates[dfRates$terminal_rate > upperThreshold, ]
# Remove from dataset, if desired.

#### SINGLE VARIABLE REGRESSION ANALYSIS ####
# Running a single variable PGLS regression analysis for each trait to determine whether significance can be detected. If so, they will be included in the multivariable regression model selection process.

# Read the tree in that you are using for PGLS analyses.
chronoTree <- read.tree("chronoTree.tre")
# Fixing the tip labels (if the species names contain underscores when they are read in).
chronoTree$tip.label <- gsub("_", " ", chronoTree$tip.label)
# First, make sure the trait data and phylo tree match (in case species were removed).
chronoTree <- drop.tip(phy = chronoTree, tip = chronoTree$tip.label[!chronoTree$tip.label %in% dfRates$species_name])
dfRates <- dfRates[match(chronoTree$tip.label, dfRates$species_name), ]
# Merge trait data and with rate data.
dfTraits <- merge(dfRates, dfTraits, by = "species_name")

### SINGLE-VARIABLE PGLS ANALYSES ###
# Use the PGLS function to perform single-variable for all of the traits. 
# We will do this by looping through all of the columns containing the trait data using lapply.
traits <- as.list(colnames(dfTraits[, 5:14]))
# Set to dataframe.
dfTraits <- as.data.frame(dfTraits)
# Start the loop.
singleVarResults <- lapply(traits, function(x) {
  # We only want the columns containing species name and dependent and independent variables.
  data <- dfTraits[, c("species_name", x, "terminal_rate")]
  # Remove NA values.
  data <- data[complete.cases(data), ]
  # Perform PGLS. The trait of interest in this case will always be the 2nd column.
  caper <- PGLS(data, chronoTree, log(terminal_rate) ~ data[, 2])
  # Take the summary of the results.
  caperSum <- summary(caper)
})
# Assign names to the list of results based on the trait of interest.
names(singleVarResults) <- traits

# Which traits have p-values 0.15 or below?
# For now, this is only taking the first p-value of the trait (I still need to change it to deal with multi-level factors).
sigVars <- lapply(singleVarResults, function(x) (x$coefficients[2,4]))
names(sigVars) <- names(singleVarResults)
# Which are below 0.15?
keepVars <- names(which(sigVars <= 0.15))
keepVars <- c("species_name", "terminal_rate", keepVars)

#### MULTIVARIABLE REGRESSION ANALYSES ####
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
len <- which(colnames(dfMultivariable) == "repro_mode")
dfMultivariableCut <- dfMultivariable[, 1:len] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariableCut <- dfMultivariableCut[complete.cases(dfMultivariableCut)]

# Now check for data variability in this subset so our tests will actually work.
# Examples:
GetTraitInfo(dfMultivariableCut$median_lat)
GetTraitInfo(dfMultivariableCut$lakes)
GetTraitInfo(dfMultivariableCut$repro_mode)

# Check for multicollinearity between variables using the variance inflation factor (vif), if desired. 
# Multicollinearity can lead to errors in the estimations of our coefficients.
fit <- lm(log(terminal_rate) ~ median_lat + lakes + repro_mode, data = dfMultivariableCut)
vif(mod = fit)

#### MODEL SELECTION ####
# In this section we will perform manual stepwise model selection. We will remove one trait at a time (usually the one with the highest p-value and examine the R squared and BIC values for each model. The model with the lowest BIC value and/or highest R squared will serve as our best-fit model.

# Now, let's perform a PGLS regression analysis using all of the variables.
# This is our "global" model.
global <- PGLS(dfMultivariableCut, chronoTree, log(terminal_rate) ~ median_lat + lakes + repro_mode)
# Check that the phylogenetic residuals are normal.
hist(global$phyres)
qqnorm(global$phyres)
qqline(global$phyres)
plot(x = fitted(global), y = global$phyres, pch = 5)
summary(global)

# Remove variables one at a time see if the fit of the model improves. 
fit1 <- PGLS(dfMultivariableCut, chronoTree, log(terminal_rate) ~ lakes + repro_mode)
# Compare the model to the global model using BIC.
BIC(global, fit1)

# Continue process until you find the model with the lowest BIC value and/or highest R squared value.

# GRAPHS #
# For example, plot terminal_rate against median_lat.
# Note: this data is uncorrected for phylogeny at the moment.
fit <- lm(log(terminal_rate) ~ median_lat, data = dfMultivariableCut)
dfMultivariableCut %>% plot_ly(x = ~median_lat) %>% 
  add_markers(y = ~log(terminal_rate)) %>% 
  add_lines(x = ~median_lat, y = fitted(fit))