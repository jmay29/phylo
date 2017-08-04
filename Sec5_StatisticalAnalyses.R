################################################################################

# # Copyright Jacqueline May (C) 2017.
# Project: MSc of Bioinformatics Thesis.
# Program Description: Multivariate analysis of environmental and biological 
#                      correlates affecting fish molecular evolution rates.

# Advisors: Dr. Sarah J. Adamowicz and Dr. Zeny Feng.

################################################################################

##### SECTION 5: STATISTICAL ANALYSES #####
# This is the section where the univariate analyses are performed as a form
# of model selection.

### PACKAGES REQUIRED ###

# For data manipulation:
#install.packages("data.table")
library(data.table)
#install.packages("dplyr")
library(dplyr)

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
library(ggplot)
#install.packages("gtools")
library(gtools)

# Load the function(s) that I designed for this script:
source("TestPhyloSig.R")

################################################################################

# First, read in your phylogenetic tree.
mainTree <- read.tree(file = "RAxML_labelledTree.EPAPercTreeAug2")
mainTree$tip.label <- gsub("_", " ", mainTree$tip.label)
# Re-name the outgroup.
mainTree$tip.label <- gsub("QUERY   ", "", mainTree$tip.label)
# Root the tree.
outgroups <- c("Acanthurus triostegus")
mainTree <- root(mainTree, outgroup = outgroups, resolve.root = TRUE)

### TRAIT: NUMBER OF NODES.
# Let's first determine the number of nodes for each species. This will be used
# as a control variable in the multiple regression analysis (to account for the
# node density effect).
number_of_nodes <- distRoot(mainTree, tips = "all", method = "nNodes")
# Make a dataframe of this information
dfNumberOfNodes <- data.frame(number_of_nodes)
dfNumberOfNodes$species_name <- row.names(dfNumberOfNodes)
# Merge with dfAllSeqsNode.
dfAllSeqsNode <- merge(dfAllSeqsWithOG, dfNumberOfNodes, by = "species_name")

### TRAIT: BRANCH LENGTHS.
# Let's calculate the sum of branch lengths now (from root to tip). These values
# will serve as our measurement of molecular evolution rate.
branch_lengths <- distRoot(mainTree, tips = "all", method = "patristic")
# Make into a dataframe.
dfBranchLengths <- data.frame(branch_lengths)
colnames(dfBranchLengths)[1] <- "branch_length"
dfBranchLengths$species_name <- names(branch_lengths)
# Check the distribution of the branch lengths.
hist(dfBranchLengths$branch_length, main = "", xlab = "Branch Length")
# Some more stats.
median(dfBranchLengths$branch_length)
mean(dfBranchLengths$branch_length)
range(dfBranchLengths$branch_length)
# Merge dfBranchLengths and dfAllSeqsNode.
dfRegression <- merge(dfBranchLengths, dfAllSeqsNode, by = "species_name")
# Merge back to the traits dataframe.
dfRegression <- merge(dfRegression, dfTraits, all.x = T, by = "bin_uri")
# Some reordering.
dfRegression <- dfRegression[, c(1:2, 4, 3, 11, 15:23)]
# Renaming for consistency purposes.
colnames(dfRegression)[2] <- "species_name"
colnames(dfRegression)[3] <- "filtered_bin_size"

### Univariate analyses ###
# Running a univariate PGLS regression analysis for each trait to determine
# whether significance can be detected. If so, they will be included in the 
# multiple regression model selection process.

# First, make sure the trait data and phylo tree are in the same order.
# Make sure the order of the data matches the tree.
mainTree <- drop.tip(phy = mainTree, 
                     tip = mainTree$tip.label[!mainTree$tip.label%in%dfRegression$species_name], 
                     rooted = T)
dfRegression <- dfRegression[match(mainTree$tip.label, 
                                   dfRegression$species_name), ]

## 1. Branch length.
# As branch length is our response variable, we will only be estimating Pagel's 
# lambda, which is a measure of phylogenetic signal.
# Use the TestPhyloSig function to estimate phylogenetic signal.
sigBL <- TestPhyloSig(dfRegression, "branch_length", mainTree, 
                      type = "continuous")

## 2. Number of nodes.
# Estimate phylogenetic signal.
sigNodes <- TestPhyloSig(dfRegression, "number_of_nodes", mainTree, 
                         type = "continuous")
# Now, make a univariate dataframe for number_of_nodes.
dfNodes <- dfRegression[, c("species_name", "branch_length", "number_of_nodes")]
# Perform a PGLS analysis using caper.
# First, we must make a comparative.data object to ensure the tree tips and data
# match.
c_data <- comparative.data(mainTree, dfNodes, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
# Note: Number of nodes will be included as a control variable in
# subsequent analyses due to its significant effect on branch length 
# (indicative of the node density effect).
caperNodes <- pgls(branch_length ~ number_of_nodes, c_data, lambda = "ML")

## 3. Median latitude.
# Estimate phylogenetic signal.
sigLat <- TestPhyloSig(dfRegression, "median_lat", mainTree, 
                       type = "continuous")
# Merge univariate dataframe to dfRegression to get branch length data.
dfLatitude <- merge(dfLatitude, dfRegression, by.x = "species_label", 
                    by.y = "species_name")
dfLatitude <- dfLatitude[, .(species_name = species_label, branch_length, 
                             number_of_nodes, median_lat = median_lat.y)]
# Dealing with ties of species_label from the merging of dataframes 
# since median_lat was determined by bin_uri and not by species name.
dup_majority_species <- which(duplicated(dfLatitude$species_name))
dfLatitude <- dfLatitude[-dup_majority_species,]
# Make sure it is of the dataframe format.
dfLatitude <- as.data.frame(dfLatitude)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfLatitude, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperLatitude <- pgls(branch_length ~ number_of_nodes + median_lat, c_data, lambda = "ML")

## 4. Maximum length.
# Estimate phylogenetic signal.
sigLength <- TestPhyloSig(dfRegression, "max_length", mainTree, 
                          type = "continuous")
# Merge univariate dataframe to dfRegression to get branch length data.
dfMaxLength <- merge(dfMaxLength, dfRegression, by = "species_name")
dfMaxLength <- dfMaxLength[, .(species_name, branch_length, number_of_nodes, 
                               max_length = max_length.y)]
# Make sure it is of the dataframe format.
dfMaxLength <- as.data.frame(dfMaxLength)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfMaxLength, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperMaxLength <- pgls(branch_length ~ number_of_nodes + max_length, c_data, lambda = "ML")

## 5. Maximum depth.
# Estimate phylogenetic signal.
sigDepth <- TestPhyloSig(dfRegression, "max_depth", mainTree, 
                         type = "continuous")
# Merge univariate dataframe to dfRegression to get branch length data.
dfMaxDepth <- merge(dfMaxDepth, dfRegression, by = "species_name")
dfMaxDepth <- dfMaxDepth[, .(species_name, branch_length, number_of_nodes, 
                             max_depth = max_depth.y)]
# Make sure it is of the dataframe format.
dfMaxDepth <- as.data.frame(dfMaxDepth)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfMaxDepth, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperMaxDepth <- pgls(branch_length ~ number_of_nodes + max_depth, c_data, lambda = "ML")

## 6. Longevity.
# Estimate phylogenetic signal.
sigLong <- TestPhyloSig(dfRegression, "longevity", mainTree, 
                        type = "continuous")
# Merge univariate dataframe to dfRegression to get branch length data.
dfLongWild <- merge(dfLongWild, dfRegression, by = "species_name")
dfLongWild <- dfLongWild[, .(species_name, branch_length, number_of_nodes, 
                             longevity = longevity.y)]
# Make sure it is of the dataframe format.
dfLongWild <- as.data.frame(dfLongWild)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfLongWild, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperLongevity <- pgls(branch_length ~ number_of_nodes + longevity, c_data, lambda = "ML")

## 7. Age at maturity.
# Estimate phylogenetic signal.
sigAge <- TestPhyloSig(dfRegression, "age_at_maturity", mainTree, 
                       type = "continuous")
# Merge univariate dataframe to dfRegression to get branch length data.
dfAgeMaturity <- merge(dfAgeMaturity, dfRegression, by = "species_name")
dfAgeMaturity <- dfAgeMaturity[, .(species_name, branch_length, number_of_nodes, 
                                   age_at_maturity = age_at_maturity.y)]
# Make sure it is of the dataframe format.
dfAgeMaturity <- as.data.frame(dfAgeMaturity)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfAgeMaturity, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperAgeMaturity <- pgls(branch_length ~ number_of_nodes + age_at_maturity, c_data, lambda = "ML")

# Discrete traits.
## 8. Streams.
# Merge univariate dataframe to dfRegression to get branch length data.
dfStreams <- merge(dfStreams, dfRegression, by = "species_name")
dfStreams <- dfStreams[, .(species_name, branch_length, number_of_nodes, 
                           streams = streams.y)]
# Make sure it is of the dataframe format.
dfStreams <- as.data.frame(dfStreams)
# Relevel the trait according to the category I want to be the reference.
dfStreams$streams <- relevel(dfStreams$streams, ref = "0")
# Also relevel the trait in dfRegression.
dfRegression$streams <- relevel(dfRegression$streams, ref = "0")
# Estimate phylogenetic signal.
sigStreams <- TestPhyloSig(dfStreams, "streams", mainTree, type = "discrete")
# Make sure it is of the dataframe format.
dfStreams <- as.data.frame(dfStreams)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfStreams, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperStreams <- pgls(branch_length ~ number_of_nodes + streams, c_data, lambda = "ML")

## 9. Lakes.
# Merge univariate dataframe to dfRegression to get branch length data.
dfLakes <- merge(dfLakes, dfRegression, by = "species_name")
dfLakes <- dfLakes[, .(species_name, branch_length, number_of_nodes, 
                       lakes = lakes.y)]
# Make sure it is of the dataframe format.
dfLakes <- as.data.frame(dfLakes)
# Relevel the trait according to the category I want to be the reference.
dfLakes$lakes <- relevel(dfLakes$lakes, ref = "0")
# Also relevel the trait in dfRegression.
dfRegression$lakes <- relevel(dfRegression$lakes, ref = "0")
# Estimate phylogenetic signal.
sigLakes <- TestPhyloSig(dfLakes, "lakes", mainTree, type = "discrete")
# Make sure it is of the dataframe format.
dfLakes <- as.data.frame(dfLakes)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfLakes, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperLakes <- pgls(branch_length ~ number_of_nodes + lakes, c_data, lambda = "ML")

## 10. Body Shape.
# Merge univariate dataframe to dfRegression to get branch length data.
dfBodyShape <- merge(dfBodyShape, dfRegression, by = "species_name")
dfBodyShape <- dfBodyShape[, .(species_name, branch_length, number_of_nodes, 
                               body_shape = body_shape.y)]
# Relevel the trait according to the category I want to be the reference.
dfBodyShape$body_shape <- relevel(dfBodyShape$body_shape, ref = "fusiform / normal")
# Also relevel the trait in dfRegression.
dfRegression$body_shape <- relevel(dfRegression$body_shape, ref = "fusiform / normal")
# Estimate phylogenetic signal.
# As there are several categories for this trait, I need to create a dataframe 
# of dummy variables first in order to estimate the D metric (it only considers
# binary traits).
dfDummyBs <- model.matrix(~ body_shape - 1, data = dfBodyShape)
dfDummyBs <- as.data.frame(dfDummyBs)
dfDummyBs$species_name <- dfBodyShape$species_name
colnames(dfDummyBs)[1] <- "elongated"
colnames(dfDummyBs)[2] <- "fusiform"
colnames(dfDummyBs)[3] <- "short"
# Now estimate phylogenetic signal separately for each category.
sigElongated <- TestPhyloSig(dfDummyBs, "elongated", mainTree, type = "discrete")
sigFusiform <- TestPhyloSig(dfDummyBs, "fusiform", mainTree, type = "discrete")
sigShort <- TestPhyloSig(dfDummyBs, "short", mainTree, type = "discrete")
# Make sure it is of the dataframe format.
dfBodyShape <- as.data.frame(dfBodyShape)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfBodyShape, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperBodyShape <- pgls(branch_length ~ number_of_nodes + body_shape, c_data, lambda = "ML")
# Since it is a multi-level factor, let's see what the p-value is for the entire trait
# by performing an ANOVA.
caperControl <- pgls(branch_length ~ number_of_nodes, c_data, lambda = "ML")
anova(caperControl, caperBodyShape)

## 11. Feeding Type.
# Merge univariate dataframe to dfRegression to get branch length data.
dfFeedingType <- merge(dfFeedingType, dfRegression, by = "species_name")
dfFeedingType <- dfFeedingType[, .(species_name, branch_length, number_of_nodes, 
                                   feeding_type = feeding_type.y)]
# Relevel the trait according to the category I want to be the reference.
dfFeedingType$feeding_type <- relevel(dfFeedingType$feeding_type, ref = "hunting macrofauna (predator)")
# Also relevel the trait in dfRegression.
dfRegression$feeding_type <- relevel(dfRegression$feeding_type, ref = "hunting macrofauna (predator)")
# Estimate phylogenetic signal.
# As there are several categories for this trait, I need to create a dataframe 
# of dummy variables first in order to estimate the D metric (it only considers
# binary traits).
dfDummyFt <- model.matrix(~ feeding_type - 1, data = dfFeedingType)
dfDummyFt <- as.data.frame(dfDummyFt)
dfDummyFt$species_name <- dfFeedingType$species_name
colnames(dfDummyFt)[1] <- "hunting"
colnames(dfDummyFt)[2] <- "browsing"
colnames(dfDummyFt)[3] <- "grazing"
colnames(dfDummyFt)[4] <- "selective"
colnames(dfDummyFt)[5] <- "variable"
# Now estimate phylogenetic signal separately for each category.
sigHunting <- TestPhyloSig(dfDummyFt, "hunting", mainTree, type = "discrete")
sigBrowsing <- TestPhyloSig(dfDummyFt, "browsing", mainTree, type = "discrete")
sigGrazing <- TestPhyloSig(dfDummyFt, "grazing", mainTree, type = "discrete")
sigSelective <- TestPhyloSig(dfDummyFt, "selective", mainTree, type = "discrete")
sigVariable <- TestPhyloSig(dfDummyFt, "variable", mainTree, type = "discrete")
# Make sure it is of the dataframe format.
dfFeedingType <- as.data.frame(dfFeedingType)
# Perform a PGLS analysis using caper.
c_data <- comparative.data(mainTree, dfFeedingType, "species_name", vcv = TRUE)
# PGLS. Does the trait have a significant effect on branch length?
caperFeedingType <- pgls(branch_length ~ number_of_nodes + feeding_type, c_data, lambda = "ML")
caperControl <- pgls(branch_length ~ number_of_nodes, c_data, lambda = "ML")
anova(caperControl, caperFeedingType)

### Multivariate (multiple regression) analyses ###
# Manual stepwise model selection. Lowest BIC value.
# First, get a dataframe of only those traits I am considering.
dfMultivariate <- dfRegression[, c("species_name", "branch_length", "number_of_nodes", 
                                   "median_lat", "longevity", "max_length", "max_depth",
                                   "age_at_maturity")]
# Set to datatable.
dfMultivariate <- setDT(dfMultivariate)
# We want the most complete dataset possible given our traits.
# First, order the columns by the amount of missing data (NA values).
dfTraitsNA <- sort(dfMultivariate[, lapply(.SD, function(x) sum(is.na(x)))])
# Reorder the original dfTraits. The columns with the least amount of NA values 
# are now first.
setcolorder(dfMultivariate, names(dfTraitsNA))
# Now, count the number of NAs in each row.
dfMultivariate[, count := rowSums(is.na(dfMultivariate))]
# Sort the rows by this order.
dfMultivariate <- dfMultivariate[order(count), ]
# Remove this count column as it is no longer needed.
dfMultivariate[, count := NULL]
# The datatable is now sorted: columns with the least amount of NAs to the left, 
# rows with with the least amount of NAs on top.
# Now I want to loop through the traits, removing one column (trait) at a time
# and count the number of complete cases. This will provide us some information
# as to which traits would provide an adequate sample size for downstream analysis.
# First, take the number of columns in dfMultivariate.
len <- ncol(dfMultivariate)
# Create a temporary variable to hold this number. This variable will hold the 
# number of subsets to check at each iteration.
tempLen <- len
# Create a vector to hold the results of the loop.
all.cc <- NULL
# Start the loop:
for (i in 1:len) {
  # Works best if you set dfMultivariate back to a dataframe.
  x <- as.data.frame(dfMultivariate)
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
names(all.cc) <- rev(colnames(dfMultivariate))
# Look at the results.
all.cc
# max_depth seems like a good cutoff point.
which(colnames(dfMultivariate) == "max_depth")
dfMultivariateCut <- dfMultivariate[, 1:6] 
# Finally, filter the original dfTraits datatable so only complete cases are kept.
dfMultivariateCut <- dfMultivariateCut %>% filter(complete.cases(.))

# Now check for data variability in this subset.
hist(dfMultivariateCut$branch_length)
range(dfMultivariateCut$branch_length)

hist(dfMultivariateCut$number_of_nodes)
range(dfMultivariateCut$number_of_nodes)

hist(dfMultivariateCut$median_lat)
range(dfMultivariateCut$median_lat)

hist(dfMultivariateCut$max_length)
range(dfMultivariateCut$max_length)

hist(dfMultivariateCut$max_depth)
range(dfMultivariateCut$max_depth)

# Check for multicollinearity between variables using the variance inflation
# factor (vif), if desired.
lm.res <- lm(branch_length ~ max_length + max_depth + median_lat, data = dfMultivariateCut)
vif(mod = lm.res)

# Backward selection using BIC.
# First, create the comparative.data object.
c_data <- comparative.data(mainTree, dfMultivariateCut, names.col = "species_name", 
                           vcv = TRUE)
# Now, let's perform a PGLS regression analysis using all of the variables.
full <- pgls(branch_length ~ number_of_nodes + max_depth
             + max_length + median_lat, data = c_data, lambda = "ML")
# Inspect for homogeneity.
hist(full$phyres)
qqnorm(full$phyres)
qqline(full$phyres)
plot(x = fitted(full), y = full$phyres, pch = 5)

# Let's remove max_depth first.
fit1 <- pgls(branch_length ~ number_of_nodes + median_lat
             + max_length, data = c_data, lambda = "ML")
# Compare the model to the full model using ANOVA and BIC.
anova(full, fit1)
BIC(full, fit1)

# Remove median_lat.
fit2 <- pgls(branch_length ~ number_of_nodes + max_length, 
             data = c_data, lambda = "ML")
# Compare the model to the full model using ANOVA and BIC.
anova(fit1, fit2)
BIC(fit1, fit2)

# Now let's plot some of the trait data against branch length using ggplot.

# Maximum length.
# First calculate a line of best fit.
lm(dfMaxLength$branch_length ~ log(dfMaxLength$max_length))
# Scatterplot.
ggplot(dfMaxLength, aes(x = log(max_length), y = branch_length)) + 
  geom_point(size = 2, shape = 19) + geom_abline(intercept = 0.54564, 
                                                 slope = -0.01343, size = 0.75, color = "red") + 
  labs(x = "Log(Maximum length (cm))", y = expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

# Age at maturity.
lm(dfAgeMaturity$branch_length ~ log(dfAgeMaturity$age_at_maturity))
ggplot(dfAgeMaturity, aes(x = log(age_at_maturity), y = branch_length)) 
+ geom_point(size = 2, shape = 19) + geom_abline(intercept = 0.48558, 
                                                 slope = -0.02061, size = 0.75, color = "red") + 
  labs(x = "Log(Maximum length (cm))", y = expression(BL[WHOLE])) + 
  theme(panel.background = element_rect(fill = "grey95", colour = 'royalblue4'))

# Save the image in your current directory.
png(filename="ageVsBranchLength.png", 
    units="in", 
    width=7, 
    height=5, 
    pointsize=12, 
    res=300)
my_sc_plot(data)
dev.off()