# As I will be estimating the phylogenetic signal of traits multiple times, 
# here is a function that will make the process more efficient. It can be used
# on either continuous or discrete (categorical) traits.
## FUNCTION: TestPhyloSig ##
# Purpose: Estimates phylogenetic signal of the trait given an input dataframe
# and phylogenetic tree.
TestPhyloSig <- function(x, y, phylo, type = c("continuous", "discrete")){
  # x = Dataframe to apply function to.
  # y = Trait of interest.
  # phylo = Tree of class phylo. This is needed to estimate phylogenetic signal.
  # type = Specify whether the trait is continuous or discrete.
  
  type <- match.arg(type)
  # If the trait is continuous...
  if(type == "continuous") {
    # Make sure x is a dataframe and not a datatable.
    x <- as.data.frame(x)
    # Subset the trait from the dataframe.
    trait <- x[, y]
    # Ensure that the order of the data matches the tree.
    x <- x[match(phylo$tip.label, x$species_name), ]
    # Name it according to the tips of the provided tree.
    names(trait) <- phylo$tip.label
    # Estimate Pagel's lambda using the phylosig function from the "phytools"
    # package. A test of significance is also performed by setting test = TRUE.
    sig <- phylosig(phylo, trait, method = "lambda", test = TRUE)
    
  } else if(type == "discrete") {
    # Make sure x is a dataframe and not a datatable.
    x <- as.data.frame(x)
    # First prune the tree so that only tips in the dataframe are present. This
    # is necessary for the phylo.d function to work properly.
    prunedPhylo <- drop.tip(phy = phylo, 
                            tip = phylo$tip.label[!phylo$tip.label%in%x$species_name], 
                            rooted = T)
    # Make sure the dataframe matches the order of the tree's tip labels.
    x <- x[match(prunedPhylo$tip.label, x$species_name), ]
    # Estimate the D metric (measure of phylogenetic signal for binary traits)
    # by Fritz and Purvis (2010) using the phylo.d function from the "caper" 
    # package.
    # Temporarily rename column as "y" so that phylo.d function works (as it 
    # take y as literal name of column).
    names(x)[names(x) == y] <- "y"
    sig <- phylo.d(x, prunedPhylo, names.col = species_name, binvar = y)
  }
  return(sig)
}