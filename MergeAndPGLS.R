MergeAndPGLS <- function(x, y, phylo, formula) {
  x <- setDT(x)
  x <- merge(x, dfRegression, by = "species_name")
  x <- x[, .(species_name, branch_length, number_of_nodes, get(y))]
  setnames(x, "V4", y)
  x <- as.data.frame(x)
  c_data <- comparative.data(phylo, x, "species_name", vcv = TRUE)
  caper <- pgls(formula, c_data)
  return(caper)
}
