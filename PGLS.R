PGLS <- function(x, phylo, formula) {
  x <- as.data.frame(x)
  c_data <- comparative.data(phylo, x, "species_name", vcv = TRUE)
  caper <- pgls(formula, c_data)
  return(caper)
}
