CountConflicts <- function(x, y) {
  # Purpose:  Counts how many BIN conflicts there for a given taxonomic level and 
  #           returns a vector of the BIN names.
  # x = Dataframe or datatable to apply function to.
  # y = Name of column of interest (the column that contains counts of order/
  #     family/genus/species per BIN).

  # Which BINs have a value greater than 1/
  conflicts <- which(x[, get(y)] > 1)
  # Subset these BINs.
  conflicts <- x[conflicts, ]
  # We only want unique BIN URIs.
  conflicts <- unique(conflicts$bin_uri)
  # How many are there?
  num_conflicts <- length(conflicts)
  # Let the user know!
  print(paste0("There are ", num_conflicts, " BIN conflicts!"))
  # Return a vector of the names of conflicted BINs.
  return(conflicts)
}
