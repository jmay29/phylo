CreateColumn <- function(x, y, z) {
  # Purpose:  Creates a column that contains information on the number of 
  #           unique taxonomic names per BIN (i.e. species, genus, etc.).
  # x = Datatable to apply function to.
  # y = Name of column of interest (i.e. species_name).
  # z = Name of new column (i.e. number_of_species).
  
  # Create a new column with the unique number of names in "y" by BIN.
  x[, new_column := length(unique(y[!is.na(y)])), keyby = bin_uri]
  setnames(x,"new_column", z)
  # Let the user know!
  print("Column created!")
  # Return the updated datatable.
  return(x)
}