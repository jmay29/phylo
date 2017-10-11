AssignLabel <- function(x, y, z) {
  # Purpose: Assigns a taxonomic label to each BIN and returns a datatable.
  # x = Dataframe of species information.
  # y = Column containing taxonomic assignments (i.e. species_name).
  # z = Name of new column.
  
  # We are only looking at rows with taxonomic information.
  dfLabel <- x[grep("[A-Z]", get(y))]
  # Count the number of names per BIN.
  dfLabel <- dfLabel[, .N, by = .(bin_uri, get(y))]
  # Select the majority name in the BIN.
  dfLabel <- dfLabel[order(-N), .(new = get[1L]), keyby = bin_uri]
  # Assign the new name of the column.
  setnames(dfLabel,"new", z)
  
  # Return the datatable of label information.
  return(dfLabel)
}
