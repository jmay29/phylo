GetTraitSpecificDataBIN <- function(x, y) {
  # Purpose: Convenience function for filtering a dataframe for only those rows that have data available for a certain trait.
  # x = Dataframe of species and trait information.
  # y = Column number containing the trait of interest.
  
  # Make sure x is in dataframe format.
  x <- as.data.frame(x)
  # Find rows without data for the trait (NA values).
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  
  # Construct the single-trait datatable.
  # If there are rows without data for column...
  if (length(noY) > 0) {
    # Remove the species without data for column.
    z <- x[-noY, ]
    # Reorganize dataframe.
    # Column 1 = bin_uri, column 10 = species_label
    z <- z[c(10, y)]
    # Remove duplicate entries per BIN.
    z <- z[!duplicated(z$species_label), ] 
    
    # If all rows have data for column...
  } else {
    # If no entries need to be removed, just rename and reorganize x.
    z <- x[c(10, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$species_label), ]  
  }
  
  rm(noY)
  return(z)
  
}
