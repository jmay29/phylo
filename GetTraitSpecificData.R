# Function: GetTraitSpecificData ##
# Purpose: Filters a dataframe for presence of trait data. This function is 
#          for the single variable analyses section.
GetTraitSpecificData <- function(x, y) {
  # Filters a dataframe for only those data related to a specified trait.
  # x = Dataframe of species and trait information.
  # y = Trait (column) name of interest.
  
  # Convert to dataframe format.
  x <- as.data.frame(x)
  # Find the rows without data for the specified trait.
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  # Construct the single variable trait dataframe.
  # If there are rows without data for the trait...
  if (length(noY) > 0) {
    # Remove the species without data for the trait from the dataframe.
    z <- x[-noY, ]
    # Reorganize the datatable.
    # 1 = species_name
    z <- z[c(1, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$species_name), ] 
    # If all rows have data for column...
  } else {
    # If no entries need to be removed, just rename and reorganize x.
    # 1 = species_name
    z <- x[c(1, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$species_name), ]  
  }
  # Return the filtered dataframe.
  return(z)
}
