## FUNCTION: GetTraitSpecificDataBIN ##
# Purpose: Filters for presence of trait data. This function is for the 
#          univariate analyses section.
GetTraitSpecificDataBIN <- function(x, y) {
  # Filters a dataframe for only those rows that have data available for a 
  # certain trait.
  # x = Dataframe of species and trait information.
  # y = Trait of interest.
  
  # Make sure x is in dataframe format.
  x <- as.data.frame(x)
  # Find rows without data for the trait (NA values).
  noY <- is.na(x[, y])
  noY <- which(noY == "TRUE")
  # Construct the univariate trait datatable. This datatable will be used in the 
  # eventual univariate analysis.
  # If there are rows without data for column...
  if (length(noY) > 0) {
    # Remove the species without data for column.
    z <- x[-noY, ]
    # Reorganize datatable.
    # Column 1 = bin_uri, column 10 = species_label, 
    # column 13 = filtered_bin_size
    z <- z[c(1, 10, 12, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$bin_uri), ] 
    # If all rows have data for column...
  } else {
    # If no entries need to be removed, just rename and reorganize x.
    z <- x[c(1, 10, 12, y)]
    # Remove duplicate entries.
    z <- z[!duplicated(z$bin_uri), ]  
  }
  rm(noY)
  return(z)
}
