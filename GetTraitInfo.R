GetTraitInfo <- function(x) {
  # Purpose: Provides descriptive information for data of a certain trait.
  # x = Vector containing trait of interest.

  # Get name of trait.
  traitName <- substitute(x)
  print(traitName)
  
  # Get class of trait.
  traitType <- class(x)
  
  # If the trait is continuous...
  if(traitType == "numeric") {

    # How many rows are there?
    print(paste0("Number of observations: ", sum(!is.na(x))))
    # What is the range of the data?
    print(paste0("Range of data: ", range(x, na.rm = T)[1], " to ", range(x, na.rm = T)[2]))
    # What is the mean?
    print(paste0("Mean: ", mean(x, na.rm = T)))
    # What is the median?
    print(paste0("Median: ", median(x, na.rm = T)))
    # Plot of histogram.
    hist(x, main = "", xlab = paste(x))
    paste("Histogram plotted!")
    
    # If the trait is discrete (categorical)...
  } else if(traitType == "character" | traitType == "factor" | traitType == "integer") {
  
    # How many rows are there?
    print(paste0("Number of observations: ", sum(!is.na(x))))
    # Table of trait information:
    print("Number of observations per category: ")
    table(x)
  }
}
