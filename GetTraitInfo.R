GetTraitInfo <- function(x, y, type = c("continuous", "discrete")) {
  # Purpose: Provides descriptive information for data of a certain trait.
  # x = Dataframe of species and trait information.
  # y = Trait (column) name of interest.
  # type = Class of trait.
  
  type <- match.arg(type)
  
  # If the trait is continuous...
  if(type == "continuous") {
    # Convert to dataframe.
    x <- as.data.frame(x)
    # How many rows are there?
    print(paste0("Number of observations: ", nrow(x)))
    # What is the range of the data?
    print(paste0("Range of data: ", range(x[, y])[1], " to ", range(x[, y])[2]))
    # What is the mean?
    print(paste0("Mean: ", mean(x[, y])))
    # What is the median?
    print(paste0("Median: ", median(x[, y])))
    # Plot of histogram.
    hist(x[, y], main = "", xlab = paste(y))
    paste("Histogram plotted!")
    
    # If the trait is discrete (categorical)...
  } else if(type == "discrete") {
    # Convert to dataframe.
    x <- as.data.frame(x)
    # How many rows are there?
    print(paste0("Number of observations: ", nrow(x)))
    # Table of trait information:
    print("Number of observations per category: ")
    table(x[, y])
  }
}
