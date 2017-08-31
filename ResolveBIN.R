# Function: ResolveBIN #
# Purpose:  Resolves BINs with taxonomic conflicts by either removing records or
#           the BIN itself from a dataframe.
ResolveBIN <- function(x, y, method = c("bin_uri","recordID")){
  # x = Dataframe to apply function to.
  # y = Either the BIN or recordID to be removed.
  # method = Specify whether the item to be removed is a BIN or a recordID.
  
  # Is it a BIN or record?
  method <- match.arg(method)
  # If it's a BIN...
  if(method == "bin_uri") {
    # Identify the indices.
    rmThisBIN <- which(X$bin_uri == Y)
    # Remove from the dataframe.
    resolved <- x[-rmThisBIN, ]
  } else if(method == "recordID") {
    # Identify the indices.
    rmThisRecord <- which(x$recordID == y)
    # Remove from the dataframe.
    resolved <- x[-rmThisRecord, ]
  }
  # Return the resolved dataframe.
  return(resolved)
}