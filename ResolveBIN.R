# Function: ResolveBIN #
# Purpose:  Resolves BINs with taxonomic conflicts by either removing records or
#           the BIN itself from a dataframe.
ResolveBIN <- function(x, y, method = c("bin_uri","recordID")){
  # x = Either the BIN or recordID to be removed.
  # y = Dataframe to apply function to.
  # method = Specify whether the item to be removed is a BIN or a recordID.
  method <- match.arg(method)
  if(method == "bin_uri") {
    rmThisBIN <- which(y$bin_uri == x)
    resolved <- y[-rmThisBIN, ]
  } else if(method == "recordID") {
    rmThisRecord <- which(y$recordID == x)
    resolved <- y[-rmThisRecord, ]
  }
  return(resolved)
}