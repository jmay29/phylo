RemoveSequences <- function(x, y) {
  # x = Dataframe of sequence information.
  # y = BIN URIs of sequences to be removed.
  if (length(y) == 0) {
    print ("There are no sequences to remove!")
  }
  else if (length(y) > 0) {
    x <- x[!x$species_name %in% y]
  }
  return(x)
}

