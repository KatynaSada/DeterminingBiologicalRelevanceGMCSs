combinations <- function(rowMCS){
  #Function that deletes NAs from a given vector and then creates all the possible pairs of the vector without repetition
  indx <- which((is.na(rowMCS)==TRUE))
  if (length(indx) >0) {
    rowMCS <- rowMCS[-indx] # Delete NA
    } 
  return(t(combn(rowMCS, 2, FUN = NULL, simplify = TRUE)))
}