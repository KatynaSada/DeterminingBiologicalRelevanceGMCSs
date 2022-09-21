#' Count Duplicated Pairs
#' @description Function that has as an input two matrices of different/equal amount of pairs and returns pairs that match
#' @param A Matrix of pairs (has 2 columns)
#' @param B Matrix of pairs (has 2 columns), it can have different amount of rows from those of B
#'
#' @return duplicatePairs: matrix of pairs that are in both A and B
countDuplicates <- function(A,B){
  na <- nrow(A)
  nb <- nrow(B)
  AB <- rbind(A, B)
  ab <- duplicated(AB)
  duplicatePairs <- AB[which(ab),] # Return pairs that match
  return(duplicatePairs)
}