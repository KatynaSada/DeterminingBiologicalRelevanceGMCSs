#' Count duplicates number
#' @description Function that has as an input two matrices of different/equal amount of pairs and returns the amount of pairs that match
#' @param A Matrix of pairs (has 2 columns)
#' @param B Matrix of pairs (has 2 columns), it can have different amount of rows from those of B
#'
#' @return Amount of pairs that match
countDuplicatesNumber <- function(A,B){
  na <- nrow(A)
  nb <- nrow(B)
  AB <- rbind(A, B)
  ab <- duplicated(AB)
  return(sum(ab)) # Return amount of pairs that match 
}