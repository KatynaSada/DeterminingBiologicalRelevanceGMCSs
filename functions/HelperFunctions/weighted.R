weighted <- function(enrichedSig){
  idsInSet <- sapply(enrichedSig$overlapId, strsplit, split=";")
  names(idsInSet) <- enrichedSig$geneSet
  minusLogP <- -log(enrichedSig$pValue)
  minusLogP[minusLogP == Inf] <- -log(.Machine$double.eps)
  setCoverNum <- 5
  nThreads <- 1
  wscRes <- weightedSetCover(idsInSet, 1 / minusLogP, setCoverNum, nThreads)
  
  return(wscRes)
}
