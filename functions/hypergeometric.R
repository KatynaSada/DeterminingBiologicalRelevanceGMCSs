#' Hypergeometric
#' @description 
#' @param allGenes Array of all metabollic genes
#' @param MCS Data frame containing all genetic minimal cut sets, each row is a different gMCS
#' @param pairsSL All the known SL pairs 
#'
#' @return Probability of an SL pair of being a gMCS pair
hypergeometric <- function(allGenes,MCS,pairsSL){
  #Make all possible pairs of genes of each MCS
  
  numberedPairsMCS <- apply(MCS, 1, combinations)
  pairsMCS <- list_df2df(numberedPairsMCS)[,2:3] # Convert list to dataframe
  #pairsMCS<- t(apply(t(apply(pairsMCS, 1, as.numeric)),1,sort,method="quick")) # Sort pairs (required for next step), aprox 10mins
  #pairsMCS <- distinct(as.data.frame(pairsMCS)) # Delete repeated pairs 
  colnames(pairsMCS) <- c("gene1","gene2")
  
  gene1 <- pairsMCS[,1]
  gene2 <- pairsMCS[,2]
  pairsMCS <- ddply(pairsMCS,.(gene1,gene2),nrow) #counts how many time a row is repeated, so it deletes duplicates
  pairsMCS[,1:2] <- t(apply(t(apply(pairsMCS[,1:2], 1, as.numeric)),1,sort)) #sort rows 
  
  m <- dim(pairsMCS)[1] #Number of white balls in the urn, amount of different MCS pairs
  
  # Count how many MCS pairs are SL
  duplicates <- countDuplicates(pairsMCS[,1:2],pairsSL) #27
   
  x <-dim(duplicates)[1] #Number of white balls drawn from the urn, how many MCS pairs are SLs
  
  # Create all the possible pairs from the metabollic genes
  pairsAllGenes <- t(combn(allGenes, 2, FUN = NULL, simplify = TRUE)) # Obtain all the possible pairs of metabolic genes
  colnames(pairsAllGenes) <- c("gene1","gene2") 
  N <- dim(pairsAllGenes)[1]
  n <- N-m #Number of black balls in the urn, all the other pairs that are not MCS pairs
  
  #Find those SL pairs in which both of their genes are metabollic 
  indxMetabolicPairsSL <- which(complete.cases(sapply(pairsSL, match,table=allGenes))) 
  #sapply + match = NA if value is not in allGenes or position in allGenes
  #complete.cases: indicates TRUE if there are no missing values in the row (NA) 
  #which: position of rows which both have matches in allGenes
  metabolicPairsSL<-pairsSL[indxMetabolicPairsSL,] #Keep those pairs that have both of their genes in allGenes
  colnames(metabolicPairsSL) <- c("gene1","gene2")
  k <- dim(metabolicPairsSL)[1] # Number of balls drawn from the urn, all SL pairs that have metabollic genes
  
  probability <- phyper(x-0.1,m,n,k,lower.tail = FALSE, log.p = FALSE) # lower.tail = FALSE, probabilities are P[X > x], substract 0.1 to take into account x (greater or equal to)
  # We use the lower tail so that we are able to see a significant value
  
  
  ## Matrices  ------------------------------------------------------------
  # Create matrix of only MCS of pairs that match 
  # genesMCS <- unique(c(pairsMCS$gene1,pairsMCS$gene2)) #Dimension of matrix is equal to all genes that are part of a MCS
  # ind1 <- match(duplicates$gene1,genesMCS) 
  # ind2 <- match(duplicates$gene2,genesMCS) 
  # matrixMatch <- sparseMatrix(i=ind1, j=ind2,x = 1,dims = c(length(genesMCS), length(genesMCS)))
  # rownames(matrixMatch) <- colnames(matrixMatch) <- genesMCS
  # matrixMatch <- matrixMatch + t(matrixMatch) ## Make it symmetric
  # image(matrixMatch)
  
  # Create matrix of all genes of pairs that match 
  #genesAll <- unique(c(pairsMCS$gene1,pairsMCS$gene2)) #Dimension of matrix is equal to all genes that are part of a MCS
  # allGenes <- as.numeric(allGenes)
  # ind12 <- match(duplicates$gene1,allGenes) 
  # ind22 <- match(duplicates$gene2,allGenes) 
  # matrixMatch <- sparseMatrix(i=ind12, j=ind22,x = 1,dims = c(length(allGenes), length(allGenes)))
  # rownames(matrixMatch) <- colnames(matrixMatch) <- allGenes
  # matrixMatch <- matrixMatch + t(matrixMatch) ## Make it symmetric
  # image(matrixMatch)
  
  return(probability)
}  