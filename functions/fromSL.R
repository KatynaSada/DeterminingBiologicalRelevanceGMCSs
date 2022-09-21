fromSL <- function(allGenes,MCS,pairsSL){
  #Find those SL pairs in which both of their genes are metabollic 
  indxMetabolicPairsSL <- which(complete.cases(sapply(pairsSL, match,table=allGenes))) 
  #sapply + match = NA if value is not in allGenes or position in allGenes
  #complete.cases: indicates TRUE if there are no missing values in the row (NA) 
  #which: position of rows which both have matches in allGenes
  metabolicPairsSL<-pairsSL[indxMetabolicPairsSL,] #Keep those pairs that have both of their genes in allGenes
  row.names(metabolicPairsSL) <- 1:nrow(metabolicPairsSL)
  colnames(metabolicPairsSL) <- c("gene1","gene2")
  
  #Make all MCS pairs
  numberedPairsMCS <- apply(MCS, 1, combinations)
  pairsMCS2 <- list_df2df(numberedPairsMCS) # Convert list to dataframe
  colnames(pairsMCS2)[2:3] <- c("gene1","gene2")
  pairsMCS2 <- pairsMCS2[,2:3] 

  #Count in how many MCS each pair of all metabolic genes is present (all of those pairs which are not part of pairsMCS are not part of the table, logicamente)
  gene1 <- pairsMCS2[,1]
  gene2 <- pairsMCS2[,2]
  
  countPairs <- ddply(pairsMCS2,.(gene1,gene2),nrow) #counts how many time a row is repeated 
  countPairs[,1:2] <- t(apply(t(apply(countPairs[,1:2], 1, as.numeric)),1,sort)) #countpairs es como pairsMCS en hypergeometric
  
  #Find which SL metabollic pairs are in countPairs
  indxSLandMCS <- row.match(countPairs[,1:2],metabolicPairsSL)
  metabollicSLin<- countPairs[which(!is.na(indxSLandMCS)),][,c(1,2)] 
  
  geneSymbols1<- sl_human$HGNC_1[match(metabollicSLin$gene1,sl_human$Entrez_1)]
  geneSymbols2<- sl_human$HGNC_2[match(metabollicSLin$gene2,sl_human$Entrez_2)]
  metabollicSLinNames <- cbind(metabollicSLin,geneSymbols1,geneSymbols2)
  
  #View all SynthLethDB Data
  pairsSAll1 <- sl_humanOutput[row.match(metabollicSLin,cbind(sl_humanOutput$Entrez_1,sl_humanOutput$Entrez_2)),]
  pairsSAll2 <- sl_humanOutput[row.match(metabollicSLin,cbind(sl_humanOutput$Entrez_2,sl_humanOutput$Entrez_1)),]
  pairsSAll1[-which(is.na(pairsSAll2[,1])),] <- pairsSAll2[-which(is.na(pairsSAll2[,1])),]
  
  }

