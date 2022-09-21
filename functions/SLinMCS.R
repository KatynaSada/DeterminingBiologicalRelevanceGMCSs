SLinMCS <- function(allGenes,MCS,pairsSL){
  # How many MCS have a SL (or more than 1)?? Find percentage of SL in each gMCS 
  # HAY COSAS QUE SE REPITEN EN fromSL y en KStest
  
  # de aqui
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
  countPairs[,1:2] <- t(apply(t(apply(countPairs[,1:2], 1, as.numeric)),1,sort))
  
  #Find which SL metabollic pairs are in countPairs
  indxSLandMCS <- row.match(countPairs[,1:2],metabolicPairsSL)
  metabollicSLin<- countPairs[which(!is.na(indxSLandMCS)),][,c(1,2)] 
  # a aquiiiiii
  
  #Make a list of the 27 genes in ascending and descending order (it takes a lot of time to order all the MCS pairs)
  bothOrdersSL <- rbind(metabollicSLin,metabollicSLin)
  #Order the last 27
  bothOrdersSL[((dim(bothOrdersSL)[1]/2)+1):dim(bothOrdersSL)[1],]<-t(apply( bothOrdersSL[((dim(bothOrdersSL)[1]/2)+1):dim(bothOrdersSL)[1],], 1, sort, decreasing=T)) 
  
  #indxMCS<-duplicated(rbind(bothOrdersSL,numberedPairsMCS[,2:3]))[(nrow(bothOrdersSL)+1):(nrow(bothOrdersSL)+nrow(numberedPairsMCS))]
  #numberedPairsMCS[which(indxMCS),]
  
  #Set column names to all the elements of the list
  numberedPairsMCS<- lapply(numberedPairsMCS,function(mat){
    colnames(mat) <- c("gene1","gene2")
    return(mat)})
  
  #Find how much SLs are in each gMCS
  SLperMCS<-lapply(numberedPairsMCS, countDuplicatesNumber,A=bothOrdersSL)
  SLperMCS<- t(as.data.frame(SLperMCS))
  
  #Verify...
  #sum(SLperMCS)
  #sum(tableWith0[,3]) #In KStest
  
  #Find which SL pairs are in each gMCS
  SLperMCSpairs<-lapply(numberedPairsMCS, countDuplicates,A=bothOrdersSL)
  # SLperMCSpairs[["38373"]]
  
  # Percentage of SLs in each gMCS
  amountPairs <- do.call(rbind,lapply(numberedPairsMCS, dim))[,1] #Amount of pairs in each gMCS
  percentageTable <- cbind(SLperMCS,amountPairs,SLperMCS/amountPairs)
  colnames(percentageTable)<-c("Total SL","Total Pairs","Percentage")
  rownames(percentageTable) <- str_remove(rownames(percentageTable),"X")
  
  percentageTable <- as.data.frame(percentageTable[order(-percentageTable[,3]),]) #Order percentage
  
  # Percentage of gMCS that each SL belong to 
  percentageTableSLs <- countPairs[which(!is.na(indxSLandMCS)),]
  percentageTableSLs <-cbind(percentageTableSLs,percentageTableSLs[,3]/nrow(MCS)) #Percentage is calculated dividing the amount of gMCSs a SL belongs to into the total amount of gMCSs (41541)
  colnames(percentageTableSLs)<-c("gene1","gene2","Total gMCSs","Percentage")
  
  percentageTableSLs <- percentageTableSLs[order(-percentageTableSLs[,3]),] #Order percentage
 
  geneSymbols1<- sl_human$HGNC_1[match(percentageTableSLs$gene1,sl_human$Entrez_1)]
  geneSymbols2<- sl_human$HGNC_2[match(percentageTableSLs$gene2,sl_human$Entrez_2)]
  
  percentageTableSLs <-cbind(geneSymbols1,geneSymbols2,percentageTableSLs)
  
  
  # # Obtain the HYPERGEOMETRIC of each gMCS
  # x <- percentageTable$`Total SL` # Amount of SL in each gMCS
  # k <- percentageTable$`Total Pairs` # Total amount of pairs in each gMCS
  # m <- dim(countPairs)[1] #  All gMCS pairs, 20.796
  # 
  # # Create all the possible pairs from the metabollic genes
  # pairsAllGenes <- t(combn(allGenes, 2, FUN = NULL, simplify = TRUE)) # Obtain all the possible pairs of metabolic genes
  # colnames(pairsAllGenes) <- c("gene1","gene2") 
  # N <- dim(pairsAllGenes)[1] # Total number of metabolic pairs
  # n <- N - m
  # individualProbabilities <- phyper(x-0.1,m,n,k,lower.tail = FALSE, log.p = FALSE) # lower.tail = FALSE, probabilities are P[X > x], substract 0.1 to take into account x (greater or equal to)
  # percentageTableANDhypergeometric <- cbind(percentageTable,individualProbabilities)
  # percentageTableANDhypergeometric <- percentageTableANDhypergeometric[order(percentageTableANDhypergeometric$individualProbabilities),] #Order probability
  # 
  # cor(percentageTableANDhypergeometric$individualProbabilities,percentageTableANDhypergeometric$Percentage)
  # 
  #CORREGIDOOOO
  
  # Obtain the HYPERGEOMETRIC of each gMCS
  x <- percentageTable$`Total SL` # Amount of SL in each gMCS
  k <- 2183
  m <- percentageTable$`Total Pairs` # Total amount of pairs in each gMCS

  # Create all the possible pairs from the metabollic genes
  pairsAllGenes <- t(combn(allGenes, 2, FUN = NULL, simplify = TRUE)) # Obtain all the possible pairs of metabolic genes
  colnames(pairsAllGenes) <- c("gene1","gene2") 
  N <- dim(pairsAllGenes)[1] # Total number of metabolic pairs
  n <- N - m
  individualProbabilities2 <- phyper(x-0.1,m,n,k,lower.tail = FALSE, log.p = FALSE) # lower.tail = FALSE, probabilities are P[X > x], substract 0.1 to take into account x (greater or equal to)
  percentageTableANDhypergeometric2 <- cbind(percentageTable,individualProbabilities2)
  percentageTableANDhypergeometric2 <- percentageTableANDhypergeometric2[order(percentageTableANDhypergeometric2$individualProbabilities),] #Order probability
  
  
}
  