KStest <- function(allGenes){
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
  indxSLandMCS <- row.match(metabolicPairsSL,countPairs[,1:2])
  metabollicSLin<- metabolicPairsSL[which(!is.na(indxSLandMCS)),]
  tableSLinMCS <- countPairs[indxSLandMCS[which(!is.na(indxSLandMCS))],] #mismas 27
  colnames(tableSLinMCS) <- c("gene1","gene2","count") 
  dim(tableSLinMCS)[1]
  
  #Add to the table the SL pairs that appear 0 times
  idxtable0 <- duplicated(rbind(metabollicSLin, metabolicPairsSL))[(nrow(metabollicSLin)+1):(nrow(metabollicSLin)+nrow(metabolicPairsSL))]
  tableWith0 <- metabolicPairsSL[-which(idxtable0),]
  tableWith0 <- cbind(tableWith0, rep(0,dim(tableWith0)[1]))
  colnames(tableWith0) <- c("gene1","gene2","count") 
  tableWith0 <- rbind(tableSLinMCS,tableWith0)
  plot(density(tableSLinMCS[,3]))
  plot(density(tableWith0[,3]))
  
  
  ## Random pairs  ------------------------------------------------------------
  #GOAL: hacer mil tiradas aleatorias y sacar mil p-values con ks-test
  #hacerlo a 1 cola, tienen que ser diferentes y uno tiene que quedar mas a la izquierda
  #nosotros mas numeros en la no aleatoria y mas ceros en la aleatoria, mas a la derecha la nuestra

  pairsAllGenes <- t(combn(allGenes, 2, FUN = NULL, simplify = TRUE)) # Obtain all the possible pairs of metabolic genes
  #pairsAllGenes<- t(apply(t(apply(pairsAllGenes, 1, as.numeric)),1,sort)) # Sort in ascending order
  colnames(pairsAllGenes) <- c("gene1","gene2")
  
  #Create list of 1000 matrixes of 2183 pairs of metabolic genes
  listFalseSL <- list()
  for (i in 1:1000) {
    
    indxRandom <- sample(1:dim(pairsAllGenes)[1],dim(metabolicPairsSL)[1],replace = F) #Choose random indexes without replacement
    
    falseSL <- pairsAllGenes[indxRandom,] #Obtain the random pairs
    falseSL <- t(apply(t(apply(falseSL, 1, as.numeric)),1,sort))
    # if( typeof(falseSL) == "list"){ #Sometimes there´s an error an falseSL becomes a list
    #   falseSL <- do.call(rbind,falseSL)
    # }
    #
    colnames(falseSL) <- c("gene1","gene2")
    listFalseSL[[i]]<-falseSL
    rm(falseSL)
  }
  
  pvaluesKS<- lapply(listFalseSL,function(falseSL){
    #Find how many times each of the pairs of falseSL appears in a MCS
    indxFalseSLandMCS <- row.match(as.data.frame(falseSL),countPairs[,1:2]) #Match pairs
    indxFalseSLandMCS<-which(!is.na(indxFalseSLandMCS)) #Indexes of pairs that match
    
    if (length(indxFalseSLandMCS)==0){
      falseTableWith0 <- cbind(falseSL, rep(0,dim(falseSL)[1]))
      colnames(falseTableWith0) <- c("gene1","gene2","count")
    }else{
      
      falseSLin<- falseSL[indxFalseSLandMCS,] #SLs that are in the MCS
      falsetableSLinMCS <- countPairs[indxFalseSLandMCS[which(!is.na(indxFalseSLandMCS))],] #Includes count
      colnames(falsetableSLinMCS) <- c("gene1","gene2","count")
      dim(falsetableSLinMCS)[1]
      
      #Add to the table the SL false pairs that appear 0 times
      falseTableWith0 <- falseSL[-indxFalseSLandMCS,] #Obtain the pairs that don't have a match
      falseTableWith0 <- cbind(falseTableWith0, rep(0,dim(falseTableWith0)[1]))
      colnames(falseTableWith0) <- c("gene1","gene2","count")
      falseTableWith0 <- rbind(falsetableSLinMCS,falseTableWith0)}
    
    #Our distribuition should be more to the right than the random one, the mean should be larger
    #pvalue<- ks.test(tableWith0[,3],falseTableWith0[,3],alternative = "less")
    pvalue<- ks.test(tableWith0[,3],falseTableWith0[,3],alternative = "greater") #null hypothesis: x is not greater than y
    return(pvalue[["p.value"]])
  })
  
  pvaluesKS <- do.call(rbind,pvaluesKS) #convert from list to vector
  hist(pvaluesKS,1000)
  
  
  plot(density(tableWith0[,3]))
  lines(density(falseTableWith0[,3]),col="red")

#Algo muy importante es que esas 27 salen en MUCHOS MCS
}