cosmic<- function(metabollicSLin,tableWith0){
   #COSMIC
   #ver el pvalue de cada pareja de gMCS
   #luego hacer lo mismo con todas las parejas que no estan en los gMCS
   #hacer un boxplot de los pvalores
   #luego en el futuro veremos que test usar 
   
   #Load pvalues of being mutually exclusive
   load("./input_data/pvals.RData") # View(as.matrix(pvals))
   #Obtain pvalues of pairs being mutually exclusive for the 27 SL pairs that are in the MCS
   #The diagonal contains the wanted values (por ahora hacer la matriz con todas las combinaciones entre esas 27 parejas  y luego coger la diagonal es lo mas rapido que he logrado hacer)
   pvalsCosmic27 <- (diag(as.matrix(pvals[match(metabollicSLin[,1],pvals@Dimnames[[1]]),match(metabollicSLin[,2],pvals@Dimnames[[2]])])))
   pvalsCosmic27Table <- cbind(metabollicSLin,exp(pvalsCosmic27))
   
   geneSymbols1<- sl_human$HGNC_1[match(pvalsCosmic27Table$gene1,sl_human$Entrez_1)]
   geneSymbols2<- sl_human$HGNC_2[match(pvalsCosmic27Table$gene2,sl_human$Entrez_2)]
   
   
   pvalsCosmic27Table <-cbind(geneSymbols1,geneSymbols2,pvalsCosmic27Table)
   pvalsCosmic27Table <- pvalsCosmic27Table[order(-pvalsCosmic27Table$`exp(pvalsCosmic27)`),]
   
   
   #Obtain pvalues of pairs being mutually exclusive for all the other SL pairs 
   pvalsCosmicAllMinus27<- log(diag(as.matrix(pvals[match(tableWith0[dim(metabollicSLin)[1]+1:dim(tableWith0)[1] ,1],pvals@Dimnames[[1]]),match(tableWith0[28:2183,2],pvals@Dimnames[[2]])])))
   pvalsCosmicAllMinus27Table <- cbind(tableWith0[28:2183,],exp(pvalsCosmicAllMinus27))
   #Delete NA values
   pvalsCosmicAllMinus27Table <- pvalsCosmicAllMinus27Table[-which(is.na(pvalsCosmicAllMinus27)),]
   
   #Plot them!
   means1 <- c(mean(pvalsCosmic27),log(mean(pvalsCosmicAllMinus27Table$`exp(pvalsCosmicAllMinus27)`)))
   boxplot(pvalsCosmic27,pvalsCosmicAllMinus27,col = c("gray","yellow"),names = c("SL", "Not SL"),ylab = "Logarithm of p-values")
   points(1:2, means1, col = "red")
   text(1:3, means1 -0.5, labels = round(means1,digits=4))
   
   means2 <- c(mean(exp(pvalsCosmic27)),(mean(pvalsCosmicAllMinus27Table$`exp(pvalsCosmicAllMinus27)`)))
   boxplot(exp(pvalsCosmic27),exp(pvalsCosmicAllMinus27),col = c("gray","yellow"),names = c("SL", "Not SL"),ylab = "P-values")
   points(1:2, means2, col = "red")
   text(1:3, means2 + 0.04, labels = round(means2,digits=4))
   
   hist(exp(pvalsCosmic27),20) 
   hist(exp(pvalsCosmicAllMinus27),1000) 
   
   #Percentage of values that have a significant p-value
   sum(exp(pvalsCosmic27)<0.5)/length(pvalsCosmic27)
   sum(pvalsCosmicAllMinus27Table$`exp(pvalsCosmicAllMinus27)`<0.5)/nrow(pvalsCosmicAllMinus27Table)
   
   # pvalsCosmic27<- apply(metabollicSLin,1,function(pareja){
   #  valor <- pvals[match(pareja[1],pvals@Dimnames[[1]]),match(pareja[2],pvals@Dimnames[[2]])]
   #  return(valor)})
   
   #COSMIC SCORES
   
   #Obtain log of pvalues of all MCS pairs
   pvalsCosmicAll<- diag(as.matrix(pvals[match(countPairs[,1],pvals@Dimnames[[1]]),match(countPairs[,2],pvals@Dimnames[[2]])]))
   pvalsCosmicAll[which(pvalsCosmicAll==0)]<-0.0000000001
   # NA values indicate one or both of the genes are not in the matrix
   
   # Add gene names 
   pvalsCosmicAll <- cbind(countPairs[,1:2],pvalsCosmicAll) 
   # which(pvalsCosmicAll==0)
   
   #Find the scores of every MCS 
   # THIS PART TAKES TIME
   scores <- lapply(numberedPairsMCS,function(AgMCSpairs){
      AgMCSpairs <-t(apply(AgMCSpairs,1,sort)) #sort pairs
      positionsPvals <- row.match(as.data.frame(AgMCSpairs), pvalsCosmicAll[,1:2]) # Find the position of the pvals of those pairs
      score<- -2*sum(log(pvalsCosmicAll[positionsPvals,3]),na.rm = TRUE) # Obtain pvals of those pairs, then their ln and finally add them
      return(score)}
      )
   # scores is a list with the sum of pvals of the pairs in each MCS
   
   
   scores <- unlist(scores)
   numberOfPairs <- unlist(lapply(numberedPairsMCS, nrow))
   chiScores <- pchisq(scores,2*numberOfPairs,lower.tail = FALSE)
   chiTable <- as.data.frame(cbind(numberOfPairs,scores,chiScores))
   chiTable <- chiTable[order(chiTable$chiScores),]
   
   #scores <- scores[order(-scores)]
   
   
   #percentageTable se crea en SLinMCS.R
   percentageANDscore <- as.data.frame(cbind(percentageTable,chiTable[match(rownames(percentageTable),rownames(chiTable)),]))
   #percentageANDscore <- as.data.frame(cbind(percentageANDscore,chiScores[match(rownames(percentageANDscore),names(chiScores))]))
   #colnames(percentageANDscore) <- c("Percentage","Score 1","Cosmic score")
   #percentageANDscore  <- percentageANDscore[order(percentageANDscore$`Cosmic score`),]

   #all <- as.data.frame(cbind(percentageANDscore,percentageTableANDhypergeometric2$individualProbabilities2[match(rownames(percentageANDscore),rownames(percentageTableANDhypergeometric2))]))
   

   summary(lm(percentageANDscore$Percentage~percentageANDscore$`Cosmic score`))
   plot(percentageANDscore$Percentage,percentageANDscore$`Cosmic score`)
 
   # pval is 2(-16), if the percentage increases, then the score increases (coefficient=0.02232)
   
   #ESTO PASA A COSMIC
   pvalsCosmicMetabolic <- diag(as.matrix(pvals[match(pairsAllGenes[1:20000,1],pvals@Dimnames[[1]]),match(pairsAllGenes[1:20000,2],pvals@Dimnames[[2]])]))
   pvalsCosmicMetabolic <- cbind(countPairs[,1:2],pvalsCosmicMetabolic)
   
   boxplot(pvalsCosmicMetabolic,pvalsCosmic27)
   }
