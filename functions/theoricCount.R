#Find which genes are part of more minimal cut sets and study their functions
theoricCount <- function(MCS,allGenesData,sl_human){
  
   sl_humanDRAFT <- sl_human[,1:4]
   colnames(sl_humanDRAFT) <- c("symbol","entrez_id","symbol","entrez_id")
   ALLSLnames <- rbind(sl_humanDRAFT[,1:2],sl_humanDRAFT[,3:4])
   
   #First delete the NAs of every MCS
   MCSnoNAS <- apply(MCS, 1,function(rowMCS){
     #Function that deletes NAs from a given vector and then creates all the possible pairs of the vector without repetition
     indx <- which((is.na(rowMCS)==TRUE))
     if (length(indx) >0) {
       rowMCS <- rowMCS[-indx] # Delete NA
     } 
     return(rowMCS)
   } )
   
   ## METABOLIC GENES COUNT IN EACH GMCS ----
   AllMCSgenes <- unlist(MCSnoNAS) #Convert list to vector, only 524 genes are part of the gMCS
   MCSgenesCount <- ddply(as.data.frame(AllMCSgenes),.(AllMCSgenes),nrow) #counts how many time a row is repeated 
   geneNames<- allGenesData$symbol[match(MCSgenesCount$AllMCSgenes,allGenesData$entrez_id)]
   MCSgenesCount <- cbind(geneNames,MCSgenesCount)
   colnames(MCSgenesCount)[3] <- "count"
   MCSgenesCount <- MCSgenesCount[order(-MCSgenesCount$count),]
   colnames(MCSgenesCount)<-c("symbol","MCSgene","count")
   #Didn't include the zero values
   
   barplot(MCSgenesCount$count[1:20], las = 2, names.arg = MCSgenesCount$symbol[1:20], col ="yellow", main ="Most frequent metabolic genes in gMCSs", ylab = "frequencies")
   
   TOPenrichMCS <- ORA("TOPenrichMCS",as.vector(MCSgenesCount[1:131,1]),NULL)
   WORSTenrichMCS <- ORA("WORSTenrichMCS",as.vector(MCSgenesCount[(nrow(MCSgenesCount)-130):nrow(MCSgenesCount),1]),NULL)
   
   # weighted set cover
   indxTOPenrichMCS <-  which(is.na(match(TOPenrichMCS$geneSet,WORSTenrichMCS$geneSet)))
   indxWORSTenrichMCS <- which(is.na(match(WORSTenrichMCS$geneSet,TOPenrichMCS$geneSet)))
   totalmatchMCS <- sum((match(TOPenrichMCS$geneSet,WORSTenrichMCS$geneSet)!="NA"),na.rm = TRUE)

   wscTOPenrichMCS <- weighted(TOPenrichMCS[indxTOPenrichMCS,])
   wscWORSTenrichMCS <- weighted(WORSTenrichMCS[indxWORSTenrichMCS,])
   
   topSetsMCS <- TOPenrichMCS[match(wscTOPenrichMCS$topSets,TOPenrichMCS$geneSet),]
   worstSetsMCS <- WORSTenrichMCS[match(wscWORSTenrichMCS$topSets,WORSTenrichMCS$geneSet),]
   
   coverTOPMCS <- length(unique(unlist(sapply(topSetsMCS$userId, strsplit, split=";"))))/131
   coverWORSTMCS <- length(unique(unlist(sapply(worstSetsMCS$userId, strsplit, split=";"))))/131
   
   
    #Obtain the count of the SL gMCS genes
   #MCSgenesCount$count[match(uniqueSLinMCS,MCSgenesCount$AllMCSgenes)]
   #5770 3849 5618  409 9037 9084    5    2    6   58   18 9430  148 5423 2274  340
   
   #Obtain if the SL gMCS genes are in the top 50 of the MCSgenesCount
   #MCSgenesCount$geneNames[match(uniqueSLinMCS,MCSgenesCount[1:50,2])] #NEU1  RRM1  SLC29A1 SLC29A2  CMPK1   CTSA     

   ## SL GENES (all genes that belong to a pair) COUNT IN EACH GMCS ----
   genesSL <- unique(c(pairsSL[,1],pairsSL[,2])) #Make array with all genes that belong to an SL pair, there are 5130 genes
   SLgenesCount <- MCSgenesCount[-which(is.na(match(MCSgenesCount$MCSgene,genesSL))),] #Take the count of the SLs from the MCS genes count
   colnames(SLgenesCount)<-c("symbol","SLgene","count")
   
   # #Add those genes that have a count of 0
   # NOTinSLgenesCount <- as.data.frame(genesSL[which(is.na(match(genesSL,MCSgenesCount$AllMCSgenes)))])
   # NOTinSLgenesCount<-cbind(NOTinSLgenesCount,rep(0,nrow(NOTinSLgenesCount)))
   # colnames(NOTinSLgenesCount) <- c("SLgene","count")
   # NOTgeneNamesSL<- ALLSLnames$symbol[match(NOTinSLgenesCount$SLgene,ALLSLnames$entrez_id)]
   # NOTinSLgenesCount <-cbind(NOTgeneNamesSL,NOTinSLgenesCount)
   # colnames(NOTinSLgenesCount)[1]<-"symbol"
   # 
   # SLgenesCount <- rbind(SLgenesCount,NOTinSLgenesCount)

   
   barplot(SLgenesCount$count[1:20], las = 2, names.arg = SLgenesCount$symbol[1:20], col ="lightblue", main ="Most frequent SL genes in gMCSs", ylab = "frequencies")
   
   TOPenrichSL <- ORA("TOPenrichSL",as.vector(SLgenesCount[1:32,1]),NULL)
   WORSTenrichSL <- ORA("WORSTenrichSL",as.vector(SLgenesCount[(nrow(SLgenesCount)-31):nrow(SLgenesCount),1]), NULL)
   
   # weighted set cover
   indxTOPenrichSL <-  which(is.na(match(TOPenrichSL$geneSet,WORSTenrichSL$geneSet)))
   indxWORSTenrichSL <- which(is.na(match(WORSTenrichSL$geneSet,TOPenrichSL$geneSet)))
   totalmatchSL <- sum((match(TOPenrichSL$geneSet,WORSTenrichSL$geneSet)!="NA"),na.rm = TRUE)
   
   wscTOPenrichSL <- weighted(TOPenrichSL[indxTOPenrichSL,])
   wscWORSTenrichSL <- weighted(WORSTenrichSL[indxWORSTenrichSL,])
   
   topSetsSL <- TOPenrichSL[match(wscTOPenrichSL$topSets,TOPenrichSL$geneSet),]
   worstSetsSL <- WORSTenrichSL[match(wscWORSTenrichSL$topSets,WORSTenrichSL$geneSet),]

   coverTOPSL <- length(unique(unlist(sapply(topSetsSL$userId, strsplit, split=";"))))/32
   coverWORSTSL <- length(unique(unlist(sapply(worstSetsSL$userId, strsplit, split=";"))))/32
   
   # ver cuantos genes coinciden entre los dos enrichments
   idsTOPmcs <- unique(unlist(sapply(topSetsMCS$userId, strsplit, split=";")))
   idsTOPsl <- unique(unlist(sapply(topSetsSL$userId, strsplit, split=";")))
   idsTOPmcs[which((match(idsTOPmcs,idsTOPsl))!="NA")]
   
   
   ## SL GENES COUNT IN EACH SL PAIR ----
   SLgenesCountinSL <- ddply(as.data.frame(c(pairsSL[,1],pairsSL[,2])),.(c(pairsSL[,1],pairsSL[,2])),nrow) #counts how many time a row is repeated 
   colnames(SLgenesCountinSL) <- c("SLgene","count")
   
   geneNamesSL<- ALLSLnames$symbol[match(SLgenesCountinSL$SLgene,ALLSLnames$entrez_id)]
   
   SLgenesCountinSL <- cbind(geneNamesSL,SLgenesCountinSL)
   SLgenesCountinSL <- SLgenesCountinSL[order(-SLgenesCountinSL$count),]
   colnames(SLgenesCountinSL)[1]<-"symbol"
   
   barplot(SLgenesCountinSL$count[1:20], las = 2, names.arg = SLgenesCountinSL$symbol[1:20], col ="pink", main ="Most frequent SL genes in SL pairs", ylab = "frequencies")
   
   TOPenrichSLinSL <- ORA("TOPenrichSLinSL",as.vector(SLgenesCountinSL[1:1283,1]),NULL)
   WORSTenrichSLinSL <- ORA("WORSTenrichSLinSL",as.vector(SLgenesCountinSL[(nrow(SLgenesCountinSL)-2031):nrow(SLgenesCountinSL),1]),NULL)
   
   AllenrichSLinSL <- ORA("AllenrichSLinSL",as.vector(SLgenesCountinSL[,1]),NULL)
   
   # weighted set cover
   indxTOPenrichSLinSL <-  which(is.na(match(TOPenrichSLinSL$geneSet,WORSTenrichSLinSL$geneSet)))
   indxWORSTenrichSLinSL <- which(is.na(match(WORSTenrichSLinSL$geneSet,TOPenrichSLinSL$geneSet)))
   totalmatchSLinSL <- sum((match(TOPenrichSLinSL$geneSet,WORSTenrichSLinSL$geneSet)!="NA"),na.rm = TRUE)
   
   wscTOPenrichSLinSL <- weighted(TOPenrichSLinSL[indxTOPenrichSLinSL,])
   wscWORSTenrichSLinSL <- weighted(WORSTenrichSLinSL[indxWORSTenrichSLinSL,])
   
   topSetsSLinSL <- TOPenrichSLinSL[match(wscTOPenrichSLinSL$topSets,TOPenrichSLinSL$geneSet),]
   worstSetsSLinSL <- WORSTenrichSLinSL[match(wscWORSTenrichSLinSL$topSets,WORSTenrichSLinSL$geneSet),]
  
   coverTOPSLinSL <- length(unique(unlist(sapply(topSetsSLinSL$userId, strsplit, split=";"))))/1283
   coverWORSTSLinSL <- length(unique(unlist(sapply(worstSetsSLinSL$userId, strsplit, split=";"))))/2032
   

}
