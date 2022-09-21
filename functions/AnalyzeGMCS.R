
number <- "2816"

rowMCS <- MCS[number,]

indx <- which((is.na(rowMCS)==TRUE))
if (length(indx) >0) {
  rowMCS <- rowMCS[-indx] # Delete NA
} 

allGenesData[match(rowMCS,allGenesData$entrez_id),2]

cbind(as.vector(allGenesData$symbol[match(SLperMCSpairs[[number]][,1],allGenesData$entrez_id)]),as.vector(allGenesData$symbol[match(SLperMCSpairs[[number]][,2],allGenesData$entrez_id)]))

allGenesData[match(unique(c(SLperMCSpairs[[number]][,1],SLperMCSpairs[[number]][,2])),allGenesData$entrez_id),2]


