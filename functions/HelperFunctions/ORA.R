#' Over Representation Analysis
#' @description With a given list of genes it returns its ORA
#' @param projectName The name of the file that is goin to have all the results. 
#' @param genes Array of genes to be analyzed 
#' @param referenceGenes Array of background genes
#'
#' @return A list with the results of the enrichment analysis
ORA<-function(projectName,genes,referenceGenes){
  
  #fileName<-"./input_data/MCSgenesTOP.txt"
  #genes<-MCSgenesCount[1:20,1]
  
  outputDirectory <- "./ORA"
  
  enrichResult <-WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                             enrichDatabase="geneontology_Biological_Process", enrichDatabaseFile=NULL, enrichDatabaseType=NULL,
                             enrichDatabaseDescriptionFile=NULL, interestGeneFile=NULL, interestGene=genes,
                             interestGeneType="genesymbol", collapseMethod="mean", 
                             referenceGeneFile=NULL, referenceGene=referenceGenes,
                             referenceGeneType="genesymbol", referenceSet="genome", 
                             minNum=2, maxNum=1000, fdrMethod="BH", sigMethod="fdr",
                             fdrThr=0.05, topThr=10, reportNum=20, perNum=1000, isOutput=TRUE, outputDirectory=outputDirectory,
                             projectName=projectName, dagColor="continuous",hostName="http://www.webgestalt.org/")
  return(enrichResult)
}