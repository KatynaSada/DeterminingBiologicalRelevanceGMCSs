investigate<- function(metabollicSLin,sl_human,allGenesData){

  #How many genes are there in those 27 pairs?
  uniqueSLinMCS <- unique(c(metabollicSLin[,1],metabollicSLin[,2])) #16 genes!!!
  length(uniqueSLinMCS)
  
  infoGenes<-allGenesData[match(uniqueSLinMCS,allGenesData[,16]),] # Information we have of those genes
  
  library(biomaRt)
  #https://www.ensembl.org/biomart/martview/12788329192943d205efa9813102a21b
  #BIOMART

  #Download data from ensembl
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  data16genes = getBM(attributes = c('entrezgene_id', 'external_gene_name', 'gene_biotype','name_1006','definition_1006','hpa_accession'), 
                filters = 'entrezgene_id', 
                values = uniqueSLinMCS, 
                mart = ensembl,
                uniqueRows=TRUE)

  
  #Gene Ontology
  enrichResult16 <- ORA("enrichResult16_0.05",as.vector(allGenesData$symbol[match(uniqueSLinMCS,allGenesData$entrez_id)]),as.vector(allGenesData$symbol))
  

  
  #EXCEL 
  #STRING 
  library(STRINGdb)
  string_db <- STRINGdb$new( version="10", species=9606, score_threshold=0, input_directory="" ) #9606 for human
  
  # CMPK1 <- string_db$mp(3177)
  # otro <- string_db$mp(6240)
  
  #Obtain ids of STRING database for the 27 pairs (and include gene symbol in table)
  string_ids <- string_db$map(infoGenes, "symbol", removeUnmappedRows = FALSE ) #we map the gene names to the STRING database identifiers using the ”map” method
  string_idsPairs<- cbind(string_ids$STRING_id[match(metabollicSLin[,1],string_ids$entrez_id)],string_ids$STRING_id[match(metabollicSLin[,2],string_ids$entrez_id)])
  symbols <-cbind(string_ids$symbol[match(metabollicSLin[,1],string_ids$entrez_id)],string_ids$symbol[match(metabollicSLin[,2],string_ids$entrez_id)])
  
  string_idsPairs <-cbind(metabollicSLin,string_idsPairs,symbols)
  string_idsPairs <- as.data.frame(string_idsPairs)
  row.names(string_idsPairs) <- 1:nrow(string_idsPairs)
  colnames(string_idsPairs) <- c("gene1","gene2","stringID1","stringID2","symbol1","symbol2")
  
  #string_id<-as.character(unlist(string_idsPairs[,c(4,5)][1,]))
  
  #string_id<-string_idsPairs[,c(3,4)][1,c(1,2)]
  #string_id<-unlist(string_id)
  
  #Obtain interactions, articles and STRING link (to obtain more info)
  allPairsInfo<- apply(string_idsPairs[,3:4],1,function(string_id){
    
    name <- paste(string_ids$symbol[match(string_id[1],string_ids$STRING_id)],string_ids$symbol[match(string_id[2],string_ids$STRING_id)])
    interactions<-string_db$get_interactions(as.character(unlist(string_id)))
    articles<-string_db$get_pubmed_interaction(as.character(string_id[[1]]),as.character(string_id[[2]])) #retrieve the pubmed identifiers of the articles that contain the name of both the proteins (if any)

    articlesNames <- lapply(articles, function(articleNum){
      
    if (substr(articleNum,1,4)=="PMID"){
    thepage <- readLines(paste("https://pubmed.ncbi.nlm.nih.gov/", gsub("PMID:", "", articleNum),"/",sep = ""))
    } else if (substr(articleNum,1,4)=="OMID") {
    thepage <- readLines(paste("https://www.omim.org/entry/", gsub("OMID:", "", articleNum),sep = ""))
    } 
      
    articleName<-thepage[25]
    articleName<-gsub("    <title>","",articleName)
    articleName<-gsub(" - PubMed</title>","",articleName)
    return(articleName)
    })
    
    articlesInfo <- cbind(articles,unlist(articlesNames))
      
    link<- string_db$get_link(string_id)

     info <- list(name,interactions,articlesInfo,links)

    #enrichment <- string_db$get_enrichment(string_id)
    #info <- list(name,interactions,articles,link,enrichment)

    return(info)}
    )
  
  names(allPairsInfo) <- as.character(interaction(string_idsPairs[,5:6],sep=" & ")) #Rename list with gene symbols names
  
  #Create Excel with genes data
  library("xlsx")
  # Write the first data set in a new workbook

  wb<-createWorkbook(type="xlsx")
  TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=24, color="deeppink3", isBold=TRUE, name="Helvetica")
  SUB_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16, color="darkgray", isBold=TRUE, underline=1,name="Helvetica")
  TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
  TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) + Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") + Border(color="black", position=c("TOP", "BOTTOM"), pen=c("BORDER_THIN", "BORDER_THICK")) 
  
   lapply(allPairsInfo, function(OnePairInfo){
    sheet <- createSheet(wb, sheetName = OnePairInfo[[1]])

    rows <- createRow(sheet, 1:10) # 10 rows
    cells <- createCell(rows, colIndex=1:16) # 16 columns
    setColumnWidth(sheet, 1:16, 15)
    
    #Interactions
    addDataFrame(OnePairInfo[[2]], sheet, startRow=2, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE, row.names = FALSE)
    
    #Hyperlink
    setCellValue(cells[[5,1]], "StringDB Link")
    setCellStyle(cells[[5,1]], SUB_STYLE)
    setCellValue(cells[[6,1]], OnePairInfo[[4]])
    addHyperlink(cells[[6,1]],OnePairInfo[[4]])
    
    #Articles
    addDataFrame(OnePairInfo[[3]], sheet, startRow=8, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE,rownamesStyle = TABLE_ROWNAMES_STYLE, row.names = FALSE)
    setCellValue(cells[[8,1]], "Pubmed Article ID")
    setCellValue(cells[[8,2]], "Pubmed Article Name")
    
    #Title
    xlsx.addTitle(sheet, rowIndex=1, title=OnePairInfo[[1]], titleStyle = TITLE_STYLE)

    # addDataFrame(OnePairInfo[[2]], sheet, startRow=1, startColumn=1)
    # addDataFrame(OnePairInfo[[4]], sheet, startRow=4, startColumn=1)
    # addDataFrame(OnePairInfo[[3]], sheet, startRow=6, startColumn=1)
    #addDataFrame(OnePairInfo[[5]], sheet, startRow=20, startColumn=1)
    

  })
   
  
  saveWorkbook(wb, "GeneData.xlsx") #Save data
  
  #llPairsInfo[["1"]][[3]][1]
  thepage = readLines(paste("https://pubmed.ncbi.nlm.nih.gov/",gsub("PMID:", "", allPairsInfo[["1"]][[3]][2]),"/",sep = ""))
  articleName<-thepage[25]
  articleName<-gsub("    <title>","",articleName)
  articleName<-gsub(" - PubMed</title>","",articleName)
  

  
  library(easyPubMed)
  my_abstracts_txt <- fetch_pubmed_data(OnePairInfo[[3]])
  
  dami_query_string <- " "
  dami_on_pubmed <- get_pubmed_ids(dami_query_string)
  dami_papers <- fetch_pubmed_data(dami_on_pubmed)
  titles <- custom_grep(dami_papers, "ArticleTitle", "char")
  print(titles)
}

