
#Load functions
source("./functions/hypergeometric.R")
source("./functions/HelperFunctions/combinations.R")
source("./functions/HelperFunctions/countDuplicates.R")

source("./functions/HelperFunctions/countDuplicatesNumber.R")
source("./functions/fromSL.R")
source("./functions/investigate.R")
source("./functions/HelperFunctions/xlsx.addTitle.R")
source("./functions/HelperFunctions/ORA.R")

source("./functions/HelperFunctions/weighted.R")

## Load Data ------------------------------------------------------------

library(Matrix)
library(readr)
library(dplyr) #distinct
library(plyr) #ddply
library(qdapTools) #list_df2df
library(prodlim)
library(stringr)
library(WebGestaltR)

#Import synthetic lethal data
sl_human <- read_delim("./input_data/sl_human",
                       "\t", escape_double = FALSE, col_names = FALSE,  col_types="cdcdccccd",
                       trim_ws = TRUE)
colnames(sl_human) <- c("HGNC_1","Entrez_1","HGNC_2","Entrez_2","PubmedID","Source","TypeSL","tissue-cell_line","Score")
namessl1 <- sl_human$Entrez_1
namessl2 <- sl_human$Entrez_2

# missnames1 <- which(is.na(match(namessl1,namesGO)))
# missnames2 <- which(is.na(match(namessl2,namesGO)))
# missnames <- unique(c(missnames1,missnames2))
# sl_humanOutput <- sl_human[-missnames,]

## SELECT THE SOURCE!
#sl_humanOutput <- sl_human[grep("Mining",sl_human$Source),]
sl_humanOutput <- sl_human
pairsSL <- sl_humanOutput[,c(2,4)]
pairsSL<- t(apply(t(apply(pairsSL, 1, as.numeric)),1,sort)) # Sort in ascending order
pairsSL <- distinct(as.data.frame(pairsSL)) # Make sure there are no repeated pairs 
colnames(pairsSL) <- c("gene1","gene2") 


# Import minimal cut sets
MCS <- as.data.frame(read.csv(file="./input_data/CHECK_structural_Recon3D_50000gMCSs_60seconds_ENTREZ_20190613.csv", header=FALSE, sep=","))
MCS <- MCS[-(which((is.na(MCS$V2)==TRUE))),] #Delete MCS which contain only one gene (can't make a pair)

# Import all metabolic genes
allGenesData <- as.data.frame(read.csv(file="./input_data/MetabolicGenes.csv", header=TRUE, sep=";"))
allGenes<- unique(allGenesData[,16]) # Only keep the entrez column

## Part one ------------------------------------------------------------
probability <- hypergeometric(allGenes,MCS,pairsSL) 



