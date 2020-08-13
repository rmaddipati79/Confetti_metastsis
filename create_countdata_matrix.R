###################################
# Author: Komal S Rathi
# Function: Create countdata matrix 
# Date: 02/26/2019
###################################

# Read in formatted data
load("data/Archive/formattedData_mouse850.RData")
countData2 <- readRDS('data/Archive/852_853.RDS')
countData2$gene_id <- gsub('[.].*','',countData2$gene_id)
countData2$gene_symbol <- NULL
countData <- merge(countData, countData2, by.x = 'row.names', by.y = 'gene_id')
rownames(countData) <- countData$Row.names
countData$Row.names <- NULL
colnames(countData) <- gsub('-','_',colnames(countData))
countData <- countData[,order(colnames(countData))]

# sample data
sampleData2 <- read.delim('data/ALL MICE Master samplesheet 072418.txt')
sampleData <- sampleData2
sampleData$Column <- gsub('-','_',sampleData$Column)
rownames(sampleData) <- sampleData$Column
sampleData <- sampleData[order(rownames(sampleData)),]

# save RData
if(identical(rownames(sampleData), colnames(countData))){
  rm(countData2, sampleData2)
  colnames(countData) <- sampleData$Non.Merge.Name
  save.image(file = 'data/all_mouse_matrix_072418.RData')
}