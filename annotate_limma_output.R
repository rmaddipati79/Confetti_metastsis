##################################################
# Author: Komal S Rathi
# Function: Annotate Limma output with MYC binding and Expression  
# Date: 07/14/2019
##################################################

setwd('~/Projects/PancModelUpdate/')
library(reshape2)
library(dplyr)

# data with all biotypes
myc.binding <- read.delim('results/MYC_binding/PrimaryMetvsPrimaryNonMet_allMice_MYCBinding_allBiotypes.txt')
limma.out <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes.txt')

covar <- function(x) {
  sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE)
}

# 1. add FPKM values from actual data
# this has already been merged for clones
load('data/all_mouse_matrix_FPKM_fullset.RData')

# just primary tumors
sampleData <- sampleData[which(sampleData$PrimCRLab != 'Met'),]
sampleData <- sampleData[,c('Merge.Name','PrimCRLab')]

load('results/Limma/PrimaryMet_vs_PrimaryNonMet/pmetvspnmet.RData')
countData_pmetvspnmet <- res[[1]]

# load('data/all_mouse_matrix_FPKM_fullset.RData')
countData <- countData[rownames(countData) %in% rownames(countData_pmetvspnmet),colnames(countData) %in% colnames(countData_pmetvspnmet)]
countData <- melt(as.matrix(countData), varnames = c('gene_id','sample'), value.name = 'FPKM')
countData <- merge(countData, sampleData, by.x = 'sample', by.y = 'row.names')
countData <- countData %>% group_by(gene_id, PrimCRLab) %>%
  summarise(mean = mean(FPKM),
            sd = sd(FPKM),
            var = covar(FPKM),
            count = n()) %>% 
  as.data.frame()
primarymet <- countData[which(countData$PrimCRLab == "Primary_Met"),]
primarynmet <- countData[which(countData$PrimCRLab == "Primary_NoMet"),]

colnames(primarymet) <- paste0('PrimaryMet_',colnames(primarymet))
colnames(primarynmet) <- paste0('PrimaryNMet_',colnames(primarynmet))


# 2. add FPKM values from black 6 cell lines
load('data/black6mouse_mm10_FPKM_merged.RData')
black6.matrix <- as.data.frame(black6.matrix)
black6.matrix <- black6.matrix[rownames(black6.matrix) %in% rownames(countData_pmetvspnmet),]
black6.matrix <- melt(as.matrix(black6.matrix), varnames = c('gene_id','sample'), value.name = 'FPKM')
black6.matrix <- merge(black6.matrix, black6.ss, by.x = 'sample', by.y = 'rep')
black6.matrix <- black6.matrix %>% group_by(gene_id, PrimCRLab) %>%
  summarise(mean = mean(FPKM),
            sd = sd(FPKM),
            var = covar(FPKM),
            count = n()) %>% 
  as.data.frame()
tcellhigh <- black6.matrix[which(black6.matrix$PrimCRLab == "T cell high"),]
tcelllow <- black6.matrix[which(black6.matrix$PrimCRLab == "T cell low"),]

colnames(tcellhigh) <- paste0('B6_CL_TcellHigh_',colnames(tcellhigh))
colnames(tcelllow) <- paste0('B6_CL_TcellLow_',colnames(tcelllow))

# 3. add FPKM values from black 6 bulk tumor
load('data/black6mouse_bulk_mm10_FPKM_merged.RData')
black6.bulk <- black6.bulk[rownames(black6.bulk) %in% rownames(countData_pmetvspnmet),]
black6.bulk <- melt(as.matrix(black6.bulk), varnames = c('gene_id','sample'), value.name = 'FPKM')
black6.bulk <- merge(black6.bulk, black6.ss, by.x = 'sample', by.y = 'Clones')
black6.bulk <- black6.bulk %>% group_by(gene_id, PrimCRLab) %>%
  summarise(mean = mean(FPKM),
            sd = sd(FPKM),
            var = covar(FPKM),
            count = n()) %>% 
  as.data.frame()
bulk.tcellhigh <- black6.bulk[which(black6.bulk$PrimCRLab == "T cell high"),]
bulk.tcelllow <- black6.bulk[which(black6.bulk$PrimCRLab == "T cell low"),]

colnames(bulk.tcellhigh) <- paste0('B6_Bulk_TcellHigh_',colnames(bulk.tcellhigh))
colnames(bulk.tcelllow) <- paste0('B6_Bulk_TcellLow_',colnames(bulk.tcelllow))

to.add <- cbind(primarymet, 
      primarynmet[,2:ncol(primarynmet)],
      tcellhigh[,2:ncol(tcellhigh)],
      tcelllow[,2:ncol(tcelllow)],
      bulk.tcellhigh[,2:ncol(bulk.tcellhigh)],
      bulk.tcelllow[,2:ncol(bulk.tcelllow)])
res <- merge(myc.binding, to.add, by.x = 'gene_id', by.y = 'PrimaryMet_gene_id')
res <- res[,!colnames(res) %in% c('AveExpr','t','B')]
write.table(res, file = 'results/MYC_binding/PrimaryMetvsPrimaryNonMet_allMice_MYCBinding_allBiotypes_withExpression.txt', quote = F, sep = "\t", row.names = F)
