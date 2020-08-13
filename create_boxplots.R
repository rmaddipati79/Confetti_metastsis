# Author: Komal S. Rathi
# Date: 05/13/2020
# Function: Create boxplots
library(tidyverse)
library(sva)
library(ggpubr)

# boxplots
source('code/Utils/pubTheme.R')
source("code/Utils/boxplot_accessory.R")

# load FPKM matrix (n = 54)
load('data/all_mouse_matrix_FPKM_fullset.RData') 

# genes of interest
goi <- read.delim('data/genesofinterest.txt', stringsAsFactors = F) # add some genes of interest
outputLimma_pmetvspnmet.volcano <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt', stringsAsFactors = F)
deGenes <- outputLimma_pmetvspnmet.volcano %>%
  .$gene_symbol
deGenes <- unique(c(deGenes, goi$gene_symbols))

# use geneid + genesym for plotting
geneids <- geneAnnot[which(geneAnnot$gene_symbol %in% deGenes),'gene_id']

# subset matrix to GOI and take log2(FPKM + 1) and then ComBat
countDataTrimSVA <- log2(countData[which(rownames(countData) %in% outputLimma_pmetvspnmet.volcano$gene_id),] + 1)
mouseID <- factor(sampleData$Mouse.ID)
countDataTrimSVA <- ComBat(as.matrix(countDataTrimSVA), batch = mouseID, par.prior = T)

# annotate p-values for multiple groups
for(i in 1:length(geneids)){
  print(i)
  id <- geneAnnot[which(geneAnnot$gene_id == geneids[i]),1]
  sym <- geneAnnot[which(geneAnnot$gene_id == geneids[i]),2]
  gene <- paste0(sym,'_',id)
  fname <- paste0("results/Boxplots/allMice/withmets/", gene, "_annotated.eps")
  boxplot.multiple.groups(myMat = countDataTrimSVA, mySampData = sampleData, myGeneID = id, myGeneSym = sym, mets = 'yes', test = 't.test')
  ggsave(filename = fname, width = 7, height = 5)
}

for(i in 1:length(geneids)){
  print(i)
  id <- geneAnnot[which(geneAnnot$gene_id == geneids[i]),1]
  sym <- geneAnnot[which(geneAnnot$gene_id == geneids[i]),2]
  gene <- paste0(sym,'_',id)
  fname <- paste0("results/Boxplots/allMice/nomets/", gene, "_annotated.eps")
  boxplot.multiple.groups(myMat = countDataTrimSVA, mySampData = sampleData, myGeneID = id, myGeneSym = sym, mets = 'no', test = 't.test')
  ggsave(filename = fname, width = 7, height = 5)
}