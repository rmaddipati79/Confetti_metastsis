###############################
# Author: Komal S Rathi
# Function: 
# Map Mouse to Human genes and 
# Create input matrix for 
# Immune signatures/Gene set enrichments etc
# Date: 12/09/2019
###############################

library(biomaRt)
library(sva)
library(tidyverse)
library(limma)
library(edgeR)

setwd('~/Projects/PancModelUpdate/')

# function to format matrix of counts
collapse.rnaseq <- function(countData, sampleData, outfile, type) {
  if(identical(colnames(countData), rownames(sampleData))){
    rownames(sampleData) <- paste0('PD',sampleData$Mouse.ID, '_', sampleData$Unique.Identifier, '_', sampleData$Sample.Type, sampleData$Organ.site, sampleData$Color,'_',sampleData$Relationship.between.primary.tumors,'',sampleData$GivesRiseToMet,sampleData$Primary.to.Metastasis.relationship)
    colnames(countData) <- rownames(sampleData)
  }
  
  # remove non-liver mets (n = 67)
  sampleData$PrimCR <- paste(sampleData$Sample.Type, sampleData$GivesRiseToMet, sep = "_")
  sampleData$PrimCRLab <- ifelse(sampleData$PrimCR == "P_0", "Primary_NoMet", ifelse(sampleData$PrimCR == "P_1", "Primary_Met", "Met"))
  sampleData$PrimCRLab <- factor(sampleData$PrimCRLab, levels = c("Primary_NoMet", "Primary_Met", "Met"))
  if(length(which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li")) == 0){
    sampleDataPrim <- sampleData
  } else {
    sampleDataPrim <- sampleData[-which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li"),]
  }
  countDataPrim <- countData[,rownames(sampleDataPrim)]
  
  # remove genes low expressing genes
  maxValueGenes <- apply(countDataPrim, FUN = max, MARGIN = 1)
  if(type == "counts"){
    expressedGenes <- names(maxValueGenes[maxValueGenes > 100])
  } else {
    expressedGenes <- names(maxValueGenes[maxValueGenes > 10])
  }
  countDataPrim <- countDataPrim[expressedGenes,]
  countDataPrim <- countDataPrim[,sort(colnames(countDataPrim))]
  
  # find all things to be merged
  sampleDataPrim$Primary_clone <- paste0(sampleDataPrim$Mouse.ID,'_',sampleDataPrim$Relationship.between.primary.tumors)
  
  # separate Mets because they are not to be merged (n = 19)
  sampleData.onlyMets <- sampleDataPrim[which(sampleDataPrim$PrimCRLab == "Met"),]
  countData.onlyMets <- countDataPrim[,rownames(sampleData.onlyMets)]
  
  # find all things to be merged
  sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
  countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  
  # sort 
  countDataPrim <- countDataPrim[,sort(colnames(countDataPrim))]
  primaryCloneCounts <- plyr::count(sampleDataPrim$Primary_clone)
  clonesToMerge <- as.character(primaryCloneCounts[primaryCloneCounts$freq > 1, 'x'])
  
  for(i in 1:length(clonesToMerge)) {
    # Merge rows in sample data
    sampsToMerge <- rownames(sampleDataPrim[sampleDataPrim[,"Primary_clone"] %in% clonesToMerge[i],])
    firstSamp <- sampsToMerge[1]
    otherSamps <- sampsToMerge[2:length(sampsToMerge)]
    sampleDataPrim <- sampleDataPrim[setdiff(colnames(countDataPrim), otherSamps),] # remove other rows
    
    # Merge columns in count data
    tmpCol <- rowMeans(countDataPrim[,sampsToMerge])
    countDataPrim[,firstSamp] <- tmpCol # Replace the first sample with mean
    countDataPrim <- countDataPrim[, setdiff(colnames(countDataPrim), otherSamps)] # Remove other sample
  }
  
  # merge merged Primary and non-merged Mets
  sampleDataPrim  <- rbind(sampleDataPrim, sampleData.onlyMets)
  if(identical(rownames(countDataPrim), rownames(countData.onlyMets))){
    countDataPrim <- cbind(countDataPrim, countData.onlyMets)
  }
  
  # normalize using voom
  # batch correction by mouse identifier
  clonRel <- factor(sampleDataPrim$PrimCRLab)
  design <- model.matrix(~clonRel)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  # only do voom normalization on count data
  if(type == "counts"){
    y <- DGEList(counts = as.matrix(countDataPrim), genes = rownames(countDataPrim))
    y <- calcNormFactors(y)
    v <- voom(counts = y, design = design, plot=F)
    voomData <- v$E
    countDataPrim <- voomData
    countDataPrim <- ComBat(dat = voomData, batch = mouseID , mod = design, par.prior = T)
  } 
  
  # function to convert data to human genes
  # also output gene type (to do)
  convertMouseGeneList <- function(x){
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesV2 = getLDS(attributes = c("mgi_symbol"), 
                     filters = "mgi_symbol", 
                     values = x, 
                     mart = mouse, 
                     attributesL = c("hgnc_symbol","gene_biotype"), 
                     martL = human, uniqueRows=T)
    return(genesV2)
  }
  
  # map human to mouse
  # geneAnnot$gene_type <- NULL
  countDataPrim <- merge(countDataPrim, geneAnnot, by.x = 'row.names', by.y = 'gene_id')
  
  # mouse matrix
  countDataPrim.mouse <- countDataPrim %>%
    filter(gene_type == "protein_coding") %>%
    dplyr::select(-c(Row.names, gene_type)) %>%
    dplyr::mutate(means = rowMeans(.[1:nrow(sampleDataPrim)])) %>%
    arrange(desc(means)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    dplyr::select(-c(means)) %>%
    column_to_rownames(var = "gene_symbol")
  
  human.genes <- convertMouseGeneList(countDataPrim$gene_symbol)
  countDataPrim <- merge(countDataPrim, human.genes, by.x = 'gene_symbol', by.y = 'MGI.symbol')
  countDataPrim <- countDataPrim %>% 
    filter(Gene.type == "protein_coding") %>%
    dplyr::select(-c(gene_symbol, gene_type, Gene.type, Row.names)) %>%
    dplyr::mutate(means = rowMeans(.[1:nrow(sampleDataPrim)])) %>% 
    arrange(desc(means)) %>% 
    distinct(HGNC.symbol, .keep_all = TRUE) %>% 
    filter(HGNC.symbol != "") %>%
    dplyr::select(-c(means)) %>%
    column_to_rownames(var = "HGNC.symbol")
  sampleDataPrim <- sampleDataPrim[,c("Non.Merge.Name", "Merge.Name", "Mouse.ID", "PrimCRLab")]
  
  save(countDataPrim.mouse, countDataPrim, sampleDataPrim, human.genes, file = outfile)
}

# count data
load('data/all_mouse_matrix_072418.RData')
collapse.rnaseq(countData = countData, 
                sampleData = sampleData, 
                outfile = 'data/collapsed_counts_to_human.RData', 
                type = "counts")

# fpkm data
load('data/all_mouse_matrix_FPKM_072418.RData')
collapse.rnaseq(countData = countData, 
                sampleData = sampleData, 
                outfile = 'data/collapsed_fpkm_to_human.RData', 
                type = "fpkm")


