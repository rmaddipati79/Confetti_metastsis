setwd('~/Projects/PancModelUpdate/')
library(ggplot2)
library(scatterplot3d)
library(pheatmap)
library(dplyr)
library(sva)

# Only protein coding genes were used. 
# Clones were merged but not for individual mice. 
# Where there is no heatmap available, the dimensions were insufficient i.e. < 3

# matrix
load('data/all_mouse_matrix_072418.RData')

# functions to create heatmaps
# primary mets, primary non-mets and mets
createSampleHeatmap <- function(mat, mets, fname){
  
  # Read in formatted data
  load("data/all_mouse_matrix_072418.RData")
  sampleData[, ] <- lapply(sampleData[, ], as.character)
  
  ######################################################
  # 1. Format more and QC Plots
  ######################################################
  if(identical(colnames(countData), rownames(sampleData))){
    rownames(sampleData) <- paste0('PD',sampleData$Mouse.ID, '_', sampleData$Unique.Identifier, '_', sampleData$Sample.Type, sampleData$Organ.site, sampleData$Color,'_',sampleData$Relationship.between.primary.tumors,'',sampleData$GivesRiseToMet,sampleData$Primary.to.Metastasis.relationship)
    colnames(countData) <- rownames(sampleData)
  }
  
  # 1.a define groups
  # remove non-liver mets
  sampleData$PrimCR <- paste(sampleData$Sample.Type, sampleData$GivesRiseToMet, sep = "_")
  sampleData$PrimCRLab <- ifelse(sampleData$PrimCR == "P_0", "Primary_NoMet", ifelse(sampleData$PrimCR == "P_1", "Primary_Met", "Met"))
  sampleData$PrimCRLab <- factor(sampleData$PrimCRLab, levels = c("Primary_NoMet", "Primary_Met", "Met"))
  if(length(which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li")) == 0) {
    sampleDataPrim <- sampleData
  } else {
    sampleDataPrim <- sampleData[-which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li"),]
  }
  countDataPrim <- countData[,rownames(sampleDataPrim)]
  
  # subset matrix (either all samples or individual mice)
  if(mat == "all"){
    print("Using all mice...")
  } else {
    print("Use the specified mouse id")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$Mouse.ID == mat),]
    countDataPrim <- countData[,rownames(sampleDataPrim)]
  }
  
  # if mets are to be merged or not
  if(mets == "yes"){
    print("Using mets as well...")
  } else {
    print("Remove mets...")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # remove genes low expressing genes
  maxValueGenes <- apply(countDataPrim, FUN = max, MARGIN = 1)
  expressedGenes <- names(maxValueGenes[maxValueGenes > 100])
  countDataPrim <- countDataPrim[expressedGenes,]
  countDataPrim <- countDataPrim[,sort(colnames(countDataPrim))]
  
  # 1.b Collapse columns & Annotations
  # Really just keeping first row for sample data
  # For count data calculating mean
  sampleDataPrim$Primary_clone <- paste0(sampleDataPrim$Mouse.ID,'_',sampleDataPrim$Relationship.between.primary.tumors)
  
  # separate Mets because they are not to be merged (n = 19)
  if(mets == "yes"){
    print("Lets separate the mets in another object...")
    sampleData.onlyMets <- sampleDataPrim[which(sampleDataPrim$PrimCRLab == "Met"),]
    countData.onlyMets <- countDataPrim[,rownames(sampleData.onlyMets)]
    
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # find all things to be merged
  # only merge clones if all samples are involved
  if(mat == "all"){
    primaryCloneCounts <- plyr::count(sampleDataPrim$Primary_clone)
    clonesToMerge <- as.character(primaryCloneCounts[primaryCloneCounts$freq > 1, 'x'])
    
    for(i in 1:length(clonesToMerge)) {
      # Merge rows in sample data
      sampsToMerge <- rownames(sampleDataPrim[sampleDataPrim$Primary_clone %in% clonesToMerge[i],])
      firstSamp <- sampsToMerge[1]
      otherSamps <- sampsToMerge[2:length(sampsToMerge)]
      sampleDataPrim <- sampleDataPrim[setdiff(colnames(countDataPrim), otherSamps),] # remove other rows
      
      # Merge columns in count data
      tmpCol <- rowMeans(countDataPrim[,sampsToMerge])
      countDataPrim[,firstSamp] <- tmpCol # Replace the first sample with mean
      countDataPrim <- countDataPrim[, setdiff(colnames(countDataPrim), otherSamps)] # Remove other sample
    }
  }
  
  # merge merged Primary and non-merged Mets
  if(mets == "yes"){
    print("Merge back the mets in...")
    sampleDataPrim  <- rbind(sampleDataPrim, sampleData.onlyMets)
    if(identical(rownames(countDataPrim), rownames(countData.onlyMets))){
      countDataPrim <- cbind(countDataPrim, countData.onlyMets)
    }
  }
  
  # final manipulation (add shorter names)
  if(identical(colnames(countDataPrim), rownames(sampleDataPrim))){
    if(mat == "all"){
      rownames(sampleDataPrim) <- sampleDataPrim$Merge.Name
      colnames(countDataPrim) <- rownames(sampleDataPrim)
    } else {
      rownames(sampleDataPrim) <- sampleDataPrim$Non.Merge.Name
      colnames(countDataPrim) <- rownames(sampleDataPrim)
    }
  }
  
  ########################################
  # Correlation plot between samples 
  ########################################
  sampleDataPrim$Group <- factor(sampleDataPrim$PrimCRLab)
  fname <- paste0('results/heatmaps_allmice/', fname, '.pdf')
  pdf(file = fname, width = 10, height = 10)
  myCor <- cor(log2(countDataPrim+1))
  myCorNames <- rownames(myCor)
  pheatmap(myCor,annotation_row = sampleDataPrim[c("Group","Mouse.ID")], main = "Sample Correlation Heatmap\n", display_numbers = TRUE)
  dev.off()
  
  # keep only protein-coding genes, remove Gm and Riken genes
  genes <- geneAnnot[which(geneAnnot$gene_type == "protein_coding"),]
  genes <- genes[grep('^Gm|Rik$', genes$gene_symbol, invert = T),]
  countDataPrim <- countDataPrim[which(rownames(countDataPrim) %in% genes$gene_id),]
  
  return(list(countDataPrim, sampleDataPrim))
}

# n is for number of variable genes
# matrix and sample sheet come from the function above
createGeneHeatmap <- function(countDataPrim, sampleDataPrim, fname, title, n = NULL){
  
  # combat adjust all protein coding genes
  countDataTrim <- log2(countDataPrim+1)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  if(length(levels(factor(mouseID))) > 1){
    print("More than 2 levels get combat adjustment...")
    countDataTrim <- ComBat(as.matrix(countDataTrim), batch=mouseID, par.prior=T)
  }
  
  # now subset depending on N
  if(is.null(n)){
    print("Use all expressed protein coding genes...")
  } else {
    print("Use n most-variable protein coding genes...")
    # calculate the n most variable genes
    var_genes <- apply(countDataTrim, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:n]
    highly_variable_lcpm <- countDataTrim[select_var,]
  }
  
  sampleDataPrim$Group <- factor(sampleDataPrim$PrimCRLab)
  fname <- paste0('results/heatmaps_allmice/', fname, '.pdf')
  pdf(file = fname, width = 10, height = 10)
  # pheatmap
  if(is.null(n)){
    # for all expressed genes
    print("All Expressed Dimensions:")
    print(dim(countDataTrim))
    if(ncol(countDataTrim) > 2){
      pheatmap(countDataTrim, 
               scale="row", 
               show_rownames=F, 
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  } else {
    # n most variable genes
    print("Most Variable Dimensions:")
    print(dim(highly_variable_lcpm))
    if(ncol(highly_variable_lcpm) > 2){
      pheatmap(highly_variable_lcpm, 
               scale="row", 
               show_rownames=F, 
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  }
  dev.off()
}

# unclustered columns
# n is for number of variable genes
# matrix and sample sheet come from the function above
# first Primary Mets and then Primary non-mets
# Sort by number
createGeneHeatmap.unclustered <- function(countDataPrim, sampleDataPrim, fname, title, n = NULL){
  
  # combat adjust all protein coding genes
  countDataTrim <- log2(countDataPrim+1)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  if(length(levels(factor(mouseID))) > 1){
    print("More than 2 levels get combat adjustment...")
    countDataTrim <- ComBat(as.matrix(countDataTrim), batch=mouseID, par.prior=T)
  }
  sampleDataPrim$PrimCRLab <- factor(sampleDataPrim$PrimCRLab, levels = c("Primary_Met","Primary_NoMet","Met"))
  sampleDataPrim <- sampleDataPrim[order(sampleDataPrim$PrimCRLab, sampleDataPrim$Merge.Name),]
  countDataTrim <- countDataTrim[,sampleDataPrim$Merge.Name]
  
  # now subset depending on N
  if(is.null(n)){
    print("Use all expressed protein coding genes...")
  } else {
    print("Use n most-variable protein coding genes...")
    # calculate the n most variable genes
    var_genes <- apply(countDataTrim, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:n]
    highly_variable_lcpm <- countDataTrim[select_var,]
  }
  
  sampleDataPrim$Group <- factor(sampleDataPrim$PrimCRLab)
  fname <- paste0('results/heatmaps_allmice/', fname, '.pdf')
  pdf(file = fname, width = 10, height = 10)
  # pheatmap
  if(is.null(n)){
    # for all expressed genes
    print("All Expressed Dimensions:")
    print(dim(countDataTrim))
    if(ncol(countDataTrim) > 2){
      pheatmap(countDataTrim, 
               scale="row", 
               show_rownames=F, cluster_cols = FALSE,
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  } else {
    # n most variable genes
    print("Most Variable Dimensions:")
    print(dim(highly_variable_lcpm))
    if(ncol(highly_variable_lcpm) > 2){
      pheatmap(highly_variable_lcpm, 
               scale="row", 
               show_rownames=F, cluster_cols = FALSE,
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  }
  dev.off()
}

### clustered
## all mice + all samples 
system('mkdir -p ~/Projects/PancModelUpdate/results/heatmaps_allmice/allmice_with_mets')
res <- createSampleHeatmap(mat = 'all', mets = 'yes', fname = 'allmice_with_mets/allmice_sample_heatmap')
countData.forHeatmap <- res[[1]]
sampleData.forHeatmap <- res[[2]]

# all expressed genes
createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                  sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                  fname = 'allmice_with_mets/allmice_expressed_heatmap', 
                  title = "Heatmap All Expressed Genes")

# n variable genes
nums <- c(250, 500, 1000, 2000, 5000)
for(j in 1:length(nums)){
  createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                    sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                    fname = paste0('allmice_with_mets/allmice_top', nums[j], '_heatmap'), 
                    title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
}

## all mice + no mets 
system('mkdir -p ~/Projects/PancModelUpdate/results/heatmaps_allmice/allmice_without_mets')
res <- createSampleHeatmap(mat = 'all', mets = 'no', fname = 'allmice_without_mets/allmice_sample_nomets_heatmap')
countData.forHeatmap <- res[[1]]
sampleData.forHeatmap <- res[[2]]

# all expressed genes
createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                  sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                  fname = 'allmice_without_mets/allmice_expressed_nomets_heatmap', 
                  title = "Heatmap All Expressed Genes")

# n variable genes
for(j in 1:length(nums)){
  createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                    sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                    fname = paste0('allmice_without_mets/allmice_top', nums[j], '_nomets_heatmap'), 
                    title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
}

# individual mice
mice <- unique(sampleData$Mouse.ID)
for(i in 1:length(mice)){
  nums <- c(250, 500, 1000, 2000, 5000)
  system(paste0('mkdir -p ~/Projects/PancModelUpdate/results/heatmaps_allmice/', mice[i]))
  # with mets
  res <- createSampleHeatmap(mat = mice[i], mets = 'yes', fname = paste0(mice[i],'/',mice[i],'_sample_heatmap'))
  countData.forHeatmap <- res[[1]]
  sampleData.forHeatmap <- res[[2]]
  
  # all expressed genes
  createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                    sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                    fname = paste0(mice[i],'/',mice[i],'_expressed_heatmap'), 
                    title = "Heatmap All Expressed Genes")
  
  # variable genes
  for(j in 1:length(nums)){
    createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                      sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                      fname = paste0(mice[i],'/',mice[i],'_top', nums[j], '_heatmap'), 
                      title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
  }
  
  # without mets
  res <- createSampleHeatmap(mat = mice[i], mets = 'no', fname = paste0(mice[i],'/',mice[i],'_sample_nomets_heatmap'))
  countData.forHeatmap <- res[[1]]
  sampleData.forHeatmap <- res[[2]]
  
  # all expressed genes
  createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                    sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                    fname = paste0(mice[i],'/',mice[i],'_expressed_nomets_heatmap'), 
                    title = "Heatmap All Expressed Genes")
  
  # n most variable
  for(j in 1:length(nums)){
    createGeneHeatmap(countDataPrim = countData.forHeatmap, 
                      sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                      fname = paste0(mice[i],'/',mice[i],'_top', nums[j], '_nomets_heatmap'), 
                      title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
  }
}

### unclustered
## all mice + all samples 
system('mkdir -p ~/Projects/PancModelUpdate/results/heatmaps_allmice/allmice_with_mets_unclustered')
res <- createSampleHeatmap(mat = 'all', mets = 'yes', fname = 'allmice_with_mets_unclustered/allmice_sample_heatmap')
countData.forHeatmap <- res[[1]]
sampleData.forHeatmap <- res[[2]]

# all expressed genes
createGeneHeatmap.unclustered(countDataPrim = countData.forHeatmap, 
                              sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                              fname = 'allmice_with_mets_unclustered/allmice_expressed_heatmap', 
                              title = "Heatmap All Expressed Genes")

# n variable genes
nums <- c(250, 500, 1000, 2000, 5000)
for(j in 1:length(nums)){
  createGeneHeatmap.unclustered(countDataPrim = countData.forHeatmap, 
                                sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                                fname = paste0('allmice_with_mets_unclustered/allmice_top', nums[j], '_heatmap'), 
                                title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
}

## all mice + no mets 
system('mkdir -p ~/Projects/PancModelUpdate/results/heatmaps_allmice/allmice_without_mets_unclustered')
res <- createSampleHeatmap(mat = 'all', mets = 'no', fname = 'allmice_without_mets_unclustered/allmice_sample_nomets_heatmap')
countData.forHeatmap <- res[[1]]
sampleData.forHeatmap <- res[[2]]

# all expressed genes
createGeneHeatmap.unclustered(countDataPrim = countData.forHeatmap, 
                              sampleDataPrim = sampleData.forHeatmap, n = NULL, 
                              fname = 'allmice_without_mets_unclustered/allmice_expressed_nomets_heatmap', 
                              title = "Heatmap All Expressed Genes")

# n variable genes
for(j in 1:length(nums)){
  createGeneHeatmap.unclustered(countDataPrim = countData.forHeatmap, 
                                sampleDataPrim = sampleData.forHeatmap, n = nums[j], 
                                fname = paste0('allmice_without_mets_unclustered/allmice_top', nums[j], '_nomets_heatmap'), 
                                title = paste0('Heatmap Top ', nums[j], ' Most Variable Genes'))
}

# for degs only, unclustered
createSampleHeatmap.allbiotypes <- function(mat, mets, fname){
  
  # Read in formatted data
  load("data/all_mouse_matrix_072418.RData")
  sampleData[, ] <- lapply(sampleData[, ], as.character)
  
  ######################################################
  # 1. Format more and QC Plots
  ######################################################
  if(identical(colnames(countData), rownames(sampleData))){
    rownames(sampleData) <- paste0('PD',sampleData$Mouse.ID, '_', sampleData$Unique.Identifier, '_', sampleData$Sample.Type, sampleData$Organ.site, sampleData$Color,'_',sampleData$Relationship.between.primary.tumors,'',sampleData$GivesRiseToMet,sampleData$Primary.to.Metastasis.relationship)
    colnames(countData) <- rownames(sampleData)
  }
  
  # 1.a define groups
  # remove non-liver mets
  sampleData$PrimCR <- paste(sampleData$Sample.Type, sampleData$GivesRiseToMet, sep = "_")
  sampleData$PrimCRLab <- ifelse(sampleData$PrimCR == "P_0", "Primary_NoMet", ifelse(sampleData$PrimCR == "P_1", "Primary_Met", "Met"))
  sampleData$PrimCRLab <- factor(sampleData$PrimCRLab, levels = c("Primary_NoMet", "Primary_Met", "Met"))
  if(length(which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li")) == 0) {
    sampleDataPrim <- sampleData
  } else {
    sampleDataPrim <- sampleData[-which(sampleData$PrimCRLab == "Met" & sampleData$Organ.site != "Li"),]
  }
  countDataPrim <- countData[,rownames(sampleDataPrim)]
  
  # subset matrix (either all samples or individual mice)
  if(mat == "all"){
    print("Using all mice...")
  } else {
    print("Use the specified mouse id")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$Mouse.ID == mat),]
    countDataPrim <- countData[,rownames(sampleDataPrim)]
  }
  
  # if mets are to be merged or not
  if(mets == "yes"){
    print("Using mets as well...")
  } else {
    print("Remove mets...")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # remove genes low expressing genes
  maxValueGenes <- apply(countDataPrim, FUN = max, MARGIN = 1)
  expressedGenes <- names(maxValueGenes[maxValueGenes > 100])
  countDataPrim <- countDataPrim[expressedGenes,]
  countDataPrim <- countDataPrim[,sort(colnames(countDataPrim))]
  
  # 1.b Collapse columns & Annotations
  # Really just keeping first row for sample data
  # For count data calculating mean
  sampleDataPrim$Primary_clone <- paste0(sampleDataPrim$Mouse.ID,'_',sampleDataPrim$Relationship.between.primary.tumors)
  
  # separate Mets because they are not to be merged (n = 19)
  if(mets == "yes"){
    print("Lets separate the mets in another object...")
    sampleData.onlyMets <- sampleDataPrim[which(sampleDataPrim$PrimCRLab == "Met"),]
    countData.onlyMets <- countDataPrim[,rownames(sampleData.onlyMets)]
    
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # find all things to be merged
  # only merge clones if all samples are involved
  if(mat == "all"){
    primaryCloneCounts <- plyr::count(sampleDataPrim$Primary_clone)
    clonesToMerge <- as.character(primaryCloneCounts[primaryCloneCounts$freq > 1, 'x'])
    
    for(i in 1:length(clonesToMerge)) {
      # Merge rows in sample data
      sampsToMerge <- rownames(sampleDataPrim[sampleDataPrim$Primary_clone %in% clonesToMerge[i],])
      firstSamp <- sampsToMerge[1]
      otherSamps <- sampsToMerge[2:length(sampsToMerge)]
      sampleDataPrim <- sampleDataPrim[setdiff(colnames(countDataPrim), otherSamps),] # remove other rows
      
      # Merge columns in count data
      tmpCol <- rowMeans(countDataPrim[,sampsToMerge])
      countDataPrim[,firstSamp] <- tmpCol # Replace the first sample with mean
      countDataPrim <- countDataPrim[, setdiff(colnames(countDataPrim), otherSamps)] # Remove other sample
    }
  }
  
  # merge merged Primary and non-merged Mets
  if(mets == "yes"){
    print("Merge back the mets in...")
    sampleDataPrim  <- rbind(sampleDataPrim, sampleData.onlyMets)
    if(identical(rownames(countDataPrim), rownames(countData.onlyMets))){
      countDataPrim <- cbind(countDataPrim, countData.onlyMets)
    }
  }
  
  # final manipulation (add shorter names)
  if(identical(colnames(countDataPrim), rownames(sampleDataPrim))){
    if(mat == "all"){
      rownames(sampleDataPrim) <- sampleDataPrim$Merge.Name
      colnames(countDataPrim) <- rownames(sampleDataPrim)
    } else {
      rownames(sampleDataPrim) <- sampleDataPrim$Non.Merge.Name
      colnames(countDataPrim) <- rownames(sampleDataPrim)
    }
  }
  
  ########################################
  # Correlation plot between samples 
  ########################################
  sampleDataPrim$Group <- factor(sampleDataPrim$PrimCRLab)
  fname <- paste0('results/heatmaps_allmice/', fname, '.pdf')
  pdf(file = fname, width = 10, height = 10)
  myCor <- cor(log2(countDataPrim+1))
  myCorNames <- rownames(myCor)
  pheatmap(myCor,annotation_row = sampleDataPrim[c("Group","Mouse.ID")], main = "Sample Correlation Heatmap\n", display_numbers = TRUE)
  dev.off()
  
  # keep all biotypes
  # keep only protein-coding genes, remove Gm and Riken genes
  # genes <- geneAnnot[which(geneAnnot$gene_type == "protein_coding"),]
  # genes <- genes[grep('^Gm|Rik$', genes$gene_symbol, invert = T),]
  # countDataPrim <- countDataPrim[which(rownames(countDataPrim) %in% genes$gene_id),]
  
  return(list(countDataPrim, sampleDataPrim))
}

createGeneHeatmap.unclustered.custom <- function(countDataPrim, sampleDataPrim, fname, title, genelist = NULL){
  
  # combat adjust all genes
  countDataTrim <- log2(countDataPrim+1)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  if(length(levels(factor(mouseID))) > 1){
    print("More than 2 levels get combat adjustment...")
    countDataTrim <- ComBat(as.matrix(countDataTrim), batch=mouseID, par.prior=T)
  }
  sampleDataPrim$PrimCRLab <- factor(sampleDataPrim$PrimCRLab, levels = c("Primary_Met","Primary_NoMet","Met"))
  sampleDataPrim <- sampleDataPrim[order(sampleDataPrim$PrimCRLab, sampleDataPrim$Merge.Name),]
  countDataTrim <- countDataTrim[,sampleDataPrim$Merge.Name]
  
  # now subset depending on genelist
  if(is.null(genelist)){
    print("Use from above")
  } else {
    print("Use custom list")
    custom <- countDataTrim[which(rownames(countDataTrim) %in% genelist$gene_id),]
  }
  
  sampleDataPrim$Group <- factor(sampleDataPrim$PrimCRLab)
  fname <- paste0('results/heatmaps_allmice/', fname, '.pdf')
  pdf(file = fname, width = 10, height = 10)
  # pheatmap
  if(is.null(genelist)){
    # for all expressed genes
    print("All Expressed Dimensions:")
    print(dim(countDataTrim))
    if(ncol(countDataTrim) > 2){
      pheatmap(countDataTrim, 
               scale="row", 
               show_rownames=F, cluster_cols = FALSE,
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  } else {
    # use genelist
    print("Using genelist:")
    print(dim(custom))
    if(ncol(custom) > 2){
      pheatmap(custom, 
               scale="row", 
               show_rownames=F, cluster_cols = FALSE,
               annotation_col= sampleDataPrim[c("Group","Mouse.ID")], 
               color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
               main=title)
    } else {
      print("Dimensions insufficient for clustering...")
    }
  }
  dev.off()
}
res <- createSampleHeatmap.allbiotypes(mat = 'all', mets = 'no', fname = 'allmice_without_mets_unclustered/allmice_sample_nomets_allbiotypes_heatmap')
countData.forHeatmap <- res[[1]]
sampleData.forHeatmap <- res[[2]]

# all expressed genes
limma <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt')
createGeneHeatmap.unclustered.custom(countDataPrim = countData.forHeatmap, 
                              sampleDataPrim = sampleData.forHeatmap, genelist = limma, 
                              fname = 'allmice_without_mets_unclustered/allmice_expressed_nomets_onlydegs_heatmap', 
                              title = "Heatmap All DE Genes")

# logFC > 1
limma.sub <- limma[which(abs(limma$logFC) > 1),]
createGeneHeatmap.unclustered.custom(countDataPrim = countData.forHeatmap, 
                                     sampleDataPrim = sampleData.forHeatmap, genelist = limma.sub, 
                                     fname = 'allmice_without_mets_unclustered/allmice_expressed_nomets_onlydegs_logfc1_heatmap', 
                                     title = "Heatmap DE Genes (Abs logFC > 1)")

# logFC > 1.5
limma.sub <- limma[which(abs(limma$logFC) > 1.5),]
createGeneHeatmap.unclustered.custom(countDataPrim = countData.forHeatmap, 
                                     sampleDataPrim = sampleData.forHeatmap, genelist = limma.sub, 
                                     fname = 'allmice_without_mets_unclustered/allmice_expressed_nomets_onlydegs_logfc15_heatmap', 
                                     title = "Heatmap DE Genes (Abs logFC > 1.5)")
