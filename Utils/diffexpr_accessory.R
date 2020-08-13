# differential expression Primary Met vs Primary Non-mets
subMat.pmetvspnmet <- function(sampleDataPrim, countData, groups, contrast){
  sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab %in% groups),]
  countDataPrim <- countData[,rownames(sampleDataPrim)]
  
  print(dim(countDataPrim))
  print(dim(sampleDataPrim))
  
  # remove genes low expressing genes
  maxValueGenes <- apply(countDataPrim, FUN = max, MARGIN = 1)
  expressedGenes <- names(maxValueGenes[maxValueGenes > 100])
  countDataPrim <- countDataPrim[expressedGenes,]
  countDataPrim <- countDataPrim[,sort(colnames(countDataPrim))]
  
  # find all things to be merged
  sampleDataPrim$Primary_clone <- paste0(sampleDataPrim$Mouse.ID,'_',sampleDataPrim$Relationship.between.primary.tumors)
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
  
  # normalize using voom
  clonRel <- factor(sampleDataPrim$PrimCRLab)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  design <- model.matrix(~clonRel)
  y <- DGEList(counts = as.matrix(countDataPrim), genes = rownames(countDataPrim))
  y <- calcNormFactors(y)
  design <- model.matrix(~clonRel)
  v <- voom(y,design,plot=TRUE)
  voomData <- v$E
  voomDataSVA <- ComBat(voomData, batch=mouseID , mod=design, par.prior=T)
  
  # create contrasts and fit model
  design <- model.matrix(~0+clonRel)
  fit <- lmFit(voomDataSVA, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  outputLimma <- topTable(fit2, coef = 1, number=Inf)
  return(list(countDataPrim, sampleDataPrim, outputLimma))
}

# differential expression Mets vs Primary non-mets
subMat.metvspnmet <- function(sampleDataPrim, countData, groups, contrast){
  sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab %in% groups),]
  countDataPrim <- countData[,rownames(sampleDataPrim)]
  
  print(dim(countDataPrim))
  print(dim(sampleDataPrim))
  
  # remove genes low expressing genes
  maxValueGenes <- apply(countDataPrim, FUN = max, MARGIN = 1)
  expressedGenes <- names(maxValueGenes[maxValueGenes > 100])
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
  clonRel <- factor(sampleDataPrim$PrimCRLab)
  mouseID <- factor(sampleDataPrim$Mouse.ID)
  design <- model.matrix(~clonRel)
  y <- DGEList(counts = as.matrix(countDataPrim), genes = rownames(countDataPrim))
  y <- calcNormFactors(y)
  design <- model.matrix(~clonRel)
  v <- voom(y,design,plot=TRUE)
  voomData <- v$E
  voomDataSVA <- ComBat(voomData, batch=mouseID , mod=design, par.prior=T)
  
  # create contrasts and fit model
  design <- model.matrix(~0+clonRel)
  fit <- lmFit(voomDataSVA, design)
  contrast.matrix <- makeContrasts(contrasts = contrast, levels=design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  outputLimma <- topTable(fit2, coef = 1, number=Inf)
  return(list(countDataPrim, sampleDataPrim, outputLimma))
}

# add annotation
addAnnot <- function(x, geneAnnot){
  x <- merge(geneAnnot, x, by.x = 'gene_id', by.y = 'row.names')
  rownames(x) <- x$gene_id
  x$DEGene <- x$adj.P.Val < 0.05 
  x$DEGAnnot <- 'Unchanged'
  x$DEGAnnot[x$adj.P.Val < 0.05 & x$logFC > 0] <- "Up"
  x$DEGAnnot[x$adj.P.Val < 0.05 & x$logFC < 0] <- "Down"
  return(x)
}

# subset to protein coding genes
subset.PC <- function(x){
  x <- x[which(x$gene_type == "protein_coding"),]
  x <- x[grep('^Gm|Rik$|^[0-9]',x$gene_symbol,invert = TRUE),]
  return(x)
}