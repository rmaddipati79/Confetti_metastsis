###################################
# Author: Komal S Rathi
# Function: Create FPKM matrix 
# Date: 03/21/2019
###################################

# get and merge all FPKM data 
dat <- readRDS('data/all_mice_FPKM.RDS')
colnames(dat) <- gsub('-','_',colnames(dat))
samplesheet <- read.delim('data/ALL MICE Master samplesheet 072418.txt')
samplesheet$Column <- gsub('-','_',samplesheet$Column)
dat$gene_id <- gsub('[.].*','', dat$gene_id)
rownames(dat) <- dat$gene_id 
dat$gene_id <- NULL
dat$gene_symbol <- NULL
dat <- dat[,which(colnames(dat) %in% samplesheet$Column)]
rownames(samplesheet) <- samplesheet$Non.Merge.Name
samplesheet <- samplesheet[order(samplesheet$Column),]
dat <- dat[,order(colnames(dat))]

if(identical(samplesheet$Column, colnames(dat))){
  # rm(countData2, sampleData2)
  # colnames(countData) <- sampleData$Non.Merge.Name
  print("Yes")
  colnames(dat) <- samplesheet$Non.Merge.Name
}

load('data/geneAnnot.RData')
sampleData <- samplesheet
countData <- dat 
save(countData, sampleData, file = 'data/all_mouse_matrix_FPKM_072418.RData')

# merge clones
# we will probably get rid of this
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
expressedGenes <- names(maxValueGenes[maxValueGenes > 10])
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

countData <- countDataPrim
sampleData <- sampleDataPrim
load('data/geneAnnot.RData')
save(countData, sampleData, geneAnnot, file = 'data/all_mouse_matrix_FPKM.RData')
