#####################################
# Author: Komal Rathi, Pichai Raman
# Function: Assign molecular subtypes
# Date: 09/18/2018
#####################################

# mySigUp <- ADEX
# exprs_pr <- exprs
# mySigDown <- NULL

# load libraries
library(pheatmap)
source('../Utils/pubTheme.R')

# default ggplot color function
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# z-score function
myZ <- function(x) {
  (x-mean(x))/sd(x)
}

# create signature based on list
createSigScore <- function(mySigUp, exprs_pr, sigName = "NoName", mySigDown = NULL){
  # Genes in intersection
  mySig <- union(mySigUp, mySigDown)
  intGenes <- intersect(rownames(exprs_pr), mySig)
  exprs_pr_cmat <- exprs_pr[intGenes,]
  exprs_pr_cmat_z <- data.frame(t(apply(exprs_pr_cmat, FUN=myZ, MARGIN=1)))
  
  # Pull out up and down genes
  discHitsGenesUp <- intersect(mySigUp, intGenes)
  discHitsGenesDown <- intersect(mySigDown, intGenes)
  
  # Generate Scores
  # myScores <- colSums(exprs_pr_cmat_z[discHitsGenesUp,]*(-1*log10(signatureList[discHitsGenesUp,"P.Value"])))-colSums(exprs_pr_cmat_z[discHitsGenesDown,]*(-1*log10(signatureList[discHitsGenesDown,"P.Value"])))
  myScores <- colSums(exprs_pr_cmat_z[discHitsGenesUp,])-colSums(exprs_pr_cmat_z[discHitsGenesDown,])
  aucDF <- data.frame(myScores)
  aucDF[,"sample"] <- rownames(aucDF)
  # classLabs <- data.frame(c(lowSamps, highSamps), c(rep(1, length(lowSamps)), rep(0, length(highSamps))))
  # rownames(classLabs) <- classLabs[,1];
  # aucDF <- aucDF[rownames(classLabs),]
  # aucDF <- cbind(aucDF[1], classLabs[2]);
  colnames(aucDF) <- c("score", "sample");
  aucDF[,"Signature"] <- sigName;
  return(list(aucDF, myScores));
}

mol.subtype <- function(signature, exprs){
  subtypes <- unique(signature$Subtype)
  for(i in 1:length(subtypes)){
    print(subtypes[i])
    dat.list <- as.character(signature[signature$Subtype == subtypes[i], 1])
    dat <- createSigScore(dat.list, exprs_pr = exprs, sigName= subtypes[i])[2]
    dat <- as.data.frame(dat)
    # normalize by length
    dat[,1] <- dat[,1]/length(dat.list)
    if(i == 1){
      total <- dat
    } else {
      total <- cbind(total, dat)
    }
  }
  colnames(total) <- subtypes
  # take max score
  total$class <- colnames(total)[max.col(total, ties.method = "first")]
  return(total)
}

# pca all mice together
create.pca.allmice <- function(signature, pca.title, file.name, dir.name, exprs){
  to.add <- sampleDataPrim[,c("PrimCRLab","Mouse.ID","Merge.Name")]
  to.add <- merge(signature, to.add, by = 'row.names')
  rownames(to.add) <- to.add$Row.names
  
  # PCA plot for all expressed genes
  prData <- prcomp(exprs)
  pca.data <- prData$rotation
  pca.data <- data.frame(pca.data)[1:3]
  
  to.add <- to.add[rownames(pca.data),]
  if(identical(rownames(pca.data), rownames(to.add))){
    pca.data <- cbind(pca.data, to.add)  
  }
  
  # create PCA
  system(paste0('mkdir -p ', dir.name))
  fname <- paste0(dir.name,'/', file.name, '_2Dpca.pdf')
  pdf(fname, width = 8, height = 10)
  p1 <- ggplot(data = pca.data, aes(PC1, PC2, group = PrimCRLab)) + 
    geom_point(aes(color = class, shape = PrimCRLab)) +
    geom_text(aes(label = Mouse.ID), size = 4, hjust = 0.5, vjust = 1.5) +
    theme_bw() 
  p2 <- ggplot(data = pca.data, aes(PC2, PC3, group = PrimCRLab)) + 
    geom_point(aes(color = class, shape = PrimCRLab)) +
    geom_text(aes(label = Mouse.ID), size = 4, hjust = 0.5, vjust = 1.5) +
    theme_bw() 
  gridExtra::grid.arrange(p1, p2, top = textGrob(pca.title, gp = gpar(fontsize=20,font=3)))
  dev.off()
  
  # Add Color
  n <- unique(pca.data$class)
  Mouse.Color <- gg_color_hue(length(n))
  cols <- data.frame(class = n, Mouse.Color)
  pca.data <- merge(pca.data, cols, by.x = 'class')
  
  # Add Shape
  pca.data$Mouse.Shape <- ifelse(pca.data$PrimCRLab == "Primary_Met", 15, ifelse(pca.data$PrimCRLab == "Primary_NoMet", 17, 16))
  
  # create 3D PCA
  fname <- paste0(dir.name, '/', file.name, '_3Dpca.pdf')
  pdf(fname, width = 10, height = 8)
  s3d <- scatterplot3d(pca.data[,"PC1"], pca.data[,"PC2"], pca.data[,"PC3"], 
                       xlab="PC1", ylab="PC2", zlab="PC3", 
                       color=as.character(pca.data[,"Mouse.Color"]), 
                       pch=pca.data[,"Mouse.Shape"], main=pca.title,
                       cex.symbols=1)
  tmpLegend <- unique(pca.data[,c("PrimCRLab", "Mouse.Color", "Mouse.Shape", "class")])
  tmpLegend[,"Mouse.ID"] <- paste(tmpLegend[,"PrimCRLab"], tmpLegend[,"class"], sep=" : ")
  legend(x="bottomright", pch=tmpLegend[,"Mouse.Shape"], col = as.character(tmpLegend[,"Mouse.Color"]), legend = tmpLegend[,"Mouse.ID"], cex = 0.5, inset=0)
  text(s3d$xyz.convert(pca.data[,"PC1"], pca.data[,"PC2"], pca.data[,"PC3"]), labels=pca.data$Mouse.ID, pos=1, cex = 0.5)
  dev.off()
}

# barplot of mice with subtypes
create.barplot <- function(signature, bar.title, file.name, dir.name){
  to.add <- sampleDataPrim[,c("PrimCRLab","Mouse.ID","Merge.Name")]
  to.add <- merge(signature, to.add, by = 'row.names')
  rownames(to.add) <- to.add$Row.names
  
  to.add <- dcast(to.add, PrimCRLab~class)
  to.add <- melt(to.add, variable.name = 'class')
  to.add <- to.add %>% group_by(PrimCRLab) %>%
    mutate(sum = sum(value))
  to.add$value <- signif(to.add$value/to.add$sum*100, digits = 2)
  to.add <- to.add[which(to.add$value != 0),]
  
  fname <- paste0(dir.name, '/', file.name, '_barplot.pdf')
  pdf(fname, width = 8, height = 8)
  p <- ggplot(to.add, aes(x = PrimCRLab, y = value, fill = factor(class), label = paste0(value,"%"))) +
    geom_bar(stat="identity", width = 0.7) +
    labs(x = "Signature", y = "Percent (%)", fill = "Category") +
    geom_text(size = 4, position = position_stack(vjust = 0.5)) +
    theme_bw() + theme_Publication() + ggtitle(bar.title)
  print(p)
  dev.off()
}

# pca for individual mice
create.pca.indmice <- function(signature, pca.title, file.name, dir.name, exprs){
  to.add.total <- sampleDataPrim[,c("PrimCRLab","Mouse.ID","Merge.Name")]
  to.add.total <- merge(signature, to.add.total, by = 'row.names')
  rownames(to.add.total) <- to.add.total$Row.names
  
  mouse.id <- unique(to.add.total$Mouse.ID)
  for(i in 1:length(mouse.id)){
    print(i)
    to.add <- to.add.total[which(to.add.total$Mouse.ID == mouse.id[i]),]
    exprs.sub <- exprs[,which(colnames(exprs) %in% to.add$Row.names)]
    
    # PCA plot for all expressed genes
    prData <- prcomp(exprs.sub)
    pca.data <- prData$rotation
    pca.data <- data.frame(pca.data)[1:3]
    
    to.add <- to.add[rownames(pca.data),]
    if(identical(rownames(pca.data), rownames(to.add))){
      pca.data <- cbind(pca.data, to.add)  
    }
    
    # create PCA
    fname <- paste0(dir.name, '/', file.name, '_', mouse.id[i], '_2Dpca.pdf')
    pdf(fname, width = 8, height = 10)
    p1 <- ggplot(data = pca.data, aes(PC1, PC2, group = PrimCRLab)) + 
      geom_point(aes(color = class, shape = PrimCRLab)) +
      geom_text(aes(label = Merge.Name), size = 4, hjust = 0.5, vjust = 1.5) +
      theme_bw() 
    p2 <- ggplot(data = pca.data, aes(PC2, PC3, group = PrimCRLab)) + 
      geom_point(aes(color = class, shape = PrimCRLab)) +
      geom_text(aes(label = Merge.Name), size = 4, hjust = 0.5, vjust = 1.5) +
      theme_bw() 
    gridExtra::grid.arrange(p1, p2, top = textGrob(pca.title, gp = gpar(fontsize=20,font=3)))
    dev.off()
    
    # Add Color
    n <- unique(pca.data$class)
    Mouse.Color <- gg_color_hue(length(n))
    cols <- data.frame(class = n, Mouse.Color)
    pca.data <- merge(pca.data, cols, by.x = 'class')
    
    # Add Shape
    pca.data$Mouse.Shape <- ifelse(pca.data$PrimCRLab == "Primary_Met", 15, ifelse(pca.data$PrimCRLab == "Primary_NoMet", 17, 16))
    
    # create 3D PCA
    fname <- paste0(dir.name, '/', file.name, '_', mouse.id[i], '_3Dpca.pdf')
    pdf(fname, width = 10, height = 8)
    s3d <- scatterplot3d(pca.data[,"PC1"], pca.data[,"PC2"], pca.data[,"PC3"], 
                         xlab="PC1", ylab="PC2", zlab="PC3", 
                         color=as.character(pca.data[,"Mouse.Color"]), 
                         pch=pca.data[,"Mouse.Shape"], main=pca.title,
                         cex.symbols=1)
    tmpLegend <- unique(pca.data[,c("PrimCRLab", "Mouse.Color", "Mouse.Shape", "class")])
    tmpLegend[,"Mouse.ID"] <- paste(tmpLegend[,"PrimCRLab"], tmpLegend[,"class"], sep=" : ")
    legend(x="bottomright", pch=tmpLegend[,"Mouse.Shape"], col = as.character(tmpLegend[,"Mouse.Color"]), legend = tmpLegend[,"Mouse.ID"], cex = 0.7, inset=0)
    text(s3d$xyz.convert(pca.data[,"PC1"], pca.data[,"PC2"], pca.data[,"PC3"]), labels=pca.data$Merge.Name, pos=1, cex = 0.7)
    dev.off()
  }
}

# heatmap all mice for all three classifications 
# heatmap primary mets and primary non-mets for all three classifications
# revised version
createSampleHeatmap <- function(countDataPrim, sampleDataPrim, signature, mets, fname, heatmap.title, cluster.col = TRUE){
  
  if(mets == "no") {
    print("removing mets")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # merge signature with sample data
  to.add <- merge(signature, sampleDataPrim, by = 'row.names')
  rownames(to.add) <- to.add$Row.names
  to.add <- to.add[colnames(countDataPrim),]
  if(identical(rownames(to.add), colnames(countDataPrim))) {
    rownames(to.add) <- to.add$Merge.Name
    colnames(countDataPrim) <- rownames(to.add)
  }
  
  if(cluster.col == FALSE){
    print("Reorder columns by PrimCRLab")
    to.add <- to.add[order(to.add$PrimCRLab),]
    countDataPrim <- countDataPrim[,rownames(to.add)]
  }
  
  # create correlation heatmap
  fname <- paste0(fname,'.pdf')
  pdf(fname, width = 14, height = 12)
  myCor <- cor(countDataPrim)
  myCorNames <- rownames(myCor)
  pheatmap(t(myCor),
           cluster_cols = cluster.col,
           cluster_rows = cluster.col,
           annotation_col = to.add[c("PrimCRLab","Mouse.ID","bailey","collison","moffitt")], 
           main = "Sample Correlation Heatmap\n", 
           display_numbers = FALSE)
  dev.off()
}

# revised version
# n is for number of variable genes
# matrix and sample sheet come from the function above
createGeneHeatmap <- function(countDataPrim, sampleDataPrim, signature, fname, title, num = NULL, mets, subset = NULL, cluster.col = TRUE){
  
  # expressed genes only
  if(is.null(subset) & is.null(num)){
    fname <- paste0(fname,'_13627_Expressed.pdf')
  }
  
  # subset using gene list
  if(!is.null(subset)){
    print("Subset using provided genelist")
    fname <- paste0(fname,'_504_DiffExpr.pdf')
    countDataPrim <- countDataPrim[rownames(countDataPrim) %in% subset,]
  }
  
  # remove mets
  if(mets == "no") {
    print("removing mets")
    sampleDataPrim <- sampleDataPrim[which(sampleDataPrim$PrimCRLab != "Met"),]
  }
  countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  
  # now subset depending on N
  if(!is.null(num)) {
    print("Use n most-variable protein coding genes...")
    fname <- paste0(fname,'_', num,'_MostVariable.pdf')
    # calculate the n most variable genes
    var_genes <- apply(countDataPrim, 1, var)
    select_var <- names(sort(var_genes, decreasing=TRUE))[1:num]
    countDataPrim <- countDataPrim[select_var,]
  }
  
  # merge with signature to add Class
  sampleDataPrim <- merge(signature, sampleDataPrim, by = 'row.names')
  rownames(sampleDataPrim) <- sampleDataPrim$Row.names
  sampleDataPrim <- sampleDataPrim[colnames(countDataPrim),]
  if(identical(rownames(sampleDataPrim), colnames(countDataPrim))) {
    rownames(sampleDataPrim) <- sampleDataPrim$Merge.Name
    colnames(countDataPrim) <- rownames(sampleDataPrim)
  }
  
  if(cluster.col == FALSE){
    print("Reorder columns by PrimCRLab")
    sampleDataPrim <- sampleDataPrim[order(sampleDataPrim$PrimCRLab),]
    countDataPrim <- countDataPrim[,rownames(sampleDataPrim)]
  }
  
  # create heatmap
  if(ncol(countDataPrim) > 2){
    pdf(fname, width = 12, height = 12)
    p <- pheatmap(countDataPrim,
             cluster_cols = cluster.col,
             scale = "row", 
             show_rownames = F, 
             annotation_col = sampleDataPrim[c("PrimCRLab","Mouse.ID","bailey","collison","moffitt")], 
             color = colorRampPalette(c("blue", "white", "firebrick3"))(50), 
             main = title)
    p
    dev.off()
  } else {
    print("Dimensions insufficient for clustering...")
  }
}
