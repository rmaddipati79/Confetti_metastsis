boxplot.multiple.groups <- function(myMat, mySampData, myGeneID, myGeneSym, mets = 'yes', test){

  myTitle <- paste0(myGeneSym,' Gene Expression\n')
  tmpData <- t(myMat[which(rownames(myMat) %in% myGeneID),,drop = F])
  tmpDF <- merge(mySampData, tmpData, by = 'row.names')
  colnames(tmpDF)[ncol(tmpDF)] <- "Gene"
  if(mets == 'no'){
    tmpDF <- tmpDF[which(tmpDF$PrimCRLab != 'Met'),]
    my_comparisons <- list( c("Primary_NoMet", "Primary_Met"))
  } else {
    my_comparisons <- list( c("Primary_NoMet", "Primary_Met"), c("Primary_Met", "Met"), c("Primary_NoMet", "Met") )
  }
  p <- ggplot(tmpDF, aes(factor(PrimCRLab), Gene)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot(lwd = 0.5, fatten = 0.7, outlier.shape = 1, width = 0.5, outlier.size = 1) +
    theme_bw() + theme_Publication() +
    xlab("\nSample Type") + ylab("Log2 FPKM\n") + ggtitle(myTitle) +
    stat_compare_means(comparisons = my_comparisons, method = test)
  return(p)
}
