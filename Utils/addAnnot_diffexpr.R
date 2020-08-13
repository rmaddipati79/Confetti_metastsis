addAnnot <- function(x, geneAnnot, foldchange){
  
  subset.PC <- function(geneAnnot){
    geneAnnot <- geneAnnot[which(geneAnnot$gene_type %in% c('protein_coding')),]
    geneAnnot <- geneAnnot[grep('Rik$|Rik[0-9]{1}$', geneAnnot$gene_symbol, invert = TRUE),]
    geneAnnot <- geneAnnot[grep('Gm[0-9]{1,5}$', geneAnnot$gene_symbol, invert = TRUE),]
    geneAnnot <- geneAnnot[,c('gene_id','gene_type','gene_symbol','logFC','adj.P.Val','DEGene','DEGAnnot')]
    return(geneAnnot)
  }
  
  x <- merge(geneAnnot, x, by.x = 'gene_id', by.y = 'row.names')
  rownames(x) <- x$gene_id
  x$DEGene <- x$adj.P.Val < 0.05
  x$DEGAnnot <- 'Unchanged'
  x$DEGAnnot[x$adj.P.Val < 0.05 & x$logFC > (foldchange)] <- "Up"
  x$DEGAnnot[x$adj.P.Val < 0.05 & x$logFC < -(foldchange)] <- "Down"
  x <- x[which(x$DEGAnnot != "Unchanged"),]
  x <- subset.PC(x)
  return(x)
}