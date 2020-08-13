# accessory for volcano plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

# volcano plot, takes in limma analysis
plotVolcano <- function(result, title = "Volcano Plot", otherCol = c("gray","blue"), lfcutoff = NULL, pvalcutoff = 0.05) {
  if(!is.null(lfcutoff)){
    result$DEGAnnot <- 'Other'
    result$DEGAnnot[result$adj.P.Val < pvalcutoff & result$logFC > lfcutoff] <- "Up"
    result$DEGAnnot[result$adj.P.Val < pvalcutoff & result$logFC < -(lfcutoff)] <- "Down"
  }
  nums <- plyr::count(result$DEGAnnot)
  result <- merge(result, nums, by.x = 'DEGAnnot', by.y = 'x')
  result$DEGAnnot <- paste0(result$DEGAnnot," (n = ", result$freq, ")")
  p <- ggplot(result, aes(x = logFC, y = adj.P.Val, color = DEGAnnot)) +
    geom_point() + scale_y_continuous(trans = reverselog_trans(10)) + 
    ggtitle(title) + scale_colour_manual(values = otherCol)
  p <- p + theme_Publication() + labs(color = "DEG")
  return(p)
}

# volcano plot text, takes in limma analysis
plotVolcanoText <- function(result, title = "Volcano Plot", otherCol = c("gray","blue"), lfcutoff = NULL, syms = NULL, pvalcutoff = 0.05){
  if(!is.null(lfcutoff)){
    result$DEGAnnot <- 'Other'
    result$DEGAnnot[result$adj.P.Val < pvalcutoff & result$logFC > lfcutoff] <- "Up"
    result$DEGAnnot[result$adj.P.Val < pvalcutoff & result$logFC < -(lfcutoff)] <- "Down"
  }
  if(!is.null(syms)){
    result$gene_label <- ifelse(result$gene_symbol %in% syms, result$gene_symbol, '')
  } else {
    result$gene_label <- result$gene_symbol
  }
  nums <- plyr::count(result$DEGAnnot)
  result <- merge(result, nums, by.x = 'DEGAnnot', by.y = 'x')
  result$DEGAnnot <- paste0(result$DEGAnnot," (n = ", result$freq, ")")
  if(is.null(syms)){
    p <- ggplot(result, aes(x = logFC, y = adj.P.Val, label = gene_label)) + 
      geom_point(aes(color = DEGAnnot), shape = 1) + 
      geom_text(aes(color = DEGAnnot), vjust = 0, nudge_y = 0.05) + 
      scale_y_continuous(trans = reverselog_trans(10)) + 
      ggtitle(title) + scale_colour_manual(values = otherCol) + 
      theme_Publication()
  } else {
    p <- ggplot(result, aes(x = logFC, y = adj.P.Val, label = gene_label)) + 
      geom_point(aes(color = DEGAnnot), shape = 1) + 
      geom_text(vjust = 0, nudge_y = 0.05) + 
      scale_y_continuous(trans = reverselog_trans(10)) + 
      ggtitle(title) + scale_colour_manual(values = otherCol) + 
      theme_Publication() 
  }
  return(p)
}

