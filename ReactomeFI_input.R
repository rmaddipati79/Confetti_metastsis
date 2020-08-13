# ReactomeFI input
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

convertMouseGeneList <- function(x){
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}
df <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes.txt', stringsAsFactors = F)


# 1. Top 50 DEG UP  
up50 <- unique(df[which(df$DEGAnnot == 'Up'),])
up50 <- up50[order(up50$adj.P.Val),]
up50 <- up50[1:50,'gene_symbol']
up50 <- convertMouseGeneList(up50)
write.table(up50$HGNC.symbol, file = 'results/ReactomeFI/top50_up.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# 2. Top 50 DEG Down (nothing below 0.29 FDR)
down50 <- unique(df[which(df$DEGAnnot == 'Down'),])
down50 <- down50[order(down50$adj.P.Val),]
down50 <- down50[1:50,'gene_symbol']
down50 <- convertMouseGeneList(down50)
write.table(down50$HGNC.symbol, file = 'results/ReactomeFI/top50_down.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# 3. genes with LogFC > 1  
upreg <- unique(df[which(df$DEGAnnot == 'Up' & df$logFC > 1),'gene_symbol'])
upreg <- convertMouseGeneList(upreg)
write.table(upreg$HGNC.symbol, file = 'results/ReactomeFI/upreg_logFC1.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# 4. genes with LogFC < -1
downreg <- unique(df[which(df$DEGAnnot == 'Down' & df$logFC < -1),'gene_symbol'])
downreg <- convertMouseGeneList(downreg)
write.table(downreg$HGNC.symbol, file = 'results/ReactomeFI/downreg_logFC1.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# 5. all down DEGs (nothing below 0.29 FDR)
alldown <- unique(df[which(df$DEGAnnot == 'Down'),'gene_symbol'])
alldown <- convertMouseGeneList(alldown)
write.table(alldown$HGNC.symbol, file = 'results/ReactomeFI/downreg.txt', quote = F, sep = "\t", row.names = F, col.names = F)

# 6. all up DEGs
allup <- unique(df[which(df$DEGAnnot == 'Up'),'gene_symbol'])
allup <- convertMouseGeneList(allup)
write.table(allup$HGNC.symbol, file = 'results/ReactomeFI/upreg.txt', quote = F, sep = "\t", row.names = F, col.names = F)
