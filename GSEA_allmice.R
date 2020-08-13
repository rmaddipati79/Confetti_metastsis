#####################################
# Author: Komal Rathi, Pichai Raman
# Function: GSEA 
# Date: 04/03/2019
#####################################

# GSEA
# GSEA on countData_PrimarymetvsPrimarynonmets (now countData)
# sampleData_PrimarymetvsPrimarynonmets (now sampleData)
library(gage)
library(MSigDB)
library(dplyr)

load('data/all_mouse_matrix_072418.RData')
load('results/Limma/PrimaryMet_vs_PrimaryNonMet/pmetvspnmet.RData')
countData <- res[[1]]
sampleData <- res[[2]]
countData <- log2(countData+1)
ref.idx = grep('Primary_NoMet', sampleData$PrimCRLab)
samp.idx = grep('Primary_Met', sampleData$PrimCRLab)

# prepare data to match MSIGDB datasets
if(symbol == 'yes'){
  dat <- merge(countData, geneAnnot, by.x = 'row.names', by.y = 'gene_id')
  convertMouseGeneList <- function(x){
    require("biomaRt")
    human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    return(genesV2)
  }
  human.genes <- convertMouseGeneList(dat$gene_symbol)
  dat <- merge(dat, human.genes, by.x = 'gene_symbol', by.y = 'MGI.symbol')
  dat$gene_symbol <- NULL
  dat$gene_type <- NULL
  dat$Row.names <- NULL
  dat <- dat %>% mutate(means = rowMeans(.[1:nrow(sampleData)])) %>% 
    arrange(desc(means)) %>% 
    distinct(HGNC.symbol, .keep_all = TRUE) %>% ungroup()
  rownames(dat) <- dat$HGNC.symbol
  dat$HGNC.symbol <- NULL
  dat$means <- NULL
} else {
  dat <- countData
}

# create input files for GSEA web version
# expression file
s <- sampleData[order(sampleData$PrimCRLab, decreasing = T),]
gct <- dat[,rownames(s)]
add <- data.frame(NAME = c("#1.2",nrow(gct),"NAME"), Description = c('', ncol(gct), "Description"))
add[,3:22] <- ''
colnames(add)[3:22] <- colnames(gct)
add[3,3:22] <- colnames(gct)
annot <- data.frame(NAME = rownames(gct), Description = 'na')
annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
add <- rbind(add, annot)
if(symbol == 'yes'){
  write.table(add, file = 'results/Enrichment/GSEA_web_version/input_files/input_mat.gct', quote = F, sep = "\t", col.names = F, row.names = F)
} else {
  write.table(add, file = 'results/Enrichment/GSEA_web_version/input_files/input_mat_ensemblids.gct', quote = F, sep = "\t", col.names = F, row.names = F)
}

# phenotype file
ph <- data.frame(A = c(ncol(gct),"#Primary_Met",''), B = c('2', 'Primary_NoMet',''), C = c(1,'',''), stringsAsFactors = F)
ph[,4:20] <- ''
ph[3,] <- s$PrimCRLab
write.table(ph, file = 'results/Enrichment/GSEA_web_version/input_files/phenotype.cls', quote = F, sep = "\t", col.names = F, row.names = F)

# chp file
chp <- data.frame('Probe Set ID' = rownames(gct), 'Gene Symbol' = rownames(gct), 'Gene Title' = 'na', check.names = F)
write.table(chp, file = 'results/Enrichment/GSEA_web_version/input_files/mapping.chip', quote = F, sep = "\t", row.names = F)

# what databases are available
names(MSigDB)

# hallmark genes
hallmark <- gage(exprs = dat, 
                 gsets = MSigDB$HALLMARK, 
                 ref = ref.idx, samp = samp.idx, compare = "unpaired")
hallmark.greater <- as.data.frame(hallmark$greater[,1:5])
hallmark.lesser <- as.data.frame(hallmark$less[,1:5])
write.table(hallmark.greater, file = 'results/Enrichment/GSEA/hallmark_upreg_allmice.txt', quote = F, sep = "\t")
write.table(hallmark.lesser, file = 'results/Enrichment/GSEA/hallmark_downreg_allmice.txt', quote = F, sep = "\t")

# c1 positional
c1.pos <- gage(exprs = dat, 
               gsets = MSigDB$C1_POSITIONAL, 
               ref = ref.idx, samp = samp.idx, compare = "unpaired")
head(c1.pos$greater)

# c2 curated
c2.curated <- gage(exprs = dat, 
                   gsets = MSigDB$C2_CURATED, 
                   ref = ref.idx, samp = samp.idx, compare = "unpaired")
head(c2.curated$greater)

# c3 motif
c3.motif <- gage(exprs = dat, 
                 gsets = MSigDB$C3_MOTIF, 
                 ref = ref.idx, samp = samp.idx, compare = "unpaired")
c3.motif.greater <- as.data.frame(c3.motif$greater[,1:5])
c3.motif.lesser <- as.data.frame(c3.motif$less[,1:5])
write.table(c3.motif.greater, file = 'results/Enrichment/GSEA/c3motif_upreg_allmice.txt', quote = F, sep = "\t")
write.table(c3.motif.lesser, file = 'results/Enrichment/GSEA/c3motif_downreg_allmice.txt', quote = F, sep = "\t")

# c4 computational
c4.comp <- gage(exprs = dat, 
                gsets = MSigDB$C4_COMPUTATIONAL, 
                ref = ref.idx, samp = samp.idx, compare = "unpaired")
head(c4.comp$greater)

# c5 GO
c5.go <- gage(exprs = dat, 
              gsets = MSigDB$C5_GENE_ONTOLOGY, 
              ref = ref.idx, samp = samp.idx, compare = "unpaired")
c5.go.greater <- as.data.frame(c5.go$greater[,1:5])
c5.go.lesser <- as.data.frame(c5.go$less[,1:5])
write.table(c5.go.greater, file = 'results/Enrichment/GSEA/c5go_upreg_allmice.txt', quote = F, sep = "\t")
write.table(c5.go.lesser, file = 'results/Enrichment/GSEA/c5go_downreg_allmice.txt', quote = F, sep = "\t")

# c6 oncogenic signatures
c6.onco <- gage(exprs = dat, 
                gsets = MSigDB$C6_ONCOGENIC_SIGNATURES, 
                ref = ref.idx, samp = samp.idx, compare = "unpaired")
head(c6.onco$greater)

# c7 immunologic signatures
c7.immuno <- gage(exprs = dat, 
                  gsets = MSigDB$C7_IMMUNOLOGIC_SIGNATURES, 
                  ref = ref.idx, samp = samp.idx, compare = "unpaired")
head(c7.immuno$greater)

# convert mouse gene symbol to human entrez ID
dat <- merge(countData, geneAnnot, by.x = 'row.names', by.y = 'gene_id')
convertMouseGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("entrezgene"), martL = human, uniqueRows=T)
  return(genesV2)
}
human.genes <- convertMouseGeneList(dat$gene_symbol)
dat <- merge(dat, human.genes, by.x = 'gene_symbol', by.y = 'MGI.symbol')
dat$gene_symbol <- NULL
dat$Row.names <- NULL
dat$gene_type <- NULL
dat <- dat %>% mutate(means = rowMeans(.[1:nrow(sampleData)])) %>% arrange(desc(means)) %>% distinct(NCBI.gene.ID, .keep_all = TRUE) %>% ungroup()
rownames(dat) <- dat$NCBI.gene.ID
dat$NCBI.gene.ID <- NULL
dat$means <- NULL

# kegg
data(kegg.gs)
kegg <- gage(exprs = dat, 
     gsets = kegg.gs, 
     ref = ref.idx, samp = samp.idx, compare = "unpaired")
kegg.greater <- as.data.frame(kegg$greater[,1:5])
kegg.lesser <- as.data.frame(kegg$less[,1:5])
write.table(kegg.greater, file = 'results/Enrichment/GSEA/kegg_upreg_allmice.txt', quote = F, sep = "\t")
write.table(kegg.lesser, file = 'results/Enrichment/GSEA/kegg_downreg_allmice.txt', quote = F, sep = "\t")

# kegg.gs.dise
data("kegg.gs.dise")
kegg.dise <- gage(exprs = dat, 
             gsets = kegg.gs.dise, 
             ref = ref.idx, samp = samp.idx, compare = "unpaired")
kegg.dise.greater <- as.data.frame(kegg.dise$greater[,1:5])
kegg.dise.lesser <- as.data.frame(kegg.dise$less[,1:5])
write.table(kegg.greater, file = 'results/Enrichment/GSEA/kegg_upreg_allmice.txt', quote = F, sep = "\t")
write.table(kegg.lesser, file = 'results/Enrichment/GSEA/kegg_downreg_allmice.txt', quote = F, sep = "\t")

# go.gs
data("go.gs")
go <- gage(exprs = dat, 
                  gsets = go.gs, 
                  ref = ref.idx, samp = samp.idx, compare = "unpaired")
go.greater <- as.data.frame(go$greater[,1:5])
go.lesser <- as.data.frame(go$less[,1:5])
write.table(go.greater, file = 'results/Enrichment/GSEA/gogs_upreg_allmice.txt', quote = F, sep = "\t")
write.table(go.lesser, file = 'results/Enrichment/GSEA/gogs_downreg_allmice.txt', quote = F, sep = "\t")

# gageData
# now, we will use gene sets that have mouse entrez ids 
dat <- merge(countData, geneAnnot, by.x = 'row.names', by.y = 'gene_id')
getMouseID <- function(x){
  require("biomaRt")
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 <- getBM(attributes = c("mgi_symbol","entrezgene"), filters = "mgi_symbol", values = x, mart = mouse)
  return(genesV2)
}
geneids <- getMouseID(dat$gene_symbol)
dat <- merge(dat, geneids, by.x = 'gene_symbol', by.y = 'mgi_symbol')
dat$gene_symbol <- NULL
dat$Row.names <- NULL
dat$gene_type <- NULL
dat <- dat[!is.na(dat$entrezgene),]
dat <- dat %>% mutate(means = rowMeans(.[1:nrow(sampleData)])) %>% arrange(desc(means)) %>% distinct(entrezgene, .keep_all = TRUE) %>% ungroup()
rownames(dat) <- dat$entrezgene
dat$entrezgene <- NULL
dat$means <- NULL


library(gageData)
data("go.sets.mm")
data("go.subs.mm")
data("kegg.sets.mm")
go.mf.mm = go.sets.mm[go.subs.mm$MF]
go.bp.mm = go.sets.mm[go.subs.mm$BP]
go.cc.mm = go.sets.mm[go.subs.mm$CC]

# molecular function
go.mf <- gage(exprs = dat, 
           gsets = go.mf.mm, 
           ref = ref.idx, samp = samp.idx, compare = "unpaired")
go.mf.greater <- as.data.frame(go.mf$greater[,1:5])
go.mf.lesser <- as.data.frame(go.mf$less[,1:5])
write.table(go.mf.greater, file = 'results/Enrichment/GSEA/gomf_upreg_allmice.txt', quote = F, sep = "\t")

# biological process
go.bp <- gage(exprs = dat, 
              gsets = go.bp.mm, 
              ref = ref.idx, samp = samp.idx, compare = "unpaired")
go.bp.greater <- as.data.frame(go.bp$greater[,1:5])
go.bp.lesser <- as.data.frame(go.bp$less[,1:5])
write.table(go.bp.greater, file = 'results/Enrichment/GSEA/gobp_upreg_allmice.txt', quote = F, sep = "\t")

# cellular component
go.cc <- gage(exprs = dat, 
              gsets = go.cc.mm, 
              ref = ref.idx, samp = samp.idx, compare = "unpaired")
go.cc.greater <- as.data.frame(go.cc$greater[,1:5])
go.cc.lesser <- as.data.frame(go.cc$less[,1:5])
write.table(go.cc.greater, file = 'results/Enrichment/GSEA/gocc_upreg_allmice.txt', quote = F, sep = "\t")

# kegg
kegg.mm <- gage(exprs = dat, 
              gsets = kegg.sets.mm, 
              ref = ref.idx, samp = samp.idx, compare = "unpaired")
kegg.mm.greater <- as.data.frame(kegg.mm$greater[,1:5])
kegg.mm.lesser <- as.data.frame(kegg.mm$less[,1:5])
write.table(kegg.mm.greater, file = 'results/Enrichment/GSEA/keggMouse_upreg_allmice.txt', quote = F, sep = "\t")
write.table(kegg.mm.lesser, file = 'results/Enrichment/GSEA/keggMouse_downreg_allmice.txt', quote = F, sep = "\t")
save.image(file = 'results/Enrichment/GSEA/gsea.RData')
