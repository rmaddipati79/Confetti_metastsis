###############################
# Author: Komal S Rathi
# Function: Molecular Subtype heatmaps
# PMet-high vs PMet-low
# Date: 10/08/2018
###############################

# load libraries
library(tidyverse)

# set working dir
setwd('~/Projects/PancModelUpdate/')

# load data
# molecular subtype signatures 
load('results/molecular_subtyping/molecular_subtypes.RData')

# merge all three signatures
if(identical(rownames(bailey.sig), rownames(collison.sig)) & identical(rownames(moffitt.sig), rownames(collison.sig))){
  sigs <- cbind(bailey = bailey.sig$class, collison = collison.sig$class, moffitt = moffitt.sig$class)
  rownames(sigs) <- rownames(bailey.sig)
}
rm(bailey.sig, collison.sig, moffitt.sig)

# load expression data and metadata
# expressed data for Pmet-high vs Pmet-low
load('data/collapsed_counts_to_human.RData')

# 13627 expressed protein coding genes
outputLimma_pmetvspnmet.volcano <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt')

# subset to protein coding expressed genes
human.genes <- human.genes %>% 
  filter(HGNC.symbol != "") %>%
  filter(MGI.symbol %in% outputLimma_pmetvspnmet.volcano$gene_symbol)
countDataPrim <- countDataPrim[rownames(countDataPrim) %in% human.genes$HGNC.symbol,]

# differentially expressed genes with adj Pval < 0.05 and abs logFC > 1
# 108 Up, 396 Down
diffexpr.pc <- outputLimma_pmetvspnmet.volcano %>% 
  filter(gene_type ==  "protein_coding" & 
           DEGAnnot != "Unchanged" & 
           adj.P.Val < 0.05 &
           abs(logFC) > 1)
diff.expr.human.genes <- human.genes %>% 
  filter(MGI.symbol %in% diffexpr.pc$gene_symbol) %>%
  .$HGNC.symbol
rm(diffexpr.pc, human.genes, outputLimma_pmetvspnmet.volcano)

# source code
source('code/molecular_subtyping/molecular_subtypes_helper.R')

create.heatmaps <- function(signature, heatmap.title, fname, mets, cluster.col){
  
  # sample correlation heatmap
  print("Sample Heatmap")
  createSampleHeatmap(countDataPrim, sampleDataPrim,
                      signature = signature, 
                      heatmap.title = "Sample Heatmap", 
                      fname = fname, mets = mets, cluster.col)
  
  # expressed protein coding genes heatmap
  print("Expressed genes")
  createGeneHeatmap(countDataPrim = countDataPrim, 
                    sampleDataPrim = sampleDataPrim, 
                    signature = signature,
                    num = NULL, fname = fname,
                    title = "Expressed Protein Coding", 
                    mets = mets, subset = NULL, cluster.col)
  
  # top 1000 most variable protein coding genes
  print("Top 1000 most variable genes")
  createGeneHeatmap(countDataPrim = countDataPrim, 
                    sampleDataPrim = sampleDataPrim, 
                    signature = signature,
                    num = 1000, fname = fname,
                    title = "Top 1000 Most Variable", mets = mets, subset = NULL,
                    cluster.col)
  
  # differentially expressed gene set heatmap
  # 396 up and 108 down (abs logFC > 1 and  p < 0.01)
  print("Diffexpr genes")
  createGeneHeatmap(countDataPrim = countDataPrim, 
                    sampleDataPrim = sampleDataPrim, 
                    signature = signature,
                    num = NULL, fname = fname,
                    title = "Diff. Expressed Genes\nAdj Pval < 0.05 and abs logFC > 1", 
                    mets = mets, subset = diff.expr.human.genes, 
                    cluster.col)
}

# without mets
# automatic clustering
create.heatmaps(signature = sigs, 
                fname = "results/molecular_subtyping/auto/molecularSubtypes_NoMets", 
                mets = "no", cluster.col = TRUE)
# manual clustering
create.heatmaps(signature = sigs, 
                fname = "results/molecular_subtyping/manual/molecularSubtypes_NoMets", 
                mets = "no", cluster.col = FALSE)

# with mets
# automatic clustering
create.heatmaps(signature = sigs, 
                fname = "results/molecular_subtyping/auto/molecularSubtypes_Mets", 
                mets = "yes", cluster.col = TRUE)

# manual clustering
create.heatmaps(signature = sigs, 
                fname = "results/molecular_subtyping/manual/molecularSubtypes_Mets", 
                mets = "yes", cluster.col = FALSE)

