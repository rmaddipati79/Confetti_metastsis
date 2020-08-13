###############################
# Author: Komal S Rathi
# Function: Limma Differential Analysis
# Date: 10/08/2018
###############################

library(edgeR)
library(sva)
library(limma)
library(ggplot2)
library(scales)
library(ggpubr)
library(tidyverse)
source("code/Utils/pubTheme.R")
source("code/Utils/diffexpr_accessory.R")
source("code/Utils/volcanoPlotAccessory.R")

# differential expression
# load matrix (mice = 7)
load('data/all_mouse_matrix_072418.RData')

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

# diffexpr results
res <- subMat.pmetvspnmet(sampleDataPrim = sampleDataPrim, 
                          countData = countData, 
                          groups = c('Primary_Met', 'Primary_NoMet'), 
                          contrast = 'clonRelPrimary_Met-clonRelPrimary_NoMet')
save(res, file = 'results/Limma/PrimaryMet_vs_PrimaryNonMet/pmetvspnmet.RData')

countData_pmetvspnmet <- res[[1]]
sampleData_pmetvspnmet <- res[[2]]
outputLimma_pmetvspnmet <- res[[3]]

# add gene annotation to diffexpr results
outputLimma_pmetvspnmet <- addAnnot(outputLimma_pmetvspnmet, geneAnnot)
write.table(outputLimma_pmetvspnmet, file = 'results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes.txt', quote = F, sep = "\t", row.names = F)

# differential expression results output
diffexpr <- outputLimma_pmetvspnmet %>%
  filter(adj.P.Val < 0.05)
write.table(diffexpr, file = 'results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt', quote = F, sep = "\t", row.names = F)

# subset to protein coding genes and remove riken genes
outputLimma_pmetvspnmet.volcano <- subset.PC(outputLimma_pmetvspnmet)
write.table(outputLimma_pmetvspnmet.volcano, file = 'results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt', quote = F, sep = "\t", row.names = F)

# upregulated genes
upreg <- outputLimma_pmetvspnmet.volcano %>%
  filter(DEGAnnot == "Up")
write.table(upreg, 
            file = 'results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_Upreg.txt', 
            quote = F, sep = "\t", row.names = F)

#########################
# Volcano plot
#########################
plotVolcano(outputLimma_pmetvspnmet.volcano, 
            title = "TMH vs. TML (Abs LogFC >1, p<0.01)", 
            otherCol = c('firebrick3','gray','forestgreen'), lfcutoff = 1, pvalcutoff = 0.01)
ggsave("results/Volcano/PrimaryMetvsPrimaryNonMetClones_VolcanoPlot_allMice_pc_logfc1_pval001.eps", width=15, height=10)

plotVolcano(outputLimma_pmetvspnmet.volcano, 
            title = "TMH vs. TML (Abs LogFC >1, p<0.05)", 
            otherCol = c('firebrick3','gray','forestgreen'), lfcutoff = 1, pvalcutoff = 0.05)
ggsave("results/Volcano/PrimaryMetvsPrimaryNonMetClones_VolcanoPlot_allMice_pc_logfc1_pval005.eps", width=15, height=10)

plotVolcanoText(outputLimma_pmetvspnmet.volcano, 
            title = "TMH vs. TML (Abs LogFC >1, p<0.01)", 
            otherCol = c('firebrick3','gray','forestgreen'), lfcutoff = 1, pvalcutoff = 0.01)
ggsave("results/Volcano/PrimaryMetvsPrimaryNonMetClones_VolcanoPlot_allMice_pc_logfc1_pval001_genesym.eps", width=15, height=10)

plotVolcanoText(outputLimma_pmetvspnmet.volcano, 
                title = "TMH vs. TML (Abs LogFC >1, p<0.05)", 
                otherCol = c('firebrick3','gray','forestgreen'), lfcutoff = 1, pvalcutoff = 0.05)
ggsave("results/Volcano/PrimaryMetvsPrimaryNonMetClones_VolcanoPlot_allMice_pc_logfc1_pval005_genesym.eps", width=15, height=10)

# differential expression
res <- subMat.metvspnmet(sampleDataPrim = sampleDataPrim, 
                         countData = countData, 
                         groups = c('Met', 'Primary_NoMet'), 
                         contrast = 'clonRelMet-clonRelPrimary_NoMet')
countData_metvspnmet <- res[[1]]
sampleData_metvspnmet <- res[[2]]
outputLimma_metvspnmet <- res[[3]]
outputLimma_metvspnmet <- addAnnot(outputLimma_metvspnmet, geneAnnot)
outputLimma_metvspnmet.volcano <- subset.PC(outputLimma_metvspnmet)
write.table(outputLimma_metvspnmet.volcano, file = 'results/Limma/Met_vs_PrimaryNonMet/MetvsPrimaryNonMet_Limma_allMice.txt', quote = F, sep = "\t", row.names = F)

#########################
# Scatter plot
#########################
one <- outputLimma_pmetvspnmet.volcano
two <- outputLimma_metvspnmet.volcano

one <- one[,c('gene_symbol','logFC','DEGene')]
one$Label <- ifelse(one$DEGene == TRUE & one$logFC < 0,'Downreg', 'Upreg')
one$Label <- ifelse(one$DEGene == FALSE,'NotSig',one$Label)
colnames(one)[2:ncol(one)] <- paste0('one_',colnames(one)[2:ncol(one)])
one <- unique(one)

two <- two[,c('gene_symbol','logFC','DEGene')]
two$Label <- ifelse(two$DEGene == TRUE & two$logFC < 0,'Downreg', 'Upreg')
two$Label <- ifelse(two$DEGene == FALSE,'NotSig',two$Label)
colnames(two)[2:ncol(two)] <- paste0('two_',colnames(two)[2:ncol(two)])
two <- unique(two)
total <- merge(one, two, by = 'gene_symbol')

total$NewLabel <- NA
total$NewLabel[total$one_Label == "Upreg" & total$two_Label == "Upreg"] <- 'Upregulated'
total$NewLabel[total$one_Label == "Downreg" & total$two_Label == "Downreg"] <- 'Downregulated'
total$NewLabel[is.na(total$NewLabel)] <- 'Others'

summ <- summary(lm(formula = one_logFC ~ two_logFC, data = total))
summ

nums <- plyr::count(total$NewLabel)
total <- merge(total, nums, by.x = 'NewLabel', by.y = 'x')
total$NewLabel <- paste0(total$NewLabel," (n = ", total$freq, ")")
ggplot(total, aes(one_logFC, two_logFC)) +
  geom_point(size = 0.5, aes(color = NewLabel)) +
  xlab('Primary Mets vs Primary Non-Mets (logFC)') + 
  ylab('Mets vs Primary Non-Mets (logFC)') + theme_bw() +
  theme_Publication() + ggtitle(paste0('Protein-coding genes (n = ',nrow(total),')')) +
  guides(color=guide_legend(title="Status")) + geom_smooth(method = "lm", se = FALSE, lwd = 0.5) +
  annotate("text", x=2.5, y=-5, label = "R^2 = 0.20\nP-value < 2.2e16") + 
  scale_color_manual(values=c("#F8766D", "gray", "#00BA38"))
ggsave("results/Scatter/Scatter_PrimaryMetvsNonMets_MetsvsPrimaryNonMets_allMice.eps", width=8, height=8)

