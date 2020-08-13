#################################################
# Code to check MYC Association with survival in PDAC
# Author: Pichai Raman
# Date: 10/10/2018
#################################################

library(tidyverse)
library(GGally)
library(biomaRt)
library(survminer)
library(RDiseaseXpress)

#Source Code
source("code/surv_analysis/PR_Utils.R")
source("code/surv_analysis/kapmGroup.R")
source("code/surv_analysis/KaplanScan.R")


# Read TCGA Data
expData <- read.delim("data/surv_analysis/data_RNA_Seq_v2_expression_median.txt", check.names = F)
clinData <- read.delim("data/surv_analysis/paad_tcga_clinical_data.txt")  # clinical data
rownames(clinData) <- clinData$Sample.ID # n = 176


# format expression data 
expData <- expData[-2]
expData <- na.omit(expData)
expData[,"maxVal"] <- rowMeans(expData[2:ncol(expData)])
expData <- expData[order(-expData[,"maxVal"]),]
expData <- expData[!duplicated(expData[,1]),]
rownames(expData) <- expData[,1]
expData <- expData[-1]
expData <- expData[-ncol(expData)]

# first use as is
create.surv.plots <- function(clinData, expData, suffix){
  # samples with expression  data
  intSamples <- intersect(rownames(clinData), colnames(expData))
  expData <- expData[,intSamples]
  clinData <- clinData[intSamples,]
  
  print(nrow(clinData))
  fname <- file.path("results/survival_analysis/", 
                     paste0("Samples_used_", suffix, ".txt"))
  write.table(x = clinData$Sample.ID, file = fname, row.names = F, quote = F, col.names = F)
  
  # overall survival
  clinData$Overall.Survival.Status <- ifelse(clinData$Overall.Survival.Status == "DECEASED", 1, 0)
  
  # 1. MYC expression
  mycGeneExpression <- cbind(clinData, t(expData["MYC",]))
  mycGeneExpression$MYC_Status <- ifelse(mycGeneExpression$MYC > median(mycGeneExpression$MYC), "High", "Low")
  
  # plot survival curve
  fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ MYC_Status,
                 data = mycGeneExpression)
  p1 <- ggsurvplot(title = paste0("Total cases:", nrow(mycGeneExpression)),
    fit, data = mycGeneExpression,  
    risk.table = T, pval = T,            
    conf.int = F, xlim = c(0,80),    
    break.time.by = 10, ggtheme = theme_minimal(),
    risk.table.y.text.col = T, risk.table.y.text = F)
  fname <- file.path("results/survival_analysis/", 
                     paste0("MYCEXP_PANC_ALLSUBTYPES_TCGA_", suffix, ".eps"))
  ggsave(filename = fname,  plot = print(p1), width = 8, height = 6)
  
  # 2. Hallmark MYC signature
  mycSigHallmark <- read.delim("data/surv_analysis/MYCSig.txt", skip=1, stringsAsFactors=F)[,1]
  
  # z-score function
  myZ <- function(x) {
    x <- log2(x+1)
    out <- (x-mean(x))/sd(x)
    return(out)
  }
  tcgaExpZ <- data.frame(t(apply(expData, FUN=myZ, MARGIN=1)))
  mycHallmarkSigScore <- colSums(tcgaExpZ[mycSigHallmark,], na.rm=T)
  mycHallmarkSigScore <- data.frame(names(mycHallmarkSigScore), mycHallmarkSigScore)
  colnames(mycHallmarkSigScore) <- c("SAMPLE_ID", "MYCHallmarkScore")
  mycHallmarkSigScore[,1] <- gsub("\\.", "-", mycHallmarkSigScore[,1] )
  mycHallmarkSigScore[,"Sample_ABV"] <- gsub("-01", "",mycHallmarkSigScore[,"SAMPLE_ID"])
  tcgaSurvSigExp <- merge(mycHallmarkSigScore, mycGeneExpression, by.x="SAMPLE_ID", by.y="Sample.ID")
  
  # create groups for High MYC and low MYC using quartiles
  tcgaSurvSigExp$MYCHallmarkScore_Status <- ifelse(tcgaSurvSigExp$MYCHallmarkScore > 0, "High", "Low")
  
  # plot survival curve
  fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ MYCHallmarkScore_Status,
                 data = tcgaSurvSigExp)
  p2 <- ggsurvplot(title = paste0("Total cases:", nrow(tcgaSurvSigExp)),
    fit, data = tcgaSurvSigExp,   
    risk.table = T, pval = T,
    conf.int = F, xlim = c(0,80), 
    break.time.by = 10, ggtheme = theme_minimal(), 
    risk.table.y.text.col = T, risk.table.y.text = F)
  fname <- file.path("results/survival_analysis/", 
                     paste0("HALLMARK_PANC_ALLSUBTYPES_TCGA_", suffix, ".eps"))
  ggsave(filename = fname,  plot = print(p2), width = 8, height = 6)
  
  # 3. PDAC MYC Signature
  mycSig <- read.delim("results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt")
  mycSig <- mycSig[mycSig[,"adj.P.Val"] < 0.05,]
  mycSig <- mycSig[abs(mycSig[,"logFC"]) > log2(1.5),]
  upGenes <- as.character(mycSig[mycSig[,"DEGAnnot"] == "Up","gene_symbol"])
  downGenes <- as.character(mycSig[mycSig[,"DEGAnnot"] == "Down","gene_symbol"])
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  upGenesHuman = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = upGenes ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  downGenesHuman = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = downGenes ,mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  upGenesHuman <- intersect(upGenesHuman[,2], rownames(expData))
  downGenesHuman <- intersect(downGenesHuman[,2], rownames(expData))
  
  sigScore <- colSums(tcgaExpZ[upGenesHuman,], na.rm=T)-colSums(tcgaExpZ[downGenesHuman,], na.rm=T)
  sigScore <- data.frame(names(sigScore), sigScore)
  colnames(sigScore) <- c("SAMPLE_ID", "SignatureScore")
  sigScore[,1] <- gsub("\\.", "-", sigScore[,1] )
  sigScore[,"Sample_ABV"] <- gsub("-01", "",sigScore[,"SAMPLE_ID"])
  
  #Merge & format
  tcgaSurvSigExp <- merge(sigScore, tcgaSurvSigExp, by = "SAMPLE_ID")
  
  #Create Groups for High MYC and low MYC using quartiles
  tcgaSurvSigExp$SignatureScore_Status <- ifelse(tcgaSurvSigExp$SignatureScore>0, "High", "Low")
  
  fit <- survfit(Surv(Overall.Survival..Months., Overall.Survival.Status) ~ SignatureScore_Status,
                 data = tcgaSurvSigExp)
  
  p3 <- ggsurvplot(title = paste0("Total cases:", nrow(tcgaSurvSigExp)),
                   fit, data = tcgaSurvSigExp,   
                   risk.table = T, pval = T,
                   conf.int = F, xlim = c(0,80), 
                   break.time.by = 10, ggtheme = theme_minimal(), 
                   risk.table.y.text.col = T, risk.table.y.text = F)
  
  fname <- file.path("results/survival_analysis/", 
            paste0("SIGNATURE_PANC_ALLSUBTYPES_TCGA_", suffix, ".eps"))
  ggsave(filename = fname,  plot = print(p3), width = 8, height = 6)
}

# full clinical
clinData.full <- clinData
create.surv.plots(clinData = clinData.full, expData = expData, suffix = "full")

# subset samples
clinData.sub <- clinData[!clinData$Sample.ID %in%
                           c("TCGA-FB-A7DR-01", "TCGA-HV-A7OP-01",
                             "TCGA-HZ-8638-01", "TCGA-IB-AAUT-01",
                             "TCGA-RB-A7B8-01", "TCGA-US-A776-01",
                             "TCGA-XN-A8T5-01"),]
create.surv.plots(clinData = clinData.sub, expData = expData, suffix = "subset")

# pichai's subset
clinData.sub <- clinData.full
clinData.sub <- clinData.sub[clinData.sub$Cancer.Type.Detailed != "Pancreatic Neuroendocrine Tumor",]
clinData.sub <- clinData.sub[clinData.sub$Neoplasm.Histologic.Type.Name %in% c("Pancreas-Adenocarcinoma Ductal Type"),]
create.surv.plots(clinData = clinData.sub, expData = expData, suffix = "old")
