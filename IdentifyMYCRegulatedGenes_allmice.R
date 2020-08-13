#########################################
# code to identify genes regulated by MYC
# Author: Pichai Raman
# Date: 12/6/2016
########################################

# Alright first we need to load all the libraries
library(seqinr)
library(biomaRt)
library(Biostrings)
library(MotifDb)
library(MotIV)
library(seqLogo)
library(motifStack)
library(Biostrings)
library(GenomicFeatures)

# Now let's create a function that takes in a gene, pulls the sequence and then searches for motifs in MotifDB
# Accessory Function to pull relevant information  from matchPWM function
matchPWMHelp <- function(motif, subject, gene, min.score) {
  out <- matchPWM(motif, subject, min.score)
  rangeOut <- as.data.frame(ranges(out))
  if(dim(rangeOut)[1]==0)
  {
    rangeOut <- data.frame(t(c(0,0,0)))
  }
  res <- data.frame(gene, rangeOut)
  colnames(res) <- c("gene", "start", "end", "width")
  return(res)
}


#Function, just enter gene and how far upstream in promoter you want to search
pullMatchMotifDB <- function(SeqObj, motifs, min.score="90%") {
  #Pull Sequence
  tmpSeq <- SeqObj[[1]]
  gene <- SeqObj[[2]]
  #Get matches
  print(gene)
  if(length(tmpSeq)>0)
  {
    tmpOut <- lapply(motifs, FUN=matchPWMHelp, subject=tmpSeq, gene=gene, min.score=min.score)
    out <- do.call("rbind", tmpOut)
    out <- data.frame(out)
  }
  if(length(tmpSeq)==0)
  {
    out <- c(1:length(motifs))
    out <- data.frame(out)
  }
  return(out)
}

#Initialize Variable
myMart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#Have to get genes first
data <- read.delim('results/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt')
myGenesEG <- unique(as.character(data[,"gene_id"]))

#Okay refactor code a bit, pull all sequences in one fell swoop
tmpSeqObj <- biomaRt::getSequence(id = myGenesEG, type = "ensembl_gene_id", seqType = "gene_flank", upstream = 500, mart = myMart)

#Now have to get motifs
MYCMotif <- query (query (MotifDb, 'Mmusculus'), 'MYC')

#Now query and find motif in strings
MotifSearch <- apply(tmpSeqObj, FUN = pullMatchMotifDB, MARGIN = 1, motifs = MYCMotif, min.score = "85%")
MYC_SEARCH <- do.call("rbind", MotifSearch)
MYC_SEARCH[,"TF"] <- rownames(MYC_SEARCH)

removeTrail <- function(x) {
  ind <- regexpr("\\.", x)[[1]]-1
  nx <- substring(x, 1, ind)
  return(nx)
}
MYC_SEARCH[,"TF"] <- as.character(lapply(as.character(MYC_SEARCH[,"TF"]), FUN = removeTrail))
justMYC <- grep("-Myc-", MYC_SEARCH[,"TF"])
MYC_SEARCH <- MYC_SEARCH[justMYC,]
MYC_SEARCH[,"HIT"] <- ifelse(MYC_SEARCH[,"width"] > 0, 1, 0)
MYC_SEARCH <- MYC_SEARCH[MYC_SEARCH["HIT"] > 0,]

#Final list of genes that have a binding motif
tfBindingSiteResults <- unique(as.character(MYC_SEARCH[,"gene"]))

data.sub <- data[as.character(data$gene_id) %in% tfBindingSiteResults,] 

data$HIT <- 'NO'
data$HIT[which(data$gene_id %in% tfBindingSiteResults)] <- "YES"
write.table(data, file = 'results/MYC_binding/PrimaryMetvsPrimaryNonMet_allMice_MYCBinding.txt', quote = F, sep = "\t", row.names = F)

# calculate proportions
x <- data
x$status <- ifelse(x$DEGene == T & x$logFC > 0, 'Upreg', 
                   ifelse(x$DEGene == T & x$logFC < 0, 'Downreg', 'NotSig'))
x <- plyr::count(x, c('status','HIT'))

# upregulated
ftest.upreg <- matrix(c(418, 2367, 1308, 8161), ncol = 2, nrow = 2, dimnames = list(c('Upreg','NDE'),c('Present','Not-Present')))
fisher.test(ftest.upreg, alternative = 'greater')

# downregulated
ftest.dnreg <- matrix(c(311, 2367, 1062, 8161), ncol = 2, nrow = 2, dimnames = list(c('Downreg','NDE'),c('Present','Not-Present')))
fisher.test(ftest.dnreg, alternative = 'greater')
