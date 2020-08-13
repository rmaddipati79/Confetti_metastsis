###############################
# Author: Komal S Rathi
# Function: Molecular Subtype signatures
# and PCA plots
# Date: 09/18/2018
###############################

# load libraries
library(scatterplot3d)
library(ggplot2)
library(grid)
library(gridExtra)
library(pheatmap)

# set working directory
setwd('~/Projects/PancModelUpdate/')

# source code
source('code/molecular_subtyping/molecular_subtypes_helper.R')

# load expression matrix 
# removed low expressing genes, merged clones, batch corrected by mouse id and mapped to human gene symbols
load('data/collapsed_counts_to_human.RData')

# Bailey signature
baileySig <- read.delim("data/molecular_subtyping/Bailey.txt", stringsAsFactors = F)
bailey.sig <- mol.subtype(signature = baileySig, exprs = countDataPrim)

# Moffitt
moffittSig <- read.delim("data/molecular_subtyping/Moffitt.txt", stringsAsFactors = F)
moffittSig <- moffittSig[which(moffittSig$Subtype %in% c("Basal-Like", "Classical")),] # only basal and classical
moffitt.sig <- mol.subtype(signature = moffittSig, exprs = countDataPrim)

# Collison 
collisonSig <- read.delim('data/molecular_subtyping/collison.txt', stringsAsFactors = F)
collison.sig <- mol.subtype(signature = collisonSig, exprs = countDataPrim)

# save here and use in heatmap generation
save(bailey.sig, collison.sig, moffitt.sig, file = 'results/molecular_subtyping/molecular_subtypes.RData')

# PCA - all mice
create.pca.allmice(signature = bailey.sig, 
                   exprs = countDataPrim,
                   pca.title = "Bailey Signature", 
                   file.name = "Bailey_Signature_AllMice", 
                   dir.name = 'results/molecular_subtyping/bailey')
create.pca.allmice(signature = collison.sig, 
                   exprs = countDataPrim,
                   pca.title = "Collison Signature", 
                   file.name = "Collison_Signature_AllMice", 
                   dir.name = 'results/molecular_subtyping/collison')
create.pca.allmice(signature = moffitt.sig, 
                   exprs = countDataPrim,
                   pca.title = "Moffitt Signature", 
                   file.name = "Moffitt_Signature_AllMice", 
                   dir.name = 'results/molecular_subtyping/moffitt')

# Stacked Barplot
create.barplot(signature = bailey.sig,
               bar.title = "Bailey Signature",
               file.name = "Bailey_Signature_summary", 
               dir.name = 'results/molecular_subtyping/bailey')
create.barplot(signature = collison.sig,
               bar.title = "Collison Signature",
               file.name = "Collison_Signature_summary", 
               dir.name = 'results/molecular_subtyping/collison')
create.barplot(signature = moffitt.sig,
               bar.title = "Moffitt Signature",
               file.name = "Moffitt_Signature_summary", 
               dir.name = 'results/molecular_subtyping/moffitt')

# PCA individual mice 
create.pca.indmice(signature = moffitt.sig, 
                   exprs = countDataPrim,
                   pca.title = "Moffitt Signature", 
                   file.name = "Moffitt_Signature_Mouse",
                   dir.name = 'results/molecular_subtyping/moffitt')
create.pca.indmice(signature = bailey.sig, 
                   exprs = countDataPrim,
                   pca.title = "Bailey Signature", 
                   file.name = "Bailey_Signature_Mouse", 
                   dir.name = 'results/molecular_subtyping/bailey')
create.pca.indmice(signature = collison.sig, 
                   exprs = countDataPrim,
                   pca.title = "Collison Signature", 
                   file.name = "Collison_Signature_Mouse", 
                   dir.name = 'results/molecular_subtyping/collison')
