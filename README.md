# Maddipati_PDAC_metastasis

1.	Expression Data:

Code:
create_countdata_matrix.R
create_fpkm_matrix.R

Output Data:
Count Matrix with ensembl gene ids:
data/all_mouse_matrix_072418.RData

FPKM Matrix with ensembl gene ids:
data/all_mouse_matrix_FPKM_072418.RData

Methods:
We first aligned the fastq files of all mice to mouse mm10 genome with STAR aligner v2.5.2b. We then quantified the gene expression in terms of expected counts and FPKM (Fragments per kilobase per million mapped reads) using gencode annotation vM12 with RSEM v1.2.28. 

The RSEM files were used to create the first level of matrices of expression data where row names are ensembl gene ids and column names are mouse sample identifiers. No filters or merging was done for this. These are the starting datasets for a lot of downstream scripts.

 


2.	Collapsed Mouse/Human Matrices:

Code:
create_mouseTohuman_matrix.R

Input Data:
Count matrix with ensembl gene ids:
Data/all_mouse_matrix_072418.RData

FPKM matrix with ensembl gene ids:
Data/all_mouse_matrix_FPKM_072418.RData

Output Data:
Collapsed count matrix with gene symbols, filtered for low expression and merged clones:
Data/collapsed_counts_to_human.RData

Collapsed fpkm matrix with gene symbols, filtered for low expression and merged clones:
Data/collapsed_fpkm_to_human.RData

Methods:
The raw count and fpkm matrices from above were used to create a second level of expression matrices where gene symbols are row names and mouse identifiers are column names. First we removed low expressing genes by taking a max value of 100 for count and max value of 10 for fpkm data. Then, we merged the primary mets and primary non-mets by taking the mean value of the clones.This way we got a total number of samples to 54. We then normalized the count data using voom and batch corrected using mouse identifier. We then collapsed the matrices to unique gene symbols by taking the maximum value of the mean across all samples.

For downstream analyses like gene set enrichment or immune profiling, we filtered the matrix by all the protein coding genes and mapped the mouse gene symbols to the human gene symbols.

 

3.	Limma: All mice differential expression

Code:
main_allmice_diffexpr.R

Input Data:
Count matrix with ensembl gene identifiers:
Data/all_mouse_matrix_072418.RData

Output Data:
Differential expression limma results, count matrix subset of merged Primary mets and Primary non-mets filtered for low expression: Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/pmetvspnmet.RData

Differentially expressed genes with Adj. p-value < 0.05 for Primary mets vs Primary non-mets:
Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt

Differentially expressed genes with Adj. p-value < 0.05 and logFC > 0 for Primary mets vs Primary non-mets results annotated by gene types: Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes.txt

Differentially expressed genes with Adj. p-value < 0.05 and logFC > 0 for Primary mets vs Primary non-mets, results annotated by gene types and subset to protein coding genes and filtered for Riken genes:
Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt

Upregulated genes with Adj. p-value < 0.05 and logFC > 0:
Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_Upreg.txt

Differentially expressed genes with Adj. p-value < 0.05 and logFC > 0 for Mets vs Primary non-mets, annotated by gene types and subset to protein coding genes and filtered for Riken genes:
Analysis_01092020/Limma/Met_vs_PrimaryNonMet/MetvsPrimaryNonMet_Limma_allMice.txt

Methods:

For differential expression analysis, we started by using expected count data from RSEM for all mice. First, we removed low expressing genes by taking a cutoff of max expected count > 100 across all samples. Then, we merged the primary mets and primary non-mets by taking the mean value of the clones. This way we got a total number of samples to 54 i.e. Primary mets  (n  = 7), Primary non-mets (n = 13) and Mets (n = 34). We then normalized the count data using voom function from R package limma followed by batch correction of mouse identifiers using the R package ComBat. We used the limma package to perform differential expression between Primary mets vs Primary non-mets and Mets vs Primary non-mets.

 

4.	Box Plots

Code:
main_allmice_diffexpr.R

Output:
Boxplots for all differentially expressed genes where y-axis is log2FPKM and x-axis is the sample type. P-value represents t-test p-value:
Analysis_01092020/Boxplots/

Methods:
We created boxplots where y-axis is log2 FPKM and x-axis is the sample type. Significance was determined using t-test p-value.

 

5.	Volcano

Code:
main_allmice_diffexpr.R

Output:


Volcano of genes with Adj.p-value <  0.05 and abs. logFC > 1:
Analysis_01092020/Volcano/PrimaryMetvsPrimaryNonMetClones_VolcanoPlot_allMice_pc_logfc1_pval005.eps

Methods:
For Volcano plots, we created two different plots one with Adj. p-value < 0.05 and abs logFC > 1 

 
 

6.	Enrichment: GSEA (Web version)

Code:
GSEA_allmice.R
 
Input Data:
Analysis_01092020/Enrichment/GSEA_web_version/input_files

Output Data:
Analysis_01092020/Enrichment/GSEA_web_version

Methods:
For gene set enrichment analysis, we used the expression matrix of expressed genes (using the cut-off of max expected count across all samples > 100) from Primary mets and Primary non-mets (total samples = 20). We mapped a total of 15629 mouse genes to 13514 human genes. We then log-normalized the count data and created input files for the GSEA web version. Then, we performed enrichment using different gene sets from MSigDB like GO, Hallmark, TF, KEGG, Immunologic signatures and Oncogenic signatures.

 

8.	Enrichment: TF Enrichment using Metacore

Input Data:
Differentially expressed genes with Adj. p-value < 0.05 for Primary mets vs Primary non-mets:
Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt

Output Data:
Analysis_01092020/Enrichment/Metacore

Methods:
We performed a Transcription factor enrichment using Metacore on all differentially expressed genes (n = 3585) from Primary mets vs Primary non-mets analysis. We found Myc as the top Transcription factor significantly enriched using this analysis.

 

9.	Myc motif enrichment

Code:
Myc motif enrichment:
IdentifyMYCRegulatedGenes_allmice.R
Annotate Limma differential output with Myc binding and expression data: annotate_limma_output.R

Input Data:
Analysis_01092020/Limma/PrimaryMet_vs_PrimaryNonMet/PrimaryMetvsPrimaryNonMet_Limma_allMice.txt

Output Data:
Analysis_01092020/MYC_binding/

Methods:
We used the MotifDb to search for a Myc binding motif 500 bp upstream in the promoter region of each expressed gene (n = 13,627). We used a cutoff of 85% match as a hit to determine the presence of Myc motif. We then compared the probability of Myc binding enrichment upstream of upregulated protein-coding genes vs all expressed genes and got a fisher-test p-value of 0.059.

10.	ReactomeFI (Reactome Pathway Enrichment)

Code:
Script to create input files for ReactomeFI plugin:
ReactomeFI_input.R

Output Data:
Analysis_01092020/ReactomeFI

Methods:
We used the above gene lists to find enrichment of Reactome Pathways with the ReactomeFI plugin of Cytoscape.

11.	IPA analysis: Bar Plot of Upstream regulators

Code:
Barplot function for upstream regulators:
barplot_upstream_regulators.R

Input Data:
Upstream regulators analysis for differentially expressed genes with logFC > 1 and Adj. p-value <0.01: Analysis_01092020/IPA/Networks/all_diffexpr_Adjpval001_logFC1/all_diffexpr_Adjpval001_logFC1_upstream_regulators.txt

Methods:
We used all genes with Adj P-value < 0.01 and logFC > 1 to identify enrichment of upstream regulators using IPA. For this script, we use the upstream regulator output from IPA to create bar plots where y-axis are the upstream regulators and x-axis is -log10 P-value.

 

12.	 IPA analysis (Network and Upstream Regulators)

Input Data:
All differentially expressed genes

Output Data:
Analysis_01092020/IPA/Networks

Methods:
Upstream regulator analysis was done on Met High, Met Low and Diff.expr genes. Network analysis was done

13.	Survival Analysis
Code:
surv_analysis/survival-pichai/analysis_V2.R

Input Data:
TCGA PAAD expression data:
Data/surv_analysis/data_RNA_Seq_v2_expression_median.txt

TCGA PAAD patient data:
Data/surv_analysis/data_bcr_clinical_data_patient.txt

TCGA PAAD sample data:
Data/surv_analysis/data_bcr_clinical_data_sample.txt

Differentially expressed genes between Primary met vs Primary non-mets with Adj. P-value < 0.05:
Data/surv_analysis//PrimaryMetvsPrimaryNonMet_Limma_allMice_allBiotypes_diffexpr.txt


Methods:
For survival analysis, we downloaded the TCGA PAAD expression and patient and sample level clinical data from cBioportal. We filtered the samples by subtype with Pancreatic Neuroendocrine Tumor and histological diagnosis with Pancreas-Adenocarcinoma Ductal Type. After filtering we obtained 147 samples.
We also used Met high signature from the differential expression analysis between Primary mets vs Primary non-mets to find correlation to survival. We took the 3585 differentially expressed genes and filtered them using abs. logFC > 1.5 to get a list of 736 highly upregulated and 1036 highly downregulated genes. We then subtracted the sum of all downregulated gene expression from the upregulated gene expression for each sample and obtained the Primary mets signature score. Similar to the hallmark signature, we divided this PMet-high signature into high and low strata by taking a cutoff of >0. We found that Pmet high signature is significantly correlated with survival at a P-value of 0.0053.

14. Heatmaps
Code:
heatmap_allmice.R

Input Data:
Expected count matrix: 
Data/all_mouse_matrix_072418.RData

Output Data:
Analysis_01092020/heatmaps_allmice
Prior to plotting the heatmaps, the expression matrix was filtered for low expression (max expected count across all samples > 100), only protein coding genes and log2 normalized.

15.	Molecular Subtyping

Code:
Classification of molecular subtypes: 
molecular_subtypes.R

Helper scripts for classification:
molecular_subtypes_helper.R

Heatmaps with predicted classification: 
molecular_subtypes_heatmaps.R

Input Data:
Voom normalized count matrix after filtering for low expressing genes, merging clones, batch correcting by mouse id and mapping to human gene symbols:
Data/collapsed_counts_to_human.RData

Output Data:
Molecular signature classification for three signatures:
Analysis_01092020/molecular_subtyping/molecular_subtypes.RData

Methods:
We used three different molecular signatures to classify the samples i.e. Bailey, Moffitt and Collison. For each sample, we calculated a score by subtracting the sum of normalized expression for genes corresponding to specific classes within a particular molecular signature. We then took the maximum score observed across each class and assigned to the samples.

After predicting the molecular subtypes, we created several heatmaps and PCA plots with and without met samples. 

For heatmaps, we performed clustering using Differentially expressed genes with Adj. P-value < 0.05 and abs. logFC > 1 (n =  504) with samples annotated with predicted molecular subtypes


