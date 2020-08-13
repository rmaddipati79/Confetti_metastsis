##########################################
#Code to pull out all the reads associated with a gene
#or transcript as well as plot it.
##########################################


#Call libraries and data
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

#Read in annotation Data
geneAnnotEntrez <- read.delim("/Users/pichairaman/Documents/Data/AnnotationFiles/Entrez_gene/Homo_sapiens.gene_info", header=F);
geneAnnotEntrez <- geneAnnotEntrez[,c("V3", "V2", "V12", "V8")]
colnames(geneAnnotEntrez) <- c("GENE_SYMBOL", "GENE_ID", "GENE_DESC", "GENE_LOC");
exonsGene <- exonsBy(txdb, by = "gene")

#Enter a gene and back a GRanges object associate with it
getExonsForGene <- function(x)
{
x.gi <- as.character(na.omit(geneAnnotEntrez[geneAnnotEntrez[,"GENE_SYMBOL"]==x,"GENE_ID"]));
exons.x <- exonsGene[[x.gi]];
return(exons.x);
}
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=getExonsForGene("EGFR"), what=what)
bamFile <- "/Users/pichairaman/Documents/CodeTools/PullDataFromBam/data/TARGET-30-PAUDDK-02A.sorted.bam";

bam <- scanBam(bamFile, param=param)
