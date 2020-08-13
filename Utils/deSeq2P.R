###########################################
#Code to run DeSeq2 on a htseq count matrix
#
#Function requires a matrix of counts and 2 or more groups
#returns back a list
#
#First object is gene, Log-FC, p-value
#Second object is the actual deseq2 output
#
###########################################

#####################
#EXAMPLES
#
#2 Groups:
#load("../data/mouse.brain.htseq.rda");
#groupsName <- list(colnames(mouse.brain.htseq)[1:5], colnames(mouse.brain.htseq)[6:10]);
#groupsNum <- list(c(1), c(10));
#deSeq2Output<- deSeq2P(mouse.brain.htseq, groupsNum);
#
#####################


#Function toDeSeq2, by default first group is control, only 2 group comparisons allowed at this point
deSeq2P <- function(mtrx=NULL, grps=NULL)
{
 #Call libraries
require("DESeq2");

#reassign vars
dataMat <- mtrx;
groups <- grps;

#Handle if its a numeric vector
if(class(groups[[1]])=="integer")
{
tmpGroups <- colnames(dataMat)[groups[[1]]];
for(i in 2:length(groups))
{
tmpGroups <- list(tmpGroups, colnames(dataMat)[groups[[i]]]);
}
groups <- tmpGroups;
}
    
#Pull out columns from dataMat
dataMat <- dataMat[,unlist(groups)];
    
#Generate targets file
targets <- data.frame(groups[[1]],paste("Group_", rep(1, length(groups[[1]])), sep=""))
colnames(targets) <- c("Sample", "Condition");
for(i in 2:length(groups))
{
tmpTargets <- data.frame(groups[[i]],paste("Group_", rep(i, length(groups[[i]])), sep=""))
colnames(tmpTargets) <- c("Sample", "Condition");
targets <- rbind(targets, tmpTargets);
}
rownames(targets) <- targets[,1];
targets <- targets[-1];

    dds <- DESeqDataSetFromMatrix(countData = dataMat, colData = targets, design = ~ Condition);
    dds <- DESeq(dds);
    res <- results(dds);
    res <- data.frame(res);
    res <- res[order(res[,"pvalue"]),];
    output <- list(res[,c("log2FoldChange", "pvalue")], res, dataMat);
    return(output);
}

