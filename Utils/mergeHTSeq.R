####################################
#Code simply takes output in a htseq folder and
#binds it together in a matrix
####################################

########################
#Example
########################


#Call to merge files, just needs a directory location
mergeHTSeq <- function(myDir=NULL)
{
    myFiles <- list.files(myDir);
    myDF<- read.delim(paste(myDir, myFiles[1], sep="/"), header=F);

    for(i in 2:length(myFiles))
    {
        tmpDF<- read.delim(paste(myDir, myFiles[i], sep="/"), header=F);
        myDF <- cbind(myDF, tmpDF[2]);
    }

    myHeader <- c("ID", gsub(".htseq.tsv", "", myFiles));
    colnames(myDF) <- myHeader;
    rownames(myDF) <- myDF[,1];
    myDF <- myDF[-1];
    return(myDF);
}
