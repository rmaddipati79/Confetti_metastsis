##############################################
#Purpose: R Utilities and often used function
#Author: Pichai Raman
#Date: 10/10/2018
##############################################



#Create DF with unique rows, when you have genes as a column and unwanted columns
getUniqueDF <- function(myDF=NULL, rowCol=NULL, remCols=NULL, chooseFun="max")
{
	#Print number of rows initially
	print(paste("Data frame has", nrow(myDF), "rows"));

	#Remove Rows
	myDF <- myDF[, setdiff(colnames(myDF), remCols)]

	#Remove duplicates
	myDF[,"AGG_VALUE_TMP"] <- apply(myDF[-1], MARGIN=1, FUN=chooseFun)
	myDF <- myDF[order(-myDF[,"AGG_VALUE_TMP"]),];
	myDF <- myDF[!duplicated(myDF[,rowCol]),]
	myDF <- myDF[, setdiff(colnames(myDF), "AGG_VALUE_TMP")]

	#Set Rownames
	myDF <- myDF[!is.na(myDF[,1]),]
	rownames(myDF) <- myDF[,1];
	myDF <- myDF[-1];

	#Print number of rows initially
	print(paste("Data frame now has", nrow(myDF), "rows"));

	return(myDF);
}
