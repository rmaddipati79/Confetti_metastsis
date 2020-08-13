###############################
#Code to compare 2 lists, coming from
#processing through voom/limma or DE, input
#should be a df with rownames as ID's or genes 
#first column is FC and second p-value
#You can compare either via a scatter, venn diagram
#or you can compare with an ROC type curve
###############################


#Comparing lists with a scatter plot
scatterCompare <- function(gl1=NULL, gl2=NULL, comp="pval")
{
#Call libraries
require(ggplot2);

if(comp=="pval")
{
log=T;
col=2;
}
if(comp=="fc")
{
log=F;
col=1;
}

#Make sure they are both numeric
gl1[,col] <- as.numeric(gl1[,col]);
gl2[,col] <- as.numeric(gl2[,col]);
gl1 <- gl1[rownames(gl2),];
gl1x <- gl1[,col];
gl2y <- gl2[,col];


if(log==T)
{
gl1x <- (-1)*log10(gl1x);
gl2y <- (-1)*log10(gl2y);
}

xLoc <- sort(runif(1000, min(na.omit(gl1x)), max(na.omit(gl1x))))[200];
yLoc <- sort(runif(1000, min(na.omit(gl2y)), max(na.omit(gl2y))))[800];

#Function to display R-square
lm_eqn = function(x, y){
    m = lm(y ~ x);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(coef(m)[1], digits = 2), 
              b = format(coef(m)[2], digits = 2), 
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
}

p <- qplot(gl1x, gl2y)+geom_point()+geom_smooth(method="lm")+theme_bw()+geom_text(aes(x=xLoc, y=yLoc, label=lm_eqn(gl1x, gl2y)), parse=T);


}



#Comparing lists with a Venn Diagram, can only compare p-values
vennCompare <- function(gl1=NULL, gl2=NULL, num=500, pval=.01, useNum=T, plot=T)
{
#Call libraries
require(VennDiagram);
require(gridExtra)

#Compare p-value column
col=2;

#Make sure they are both numeric
gl1 <- na.omit(gl1);
gl2 <- na.omit(gl2);
gl1[,col] <- as.numeric(gl1[,col]);
gl2[,col] <- as.numeric(gl2[,col]);

if(useNum==T)
{
gl1 <- gl1[order(gl1[,col]),];
gl2 <- gl2[order(gl2[,col]),];
myX <- rownames(gl1)[1:num];
myY <- rownames(gl2)[1:num];
}
if(useNum==F)
{
myX <- gl1[gl1[,2]<pval,];
myY <- gl2[gl2[,2]<pval,];
myX <- rownames(myX);
myY <- rownames(myY);
}

#Functions to run functional enrichment on any gene set 
runHypGeom <- function(set, genes,n=20000)
{
#number of white balls
x <- length(intersect(genes, set));

#white balls
m <- length(genes);

#black balls
n2 <- n-m; 

#balls drawn from the urn 
k <- length(set);


out <- phyper(x-1, m, n2, k, lower.tail=F);
setSize <- k;
overLap <- x;
numGenes <- m;

myRet <- c(setSize, numGenes, overLap, out); 
return(myRet);

}

myP <- round(runHypGeom(myX, myY, n=length(rownames(gl1)))[4],3);

if(plot==T)
{
p <- draw.pairwise.venn(length(myX), length(myY), length(intersect(myX, myY)), c("gl1", "gl2"),  fill = c("blue", "red"), alpha=.5, label.col="black", cex = 2,cat.fontface = 1, cat.cex=2, cat.dist=0.05);
grid.arrange(gTree(children=p), main=paste("P-value of Intersection :", myP))
}
return(length(intersect(myX, myY))/length(union(myX, myY)));
}


gl1 <- gfpmchComp_DS2[[1]];
gl2 <- gfpmchComp_VL[[1]];

vennCompare(gl1, gl2, plot=T);





#Comparing lists with an ROC curve
rocCompare <- function(gl1=NULL, gl2=NULL, useNum=T, dens=100)
{
#Call libraries
require(ggplot2);

#Compare p-value column
col=2;

#Make sure they are both numeric 
gl1[,col] <- as.numeric(gl1[,col]);
gl2[,col] <- as.numeric(gl2[,col]);

for(i in 1
vennCompare(gl1, gl2, plot=T);





}


































