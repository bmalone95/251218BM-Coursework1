#inQuestion 1 - Read in the differential expression table and produce a data.frame of all results. How many genes have a padj < 0.05.
DEgenes <- read.delim("DE_Genes/GM12878_Minus_HeLa_DEG.csv",sep= "\t",header=TRUE)
DEgenes <- read.delim("DE_Genes/GM12878_Minus_HeLa_DEG.csv",sep= ",",header=TRUE)
DEgenespadj <- DEgenes[DEgenes$padj <0.05,]
count(DEgenespadj)
4719
#Question 2 - Now with these genes with a padj < 0.05, create a scatter plot (as seen below) of -log10 pvalues on Y axis and log2FoldChange on X axis using ggplot2.
log10pval <- -log10(DEgenespadj$pvalue)
DEgenespadj2 <- cbind(DEgenespadj,log10pval)
DEgenes_scatter <- ggplot(data = DEgenespadj2,mapping = aes(x=log2FoldChange,y=log10pval)) + geom_point()
DEgenes_scatter
DEgenes_scatter + ggtitle("Volcano plot of GM12878 Minus HeLa") + ylab("-log10 of P value")+ theme_minimal()
ggsave(DEgenes_scatter,filename="Exercise1pt2_scatter.png",width=15,height=15,units="cm")
#Question 3 - Read in the absolute expression table, add 1 to every value in table and make a boxplot of log10 expression values for all samples.
ABSgenes <- read.delim("DE_Genes/Expression.csv",sep=",",header = TRUE)
ABSgenesadd1 <- ABSgenes[,2:5]+1.0000001
Expressionlog10 <- log10(ABSgenesadd1[,1:4])
LogABSgenes <- cbind(ABSgenesadd1,Expressionlog10)
boxplot(LogABSgenes[,5:8])
boxplot(LogABSgenes[,5:8],boxwex=0.5,names = c("GM12878_1","GM12878_2","HeLa_1","HeLa_2"),las=2)
##Question 4 - Now create a similar boxplot with just genes that have a padj < 0.05 and a log2FoldChange > 1
ABSDEmerge <- merge(ABSgenes,DEgenes,by="ID")
padjbox <- ABSDEmerge[ABSDEmerge$padj <0.05, ]
padjlog2box <- padjbox[padjbox$log2FoldChange>1,]
padjlog2boxadd1 <- padjlog2box[,2:5]+1.000001
log10exp <- log10(padjlog2boxadd1[,1:4])
boxplot(log10exp[,1:4],boxwex=0.5,names = c("GM12878_1","GM12878_2","HeLa_1","HeLa_2"),las=2)
##Question 5 - Using the absolute expression table, identify the genes whose expression is in the top 60%. Filter the results from the differential expression table to these results and plot the log2 basemean on X and log2FoldChange on Y. Highlight genes who have padj < 0.05
averageegeneexpression <- rowMeans(ABSgenes[1:23434,2:5])
ABSgenesave <- cbind(ABSgenes,averageegeneexpression)
ABSgenes60 <- subset(ABSgenesave,averageegeneexpression >= quantile(averageegeneexpression,0.4))
DEABS60merge <- merge(ABSgenes60,DEgenes)
DEscatter <- ggplot(DEABS60merge,mapping=aes(x=log2(baseMean),y=log2FoldChange,colour=padj)) + geom_point() +theme_minimal() + scale_colour_gradientn(colours=c("cyan","red","red"),values = c(0,0.05,1),labels=c("NoSig","Sig","Sig"),breaks=c(0,0.05,1))
DEscatter
##Part 2
##Question 1 - Read in H3K27Ac_Limb_1.txt file and report the number of genomic locations listed in file.
HOMER <- read.table("HOMER_peaks/H3K27Ac_Limb_1.txt")
summary(HOMER)
HOMERGranges <- GRanges(HOMER$V2,IRanges(HOMER$V3,HOMER$V4),RegSize=HOMER$V7)
### 20447 peaks listed
###Question 2 - Make a histogram of the log10 of regions sizes as shown below using base graphics.
hist(log10(HOMER[,7]),main="Histogram of log10 of region sizes",xlab = "region sizes (log10)")
### Question 3 - Make a density plot of the log10 of regions sizes as shown below using ggplot graphics.
HOMERdensity <- ggplot(HOMER,mapping=aes(x=log10(HOMER[,7]))) +geom_density(colour="dark green",fill="dark green")+theme_minimal()+labs(x="region sizes (log 10)",title = "Histogram log10 of region sizes")
###Question 4 - Make a density plot for each chromosome of the log10 of regions sizes as shown below using ggplot graphics
HOMERdensitywrap <- ggplot(HOMER,mapping=aes(x=log10(HOMER[,7]))) +geom_density(colour="dark green",fill="dark green")+theme_minimal()+labs(x="region sizes (log 10)",title = "Histogram log10 of region sizes")+facet_wrap(HOMER$V2)
###Question 5 - Make a boxplot plot of the log10 of findPeaks.Score for each chromosome as shown below using ggplot graphics.
HOMERboxwrap <- ggplot(HOMER,mapping=aes(y=log10(HOMER[,9]),x=HOMER$V2,fill=HOMER$V2)) +geom_boxplot()+theme_minimal()+labs(y="region sizes",x="region size(log 10)",title = "Histogram log10 of region sizes")+ scale_fill_hue(l=60,c=80)+coord_flip()
##Question 6 - Export the Homer genomic regions as a BED3 file.
library(nucleR)
export.bed(HOMERGranges,name = "HomerGrange.bed")

