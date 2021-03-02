library(ggplot2)
library(DESeq2)
library(reshape)
library(ggfortify)


readCSVcountTable<-function(fname){
    ctable<-read.csv(fname, header=TRUE, sep=",")
    rownames(ctable)<-ctable[,1]
    ctable<- ctable[, -1]
    return(ctable)
}

DESeqnormCs <- function(ctable){
    dummycons<-data.frame(condition = gsub("_.+", "", names(ctable)), row.names=names(ctable))
    DESeq.dat<-DESeqDataSetFromMatrix( countData = ctable, colData = dummycons, design = ~ condition)
    DESeq.dat<- DESeq.dat[rowSums(counts(DESeq.dat)) > 0,]
    DESeq.rlog <- rlog(DESeq.dat, blind =TRUE)
    return(DESeq.rlog)
    }

RundDGE<-function(ctable, conditions, baseCondition){
    DESeq.dat<-DESeqDataSetFromMatrix( countData = ctable, colData = conditions, design = ~ condition)
    colData(DESeq.dat)$condition <- relevel(colData(DESeq.dat)$condition, baseCondition)
    De<-DESeq(DESeq.dat)
    R<-results(De, independentFiltering=TRUE, alpha=0.05)
    return(R)
}





#### EXAMPLE BOXPLOT
### melt the dataframe

## Dm<-melt(D)

###Add a column

##a<-gsub("_.+", "", s$X2)
##Dm<-cbind(Dm,a)

##ggplot(nb,aes(x=nb$X2, y=value, fill=app )) + geom_boxplot() + geom_jitter(alpha = 0.001)


###EXAMPLE principal component analysis
### plotPCA function does not uses ALL variables (genes) for PCA, by default using only the top 500 genes with highest count (row) variance. 
##plotPCA(N,ntop= 100)

###MAplot
##plotMA(R, alpha=0.05, ylim=c(-10,10), colSig="red3")
