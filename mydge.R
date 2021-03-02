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
    
    normc <- assay(DESeq.rlog)
    return(normc)
    }






#### EXAMPLE BOXPLOT
### melt the dataframe

## Dm<-melt(D)

###Add a column

##a<-gsub("_.+", "", s$X2)
##Dm<-cbind(Dm,a)

##ggplot(nb,aes(x=nb$X2, y=value, fill=app )) + geom_boxplot() + geom_jitter(alpha = 0.001)


#EXAMPLE principal component analysis
#pcaN<-prcomp(t(N))
