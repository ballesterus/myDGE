library(ggplot2)
library(DESeq2)
library(reshape)



readCSVcountTable<-function(fname){
    ctable<-read.csv(fname, header=TRUE, sep=",")
    rownames(ctable)<-ctable[,1]
    ctable<- ctable[, -1]
    return(ctable)
}

normalizeRLS <- function(ctable){
    dummycons<-data.frame(condition = gsub("_.+", "", names(ctable)), row.names=names(ctable))
    DESeq.dat<-DESeqDataSetFromMatrix( countData = ctable, colData = dummycons, design = ~ condition)
    DESeq.dat<- DESeq.dat[rowSums(counts(DESeq.dat)) > 0,]
    DESeq.rlog <- rlog(DESeq.dat, blind =TRUE)
    return(DESeq.rlog)
    }

rundDESeq<-function(ctable, expdesign){
    DESeq.dat<-DESeqDataSetFromMatrix( countData = ctable, colData = conditions, design = ~ condition)
    colData(DESeq.dat)$condition <- relevel(colData(DESeq.dat)$condition, baseCondition)
    De<-DESeq(DESeq.dat)

    
    R<-results(De, independentFiltering=TRUE, alpha=0.05, contrast)
    return(R)
}


conditionsfromTable<-function(ctable,exp){
    r<-data.frame(row.names=names(ctable),condition=gsub(exp,"",names(ctable)))
    return(r)
    
}

pairwiseConstrastlist <- function(conditionTable){
    e<-unique(conditionTable$condition)
    r<-data.frame(matrix(ncol=3,nrow=0))
    for (i in e){
        for (j in e){
            if (i != j ){
                v <- c("condition",i,j)
                r <- rbind(r,v)
                }
            }
        }
    names(r)<-c("fact", "num", "denom")
    return(r)
    }

subsetTreatment <- function(ctable, ...){
    l <- list(...)
    r<-data.frame(row.names=rownames(ctable))
    for (p in l){
        mc<-C[,grep(p, names(C))]
        r<- cbind(r,mc)
    }
    return(r)
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
