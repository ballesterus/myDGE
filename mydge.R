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

runDESeq<-function(ctable, expdesign){
    DESeq.dat <- DESeqDataSetFromMatrix( countData = ctable, colData = expdesign, design = ~ condition)
    De <-DESeq(DESeq.dat)
    pwcomps<- getpairwiseContrasts(expdesign)
    for(i in 1:nrow(pwcomps)){
        oname<- paste(pwcomps[i,2], "vs", pwcomps[i,3], sep= "")
        dir.create(oname)
        r<- results(De, independentFiltering =TRUE, alpha=0.05, contrast = c(pwcomps[i,1], pwcomps[i,2], pwcomps[i,3]))
        print(summary(r))
        write.csv(r,file=paste(oname,"/DESeq_", oname, "_results.csv", sep=""))
        pdf(paste(oname, "/MAplot_", oname, ".pdf", sep =""))
        plotMA (r, alpha = 0.05 , main = oname , ylim = c ( -10 ,10) )
        dev.off()
 
    }
}

conditionsfromTable<-function(ctable,exp){
    r<-data.frame(row.names=names(ctable),condition=gsub(exp,"",names(ctable)))
    return(r)    
}

getpairwiseContrasts <- function(conditionTable){
    e<-unique(conditionTable$condition)
    r<-data.frame(matrix(ncol=3,nrow=0))
    f<-names(conditionTable)
    for (i in e){
        for (j in e){
            if (i != j ){
                v <- c(f,i,j)
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
