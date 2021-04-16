library(ggplot2)
library(DESeq2)
library(reshape)
library(edgeR)
library(NMF)

myplotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
                 factor(apply(intgroup.df, 1, paste, collapse = ":"))
             }
        else {
                 colData(object)[[intgroup]]
             }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
        ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance"))
    + coord_fixed()
}


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


runEdgeR<-function(ctable,contable,refgroup){
    g <- relevel(factor(contable$condition), ref=refgroup)
    s <- DGEList(counts=ctable, group= g)
    f <- filterByExpr(s)
    s <- s[f,]
    s$samples$lib.size<-colSums(s$counts)
    design <- model.matrix(~group, data=y$samples)
    y <-estimateDisp(y,design)
    fit <-glmQLFit(y,design)
    for(i in 2:length(levels(g))){
        oname <- paste("edgeR_", levels(g)[i], "vs", refgroup, sep="")
        dir.create(oname)
        glf<-glmQLFTest(fit, coef=i)
        o <- topTags(glf, n= Inf, adjust.method = "BH")
        write.csv(o,file=paste(oname, "/", oname, "result.csv", sep=""))
    }
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

normalizeTMM <- function(ctab){
    dgl<- DGEList(counts=ctab, genes=rownames(ctab))
    dgl<-calcNormFactors(dgl, method= "TMM")
    dgl <- cpm(dgl)
    return(dgl)
}


getNormCounts<-function(gidlist,normCountsdf){
    D<- normCountsdf[gidlist,]
    return(D)
}


plotHM<-function(adf){
    aheatmap(adf, Rowv=TRUE, Colv=TRUE, distfun="euclidean", hclustfun="average", scale ="row", color="RdYlGn")
}


plotVolcano<-function(data, pcol,fccol){
    s<-data$pcol
    s[s < 0.05]<- "de"
    s[s >= 0.05]<- "non-de"
    D<-cbind(data, s)
    print(nrow(D))
    ggplot(D) +
        geom_point(aes(x = fccol, y=pcol, colour=s))+
        xlab("log2 fold change") +
        ylab("-log10 adjusted p-val")
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


##plotMA(R, alpha=0.05, ylim=c(-10,10), colSig="red3"s
