# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting
getPCA = function (x, intgroup = "condition", ntop = 500, pcs=3, cdata=NULL) 
    {     
        library(genefilter)
        library(GenomicRanges)
        
        if(class(x) == 'SeqExpressionSet'){
            cdata = pData(x)
            x = counts(x)
            
        }
            
        if(class(x) %in% c('SummarizedExperiment','DESeqDataSet')){
            cdata = colData(x)
            x = assays(x)[[1]]
            
        }
        if(class(x) == 'matrix'){
            cdata = cdata
            x = x
            
        }
            
        
        rv <- rowVars(x)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                           length(rv)))]
        library(pcaMethods)
        pca <- prcomp(t(x[select, ]))
        percentVar <- pca$sdev^2/sum(pca$sdev^2)
        if (!all(intgroup %in% names(cdata))){
            stop("the argument 'intgroup' should specify columns of cdata")
        }         
        intgroup.df <- as.data.frame(cdata[, intgroup, drop = FALSE])
        group <- factor(apply(intgroup.df, 1, paste, collapse = " : "))
        pc = data.frame(pca$x[,1:min(pcs, ncol(pca$x))])
        d <- data.frame(pc, group = group, 
                        intgroup.df, names = colnames(x))
        attr(d, "percentVar") <- percentVar[1:pcs]
        rownames(d) = colnames(x)
        return(d) 
    }


# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting

PlotPCA = function(tab, cdata, ntop=500, intgroup, title='PCA', shape=1, color='black'){
    
    data = getPCA(tab,intgroup=intgroup, cdata=cdata, pcs=3, ntop=ntop)
    data$color = color
    data$shape = shape
    
    percentVar = round(100 * attr(data, "percentVar"))
    
    print(ggplot(data, aes(PC1, PC2, color=color, shape=shape)) +
              geom_point(size=6) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance")) +
              scale_size_manual(values=seq(5,10,length.out=3))+ ggtitle(title)) 
    
    print(ggplot(data, aes(PC1, PC3, color=color, shape=shape)) +
              geom_point(size=6) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC3: ",percentVar[3],"% variance")) +
              scale_size_manual(values=seq(5,10,length.out=3))+ ggtitle(title))
    
    print(ggplot(data, aes(PC2, PC3, color=color, shape=shape)) +
              geom_point(size=6) +
              xlab(paste0("PC2: ",percentVar[2],"% variance")) +
              ylab(paste0("PC3: ",percentVar[3],"% variance")) +
              scale_size_manual(values=seq(5,10,length.out=3))+ ggtitle(title)) 
    
}


# ------------------------------------------------------------------- #
# a function that gets the first n principal components for plotting