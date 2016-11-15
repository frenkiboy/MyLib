# ------------------------------------------------------------------------------ #
# functions for selecting differentially expressed genes
diffNum = function(res, nomark='None'){

    require(stringr)
    diff = res[,str_detect(colnames(res),'diff'),drop=FALSE]
    sum(rowSums(diff != nomark)>0)
}

diffTab = function(res){

    require(stringr)
    diff = res[,str_detect(colnames(res),'diff'),drop=FALSE]
    levs = sort(unique(unlist(lapply(diff,unique))))
    sapply(diff, function(x)table(factor(x,levels=levs)))
}

diffSel = function(res, nomark='None'){

    require(stringr)
    diff = res[,str_detect(colnames(res),'diff'),drop=FALSE]
    res[rowSums(diff != nomark)>0,]
}

diffAll = function(res, nomark='None'){

    require(stringr)
    diff = res[,str_detect(colnames(res),'diff'),drop=FALSE]
    res[rowSums(diff != nomark) == ncol(diff),]
}

diffCol = function(res,i, nomark='None'){

    require(stringr)
    res[res[,i] != nomark,]
}

diffMark = function(res, lfc, pval, log.col=NULL, pval.col=NULL, nomark='No'){

	diff = rep(nomark, nrow(res))

  if(is.null(log.col))
    log.col = which(colnames(res) == 'log2FoldChange')

  if(is.null(pval.col))
    pval.col = which(colnames(res) == 'padj')

	lfc.u  = which(res[,log.col] > lfc)
	lfc.d  = which(res[,log.col] < -lfc)
	pval.w = which(res[,pval.col] < pval)
	diff[intersect(lfc.u, pval.w)] = 'Up'
	diff[intersect(lfc.d, pval.w)] = 'Down'
	return(diff)
}

getResults = function(des, contrasts, lfc, pval, independentFiltering=FALSE){

    library(data.table)
    lres = list()
    for(i in 1:length(contrasts)){

        name = names(contrasts)[i]
        print(name)
        res = results(des, contrasts[[name]], independentFiltering=independentFiltering)
        res = res[,c('log2FoldChange','padj')]
        res$diff = diffMark(res, lfc, pval)
        colnames(res) = paste(colnames(res),name, sep='.')
        if(is.null(rownames(res)))
            rownames(res) = as.character(1:nrow(res))
        res = data.frame(id=rownames(res),res)
        lres[[name]] = data.table(res)
    }
    mres = MergeDataTable(lres, key='id',all=TRUE)
    return(mres)
}

# ------------------------------------------------------------------------------ #
# functions for making contrasts
makeBinaryContrasts = function(data,column='sample'){

    if(class(data) == 'data.frame')
        contrasts = expand.grid(sort(unique(data[[column]])), sort(unique(data[[column]])))

    if(class(data) == 'character')
        contrasts = expand.grid(sort(unique(data)), sort(unique(data)))

    if(class(data) == 'factor')
        contrasts = expand.grid(sort(levels(data)), sort(levels(data)))

    contrasts = subset(contrasts, Var1 != Var2)
    contrasts = contrasts[order(contrasts$Var1),]
    contrasts = contrasts[!duplicated(apply(contrasts, 1, function(x)paste(sort(x), collapse='-'))),]
    contrasts = contrasts[order(contrasts[,1]),]
    contrasts = with(contrasts, paste(Var1, Var2, sep='-'))
    return(contrasts)
}



# ------------------------------------------------------------------------------ #
getMeans = function(mat, factors, unique=TRUE){

    if(unique==TRUE){
        mean.mat = lapply(factors, function(x){
                                col_ind=which(str_detect(colnames(mat),x))
                                if(length(col_ind) == 1){
                                    return(mat[,col_ind])
                                }else{
                                    return(rowMeans(mat[,col_ind]))
                                }})
    }else{
        mean.mat = lapply(unique(factors), function(x){
            col_ind=which(str_detect(factors,x))
            if(length(col_ind) == 1){
                return(mat[,col_ind])
            }else{
                return(rowMeans(mat[,col_ind]))
            }})
        factors = unique(factors)

    }
    mean.mat = data.frame(mean.mat)
    colnames(mean.mat) = paste('mean', factors,sep='.')
    return(mean.mat)
}

getMeans.DESeqDataSet = function(dds, login=FALSE, logout=FALSE){

    if(login==FALSE)
      mat = log2(counts(dds, normalized=TRUE)+1)

    colnames(mat) = as.character(dds$Factor)
    dmeans = getMeans(mat, levels(dds$Factor))
    if(logout == FALSE)
        dmeans = data.frame(2^dmeans)

    dmeans = round(dmeans, 2)
    dmeans$id = rownames(dmeans)
    return(dmeans)
}

getMeans.VST = function(vst, logout=FALSE){

  mat = assays(vst)[[1]]
  colnames(mat) = as.character(vst$Factor)
  dmeans = getMeans(mat, levels(vst$Factor))
  if(logout == FALSE)
      dmeans = data.frame(2^dmeans)

  dmeans$id = rownames(dmeans)
  return(dmeans)
}


# ---------------------------------------------------------------------------- #
plotDESeqDiagnostics = function(dds, contrasts, outpath, name){

    require(ggplot2)
    d = data.frame('name'=paste(dds$Factor, dds$replicate),'sizeFactors'=sizeFactors(dds), Factor=dds$Factor)
    vsd = rlog(dds)
    message('Starting...')
    pdf(file.path(outpath,DateNamer(paste(name, 'DESeq.Diagnostics.pdf',sep='_'))), width=6, height=6)
        print(qplot(data=d, x=name, y=sizeFactors, color=Factor) + ggtitle('sizeFactors') + theme(axis.text.x = element_text(angle = 90, hjust = 1)))
        plotDispEsts(dds)
        for(i in contrasts){
            message('Doing...')
            plotMA(results(dds, contrast=i))
        }
        plotPCA(vsd, intgroup='Factor')
        plotSparsity(dds)
    dev.off()
    message('Finished!')
}


# ---------------------------------------------------------------------------- #

getResults_limma = function(fit, contrasts, lfc=1, pval=0.05){
    require(data.table)
    message('Results... ')
    ltop = lapply(contrasts, function(x){
      top = topTable(fit, coef=x, number=100000)
      top = top[,c(1,5)]
      top$diff = diffMark(top, lfc, pval, 1, 2)
      colnames(top)=paste(x,colnames(top),sep='.')
      top$id = rownames(top)
      data.table(top)
    })
  message('Merging... ')
  results = MergeDataTable(ltop, key='id')
  colnames(results) = str_replace(colnames(results),'-','_')
  return(results)
}

# ---------------------------------------------------------------------------- #
get_limma = function(eset, samps){

  message('Contrasts... ')
  cont=makeBinaryContrasts(samps)
  contrast.matrix = makeContrasts(contrasts=cont,
                                  levels=unique(samps))

  message('Design... ')
  design = model.matrix(~0+samps)
  colnames(design) = str_replace(colnames(design),'samps','')

  message('Fit... ')
  fit  = lmFit(eset, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  message('eBayes... ')
  fit2 = eBayes(fit2, robust=TRUE)
  return(fit2)
}
