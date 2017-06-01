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
get_DifferentialExpression = function(
    trans,
    bamfiles,
    coldata,
    design=NULL,
    nreads=5,
    nsamp=3,
    contlist=NULL,
    ignore.strand=FALSE,
    independent.filtering=TRUE,
    betaPrior=TRUE,
    preprocess.reads=NULL,
    singleEnd=TRUE,
    invertStrand=FALSE,
    merge_id = 'transcript_id',
    annotation=NULL,
    outpath,
    name,
    cnts.name=NULL,
    load=FALSE,
    lfc=1,
    padj=0.01
	){

    library(GenomicRanges)
    library(GenomicAlignments)
    library(DESeq2)
    library(sva)
    library(dplyr)
    source(file.path(lib.path, 'Annotate_Functions.R'), local=TRUE)
    source(file.path(lib.path, 'BamWorkers.R'), local=TRUE)
    source(file.path(lib.path, 'DifferentialExpression.R'), local=TRUE)
    source(file.path(lib.path, 'ScanLib.R'), local=TRUE)
    if(is.null(contlist))
        stop('Please specify the contrast list')

    # if(class(trans) == 'GRangesList'){
    #   utrans = unlist(trans)
    # }else{
    #   utrans = trans
    # }
    #
    # if(!any(id.col %in% colnames(values(utrans))))
    #     stop('id column is invalid')

    if(is.null(design))
        design = formula('~Factor')

    outfile=file.path(outpath, paste(name, 'Cnts', 'rds', sep='.'))
    if(load){
      message('Loading...')
      txhits=readRDS(outfile)
    }else{
      message('Summarize...')
      ranges=trans
      if(invertStrand)
	ranges = invertStrand(ranges)
      txhits = summarizeOverlaps(ranges, BamFileList(bamfiles),
                                 ignore.strand=ignore.strand,
                                 param=ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE)),
                                 preprocess.reads=preprocess.reads,
				 singleEnd=singleEnd)
      message('Saving...')
      saveRDS(txhits, outfile)
    }


    message('DES...')
    colData(txhits) = DataFrame(coldata)
    ass = assays(txhits)[[1]]
    ass = ass[rowSums(ass > nreads)>nsamp,]
    
    if(!is.null(cnts.name)){
        colnames(ass) = coldata[[cnts.name]]
    }else{
	colnames(ass) = BamName(bamfiles)
    } 
    dds = DESeqDataSetFromMatrix(ass, colData=coldata, design=design)
    des = DESeq(dds, parallel=FALSE, betaPrior=betaPrior)
    colnames(des) = colnames(ass)
    vsd = varianceStabilizingTransformation(des)
    cnts = as.data.frame(counts(des, normalized=TRUE))
    cnts$id = rownames(cnts)

    message('Results...')
    res = getResults(des, contlist, lfc=lfc, pval=padj,
                     independentFiltering=independent.filtering)
    means = getMeans.DESeqDataSet(des)
    message('Dat...')
    if(is.null(annotation)){
   	ann = Get_Annotation(trans)
    }else{
	ann = annotation
    }

    ann$id = ann[[merge_id]]

    dat = merge(res, means, by='id')
    dat = merge(dat, cnts, by='id')
    dat = merge(ann, dat, by='id')%>%
        mutate(id = NULL)



    return(list(trans=trans, txhits = txhits, des = des, vsd=vsd, res=res, dat=dat))
}
#
# # ---------------------------------------------------------------------------- #
# getAnnotation_GrangesList = function(gl){
#
#   ran = range(gl)
#   glu = unlist(gl)
#   tab = as.data.frame(unlist(ran))
#   tab$transcript_id = names(ran)
#   tab = merge(tab, unique(as.data.frame(values(glu))[,c('gene_id','transcript_id','gene_biotype')]), by='transcript_id')
#   tab$twidth = sum(width(gl))[tab$transcript_id]
#   return(tab)
#
# }


# ---------------------------------------------------------------------------- #

getResults_limma = function(fit, contrasts, lfc=1, pval=0.05, nres=1000000){
    require(data.table)
    message('Results... ')
    ltop = lapply(contrasts, function(x){
      top = topTable(fit, coef=x, number=nres)
      top = top[,c('logFC','adj.P.Val')]
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
  design = design[,match(rownames(contrast.matrix), colnames(design))]


  message('Fit... ')
  fit  = lmFit(eset, design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  message('eBayes... ')
  fit2 = eBayes(fit2, robust=TRUE)
  return(fit2)
}

# ---------------------------------------------------------------------------- #
get_limma_tab = function(expr, samps, lfc=1, padj=0.05){

  library(limma)
  lm = get_limma(expr, samps)
  cont = makeBinaryContrasts(unique(samps))
  res = getResults_limma(lm, cont, lfc=lfc, padj=padj)

	if(class(expr) == 'expressionSet'){

		dat = as(featureData(expr),'data.frame') %>%
		dplyr::select(1,9,10,11)
		colnames(dat) = str_replace(colnames(dat),' ','_')
		colnames(dat) = tolower(colnames(dat))
		means = getMeans(exprs(expr), samps, unique=FALSE)
		means$id = rownames(means)
  		tab = merge(dat, res, by='id')
	}else{
		means = getMeans(expr, samps, unique=FALSE)
		means$id = rownames(means)
		tab = res
	}

  tab = merge(tab, means, by='id')
  return(tab)
}
