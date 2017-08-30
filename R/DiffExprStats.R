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

