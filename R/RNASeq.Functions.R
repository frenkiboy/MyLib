# INFO: function for analysis of rnaseq data
# DATE: 20.05.2011.
# AUTH: v+

# ---------------------------------------- #
MakeBinaryContrasts = function(v){
    
    v = unique(v)
    estatus = expand.grid(Var1=v, Var2=v)
    estatus = estatus[!duplicated(apply(estatus, 1, function(x)paste(sort(x), collapse=":"))),]
    estatus = subset(estatus, Var1 != Var2)
    estatus = Reduce(function(...)paste(..., sep='-'), estatus)
    return(estatus)
}    
    
	
	# {{2}}
	# does DESeq differential expression on a table of counts
	# for no replicates
# 	DiffExpNoRep = function(expr){
# 	
# 		library(DESeq)
# 		cds = newCountDataSet(expr, c('T','N'))
# 		cds = estimateSizeFactors(cds)
# 		cds = estimateVarianceFunctions(cds, method="blind" )
# 		d = nbinomTest(cds, 'T', 'N')
# 		d$padj[is.na(d$padj)] = 1
# 		return(d)
# 	}
# 	
	# one sample has replicates
# 	DiffExpOneRep = function(expr){
# 	
# 		library(DESeq)
# 		cds = newCountDataSet(expr, c('T','N'))
# 		cds = estimateSizeFactors(cds)
# 		cds = estimateVarianceFunctions(cds, method="blind" )
# 		d = nbinomTest(cds, 'T', 'N')
# 		d$padj[is.na(d$padj)] = 1
# 		return(d)
# 	}
	#/{{2}}
	
#/{1} FUNCTIONS
