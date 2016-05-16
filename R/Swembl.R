# -------------------------------------------------------------------------------- #
# P - Penalty per base pair applied after gap specifed above
# p - Penalty applied to a gap of one base pair
# d - Gap at which penalty per base pair increases
# x - Penalty factor for reference reads
# f - Fragment length for extending fragments
# t - threshold value for number of sample reads starting at one base pair
# i - input file
# o - output file


run_Swembl = function(name, param.tab, remove=TRUE){
    
    require(data.table)
    require(GenomicRanges)
    swembl = '~/bin/Software/PeakCalling/SWEMBL.3.3.1/SWEMBL'
    outname = paste(name,with(subset(param.tab, !param %in% c('i','o')), paste(paste0(param,value),collapse='.')),'swembl', sep='.')
    param.tab$value[param.tab$param == 'o'] = file.path(param.tab$value[param.tab$param == 'o'], outname)
    
    params = with(param.tab, paste(paste0('-',param), value, collapse=' '))
    command = paste(swembl, params,'-C','2> /dev/null')
    message('Running swembl ...')
    system(command, wait=TRUE, ignore.stdout=TRUE)
    r = read.table(subset(param.tab, param=='o')$value, header=TRUE, sep='\t')
    colnames(r)[1:3] = c('chr','start','end')
	
	if(remove)
		file.remove(subset(param.tab, param=='o')$value)
		
    if(nrow(r) > 0){
		message('Returning ...')
        g = makeGRangesFromDataFrame(r, keep.extra.columns=TRUE)
		return(g)
    }
}


test_Swembl = function(name, infile, outpath, params, nthreads=16, remove=TRUE){
    
    options(scipen=100)
    library(doMC)
    registerDoMC(nthreads)
    param.exp = do.call(expand.grid,params)
    lg = list()
    lg = foreach(i = 1:nrow(param.exp))%dopar%{
    
        message(round(i/nrow(param.exp),2))
        param.tab = data.frame(param=colnames(param.exp),
                               value=as.character(param.exp[i,]))
        param.tab = rbind(param.tab, data.frame(param=c('i','o'), value=c(infile, outpath)))
        g = run_Swembl(name, param.tab, remove=remove)
        return(g)
    }
	param.n = Reduce(function(x,y)paste(x,y,sep='.'),param.exp)
	names(lg) = param.n
	return(lg)
}
