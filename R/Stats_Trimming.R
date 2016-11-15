
# ------------------------------------------------------------------------ #
# for BBduk
TrimmingStats_BBduk = function(path){

    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\n', quiet=TRUE)
    sind = grep('Input:',s)
    s = s[sind:(sind+3)]
    s = str_replace(s,'\t\t','\t')
    s = str_replace(s,':.+?\t','\t')
    s = str_replace(s,' reads.+?\t','\t')
    s = str_replace(s,' bases.+','')
    s = data.frame(do.call(rbind, strsplit(s,'\t')))
    colnames(s) = c('step','reads','bases')
	  s[,1] = c('Total','QTrimmed','ATrimmed','Result')
    ms = melt(s, id.vars='step')
    colnames(ms)[3] = 'cnts'
    ms$cnts = as.numeric(ms$cnts)
    ms$freq = ms$cnts
    ms$freq[1:4] = round(ms$freq[1:4]/ms$freq[1],2)
    ms$freq[5:8] = round(ms$freq[5:8]/ms$freq[5],2)

    sname = basename(path)
    sname = str_replace(sname,'.log','')
    sname = str_replace(sname,'.stat.+','')
    ms = data.table(sample = sname, ms)
    
    return(ms)
}

# ------------------------------------------------------------------------ #
# general function for getting mapping statistics
# output = data.table: sample, mapped, cnts
GetTrimmingStats = function(path, which.stats=NULL){

    require(stringr)
    require(data.table)


	if(is.null(which.stats))
		stop('Specify the trimmer')

	pattern = list(bbduk='log$')[[which.stats]]

	if(length(path) == 0)
            stop('There are no files')

	if(length(path) == 1){

    if(!file.exists(path))
      stop('The input directory does not exist')
    
	if(which.stats=='bbduk'){    
	    files = list.files(path, full.names=TRUE, pattern='log$', recursive=TRUE)
	    sl = sapply(files, function(x)length(scan(x, what='character',sep='\n', quiet=TRUE)))
	    if(!all(sl==16) & any(sl==16))
            stop('stat files got mixed')
	    if(all(sl == 16))
	        files = list.files(path, full.names=TRUE, pattern='stat', recursive=TRUE)
	}
		
	}else{
		files=path
	}
    # ---------------------------------- #
    message('Reading the annotation...')
    if(which.stats == 'bbduk')
        lstat = lapply(files, TrimmingStats_BBduk)

    dstat = rbindlist(lstat)
    return(dstat)
}



# ------------------------------------------------------------------------ #
# StatsPlotting function
plot_TrimmingStats = function(dstat, outpath=NULL, name='Trimming_Mapping', width=6, height=6){

    require(stringr)
    require(ggplot2)
    gsub = subset(dstat, step %in% c('Total','Result') & variable == 'reads')
    gsub$step = factor(gsub$step, levels=c('Total','Result'), ordered=TRUE)
    g.tot = ggplot(gsub, aes(x=sample, y=cnts, fill=step)) +
            geom_bar(stat='identity', position="dodge") +
            xlab('sample') + ylab('Number of reads') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Number of reads per sample')

    g.map = ggplot(subset(dstat, variable == 'reads' & step=='Result'), aes(x=sample, y=freq)) +
            geom_bar(stat='identity') +
            xlab('sample') + ylab('Percentage') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Percentage of remaining reads after filtering') + ylim(c(0,1))


    if(!is.null(outpath)){
        pdf(file.path(outpath, DateNamer(paste(name, 'pdf', sep='.'))), width=width, height=height)
            print(g.tot)
            print(g.map)
        dev.off()
    }
}
