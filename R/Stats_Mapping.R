
# ------------------------------------------------------------------------ #
# for bowtie
MappingStats_Bowtie = function(path){

    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\t', quiet=TRUE)
    s = str_replace(s,'^.+: ','')
    s = str_replace(s,'Reported ','')
    s = str_replace(s,' .+','')
    s = as.numeric(s)
    d = data.table(sample= str_replace(basename(path),'.log',''),
                   mapped= c('reads.total','map.uniq','map.mult','map.disc','map.total'),
                   cnts  = c(s[1:4],s[2]+s[3]))

    d[,freq := round(cnts/cnts[1],3)]
    return(d)
}

# ------------------------------------------------------------------------ #
# for bowtie
MappingStats_STAR = function(path){

    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\n', quiet=TRUE)
    s = as.numeric(str_replace(s[c(5,8,23,25)],'^.+\\t',''))
    d = data.table(sample= str_replace(basename(path),'Log.final.out',''),
                   mapped= c('reads.total','map.uniq','map.mult','map.disc','map.total'),
                   cnts  = c(s,s[2]+s[3]))

    d[,freq := round(cnts/cnts[1],3)]
    return(d)
}


# ------------------------------------------------------------------------ #
# for bowtie
MappingStats_Salmon = function(path){

    require(stringr)
    require(data.table)
    s = scan(path, what='character', sep='\t', quiet=TRUE)
    s = str_replace(s,'^.+: ','')
    s = str_replace(s,'Reported ','')
    s = str_replace(s,' .+','')
    s = as.numeric(s)
    d = data.table(sample= str_replace(basename(path),'.log',''),
                   mapped= c('reads.total','map.uniq','map.mult','map.disc','map.total'),
                   cnts  = c(s[1:4],s[2]+s[3]))

    d[,freq := round(cnts/cnts[1],3)]
    return(d)
}


# ------------------------------------------------------------------------ #
# general function for getting mapping statistics
# output = data.table: sample, mapped, cnts
GetMappingStats = function(path, which.stats=NULL){

    require(stringr)
    require(data.table)

	if(is.null(which.stats))
		stop('Specify the mapper')

  which.list = list(bowtie='log$', star='Log.final.out$', STAR='Log.final.out$')

  if(!which.stats %in% names(which.list))
    stop('which.stats not correct')
  pattern = which.list[[which.stats]]


	if(length(path) == 0)
            stop('There are no files')

	if(length(path) == 1){
		files = list.files(path, full.names=TRUE, pattern=pattern, recursive=TRUE)
	}else{
		files=path
	}
    # ---------------------------------- #
    message('Reading the annotation...')
    if(which.stats == 'bowtie')
        lstat = lapply(files, MappingStats_Bowtie)

    if(which.stats %in% c('star','STAR'))
        lstat = lapply(files, MappingStats_STAR)

    dstat = rbindlist(lstat)
    return(dstat)
}


# ------------------------------------------------------------------------ #
# StatsPlotting function
plot_MappingStats = function(dstat, outpath=NULL, name='Stats_Mapping', width=6, height=6){

    require(stringr)
    require(ggplot2)
    gsub = subset(dstat, str_detect(mapped, 'total'))
    gsub$mapped = factor(gsub$mapped, levels=c('reads.total','map.total'), ordered=TRUE)
    g.tot = ggplot(gsub, aes(x=sample, y=cnts, fill=mapped)) +
            geom_bar(stat='identity', position="dodge") +
            xlab('sample') + ylab('Number of reads') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Number of reads per sample')

    g.map = ggplot(subset(dstat, !str_detect(mapped, 'total')  ), aes(x=sample, y=freq, fill=mapped)) +
            geom_bar(stat='identity') +
            xlab('sample') + ylab('Number of reads') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle('Number of reads per sample')


    if(!is.null(outpath)){
        pdf(file.path(outpath, DateNamer(paste(name, 'pdf', sep='.'))), width=width, height=height)
            print(g.tot)
            print(g.map)
        dev.off()
    }
}


# ------------------------------------------------------------------------ #
list_star.logfiles = function(path){

	list.files(path, full.names=TRUE, recursive=TRUE, pattern='Log.final.out')
}
