# --------------------------------------------------------------- #
parseSalmonLog = function(file){
    
    s = scan(file, what='character', sep='\n',quiet=TRUE)
    name = basename(dirname(file))
    ind = which(grepl('jointLog',s))[4]-1
    s = s[ind]
    s = unlist(str_split(s,' '))
    s = s[c(5,12)]
    s = as.numeric(s)
    d = data.frame(name=name,
                   total =s[2],
                   mapped=s[1],
                   mapped.perc=round(s[1]/s[2],2))
    d
}

# --------------------------------------------------------------- #
parseSalmonStat = function(file){
    
    s = scan(file, what='character', sep='\n',quiet=TRUE)
    name = basename(dirname(file))
    message(name)
    s = s[4:5]
    s = sub('^.+: ','',s)
    s = as.numeric(s)
    d = data.frame(name = name,
                   strand.correct.percentage=round(s[1]/(s[1]+s[2]),3))
    d
}


# --------------------------------------------------------------- #
stats_Salmon = function(inpath){
    
    log.files  = list.files(inpath, full.names=TRUE, recursive=TRUE, pattern='log')
    log.files  = log.files[!grepl('quant',log.files)]
    stat.files = list.files(inpath, full.names=TRUE, recursive=TRUE, pattern='libFormatCounts.txt')

    logd = do.call(rbind, lapply(log.files, parseSalmonLog))
    statd = do.call(rbind, lapply(stat.files, parseSalmonStat))
    mstat = merge(logd, statd, by='name')
    return(mstat)
}