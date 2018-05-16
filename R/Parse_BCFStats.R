path = 'Proseq_HEK_KO_0m_br1.filt.snps.stats.txt'

# ---------------------------------------------------------------------------- #
parse_BCFstats = function(path){

    library(stringr)
    sname = str_replace(basename(path),'.filt.snps.stats.txt','')
    s = scan(path, what='character', sep='\n')
    ind  = unique(str_replace(s, '\\t.+', ''))
    ind  = ind[!str_detect(ind,'#')]
    ind  = setdiff(ind, 'ID')
    wind = lapply(setNames(ind, ind), function(x)range(which(str_detect(s,paste0(x,'\\t')))))
    lout = list()
    for(i in names(wind)){

        message(i)
        d = s[Reduce(':',wind[[i]])]
        d = str_replace(d,'^\\s+#','')
        d = do.call(rbind, strsplit(d,'\\t'))
        d = data.frame(d)
        colnames(d) = str_replace(d[1,],'\\[.+\\]','')
        d = d[-1,-1]
        d = data.frame(d)

        if(i == 'SN')
            d$value = as.numeric(d$value)
        if(i == 'ST')
            d$count = as.numeric(d$count)

        if(i %in% setdiff(names(wind),c('SN','ST')))
            d = suppressWarnings(apply(d, 2, as.numeric))

        d = data.frame(sample = sname, as.data.frame(d))
        lout[[i]] = d
    }
    return(lout)
}

# ---------------------------------------------------------------------------- #
parse_BCFstats_Files = function(
    path,
    suffix = '.filt.snps.stats.txt'
){

    infiles = list.files(path, recursive=TRUE, pattern=suffix, full.names=TRUE)
    if(length(infiles) == 0)
        stop('There are no input files')

    lout   = suppressMessages(lapply(infiles, parse_BCFstats))
    snames = names(lout[[1]])
    dout   = lapply(setNames(snames,snames), function(x)
        do.call(rbind, lapply(lout, '[[', x))
    )
    return(dout)
}
