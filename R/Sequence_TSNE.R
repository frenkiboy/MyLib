
# ---------------------------------------------------------------------------- #
# does tsne dimensionallity reduction on sequence kmer

calculate_kmers = function(seq.all, width=6, normalize=TRUE, reverse.complement=TRUE, smooth=FALSE, alpha=.3){
    
    if(class(seq.all) != 'DNAStringSet')
        stop('ds needs to be DNAStringSet')
    
    
    ss = seq(1, length(seq.all), by=10000)
    ds = data.frame(x=c(1,ss[-1]+1), y=c(tail(ss,-1),length(seq.all)))
    set.seed(1)
    
    require(data.table)
    # ------------------------------------------------- #
    message('Kmers...')
    lk = list()
    for(i in 1:nrow(ds)){
        print(i)
        lk[[as.character(i)]] = data.table(oligonucleotideFrequency(seq.all[ds[i,1]:ds[i,2]], width=width, as.prob=FALSE, step=1))
    }
    kmers = rbindlist(lk)
    
    message('Revcomp...')
    if(reverse.complement){
        dk = data.frame(k1=colnames(kmers), k2=as.character(reverseComplement(DNAStringSet(colnames(kmers)))))
        dk = dk[!duplicated(apply(dk,1, function(x)paste(sort(x),collapse=':'))),]
        lk = list()
        for(i in 1:nrow(dk)){
            if(dk[i,1] != dk[i,2]){
                lk[[dk[i,1]]] = rowSums(kmers[,colnames(kmers) %in% dk[i,,drop=T],with=FALSE])
            }else{
                lk[[dk[i,1]]] = kmers[,colnames(kmers) == dk[i,1],with=FALSE]
            }
        }
        kmers = data.table(data.frame(lk))
    }
    
    message('Normalization...')
    kmers = kmers/rowSums(kmers)
    
    if(smooth){
        message('Smoothing...')
        if(width > 6)
            stop('with > 6 is not supported')
        
        if(reverse.complement)
            stop('revcomp is not supported')
        require(Biostrings)
        str = mkAllStrings(c('A','C','G','T'), width)
        A = as.matrix(stringDist(str, method='hamming'))
        A[A > 1] = 0
        A = A/colSums(A)
        
        kmers.0 = as.matrix(kmers)
        kmers.smooth = kmers.0
        for(i in 1:10){
            kmers.smooth = alpha*(kmers.smooth%*%A) + (1-alpha)*(kmers.0)
        }
        kmers = kmers.smooth
    }
    
    kmers = as.data.table(kmers)
    if(normalize)
        kmers = log((kmers+1)/exp(rowMeans(log(kmers+1))))
    
    message('Returning...')
    return(kmers)
    
}


# ---------------------------------------------------------------------------- #
plot_tsne = function(kmers, coldata, color=NULL, shape=NULL, outname, per=30, size=.8){
    
    require(cowplot)
    require(Rtsne)
    coldata = coldata[!duplicated(kmers),]
    kmers = kmers[!duplicated(kmers),]
    
    message('tsne...')
    ktsne = Rtsne(kmers, perplexity=per, dims=3)
    dat = data.frame(ktsne$Y, coldata)
    colnames(dat)[1:3] = c('X','Y','Z')
    
    message('Plotting...')
    # pdf(file.path(path.out, DateNamer(paste(outname,per,'pdf', sep='.'))),width=25, height=9)
    gl=list(
        g1=ggplot(dat, aes_string('X','Y', color=color, shape=shape)) + geom_point(size=size) +
            ggtitle(paste('X - Y')),
        g2=ggplot(dat, aes_string('Z','Y', color=color, shape=shape)) + geom_point(size=size) +
            ggtitle(paste('Z - Y')),
        g3=ggplot(dat, aes_string('X','Z', color=color, shape=shape)) + geom_point(size=size) +
            ggtitle(paste('X - Z'))
    )
    grobs = ggplotGrob(gl[[1]])$grobs
    gl = lapply(gl, function(x)x + theme(legend.position="none"))
    legend = grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]

    pg = plot_grid(plotlist=gl, nrow=1, ncol=3) 
    pg = plot_grid( pg, legend, rel_widths = c(9, .3))
    title = ggdraw() + draw_label(outname, fontface='bold')
    pg = plot_grid(title, pg, ncol=1, rel_heights=c(0.1, 1))
    print(pg)
    # dev.off()
    
}


# ---------------------------------------------------------------------------- #
run_kmerTSNE = function(seq, 
                        coldata,  
                        kmer.size = 4:7,
                        reverse.complement=TRUE,
                        normalize=TRUE,
                        smooth=FALSE,
                        perplexities = c(5,10,30, 50), 
                        alpha = .3,
                        path.out, 
                        outname, 
                        ncores=8, 
                        color=NULL, 
                        shape=NULL,
                        size=0.8){
    
 
    require(doMC)
    registerDoMC(ncores)
    foreach(kmer.width = kmer.size)%dopar%{
        
       
        
        if(normalize)
            outname=paste(outname, 'gnorm', sep='.')
        
        if(reverse.complement)
            outname=paste(outname, 'rc', sep='.')
        
        if(smooth)
            outname=paste(outname, 'sm', sep='.')
        
        kmers = calculate_kmers(seq, 
                                    kmer.width, 
                                    normalize=normalize, 
                                    reverse.complement=reverse.complement,
                                    smooth=smooth,
                                    alpha=alpha)
        var = apply(kmers, 2, var)
        dif = apply(kmers, 2, function(x)length(unique(x)))
        kmer.list = list(k.all     = kmers,
                         k.sel.var = kmers[, var > median(var), with=FALSE],
                         k.sel.10  = kmers[, order(-var),        with=FALSE][,1:min(10, ncol(kmers)), with=FALSE],
                         k.sel.100 = kmers[, order(-var),        with=FALSE][,1:min(100,ncol(kmers)), with=FALSE],
                         k.sel.500 = kmers[, order(-var),        with=FALSE][,1:min(500,ncol(kmers)), with=FALSE])
        
        if(!all(c(color, shape) %in% colnames(coldata)))
            stop('coldata is missing columns')
        
        pdf(file.path(path.out, DateNamer(paste(outname, kmer.width,'pdf', sep='.'))),width=25, height=8)
            for(per in perplexities){
                print(per)
                lapply(names(kmer.list), function(x){
                    plot_tsne(kmer.list[[x]],
                              coldata,
                              paste(outname, x,'tsne', kmer.width, per, sep='.'),
                              per=per, 
                              color = color,
                              shape = shape, 
                              size = size)
                    })
            }
        dev.off()
    }
}

