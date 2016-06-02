
# ----------------------------------------------------------------------------- #
# does tsne dimensionallity reduction on sequence kmer

calculate_kmers = function(seq.all, width=6){
    
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
    kmers.red = data.table(data.frame(lk))+1
    
    message('Normalization...')
    kmers.norm = kmers.red/rowSums(kmers.red)
    kmers.geo = log((kmers.norm)/exp(rowMeans(log(kmers.norm))))
    
    message('Returning...')
    return(kmers.geo)
    
}



