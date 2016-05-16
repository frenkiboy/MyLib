
path.deepbind = '/home/vfranke/bin/Software/MotifDiscovery/deepbind'
run_deepbind = function(seqs, type='TF'){

    lseq = lapply(1:length(seqs), function(x){
      d = DNAStringSet(as.character(Views(seqs[[x]], IRanges(seq(1, length(seqs[[x]]), 50), width=50))))
      names(d) = paste('s',x,1:length(d),sep='_')
      d
    })
    dseq = do.call(c, lseq)

    require(Biostrings)
    id.tab = read.table(file.path(path.deepbind,'db/db.tsv'), header=TRUE, sep='\t')
    id.tab = id.tab[id.tab$Type == type & id.tab$Labels != 'deprecated',]
    ids.file = file.path('~/Tmp/idfile.ids')
    write.table(id.tab$ID, ids.file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')


    rnumb = paste(paste(sample(1:10), collapse=''), 'fa', sep='.')
    seqs.file = file.path('~/Tmp',rnumb)
    writeXStringSet(dseq, seqs.file)

    deepbind = file.path(path.deepbind, 'deepbind')
    command = paste(deepbind,'--echo',
                    ids.file,
                    seqs.file)
    out = system(command, intern=TRUE)
    out = do.call(rbind, strsplit(out, '\\t'))
    dout = data.frame(apply(out[-1,], 2, as.numeric))
    colnames(dout) = out[1,]
    return(list(hits=dout, annot=id.tab))


}

annotate_deepbind = function(deep){

    require(reshape)
    mhits = melt(data.frame(seq = 1:nrow(deep$hits),deep$hits), id='seq')
    mhits = subset(mhits, value > 0)
    mannot = merge(deep$annot[,1:4], mhits, by.x='ID', by.y='variable')
    mannot = mannot[order(-mannot$value),]
    mannot = mannot[!duplicated(with(mannot, paste(Protein, seq))),]
    return(mannot)
}
