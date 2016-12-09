# ---------------------------------------------------------------------------- #
define_Clusters = function(bamfiles, path.rdata, mirna, gtf, width=100,
                           nworkers=32,
                           cnt.ind=5,
                           samp.ind=5,
                           clust.pattern=NULL,
                           new=FALSE){
    
    
    bamfiles = select_snRNA_Bamfiles()
    if(is.null(clust.pattern))
        clust.pattern = paste('si',samp.ind,'ci',cnt.ind, sep='.')
    
    if(new){
        message('Prepareing reads...')
        prepare_snRNA_Reads(bamfiles=bamfiles,
                            path.rdata=path.rdata,
                            mirna=mirna,
                            width=width,
                            nworkers=nworkers)
        
        
        message('Creating clusters...')
        define_snRNA_Clusters(path.rdata=path.rdata,
                              width=width,
                              cnt.ind=cnt.ind,
                              samp.ind=samp.ind)
    }
    
    clust.file = sort(list.files(path.rdata, pattern=clust.pattern, full.names=TRUE))
    if(!file.exists(clust.file))
        stop('The cluster file does not exist')
    clust = readRDS(clust.file)[1]
    
    # clust$gene_annot = suppressWarnings(AnnotateRanges(clust, lannot, type='precedence'))
    clust$rclust = annotate_Antisense(clust$rclust, gtf)
    clust$rclust$mirbase = GetAnnotOverlaps(clust$rclust, mirna, 'Name')
    return(clust)
    
}

# ---------------------------------------------------------------------------- #
prepare_snRNA_Reads = function(bamfiles,
                               path.rdata,
                               mirna,
                               width=100,
                               clust.width.bord = c(21,26,30),
                               nworkers
                               
){
    
    donefiles = list.files(path.rdata, pattern='Reads.rds')
    require(doMC)
    registerDoMC(nworkers)
    foreach(i = 1:length(bamfiles))%dopar%{
        
        bamfile = bamfiles[i]
        name = BamName(bamfile)
        if(name %in% donefiles)
            return()
        
        name = str_replace(name, '_genome', '')
        print(name)
        
        reads = readGAlignments(bamfile, use.names=TRUE)
        seqlevels(reads, force=TRUE) = setdiff(seqlevels(reads), 'chrM')
        g = granges(reads)
        g$name = paste(name, names(g), sep='.')
        tab = data.table(name=names(reads))
        tab = tab[,.N, by=name]
        tab$width = qwidth(reads)[match(tab$name, names(reads))]
        tab$short = cut(tab$width, breaks=c(min(tab$width), clust.width.bord, max(tab$width)), include.lowest=TRUE, right=FALSE)
        tab$short = factor(tab$short)
        tab$uniq = ifelse(tab$N==1, 'Uniq', 'Mult')
        tab$set = name
        tab[,rname := paste(set, name, sep='.')]
        wtab = tab[,.N,by=c('set','width')]
        
        clust = reduce(resize(g, width=width(g)+width, fix='center'))
        clust = resize(clust, width= width(clust)-width, fix='center')
        
        fo = data.table(as.matrix(findOverlaps(clust, g)))
        fo$rname = g$name[fo$subjectHits]
        fo$set    = tab$set[match(fo$rname, tab$rname)]
        fo$uniq   = tab$uniq[match(fo$rname, tab$rname)]
        fo$short   = tab$short[match(fo$rname, tab$rname)]
        fo$width   = tab$width[match(fo$rname, tab$rname)]
        
        foc = fo[,length(unique(subjectHits)),by=c('queryHits','set','uniq','short')]
        fos = dcast.data.table(foc, queryHits~set+uniq+short,  value.var='V1', fun=sum)
        fos[is.na(fos)] = 0
        
        mirfo = data.table(as.matrix(findOverlaps(mirna, g)))
        mirfo$rname = g$name[mirfo$subjectHits]
        mirfo$set    = tab$set[match(mirfo$rname, tab$rname)]
        mirfo$uniq   = tab$uniq[match(mirfo$rname, tab$rname)]
        mirfo$short   = tab$short[match(mirfo$rname, tab$rname)]
        mirfo$width   = tab$width[match(mirfo$rname, tab$rname)]
        
        mirfoc = mirfo[,length(unique(subjectHits)),by=c('queryHits','set','uniq','short')]
        mirfos = dcast.data.table(mirfoc, queryHits~set+uniq+short,  value.var='V1', fun=sum)
        mirfos[is.na(mirfos)] = 0
        
        
        l=list(fos=fos, clust=clust, mirfos=mirfos, wtab=wtab)
        saveRDS(l, file = file.path(path.rdata, paste(name, 'Reads.rds', sep='.')))
    }
}


# ---------------------------------------------------------------------------- #
define_Clusters  = function(path.rdata, width=100, cnt.ind, samp.ind){
    
    rdata.files = list.files(path.rdata, full.name=TRUE, pattern='Reads.rds')
    lbam = lapply(rdata.files, readRDS)
    names(lbam) = str_replace(rdata.files,'.Reads.rds','')
    
    r = reduce(unlist(GRangesList(lapply(lbam, '[[', 'clust'))))
    rclust = suppressWarnings(reduce(resize(r, width=width(r)+width, fix='center')))
    rclust = resize(rclust, width= width(rclust)-width, fix='center')
    
    lr = list()
    lr = foreach(i = 1:length(lbam), .inorder=FALSE)%dopar%{
        
        print(i)
        set = lbam[[i]]
        clust = set$clust
        fos   = set$fos
        
        fo = data.table(as.matrix(findOverlaps(clust,rclust)))
        fo = merge(fo, fos, by='queryHits',all=TRUE)
        fo$queryHits=NULL
        fo = fo[,lapply(.SD, function(x)sum(as.numeric(x))), by='subjectHits']
        return(fo)
    }
    m = MergeDataTable(lr, 'subjectHits')
    m = m[, lapply(.SD, function(x){x[is.na(x)] = 0;x})]
    m = format_Colnames(m)
    m.ind = filterCluster(m, cnt.ind, samp.ind)
    
    saveRDS(list(rclust=rclust,m=m), file=file.path(path.rdata, DateNamer(paste('snRNA_Cluster_All','rds',sep='.'))))
    saveRDS(list(rclust=rclust[m.ind], m=m[m.ind,]), file=file.path(path.rdata, DateNamer(paste('snRNA_Cluster','si',samp.ind,'ci',cnt.ind,'rds',sep='.'))))
}
# ---------------------------------------------------------------------------- #
filterCluster = function(m, cnt.ind=5, samp.ind=5){
    
    m.uniq = dplyr::select(m, contains('Uniq'))
    r.size = unique(str_replace(colnames(m.uniq),'^.+_',''))
    m.ind = lapply(r.size, function(x)rowSums(dplyr::select(m.uniq, contains(x))>cnt.ind)>samp.ind)
    m.ind = Reduce('|', m.ind)
    return(m.ind)
}

