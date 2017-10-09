# ----------------------------------------------------------------------------- #
findRegions = function(v, sa=NULL, nw=-1, cw=1, gap=2000, gval=-100){


    require(dplyr)
    if(is.null(sa))
        sa = seq_along(v)

    if(gval > 0)
        stop('gval needs to be a negative value')
    if(nw > 0)
        stop('nw needs to be a negative value')
    if(cw <= 0)
        stop('cw needs to be a positive value')


    counter = 0
    pcnt = 0
    message('Starting...')
    rpos = list()
    for(s in sa){

        cval = as.vector(v[s])
        if(counter <1 & cval == 0){
            next()
        }
        pval = as.vector(v[s-1])

        if(counter == 0 & pval == 0 & cval > 0){
            lpos = list()
            lpos$start = s
            pcnt = 0
            ncnt = 0
            cw = cw
            nw = nw
            gap = gap
            gval = gval
        }
        if(cval == 0){

            ncnt = ncnt + 1
            if(ncnt >= gap)
                nw = gval

            counter = counter + nw
        }else if(cval > 0){
            counter = counter + cval*cw
            ncnt = 0
        }

        if(cval > 0){
            pcnt = pcnt + cval
        }

        if(counter <= 0){
            lpos$end = s-1
            lpos$score = pcnt
            rpos = c(rpos, list(data.frame(lpos)))
            counter=0
        }
        # message(paste('s',s,'cval',cval,'pval',pval,'cnt',counter,'pcnt',pcnt))
    }
    message('Returning...')
    dpos = do.call(rbind, rpos)
    return(dpos)
}

# ----------------------------------------------------------------------------- #
findRegionsRLE = function(v, nw=-1, cw=1, gap=2000, gval=-100){


    require(dplyr)


    if(gval > 0)
        stop('gval needs to be a negative value')
    if(nw > 0)
        stop('nw needs to be a negative value')
    if(cw <= 0)
        stop('cw needs to be a positive value')

    message('Starting...')


    vals = runValue(v)
    lens = runLength(v)
    if(vals[1] == 0){
        start = 2
    }else{
        start = 1
    }
    len = length(vals)
    sa = start:len

    counter = 0
    pcnt = 0
    rpos = list()
    for(s in sa){
        cat(round(s/len,2),'\r')

        cval = vals[s]
        clen = lens[s]
        pval = vals[s-1]

        if(counter == 0 & pval == 0 & cval > 0){
            lpos = list()
            lpos$start = sum(lens[1:(s-1)])+1
            pcnt = 0
            ncnt = 0
            nwu = nw
        }
        if(cval == 0){

            if(clen >=gap){
                nvec = rep(c(nwu,gval),times=c(gap,clen-gap))
            }else{
                nvec = rep(nwu,times=clen)
            }
            snvec = cumsum(nvec)
            sndiff = counter + snvec
            if(any(sndiff <= 0)){
                ecoord = min(which(sndiff <= 0))-1
                lpos$end = sum(lens[1:(s-1)])+ecoord
                lpos$score = pcnt
                rpos = c(rpos, list(data.frame(lpos)))
                counter=0
            }else{
                counter = counter + sum(nvec)
            }
        }else if(cval > 0){
            counter = counter + (cval*clen)*cw
            pcnt = pcnt + cval*clen
            ncnt = 0
        }

        # message(paste('s',s,'cval',cval,'pval',pval,'cnt',counter,'pcnt',pcnt))
    }
    message('Returning...')
    dpos = do.call(rbind, rpos)
    return(dpos)
}

findRegionsGenome = function(r, strand='+', nw=-1, cw=1, gap=2000, gval=-100){

    l.chr = list()
    chrs = setdiff(names(r),c('chrM','chrY','chrHSV1'))
    r = r[chrs]
    l.chr = foreach(v=r)%dopar%{
        if(strand == '-')
            v = rev(v)
        regs = findRegionsRLE(v, nw=nw, cw=cw, gap=gap, gval=gval)

        regs = data.frame(regs, strand=strand)
        if(strand == '-'){
            rs = length(v) - regs$end   + 1
            re = length(v) - regs$start + 1
            regs$start = rs
            regs$end   = re
        }
        regs
    }
    dregs = data.frame(chr=rep(chrs, times=sapply(l.chr, nrow)),do.call(rbind, l.chr))

    gpos = makeGRangesFromDataFrame(dregs, keep.extra.columns = TRUE)
    return(gpos)
}

# ----------------------------------------------------------------------------- #
# loops through bw files and finds antisense regions
Sample_FindRegion = function(bw.files, gtf, param=NULL, outpath, export.bw=TRUE, export.bed=TRUE, mpat='m.bw', 
                             normalize=TRUE){

    if(is.null(param))
        stop('please spcify the parameters')

    lregs=list()
    for(i in 1:length(bw.files)){

        bw.file = bw.files[i]
        bwname = str_replace(basename(bw.file),'.bw','')
        print(bwname)
        bw = import.bw(bw.file, as='GRanges')
        
        if(normalize){
            total = sum(as.numeric(bw$score))
            norm.fac = (10^(nchar(as.character(total))-1))/total
            bw$score = round(bw$score*(norm.fac),3)
        }

        strand = ifelse(str_detect(basename(bw.file),mpat),'-','+')

        bw = suppressWarnings(bw[countOverlaps(bw,gtf[strand(gtf) == strand])==0])
        cov = coverage(bw, weight=bw$score)


        regs = findRegionsGenome(cov,
                                 strand,
                                 nw=param$nw,
                                 cw=param$cw,
                                 gap=param$gap,
                                 gval=param$gval)
        lregs[[bwname]] = regs
        
        if(export.bw)
            export.bw(cov, file.path(outpath, DateNamer(paste(bwname, 'sub', 'bw', sep='.'))))
        
        if(export.bed){
            param.name = paste(names(param), param[1,], sep='.', collapse='_')
            write.table(as.data.frame(sort(regs))[,1:3], file.path(outpath, DateNamer(paste(bwname,param.name,'bed',sep='.'))),row.names=F,col.names=F,quote=F, sep='\t')
        }
    }
    return(lregs)

}


# ----------------------------------------------------------------------------- #
regionError = function(gpos, gtf.anti){

    fog = as.data.table(findOverlaps(gpos, gtf.anti))
    fog$ns = start(gpos)[fog$queryHits]
    fog$ne = end(gpos)[fog$queryHits]
    fog$ks = start(gtf.anti)[fog$subjectHits]
    fog$ke = end(gtf.anti)[fog$subjectHits]
    fog[,sdiff:=ks-ns]
    fog[,ediff:=ke-ne]
    fog[,adiff:= abs(sdiff)+abs(ediff)]
    fog[,jacc := round((pmin(ne,ke)-pmax(ns,ks))/ (pmax(ne,ke)-pmin(ns,ks)),2)]
    fog = arrange(fog, desc(adiff)) %>% filter(!duplicated(queryHits))
    return(fog)
}
