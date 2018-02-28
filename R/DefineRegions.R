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

# ---------------------------------------------------------------------------- #
# loops through bw files and finds antisense regions
Sample_FindRegion = function(
  bw.files,
  gtf        = NULL,
  param      = NULL,
  outpath,
  export.bw  = TRUE,
  export.bed = TRUE,
  mpat       = 'm.bw',
  normalize  = TRUE,
  subset.chr = NULL,
  strand     = TRUE,
  lower      = 0,
  upper      = 'max'
){

    source(file.path(lib.path, 'ScanLib.R'), local=TRUE)
    library(rtracklayer)
    if(is.null(param))
        stop('please specify the parameters')

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
        if(!is.null(subset.chr)){
          seqlevels(bw,  pruning.mode='coarse') = intersect(seqlevels(bw), subset.chr)
        }

        file.strand = ifelse(str_detect(basename(bw.file),mpat),'-','+')

        if(!is.null(gtf))
          bw = suppressWarnings(bw[countOverlaps(bw,gtf[strand(gtf) == file.strand])==0])

        cov = coverage(bw, weight=bw$score)

        message('Finding Regions ...')
            regs_raw = findRegionsGenome(
                                     cov,
                                     file.strand,
                                     nw   = param$nw,
                                     cw   = param$cw,
                                     gap  = param$gap,
                                     gval = param$gval)

        message('Defining Regions ...')
            regs_def = DefineRegionBorders(
                                regs_raw,
                                cov,
                                down   = param$down,
                                up     = param$up,
                                strand = strand,
                                lower  = lower,
                                upper  = upper)


        lregs[[bwname]]$regs_raw = regs_raw
        lregs[[bwname]]$regs_def = regs_def

        if(export.bw)
            export.bw(cov, file.path(outpath, DateNamer(paste(bwname, 'sub', 'bw', sep='.'))))

        if(export.bed){
            message('Export BED ...')
                param.name = paste(names(param), param[1,], sep='.', collapse='_')
                write.table(as.data.frame(sort(regs_raw))[,1:3], file.path(outpath, DateNamer(paste(bwname,param.name,'raw.bed',sep='.'))),row.names=F,col.names=F,quote=F, sep='\t')
                write.table(as.data.frame(sort(regs_def))[,1:3], file.path(outpath, DateNamer(paste(bwname,param.name,'def.bed',sep='.'))),row.names=F,col.names=F,quote=F, sep='\t')
        }
    }
    return(lregs)

}

# ---------------------------------------------------------------------------- #
# Takes a GRanges and a RleList and defines the regions that contain the bulk of the coverage
# RLEs are not stranded so you have to run the function for the + and minus strand separately
DefineRegionBorders = function(g, r, down=0.1, up=0.9, strand=FALSE, lower=0, upper='max'){

    # gets the chromosome names
    g$'_ind' = 1:length(g)
    g = sort(g)
    chrs = unique(as.character(seqnames(g)))


    # whether the region reduction should be strand oriented
    if(!strand){
        gsrl = as(g, 'RangesList')
        lregs = list()
        lregs = foreach(chr = chrs)%dopar%{
            v = Views(r[chr], gsrl[chr])
            va = viewApply(v[[chr]], function(x)GetRegs(x, down=down, up=up, strand='*', lower=lower, upper=upper), simplify=FALSE)
            regs = as.data.frame(do.call(rbind, va))
            return(regs)
        }

        dregs = do.call(rbind, lregs)
        setnames(dregs, 1:2, c('rstart','rend'))
        dregs$width=with(dregs, rend-rstart+1)
        dregs$gwidth = width(g)

        if(with(dregs, sum(width > gwidth)) != 0)
            stop('Final widths larger then starting widths')

        if(any(dregs$width < 0))
            stop('Some regions have width 0')

        end(g)   = start(g) + dregs$rend
        start(g) = start(g) + dregs$rstart
        ds = g

    }else{
        ls = list()

        strand = as.character(unique(strand(g)))
        for(s in strand){
            message(s)
            gs = g[strand(g) == s]
            if(any(width(gs) == 1))
                warning('Some regions have width 1')

            gsrl = as(gs, 'RangesList')
            lregs = foreach(chr = chrs)%dopar%{
                v = Views(r[chr], gsrl[chr])
                va = viewApply(v[[chr]], function(x)GetRegs(x, down=down, up=up, strand=s, lower=lower, upper=upper), simplify=FALSE)
                regs = as.data.frame(do.call(rbind, va))

                return(regs)
            }

            dregs = do.call(rbind, lregs)
            setnames(dregs, 1:2, c('rstart','rend'))
            dregs$width=with(dregs, rend-rstart+1)
            dregs$gwidth = width(gs)
            inf.ind = is.infinite(dregs$rstart) | is.infinite(dregs$rend)
            dregs$rstart[inf.ind] = 1
            dregs$rend[inf.ind] = 1
            dregs$width[inf.ind] = 1

            if(with(dregs, sum(width > gwidth)) != 0)
                stop('Final widths larger then starting widths')

            if(any(dregs$width < 0))
                stop('Some regions have width 0')

            end(gs)   = start(gs) + dregs$rend
            start(gs) = start(gs) + dregs$rstart
            gs = gs[dregs$width > 1]
            ls[[s]] = gs
        }
        ds = unlist(GRangesList(ls), use.names=FALSE)

    }
    ds = ds[order(ds$'_ind')]
    ds$'_ind' = NULL
    return(ds)
}


GetRegs = function(x, down=0.1, up=0.9, strand='*', lower=0, upper='max'){

    if(!strand %in% c('+','-','*'))
        stop('strand is not an allowed value')

    v = as.vector(x)
    if(lower > 0)
        v[v < lower] = 0

    if(upper != 'max')
        v[v > max] = max

    if(strand == '+' | strand == '*'){
        cs = cumsum(v)
        cs = cs/max(cs)
        reg = c(min(which(cs >= down)), max(which(cs <= up)))
    }
    if(strand == '-'){
        cs = cumsum(rev(v))
        cs = cs/max(cs)
        reg = c(min(which(cs >= down)), max(which(cs <= up)))
        reg = length(v) - rev(reg) + 1
    }

    if(all(v == 0))
        reg=c(0,0)

    return(reg)
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
