
# ---------------------------------------------------------------------------- #
### garbage collect function
collect.garbage = function()
{
  while(gc()[2,4] != gc()[2,4]){}
}


### makes a tilling window for the given chromosomes
MakeTillingWindows = function(w, chrlen){

	wins = lapply(chrlen, function(x)seq(1,x-w,w))
	g.l = lapply(seq(along=wins), function(x)GRanges(seqnames = rep(names(chrlen)[x], times=length(wins[[x]])), IRanges(start=wins[[x]], width=w), seqlengths=chrlen))
	g.l = GRangesList(g.l)
	names(g.l) = names(wins)
	g.l
}

### makes a tilling window for the given chromosomes
MakeSlidingWindows = function(w, s, chrlen){

	wins = lapply(chrlen, function(x)seq(1,x-w,s))


	g.l = lapply(seq(along=wins), function(x)GRanges(seqnames = rep(names(chrlen)[x], times=length(wins[[x]])), IRanges(start=wins[[x]], width=w), seqlengths=chrlen))
	g.l = GRangesList(g.l)
	g.l
}


# ---------------------------------------------- #
### removes random chromosomes
RandomRemover = function(x){
	if (length(grep(pattern = "random", x[,1])) > 0){
		x = x[-c(grep(pattern = "random", x[,1])),]
		return(x)
	}else{
		return(x)
	}
}

### correts the ranges of regions that fall of the chromosomes
RegionCorrector = function(bed, seqlen){
	if(!identical(sort(unique(as.character(bed[,1]))), sort(unique(names(seqlen)))))
		stop('Provided seqlengths do not match with the names in the bed file')

	cat('Correcting bed...\n')
	for(chr in names(seqlen)){
		ind.l = bed[,1] == chr & bed[,3] > seqlen[chr]
		ind.s = bed[,1] == chr & bed[,2] < 1
		cat(chr,'l:', sum(ind.l),'s:',sum(ind.s),'\n')
		bed[ind.l,3] = seqlen[chr]
		bed[ind.s,2] = 1

	}
	return(bed)
}

### removes regions that fall out of chromosomes
RegionCleaner = function(tab, genome){

	l.tab=list()
	for(i in seqnames(genome)){
		print(i)
		l.tab[[i]] = tab[tab[,1] == i & tab[,3] <= seqlengths(genome)[i],]

	}
	return(do.call(rbind, l.tab))
}


# ---------------------------------------------------------------------------- #
### does a sliding window over a vector
TillingWindow = function(x, windowsize=3, step=2, floor=T){

	idx1 = seq(1,length(x),by=step)
	idx2 = idx1 + windowsize
	idx2[idx2>(length(x)+1)] = length(x)+1
	cx = c(0,cumsum(x))
	r=(cx[idx2]-cx[idx1])/windowsize
	if(floor == T)
		floor(r)
	return(r)
}

### makes a tilling window of maximum values
TillingWindowMax = function(x, windowsize){

	dif = length(x) %% windowsize
	m = matrix(rev(rev(x)[-(0:dif)]), ncol=windowsize)
	maxs = apply(m, 1, max)
	return(as.vector(maxs))
}

# ---------------------------------------------------------------------------------------- #
### construct the current date
Dater = function(){
	format(Sys.Date(), format="%y%m%d")

}

### construct a filename using the current date
DateNamer = function(name, sep='.'){
    return(paste(Dater(), name, sep=sep))

}

### gets colors for a factor variable
GetColors = function(n) {

	black = "#000000"
	if (n <= 9) {
		library(RColorBrewer)
		c(black,brewer.pal(n-1, "Set2"))
	} else {
		c(black,hcl(h=seq(0,(n-2)/(n-1),length=n-1)*360,c=100,l=65,fixup=TRUE))
	}
}


### sets the colnames of the first three columns of a data frame to be chr start end
Colnamer = function(d) {

	colnames(d)[1:3] = c('chr','start','end')
	return(d)
}

# ---------------------------------------------------------------------------- #
### Merges a list of data tables given their ids
MergeDataTable = function(l, key=NULL, all=TRUE, fill=NULL){


	if(is.null(key))
		stop('You have to give a key')

	l = lapply(l, function(x)setkeyv(x, cols=key))
	r = Reduce(function(x,y)merge(x,y, all=all,by=key), l)
    if(!is.null(fill))
        r[is.na(r)] = fill
	return(r)
}

# ---------------------------------------------------------------------------- #
### given two ranges - gets the overlapping data table
GetOverlaps = function(reg1,reg2, colname=NULL, ignore.strand=FALSE, funct='sum'){

    fun = match.fun(funct)
	if(is.null(colname)){
		cat('Using the default colname...\n')
		values(reg2)$weight=1
		colname=weight
	}
	require(data.table)
	fo = dtfindOverlaps(reg1,reg2, ignore.strand=ignore.strand)
	fo$weight = values(reg2)[[colname]][fo$subjectHits]
	fo = fo[,fun(weight), by=queryHits]
	v = rep(0, length(reg1))
	v[fo$queryHits] = fo$V1
	return(v)
}

# ---------------------------------------------------------------------------- #
GetAnnotOverlaps = function(reg1,reg2, colname=NULL, ignore.strand=FALSE, null.fac='None', sep=':'){

    if(is.null(values(reg2)[[colname]]))
        stop('Please specify a colname')

    fo = dtfindOverlaps(reg1,reg2, ignore.strand=ignore.strand)
    fo$colname = values(reg2)[[colname]][fo$subjectHits]
    fo = fo[,paste(unique(colname), collapse=sep), by=queryHits]
    v = rep(null.fac, length(reg1))
    v[fo$queryHits] = fo$V1
    return(v)

}

# ---------------------------------------------------------------------------- #
dtfindOverlaps = function(reg1, reg2, ignore.strand=FALSE){
		require(data.table)
		require(GenomicRanges)
		as.data.table(findOverlaps(reg1,reg2, ignore.strand=ignore.strand))
}


# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
GCoords = function(x, cor.sep='-', chr.sep=':'){
		paste(as.character(seqnames(x)), paste(start(x), end(x), sep=cor.sep), sep=chr.sep)
}



# ---------------------------------------------------------------------------- #
WriteRegs = function(g, outpath, name=NULL, strand=TRUE, sname=c('-'='m','+'='p','*'='n')){

    if(is.null(name))
    	stop('Please specify the name')


    names(g) = NULL
    g = sort(g)
    d = as.data.frame(g)
    if(!strand)
    	d$strand = '*'

    for(i in unique(d$strand)){
    	s = sname[i]
    	print(s)
    	write.table(d[d$strand == i,1:3], file.path(outpath, DateNamer(paste(name, s, 'bed', sep='.'))), row.names=F, col.names=F, quote=F, sep='\t')
    }
}

# ---------------------------------------------------------------------------- #
flipStrand = function(g){
	levels(strand(g)) = c('-','+','*')
	return(g)
}


# ---------------------------------------------------------------------------- #
chrSeqlevels = function(g){
	seqlevels(g) = paste0('chr',seqlevels(g))
	return(g)
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

				gs = g[strand(g) == s]
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


# ---------------------------------------------------------------------------- #
# takes the log of of the difference of two variables
difflog = function(x, logbase=2, zinf=TRUE){

    mind = x<0
    x = log(abs(x),logbase)
    x[which(mind)] = x[which(mind)]*-1
    if(zinf)
        x[is.infinite(x)] = 0
    x
}

# ---------------------------------------------------------------------------- #
# splits a character vector and returns a data frame with numbered elements
split.numlist = function(v, sep=';', id.name='id', num.name='num'){

    numlist(strsplit(v, split=sep, id.name, num.name))
}


#given a list of elements returns a data frame with numbered elements
numlist = function(l, id.name='id', num.name='num'){

    d = data.frame(unlist(l), rep(1:length(l) ,times=sapply(l, length)))
    names(d) = c(id.name, num.name)
    return(d)
}


# ---------------------------------------------------------------------------- #
# orders the Ensembl type seqnames
orderSeqnames = function(g){

    g = chrSeqlevels(g)
    g = keepStandardChromosomes(g)
    g = unlist(g)
    values(g) = NULL
    g

}

# ---------------------------------------------------------------------------- #
# gets the uniqe names given a character vector
Uniquer = function(v){

    require(data.table)
    d = data.table(id=v)
    d[,key:=1:.N,by=id]
    d[,key:=paste(id,key,sep='.')]
    d$key

}


# ---------------------------------------------------------------------------- #
markExonNumber = function(g, id='transcript_id'){

  library(data.table)
  d = data.table(strand=as.character(strand(g)), id = values(g)[[id]])
  d[d$strand == '+' , `:=`( COUNT = .N , IDX = 1:.N ) , by = id[strand == '+']]
  d[d$strand == '-' , `:=`( COUNT = .N , IDX = .N:.1) , by = id[strand == '-']]
  g$exon_number= d$IDX
  g$exon_count = d$COUNT
  g$exon_first = g$exon_number == 1
  g$exon_last  = g$exon_number == g$exon_count
  return(g)
}



# ---------------------------------------------------------------------------- #
GRangesTodata.frame = function(g){

  sind = which(sapply(values(g), class) == 'list')
  if(length(sind) > 0)
  	g = g[,-sind]
  dind = which(sapply(values(g), class) == 'DataFrame')
  if(length(dind)>0){
    df = do.call(cbind, lapply(dind, function(x)as.data.frame(values(g)[,x])))
    ann = cbind(as.data.frame(g[,-dind]), df)
  }else{
    ann = as.data.frame(g)
  }
  return(ann)
}

# ---------------------------------------------------------------------------- #
# show method for the DataFrame as DataFrame object
# suppressPackageStartupMessages(library(S4Vectors))
# setMethod("showAsCell",signature("DataFrame"),
#           function(object){
#               cnams = paste(colnames(object),collapse=':')
#               if(nchar(cnams) > 10)
#                   cnams = paste0(substring(cnams, 1,10),'...')
#               rep(cnams, nrow(object))
#           })


#####--------------------/FUNCTIONS/---------------------------#####
