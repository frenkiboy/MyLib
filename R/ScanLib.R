
#####-------------------- INTERNALS ---------------------------#####
### finds the overlaps between two bed formatted files
.OverlappingIntervals = function(bed1, bed2, minoverlap=1){

	require(IRanges)
	### converts bed do IRanges
	bed1.range = IRanges(bed1[,2], bed1[,3])
	bed2.range = IRanges(bed2[,2], bed2[,3])
	
	### finds overlaps
	ovlap = as.matrix(findOverlaps(bed1.range, bed2.range, minoverlap = minoverlap))
	
	###creates the output table and returns the value
	bed.ovlap = cbind(bed1[ovlap[,1],], bed2[ovlap[,2],])
	return(bed.ovlap)
}

.FindNeighbour = function(bed1, bed2){
	
	require(IRanges)
	### converts bed do IRanges
	bed1.range = IRanges(bed1[,2], bed1[,3])
	bed2.range = IRanges(bed2[,2], bed2[,3])
	
	### finds the nearest neighbour of bed1 in bed2
	near = nearest(bed1.range, bed2.range)
	
	bed = cbind(bed1, bed2[near,])
	return(bed)
}
#####-------------------/ INTERNALS /--------------------------#####



#####-------------------- FUNCTIONS ---------------------------#####
### garbage collect function
collect.garbage = function()
{
  while(gc()[2,4] != gc()[2,4]){}
}

### makes a unique key for each row in a bed file
BedKey = function(bed, sep="_"){
	bed$key = paste(bed[,1], bed[,2], bed[,3], sep=sep)
	return(bed)
}

### takes two bed files and returns the overlap statistics
BedSet = function(bed1, bed2, r=4){

	stopifnot(!is.null(bed1) | !is.null(bed2))
	bed1$uniq.id = 1:nrow(bed1)
	bed2$uniq.id = 1:nrow(bed2)
		
	o = FeatureOverlapBED(bed1, bed2)
	
	l=list()
	d = data.frame(bed1	 = nrow(bed1),
						 bed1.o  = sum(bed1$uniq.id %in% o$uniq.id.a),
						 bed1.pc = sum(bed1$uniq.id %in% o$uniq.id.a) / nrow(bed1), 
						 bed2	 = nrow(bed2),
						 bed2.o  = sum(bed2$uniq.id %in% o$uniq.id.b),
						 bed2.pc = sum(bed2$uniq.id %in% o$uniq.id.b) / nrow(bed2)
						 )
						 
	return(round(d, r))
}


### annotates bed1 by overlaps with bed2
DesignateSetBED = function(bed1, bed2, set=deparse(substitute(bed2)), index=c('No','Yes')){
	
		stopifnot(!is.null(bed1) | !is.null(bed2))
		bed1$uniq.id = 1:nrow(bed1)
		
		o = FeatureOverlapBED(bed1, bed2)
		bed1[[set]] = index[1]
		bed1[[set]][bed1$uniq.id %in% o$uniq.id.a] = index[2]
		return(subset(bed1, select=-uniq.id))
}


### takes a dataframe of ids, and makes one vector of sorted ids
IDsorter = function(ids, abbrev = ""){
		tmp.id = as.character(unlist(ids))
		tmp.id = unique(unlist(strsplit(tmp.id, split = ":")))
		if(abbrev != ""){
			tmp.id = paste(abbrev, ".", tmp.id, sep = "")
		}
		new.id = paste(sort(tmp.id), collapse = ":")
		return(new.id)
}

### Takes a list of id's and designates them with an abbreviation
 ## calls IDsorter
IDchanger = function(data, abbrev = ""){
	id.col = grep("id", colnames(data))
	
	if(id.col == 0){
		cat("No colunms designated as id!!!\n")
		return(NULL)
	}
	ids = as.vector(data[,id.col])
	data[,id.col] = sapply(ids, IDsorter, abbrev)
	return(data)
}

### Checks if the bed file has an column name containing 'id', if it doesn't then it creates 1:nrow(bed) as ids
IDChecker = function(bed){

	if(sum(grep('id', names(bed))) ==0 ){
		bed$id = as.character(1:nrow(bed))
	}else{
		names(bed)[grep('id', names(bed))] = 'id'
	}
	return(bed)

}

### makes the crossection of two sets
FeatureCrossSection = function(bed1, bed2){
	
	require(IRanges)
	### removes negative length ranges
	bed1 = bed1[bed1[,2] < bed1[,3],]
	bed2 = bed2[bed2[,2] < bed2[,3],]
	
	### checks if the bed files have an id column, if they don't creates 1:nrow(bed), as id
	bed1 = IDChecker(bed1)
	bed2 = IDChecker(bed2)
	
	### reduces both bed files to nonoverlaping region sets
	bed1 = FeatureUnion(bed1)
	bed2 = FeatureUnion(bed2)
	
	### data.frame for holding the results
	l = list()
	
	chrs = unique(as.vector(bed1[,1]))
	chrs = chrs[chrs %in% unique(as.vector(bed2[,1]))]
	for(i in seq(along = chrs)){
		a = bed1[bed1[,1] == chrs[i],]
		b = bed2[bed2[,1] == chrs[i],]
		a.ids = as.vector(a[,grep('id', names(a))])
		b.ids = as.vector(b[,grep('id', names(b))])
		
		### finds overlapping intervals
		a.range = IRanges(a[,2], a[,3])
		b.range = IRanges(b[,2], b[,3])
		p = intersect(a.range, b.range)
		
		### concatenates the ids
		a.ovlap = as.matrix(findOverlaps(a.range, p, multiple = T))
		b.ovlap = as.matrix(findOverlaps(b.range, p, multiple = T))
		id.ind = merge(a.ovlap, b.ovlap, by.x = 'subject', by.y = 'subject')
		ids = cbind(a.ids[id.ind[,2]], b.ids[id.ind[,3]])
		### id takes uniq ids for each intersect and concatenates them with : as separator - don't ask me how it works
		if(nrow(ids) != 0){
			ids.con = tapply(1:nrow(ids), id.ind[,1], FUN = function(x){id=IDsorter(ids[as.vector(x),1:2]); return(id)})
		
			d = data.frame(chrs[i], start(p), end(p), ids.con)
			names(d) = c('chr','start','end','id')
			l[[i]] = d
		}
	}
	data.intersect = do.call(rbind,l)
	return(data.intersect)
}

### takes two bed files and returns overlapped rows between them
FeatureOverlapBED<-function(bed1,bed2, minoverlap=1)
{
    #get unique shared chromosomes between 2 files
    chrs=as.vector(unique(bed1[,1]))
    chrs=chrs[chrs %in% unique(bed2[,1])]
	l = list()
	### designates bed1 colnames as x.a, bed 2 colnames as x.b
	names(bed1) = paste(names(bed1), '.a', sep = "")
	names(bed2) = paste(names(bed2), '.b', sep = "")
    for (i in 1:length(chrs))
    {
		cat('Overlapping chr:',chrs[i],"\n")
        overlap=.OverlappingIntervals(
        bed1[bed1[,1]==as.character(chrs[i]),],
        bed2[bed2[,1]==as.character(chrs[i]),],
		minoverlap=minoverlap)
        
        l[[i]]= overlap
    }
	cat('Overlapping done!!!\n')
	cat('Creating the results\n\n')
	result = do.call(rbind, l)
    if(nrow(result)==0){return(NULL)}
    return(result)
}

### takes two bed files and outputs the overlaping rows between them, that are on the same strand
FeatureOverlapBED.strand<-function(bed1,bed2, minoverlap=1)
{
    ### get unique shared chromosomes between 2 files
    chrs=unique(as.vector(bed1[,1]))
    chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
    l = list()
	### renames the colnames of bed1 and bed2
	names(bed1) = paste(names(bed1), '.a', sep = "")
	names(bed2) = paste(names(bed2), '.b', sep = "")
	
	### finds the colnames named strand
	bed1.strand = grep('strand', names(bed1))
	bed2.strand = grep('strand', names(bed2))
    for (i in 1:length(chrs) )
    {
		cat('Overlapping chr:',chrs[i],"\n")
        overlap.plus=.OverlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,bed1.strand]=="+",],
        bed2[bed2[,1]==chrs[i] & bed2[,bed2.strand]=="+",],
		minoverlap=minoverlap)

        overlap.minus=.OverlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,bed1.strand]=="-",],
        bed2[bed2[,1]==chrs[i] & bed2[,bed2.strand]=="-",],
		minoverlap=minoverlap)
        
        l[[i]] = rbind(overlap.plus,overlap.minus)
    }
	cat('Overlapping done!!!\n')
	cat('Creating the results\n\n')
	result = do.call(rbind,l)
	cat('Returning the results\n\n')
    if(nrow(result)==0){return(NULL)}
    return(result)
}

### takes a vector of chromosome names and returns a logical vector of full chromosomes
FullChromosomes = function(chrs, y=TRUE, m=FALSE){

	if(!is.character(chrs))
		stop('chrs needs to be a character vector\n')
	require(stringr)
	ind = (!str_detect(chrs, 'random')) & (!str_detect(chrs, 'NT'))
	if(y==FALSE)
		ind = ind & !str_detect(chrs, 'chrY')
	if(m==FALSE)
		ind = ind & !str_detect(chrs, 'chrM')
	return(ind)
}

### finds the nearest neighbour of two bed strands
NearestNeighbourBED<-function(bed1,bed2)
{
	### get unique shared chromosomes between 2 files
	chrs=unique(as.vector(bed1[,1]))
	chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
	l = list()
	### renames the colnames of bed1 and bed2
	names(bed1) = paste(names(bed1), '.a', sep = "")
	names(bed2) = paste(names(bed2), '.b', sep = "")
	
	for (i in 1:length(chrs) )
	{
		cat('Finding neighbours chr:',chrs[i],"\n")
	    nn =.FindNeighbour(
		bed1[bed1[,1]==chrs[i],],
		bed2[bed2[,1]==chrs[i],]) 
	    l[[i]] = nn
	}
	cat('Finding neighbours done!!!\n')
	cat('Creating the results\n\n')
	result = do.call(rbind,l)
	cat('Returning the results\n\n')
	if(nrow(result)==0){return(NULL)}
	return(result)
}


### finds the nearest neighbour of two bed strands
NearestNeighbourBED.strand<-function(bed1,bed2)
{
    ### get unique shared chromosomes between 2 files
    chrs=unique(as.vector(bed1[,1]))
    chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
    l = list()
	### renames the colnames of bed1 and bed2
	names(bed1) = paste(names(bed1), '.a', sep = "")
	names(bed2) = paste(names(bed2), '.b', sep = "")
	
	### finds the colnames named strand
	bed1.strand = grep('strand', names(bed1))
	bed2.strand = grep('strand', names(bed2))
    for (i in 1:length(chrs) )
    {
		cat('Finding neighbours chr:',chrs[i],"\n")
        nn.plus=.FindNeighbour(
        bed1[bed1[,1]==chrs[i] & bed1[,bed1.strand]=="+",],
        bed2[bed2[,1]==chrs[i] & bed2[,bed2.strand]=="+",])

        nn.minus=.FindNeighbour(
        bed1[bed1[,1]==chrs[i] & bed1[,bed1.strand]=="-",],
        bed2[bed2[,1]==chrs[i] & bed2[,bed2.strand]=="-",],)
        
        l[[i]] = rbind(nn.plus, nn.minus)
    }
	cat('Finding neighbours done!!!\n')
	result = do.call(rbind,l)
	cat('Creating the results\n\n')
    if(nrow(result)==0){return(NULL)}
    return(result)
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



### makes the union of two sets
FeatureUnion = function(bed){

	require(IRanges)
	
	### prepares the data
	## removes start > end rows
	data = bed[bed[,2] < bed[,3],]
	### sets the id's
	data = IDChecker(data)
	data = data[order(data[,1], data[,2]),]
	
	###data.frame for the results
	l = list()
	chr = as.vector(unique(bed[,1]))
	cat('Collapsing the regions...\n')
	for(i in seq(along = chr)){
		
		cat('Collapsing chr:', chr[i],"\n")
		### reduces the regions to a nonoverlapping set
		a = data[data[,1] == chr[i],]
		a.range = IRanges(a[,2], a[,3])
		a.union = reduce(a.range)
		### find the overlaps between the reduced regions and the old set
		d = as.matrix(findOverlaps(a.range, a.union, multiple = T))
		### creates new ids
		ids = tapply(a$id, d[,2], IDsorter)
		
		if(length(a.union) != 0){
			u = data.frame(chr[i], start(a.union), end(a.union), ids, stringsAsFactors = F)
			names(u) = c('chr','start','end','id')
			l[[i]] = u
		}
	}
	reg.union = do.call(rbind, l)
	cat("intervals:", nrow(data), ",collapsed to:", nrow(reg.union), "\n")
	return(reg.union)
}

### takes the cross section of the desired regions with the filter
RegionFilter = function(regions,filter){
	
	### filters given intervals (takes cross section of regions & filter )
	reg.filter = FeatureCrossSection(regions, filter)
	cat('Regions filtered!\nRegions left: ', nrow(reg.filter), '\n')
	
	### checks if there was any overlap
	if(length(reg.filter) == 0 ){
		reg.filter = data.frame(matrix(nrow=0, ncol = length(c(names(regions), names(filter)))))
		names(reg.filter) = c(names(regions), names(filter))
		return(reg.filter)		
	}	
	return(reg.filter)
}


### intersects a bed file with a given rle file given a rle cutoff
## gives as output only intersect lengths that are greater then the 1 quartile
RleBedIntersect = function(bed, rle, cutoff = 0, first.quartille = F){

	require(IRanges)
	### checks for the id column
	bed = IDChecker(bed)
	### reduces the bed file to a non overlapping set
	bed = FeatureUnion(bed)

	
	### list for results
	l = list()
	
	### takes the chromosomes
	chr = as.vector(sort(unique(bed[,1])))
	chr = chr[chr %in% names(rle)]
	### loops the each chromosome and does the intersect
	cat('Looping through the chrs and intersecting..\n')
	for(i in seq(along = chr)){
		
		cat('Intersecting chr:',chr[i], "\n")
		bed.chr = bed[bed[,1] == chr[i],]
		### collapses the bed file to a nonoverlapping regions set 
		rle.chr = rle[[chr[i]]]
		
		### converts the stuff to IRanges objects
		iranges.bed = IRanges(bed.chr[,2], bed.chr[,3])
		iranges.rle = IRanges(rle.chr > cutoff)
		
		### takes the intersect
		int.ranges = intersect(iranges.rle, iranges.bed)
		### finds the id's from the bed for each overlapping patch
		int.pairs = as.matrix(findOverlaps(int.ranges, iranges.bed, multiple = T))
		bed.id = tapply(int.pairs[,2], int.pairs[,1], FUN = function(x){id=IDsorter(as.vector(bed.chr$id[x])); return(id)})
		
		### checks if there are any overlaps
		if(length(int.ranges) != 0){
			d = data.frame(chr[i], start(int.ranges), end(int.ranges), bed.id, stringsAsFactors = F)
			names(d) = c('chr','start','end','id')
			l[[i]] = d
		}
		
		
	}
	cat('Intersecting done\n')
	reg.intersect = do.call(rbind, l)
	### removes <= 0 length coordinates
	reg.intersect = reg.intersect[(reg.intersect[,3] > reg.intersect[,2]),]
	### removes all lengths < 1 quartille
	if (first.quartille == T){
			cat('Removing all regions whose lengths < 1\'st quartille\n')
			reg.intersect = reg.intersect[(reg.intersect[,3] - reg.intersect[,2]) > summary((reg.intersect[,3] - reg.intersect[,2]))[2],]
	}
	cat('Returning the intersected regions\n')
	return(reg.intersect)
}


### reads in and parses JASPAR formated PFM file; input - fasta like matrix file; output list
PFMparser = function(path){
	
	###opens the connection and reads in the file
	cat('Reading in the PFM matrix\n')
	connection = file(path)
	s = scan(connection,  what = 'character', sep = '>')
	close(connection)
	
	### removes brackets
	s = sub("\\[", "", s)
	s = sub("\\]", "", s)
	s = strsplit(s, split = "\\s+")
	a = seq(2, length(s), 6)
	l = vector(mode = 'list', length = length(a))
	

	### loops through each matrix and parses it into a list
	cat('Parsing the matrices...\n')
	for(i in seq(along=a)){
		
		m = data.frame(do.call(rbind, s[1:4+a[i]]), stringsAsFactors  = F)
		rownames(m) = m[,1]
		m = m[,-1]
		m = data.matrix(m)
		l[[i]]$ID = unlist(s[[a[i]]])
		l[[i]]$PFM = m
		
	}
	cat('Parsing done, returning the values.\n')
	return(l)
}

### writes a PWM/PFM in a jaspar formatted style - can be parsed by PFM parser
writePWM = function(pwm, file='./'){

	str = apply(pwm, 1, paste, collapse=' ')
	str = paste('[', str, ']')
	names(str) = rownames(pwm)
	### prints the header for the matrix
	cat('>', basename(file), '\n', file=file, append=T, sep='')
	### prints the rest of the matrix
	### when you get time write this properly
	for(i in 1:length(str)){
		cat(names(str)[i], str[i], '\n', file=file, append=T, sep='')
	}	
}


### takes a list of PFMs and converts them to PWMs
PFMtoPWM = function(PFM, p = rep(0.25, 4)){
	
	for (i in 1:length(PFM)){
		
		pfm = PFM[[i]]$PFM
		PWM = matrix(ncol = ncol(pfm), nrow = nrow(pfm))
		rownames(PWM) = rownames(pfm)
		N = colSums(pfm)
		### goes through each frequency element and calculates the probability
		for(j in 1:nrow(pfm)){
		
			curr.p = p[names(p) == rownames(pfm)[j]]
			for(k in 1:ncol(pfm)){
				### Wasserman & Sandelin NRG 2004
				PWM[j,k] = log2((pfm[j,k] + sqrt(N[k]) * curr.p) / (N[k] + sqrt(N[k])) / curr.p)			
			}
		}
		PFM[[i]]$PWM = PWM
	}
	return(PFM)
}


### checks if the file has a header
HeaderChecker = function(path){
	file = read.table(path, header = F, nrows = 1, stringsAsFactors = F)
	if(sum(sapply(file[1,], is.character)) == ncol(file)){
		header = T
	}else{
		header = F
	}
	return(header)
	
}


### takes a bed formatted data.frame and scans the sequences - use for more than 1500 seq.
MatchPWM.Large = function(pwm, data, min.score){

	chr = unique(data[,1])
	match = list()
	for(i in seq(along=chr)){
		data.tmp = data[data[,1] == chr[i],]
		v = Views(unmasked(Mmusculus[[chr[i]]]), start=data.tmp[,2], end=data.tmp[,3])
		match[[i]] = countPWM(pwm=pwm, v, min.score=min.score)
	}
	return(match)
}

# ---------------------------------------------- #
# pastes the date to the name of the outdir
OutpathDate = function(name){
	paste(strsplit(as.character(Sys.time()), '\\s')[[1]][1], name, sep='_')
}

# creates the outpath directory with a given name and date
CreateOutputDir = function(outpath, name){
	outpath = file.path(outpath, OutpathDate(name))
		dir.create(outpath, showWarnings=F)
	return(outpath)
}


# ---------------------------------------------- #
### creates the output name
OutNameMaker = function(path, name){
	file = unlist(strsplit(path, split = "/"))
	file = paste(name, file[length(file)], sep = ":")
	file = sub("....$", "", file)
	return(file)
}

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


### extends the coordinates of a bed file
RegionExtender = function(bed, downstream = 1000, upstream = 1000, strand = F){
	
	if(strand == T){
		
		### checks if the strand column exists
		if (length(grep("strand", names(bed))) == 0){
			cat('Please specify the strand column!\n')
			return(NULL)
		}
		
		bed.m = bed$strand == '-'
		bed[bed.m,2] = bed[bed.m,2] - downstream
		bed[bed.m,3] = bed[bed.m,3] + upstream
		
		bed.p = bed$strand == '+'
		bed[bed.p,2] = bed[bed.p,2] - upstream
		bed[bed.p,3] = bed[bed.p,3] + downstream
		colnames(bed)[2:3] = c('start','end')
		return(bed)
	}
	else if(strand == F){
		bed[,2] = bed[,2] - upstream
		bed[,3] = bed[,3] + downstream
		colnames(bed)[2:3] = c('start','end')
		return(bed)
	}
}


### takes a bed file and designates start and end as same coordinates, based on the strand
StartStrandMaker = function(bed){
	
	### for + strand, end = start
	strand = grep('strand', names(bed))
	bed[bed[, strand] == '+', 3] = bed[bed[,strand] == '+', 2]
	### for - strand, start = end
	bed[bed[,strand] == '-', 2] = bed[bed[,strand] == '-', 3]
	return(bed)
}

### takes a bed file and dsignates the promoter region, but truncates the downstream region if the transcipts is shorter than the designated region
TruncPromDesignator = function(bed, downstream=500, upstream=500){

	library(stringr)
	w.ind = (bed$end - bed$start) > downstream
	s.ind = bed[str_detect(names(bed), 'strand')] == '+'
	prom = StartStrandMaker(bed)
	prom$start[s.ind] = prom$start[s.ind] - upstream
	prom$end[!s.ind]  = prom$end[!s.ind]  + upstream
	
	prom$end[s.ind & w.ind] = prom$end[s.ind & w.ind] + downstream
	prom$end[s.ind & !w.ind] = bed$end[s.ind & !w.ind]
	
	prom$start[!s.ind & w.ind] = prom$start[!s.ind & w.ind] - downstream
	prom$start[!s.ind & !w.ind] = bed$start[!s.ind & !w.ind]
	return(prom)
}

### Takes two bed files and returns x1 that doesnt overlap x2
TakeUniqRegions = function(x1, x2){

	require(GenomicRanges)
	gp = GRanges(x1[,1], IRanges(x1[,2], x1[,3]))
	gs = GRanges(x2[,1], IRanges(x2[,2], x2[,3]))
	o = countOverlaps(gp, gs)
	x3 = x1[o == 1,]
	x3
}



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


### selects unique elements from a granges object based on start strand and id
UniqueRanges = function(ranges){

	cat('Taking unique ranges...\n')
	id = paste(as.character(seqnames(ranges)), start(ranges), as.character(strand(ranges)))
	ranges.uniq = ranges[!duplicated(id)]
	return(ranges)
}

### takes a bed file and designates start and end as as the end coordinates
EndStrandMaker = function(bed){

	### for + strand, end = start
	strand = grep('strand', names(bed))
	bed[bed[, strand] == '+', 2] = bed[bed[,strand] == '+', 3]
	### for - strand, start = end
	bed[bed[,strand] == '-', 3] = bed[bed[,strand] == '-', 2]
	return(bed)
}


###calculates TFBS z score
Zscore = function(filtered.hit, bck.hit, filtered.size, bck.size){
	exp.hit = bck.hit * (filtered.size/bck.size)
	p = bck.hit/bck.size
	sd.hit = sqrt(p*(1-p)*filtered.size)
	
	if(sd.hit == 0){
		return(0)
	}else{
		z.score = (filtered.hit - exp.hit - 0.5) / sd.hit
		return(z.score)
	}
}

### converts Xstring set to data frame
XStringToDataFrame = function(x){
	names(x) = paste(seq(along=x), '.', sep='')
	start = unlist(lapply(x, start))
	end = unlist(lapply(x, end))
	ids = sapply(names(start), function(a)unlist(strsplit(a, split='\\.'))[1])
	data = data.frame(start=start, end=end, id=ids, stringsAsFactors=F)
	return(data)
}

### prints the index of the loop
WorkPrinter = function(x){
	cat('Working on:', x, '\n')
}


# ---------------------------------------------------------------------------------------- #
### annotates the ranges with the corresponding list
setGeneric("AnnotateRanges",
           function(region, annotation,
                    ignore.strand=FALSE,
                    type='precedence',
                    null.fact='None',
                    collapse.char=':',
                    precedence=NULL,
					id.col=NULL)
               standardGeneric("AnnotateRanges") )

setMethod("AnnotateRanges",signature("GRanges","GRangesList"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':'){
    
    if(! class(region) == 'GRanges')
        stop('Ranges to be annotated need to be GRanges')
    
    if(! all(sapply(annotation, class) == 'GRanges'))
        stop('Annotating ranges need to be GRanges')
    
    if(!type %in% c('precedence','all'))
        stop('type may only be precedence and all')
    
    require(data.table)
    require(GenomicRanges)
    cat('Overlapping...\n')
    if(any(names(is.null(annotation))))
        stop('All annotations need to have names')
        
    if(class(annotation) != 'GRangesList')
        annotation = GRangesList(lapply(annotation, function(x){values(x)=NULL;x}))
    
    a = suppressWarnings(data.table(as.matrix(findOverlaps(region, annotation, ignore.strand=ignore.strand))))
    a$id = names(annotation)[a$subjectHits]
    a$precedence = match(a$id,names(annotation))
    a = a[order(a$precedence)]
    
    if(type == 'precedence'){
        cat('precedence...\n')
        a = a[!duplicated(a$queryHits)]
        annot = rep(null.fact, length(region))
        annot[a$queryHits] = a$id
    }
    if(type == 'all'){
        cat('all...\n')
        a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
        annot = rep(null.fact, length(region))
        annot[a$queryHits] = a$id
        
    }
    return(annot)
    
})

setMethod("AnnotateRanges",signature("GRanges","list"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':'){
		  
		  if(!all(unlist(lapply(annotation, 'class')) == 'GRanges'))
			stop('all elements of annotation need to be GRanges objects')
			
			annotation = GRangesList(lapply(annotation, function(x){values(x)=NULL;x}))
			AnnotateRanges(region, annotation, ignore.strand, type, null.fact, collapse.char)
		  
		  })
		  
		  

setMethod("AnnotateRanges",signature("GRanges","GRanges"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':', precedence=NULL, id.col=NULL){
                    
				if(is.null(id.col))
					stop('Annotation needs to have a specified id')
			
				if(is.null(values(annotation)[[id.col]]))
					stop('id.col is not a valid column')
				 
            
				if(!is.null(precedence)){
                if(!all(precedence %in% annotation$id))
                    stop('all precednce leveles have to be in the annotation id')
				}else{
					type='all'
					message('type set to all when precedence is not defined')
				}
                 
				if(!type %in% c('precedence','all'))
					stop('type may only be precedence and all')
              
				require(data.table)
				require(GenomicRanges)
				cat('Overlapping...\n')
				if(any(names(is.null(annotation))))
					stop('All annotations need to have names')
              
				a = suppressWarnings(data.table(as.matrix(findOverlaps(region, annotation, ignore.strand=ignore.strand))))
				a$id = values(annotation)[[id.col]][a$subjectHits]
              
              
              if(type == 'precedence'){
                  cat('precedence...\n')
                  a$precedence = match(a$id,precedence)[a$subjectHits]
                  a = a[order(a$precedence)]
                  a = a[!duplicated(a$queryHits)]
                  
              }
              if(type == 'all'){
                  cat('all...\n')
                  a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
              }
			  annot = rep(null.fact, length(region))
              annot[a$queryHits] = a$id
              return(annot)
              
          })

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


### given two ranges - gets the overlapping data table 
GetOverlaps = function(reg1,reg2, colname=NULL, ignore.strand=FALSE, funct='sum'){
	
    fun = match.fun(funct)
	if(is.null(colname)){
		cat('Using the default colname...\n')
		values(reg2)$weight=1
		colname=weight
	}
	require(data.table)
	fo = data.table(as.matrix(findOverlaps(reg1,reg2, ignore.strand=ignore.strand)))
	fo$weight = values(reg2)[[colname]][fo$subjectHits]
	fo = fo[,fun(weight), by=queryHits]
	v = rep(0, length(reg1))
	v[fo$queryHits] = fo$V1
	return(v)
}
    
    # ---------------------------------------------------------------------------------------- #
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

    # ---------------------------------------------------------------------------------------- #
    dtfindOverlaps = function(reg1, reg2, ignore.strand=FALSE){
        require(data.table)
        require(GenomicRanges)
        data.table(as.matrix(findOverlaps(reg1,reg2, ignore.strand=ignore.strand)))
    }
    
    GCoords = function(x, cor.sep='-', chr.sep=':'){
        paste(as.character(seqnames(x)), paste(start(x), end(x), sep=cor.sep), sep=chr.sep)
    }



	# ---------------------------------------------------------------------------------------- #
	WriteRegs = function(g, outpath, name=NULL, strand=TRUE, sname=c('-'='m','+'='p','*'='n')){
	
		if(is.null(name))
			stop('Please specify the name')
		
	
		
		d = as.data.frame(g)
		if(!strand)
			d$strand = '*'
		
		for(i in unique(d$strand)){
			s = sname[i]
			print(s)
			write.table(d[d$strand == i,1:3], file.path(outpath, DateNamer(paste(name, s, 'bed', sep='.'))), row.names=F, col.names=F, quote=F, sep='\t')
		}
	}
	
	# ---------------------------------------------------------------------------------------- #
	flipStrand = function(g){
	    levels(strand(g)) = c('-','+','*')
	    return(g)
	}
	
	
	# ---------------------------------------------------------------------------------------- #
	chrSeqlevels = function(g){
	    seqlevels(g) = paste0('chr',seqlevels(g))
	    return(g)
	}
	
	

#####--------------------/FUNCTIONS/---------------------------#####

