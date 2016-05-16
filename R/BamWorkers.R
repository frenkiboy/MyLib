### INFO: Functions for working with bam files in R
### DATE: 27.03.2010.
### AUTH: Vedran Franke

# {1} LIBRARIES
# library(Rsamtools)
#/{1} LIBRARIES



# {2} CODE

	# ------------------------------------------------------------------------------- #
	# {{1}}
	### reads in a bam file - it needs a bam index in the same folder as the bam file
	BamReader = function(bam.path=NULL, what=c('flag', 'rname', 'strand', 'pos','seq', 'qual', 'mapq'), which.ranges=NULL, ind.path=bam.path, filter=TRUE, mapq.filter=TRUE, tag=character()){
		
		cat('Constructing the parameters...\n')
		if(is.null(which.ranges)){
			param = ScanBamParam(what=what)
		}else{
			param = ScanBamParam(what=what, which=which.ranges, tag=tag)
		}

		cat('Reading in the bam file...\n')
		bam = scanBam(file=bam.path, param=param, index=ind.path)[[1]]
		### checks if something was read in
		if(sum(unlist(lapply(bam, length))) == 0){
			cat('None of the reads were read!\n')
			return(NULL)
		}
		
		### removes the bad quality reads
		if(!is.null(bam$qual) & filter == TRUE){
			bam = ReadQualFilter(bam=bam)
		}else if(is.null(bam$qual) & filter == TRUE){
			stop('Qualities need to be read in if you want to filter the data.')
		}
		if(!is.null(bam$mapq) & mapq.filter == TRUE){
			bam = MapQualFilter(bam=bam)
		}else if(is.null(bam$mapq) & mapq.filter == TRUE){
			stop('Mapping quality does not exist!')
		}
		
		
		### removes nonmapped reads
		cat('Removing non mapped reads...\n')
		ind = !is.na(bam$pos)
		bam = lapply(bam, '[', ind)
		
		cat('Returning the data...\n\n')
		return(bam)
	
	}
	#/{{1}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{1.25}}
	# reads a one chromosome from the bam file
	BamReaderChr = function(bam.path, chr, chrs, filter=F, mapq.filter=T, values=T){
	
		if(is.null(bam.path) | is.null(chr) | is.null(chrs))
			stop('Not all parameters are assigned')
		chr.len=chrs$chr.len[chrs$chr == chr]
		names(chr.len) = chr
		which.ranges = GRanges(chr, IRanges(1, chr.len))
		bam = BamReader(bam.path, which.ranges = which.ranges, filter=filter, mapq.filter=mapq.filter)
		ranges = BamToGRanges(bam, chr.len, values=values)
		return(ranges)	
	}
	#/{{1}.25}
	
	
	
	# ------------------------------------------------------------------------------- #
	# {{1.5}}
	### reads in a bam file - it needs a bam index in the same folder as the bam file
	GappedBamReader = function(bam.path=NULL, ind.path=bam.path){
		
		cat('Constructing the parameters...\n')
		param = ScanBamParam(what=c("qual"), flag=scanBamFlag(isDuplicate=NA, isUnmappedQuery=FALSE), tag="NM")
		
		cat('Reading in the gapped bam file...\n')
		gbam = readBamGappedAlignments(bam.path, index=ind.path, param=param)
		
		### checks if something was read inq
		if(length(gbam) == 0){
			return(NULL)
		}
		
		good.quals=ReadQualFinder(values(gbam)$qual)
		gbam = gbam[good.quals]
		
		cat('Returning the data...\n\n')
		return(gbam)
	}
	#/{{1.5}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{2}}
	### filters the bam reads by quality
	### filters all reads that have 3 bases under 20
	ReadQualFilter = function(bam, qual=20){
		
		cat('Filtering the bam file...\n')
		good.qual.reads = ReadQualFinder(bam$qual, qual=qual)
		bam.subset = lapply(bam, function(x)x[good.qual.reads])
		return(bam)
	}
	#/{{2}}
	
	# ------------------------------------------------------------------------------- #
	# {{2.5}}
	### filters the bam reads by quality
	### filters all reads that have 3 bases under 20
	ReadQualFinder = function(read.quals, qual=20){
		
		library(ShortRead)
		cat('Calculating the quality scores...\n')
		q.mat = as(FastqQuality(SolexaQuality(read.quals)),'matrix')
		good.qual.reads = rowSums(q.mat < 20) >= 3
		return(good.qual.reads)
	}
	#/{{2.5}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{2.75}}
	### filters reads by mapping quality
	### takes all the reads with mapping quality greater than qual
	MapQualFilter = function(bam, qual=0){
		
		cat('Filtering by mapping quality...\n')
		mapq = bam$mapq
		mapq[is.na(mapq)] = 0
		ind = mapq > qual
		bam = lapply(bam, '[', ind)
		return(bam)
	}
	
	
	#/{{2.75}}
	
	# ------------------------------------------------------------------------------- #
	# {{3}}
	### removes all ranges that are nearer the edge than a given distance
	NearEdgeReadRemover = function(ranges, distance=200){
	
		cat('Removing reads close to the chromosome edges...\n')
		ranges = ranges[as.vector(!(strand(ranges) == '-' & end(ranges) < distance))]
		ranges = ranges[as.vector(!(strand(ranges) == '+' & seqlengths(ranges)[1] - start(ranges) < distance))]
		return(ranges)
	}
	#/{{3}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{4}}
	### finds the chromosome names from the bam file - uses samtools
# 	chrFinder = function(bam.path, filter=FALSE, output='data.frame'){
# 	
# 		chr.stat = strsplit(system(paste('samtools idxstats ', bam.path, sep=''), intern=T), split='\\t')
# 		chr = data.frame(do.call(rbind, lapply(chr.stat, '[', c(1,2))), stringsAsFactors=F)
# 		names(chr)=c('chr','chr.len')
# 		chr = chr[nchar(chr$chr) > 1,]
# 		chr$chr.len = as.numeric(as.character(chr$chr.len))
# 		if(filter == TRUE){
# 			cat('Filtering random chromosomes...\n')
# 			require(stringr)
# 			chr = chr[!str_detect(chr$chr, 'random'),]
# 		}
# 		if(output == 'vector'){
# 			cat('Outputting a named vector...\n')
# 			v = chr$chr.len
# 			names(v) = chr$chr
# 			return(v)
# 		}else{
# 			cat('Outputting a data.frame...\n')
# 			return(chr)
# 		}
# 	}
    
    chrFinder = function(bam.path, filter=FALSE, output='data.frame'){
    
        s = scanBamHeader(bam.path)
        st = s[[1]]$text
        st = do.call(rbind,st[names(st) == '@SQ'])
        st[,1] = str_replace(st[,1],'SN:','')
        st[,2] = str_replace(st[,2],'LN:','')
        if(filter == TRUE){
            st = st[!str_detect(st[,1], 'random')]
        }
        
        if(output == 'data.frame'){
            vst = data.frame(chr=st[,1], chrlen=as.numeric(st[,2]))    
        
        }else{
            vst = as.numeric(st[,2])
            names(vst) = st[,1]
        }
        return(vst)
    }
    
	#/{{4}}
	

	# ------------------------------------------------------------------------------- #
	# {{4}}
	### converts a bam file to a Granges object
	BamToGRanges = function(bam, chr.len=NULL, values=F, chr.remove=T){
		
		if(is.null(bam))
			return(NULL)
			
		if(is.null(chr.len))
			stop('Chromosome lengths are missing')
			
		### constructs a named vector for 
		if(is.data.frame(chr.len)){
			v = chr.len[,2]
			names(v) = chr.len[,1]
			chr.len = v
		}
		
		if(chr.remove == T){
			cat('Filtering the chromosomes...\n')
			bam = ChrFilter(bam, names(chr.len))
		}
		
		### checks whether the chromosome names are equal
		cat('Converting bam to GRanges...\n')
		if(!all(bam$rname %in% names(chr.len))){
			stop('Bam and chrlen do not have matching chromosome names!')
		}
		cat('Finding off chromosome reads...\n')
		chrlens = chr.len[match(bam$rname, names(chr.len))]
		ind = rep(FALSE, length(bam$rname))
		ind[bam$pos + width(bam$seq) < chrlens] = TRUE
		cat('Selecting  only on chromosome reads...\n')
		bam = lapply(bam, '[', ind)
		
		if(length(unique(unlist(lapply(bam, length)))) != 1)
			stop('All bam elements do not have the same length!')
		
		ranges = GRanges(seqnames=as.character(bam$rname), ranges=IRanges(start=bam$pos, width=width(bam$seq)), strand=bam$strand, seqlengths=chr.len)
		
		### fills up the element metadata
		if(values == T){
			### finds the metadata column names
			cols = c('rname','pos','strand')
			vals = setdiff(names(bam), cols)
			for(i in vals)
				values(ranges)[[i]] = bam[[i]]	
		}
			
		return(ranges)
	}
	#/{{4}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{5}}
	### removes wierd chromosomes
	ChrFilter = function(bam, chr.names=NULL){
	
		if(is.data.frame(chr.names)){
			chr.names = chr.names[,1]
		}
		
		id = as.character(bam$rname) %in% chr.names
		bam = lapply(bam, '[', id)
		return(bam)
	
	}
	#/{{5}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{6}}
	BamToBedGraph = function(bam.path, chrs=NULL, outfile=NULL, extend=0, colnum, unique=F){
		
		source('/home/members/vfranke/MyLib/ChipSeq.Functions.R')
		if(is.null(chrs))
			chrs = chrFinder(bam.path)
			
		if(file.exists(outfile))
			file.remove(outfile)
		
		l.cov = list()
		for(chr in sort(as.character(chrs[,1]))){
			print(chr)
			chr.len = chrs[chrs[,1] == chr, 2]
			names(chr.len) = chr
			which.ranges = GRanges(chr, IRanges(1, chr.len))
			bam = BamReader(bam.path, which.ranges = which.ranges)
			ranges = BamToGRanges(bam, chr.len=chr.len)
			
			if(unique==T){
				message('Removing duplicated reads')
				ids = paste(as.vector(seqnames(ranges)), start(ranges), as.vector(strand(ranges)))
				ranges = ranges[!duplicated(ids)]
			}
			seqlengths(ranges) = chr.len
			ranges = NearEdgeReadRemover(ranges, extend)
			l.cov[[chr]] = coverage(resize(ranges, width=extend))[[1]]
		}
		ExportCoverage(l.cov, outfile, colnum=colnum)
	}
	#/{{6}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{7}}
	### selects a reads overlapping a given bed file
	BamSelectBed = function(bam.file, bed.file, outpath=dirname(bam.file), bam.name=sub('.bam','.ss.bam',basename(bam.file))){
	
		bam.outfile=file.path(outpath, bam.name)
		command = paste('samtools view -b -L', bed.file, bam.file,'>',bam.outfile, sep=' ')
		message(command)
		system(command)
		return(bam.outfile)
	
	}
	#/{{7}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{8}}
	### calculates read statistics for a bam file
	Percentager = function(d){
	
		d = d/d[,1]
		names(d) = paste(names(d), '%', sep='.')
		return(d[,-1])
	}
	
	# ------------------------------------------------------------------------------- #
	chrStats = function(bam.file){
	
		l = list()
		chr.len = chrFinder(bam.file)
		for(chr in chr.len$chr){
			
			print(chr)
			which.ranges = GRanges(seqnames=chr, IRanges(1, chr.len$chr.len[chr.len$chr == chr]))
			bam = BamReader(bam.file, which.ranges=which.ranges)
			names(chr.len) = chr
			ranges = BamToGRanges(bam, chr.len=chr.len)
			ranges.uniq = ranges[!duplicated(paste(seqnames(ranges), start(ranges), strand(ranges)))]
			l[[chr]] = data.frame(mapped = length(ranges), 
								  mapped.m = length(ranges[strand(ranges)=='-']), 
								  mapped.p = length(ranges[strand(ranges)=='+']), 
								  mapped.uniq = length(ranges.uniq),
								  mapped.uniq.m = length(ranges.uniq[strand(ranges.uniq)=='-']),
								  mapped.uniq.p = length(ranges.uniq[strand(ranges.uniq)=='+']))
		}
		d = do.call(rbind, l)
		d = cbind(d, Percentager(d[,1:3]), 'mapped.uniq.%'=d$mapped.uniq/d$mapped, Percentager(d[,4:6]))
		return(d)
	}
	
	
	#/{{8}}
	
	
	# ------------------------------------------------------------------------------- #
	#/{{9}}
	BamIndex = function(bam.path){
		bam.index.path = paste(bam.path, 'bai', sep='.')
		command = paste('samtools index', bam.path, bam.index.path)
		message('Running: ', command)
		system(command)
		return(bam.index.path)
	}
	#/{{9}}
	
	
	# ------------------------------------------------------------------------------- #
	#/{{10}}
	BamSort = function(bam.path){
		command = paste('samtools sort', bam.path, sub('.bam','',bam.path))
		message('Running: ', command)
		system(command)
	}
	
	#/{{10}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{11}}
	### takes a bam path file and gets the coverage by window
	BamWindowerChr = function(bam.file, chr, chr.len, window = 200, resize=1){
	
		which.ranges = GRanges(chr, IRanges(1, chr.len))
		bam = BamReader(bam.file, which.ranges=which.ranges)
		ranges = BamToGRanges(bam)
		seqnames(ranges) = chr
		seqlengths(ranges) = chr.len
		
		message('Calculating coverage...')
		cov = coverage(resize(ranges, width=1))[[1]]
		
		message('Getting the views...')
		v = Views(cov, start = seq(1, chr.len, window), width=window)
		d = data.frame(start(v), viewSums(v))
		return(d)
	}
	#/{{11}}
	
	# ------------------------------------------------------------------------------- #
	# {{12}}
	### takes a bam file and returns the name
	BamName = function(bam.file){
		sub('.bam','',basename(bam.file))
	}
	#/{{12}}
	
	# ------------------------------------------------------------------------------- #
	# {{13}}
	### counts the overlaps between two sets of regions with an included weight
	GetOverlaps = function(reg1,reg2, colname=NULL){
	
		library(data.table)
		if(is.null(colname))
			stop('Colname needs to be defined')
		fo = data.table(as.matrix(findOverlaps(reg1,reg2)), ignore.strand=T)
		fo$weight = values(reg2)[[colname]][fo$subjectHits]
		fo = fo[,sum(weight), by=queryHits]
		v = rep(0, length(reg1))
		v[fo$queryHits] = fo$V1
		return(v)
	}
	# {{13}}
    
    # ------------------------------------------------------------------------------- #
	# {{13}}
	### counts the overlaps between two sets of regions
	CountOverlaps = function(reg1,reg2, colname=NULL, setname=NULL, ignore.strand=FALSE){
	
		library(data.table)
		library(GenomicRanges)
		fo = data.table(as.matrix(findOverlaps(reg1,reg2, ignore.strand=ignore.strand)))
		if(is.null(colname)){
			fo$name = names(reg1)[fo$queryHits]
		}else{
			fo$name = values(reg1)[[colname]][fo$queryHits]
		}
		fo = fo[,length(subjectHits),by=name]
		
		if(!is.null(setname))
			setnames(fo,2,setname)
		return(fo)
	}
	# {{13}}
	
	
	# ------------------------------------------------------------------------------- #
	# {{4}}
	list.bamfiles=function(x){
	    return(list.files(x, full.names=TRUE, pattern='bam$', recursive=TRUE))
	}
	
#/{2} CODE

