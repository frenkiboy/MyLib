### INFO: Converts between genomic formats
### DATE: 31.03.2011
### AUTHOR: v+


# {1} LIBRARIES
#/{1} LIBRARIES


# {2} CODE
	# ------------------------------------------------------------------------------ #
	# {{ 1 }}
	## converts from bed file to GRanges object
	BedToGRanges = function(bed, values=FALSE, seqlen=NULL){
  
		library(stringr)
		library(GenomicRanges)
		bed = bed[order(bed[,1], bed[,2]),]
  
		if(sum(str_detect(names(bed), 'strand')) == 1){
			strand=as.character(bed[,str_detect(names(bed), 'strand')])
    
			if(!all(strand %in% c('+','-','*')))
				stop('unallowed strand character found')
		}else{
			strand='*'
		}
  
  
		ranges = GRanges(
					seqnames=as.character(bed[,1]), 
					ranges=IRanges(start=as.numeric(bed[,2]), end=as.numeric(bed[,3])),
					strand=strand,
					.id=1:nrow(bed))
 
		if(!is.null(seqlen)){
			gseqlev=GRanges(names(seqlen), IRanges(1, seqlen))
			ranges = subsetByOverlaps(ranges,gseqlev, type='within')
			seqlevels(ranges) = names(seqlen)
			seqlengths(ranges) = seqlen
		}
  
		if(values == TRUE && ncol(bed) > 3){
			col.ids = setdiff(names(bed), c( "seqnames", "ranges", "strand", "seqlengths", "start", "end", "width",  "element", "chr"))
			values(ranges) = bed[values(ranges)$'.id',col.ids]
		}
		values(ranges)$'.id' = NULL
  
		return(ranges)
	}
	
	
	# ------------------------------------------------------------------------------ #
	# {{ 2 }}
	## converts from granges to bed
	GRangesToBed = function(ranges, values=F){
	

		bed = data.frame(
						 chr=as.character(seqnames(ranges)), 
						 start=start(ranges),
						 end=end(ranges), 
						 strand=as.character(strand(ranges))
						 )
		
		if(values == T){
			print('This is going to be slooow...')
			v = as.list(values(ranges))
			v = lapply(v, as.character)
			v = do.call(cbind, v)
			bed = cbind(bed, v)
		}
		attr(bed, 'seqlengts') = seqlengths(ranges)
		return(bed)
	}

	
	
	# ------------------------------------------------------------------------------ #
	# {{ 3 }}
	# reads a gff file to a GRanges object
	GffToGRanges = function(gff, filter=NULL){
	
		library(plyr)
		if(ncol(gff) != 9)
			stop("Number of columns does not match gff format")
		
		if(any(gff[,5] < gff[,4])){
			warning("gff file contains ranges with negative widths...")
			gff = gff[gff[,5] > gff[,4],]
		}

		if(!is.null(filter)){		
			if(filter %in% gff[,3]){
				cat("Filtering", filter, "features...\n")
				gff = gff[gff[,3] == filter,]
			}else{
				stop("The given feature is not present in the gff file")
			}
		}
		
			
		cat('Getting the feature ids...\n')
		s = strsplit(gff$V9, split=';')
		z = sapply(s, length)
		a = split(s, z)
		gff = gff[order(z),]
		l = lapply(a, function(x){
						d = sub('^ ','', unlist(x, use.names=F))
						d = sub('^.+? ','',d)
						m = matrix(d, ncol = length(x[[1]]), byrow=T)
						colnames(m) = sub(' .+$','',sub('^ ','', x[[1]]))
						m})
		ids = rbind.fill(lapply(l, data.frame))
				
		cat('Constructing the granges...\n')
		gff$V7[!gff$V7 %in% c('+','-')] = '*'
		granges = GRanges(seqnames = gff[,1],
						  IRanges(gff[,4],gff[,5]),
						  strand = gff[,7],
						  frame = gff[,8],
						  feature.type = gff[,3],
						  .id = 1:nrow(gff))

		values(granges) = cbind(values(granges), DataFrame(ids)[granges$.id,])
		values(granges)$.id = NULL
		return(granges)			
	}
	
	# ------------------------------------------------------------------------------ #
	# converts a list of characters to a dataframe
	CharList2DF = function(cl, sep="\t"){
	
		l = strsplit(cl, sep)
		if(! length(unique(sapply(l, length))) == 1)
			stop("elements do not have the same length")
		data.frame(do.call(rbind, l))
	}
	
	
	# ------------------------------------------------------------------------------ #
	#formats the ensembl gene files
	FormatEnsembl = function(file, outpath=NULL, outname=NULL, add.chr=TRUE){
	
		header=c('Chromosome Name'	= 'chr',
				 'Gene Start (bp)'	= 'start',
				 'Gene End (bp)'	= 'end',
				 'Exon Chr Start (bp)' = 'ex.start',
				 'Exon Chr End (bp)' = 'ex.end',
				 'Constitutive Exon' = 'const.ex',
				 'Exon Rank in Transcript' = 'ex.rank',
				 'phase' = 'phase',
				 'cDNA coding start' = 'cdna.cod.start',
				 'cDNA coding end' = 'cdna.cod.end',
				 'Genomic coding start' = 'gen.cod.start',
				 'Genomic coding end' = 'gen.cod.end',
				 'Ensembl Exon ID' = 'ens.ex.id',
				 'Strand'			= 'strand',
				 'Ensembl Gene ID'	= 'ens.gene.id',
				 'Ensembl Transcript ID'	= 'ens.trans.id',
				 'Associated Gene Name'		= 'gene.name',
				 'Gene Biotype'				= 'gene.biotype',
				 'Transcript Biotype'		= 'trans.biotype')
		
		cat('Reading in the table...\n')
		if(file.exists(file)){
			d = read.table(file, header=F,sep='\t')
		}else{
			stop('the designated file does not exist')
		}
		if(any(!d[1,] %in% names(header)))
			stop('some of the header names are not present in the lib')
		
		cat('Formatin...\n')
		colnames(d) = header[unlist(d[1,])]
		d = d[-1,]
		if(any(d$chr == 'MT'))
			d$chr[d$chr == 'MT'] = 'M'
		if(add.chr == TRUE){
			ind = str_detect(d$chr, '^\\d+') | d$chr %in% c('X','Y','M')
			d$chr = paste('chr',d$chr,sep='')
		}
		if(sum(str_detect(colnames(d), 'strand')) == 1)
			d$strand = ifelse(d$strand == 1, '+', '-')
		if(sum(str_detect(colnames(d), 'strand')) > 1)
			stop('there is more than 1 strand column')
		
		### sets the start/end colnames if there are no appropriate columns
		### gene -> transcript -> exon priorities
		if(all(colnames(d) != 'start' & colnames(d) != 'end')){
			
			w = c(which(colnames(d) == 'tr.start'),
				  which(colnames(d) == 'tr.end'))
			if(length(w) == 0){
				w = c(which(colnames(d) == 'ex.start'),
					  which(colnames(d) == 'ex.end'))
			}
			if(length(w) == 0)
				stop('there are no appropriate start and end columns')
			
			colnames(d)[w] = c('start','end')
		}
		
		pos.ind = na.omit(match(c('chr','start','end','strand'),colnames(d)))
		id.ind = na.omit(match(c('ens.gene.id','ens.trans.id','ens.ex.id'),colnames(d)))
		d.new = data.frame(d[,pos.ind], d[,id.ind], d[,-c(pos.ind,id.ind)])
		
		cat('Savin...\n')
		if(is.null(outpath) & is.null(outname))
			cat('Will overwrite the original file...\n')
		if(is.null(outpath))
			outpath=dirname(file)
		if(is.null(outname))
			outname=basename(file)
		
		
		outfile=file.path(outpath, outname)
		write.table(d.new, outfile, row.names=F, col.names=T, quote=F, sep='\t')
		cat('\n')
	}
	
	
	# ------------------------------------------------------------------------------ #
	# Gets the
	GRangesID = function(x, sep='_', strand=TRUE){
		if(strand == TRUE){
			paste(as.character(seqnames(x)), 
				  start(x), 
				  end(x), 
				  as.character(strand(x)), sep=sep)
		}else{
			paste(as.character(seqnames(x)), 
				  start(x), 
				  end(x), 
				  sep=sep)
		}
	}
	
	
#/{2} CODE




