### loads an RData file and assigns it to the name variable
Assigner=function(`.path`, `.name`){
	if(! is.character(`.path`) | !is.character(`.name`))
		stop('Both arguments should be characters!')

	load(`.path`)
	assign(`.name`, get(ls()[1]), parent.frame())
}

#####----------------------------------------------------------------------------------------------#####


#####----------------------------------------------------------------------------------------------#####
# loads the corresponding genome in the memory
GenomeLoader = function(genome){

	require(genome, character.only=T)
	genome.name = unlist(strsplit(genome, split='\\.'))
	return(get(genome.name[2]))
}
#####----------------------------------------------------------------------------------------------#####


#####----------------------------------------------------------------------------------------------#####
### Loads in large files
HugeFileLoader = function(path, sep = "\t", skip = 0, header = T, nrows = 10, quote=''){

	### counts the number of lines using shell wc command, and converts the output to numeric
	line.count = paste("wc -l ", path, sep = "")
	row.count = as.numeric(strsplit(system(line.count, intern = T), split = " ")[[1]][1]) - skip

	### reads in first five lines of the file and determines the type of each column
	first5rows = read.table(path, header = TRUE, nrows = nrows, skip = skip, sep = sep, stringsAsFactors=F, quote=quote)
	tab.classes = sapply(first5rows, class)

	### reads in the data
	tab = read.table(path, header = header, colClasses = tab.classes, comment.char = "#", nrows = row.count, skip = skip, sep = sep, stringsAsFactors=FALSE)
	return(tab)
}
#####----------------------------------------------------------------------------------------------#####




#####----------------------------------------------------------------------------------------------#####
### reads data downloaded from the ucsc genome browser - all colum tables
# clean=T parameter designates removal of the random and Un chromosomes

UCSCEnsReader = function(path, clean=T){

	# path='/common/USERS/vfranke/Base/Ensembl/mm9/mm9.ens.genes'
	tab = read.table(path, header=HeaderChecker(path), quote="", sep='\t', stringsAsFactors=F)
	ex.start = strsplit(as.character(tab[,10]), split=',')
	ex.end  = strsplit(as.character(tab[,11]), split=',')
	len = unlist(lapply(ex.start, length))
	l.data=list()
	l.data$exons = data.frame(
							  chr 	   = rep(tab[,3], times=len),
							  ex.start = as.numeric(as.character(unlist(ex.start)))+1,
							  ex.end   = as.numeric(as.character(unlist(ex.end))),
							  strand   = rep(tab[,4],times=len),
							  ens.trans.id = rep(tab[,2],times=len),
							  rank = unlist(sapply(len, function(x)1:x)),
							  stringsAsFactors=F
							)


	l.data$genes = data.frame(
							  chr 		= tab[,3],
							  tx.start  = as.numeric(as.character(tab[,5]))+1,
							  tx.end 	= as.numeric(as.character(tab[,6])),
							  strand	= as.character(tab[,4]),
							  cds.start = as.numeric(as.character(tab[,7]))+1,
							  cds.end 	= as.numeric(as.character(tab[,8])),
							  exon.num = len,
							  ens.trans.id = tab[,2],
							  ens.gene.id  = tab[,13],
							  stringsAsFactors=F
							 )
	l.data$exons$ex.id = 1:nrow(l.data$exons)
	l.data$exons$ex.width = l.data$exons$ex.end - l.data$exons$ex.start
	l.data$genes$tx.width = l.data$genes$tx.end - l.data$genes$tx.start

	if(clean == T)
		l.data = UCSC.hg.ChrCleaner(l.data)

	attr(l.data$exons, 'coord') = 1
	attr(l.data$genes, 'coord') = 1

	return(l.data)
}

# head.ens = function(x){
	# lapply(x, head)
# }

### reads in ucsc refseq files
UCSCRefReader = function(path, clean=T){

	# path='/common/USERS/vfranke/Base/Ensembl/mm9/mm9.ens.genes'
	tab = HugeFileLoader(path, nrows=100, header=F)
	ex.start = strsplit(as.character(tab[,10]), split=',')
	ex.end  = strsplit(as.character(tab[,11]), split=',')
	len = unlist(lapply(ex.start, length))
	l.data=list()
	l.data$exons = data.frame(
							  chr 	   = rep(tab[,3], times=len),
							  ex.start = as.numeric(as.character(unlist(ex.start)))+1,
							  ex.end   = as.numeric(as.character(unlist(ex.end))),
							  strand   = rep(tab[,4],times=len),
							  ref.trans.id = rep(tab[,2],times=len),
  							  rank = unlist(sapply(len, function(x)1:x)),
							  stringsAsFactors=F
							)

	l.data$genes = data.frame(
							  chr 		= tab[,3],
							  tx.start  = as.numeric(as.character(tab[,5]))+1,
							  tx.end 	= as.numeric(as.character(tab[,6])),
							  strand	= as.character(tab[,4]),
							  cds.start = as.numeric(as.character(tab[,7]))+1,
							  cds.end 	= as.numeric(as.character(tab[,8])),
							  exon.num = len,
							  ref.trans.id = tab[,2],
							  ref.gene.id  = tab[,13],
							  stringsAsFactors=F
							 )
	l.data$exons$ex.id = 1:nrow(l.data$exons)
	l.data$exons$width = l.data$exons$ex.end - l.data$exons$ex.start
	l.data$genes$width = l.data$genes$tx.end - l.data$genes$tx.start


	if(clean == T)
		l.data = UCSC.hg.ChrCleaner(l.data)

	return(l.data)
}

###reads ucsc genome annotation
### reads in ucsc refseq files
UCSC2Reader = function(path, clean=T){

	tab = HugeFileLoader(path, nrows=100, header=F)
	tab[,11][tab[,11] == ''] = NA
	ex.start = strsplit(as.character(tab[,9]), split=',')
	ex.end  = strsplit(as.character(tab[,10]), split=',')
	len = unlist(lapply(ex.start, length))
	l.data=list()
	l.data$exons = data.frame(
							  chr 	   = rep(tab[,2], times=len),
							  ex.start = as.numeric(as.character(unlist(ex.start)))+1,
							  ex.end   = as.numeric(as.character(unlist(ex.end))),
							  strand   = rep(tab[,3],times=len),
							  ucsc.trans.id = rep(tab[,12],times=len),
  							  rank = unlist(sapply(len, function(x)1:x)),
							  stringsAsFactors=F
							)

	l.data$genes = data.frame(
							  chr 		= tab[,2],
							  tx.start  = as.numeric(as.character(tab[,4]))+1,
							  tx.end 	= as.numeric(as.character(tab[,5])),
							  strand	= as.character(tab[,3]),
							  cds.start = as.numeric(as.character(tab[,6]))+1,
							  cds.end 	= as.numeric(as.character(tab[,7])),
							  exon.num = len,
							  ref.trans.id = tab[,12],
							  ref.gene.id  = tab[,11],
							  stringsAsFactors=F
							 )
	l.data$exons$width = l.data$exons$ex.end - l.data$exons$ex.start
	l.data$genes$width = l.data$genes$tx.end - l.data$genes$tx.start


	if(clean == T)
		l.data = UCSC.hg.ChrCleaner(l.data)

	return(l.data)
}


### cleans the Un and Random chromosomes
UCSC.hg.ChrCleaner = function(l){

	l = lapply(l, function(x)x[grep('chr(\\d+|X|Y)$',x$chr),])
	return(l)
}


#####----------------------------------------------------------------------------------------------#####


#####----------------------------------------------------------------------------------------------#####
### reads a gff file into a bed formatted data frame
GffToBed = function(gff.path){

	gff = read.table(gff.path, stringsAsFactors=F)
	peak = unlist(lapply(strsplit(gff[,ncol(gff)], split=':'), '[', 2))
	bed = data.frame(
					 chr    = as.character(gff[,1]),
					 start  = as.numeric(gff[,4]),
					 end    = as.numeric(gff[,5]),
					 peak   = trunc(as.numeric(peak)),
					 score  = as.numeric(gff[,6]),
					 strand = gff[,7],
					 stringsAsFactors=F
					)
	return(bed)

}
#####----------------------------------------------------------------------------------------------#####



#####----------------------------------------------------------------------------------------------#####
### reads a gff file into a bed formatted data frame


BigToRleLoader = function(bed.file, genome, chr.remove=TRUE, one.start=TRUE){

	source('/home/members/vfranke/MyLib/ChipSeq.Functions.R')
	library(stringr)
	message('Reading in the wig file...')
	h = HugeFileLoader(bed.file, header =F)

	message('Getting chromosome lengths...')
	genome = GenomeLoader(genome)
	chrlen = seqlengths(genome)

	chrs = unique(h[,1])
	if(chr.remove == TRUE){
		message('Removing random chromosomes...')
		chrs = chrs[!str_detect(chrs,'random')]
	}
	if(one.start == TRUE)
		h[,2] = h[,2] + 1


	l.r = list()
	# goes chr by chr and constructs a RleList of values
	for(chr in chrs){

		message(chr)
		r = Rle(0, chrlen[names(chrlen) == chr])
		ind = h[,1] == chr
		ranges = IRanges(h[ind,2], h[ind,3])
		ranges = ranges[end(ranges) < length(r)]
		r[ranges] = rep(h[ind,4], times = width(ranges))
		l.r[[chr]] = r
	}
	l.r = RleList(l.r)
	return(l.r)
}
#####----------------------------------------------------------------------------------------------#####



#####----------------------------------------------------------------------------------------------#####
## reads swembl output to R
SwemblReader = function(path){

		tab = read.table(path, header=T, sep='\t')
		names(tab) = c('chr','start','end','chip','length','unique','score','cont','max.coverage','peak')
		tab$peak = trunc(tab$peak)
		return(tab)
	}
#####----------------------------------------------------------------------------------------------#####


#####----------------------------------------------------------------------------------------------#####
### reads a repeat masker dataframe into a GRangesObject
RepeatMaskerToGRanges = function(path){
	require(GenomicRanges)
	tab = HugeFileLoader(path, header=T)
	g = GRanges(seqnames   = tab$genoName,
				ranges     = IRanges(tab$genoStart+1, tab$genoEnd),
				strand     = tab$strand,
				rep.class  = tab$repClass,
				rep.family = tab$repFamily)
	return(g)
}
#####----------------------------------------------------------------------------------------------#####


#####----------------------------------------------------------------------------------------------#####
# cisgenome reader
CisgenomeReader = function(path){

		header = read.table(path, header=F, sep='\t', nrows=1, comment.char='')
		header[2] = 'chr'
		tab = read.table(path, header=F, sep='\t')
		colnames(tab) = header
		tab = tab[,-1]
		tab
}
#####----------------------------------------------------------------------------------------------#####


# ---------------------------------------------------------------------------- #
read_SJ = function(
  path
){
  suppressPackageStartupMessages(library(data.table))
  files = list.files(path, full.names=TRUE, pattern='SJ', recursive=TRUE)
  cnams = c('chr','start','end','strand','motif','annot','uniq','mult','overhang')
  lf = list()
  for(i in files){
    name = sub('SJ.out.tab','',basename(i))
    name = sub('\\W$','',name)
    message(name)
    tab  = data.table::fread(i)
    setnames(tab, cnams)
    tab[,sample := name]
    lf[[name]] = tab
  }
  dtab = rbindlist(lf)
  return(dtab)
}

