#{DESCRIPTION}
	#{INFO}: Collection of functions for analysis of ChipSeq data
	#		  goes from raw chip seq .export files to bed formated diff.espressed peaks
	#		 based on chipseq functions from the CTCF project
	#{DATE}: 07.01.2011
	#{BY}:	 v+
#{/DESCRIPTION}



#{LIBRARIES}
# library(rtracklayer)
# library(biomaRt)
# library(ShortRead)
#{/LIBRARIES}

#{1}
#{ LOADER }
	#{{INFO}}: loads the raw chipseq data and saves them in a Rdata object	
	#		   loading from Rdata object is 7 times faster than from *export.txt

	#{{1}}
	# takes a raw chip data and saves the unique and duplicated reads in an Rdata object
	ExportReader = function(inpath, outpath, filter=15){
		
		### checks if the input file has the .txt suffix
		### dynamically creates variable names out of input file names
		cat('Creating the output name...\n')
		outfile = basename(inpath)
		outfile = sub('control', 'cont', outfile)
		if(length(grep('.txt', outfile)) != 0){
			desc.file 		 = file.path(outpath, sub('.txt','.desc.txt', outfile))
			outfile.uniq 	 = sub('.txt', paste('.qfilt.', filter, '.AlignRead.uniq.RData', sep=''), outfile)
			outfile.dupl 	 = sub('.txt', paste('.qfilt.', filter, '.AlignRead.dupl.RData', sep=''), outfile)
			outfile.obj.uniq = sub('.txt', '.uniq', outfile)
			outfile.obj.dupl = sub('.txt', '.dupl', outfile)
		}else{
			desc.file 		 = file.path(outpath, paste(outfile, '.desc.txt', sep=''))
			outfile.uniq 	 = paste(outfile,'.qfilt.', filter, '.AlignRead.uniq.RData', sep='')
			outfile.dupl 	 = paste(outfile,'.qfilt.', filter, '.AlignRead.dupl.RData', sep='')
			outfile.obj.uniq = paste(outfile, '.uniq')
			outfile.obj.dupl = paste(outfile, '.dupl')
		}
		cat('Output name created\n\n')
					
		### reads in the file
		cat('Reading the file...', inpath,"\n")
		cat('Quality filter:', filter, "\n")
		align.filter = alignQualityFilter(filter)	
		chip.reads = compact(readAligned(inpath, filter=align.filter, type="SolexaExport"))
		cat('Chip reads: ', length(chip.reads), "\n", file=desc.file, append=T)
		cat('Reading file done!\n\n')
		
		### finds the duplicated reads
		cat('Finding duplicated reads in chip...\n')
		dupl = duplicated(paste(position(chip.reads), strand(chip.reads), chromosome(chip.reads),sep=";"))
		uniq.chip.reads = chip.reads[!dupl]
		dupl.chip.reads = chip.reads[dupl]
		
		rm('chip.reads')
		gc()
		
		cat('Uniq reads: ', length(uniq.chip.reads), "\n", file=desc.file, append=T)
		cat('Dupl reads: ', length(dupl.chip.reads), "\n", file=desc.file, append=T)
		cat('Finding duplicated reads in chip done!\n\n')
		
		### saves the data in an Rdata.object
		cat('Saving aligned read objects', "\n")
		### saves the variable names to an .RData file
		assign(outfile.obj.uniq, uniq.chip.reads)
		assign(outfile.obj.dupl, dupl.chip.reads)
		save(list=as.character(outfile.obj.uniq), file=file.path(outpath, outfile.uniq))
		save(list=as.character(outfile.obj.dupl), file=file.path(outpath, outfile.dupl))
		cat('Aligned read objects saved!\n')
		gc()
	}

	
	#{{2}}
	### takes eland export file and saves it as an Rdata object
	ExportToRdata = function(inpath, outpath, filter=15){
		
		### reads in the file
		cat('Reading the file...', inpath,"\n")
		cat('Quality filter:', filter, "\n")
		align.filter = alignQualityFilter(filter)	
		chip.reads <- compact(readAligned(inpath, filter=align.filter, type="SolexaExport"))
		cat('Reading file done!\n\n')
		
		### saves the data in an Rdata.object
		cat('Saving aligned read objects', "\t")
		outfile = basename(inpath)
		outfile = sub('control', 'cont', outfile)
		outfile = sub('.txt', paste('.qfilt:', filter, '.AlignRead.RData', sep=''), outfile)
			
		### dynamically creates variable names out of input file names and saves them to an .RData file
		outfile.obj = sub('.fixed.txt', '', outfile)
		assign(outfile.obj, chip.reads)
		save(list=as.character(outfile.obj), file=file.path(outpath, outfile))
		cat('Aligned read objects saved!')
		gc()
	}

	
	#{{3}}
	### exports the coverage for each coverage file in a .bed format
	ExportCoverage = function(cov, outfile, colnum=5, header=FALSE, append=T){
		
		### checks if header exists
		if(append == T){
			if(file.exists(outfile))
				file.remove(outfile)
		}
		### prints the header fo the bedgraph
		if(header != FALSE){
		
			if(class(header) == 'character'){
				cat(paste(header, '\n', sep=''), file=outfile)
			}else{
				cat('Using the default header...\n')
				header = paste('track type=bedGraph', paste('name=',basename(outfile),sep=''))
			}
		}
		### prints the data
		for(chr in sort(names(cov))){
			cat('Exporting:', chr, '\n')
			.ExportCoverageChr(cov[[chr]], outfile, chr, colnum=colnum, append=append)
		}
	}
	
	.ExportCoverageChr = function(cov, outfile, chr, colnum=5, append=T){
	
	
		if(!colnum %in% 4:6)
			stop('Allowed column numbers are: 4,5,6')
			
		### prints the data
        nz = runValue(cov) != 0
		
		
		if(colnum == 6){
			data = cbind(chr,
						  format(start(cov)[nz] - 1, scientific=F),
                          format((end(cov)), scientific=F)[nz],
						  '.',
						  '.',
						  sprintf('%0.3f', runValue(cov)[nz])
						 )
		}
		if(colnum == 5){
			data = cbind(chr,
						  format(start(cov)[nz] - 1, scientific=F),
                          format((end(cov)), scientific=F)[nz],
						  '.',
						  sprintf('%0.3f', runValue(cov)[nz])
						 )
		}else if(colnum == 4){
			data = cbind(chr,
						 format(start(cov)[nz] - 1, scientific=F),
						 format((end(cov)), scientific=F)[nz],
						 sprintf('%0.3f', runValue(cov)[nz])
						)
		}
		
        write.table(data, file=outfile,	append=append, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
	#{{3}}
	
	#{{4}}
	### loads the input files from the standard chipseq analysis pipeline into one big table
	StandardChipseqLoader = function(inpath, ext=0){
		
			source('/export/biggles/home/vedran/R/MyLib/ScanLib.R')
			files = list.files(inpath, full.names=T)
			l.tab = list()
			for(i in seq(along=files)){
				file = files[i]
				name = sub('.chip.uniq.+', '', basename(files[i]))
				name = sub('CHIP.', '', name)
				WorkPrinter(name)
				tab = read.table(file, header=T, comment.char='#')
				tab$start = tab$peak - ext
				tab$end = tab$peak + ext
				tab$rank = (1:nrow(tab)/nrow(tab)) * 100
				tab$set = name
				l.tab[[name]] = tab
			}
			tab = do.call(rbind, l.tab)
			return(tab)
		}

#{/LOADER}



#{2}
#{ NORMALIZATION and PEAK CALLING }
	
	#{{1}}
	# Fills the list with data loaded from r objects
	.ListFiller = function(l, path, type){
		
		### gets the name of the object
		load(path, envir=sys.frame(sys.parent(n=0)))
		obj.name = ls(pattern=type)
		### puts the object into a list
		l[[obj.name]] = get(obj.name)
		### removes the object from the envir frame
		### cleans the memory
		gc(verbose=F)
		return(l)
	}

	
	#{{2}}
	# Splits big peaks into several smaller ones
	.BigPeakSplitter = function(chip.views, chip.cov, cutoff=0.1){

		### loops through all of the big regions
		for(idx in which(width(chip.views) >  2000)){
			### extracts the regions
			region.overlap = window(chip.cov, start(chip.views[idx]),end(chip.views[idx]))
			### slices it and generates smaller regions
			region.views = slice(region.overlap, lower=max(region.overlap)*cutoff)
			if(length(region.views) > 1){
			
				new.views = Views( chip.cov, 
								   start(chip.views[idx]) + start(region.views), 
								   start(chip.views[idx]) + end(region.views)
								  )
				i = rep(TRUE,length(chip.views))
				i[idx] = FALSE
				chip.views = c(chip.views[i], new.views)
			}
		}
		return(chip.views)
	}

	
	#{{3}}
	# normalizes chip seq data by counts - takes two samples and reduces the bigger one to the size of the smaller one; works on whole genomes
	.CountNormGenome = function(chip1, chip2){
		
		require(IRanges)
		### gets the lengths of the data
		chip1.count = length(chip1)
		chip2.count = length(chip2)
		### does the normalization
		cat('Normalizing unique read counts ',paste('(chip1: ',chip1.count, ' chip2: ',chip2.count,')',sep=''), "\n")
		if(chip1.count >= chip2.count){
			uniq.chip1.reads = chip1[sample(chip1.count, chip2.count)]
			uniq.chip2.reads = chip2
			
		}else if(chip1.count < chip2.count){
			uniq.chip1.reads = chip1
			uniq.chip2.reads = chip2[sample(chip2.count, chip1.count)]
		}
		cat('Chip count normalized!\n\n')
		
		
		### returns the list object
		l = list()
		l[[1]] = uniq.chip1.reads	
		l[[2]] = uniq.chip2.reads
		return(l)
	}
	
	
	#{{5}}
	# Finds the centers on flat peaks
	.FlatPeaksCenter = function(chip.peak.loc, chip.views){
		flat = width(viewRangeMaxs(chip.views)) > 1
		chip.peak.loc[flat] = chip.peak.loc[flat] + floor((width(viewRangeMaxs(chip.views))[flat])/2)
		return(chip.peak.loc)		
	}

	
	#{{6}}
	# finds differentially expressed peaks in chip vs control
	# main chipseq analysis function
	# USAGE:
	# export files using the ExportReader to a validly named RData object - type argument used to fish out objects from enviroment to a list
	# the chip and control files need to have tags chip/cont and uniq/dupl in their filenames
	# for each analysis check the corresponding fragment lengths
	PeakFinder = function(chip.uniq.path, cont.uniq.path, outpath= './', lower=4, frag.len=200, genome.name='BSgenome.Mmusculus.UCSC.mm9', type=NA){

		# loads in the genome of interest
		genome = GenomeLoader(genome.name)
			
		# loads the chip uniqe reads
		cat('Loading uniq chip data...\n')
		chip.uniq = list()
		chip.uniq = .ListFiller(chip.uniq, chip.uniq.path, type=type)
		cat('Dataset one loaded!\n\n')

		# loads the uniq cont reads	
		cat('Loading uniq cont data...\n')
		cont.uniq = list()
		cont.uniq = .ListFiller(cont.uniq, cont.uniq.path, type=type)
		cat('Dataset two loaded!\n\n')
		
		# loads the chip dupl reads	
		cat('Loading dupl chip data...\n')
		chip.dupl.path = sub('uniq','dupl', chip.uniq.path)
		chip.dupl = list()
		chip.dupl = .ListFiller(chip.dupl, chip.dupl.path, type=type)
		cat('Dupl dataset loaded!\n\n')
		
		# normalizes the chip data by min. count
		norm.chip.list = .CountNormGenome(chip.uniq[[1]], cont.uniq[[1]])
		names(norm.chip.list) = c(names(chip.uniq), names(cont.uniq))
		
		# calculates the coverage on the normalized data
		norm.chip.list.cov = lapply(norm.chip.list, coverage, width=seqlengths(genome), extend=164L)
		chip.uniq.cov = norm.chip.list.cov[[1]]
		cont.uniq.cov = norm.chip.list.cov[[2]]
		
		# finds the coverage for each strand sepparately
		chip.p = coverage(chip.uniq[[1]][strand(chip.uniq[[1]]) == '+'], width=seqlengths(genome), extend=-35L)
		chip.m = coverage(chip.uniq[[1]][strand(chip.uniq[[1]]) == '-'], width=seqlengths(genome), extend=-35L)
		chip.dupl.cov = coverage(chip.dupl[[1]], width=seqlengths(genome), extend=164L)

		
		### stores the peak data
		peaks = data.frame()
		### stores duplication rates for the peaks
		duplication.rates = vector()
		for (chr in names(chip.uniq.cov)){
			
			cat('Working on chr:', chr,"\n")
			### Identifying enriched regions
			chip.views = slice(chip.uniq.cov[[chr]], lower=lower)
			
			### Splitting large enrichment regions
			chip.views = .BigPeakSplitter(chip.views=chip.views, chip.cov=chip.uniq.cov[[chr]])
			
			# views in the control dataset
			control.views = Views(cont.uniq.cov[[chr]], start=start(chip.views),end=end(chip.views))
			
			### Finding peaks
			chip.peak.loc		= viewWhichMaxs(chip.views) 
			chip.peak.val 		= as.numeric(chip.uniq.cov[[chr]][chip.peak.loc])
			control.peak.val 	= as.numeric(cont.uniq.cov[[chr]][chip.peak.loc])
			enriched 			= chip.peak.val > control.peak.val
			chip.peak.loc 		= chip.peak.loc[enriched]
			chip.peak.val 		= chip.peak.val[enriched]
			control.peak.val 	= control.peak.val[enriched]
			chip.views 			= chip.views[enriched]
			
			if(length(chip.peak.loc)>0){
				
				### Filtering peaks based on strand distribution
				region.fs = Views(chip.p[[chr]], start=chip.peak.loc-frag.len, width=frag.len)
				region.rs = Views(chip.m[[chr]], start=chip.peak.loc, width=frag.len)
				region.fs.rc = viewSums(region.fs)
				region.rs.rc = viewSums(region.rs)
				diff = region.fs.rc>(region.rs.rc*10) | region.rs.rc>(region.fs.rc*10) 
				strand.flagged = rep(FALSE,length(chip.peak.loc))
				strand.flagged[diff] = TRUE
				
				### Choosing center on flat peaks 
				flat = width(viewRangeMaxs(chip.views)) > 1
				chip.peak.loc[flat] = chip.peak.loc[flat] + floor((width(viewRangeMaxs(chip.views))[flat])/2)
				
				### Registering duplication level
				ucov = chip.uniq.cov[[chr]]
				dcov = chip.dupl.cov[[chr]]
				duplication.rates = c(duplication.rates, (as.integer(dcov[chip.peak.loc]) + 1) / (as.integer(ucov[chip.peak.loc]) + 1))
				
				### exp. diff. calculation
				control.100.views	= Views(cont.uniq.cov[[chr]], chip.peak.loc-100,	chip.peak.loc+100)
				control.200.views	= Views(cont.uniq.cov[[chr]], chip.peak.loc-200,	chip.peak.loc+200)
				control.500.views 	= Views(cont.uniq.cov[[chr]], chip.peak.loc-500,	chip.peak.loc+500)
				control.1000.views 	= Views(cont.uniq.cov[[chr]], chip.peak.loc-1000,	chip.peak.loc+1000)
				
				null.chromosome = sum(cont.uniq.cov[[chr]])/seqlengths(genome)[chr]
				null.peak.loc = control.peak.val
				null.100 	= viewMaxs(control.100.views)
				null.200 	= viewMaxs(control.200.views,na.rm=TRUE)
				null.500 	= viewSums(control.500.views,na.rm=TRUE) /unique(width(control.500.views))
				null.1000 	= viewSums(control.1000.views,na.rm=TRUE)/unique(width(control.1000.views))
				lambda 		= pmax(null.chromosome, null.peak.loc , null.100, null.200, null.500, null.1000)
				
				p.values = ppois(chip.peak.val,lambda, lower.tail = FALSE)
				fdr = p.adjust(p.values, method="fdr")
							
				# Fold change calculation
				fold.change = (chip.peak.val+1)/(control.peak.val+1)
				
				# Output table construction
				peaks = rbind(peaks, 
									data.frame(chr=chr, 
												start		= as.integer(start(chip.views)), 
												end			= as.integer(end(chip.views)), 
												peak		= as.integer(chip.peak.loc), 
												chip		= chip.peak.val, 
												control		= control.peak.val, 
												p.value		= p.values,
												fdr			= fdr,
												fold.change	= fold.change, 
												duplicates	= as.integer(dcov[chip.peak.loc]),
												strand.bias	= strand.flagged
											   )
							 )
			}
		}
		peaks = peaks[order(peaks$p.value),]
		peaks$duplication.bias = duplication.rates > (mean(duplication.rates)+(sd(duplication.rates)*5))
	
		cat('Writing the output table...\n')
		outname = paste('CHIP.', names(chip.uniq),'_CONTROL.', names(cont.uniq), '.txt', sep='')
		
		cat('Output name:', outname, "\n")
		write.table(peaks, file.path(outpath, outname), row.names=F, col.names=T, sep="\t", quote=F)
		
		cat('Output table written...bye!!!\n')	
	}
#{/ NORMALIZATION and PEAK CALLING}


#{3}
#{MISCHELLANEOUS}


	#{{1}}
	# expands a region around a peak and takes a subset of the data
	PeakRegionMaker = function(tab, n.start=1, n.end=0, region.expand=100, name='outname'){

		### generates new path
		cat('Creating the new path...\n')
		new.path = paste(name, '.ns.', n.start, '.ne.', n.end, '.reg.exp.', region.expand, '.txt', sep='')
		
		### expands a region around the peaks
		cat('Expanding the regions:', region.expand,'\n')
		tab$start	= tab$peak - region.expand
		tab$end		= tab$peak + region.expand
		
		### checks if the wanted start is larger than the number of rows
		if(n.start >= nrow(tab))(stop('Starting row is bigger than the number of rows!\n'))
		### checks if the wanted size is bigger than the file
		if(n.end > nrow(tab)){n.end = nrow(tab)}
		### if end smaller than start then end = start + end
		if(n.end < n.start & n.end != 0){n.end = n.start + n.end}
		
		cat('Taking the subset:\n\tstart: ',n.start,'\n\tend:', n.end, '\n')
		### takes a subset of the table
		if(n.end == 0){
			tab=tab[n.start:nrow(tab),]
		}else{
			tab=tab[n.start:n.end,]
		}
		
		
		### list for the results
		cat('Returning the values...\n \tlist$path\n \tlist$table\n\n ')
		l = list()
		l$path = new.path
		l$table = tab
		return(l)
	}

	
	#{{2}}
	# wrapper for peak region maker 
	# automaticly downloads the sequences
	PeakRegionsToFasta = function(path='./', pattern='.txt', n.start=1, n.end=0, region.expand=100, genome, sort.by=NULL, sort.ord='up'){

		suppressMessages(source('/export/biggles/home/vedran/R/MyLib/FileLoader.R'))
		genome = GenomeLoader(genome)
		library(stringr)
		if(str_detect(path, '\\....$')){
			region.files = path
		}else{
			region.files = list.files(path, recursive=F, full.names=T, pattern=pattern)
		}
		l.reg = list()
		for(i in seq(along=region.files)){
			reg.name = region.files[i]
			cat('Working on:', basename(reg.name), '\n')
			tab = HugeFileLoader(reg.name, header=HeaderChecker(reg.name), nrows=1000)
			### sorts the regions by a given column
			if(!is.null(sort.by)){
				sort.ord = ifelse(sort.ord == 'up', 'FALSE', 'TRUE')
				ord.column = which(names(tab)) == 'sort.by'
				tab = tab[order(tab[,ord.column], decreasing=sort.ord), ]
			}
			
			
			### takes the subset and expands the regions
			l.reg[[i]] = PeakRegionMaker(
										 tab=tab,  
										 n.start=n.start, 
										 n.end=n.end, 
										 region.expand=region.expand,
										 name = basename(reg.name)
										 )
			names(l.reg)[i] = sub('....$','',basename(reg.name))
			cat('Downloading the sequences..\n\n')
			### downs the sequences
			l.reg[[i]]$seq = DNAStringSet(getSeq(x=genome, names=l.reg[[i]]$table[,1], start=l.reg[[i]]$table[,2], end=l.reg[[i]]$table[,3]))
			### prints the names of the object returned
			names(l.reg[[i]]$seq) = 1:length(l.reg[[i]]$seq)
		}
		cat('Returning the values...\n')
		return(l.reg)
	}

	
	#{{3}}
	# takes a dataframe of overlapping peaks and find she mean rank for each row
	PeakScoreToRank = function(data){
		
		for(i in 1:ncol(data)){
			data = data[order(data[,i], decreasing=T), ]
			data[,i] = 100*(1:nrow(data))/nrow(data)
		}
		score = rowMeans(data)
		return(score)
	}

	
	#{{4}}
	# takes the read aligned object and filters it based on bed ranges
	ReadAlignedBedFilter = function(reads, bed){

		# constructs a data.frame from ReadAligned object
		cat("Constructing a data.frame from ReadAligned object...\n")
		read.bed = data.frame(
								chr=chromosome(reads),
								start=position(reads), 
								end=position(reads)
							)
		# defines the unique key for each read
		read.bed = BedKey(read.bed)
		cat("Data frame constructed...\n\n")
		
		# reduces the number of ranges in the bed set
		cat("Reducing the .bed data frame...\n")
		bed = bed[bed[,1] %in% read.bed[,1],]
		bed.red = reduce(GRanges(IRanges(start=bed[,2], end=bed[,3]), seqnames=bed[,1]))
		bed.red = as.data.frame(bed.red)
		cat(".bed reduced\n\n")
		
		# constructs the overlaps between the sets
		o = FeatureOverlapBED(read.bed, bed)
		read.ind = which(!(read.bed$key %in% o$key.a))
		read.subset = reads[read.ind]
		
		return(read.subset)
	}

	
	#{{5}}
	### takes a list of bed file regions and finds overlapping regions between them
	#l = list of bed files
	#region.expand = region to expand around peaks for overlaps 
	# attrib = column taken as the representative
	SetOverlapFinder = function(l, region.expand=100, attrib){
		
		require(plyr)
		source('/export/biggles/home/vedran/R/MyLib/ScanLib.R')

		for(i in seq(along=l)){	
			data = l[[i]]
			data$start = data$peak - region.expand
			data$end = data$peak + region.expand
			name = names(l)[i]
			
			if(i == 1){
				peaks.ovlap = data[,1:4]
				peaks.ovlap[[name]] = data[[attrib]]
				
			}else{
				### creates the key for constructed peaks
				peaks.ovlap$key = paste(peaks.ovlap[,1], peaks.ovlap[,2], peaks.ovlap[,3])
				### creates the key for new regions
				data$key = paste(data[,1], data[,2], data[,3])
				### does the overlap
				ovlap = FeatureOverlapBED(peaks.ovlap, data)
				### removes from the new set overlapped peaks
				tmp.ind = !data$key %in% ovlap$key.b
				peaks.ovlap.tmp = data[tmp.ind, 1:4]
				peaks.ovlap.tmp[[name]] = data[tmp.ind, attrib]
				### orders the peaks
				ovlap = ovlap[order(ovlap[[paste(attrib, '.b', sep='')]]),]
				### takes the highest peak overlapped in new set
				ovlap = ovlap[!duplicated(ovlap$key.a),]
				### get the corresponding ranks from the new set
				peaks.ovlap[[name]] = 0
				ind = match(ovlap$key.a, peaks.ovlap$key)
				peaks.ovlap[[name]][ind] = ovlap[[paste(attrib, '.b', sep='')]]
				### takes the new peak position as peaks of overlapped regions/2
				peaks.ovlap$peak[ind] = ceiling((peaks.ovlap$peak[ind] + ovlap$peak.b)/2)
				### joins the new and the old set
				peaks.ovlap = rbind.fill(peaks.ovlap, peaks.ovlap.tmp)
				### removes the key
				peaks.ovlap = subset(peaks.ovlap, select=-key)
				### sets all NA values to 0
				peaks.ovlap[is.na(peaks.ovlap)] = 0
				###constructs new peak positions
				peaks.ovlap$start = peaks.ovlap$peak - region.expand
				peaks.ovlap$end = peaks.ovlap$peak + region.expand
			}
		}
		return(peaks.ovlap)
	}
	
	
	#{{6}}
	# loads the corresponding genome in the memory
	GenomeLoader = function(genome){
		
		require(genome, character.only=T)
		genome.name = unlist(strsplit(genome, split='\\.')) 
		return(get(genome.name[2]))		
	}
	
	
	#{{7}}
	# takes a bed and coverage rle list, and returns a dataframe of counts around the bed regions
	# it implies that all ranges have the same size
	BEDCoverageToDataFrame = function(bed = NULL, cov.rle = NULL){
		
		library(stringr)
		### chscks if the data is passed to the function
		if(is.null(bed) | is.null(cov.rle)){ stop('Please specify the proper input data!\n') }
			
		### chromosome names
		chrs = unique(as.character(bed[,1]))
		chrs = chrs[chrs %in% names(cov.rle)]
		### proper colnames for the first three columns of the bed file
		names(bed)[1:3] = c('chr','start','end')
		### constructs the id column
		bed$id = 1:nrow(bed)
				
		### gets the region size - has to be equal in all ranges
		region = unique((bed$end - bed$start) / 2)
		if(length(region) > 1){ stop('Regions are of differing sizes!\n') }
		region = -region:region
		region.size = length(region)
		
		l.d = list()
		for(i in seq(along=chrs)){
			chr = chrs[i]
			bed.chr = bed[bed$chr == chr,]
			cat('\t', chr, ':', nrow(bed.chr), '\n')
			
			### gets the coverage vector for the chromosome
			cov.chr = cov.rle[[chr]]
			
			### get the values around the peaks
			v = Views(cov.chr, start=bed.chr$start, end=bed.chr$end)
			a = as.vector(viewApply(v, as.vector))
			### constructs the coordinates
			pos = rep(region, times=nrow(bed.chr))
			
			### gets the peak index
			peak.id = rep(bed.chr$id, each=length(region))
			
			positive.value = a > 0
			d = data.frame(pos=pos[positive.value], ind=peak.id[positive.value], value=a[positive.value])
			l.d[[chr]] = d
		}
		cat('Constructing the final data.frame...\n')
		data = do.call(rbind, l.d)
		
		cat('Returning the data\n')
		return(data)
	
	}
	
	
	#{{8}}
	### takes a bed file of chipseq peaks and constructs a nonoverlapping set by taking the region that has the highest -10log10 pvalue, out of overlapping ones
	ChipSeqRangeUnion = function(bed, value='p.value'){
		
		cat('Disjoining the set...\n')
		### orders the data by chromosome and start
		bed = bed[order(bed[,1], bed[,2]),]
		
		r = GRanges(seqnames=bed[,1], ranges=IRanges(start=as.integer(bed[,2]), end=as.integer(bed[,3])))
		red = IRanges::reduce(r)
		
		### checks if there were any overlapping regions
		if(length(r) != length(red)){
			
			o = as.matrix(findOverlaps(red, r))
			tab = table(o[,1])
			### finds regions that are overlapping
			o.red = o[o[,1] %in% as.numeric(names(tab))[tab > 1],]
			bed.red = bed[o.red[,2],]
			bed.red$id = o.red[,1]
			### selects the region with the highest score (-10log10 p.value)
			bed.red = bed.red[order(-bed.red[[value]]),]
			bed.red = bed.red[!duplicated(bed.red$id), -ncol(bed.red)]
			return(rbind(bed[o[!(o[,2] %in% o.red[,2]),2],], bed.red))
		}else{
			return(bed)
		}
	}
	


    #{{9}}

   CombineCoverage= function(x,y,fun='-'){
       
       
   }

#{/ MISCHELLANEOUS}
