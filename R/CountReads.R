### INFO: R Script
### DATE: 19.11.2014
### AUTHOR: Vedran Franke
rm(list=ls());gc()


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
# source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(stringr)
library(doMC)
library(GenomicAlignments)

#/{1} LIBRARIES


# {2} CODE
	# {{1}} FUNCTIONS
	#/{{1}} FUNCTIONS
	
	
	# {{2}} INPUT VARIABLES 
	
		# {{{1}}} PATH VARIABLES
		inpath = '/data/bioinformatics/Projects/EWyler_Herpes/Data/Mapped/Tophat2'
		
		outpath = '/home/vfranke/Projects/CBirchmeier_enhancer/Documentation/Stats'
		
		#/{{{1}}} PATH VARIABLES
		
		# {{{2}}} SCRIPT PARAMS
		registerDoMC(21)
		#/{{{2}}} SCRIPT PARAMS
		
	#/{{2}} INPUT VARIABLES

	
	# {{3}} MAIN
	
	files = list.files(inpath, recursive=TRUE, pattern='bam$', full.names=TRUE)
	files = files[str_detect(files, 'sorted')]
	files = files[!str_detect(files, 'corrbt')]
	
	l.samp = list()
	for(i in 1:length(files)){
	
		file = files[i]
		name = BamName(file)
		print(name)
		
		chrs = chrFinder(file)	
		l.tab = foreach(chr = chrs$chr)%dopar%{
		
			print(chr)
			w = GRanges(chr, IRanges(1, chrs$chr.len[chrs$chr == chr]))
			g = readGAlignmentsFromBam(file, param = ScanBamParam(which=w, tag='NH'), use.names=TRUE)
			tab = data.table(name=names(g), NH=values(g)$NH)	
		}
		d = unique(rbindlist(l.tab))
		l.samp[[name]] = data.frame(mapped = length(unique(d$name)),
									single = sum(d$NH==1),
									multi = sum(d$NH > 1))
	}
	d.samp = do.call(rbind, l.samp)
	
	write.table(file.path(outpath, 'MappedReads.txt'), row.names=F, col.names=T)
	
	#/{{3}} MAIN
#/{2} CODE



