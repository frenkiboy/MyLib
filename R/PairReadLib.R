# ------------------------------------------------------- #
# reads in the pair end bam file
GetReads = function(bam.file, regions=NULL){
	
		chrs = chrFinder(bam.file)
		chrs = chrs[!chrs$chr %in% c('chrM', 'chrY'),]
		l.reads = foreach(chr = chrs$chr)%dopar%{
	
			l.stats = list()
			print(chr)
			chrlen = chrs[chrs$chr == chr,]
			which=GRanges(chrlen$chr, IRanges(1, chrlen$chr.len))
			cat('Loading in the bam file...\n')
			
			param = ScanBamParam(flag = scanBamFlag(isPaired=TRUE, isProperPair=TRUE), tag='NH', which=which)
			bam = readGAlignmentsFromBam(bam.file, use.names=T, param = param)
			if(!is.null(regions)){
				n = unique(names(bam[countOverlaps(bam, regions)>0]))
				bam = bam[names(bam)%in%n]
			}
			bam = bam[order(names(bam))]
			s = seq(1, length(bam),2)
			g1=bam[s]
			g2=bam[-s]
			if(!all(names(g1) ==names(g2)))
				stop('names not equal')
			return(list(g1= g1, g2 = g2))
		}
		g1 = do.call(c, lapply(l.reads, function(x)as(x$g1, 'GRangesList')))
		g2 = do.call(c, lapply(l.reads, function(x)as(x$g2, 'GRangesList')))
		list(g1=g1,g2=g2)
}


# ------------------------------------------------------- #
# Counts the number of transcripts per exon