### calculation works, but the normalization does not


### calculates the normalization score based on hexamer sequences
	HexamerCorrectionCalculate = function(bam.path){
	
		source('/home/members/vfranke/MyLib/BamWorkers.R')
		library(GenomicRanges)
		library(ShortRead)
		chrs = chrFinder(bam.path)
		
		### construct the heptamer space
		p.base = rep(0, 4^7)
		names(p.base) = mkAllStrings(c("A", "C", "G", "T"), width = 7)
		p.compute = function(i) {
			tab = tables(subseq(bam$seq, start = i, width = 7), n = NULL)$top
			out = p.base
			common.names = intersect(names(out), names(tab))
			out[common.names] = tab[common.names]
			return(out)
		}
		
		### goes throught the chromosomes and counts the number of sequences on given positions in the reads
		p.l = list()
		for(i in chrs[,1]){
			
			print(i)
			g = GRanges(seqnames=i, IRanges(1, chrs[chrs[,1]==i,2]))
			bam = BamReader(bam.path, what=c('rname', 'seq', 'qual'), which.ranges=g)
			p.l[['biased']][[i]] = rowSums(do.call(cbind, lapply(1:2, p.compute)))
			p.l[['unbiased']][[i]] = rowSums(do.call(cbind, lapply(24:29, p.compute)))
			
		}
		### calculates the correction score and returns the values
		cat('Constructing the final data...\n')
		p.biased = Reduce('+', p.l[['biased']])
		p.unbiased = Reduce('+', p.l[['unbiased']])
		p = round(p.unbiased/p.biased,1)
		return(p)
	}
		
	
	WeightGenome = (d, cov, genome, chr){
	
		library(Biostrings)
		dict = PDict(names(d))
		mind = matchPDict(dict, genome[[chr]])
		s = startIndex(mind)
		l = unlist(lapply(s, length))
		v = rep(d, times=l)
		ind = unlist(s)
		cc = cov[[chr]]
		cc = cc[ind] * v
		return(cov)
	}
	
	d = HexamerCorrectionCalculate(file)	
	