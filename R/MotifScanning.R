### Info: Functions for motif Scanning - good part imported from ScanLib
### Date: 05.07.2011.
###	Auth: Vedran Franke


# {1}
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

#/{1}


# {2}
	### writes a PWM/PFM in a jaspar formatted style - can be parsed by PFM parser
	writePWM = function(pwm, file='./', append=T, name=NULL){

		str = apply(pwm, 1, paste, collapse=' ')
		str = paste(rownames(pwm),'[', str, ']')
		### prints the header for the matrix
        if(is.null(name))
            name=basename(file)
		cat('>', name, '\n', file=file, append=T, sep='')
		### prints the rest of the matrix
		### when you get time write this properly
		for(i in 1:length(str)){
			cat(str[i], '\n', file=file, append=append, sep='')
		}	
	}
	
#/{2}



# {3}
	### loops over a list of PFM-s and converts them to a PSSM using the Wasserman & Sandelin formula
	PFMtoPSSM = function(l.PFM, p = NULL){
	
		l.PSSM = list()
		for(i in seq(along=l.PFM)){
		
			PFM = l.PFM[[i]]$PFM
			ID = l.PFM[[i]]$ID[2]
			
			l.PSSM[[ID]] = PSSM.Wasserman.Sandelin(PFM, p=c(A=0.25, C=0.25, G=0.25, T=0.25))
		}
		return(l.PSSM)
	}
#/{3}


# {4}
	### converts pfm to pwm 
	# log2 (n(b,i) + sqrt(N)/4)/(N+sqrt(N)) /(p)
	PSSM.Wasserman.Sandelin = function(pfm, p=c(A=0.25, C=0.25, G=0.25, T=0.25)){

	
		stopifnot(!is.null(p) | !is.null(names(pfm)))
		if(!is.matrix(pfm))
			stop('pfm needs to be a matrix')
		
		### Wasserman & Sandelin NRG 2004
		Nmat = matrix(rep(colSums(pfm), nrow(pfm)), nrow=nrow(pfm), byrow=T)
		p.ind = match(names(p), rownames(pfm))

		### numerator
		m1 = pfm + (sqrt(Nmat)*(1/length(p)))
		### denominator
		m2 = (Nmat + sqrt(Nmat)) * p[p.ind]
		pwm = log2(m1/m2)
		return(pwm)
	}
#/{4}

# {5}
	### calculates the realtive scores for sequences from a pwm and a data.frame of sequences
	PWMScore = function(pwm, seqs){

		seqs = strsplit(seqs, split='')
		if(length(seqs[[1]]) != ncol(pwm)){stop('The length of the sequences does not correspond to the length of the pwm\n')}
		nucleotide.order = rownames(pwm)
		pwm.col.ind = ((1:ncol(pwm))-1)*4
		score.abs = unlist(lapply(seqs, function(x){pwm.row.ind=match(x, nucleotide.order);sum(pwm[pwm.row.ind + pwm.col.ind])}))
		score.max = sum(apply(pwm, 2, max))
		score.rel = score.abs/score.max
		return(data.frame(score.abs=score.abs, score.rel=score.rel))
	}
#/{5}



# {6}
	### scans a set of sequences with a given pwm
	ScanPWM = function(seqs = NULL, pwm = NULL, min.score='80%', get.score=FALSE){

		if(is.null(seqs) | is.null(pwm)){ stop('Please specify the proper input data!\n') }
		library(Biostrings)
		pwm.p = as.matrix(pwm)
		pwm.m = reverseComplement(pwm.p)
		
		### scans the sequences
		# cat('Scanning plus...\n')
		match.p = sapply(seqs, function(x)matchPWM(pwm.p, x, min.score=min.score))
		# cat('Scanning minus...\n')
		match.m = sapply(seqs, function(x)matchPWM(pwm.m, x, min.score=min.score))
		# cat('Scanning done!\n')
		### converts the hits to a data frame and 
		hit.p = XStringToDataFrame(match.p)
		hit.m = XStringToDataFrame(match.m)
		if(get.score == TRUE){
			hit.p = cbind(hit.p, PWMScore(pwm, hit.p$seqs))
			hit.m = cbind(hit.m, PWMScore(pwm, hit.m$seqs))
		}
		if(nrow(hit.p) != 0){
			hit.p$strand = 'p'
		}else{
			hit.p = NULL
		}
		
		
		if(nrow(hit.m) != 0){
			hit.m$strand = 'm'
		}else{
			hit.m = NULL
		}
		
		cat('Returning the data...\n\n')
		return(rbind(hit.p, hit.m))
	}
#/{6}



# {7}
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
#/{7}


# {8}
### takes a set of regions and annotates them by the pwm binding
AnnotateRegionPWM = function(genome, regions, pwm, min.score='80%'){

	chrs = unique(as.character(seqnames(regions)))
	if(!all(chrs %in% seqnames(genome)))
		stop('Not all chromosomes are present in the genome')
	if(!class(regions) == 'GRanges')
		stop('regions must be a valid GRanges object')
	
	pwms = list(p = pwm, m=reverseComplement(pwm))
	ind.match = rep('No', length=length(regions))
	for(i in chrs){
	
		cat(i,'\r')
		ind = which(seqnames(regions) == i)
		regions.chr = ranges(regions[ind])
		v = Views(genome[[i]], regions.chr)
		l.match = lapply(pwms, function(x)matchPWM(x, v, min.score=min.score))
		l.o = lapply(l.match, function(x)as.matrix(findOverlaps(regions.chr, x)))
		reg.hit = unique(unlist(lapply(l.o, function(x)x[,1])))
		ind.match[ind[reg.hit]] = 'Yes'
	}
	cat('\n')
	return(ind.match)
}
#/{8}


# {9}
# takes a meme class object and returns a PFM
BlocksToPFM = function(char){

	if(!class(char) == 'meme')
		stop('char needs to be a character vector')
		
	d = apply(do.call(rbind, strsplit(char, split='')),2, function(x)table(factor(x, levels=c('A','C','G','T'))))
	return(d)
}
#/{9}
