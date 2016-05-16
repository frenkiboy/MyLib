
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
