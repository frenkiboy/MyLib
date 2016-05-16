#{DESCRIPTION}
# INFO: Functions for working with sequences and motifs
# DATE: 1.12.2010
# AUTHOR: v+
#{/DESCRIPTION}


#{CODE}
    # ------------------------------------------------------------------------------------ #
	#{{1}}
	# first order markov model
	# input: character vector or DNAString set
	# output: DNAString set
	MarkovBackground = function(foreground.seq = NULL){
	
		if(is.null(foreground.seq)){stop('Please specify the foreground sequences for: ', match.call()[[1]], " function\n")}
		require(Biostrings)
		cat('Constructing the background sequences...\n')
		background.seq = vector()
		mx = matrix(colSums(oligonucleotideFrequency(x=DNAStringSet(foreground.seq), width=2)), nrow=4, byrow=T)
		colnames(mx) = c('A','C','G','T')
		rownames(mx) = c('A','C','G','T')
		mx = mx/rowSums(mx)
		for (i in seq(along=foreground.seq)){
			cat('seq:',i,'\r')
			seq = foreground.seq[i]
			seqind = rep(NA,nchar(seq))
			seqind[1] = sample(1:4, 1, prob=colSums(mx)/4)
			for (j in 2:nchar(seq)){
				seqind[j] = sample(1:4, 1, prob = mx[seqind[j-1],])
				
			}
			background.seq = c(background.seq, paste(c("A","C","G","T")[seqind], collapse=""))
		}
		cat('Background sequences constructed...\n\n')
		return(DNAStringSet(background.seq))
	}
	#{{/1}}
#{/CODE}

    
	# ------------------------------------------------------------------------------------ #
	#{{2}}
	nucleotideShuffle <- function (sequences) {
	    # check argument
	    if(!inherits(sequences, "DNAStringSet")) stop("sequences not a DNAStringSet")
	    # get widths of sequences
	    widths <- width(sequences)
	    if (min(widths)<max(widths)) stop ("sequences must be of equal width")
	    # turn sequences into a character matrix
	    seq.mx <- matrix(strsplit(paste0(as.character(sequences),collapse=""), "")[[1]], 
	                     byrow = T, ncol = widths[1])
	    shuffled.mx <- apply(seq.mx, 2, sample)
	    return(DNAStringSet(apply(shuffled.mx, 1, paste0, collapse="")))
	    
	}
	
	dinucleotideShuffle <- function (sequences) {
	    # check argument
	    if(!inherits(sequences, "DNAStringSet")) stop("sequences not a DNAStringSet")
	    # get widths of sequences
	    widths <- width(sequences)
	    if (min(widths)<max(widths)) stop ("sequences must be of equal width")
	    # turn sequences into a character matrix
	    seq.mx <- matrix(strsplit(paste0(as.character(sequences),collapse=""), "")[[1]], 
	                     byrow = T, ncol = widths[1])
	    shuffled.mx <- matrix(NA, ncol=ncol(seq.mx), nrow=nrow(seq.mx))
	    shuffled.mx[,1] <- seq.mx[,1]
	    for(i in 2:(ncol(shuffled.mx))) {
	        cat(i, "\n")
	        for (nucl in c("A","C","G","T","N")) {
	            next.nucl <- seq.mx[seq.mx[,i-1]==nucl,i]
	            if(is.vector(next.nucl)) {
	                shuffled.mx[shuffled.mx[,i-1]==nucl,i] <- sample(next.nucl)
	            }
	        }
	    } 
	    return(DNAStringSet(apply(shuffled.mx, 1, paste0, collapse="")))
	    
	}
	
	prepare.mx <- function (pfm) {
	    if(!inherits(pfm, "PFMatrix")) stop("argument not a PFMatrix")
	    pwm.mx <- as.matrix(toPWM(pfm))
	    pwm.mx<- pwm.mx - apply(pwm.mx,2,min)
	    pwm.mx <- pwm.mx / apply(pwm.mx,2,max)
	    return(pwm.mx)
	}
	
	prepare.revcom.mx <- function (pfm) {
	    return(reverseComplement(prepare.mx(pfm)))
	}
	# ------------------------------------------------------------------------------------ #

    # motif clustering
    ## not done!!!
# 	n = length(motifs)
# 	distmat = matrix(0, ncol=n, nrow=n)
# 	strandmat = matrix('', ncol=n, nrow=n)
# 	head(which(distmat != t(distmat),arr.ind=TRUE))
# 	for(m1 in 1:n){
# 	    for(m2 in 1:n){
# 	        
# 	        cat(paste(m1, m2),'\r')
# 	        lmo = list(
# 	            mo1 = motifs[[max(m1,m2)]],
# 	            mo2 = motifs[[min(m1,m2)]])
# 	        
# 	        wind = which.max(c(ncol(lmo$mo1), ncol(lmo$mo2)))
# 	        if(wind == 2)
# 	            lmo = rev(lmo)
# 	        nc = floor(ncol(lmo[[2]])/3)
# 	        lmo[[1]] = cbind(matrix(0.25, ncol=nc, nrow=4),lmo[[1]],matrix(0.25, ncol=nc, nrow=4))
# 	        
# 	        lmo$mo3 = data.frame(reverseComplement(as.matrix(lmo[[2]])))
# 	        
# 	        vm = list()
# 	        for(offset in 1:(ncol(lmo[[1]])-ncol(lmo[[2]])+1)){
# 	            
# 	            hd.p = vector()
# 	            hd.m = vector()
# 	            for(pos in 1:ncol(lmo[[2]])){
# 	                
# 	                hd.p = c(hd.p, HellingerDist(lmo[[2]][,pos], lmo[[1]][,pos+offset-1]))
# 	                hd.m = c(hd.p, HellingerDist(lmo[[3]][,pos], lmo[[1]][,pos+offset-1]))
# 	            }
# 	            fun = function(x)fivenum(x)[3]
# 	            strand = ifelse(fun(hd.p) > fun(hd.m), 'c', 'rc')
# 	            vm$score = c(vm$score, max(c(fun(hd.p), fun(hd.m))))
# 	            vm$strand = c(vm$strand, strand)
# 	        }
# 	        distmat[m1, m2] = 1 - max(vm$score)
# 	        strandmat[m1, m2] = vm$strand[which.max(vm$score)]
# 	    }
# 	}