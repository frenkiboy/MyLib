### INFO: Functions for calling external apps from R
### DATE: 22.03.2013.
### AUTH: Vedran Franke

# {1} LIBRARIES
library(Rsamtools)
#/{1} LIBRARIES



# {2} CODE

	# ------------------------------------------------------------------------------- #
	# {{1}} Blat - calls blat from R
	# genome - path to the genome sequence in fasta or 2 bit format
	# seq - path to the sequence pattern file
	# outpath - where to ouput the blat file
	# name - name of the output file
	# blat - path to the blat binary
	# delete.temp - boolean, to delete the temporary directory
	
	# returns a data frame with the corresponding regions
	
	Blat = function(genome, seq, outpath, name, blat=NULL, delete.temp=FALSE, extra.args=''){
	
		if(is.null(blat))
			blat = '/common/USERS/vfranke/bin/UCSC_tools/BLAT/blat'
		
		outfile = file.path(outpath, name)
		command = paste(blat, extra.args, genome, seq, outfile)
		cat('Blatting ...\n')
		system(command, wait=TRUE)
		s = scan(outfile, what='character', sep='\n')
		header =  c('match', 'mismatch','rep.match', 'Ns', 'Q.gap.count','Q.gap.bases','T.gap.count','T.gap.bases','strand','Q.name', 'Q.size', 'Q.start', 'Q.end','T.name','T.size','T.start','T.end','block.count','blockSizes','qStarts','tStarts')
		d = data.frame(do.call(rbind, strsplit(tail(s,-4),'\t')))
		if(nrow(d) == 0){
			cat('0 hits found\n')
			return()
		}
		colnames(d) = header
		d[,-c(9,10,14,19,20,21)] = apply(d[,-c(9,10,14,19,20,21)], 2, as.numeric)
		
		if(delete.temp)
			file.remove(outfile)
		return(d)
	}
	
	# ------------------------------------------------------------------------------- #
	# {{2}} Clustal Omega
	ClustalO = function(infile, outfile, ClustalO=NULL,threads=12, what='AA', force=T, format='fa'){
	
		if(!what %in% c('AA','DNA'))
			stop('can only align DNA or AA')
		if(is.null(ClustalO))
			ClustalO='/common/USERS/vfranke/bin/ClustalO/clustalo-1.1.0-linux-64'
		
		### checks whether the file exists and whether to force the outfile
		### if the file does exist and the force is off he reads the file
		if((file.exists(outfile) & force==T) | !file.exists(outfile)){
		
			cat('Running the alignmnent...\n')
			threads= paste('--threads=', threads,sep='')
			format = paste('--outfmt=', format,sep='')
			command = paste(ClustalO, '-i', infile, '-o', outfile,'--force' ,threads, format)
			system(command)
		}
		
		cat('Returning the results...\n')
		if(what == 'AA')
			a = readAAStringSet(outfile, format='fasta')
		if(what == 'DNA')
			a = readDNAStringSet(outfile, format='fasta')
		return(a)
	}
	
#/{2} CODE


