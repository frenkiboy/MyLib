# ---------------------------------------------------------- #
ParseGenbank = function(gbfile, genome.name=NULL){

	if(!file.exists(gbfile))
		stop('File has to exist')
	
	library(stringr)
	library(GenomicRanges)
	gb = scan(gbfile, what='character', sep='\n')
	
	chr.s = which(str_detect(gb, 'LOCUS'))
	chr.e = which(str_detect(gb, '//'))

	chrnames = do.call(rbind, strsplit(gb[chr.s],'\\s+'))
	chrnames = data.frame(do.call(rbind, strsplit(chrnames[,2], ':')), chrnames[,-c(1:2)])
	chrs = sapply(strsplit(chrnames[,3], '_') ,'[[', 1)
	if(!is.null(genome.name))
		chrs = str_replace(chrs, genome.name, '')
	chrnames[,3] = chrs
	seqinfo = Seqinfo(seqnames   = chrnames[,3],
					  seqlengths = as.numeric(chrnames[,5]),
					  isCircular = chrnames[,9] != 'linear',
					  genome     = ifelse(!is.null(genome.name), genome.name, NULL))
	
	lf = list()
	ls = list()
	for(i in 1:length(chr.s)){
	
		chrname = seqnames(seqinfo)[i]
		print(chrname)
		chr = gb[chr.s[i]:chr.e[i]]
		feat = ParseFeatures(chr, chrname)
		seq = ParseSequence(chr)
		lf[[chrname]] = feat
		ls[[chrname]] = seq
	}
	
	fnames = unique(unlist(lapply(lf, names)))
	feat = list()
	for(i in 1:length(fnames)){
	
		fname = fnames[i]
		print(fname)
		a = lapply(lf, '[[',fname)
		a = a[sapply(a, length) > 0]
		feat[[fname]] = unlist(GRangesList(a))
	}
	
	seq = DNAStringSet(unlist(ls))
	return(list(feat=feat, seq=ls))	
}


# ---------------------------------------------------------- #
### extracts the genes, metadata and protein sequences
### ALERT! - the function has a problem because if a gene has 2 cds it will only count one 

ParseFeatures = function(s, chrname=NULL){
  
	###loads the libraries
	library(stringr)
	library(Biostrings)
	library(doMC)
	
	### selects the range for the features
	gs = which(str_detect(s, '^FEATURES'))
	ge = which(str_detect(s,'ORIGIN'))
  
	### find the feature headers
	ss = strsplit(s, ' ')
	ws = which(sapply(ss, function(x)nchar(x[6])>0))
	ws = ws[ws >= gs & ws < ge]
	ws = c(ws, ge-1)
	
	sa = str_replace(s,'^\\s+','')
	
	### loops throught the genomic features and pulls out the data
	lt = list()
	for(i in 1:(length(ws)-1)){
    
		cat(i, '\r')
		lf = list()
		
		### loops through heach feature
		for(j in (ws[i]+1):(ws[i+1]-1)){
			
			line = s[j]
			line = str_replace(line, '^\\s+', '')
			if(str_detect(line, '^/')){
				lf = c(lf, line)
			}else{
				lf[length(lf)] = paste(lf[length(lf)], line, sep='')
			}
			
		}
		lf = sapply(lf, strsplit, '=')
		names(lf) = str_replace(sapply(lf, '[', 1), '/', '')
		lf = lapply(lf, function(x)str_replace_all(x[2],'"',''))
		
		### gets the coordinates as GRanges
		name = ss[[ws[i]]]
		name = name[nchar(name)>0]
		strand = '+'
		if(str_detect(name[2], 'complement')){
			strand = '-'
			name[2] = str_replace(name[2], 'complement', '')
			name[2] = str_replace(name[2], '\\(', '')
			name[2] = str_replace(name[2], '\\)', '')
		}
		if(str_detect(name[2], 'join')){
			name[2] = str_replace(name[2], 'join', '')
			name[2] = str_replace(name[2], '\\(', '')
			name[2] = str_replace(name[2], '\\)', '')
			ns = unlist(strsplit(name[2], ','))
			ns = strsplit(ns, '\\.\\.')
			cn = which(sapply(ns, length) == 1)
			ns[cn] = lapply(ns[cn], function(x)rep(x, 2))
			name = c(name[1], unlist(ns))
		}else{
			name = c(name[1], unlist(strsplit(name[2], '\\W+')))
		}
		ranges = GRanges(chrname, 
						 IRanges(as.numeric(name[seq(2, length(name),2)]), 
								 as.numeric(name[seq(3, length(name),2)])),
						 strand = strand)
		lf$type = name[1]
		lf$ranges = ranges
		### if the id does not exist, creates it from the coordinates
		if(is.null(lf$id)){
		
			ran = range(lf$ranges)
			lf$id = paste(lf$type, unique(seqnames(lf$ranges)), start(ran), end(ran), sep='_')
		}
		
		### fills in the list
		lt[[lf$id]][[lf$type]] = c(lt[[lf$id]][[lf$type]], lf)
	}
	
	### gets the feature types, and constructs the ranges for each feature type
	types = unique(unlist(sapply(lt, function(x)unlist(sapply(x, '[', 'type')))))
	registerDoMC(length(types))
	ltype = list()
	ltype = foreach(i = 1:length(types))%dopar%{
	
		type = types[i]
		print(type)
		tl = lapply(lt, function(x)x[sapply(x, '[', 'type') == type])
		tl = tl[sapply(tl, length) > 0]
		ranges = GRangesList(lapply(tl, function(x)x[[type]]$ranges))
		ranges = unlist(ranges)
		values(ranges)$names = names(ranges)
		
		if(!str_detect(type, 'miRNA')){
			locus = lapply(tl, function(x)x[[type]]$locus_tag)
			gene = lapply(tl, function(x)ifelse(!is.null(x[[type]]$gene), x[[type]]$gene, ''))

			values(ranges)$locus_tag = unlist(locus[names(ranges)])
			values(ranges)$gene      = unlist(gene[names(ranges)])
		
		}else{
			hairpin = lapply(tl, function(x)ifelse(!is.null(x[[type]]$HairpinSequence), x[[type]]$HairpinSequence, ''))
			values(ranges)$hairpin      = DNAStringSet(unlist(hairpin[names(ranges)]))
			
		}
		return(ranges)
	}
	names(ltype) = types
	return(ltype)
}

# ---------------------------------------------------------- #
### extracts the genome sequence
ParseSequence = function(s){
  
  ss = which(str_detect(s,'ORIGIN'))+1
  se = which(str_detect(s,'//'))-1
  seq = str_replace_all(s[ss:se],'\\d','')
  seq = toupper(str_replace_all(seq,'\\s',''))
  return(DNAString(paste(seq, collapse='')))
}
