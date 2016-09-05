


# ---------------------------------------------------------------------------- #
### Parses meme output
MemeParser = function(path, what='pwm'){

	### reading and parsing the meme output file
	print('Reading in the file...')
	connection = file(path, open='')
	file = scan(connection, what='list', sep='\n')
	close(connection)
	### gets the start and end positions of all motifs
	motif.ind.start = grep('position-specific probability matrix', file) + 3
	motif.ind.end = grep('regular expression', file) - 3

	### list for holding the pwms
	l.mat = vector(mode='list', length=length(motif.ind.start))
	print('Parsing the motifs...')
	for (i in seq(along=motif.ind.start)){
		mat = file[motif.ind.start[i]:motif.ind.end[i]]
		mat = do.call(rbind, strsplit(mat ,split='\\s+', perl=T))[,-1]
		mat = data.frame(apply(mat, 1, as.numeric))
		rownames(mat) = c('A','C','G','T')
		colnames(mat) = 1:ncol(mat)
		attrs = strsplit(file[motif.ind.start[i]-2],' ')[[1]]
		attr(mat,'nsites') = as.numeric(attrs[8])
		attr(mat,'Eval') = as.numeric(attrs[10])
		l.mat[[i]] = mat
	}
	print('Returning the pwm list...')
	return(l.mat)
}


# ---------------------------------------------------------------------------- #
### Parses meme output
read_MEME = function(path){
    
    motif.paths = list.files(path, full.names=TRUE, recursive=TRUE, pattern='txt')
    memel = lapply(motif.paths, MemeParser)
    names(memel) = basename(dirname(motif.paths))
    return(memel)
}

# ---------------------------------------------------------------------------- #
# Gets blocks information from meme output
MemeBlocks = function(path, sampname=NULL){
	print('Reading in the file...')
	connection = file(path, open='')
	file = scan(connection, what='list', sep='\n')
	close(connection)
	### gets the start and end positions of all motifs
	blocks.start = grep('BLOCKS', file) + 3
	blocks.end   = grep('^//', file) - 1

	if(length(blocks.start) == 0)
	    return(list(blocks=NULL, header=NULL))
	
	vals = apply(do.call(rbind, lapply(strsplit(file[grepl('E.value', file)][-1], split='\\s+'), '[', c(9,15))),2,as.numeric)

	l.blocks=list()
	for(i in seq(along=blocks.start)){

		block = file[blocks.start[i]:blocks.end[i]]
		l.blocks[[i]] = structure(unlist(lapply(strsplit(block, split='\\s+'), '[', 4)),
								  e.val = vals[i,2],
                                  nsites = vals[i,1],
								  class = 'meme',
								  pfm=matrix(),
								  pwm=matrix(),
								  prior=vector())
	}
	
	header.ind = grep('^MOTIF', file)
	header = data.frame(do.call(rbind, strsplit(file[header.ind],'\\s+'))[,c(2,6,9,15)])
	colnames(header) = c('name','width','hits','Eval')
	options(scipen=10)
	header[,c(1,3:4)] = apply(header[,c(1,3:4)], 2, as.numeric)
	header$name = paste('M', sprintf('%05d', header$name),sep='')
	names(l.blocks) = header$name
	if(is.null(sampname))
	    header = data.frame(sampname = basename(dirname(path)),header)
	
	header = subset(header, Eval < 0.05)
	if(nrow(header) == 0)
        return(list(blocks=NULL, header=NULL))
	
	
	return(list(blocks = l.blocks[header$name], header=header))
}


# ---------------------------------------------------------------------------- #
read_MEME_Output = function(path, out='PWM', type='prob'){

  require(TFBSTools)

	if(length(path) == 1){
   	    meme.files = list.files(path, recursive=TRUE, pattern='.txt', full.names=TRUE)
	}else{
		meme.files = path
	}
    if(length(meme.files) == 0)
        stop('There are no imput.files')

  blocks = lapply(meme.files, MemeBlocks)
  pwml = list()
  for(i in 1:length(blocks)){
      
      bl = blocks[[i]]
      if(is.null(bl$header))
          next()
      bl$pfm = lapply(bl$blocks, BlocksToPFM)
      sampname = unique(bl$header$sampname)
      message(sampname)
      pwms = lapply(bl$header$name, function(x)
          toPWM(PFMatrix(profileMatrix=as.matrix(bl$pfm[[x]]),
                   ID  = paste(sampname,x,sep='.'), 
                   name= paste(sampname,x,sep='.'), 
                   strand='+'), type='prob'))
      pwml[[sampname]] = pwms
  }
  pwmu = unlist(pwml)

  return(pwmu)

}

# ---------------------------------------------------------------------------- #
find_Motif_Hits = function(motifs, seqs, min.score='85%', ncores=24, strand='*'){
    
    require(doMC)
    registerDoMC(ncores)
    lmot = list()
    lmot = foreach(i = 1:length(motifs))%dopar%{
        motname = names(motifs)[i]
        message(motname)
        hits = searchSeq(motifs[[motname]], seqs, strand='*', min.score=min.score)
        hits = data.table(as.data.frame(hits))
        return(hits)
    }
    dmot = rbindlist(lmot)
    return(dmot)
}

# ---------------------------------------------------------------------------- #
Summarize_Motif_Hits = function(hits, outpath, name){
    
}

    
# ---------------------------------------------------------------------------- #
# takes a meme class object and returns a PFM
BlocksToPFM = function(char){

	d = apply(do.call(rbind, strsplit(char, split='')),2, function(x)table(factor(x, levels=c('A','C','G','T'))))
	return(d)
}
