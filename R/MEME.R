


#####----------------------------------------------------------------------------------------------#####
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


# ------------------------------------------------------ #
# Gets blocks information from meme output
MemeBlocks = function(path){
	print('Reading in the file...')
	connection = file(path, open='')
	file = scan(connection, what='list', sep='\n')
	close(connection)
	### gets the start and end positions of all motifs
	blocks.start = grep('BLOCKS', file) + 3
	blocks.end = grep('^//', file) - 1

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
	names(l.blocks) = paste('motif', 1:length(l.blocks), sep='.')
	return(l.blocks)
}


# ------------------------------------------------------ #
read_MEME_Blocks = function(path, out='PWM', type='prob'){

  require(TFBSTools)

	if(length(path) == 1){
  	meme.files = list.files(path, recursive=TRUE, pattern='.txt', full.names=TRUE)
	}else{
		meme.files = path
	}

  meme.names = basename(dirname(meme.files))
  blocks = lapply(meme.files, MemeBlocks)
  blocks = lapply(blocks, lapply, BlocksToPFM)
  blocks = lapply(1:length(blocks), function(x){
																		ml=blocks[[x]];
																		names(ml) = paste(meme.names[1],1:length(ml),sep='.')
																		ml
																		})

	blocks = unlist(blocks, recursive=FALSE)


  if(out == 'PWM')
      memes = lapply(blocks, function(x)toPWM(as.matrix(x), type=type))

  return(memes)

}

# ------------------------------------------------------ #
# takes a meme class object and returns a PFM
BlocksToPFM = function(char){

	d = apply(do.call(rbind, strsplit(char, split='')),2, function(x)table(factor(x, levels=c('A','C','G','T'))))
	return(d)
}
