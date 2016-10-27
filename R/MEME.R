# ---------------------------------------------------------------------------- #
# given a list of DNSStringSets run MEME
run_MEME = function(lseq, path.out.meme,  minw=8, maxw=15,nmotifs=5, revcomp=FALSE){
    
    suppressPackageStartupMessages(require(TFBSTools))
    suppressPackageStartupMessages(require(Biostrings))
    path.meme = '~/bin/Software/MotifDiscovery/meme/bin/meme'
    options(scipen = 999)

    message('Running...')
   for(i in 1:length(lseq)){
            
        setname = names(lseq)[i]
        print(setname)
        path.out.meme.samp = file.path(path.out.meme, setname)
        dir.create(path.out.meme.samp, showWarnings=FALSE, recursive=TRUE)
        
        seq.set = lseq[[i]]
        seq.set = seq.set[width(seq.set) > 10]
        path.out.meme.samp.seq = file.path(path.out.meme.samp, paste(setname, 'fa', sep='.'))
        writeXStringSet(seq.set, path.out.meme.samp.seq)
            
        arglist = list('-oc'= path.out.meme.samp,
                       '-nmotifs'= nmotifs,
                       '-mod' = 'zoops',
                       '-minw' = minw,
                       '-maxw' = maxw,
                       '-p' = 32,
                       '-maxsize' = sprintf('%s',10000000))
        if(revcomp)
            arglist[['-revcomp']] = ''
        
                                      
        command = paste(path.meme,
                        path.out.meme.samp.seq,
                        '-dna',
                        paste(names(arglist), unlist(arglist), collapse=' '))
            
        system(command, wait=TRUE)
   }
    
    message('Done...')
}

# ---------------------------------------------------------------------------- #
get_MEME_Output = function(lseq, path.out.meme, pattern='',motif.base=NULL){

    library(dplyr)
    source(file.path(lib.path, 'ScanLib.R'))
    source(file.path(lib.path, 'MEME.R'))
    source(file.path(lib.path, 'App_tomtom.R'))
    if(is.null(motif.base))
        motif.base='/home/vfranke/bin/Software/MotifDiscovery/meme/db/motif_databases/JASPAR/JASPAR_CORE_2014_vertebrates.meme'

    indirs = list.files(path.out.meme, full.names=TRUE)
    indirs = indirs[!str_detect(indirs, 'Tomtom')]
    lhits = list()
    for(i in 1:length(indirs)){
        
        sname = names(lseq)[i]
        message(sname)
        seq.set = lseq[[sname]]
  
        memel  = read_MEME_Output(indirs[i])
        tomtom = annotate_tomtom(indirs[i], motif.base, path.out.meme)
        tom.an = summarize_tomtom(tomtom)
        tom.an = na.omit(tom.an)
    
        consl = lapply(memel, function(x)paste(c('A','C','G','T')[apply(Matrix(x),2,which.max)], collapse=''))
        names(consl) = unlist(lapply(memel, function(x)x@ID))
        hits  = find_Motif_Hits(memel, seq.set, strand='+', min.score='80%')
        hits = hits %>% rename(transcript_id = seqnames)
        gsets = data.frame(variable=rep(names(lseq), times=sapply(lseq,length)), 
                           transcript_id=unlist(lapply(lseq, names)),
                           width=unlist(lapply(lseq, width)))
        hits=merge(hits, gsets, by='transcript_id',allow.cartesian=TRUE)
    

        hits.stat = hits[,list(cnts=length(unique(transcript_id))), by=c('TF', 'variable')]
        hits.stat = merge(hits.stat, data.table(gsets)[,.N,by='variable'], by='variable')
        hits.stat = merge(hits.stat, tom.an, by.x='TF', by.y='motname', all.x=TRUE)
        
        hits.stat$freq = with(hits.stat, round(cnts/N,3))
        hits.stat$freq.pseud = with(hits.stat, (cnts+1)/N)
        hits.stat.cnts = dcast(hits.stat, TF~variable, value.var='cnts', fill=0)
        hits.stat.freq = dcast(hits.stat, TF~variable, value.var='freq', fill=0)
        hits.stat.freq.pseud = dcast(hits.stat, TF~variable, value.var='freq.pseud', fill=0)
        lhits[[sname]] = list(hits = hits,
                              hits.stat.freq = hits.stat.freq,
                              hits.stat.cnts = hits.stat.cnts,
                              hits.stat.freq.pseud = hits.stat.freq.pseud,
                              hits=hits, hits.stat=hits.stat,
                              motif.annotation = tom.an,
                              consl = consl)
    }
    return(lhits)        
}



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
