# ---------------------------------------------------------------------------- #
run_dreme = function(fseq, bseq, outpath, outname, nmotifs=10, strand=TRUE, base=NULL){

  require(Biostrings)
  require(XML)
  path.out.set = file.path(outpath, outname)
    dir.create(path.out.set)
  if(!file.exists(path.out.set))
    stop('Unable to create the output directory')

  ffile = file.path(path.out.set, paste(outname, 'fwd', 'fa', sep='.'))
  writeXStringSet(check_DNAStringSet_Name(fseq), ffile)

  bfile = file.path(path.out.set, paste(outname, 'bck', 'fa', sep='.'))
  writeXStringSet(check_DNAStringSet_Name(bseq), bfile)

  message('dreme...')
  command = paste('dreme',
            '-oc', path.out.set,
            '-p', ffile,
            '-n', bfile,
            '-m', nmotifs,
            '-png')

  if(strand)
   command = paste(command, '-norc')
  system(command, wait=TRUE)

  if(!is.null(base)){
    message('tomtom...')
    path.out.set.tomtom = file.path(path.out.set, 'Tomtom')
      dir.create(path.out.set.tomtom)

      motfile = list.files(path.out.set, pattern='txt', full.names=TRUE)
      tomtom = run_tomtom(motfile, base, path.out.set.tomtom)
      if(!is.null(tomtom))
        tomtom = data.frame(sample=outname, tomtom)
  }else{
      tomtom=NULL
  }

      message('parsing...')
      res = parse_dreme(path.out.set)
      res$tomtom = tomtom

  return(res)
}


# ---------------------------------------------------------------------------- #
run_tomtom = function(motfile, base, outpath, dist='pearson'){

  require(XML)
  if(is.null(base))
    stop('motif database is not specified')
  mots = scan(motfile, what='character', sep='\n',quiet=TRUE)
  if(length(mots) == 14)
    return(NULL)

  command = paste('tomtom',
            '-oc', outpath,
            '-dist', dist,
            '-min-overlap', 2,
            '-png',
            '-evalue',
            '-thresh',2,
            motfile, base)
  system(command, wait=TRUE)
  infile = dir(outpath, pattern='xml', full.names=TRUE)
  data = xmlParse(infile)
  data = xmlToList(data)
  hits = lapply(head(data$queries[[1]],-1), function(x)do.call(rbind, x[-1]))
  quer = lapply(head(data$queries[[1]],-1), function(x)x[1]$motif$.attrs[1])


  hits.ind = sapply(hits, is.null)
  if(all(unlist(hits.ind))){

    hits = data.frame(qid=unlist(quer), gene_name=NA)

    return(hits)
  }
  hits = hits[!hits.ind]
  quer = quer[!hits.ind]

  hits = suppressWarnings(data.frame(
    qnum = rep(1:length(quer), times=sapply(hits, nrow)),
    qid  = rep(unlist(quer),   times=sapply(hits, nrow)),
    do.call(rbind, hits)))
  attrs = lapply(head(data[[2]]$target_file,-1), function(x)x$.attrs)
  attrs = suppressWarnings(data.frame(do.call(rbind, attrs)))
  hits = merge(attrs[,c('id','alt')], hits, by.x='id', by.y='target')
  hits$id=NULL
  hits[,1:3] = hits[,c(2,3,1)]
  colnames(hits)[1:3] = c('qnum','qid','gene_name')
  hits[,-c(1:4)] = apply(hits[,-c(1:4)], 2, as.numeric)
  hits$qid = sub('^q_','',hits$qid)
  hits = hits[order(hits$qnum),]
  return(hits)
}

# ---------------------------------------------------------------------------- #
parse_dreme = function(inpath){


  require(XML)
  infile = dir(inpath, pattern='xml', full.names=TRUE)
  if(length(infile) == 0)
    return(NULL)
  data = xmlParse(infile)
  data = xmlToList(data)
  if(is.null(data$motifs))
    return(NULL)
  lmot = list()
  motifs = data$motifs
  for(i in 1:length(motifs)){

    message(i)
    mot = motifs[[i]]
    motif = t(do.call(rbind, mot[names(mot) == 'pos']))
    motif = apply(motif[-1,], 2, as.numeric)
    colnames(motif) = paste0('p', 1:ncol(motif))
    rownames(motif) = c('A','C','G','T')

    d = as.data.frame(t(as.data.frame(mot$.attrs)))
    lmot[[d$id]] = list(mat=motif, d = d)

  }
  dall = do.call(rbind, lapply(lmot, '[[', 'd'))
  dall[,3:ncol(dall)] = apply(dall[,3:ncol(dall)], 2, as.numeric)
  lmat = lapply(lmot,'[[','mat')
  return(list(lmat=lmat, motifs=dall))
}

# ---------------------------------------------------------------------------- #
check_DNAStringSet_Name = function(set){

  if(any(names(set) == NULL))
    names(set) = paste('s', 1:length(set), sep='')

  return(set)
}
