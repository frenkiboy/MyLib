# ---------------------------------------------------------------------------- #
summarize_tomtom = function(tomtom, n=5){
    
    library(data.table)
    tom = rbindlist(tomtom)
    ltom = lapply(split(tom, tom$motname), function(x){
        d = x[order(x$evalue),]
        d = subset(d, evalue < 0.05)
        d = d[1:n,]
        if(nrow(d) == 0)
            return(data.table(motname=unique(tom$motname)), gene_name='')
        
        d[,list(gene_name = paste(unique(gene_name), collapse='|')), by=motname]
    })
    dtom = na.omit(rbindlist(ltom))
    return(dtom)
}


# ---------------------------------------------------------------------------- #
annotate_tomtom = function(path, base=NULL, outpath, dist='pearson'){
    
    if(is.null(base))
        stop('database is not specified')
    
    if(length(path) == 1){
        motfiles = list.files(path, recursive=TRUE, pattern='.txt', full.names=TRUE)
    }else{
        motfiles = path
    }
    
   message('tomtom...')
   path.tomtom = file.path(outpath, 'Tomtom')
        dir.create(path.tomtom)
    
    lmot = list()
    for(i in 1:length(motfiles)){
        motfile = motfiles[i]
        motname = basename(dirname(motfile))
        tomtom = run_tomtom(motfile, base, path.tomtom, dist=dist)
        if(!is.null(tomtom))
            tomtom = data.frame(sample=motname, tomtom)
       
       tomtom$motname = with(tomtom, paste(sample, paste('M', sprintf('%05d',qnum),sep=''), sep='.'))
       lmot[[motname]] = tomtom
    }
    return(lmot)
    
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
check_DNAStringSet_Name = function(set){

  if(any(names(set) == NULL))
    names(set) = paste('s', 1:length(set), sep='')

  return(set)
}
