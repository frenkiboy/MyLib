# ---------------------------------------------------------------------------- #
parse_MEMEdb = function(infile, type='log2probratio', bg=c(A=.25,C=.25,G=.25,G=.25)){
    
    library(stringr)
    library(TFBSTools)
    r = readLines(infile)
    mind = which(str_detect(r, 'MOTIF'))
    
    matlist = list()
    for(i in 1:length(mind)){
        
        cat(i,'\r')
        ind = mind[i]
        name = unlist(strsplit(r[ind],' '))
        mstart = ind + 3
        mend = mstart + (as.numeric(unlist(strsplit(r[ind+2],' '))[6])-1)
        mat = r[mstart:mend]
        mat = str_replace_all(mat,'\\t','')
        mat = do.call(rbind, strsplit(mat,'\\s+'))[,-1]
        mat = apply(mat, 1, as.numeric)
        mat = round(mat*100)
        mat = apply(mat,2,as.integer)
        rownames(mat) = c('A','C','G','T')
        
        
        pfm = PFMatrix(profileMatrix=mat, ID=name[2], name=name[3],
                       matrixClass='transcription factor')
        pwm = toPWM(pfm, type=type, bg=bg)
        matlist[[name[2]]] = pwm
        
    }
    return(matlist)
}
# ---------------------------------------------------------------------------- #
extract_Jaspar_Motifs = function(organism, database='JASPAR2017'){

  library(TFBSTools)

  matlist = getMatrixSet(database, opts=list(species=organism))
  names(matlist) = sapply(matlist, name)
  matlist = lapply(matlist, toPWM)
  return(matlist)
}


# ---------------------------------------------------------------------------- #
scan_Genome = function(matlist, genome.name, outpath, min.score='80%', ncores=16, remove=FALSE){
    
    library(doMC)
    library(TFBSTools)
    library(stringr)
    gname = str_replace(genome.name,'^.+\\.','')
    path_out_scan_genome = file.path(outpath, gname)
    dir.create(path_out_scan_genome, showWarnings=FALSE)
    
    genome = GenomeLoader(genome.name)
    chrs = names(genome)
    chrs = chrs[!str_detect(chrs, 'Un')]
    chrs = chrs[!str_detect(chrs, 'random')]
    chrs = chrs[!str_detect(chrs, 'hap')]
    chrs = setdiff(chrs, c('chrM','chrY'))
    
    gl = lapply(chrs, function(x)genome[[x]])
    names(gl) = chrs
    gl = DNAStringSet(gl)
    
    registerDoMC(ncores)
    donefiles = as.character(list.files(path_out_scan_genome))
    if(remove)
        donefiles = vector()
    
    
    foreach(m = 1:length(matlist), .errorhandling='remove')%dopar%{
        print(m)
        mat = matlist[[m]]
        name = ID(mat)
        outname = paste(name, 'ms',str_replace(min.score,'%',''), 'rds', sep='.')
        if((length(donefiles)==0) || (!outname %in% donefiles)){
            
            print(name)
            hits.p  = searchSeq(mat, gl, strand='+', min.score=min.score)
            hits.m  = searchSeq(mat, gl, strand='-', min.score=min.score)
            ghits.p = try(unlist(GRangesList( lapply(hits, function(x)as(x,'GRanges')))))
            ghits.m = try(unlist(GRangesList( lapply(hits, function(x)as(x,'GRanges')))))
            ghits = GRanges()
            if(!is.na(ghits.p))
                ghits = c(ghits, ghits.p)
            if(!is.na(ghits.m))
                ghits = c(ghits, ghits.m)
            
            if(length(ghits) > 1)
                saveRDS(ghits, file.path(path_out_scan_genome, outname))
        }
    }
}



# ---------------------------------------------------------------------------- #
parse_ScanGenome = function(inpath, regions, ncores=16){

    library(doMC)
    if(length(inpath) == 1){
        scanfiles = list.files(inpath , full.names=TRUE, pattern='rds',recursive=TRUE)
    }else{
        scanfiles = inpath
    }

  registerDoMC(ncores)
  scan.list = list()
  scan.list = foreach(i = 1:length(scanfiles), .inorder=FALSE)%dopar%{

      scanfile = scanfiles[i]
      name = basename(scanfile)
      message(name)

      scan = readRDS(scanfile)
      scan.sel = scan[countOverlaps(scan,regions, ignore.strand=TRUE) > 0]
      message(paste('left:',round(length(scan.sel)/length(scan),3)))
      return(scan.sel)
  }
  name = basename(scanfiles)
  name = str_replace(name,'.ms.+','')
  names(scan.list) = name
  return(scan.list)
}


# ---------------------------------------------------------------------------- #
# Functions for parsing MEME output
