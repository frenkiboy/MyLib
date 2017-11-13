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
#' scan_Genome - given a set of PWMs and a BSgenome object scans the genome
#' and saves the output as a BAM file
#'
#' @param matlist - list of PWMs
#' @param genome.name - name of the used genome
#' @param genome - BSgenome object
#' @param chrs - which chromosomes to use
#' @param outpath - location of output bam files
#' @param min.score - minimal scanning score
#' @param ncores - number of cores
#' @param remove - whether to repeat the scanning procedure
#'
#' @return
#' @export
#'
#' @examples
scan_Genome = function(
  matlist,
  genome.name,
  genome = NULL,
  chrs   = NULL,
  outpath,
  min.score = '80%',
  ncores    = 16,
  remove    = FALSE
){

    library(doMC)
    suppressPackageStartupMessages(library(GenomicAlignments))
    library(TFBSTools)
    library(stringr)
    
    if(is.null(genome))
        stop('Genome object is not assigned')
    
    if(is.null(chrs))
        chrs = names(genome)
    
    source(file.path(lib.path, 'FileLoader.R'), local=TRUE)
    gname = str_replace(genome.name,'^.+\\.','')
    path_out_scan_genome = file.path(outpath, gname)
    dir.create(path_out_scan_genome, showWarnings=FALSE)

 
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
        outname = paste(name, 'ms',str_replace(min.score,'%',''), 'bam', sep='.')
        if((length(donefiles)==0) || (!outname %in% donefiles)){

            print(name)
            hits.p  = searchSeq(mat, gl, strand='+', min.score=min.score)
            hits.m  = searchSeq(mat, gl, strand='-', min.score=min.score)
            ghits.p = try(unlist(GRangesList( lapply(hits.p, function(x)as(x,'GRanges')))))
            ghits.m = try(unlist(GRangesList( lapply(hits.m, function(x)as(x,'GRanges')))))
            ghits = GRanges()
            if(!is.null(ghits.p))
                ghits = c(ghits, ghits.p)
            if(!is.null(ghits.m))
                ghits = c(ghits, ghits.m)

            seqlengths(ghits) = seqlengths(genome)[seqlevels(ghits)]
            genome(ghits) = genome.name


            ga = as(ghits, 'GAlignments')
            colnames(values(ga))[3] = 'AS'
            colnames(values(ga))[4] = 'RS'
            colnames(values(ga))[8] = 'seq'
            rtracklayer::export(ga, con=BamFile(file.path(path_out_scan_genome, outname)))
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
      if(length(scan)>0){
        scan.sel = scan[countOverlaps(scan,regions, ignore.strand=TRUE) > 0]
        message(paste('left:',round(length(scan.sel)/length(scan),3)))
        return(scan.sel)
      }
  }
  name = basename(scanfiles)
  name = str_replace(name,'.ms.+','')
  names(scan.list) = name
  return(scan.list)
}


# ---------------------------------------------------------------------------- #
# Functions for parsing MEME output
