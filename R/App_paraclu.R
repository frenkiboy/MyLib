# ---------------------------------------------------------------------------- #
# runs paraclu for read clustering

run_Paraclu = function(infile, outfile, min.value){
    
    outfile_paraclu = paste(outfile,'mv',min.value,'txt', sep='.')
    command = paste(path_paraclu,min.value,infile,'>'outfile_paraclu)
    system(command, wait=TRUE)
    
    outfile_paracut = paste(outfile,'mv',min.value,'cut.txt', sep='.')
    command = paste(paste(path_paraclu,'sh',sep='.'), outfile_paraclu,'>',outfile_paracut)
    system(command, wait=TRUE)
    
}


# ---------------------------------------------------------------------------- #
Paraclu = function(bamfile, outpath,, min.value=10, ncores=12){
    
    library(GenomicAlignments)
    library(GenomicRanges)
    library(doMC)
    source(file.path(lib.path, 'BamWorkers.R'))
    
    message('Reading file...')
    param = ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE))
    ga = readGAlignments(bamfile, 
                         param=param)
    ga = resize(ga, width=1, fix='end')
    gl = lapply(c('+','-'), function(x)as(coverage(ga[strand(ga)==x])),'GRanges')
    
    message('Export infile...')
    bamname = BamName(bamfile)
    da = lapply(gl, function(x)as.data.frame(gl[[x]])[,c(1,5,2,6)])
    da = do.call(rbind(da))
    infile = file.path(outpath, paste(bamname, 'tab', sep='.'))
    
    outfile = file.path(outpath, bamname)
    
    registerDoMC(ncores)
    message('Run paraclu...')
    foreach(i = 1:length(min.value))%dopar%{
        message(i)
        run_Paraclu(infile, outfile, min.value)
    }
    message('Done!')
    
}

# ---------------------------------------------------------------------------- #
