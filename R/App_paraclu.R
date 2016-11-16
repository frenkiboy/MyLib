# ---------------------------------------------------------------------------- #
# runs paraclu for read clustering

run_Paraclu = function(infile, outfile, min.value){
    
    outfile_paraclu = paste(outfile,'mv',min.value,'txt', sep='.')
    command = paste(path_paraclu,min.value,infile,'>',outfile_paraclu)
    system(command, wait=TRUE)
    
    outfile_paracut = paste(outfile,'mv',min.value,'cut.txt', sep='.')
    command = paste(paste(path_paraclu,'-cut.sh',sep=''), outfile_paraclu,'>',outfile_paracut)
    system(command, wait=TRUE)
    
    outfile_bed = paste(outfile,'mv',min.value,'bed', sep='.')
    command = paste('cut -f1,3,3', outfile_paracut, '>',outfile_bed)
    system(command, wait=TRUE)
    
}


# ---------------------------------------------------------------------------- #
Paraclu = function(bamfile, outpath, min.value=10, ncores=12){
    
    library(GenomicAlignments)
    library(GenomicRanges)
    library(doMC)
    source(file.path(lib.path, 'BamWorkers.R'))
    
    message('Reading file...')
    param = ScanBamParam(flag=scanBamFlag(isSecondaryAlignment=FALSE))
    ga = readGAlignments(bamfile, 
                         param=param)
    gl = grglist(ga)
    gl = gl[elementNROWS(gl) == 1]
    gl = resize(gl, width=1, fix='end')
    gl = unlist(gl)
    gc = lapply(c('+','-'), function(x){
        g = as(coverage(gl[strand(gl)==x]),'GRanges')
        strand(g) = x
        g
    })
    
    message('Export infile...')
    bamname = BamName(bamfile)
    da = lapply(gc, function(x)as.data.frame(subset(x, score>0))[,c(1,5,2,6)])
    da = do.call(rbind,da)
    infile = file.path(outpath, paste(bamname, 'tab', sep='.'))
    write.table(da, infile, col.names=F, row.names=F, quote=F, sep='\t')
    
    outfile = file.path(outpath, bamname)
    
    registerDoMC(ncores)
    message('Run paraclu...')
    foreach(i = 1:length(min.value))%dopar%{
        message(i)
        run_Paraclu(infile, outfile, min.value[i])
    }
    message('Done!')
    
}

# ---------------------------------------------------------------------------- #
