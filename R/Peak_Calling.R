# ---------------------------------------------------------------------------- #
call_MACS2 = function(chipfile=NULL, contfile=NULL, outpath, name, genome='hs'){
    
    outdir = file.path(outpath, name)
        dir.create(outdir, showWarnings=FALSE)
        
    if(is.null(chipfile))
        stop('Please specify the ChIP file')
    
    command = paste('macs2 callpeak',
                    '-t', chipfile,
                    '-g', genome,
                    '-n', name,
                    '-f', 'BAM',
                    '--call-summits',
                    '--outdir', outdir)
    if(!is.null(contfile))
        command = paste(command, '-t', contfile)
    system(command)
}

# ---------------------------------------------------------------------------- #

