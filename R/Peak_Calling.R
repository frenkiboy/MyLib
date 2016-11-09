callPeaks_MACS2 = function(chipfile=NULL, contfile=NULL, outpath, name, genome='hs'){
    
    if(is.null(chipfile))
        stop('Please specify the ChIP file')
    
    command = paste('MACS2 callpeaks',
                    '-t', chipfile,
                    '-g', genome,
                    -'n', name,
                    '-f', 'BAM',
                    '--call-summits',
                    '--outdir', outpath)
    if(!is.null(contfile))
        command = paste(command, '-t' contfile)
    system(command)
}