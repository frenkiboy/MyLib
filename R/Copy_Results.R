# ---------------------------------------------------------------------------- #
file_Copy = function(infile, outpath, dirname){

    outdir=file.path(outpath, dirname)
        dir.create(outdir)

    file.copy(infile, outdir, overwrite=TRUE)
}

copy_Public = function(infile, dirname){

    file_Copy(infile, path_public, dirname)
}


copy_Private = function(infile, dirname){

    file_Copy(infile, path_private, dirname)
}

# ---------------------------------------------------------------------------- #

copy_Public_DT = function(dt, outpath, name, dirname=NULL){

  if(is.null(dirname))
    stop('Please specify the dirname')

  message('Writing...')
  infile = file.path(outpath, DateNamer(paste(name,'txt', sep='.')))

  message('Copying...')
  write.table(dt, infile, row.names=F, col.names=T, quote=F, sep='\t')

  if(!file.exists(infile))
    stop('infile does not exist')

  copy_Public(infile, dirname)

}


# ---------------------------------------------------------------------------- #

copy_Public_rds = function(l, outpath, name, dirname=name){
    
    if(is.null(dirname))
        stop('Please specify the dirname')
    message('Writing...')
    infile = file.path(outpath, DateNamer(paste(name,'rds', sep='.')))
    saveRDS(l, infile)
    
    if(!file.exists(infile))
        stop('infile does not exist')
    
    message('Copying...')
    copy_Public(infile, dirname)
    
}

