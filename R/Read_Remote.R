# ---------------------------------------------------------------------- #
read_xlsx_remote = function(path, sheets=NA){
    
    require(readxl)
    
    destfile = basename(path)
    download.file(path, destfile = destfile, mode="wb")
    
    
    if(is.na(sheets))
        sheets <- readxl::excel_sheets(destfile)
    
    x = lapply(sheets, function(X) readxl::read_excel(destfile, sheet = X))
    names(x) = sheets
    x
    
}

# ---------------------------------------------------------------------- #
read_rds_remote = function(inpath){
    
    destfile = basename(inpath)
    download.file(inpath, destfile = destfile, mode="wb")
    r = readRDS(destfile)
    return(r)
    
}

# ---------------------------------------------------------------------- #
read_rds_remote_file = function(date, dirname, file.name, path.remote, new=FALSE){
    
    
    if(new){
        read_rds_remote(file.path(path.remote, dirname, paste(date, file.name, 'rds', sep='.')))
    }else{
        message('Loading cached version...')
        outfile=file.path('./',paste(date, file.name, 'rds', sep='.'))
        readRDS(outfile)
    }
    
}
