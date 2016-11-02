# --------------------------------------------------------------------------- # 
readRDS_Integration = function(date, name){
    
    infile = file.path(integration.path,name, paste(date, name, 'rds', sep='.'))
    if(!file.exists(infile))
        stop('The input file does not exist')
    return(readRDS(infile))
    
}
