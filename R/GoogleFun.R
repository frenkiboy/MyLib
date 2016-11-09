# --------------------------------------------------------------- #
### fetches annotation from google tables 
fetchAnnotation = function(key, filter=TRUE, sheet=1){
    
    library(googlesheets)
    library(stringr)
    gap = gs_key(key)
    samps = gs_read(gap, sheet)
    
    if(filter & any(colnames(samps) == 'Approved'))
        samps = subset(samps, Approved == 'Yes')
    
    return(samps)
}


# ---------------------------------------------------------------------------- #
fetchGoogleSheet = function(key, prefix=NULL){
    
    library(googlesheets)

    
    if(is.null(prefix))
        prefix='~/'
    
    token_path = file.path(prefix,'/googlesheets_token.rds')
    if(!file.exists(token_path)){
        token <- gs_auth(cache = FALSE, key=key)
        saveRDS(token, file = token_path)
    }
    suppressMessages(gs_auth(token = token_path, verbose = FALSE))
    gap = gs_key(key)
    return(gap)
}

# ---------------------------------------------------------------------------- #
# fetches annotation from google tables 
fetch_SampleSheet = function(){
    
    library(googlesheets)
    require(stringr)
    require(readr)
    suppressMessages(require(dplyr))
    
    tab = gs_read(fetchGoogleSheet(sample_sheet_key), 'Annotation')
    return(tab)
}