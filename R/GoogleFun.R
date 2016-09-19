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
fetchGoogleSheet = function(key){
    
    library(googlesheets)

    token_path = '~/googlesheets_token.rds'
    if(!file.exists(token_path)){
        token <- gs_auth(cache = FALSE, key=key)
        saveRDS(token, file = token_path)
    }
    suppressMessages(gs_auth(token = token_path, verbose = FALSE))
    gap = gs_key(key)
    return(gap)
}
