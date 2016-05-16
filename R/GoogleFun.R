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
