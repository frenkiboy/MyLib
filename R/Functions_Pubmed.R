# ---------------------------------------------------------------------------- #
#' fetch_Abstract
#'
#' @param pids a vector of Pubmed ids
#'
#' @return a data frame with publication abstracts
fetch_Abstract = function(pids = NULL){

    source(file.path(lib.path, 'ScanLib.R'), local=TRUE)
    library(data.table)
    library(stringr)
    library(openxlsx)
    library(RISmed)
    library(dplyr)
    library(doMC)
    
    if(is.null(pids))
        stop('Please specify the pubmed ids')
    
    if(any(is.na(as.numeric(pids))) & all(nchar(pids) != 8))
       stop('PIDs are not valid')
        
    
    lpub = list()
    for(pid in pids){
        
        message(pid)
        search_topic = pid
        search_query = EUtilsSummary(search_topic, retmax=1000)
        
        records = EUtilsGet(search_query)
        pubmed_data = data.frame(pid = PMID(records),
                                 year = YearAccepted(records),
                                 'Title'=ArticleTitle(records),
                                 'Abstract'=AbstractText(records))
        lpub[[pid]] = pubmed_data
        
    }
    dpub = do.call(rbind, lpub)

    dpub = dpub %>%
        arrange(desc(year)) %>%
        filter(!duplicated(pid))
    
    return(dpub)
}

# ---------------------------------------------------------------------------- #
WriteAbstract = function(dpub, outname, path.out){
    
    outfile = file.path(path.out, outname)
    if(file.exists(outfile))
        file.remove(outfile)
    
    for(i in 1:nrow(dpub)){
        
        paper = dpub[i,]
        header=paste(colnames(paper)[1:2], paper[1:2], sep=':', collapse='\t')
        title = paper[1,3]
        abst  = paper[1,4]
        cat(header, "\n", file=outfile, append=TRUE)
        cat(title , "\n\n", file=outfile, append=TRUE)
        cat(abst  , "\n\n", file=outfile, append=TRUE)
        cat("\n\n", file=outfile, append=TRUE)
    }
}

# ---------------------------------------------------------------------------- #
.test_Functions_Pubmed.R = function(){
    
    pids = c('23872065','23824327','22100165')
    fetch_Abstract(pids)
    
}