# ---------------------------------------------------------------------------- #
read_CellCycle = cacheFile(path_RDS) %@% function(
    path_cc  = NULL,
    path_gtf = NULL
){

    if(is.null(path_cc) || !file.exists(path_cc))
        stop('Cell cycle genes are not specified')

    if(is.null(path_gtf) || !file.exists(path_gtf))
        stop('GTF genes are not specified')

    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(dplyr)
    })
    gtf = RCAS::importGtf(path_gtf, keepStandardChr=FALSE)
    annot = unique(as.data.frame(values(gtf)[,c('gene_id','gene_name')])) %>%
        subset(!is.na(gene_id))%>%
        subset(!is.na(gene_name))

    cc = read.table(path_cc) %>%
        mutate(cycle = rep(c('S','G2M'), times=c(43,55))) %>%
        `colnames<-`(c('gene_name','cycle')) %>%
        merge(annot, by='gene_name') %>%
        arrange(cycle)

    return(cc)
}



