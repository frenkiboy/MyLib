source.lib('Decorate.R')
source.lib('Decorators.R')

# ---------------------------------------------------------------------------- #
read_Gene_Annotation = cacheFile(path_RDS) %@% function(
    path_gtf = NULL
){


  if(is.null(path_gtf) || !file.exists(path_gtf))
      stop('GTF genes are not specified')

    suppressPackageStartupMessages({
        library(rtracklayer)
        library(GenomicRanges)
        library(dplyr)
    })

    if(class(path_gtf) == 'character'){
        gtf = import.gff(path_gtf)
    else{
        gtf = path_gtf
    }
    gtf = subset(gtf, type == 'exon')

    gl = range(split(gtf, gtf$gene_id))

    g = unique(as.data.frame(values(gtf))) %>%
        dplyr::select(gene_id, gene_name, gene_biotype) %>%
	distinct() %>%
        mutate(id = as.character(unlist(gl[gene_id])))
    return(g)
}


# ---------------------------------------------------------------------------- #
read_Annotation = function(
  path_gtf = NULL,
  col1 = 'gene_id',
  col2 = 'gene_name'
){
  if(is.null(path_gtf) || !file.exists(path_gtf))
      stop('GTF genes are not specified')

    suppressPackageStartupMessages({
        library(rtracklayer)
        library(GenomicRanges)
        library(dplyr)
    })

    if(class(path_gtf) == 'character'){
        gtf = import.gff(path_gtf)
    else{
        gtf = path_gtf
    }
    gtf = subset(gtf, type == 'exon')
    d = as.data.frame(values(gtf)) %>%
      dplyr::select_(col1, col2) %>%
      distinct()
    return(d)
}
