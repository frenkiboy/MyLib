source.lib('Decorate.R')
source.lib('Decorators.R')


# --------------------------------------------------------------------------------------------------------- #
# Reads the gtf annotation
ReadGTFAnnotation = function(gtf.path, which.regions='exon', ensembl=FALSE){

    source(file.path(lib.path, 'ScanLib.R'), local=TRUE)
    source(file.path(lib.path, 'GeneFunctions.R'), local=TRUE)
    require(data.table)
    require(genomation)
    require(GenomicRanges)
    require(stringr)


    if(!file.exists(gtf.path))
        stop('the gtf file does not exist')

    rds.path = str_replace(gtf.path, 'gtf$','rds')

    message('Importing gtf...')
    gtf = RCAS::importGtf(gtf.path)
    gtf.exon = subset(gtf, type == 'exon')
    gtf.exon.ge = split(gtf.exon, gtf.exon$gene_id)
    gtf.range.ge = unlist(range(gtf.exon.ge))
    gtf.exon.tr = split(gtf.exon, gtf.exon$transcript_id)

    gtf.trans = NULL
    if(any(gtf$type == 'transcript'))
        gtf.trans = subset(gtf, type == 'transcript')

    message('Selecting transcripts...')
    gtf.selected = GtfSelectTranscript(gtf.exon)
    gtf.selected = split(gtf.selected, gtf.selected$transcript_id)

    message('Making union...')
    gtf.union = reduce(gtf.exon.ge)

    message('Constructing annotation...')
    gtf.annot = unique(as.data.frame(DataFrame(gtf.exon))[,c('gene_id','transcript_id','gene_name','gene_biotype')])

    gtf.annot$gcoord  = GCoords(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
    gtf.annot$gwidth = width(gtf.range.ge[match(gtf.annot$gene_id, names(gtf.range.ge))])
    gtf.annot$tcoord  = GCoords(unlist(range(gtf.exon.tr))[match(gtf.annot$transcript_id, names(gtf.exon.tr))])
    gtf.annot$twidth = sum(width(gtf.exon.tr))[gtf.annot$transcript_id]


    return(list(gtf       = gtf.exon,
                gtf.sel   = gtf.selected,
                gtf.union = gtf.union,
                gtf.trans = gtf.trans,
                annot     = gtf.annot))
}



# ---------------------------------------------------------------------------- #
# reads the complete gene annotation from the gtf file
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
    }else{
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
# reads two columns from the gtf annotation
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
    }else{
        gtf = path_gtf
    }
    gtf = subset(gtf, type == 'exon')
    d = as.data.frame(values(gtf)) %>%
      dplyr::select_(col1, col2) %>%
      distinct()
    return(d)
}
