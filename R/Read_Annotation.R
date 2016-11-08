
# --------------------------------------------------------------------------------------------------------- #
# Reads the gtf annotation
ReadGTFAnnotation = function(gtf.path, which.regions='exon', ensembl=FALSE){

    lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
    source(file.path(lib.path, 'ScanLib.R'), local=TRUE)
    source(file.path(lib.path, 'GeneFunctions.R'), local=TRUE)
    require(data.table)
    require(genomation)
    require(GenomicRanges)
    require(stringr)


    if(!file.exists(gtf.path))
        stop('the gtf file does not exist')

    rds.path = str_replace(gtf.path, 'gtf$','rds')
    gtf = vector()
    if(file.exists(rds.path)){
        message('Reading rds file...')
        gtf = readRDS(rds.path)
        
    }else{
        message('Reading gtf file...')
        gtf = gffToGRanges(gtf.path, ensembl=ensembl)
    }
    
    seqlevels(gtf, force=TRUE) = seqlevels(gtf)[!str_detect(seqlevels(gtf),'NT')]
    seqlevels(gtf)[seqlevels(gtf) == 'chrMT'] = 'chrM'
    
    if(which.regions != 'all')
        gtf = gtf[gtf$type %in% which.regions]


    gtf.exon = gtf[gtf$type == 'exon']
    gtf.exon.ge = split(gtf.exon, gtf.exon$gene_id)
    gtf.range.ge = unlist(range(gtf.exon.ge))
    gtf.exon.tr = split(gtf.exon, gtf.exon$transcript_id)


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


    return(list(gtf=gtf,
                gtf.sel = gtf.selected,
                gtf.union=gtf.union,
                annot=gtf.annot))
}
