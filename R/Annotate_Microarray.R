
# ---------------------------------------------------------------------------- #

# cstrand - use it to combine coordinate and strand, to be able to define promoters
Get_Affy_Annotation = function(ids, organism='mmusculus', attribs=NULL, platform, cstrand=FALSE, mart_name = 'ensembl'){

  require(stringr)
  require("biomaRt")

  if(is.null(attribs))
    stop('Please specify the attributes')
  if(is.null(platform))
    stop('Please specify the platform')

  # annotates the results


  message('Fetching biomart annotation ...')
  mart = useMart(mart_name,
                 dataset=paste(organism,"gene_ensembl",sep='_'))

  filters = platform
  bm = getBM(attributes = attribs,
             filter = filters,
             values=ids,
             mart=mart)

  return(data.table(bm))
}



Parse_MogeneGEO = function(platform.id, organism='mmusculus'){

  require(GEOquery)
  require(stringr)
  require(biomaRt)
  require(data.table)

  geo = getGEO(platform.id, destdir='D:/Tmp')
  plat = Table(geo)
  plat = plat[,c('ID','probeset_id','mrna_assignment')]
  plat = data.frame(apply(plat,2,as.character), stringsAsFactors = FALSE)
  nams = lapply(strsplit(plat[,"mrna_assignment"], '\\s+'), function(x)x[str_detect(x,'ENSMUSG')])
  nams = lapply(nams, function(x)if(length(nams)==0){return('')}else{x})
  nams = lapply(nams, function(x)unique(str_replace(unlist(x),'gene:','')))
  nams = lapply(nams, function(x)paste(x, collapse=':'))
  nams = unlist(nams)
  plat[,ncol(plat)] = nams

  message('Fetching biomart annotation ...')
  mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=paste(organism,"gene_ensembl",sep='_'),
                 host="www.ensembl.org")
  attribs = c('affy_mogene_2_1_st_v1','ensembl_gene_id','chromosome_name','start_position','end_position',
              'strand','external_gene_name','gene_biotype')
  filters = 'ensembl_gene_id'
  bm = getBM(attributes = attribs,
             filter = filters,
             values=plat$mrna_assignment,
             mart=mart)
  bmd = data.table(id = as.character(bm[,1]),
                   gene_id = bm[['ensembl_gene_id']],
                   coord=paste(paste0('chr',bm[['chromosome_name']]),
                               paste(bm[['start_position']],
                                     bm[['end_position']],sep='-'),sep=':'),
                   strand=ifelse(bm[['strand']]==1,'+','-'),
                   gene_name = bm[['external_gene_name']],
                   gene_biotype = bm[['gene_biotype']])

  message('Collapsing IDs ...')
  bmd = bmd[,lapply(.SD, function(x)paste(unique(x),collapse=',')),by=id]
  return(bmd)


}


# ---------------------------------------------------------------------------- #
annotate_Refseq = function(ref.ids, organism='mmusculus'){

  library(biomaRt)
  message(organism)
  if(organism == 'mmusculus')
      symbol = 'mgi_symbol'
  if(organism == 'hsapiens')
      symbol = 'hgnc_symbol'
  mart = useMart("ensembl", dataset=paste(organism,"gene_ensembl", sep='_'))
  bmm =  getBM(attributes=c('ensembl_gene_id',symbol,'refseq_mrna'),
              filters='refseq_mrna',  values = ref.ids, mart = mart)
  bmn =  getBM(attributes=c('ensembl_gene_id',symbol,'refseq_ncrna'),
              filters='refseq_ncrna',  values = ref.ids, mart = mart)
  colnames(bmn)[3] = 'refseq_mrna'

  bm = rbind(bmm, bmn)
  colnames(bm) = c('gene_id','gene_name','refseq_id')
  return(bm)
}
