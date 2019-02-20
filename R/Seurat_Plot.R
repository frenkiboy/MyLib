# ---------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})
source.lib('ScanLib.R')
# ---------------------------------------------------------------------------- #


plotExpression = function(
    object    = NULL,
    id        = NULL,
    id_type   = 'gene_id',
    expr_type = 'scale_data',
    dr_type   = 'tsne',
    annot
){
  if(is.null(object))
    stop('Please provide a Seurat object')

  if(is.null(id))
    stop('Please provide a gene id')

  if(is.null(annot))
    stop('Please provide an annotation file')

  if(!id_type %in% c('gene_name','gene_id'))
    stop('id_type should be gene_name or gene_id')


  if(id_type == 'gene_id'){
      gene_id   = id
      gene_name = subset(annot, gene_id == id)$gene_name
  }
  if(id_type == 'gene_name'){
      gene_name = id
      gene_id   = subset(annot, gene_name == id)$gene_id
  }

  # checks whether the gene id exists
  if(!any(gene_id %in% rownames(object@raw.data)))
    return(NULL)

  g1 = object@dr[[dr_type]]@cell.embeddings %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('X1','X2')) %>%
    mutate(expr = slot(object, expr_type)[gene_id,]) %>%
    ggplot(data = .,aes(X1 , X2, color=expr)) +
      geom_point(size=.5) +
      scale_color_gradient2() +
      ggtitle(gene_name) +
      xlab(paste0(dr_type,'1')) +
      ylab(paste0(dr_type,'2'))

  return(g1)

}


# ---------------------------------------------------------------------------- #
# method when title not set
plotMetaColumn = function(
    seu         = NULL,
    column_name = 'nGene',
    title       = NULL,
    dr_type     = 'tsne',
    size        = .5
){
  if(is.null(seu))
    stop('Please provide a Seurat object')

  if(is.null(column_name))
    stop('Please provide a gene_id')

  # checks whether the gene id exists
  if(!any(column_name %in% colnames(seu@meta.data)))
    stop('Meta column does not exist')

  g1 = seu@dr[[dr_type]]@cell.embeddings %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('X1','X2')) %>%
    mutate(meta = seu@meta.data[[column_name]]) %>%
    ggplot(data = ., aes(X1 , X2, color=meta)) +
      geom_point(size=size) +
      ggtitle(column_name) +
      xlab(paste0(dr_type,'1')) +
      ylab(paste0(dr_type,'2'))

  if(!is.numeric(seu@meta.data[[column_name]])){
    cols = ggplotColors(length(unique(seu@meta.data[[column_name]])))
    g1 = g1 + scale_color_manual(values=cols)

  }else{
    g1 = g1 + scale_color_gradient2()
  }


  return(g1)

}


# ---------------------------------------------------------------------------- #
plotPCS = function(
  seu         = NULL,
  column_name = 'nGene',
  title       = NULL,
  size        = .5
){
  if(is.null(seu))
    stop('Please provide a Seurat object')

  if(is.null(column_name))
    stop('Please provide a gene_id')

  if(!any(column_name %in% colnames(seu@meta.data)))
    stop('Meta column does not exist')

  dat = seu@dr[[dr_type]]@cell.embeddings %>%
    as.data.frame() %>%
    mutate(meta_column = seu@meta.data[[column_name]])
  gl = lapply(2:(ncol(dat) - 1), function(x){
    pc = paste0('PC',x)
    message(pc)
    g =  dat %>%
      ggplot(aes_string('PC1', pc, color='meta_column')) +
        geom_point(size=size) +
        ggtitle(column_name) +
        xlab('PC1') +
        ylab(pc)

      if(!is.numeric(seu@meta.data[[column_name]])){
        cols = ggplotColors(length(unique(seu@meta.data[[column_name]])))
        g = g + scale_color_manual(values=cols)
      }else{
        g = g + scale_color_gradient2()
      }


      return(g)
    })
    return(gl)

}
