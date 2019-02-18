# ---------------------------------------------------------------------------- #
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot)
  library(dplyr)
})

# ---------------------------------------------------------------------------- #
setGeneric("plotExpression", function(object) {
  standardGeneric("plotExpression")
})

setMethod("plotExpression", signature("Seurat", "character","character","character","data.frame"),
  function(
    seu       = NULL,
    id        = NULL,
    id_type   = 'gene_id',
    expr_type = 'scale_data',
    dr_type   = 'tsne',
    annot
){
  if(is.null(seu))
    stop('Please provide a Seurat object')

  if(is.null(gene_id))
    stop('Please provide a gene_id')

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
  if(!any(gene_id %in% rownames(seu@raw.data)))
    return(NULL)

  g1 = seu@dr[[tsne]]$cell.embeddings %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('X1','X2')) %>%
    mutate(expr = slot(seu, expr_type)[gene_id,]) %>%
    ggplot(data = .,aes(X1 , X2, color=expr)) %>%
      geom_point(size=.5) +
      scale_color_gradient2() +
      ggtitle(gene_name) +
      xlab(paste0(dr_type,'1')) +
      ylab(paste0(dr_type,'2'))

  return(g1)

})


# ---------------------------------------------------------------------------- #
setGeneric("plotMetaColumn", function(object) {
  standardGeneric("plotMetaColumn")
})

# method when title not set
setMethod("plotMetaColumn", signature("Seurat", "character",NULL,"character"),
  function(
    seu         = NULL,
    column_name = 'scale_data',
    title       = NULL,
    dr_type     = 'tsne'
){
  plotMetaColumn(
    seu,
    column_name,
    column_name,
    dr_type
    )
}

setMethod("plotMetaColumn", signature("Seurat", "character","character","character"),
  function(
    seu         = NULL,
    column_name = 'scale_data',
    title       = NULL,
    dr_type     = 'tsne'
){
  if(is.null(seu))
    stop('Please provide a Seurat object')

  if(is.null(column_name))
    stop('Please provide a gene_id')

  if(is.null(annot))
    stop('Please provide an annotation file')

  # checks whether the gene id exists
  if(!any(column_name %in% colnames(seu@meta.data)))
    return(NULL)

  g1 = seu@dr[[tsne]]$cell.embeddings %>%
    as.data.frame() %>%
    magrittr::set_colnames(c('X1','X2')) %>%
    mutate(meta = seu@meta.data[[column_name]]) %>%
    ggplot(data = .,aes(X1 , X2, color=meta)) %>%
      geom_point(size=.5) +
      scale_color_brewer() +
      ggtitle(column_name) +
      xlab(paste0(dr_type,'1')) +
      ylab(paste0(dr_type,'2'))

  return(g1)

})
