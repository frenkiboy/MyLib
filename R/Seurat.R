library(Seurat)

# ---------------------------------------------------------------------------- #
#' Process_Seurat - normalizes expression and calculates PCA + tSNE
#'
#' @param seu - Seurat object
#' @param regress - which variables to regress from expression (i.e. cell cycle)
#' @param tsne - whether to calculate tSNE (slows things down)
#' @param pcs_for_tsne - how many PCs to use for tSNE

Process_Seurat = function(
  seu,
  regress      = NULL,
  tsne         = TRUE,
  pcs_for_tsne = 1:9,
  genes_for_tsne = NULL
  ){
  source(file.path(pro.path, 'Seurat.R'))
  library(Seurat)

  message('Normalize ...')
    seu = NormalizeData(seu)

    if(is.null(regress)){
      message('Scale ...')
      seu = ScaleData(object = seu)
    }

    if(!is.null(regress) && all(is.character(regress)) && all(regress %in% colnames(seu@meta.data))){
      message('Scale Regress...')
      seu = ScaleData(object = seu, vars.to.regress = regress)
    }

  message('Variable genes ...')
    seu = FindVariableGenes(object = seu, do.plot = FALSE)

  message('PCA ...')
    seu = RunPCA(seu, do.print=FALSE)

  if(tsne){
    message('TSNE ...')
    if(is.null(genes_for_tsne)){
      seu = RunTSNE(object = seu, dims.use = pcs_for_tsne, do.fast = TRUE, check_duplicates = FALSE)
    }else{
      seu = RunTSNE(object = seu, genes.use = genes_for_tsne, do.fast = TRUE, check_duplicates = FALSE)
    }
  }
  Check_Seurat(seu)

  return(seu)

}

# ---------------------------------------------------------------------------- #
Check_Seurat = function(
  seu
){
  message = vector()
  if(nrow(seu@meta.data) == 0)
    message = c(message, 'meta.data has zero rows')

  if(!all(rownames(seu@meta.data) == colnames(seu@raw.data)))
    message = c(message, 'Cells do not correspond to raw.data')

  if(!all(rownames(seu@meta.data) == colnames(seu@data)))
    message = c(message, 'Cells do not correspond to data')

  if(!all(rownames(seu@meta.data) == colnames(seu@scale.data)))
    message = c(message, 'Cells do not correspond to scale.data')

  if(length(message) > 0)
    stop(paste(message, collapse='\n'))

}


# ---------------------------------------------------------------------------- #
#' Subset_Seurat - subsets seurat object by cells
#'
#' @param seu a seurat object
#' @param ind a named vector containing the indices of which cells to keep
#'
Subset_Seurat = function(
  seu = NULL,
  ind = TRUE
){
  if(is.null(seu))
    stop('Seurat object is not supplied')

  if(any(FALSE %in% ind))
    ind = ind[ind]

  if(is.null(names(ind)))
    stop('ind must be a named vector')

  if(any(!names(ind) %in% rownames(seu@meta.data)))
    stop('ind contains unkown cell types')

  meta.data = seu@meta.data
  meta.data = meta.data[rownames(meta.data) %in% names(ind),]
  seu@meta.data  = meta.data
  seu@cell.names = rownames(seu@meta.data)

  seu@raw.data = seu@raw.data[,match(rownames(seu@meta.data), colnames(seu@raw.data))]
  seu.sub = CreateSeuratObject(
    raw.data  = seu@raw.data,
    meta.data = seu@meta.data)

  return(seu.sub)
}


# ---------------------------------------------------------------------------- #
#' Subset_Genes_Seurat - subsets seurat object by genes
#'
#' @param seu - a seurat object
#' @param genes - ENSEMBL gene id's of genes which to subset
#' @param what- what/remove - whether to keep the genes or remove the genes
#' @param process - TRUE/FALSE, whether to normalize the data after selection and
#' calculate PCA and tSNE (default TRUE)
#' @param regress - variables to regress during normalization
#' @param pcs_for_tsne - how many PCs to use for tSNE
#'
#' @examples
Subset_Genes_Seurat = function(
  seu          = NULL,
  genes        = NULL,
  what         = 'remove',
  process      = TRUE,
  regress      = NULL,
  pcs_for_tsne = 1:9
){
  if(is.null(seu))
    stop('Seurat object is not supplied')

  if(is.null(genes))
    stop('genes are not supplied')

  if(!is.character(genes))
    stop('genes have to be a character vector')

  if(any(!genes %in% rownames(seu@raw.data)))
    stop('unknown genes are supplied')

  if(what == 'remove'){
    genes_ind = !rownames(seu@raw.data) %in% genes
  }else{
    genes_ind = rownames(seu@raw.data) %in% genes
  }
  seu@raw.data   = seu@raw.data[genes_ind, ]
  seu@data       = seu@data[genes_ind, ]
  seu@scale.data = seu@scale.data[genes_ind, ]


  if(process){
    seu = Process_Seurat(seu, regress = regress, pcs_for_tsne=pcs_for_tsne)
  }else{
    warning('The Seurat object is not processed - matrices contain wrongly normalized expression values')
  }
  return(seu)
}
