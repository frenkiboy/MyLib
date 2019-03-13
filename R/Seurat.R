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
  regress        = NULL,
  tsne           = TRUE,
  pcs_for_tsne   = 1:9,
  genes_for_tsne = NULL,
  scnorm         = FALSE,
  scnorm_cond    = NULL,
  scnorm_ncores  = 8
  ){
  source(file.path(lib.path, 'Seurat.R'))
  library(Seurat)

  # -------------------------------------------------------------------------- #
  if(class(rownames(seu@raw.data)) == 'array'){
    rownames(seu@raw.data) = as.character(rownames(seu@raw.data))
  }

  if(class(rownames(seu@data)) == 'array'){
    rownames(seu@data) = as.character(rownames(seu@data))
  }

  # -------------------------------------------------------------------------- #
  message('Normalize ...')
  if(scnorm == FALSE){
    seu = NormalizeData(seu)
  }else{
    library(SCnorm)
    if(is.null(scnorm_cond))
      scnorm_cond = rep(1, nrow(seu@meta.data))

    mat = as.matrix(seu@raw.data)
    colnames(mat) = as.character(1:ncol(mat))
    names(scnorm_cond) = colnames(mat)
    sce = SCnorm(Data = mat, Conditions=scnorm_cond, NCores=scnorm_ncores)
    ass = assays(a)$normcounts
    colnames(ass) = colnames(seu@raw.data)
    seu@data = ass
  }

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
  seu  = NULL,
  ind  = TRUE,
  gind = NULL
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

  if(!is.null(gind) && all(gind %in% rownames(seu@raw.data)))
    seu@raw.data = seu@raw.data[gind,]

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


# ---------------------------------------------------------------------------- #
# ported from Seurat - they have mistake in if else = instead of <=
Cell_Cycle_Scoring = function(object, g2m.genes, s.genes, set.ident = FALSE)
{
  enrich.name <- "CellCycle"
  genes.list <- list(S.Score = s.genes, G2M.Score = g2m.genes)
  object.cc <- AddModuleScore(object = object, genes.list = genes.list,
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list,
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "S", second = "G2M", null = "G1") {
    if (all(scores <= 0)) {
      return(null)
    }
    else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score",
                               "Phase")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
  object <- AddMetaData(object = object, metadata = cc.scores)
  if (set.ident) {
    object <- StashIdent(object = object, save.name = "old.ident")
    object <- SetAllIdent(object = object, id = "Phase")
  }
  return(object)
}

# ---------------------------------------------------------------------------- #
Imprint_Scoring = function(object, paternal.genes, maternal.genes)
{
  enrich.name <- "Imprinted"
  genes.list <- list(Paternal.Score = paternal.genes, Maternal.score = maternal.genes)
  object.cc <- AddModuleScore(object = object, genes.list = genes.list,
                              enrich.name = enrich.name, ctrl.size = min(vapply(X = genes.list,
                                                                                FUN = length, FUN.VALUE = numeric(1))))
  cc.columns <- grep(pattern = enrich.name, x = colnames(x = object.cc@meta.data))
  cc.scores <- object.cc@meta.data[, cc.columns]
  assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                                 first = "P", second = "M", null = "X") {
    if (all(scores <= 0)) {
      return(null)
    }
    else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  })
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments),
                     by = 0)
  colnames(x = cc.scores) <- c("rownames", "Paternal.Score", "Maternal.Score",
                               "Imprint")
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c("Paternal.Score", "Maternal.Score","Imprint")]
  object <- AddMetaData(object = object, metadata = cc.scores)

  return(object)
}

# ---------------------------------------------------------------------------- #
#' Seurat_Meta_Counts
#'
#' @param seu Seurat object
#' @param col Meta data column
#' @param type type of the matrix (data, raw.data)
#' @param norm whether to normalize the matrix
#' @param log boolean, whether to log2(x+1)
Seurat_Meta_Counts = function(
    seu,
    col,
    data.type = 'raw.data',
    norm = TRUE,
    log  = TRUE
){
    if(class(seu) != 'Seurat')
        stop('seu is not a Seurat object')

    if(!any(col %in% colnames(seu@meta.data)))
        stop('column is not in the meta.data')

    mat = slot(seu, data.type)
    cnams = unique(seu@meta.data[[col]])
    lmat = lapply(setNames(cnams, cnams), function(x){
        message(x)
        Matrix::rowSums(mat[,seu@meta.data[[col]]==x, drop=FALSE])
    })
    dmat = data.frame(lmat)
    if(norm)
        dmat = t(t(dmat)*(1e5/colSums(dmat)))

    if(log)
        dmat = log2(dmat+1)
    return(dmat)


}


# ---------------------------------------------------------------------------- #
# converts the seurat to single cell experiment
SeurateToSingleCellExperiment = function(
  seu,
  annot,
  gene_subset = NULL
){

    if(is.null(gene_subset))
        gene_subset = rownames(seu@data)

    rowData = S4Vectors::DataFrame(
        gene_id   = gene_subset,
        gene_name = annot[match(gene_subset, annot$gene_id),]$gene_name)
     rownames(rowData) = rowData$gene_id
     colData = S4Vectors::DataFrame(seu@meta.data)
     counts  = seu@raw.data[gene_subset, ]
     colnames(counts) = as.character(colnames(counts) )
     rownames(counts) = as.character(rownames(counts) )
     logcounts = seu@data[gene_subset, ]
     colnames(logcounts)   = as.character(colnames(logcounts) )
     rownames(logcounts)   = as.character(rownames(logcounts) )
     scale_data = seu@scale.data[gene_subset, ]
     colnames(scale_data)   = as.character(colnames(scale_data) )
     rownames(scale_data)   = as.character(rownames(scale_data) )

  sce =  SingleCellExperiment::SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(
      counts    = counts,
      logcounts = logcounts,
      scaleData = scale_data
    ))


  for (dr in names(seu@dr)) {
    SingleCellExperiment::reducedDim(sce, toupper(x = dr)) = slot(
      object = slot(object = seu, name = "dr")[[dr]],
      name = "cell.embeddings"
    )
  }
  return(sce)
}


# ---------------------------------------------------------------------------- #
# Should locate all cluster specific genes
FindAllMarkersIntersection = function(
  object,
  genes.use = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  print.bar = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  return.thresh = 1e-2,
  do.print = FALSE,
  random.seed = 1,
  min.cells.gene = 3,
  min.cells.group = 3,
  latent.vars = NULL,
  assay.type = "RNA",
  ...
) {
  data.1 <- GetAssayData(object = object,assay.type = assay.type,slot = "data")
  idents.all <- sort(x = unique(x = object@ident))

  genes.all <- list()
  for (i in 1:length(x = idents.all)) {

    ident_a = idents.all[i]
    set_reduced = setdiff(idents.all, ident_a)

    genes.de <- list()
    for(j in 1:length(set_reduced)){
      ident_b = set_reduced[j]

      genes.de[[j]] <- tryCatch(
        {
          FindMarkers(
            object = object,
            assay.type = assay.type,
            ident.1 = ident_a,
            ident.2 = ident_b,
            genes.use = genes.use,
            logfc.threshold = logfc.threshold,
            test.use = test.use,
            min.pct = min.pct,
            min.diff.pct = min.diff.pct,
            print.bar = print.bar,
            min.cells.gene = min.cells.gene,
            min.cells.group = min.cells.group,
            latent.vars = latent.vars,
            max.cells.per.ident = max.cells.per.ident,
            ...
          )
        },
        error = function(cond){
          return(NULL)
        }
      )
      if (do.print) {
        message(paste("Calculating cluster", idents.all[i]))
      }
    }
    genes.de = lapply(genes.de, function(x)subset(x, avg_logFC > 0))
    genes.de = genes.de[sapply(genes.de, nrow) > 0]
    gnams = Reduce(intersect, lapply(genes.de, function(x)rownames(x)))
    if(!is.null(gnams)){
      tab = genes.de[[1]]
      tab$cluster = ident_a
      tab = tab[order(tab$p_val, -tab$avg_logFC), ]
      tab$gene_id = rownames(tab)
      genes.all[[i]] = subset(tab, gene_id %in% gnams)
    }

  }
  gde.all = do.call(rbind, genes.all)
  return(gde.all)
}
