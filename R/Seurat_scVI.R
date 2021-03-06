# ---------------------------------------------------------------------------- #
# run scVI imputation on a seurat data set
run_scVI = function(
  seu           = NULL,
  conda         = '/home/vfranke/bin/Software/miniconda3/bin/conda',
  n_epochs      = 100,
  lr            = 1e-3,
  n_batches     = 0,
  batch_indices = NULL,
  use_cuda      = FALSE,
  n_latent      = 20,
  n_layers      = 2,
  n_hidden      = 128,
  dispersion    = 'gene',
  reconstruction_loss = 'zinb'
){

  if(is.null(seu) || class(seu) != 'Seurat')
    stop('Please provide a valid seurat object')

  suppressPackageStartupMessages({
    library(reticulate)
  })

  message('Writing counts csv ...')
    tmpfile = tempfile()
    save_path = paste0(dirname(tmpfile),'/')
    cnts = as.matrix(GetAssayData(seu,'counts'))
    write.csv(cnts, tmpfile)
    system(paste('gzip', tmpfile))

    batch_name = NULL
    if(!is.null(batch_indices)){
        message('Writing batches csv ...')
        batch_file = paste(tmpfile,'batch', sep='_')
        batch_name = basename(batch_file)
        write.csv(batch_indices, batch_file)
    }
  message('Import scVI ...')
    use_condaenv(condaenv='scvi', conda=conda)
    scvi = import('scvi')

   message('Import Data ...')
     filename = paste0(basename(tmpfile),'.gz')
     local_csv_dataset = scvi$dataset$CsvDataset(
       filename       = filename,
       save_path      = save_path,
       compression    = 'gzip',
       new_n_genes    = nrow(cnts),
       batch_ids_file = batch_name
     )

   message('Run scVI ...')
     vae = scvi$models$VAE(
               n_input  = as.integer(nrow(cnts)),
               n_batch  = as.integer(n_batches),
               n_latent = as.integer(n_latent),
               n_hidden = as.integer(n_hidden),
               n_layers = as.integer(n_layers),
               dispersion = dispersion,
               reconstruction_loss = reconstruction_loss)

     trainer = scvi$inference$UnsupervisedTrainer(
       vae,
       local_csv_dataset,
       train_size = 0.75,
       use_cuda   = use_cuda,
       frequency  = 5)

     trainer$train(n_epochs=as.integer(n_epochs), lr=lr)

   message('Return scVI ...')
     latent = trainer$get_all_latent_and_imputed_values()
     latent$imputed_values = t(latent$imputed_values)
     rownames(latent$latent) = colnames(GetAssayData(seu,'counts'))
     latent$latent = as.data.frame(latent$latent)
     colnames(latent$latent) = paste('vae', colnames(latent$latent), sep='_')

     rownames(latent$imputed) = rownames(GetAssayData(seu,'counts'))
     colnames(latent$imputed) = colnames(GetAssayData(seu,'counts'))
     return(latent)
}


# ---------------------------------------------------------------------------- #
# runs the Seurat analysis with scVIS imputed datasets
# This function should be used within other functions which
# provide the corresponding seurat input
# dispersion: gene, gene-batch, gene-lael, gene-cell
Run_Seurat_scVI = cacheFile(path_RDS) %@% function(
  seu,
  algorithm           = 2,
  k.param             = 15,
  dims.use            = 10,
  n_batches           = 0,
  n_latent            = 20,
  subsets_clustering  = FALSE,
  batch_indices       = NULL,
  winsor              = 2,
  n_layers            = 2,
  n_hidden            = 128,
  variable_selection_method = 'vst',
  dispersion          = 'gene',
  reconstruction_loss = 'zinb',
  nfeatures           = 1000
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(data.table)
    library(Rtsne)
  })
  source.lib('Seurat.R')
  source.lib('ScanLib.R')

  message('scVI ...')
  latent = run_scVI(
    seu,
    batch_indices = batch_indices,
    n_batch       = n_batches,
    n_latent      = n_latent,
    n_layers      = n_layers,
    n_hidden      = n_hidden,
    dispersion    = 'gene',
    reconstruction_loss = 'zinb'
  )

  lat = latent$latent
  imp = latent$imputed
  seu@meta.data = cbind(seu@meta.data, lat)
  message('Scaling imputed ...')
  quant  =  quantile(as.matrix(imp), seq(0,1,0.01))
  winsor =  head(tail(quant,winsor),1)
  imp[imp > winsor] = winsor
  colnames(imp) = as.character(colnames(imp))
  rownames(imp) = as.character(rownames(imp))
  seu[['scvi']] = CreateAssayObject(data = imp)
  DefaultAssay(seu) = 'scvi'
  seu = ScaleData(object = seu, assay='scvi')
  seu = FindVariableFeatures(
      seu,
      do.plot=FALSE,
      selection.method = variable_selection_method,
      nfeatures = nfeatures
  )

  message('PCA ...')
  seu = RunPCA(seu, do.print=FALSE, dims = dims.use, verbose=FALSE, assay='scvi')

  message('VAE ...')
  pcmat = seu@meta.data %>% dplyr::select(contains('vae'),-contains('tsne'))
  colnames(pcmat) = paste0('PC', 1:ncol(pcmat))
  seu[['vae']] = CreateDimReducObject(
    embeddings = as.matrix(pcmat),
    assay      = 'scvi',
    key        = 'VAE_'
    )


  message('TSNE ...')
    mdims = min(dims.use, ncol(pcmat))
    seu = RunTSNE(seu,do.print=FALSE, dims = 1:mdims, reduction.use='vae', verbose=FALSE)

  message('UMAP ...')
    seu = RunUMAP(seu, do.print=FALSE, dims = 1:mdims, reduction.use='vae', do.plot=FALSE, verbose=FALSE)

  seu@meta.data = seu@meta.data  %>%
    cbind(seu$tsne@cell.embeddings) %>%
    cbind(seu$umap@cell.embeddings)

  seu = FindNeighbors(object  = seu, dims = 1:mdims, verbose = FALSE)
  res = seq(0.2, 2, 0.2)
  seu = FindClusters(
    seu,
    resolution = res,
    k.param    = 10,
    prune.SNN  = 0.1,
    dims       = 1:mdims,
    algorithm  = algorithm,
    reduction.type='vae'
  )

  return(seu)
}
