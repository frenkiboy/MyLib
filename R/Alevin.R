# ------------------------------------------------------------------------------------ #
# takes the alevin output folder and outputs a matrix
# taken from somewhere else, but forgot the source
ReadAlevin = function( base.path = NULL ){
    if (! dir.exists(base.path )){
      stop("Directory provided does not exist")
    }

    barcode.loc = file.path( base.path, "alevin/quants_mat_rows.txt" )
    gene.loc    = file.path( base.path, "alevin/quants_mat_cols.txt" )
    matrix.loc  = file.path( base.path, "alevin/quants_mat.csv" )
    if (!file.exists( barcode.loc )){
      stop("Barcode file missing")
    }
    if (! file.exists(gene.loc) ){
      stop("Gene name file missing")
    }
    if (! file.exists(matrix.loc )){
      stop("Expression matrix file missing")
    }
    matrix = as.matrix(read.csv( matrix.loc, header=FALSE))
    matrix = t(matrix[,1:ncol(matrix)-1])

    cell.names = readLines( barcode.loc )
    gene.names = readLines( gene.loc )

    colnames(matrix) = cell.names
    rownames(matrix) = gene.names
    matrix[is.na(matrix)] <- 0
    return(matrix)
}

# ------------------------------------------------------------------------------------ #
# takes the alevin output folders and merges them into a Seurat object
Alevin_To_Seurat = function(
  path         = NULL,
  min.cells    = 3, 
  min.genes    = 10, 
  project_name = 'SingleCell'
){

   
  if(is.null(path))
    stop('Please provide a path')
    
  library(Seurat)
  library(dplyr)
  library(stringr)
  dirs = list.files(path)
  lseu = list()
  for(i in 1:length(dirs)){
    name = dirs[i]
    cname = str_replace(name,'-5','')
    message(cname)
    alv.data = ReadAlevin(file.path(path, name))
    colnames(alv.data) = paste(cname, colnames(alv.data), sep='_')
    seu = CreateSeuratObject(raw.data = alv.data, min.cells = min.cells, min.features = min.genes, project = project_name)
    lseu[[name]] = seu
  }
  sem = Reduce(function(x,y)MergeSeurat(x,y), lseu)
  return(sem)
}
