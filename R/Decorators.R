# ------------------------------------------------------------------ #
#' checkLoad - function decorator to save function output
#' checkLoad is a function decorator which saves the output of a function to an RDS object
#' TODO
#' [] add hash sum check for object and parameters
#'
#'
#' @param f function to be executed
#' @param inpath location where to save the function output
#' @param load boolean, whether to load the saved output or force execution

#' @return f function output
#' @export
source(file.path(lib.path, 'Decorate.R'))
cacheFile = decorator %@% function(f, inpath='./'){

    library(digest)
    fname = deparse(substitute(f))
    function(..., load=TRUE){
        print(match.call())
        outfile = file.path(inpath,paste(fname,'rds',sep='.'))
        if(load && file.exists(outfile)){
            message('Returning loaded data ...')
            readRDS(outfile)
       }else{
           message('Running function ...')
            dat = f(...)

            saveRDS(dat, outfile)
            dat
        }
    }
}
