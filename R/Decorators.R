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
#' @examples
#' source(file.path(lib.path, 'Decorate.R'))
#' f = function(x=1,y=2)x
#' g = cacheFile('./') %@% f
#' g(2,3)
#' g(2,3)
#'
#' g(x=2,y=3)
#' g(x=2,y=3)

#' f = function(x=1,y=2){x;print(y);x+y}
#' g = cacheFile('./') %@% f
#' g(2,3)


source(file.path(lib.path, 'Decorate.R'))
cacheFile = function(inpath)decorator %@% function(f){

    library(digest)
    argnames = head(as.list(args(as.list(environment())[[1]])),-1)
    body = as.list(body(f))
    function(..., load=TRUE, .anames = argnames, .fbody=body){

        fcall = as.list(match.call())
        fname = fcall[[1]]
        args  = head(fcall[-1],-1)

        # finds out the function arguments and creates the argument hash
        if(!is.null(names(args))){
            named_args = setdiff(names(args),'')
            if(!is.null(named_args))
                for(i in named_args)
                    .anames[[i]] = args[[i]]

            pos_args = which(names(args) == '')
            if(length(pos_args) > 0)
                for(i in pos_args)
                    .anames[[i]] = args[[i]]
        }else{
            for(i in seq_along(args))
                .anames[[i]] = args[[i]]
        }
        hashlist = list(anames = .anames, body=.fbody)
        args_hash = digest(hashlist)

        outfile = file.path(inpath,paste(fname, args_hash, 'rds', sep='.'))
        if(load && file.exists(outfile)){
            print(paste0(fname,': Returning loaded data ...'))
            print(outfile)
            readRDS(outfile)$dat
       }else{
            print(paste0(fname,': Running function ...'))
            dat = f(...)

            saveRDS(list(dat=dat, args=.anames, body=.fbody), outfile)
            dat
        }
    }
}
