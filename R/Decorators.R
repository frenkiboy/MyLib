# ---------------------------------------------------------------------------- #
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
#' lib.path = "https://raw.githubusercontent.com/frenkiboy/MyLib/master/R"
#' source(file.path(lib.path, 'Decorate.R'))
#' f = function(x=1,y=2)x
#' g = cacheFile('./') %@% f
#' g(2,3)
#' g(2,3)
#'
#' g(x=2,y=3)
#' g(x=2,y=3)

# All should give the same hash
#' a = 3
#' f = function(x=1,y=a)x
#' g = cacheFile('./') %@% f
#' g(2)
#' g(2,3)
#' g(x=2,y=3)
#' g(x=2,3)
#' g(2,y=3)


#' f = function(x=1,y=2){x;print(y);x+y}
#' g = cacheFile('./') %@% f
#' g(2,3)

source(file.path(lib.path, 'Decorate.R'))
cacheFile = function(inpath)decorator %@% function(f){
    
    library(digest)
    argnames = head(as.list(args(as.list(environment())[[1]])),-1)
    body = lapply(as.list(body(f)), as.character)
    function(..., .load=TRUE, .anames = argnames, .fbody=body){

        # -------------------------------------------------------------------- #
        fcall = as.list(match.call())
        
        # extracts the function name
        fname = fcall[[1]]
        
        # removes the funciton name from the call
        args  = fcall[-1]
        if(!is.null(names(args)) && any(names(args) == '.load'))
            args  = args[names(args) != '.load']

        # replaces default arguments with set arguments
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
        
        # evaluates global variables from .anames
        .anames = lapply(.anames, function(x){
            if(is.name(x)){
                eval(x)
            }else{
                x
            }
        })

        # -------------------------------------------------------------------- #
        # creates the argument hash
        hashlist = list(anames = .anames, body = .fbody)
        args_hash = digest(hashlist, algo='md5')
        message(args_hash)
        
        # -------------------------------------------------------------------- #
        outfile = file.path(inpath,paste(fname, args_hash, 'rds', sep='.'))
        if(.load && file.exists(outfile)){
            message(paste0(fname,': Returning loaded data ...'))
            message(outfile)
            readRDS(outfile)$dat
       }else{
           message(paste0(fname,': Running function ...'))
            dat = f(...)

            saveRDS(list(dat=dat, args=.anames, body=.fbody), outfile)
            dat
        }
    }
}

# ---------------------------------------------------------------------------- #
# deorator tests
# testthat::test_that('cacheDecorator'){
#     
#     # 1. test that the cache works 
#     f = function(x=1,y=2)x
#     g = cacheFile(tempdir()) %@% f
#     t1 = g(2,3, .load=FALSE)
#     testthat::expect_message(g(2,3, .load=FALSE), "g: Running function ...")
#     
#     t2 = g(2,3, .load=TRUE)
#     testthat::expect_message(g(2,3, .load=TRUE), "g: Returning loaded data ...")
#     
#     testthat::expect_equal(t1, t2)
#     
#     # 2. test that named arguments don't cause a new hash
#     t3 = g(x=2,y=3, .load=TRUE)
#     testthat::expect_equal(t1, t3)
#     
#     # 3. no aruments should not start a new hash
#     t4.1 = g(1,2)
#     t4.2 = g(.load=TRUE)
#     testthat::expect_message(g(.load=TRUE), "g: Returning loaded data ...")
#     testthat::expect_equal(t4.1, t4.2)

      # 4. different calling arguments
      # 5. different function body

# }

