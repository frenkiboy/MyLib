# blatantly stolen from klmr/decorator
`%@%` = function (decorator, f) UseMethod('%@%')

`%@%.default` = function (decorator, f)
    stop(deparse(substitute(decorator)), ' is not a decorator')

`%@%.decorator` = function (decorator, f) {
    pretty_decorators = as.list(match.call())[-1]
    # Patch delayed decorator.
    if (! is.null({pretty_patched = attr(decorator, 'calls')}))
        pretty_decorators = c(pretty_patched, pretty_decorators[-1])

    # Handle operator precedence so that we can chain decorators.
    if (inherits(f, 'decorator'))
        .delayed_decorate(decorator, f, pretty_decorators)
    else
        prettify(decorator(f), f, pretty_decorators[-length(pretty_decorators)])
}

decorator = function (f)
    structure(f, class = 'decorator')

decorator = decorator(decorator)

print.decorated = function (x, useSource = TRUE, ...) {
    bare = function (f) {
        bare = unclass(f)
        attr(bare, 'decorators') = NULL
        bare
    }

    fun_def = capture.output(print.function(bare(x), useSource = useSource, ...))
    for (decorator in attr(x, 'decorators'))
        cat(deparse(decorator), '%@%\n')
    cat(fun_def, sep = '\n')
    invisible(x)
}

modules:::register_S3_method('print', 'decorated', print.decorated)

prettify = function (f, original, decorator_calls) {
    attr(f, 'srcref') = pretty_code(original)
    attr(f, 'decorators') = decorator_calls
    class(f) = c(class(f), 'decorated')
    f
}

pretty_code = function (f) {
    srcref = attr(f, 'srcref')
    if (is.null(srcref)) body(f) else srcref
}

.delayed_decorate = function (d1, d2, decorator_calls)
    structure(decorator(function (f) d1(d2(f))), calls = decorator_calls)

#
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
#'
#'
#' l = list(x=1:10, y=5:10)
#' f = cacheFile('./') %@% function(x)print(x)
#' lapply(l, f)
#'

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
        
        # checks whether there are any extra ... arguments
        # converts them into a list
        .dotind = names(.anames) == '...'
        if(any(.dotind)){
          .anames = .anames[!.dotind]
        }

        # evaluates global variables from .anames
        if(length(args) > 0){
            for(i in 1:length(.anames)){
                if(is.call(.anames[[i]]) | is.name(.anames[[i]])){
                  .eval = eval(.anames[[i]], envir=parent.frame())
                  if(is.null(.eval))
                    .eval = list(NULL)
                  .anames[[i]] = .eval
                }
            }
        }
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
