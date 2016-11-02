# INFO: functions that should be defined when R starts
# DATE: 05012015
# AUTH: vedran franke; vedran.franke@gmail.com
# fuctions

# -------------------------------------------------- #
# general options
options(stringsAsFactors=FALSE)



.First = function(){
	# -------------------------------------------------- #
	# source
	
	source("http://bioconductor.org/biocLite.R")
    
    varfile=list.files(getwd(), pattern='Config', full.names=TRUE)
    if(file.exists(varfile) && length(varfile) == 1)
	    eval(source(varfile), envir=parent.frame())
    
    try(source(file.path(Sys.getenv("MYLIB"),"RFun/Startup.R"))) 
}
