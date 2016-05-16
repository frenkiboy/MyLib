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
	source('~/.Renviron')
    	source("http://bioconductor.org/biocLite.R")
	source(file.path(Sys.getenv("MYLIB"),"RFun/Startup.R"))
}
