#'---
#'author: "Vedran Franke"
#'output:
#'  knitrBootstrap::bootstrap_document:
#'    title: "%TITLE"
#'    theme.chooser: TRUE
#'    highlight.chooser: TRUE
#'---



#+ global_options, include=FALSE
# needed for knitrBootstrap to hide code and messages with toggle on/off option
knitr::opts_chunk$set(bootstrap.show.code=TRUE,bootstrap.show.message=FALSE, fig.width=12, fig.height=16)

#+ load_local, include=FALSE
lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
source(file.path(lib.path, 'FileLoader.R'))

#+ load_global, include=FALSE
library(data.table)
library(stringr)
library(doMC)