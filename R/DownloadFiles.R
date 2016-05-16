### INFO: R Script
### DATE: 29.04.2013
### AUTHOR: Vedran Franke
rm(list=ls())


# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
source(file.path(lib.path, 'FileLoader.R'))
source(file.path(lib.path, 'FormatConverters.R'))
source(file.path(lib.path, 'BamWorkers.R'))
source(file.path(lib.path, 'ScanLib.R'))
library(data.table)
library(biomaRt)
library(GEOquery)
library(doMC)
library(sra)
library(SRAdb)

#/{1} LIBRARIES


# {2} CODE
    # {{1}} FUNCTIONS
    #/{{1}} FUNCTIONS
    
    
    # {{2}} INPUT VARIABLES 
    
    # {{{1}}} PATH VARIABLES
    
    
     outdir='/data/akalin/Base/AccessoryData/'
    
    sradb = '~/Tmp/SRAmetadb.sqlite'
    
    
    #/{{{1}}} PATH VARIABLES
    
    # {{{2}}} SCRIPT PARAMS
     outname='GSE41091_Nekrasov_2012_NSMB_CTCF_ChipSeq'
    
     registerDoMC(20)
    
    #/{{{2}}} SCRIPT PARAMS
    
    #/{{2}} INPUT VARIABLES
    
    
    # {{3}} MAIN
    
        outpath = file.path(outdir, outname, 'Raw','SRA')
        dir.create(outpath, showWarnings=FALSE, recursive=TRUE)
        
        id = unlist(strsplit(outname,'_'))[1]
        geo = getGEO(exps)[[1]]
        p = pData(geo)
        d = data.frame(sample=p$title, id=basename(p$supplementary_file_1))
        
        if(!file.exists(sradb))
            getSRAdbFile(destdir=dirname(sradb))
        
        con = dbConnect(dbDriver("SQLite"),sradb)

        sra = getSRAinfo(in_acc=d$id, sra_con=con)
        sra = sra[c('experiment','ftp')]
        m = merge(d, sra, by.x='id', by.y='experiment')
        
        foreach(i = 1:nrow(m))%dopar%{
            
            print(i)
            file=m$ftp[i]
            download.file(file, file.path(outpath, paste(m$sample[i], 'sra', sep='.')))
        }
    
    #/{{3}} MAIN
#/{2} CODE



