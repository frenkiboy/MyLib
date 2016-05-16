### INFO: Parses FastQC output
### DATE: %DATE
### AUTHOR: %AUTH

# -------------------------------------------------------- #
parse_FastQCfile = function(path=NULL){
    # function that reads output from fastqc program
    # input: a zipped file e.g. SRR1731137_fastqc.zip
    # output: a character vector of pass/fail/warn
    
    require(stringr)
    require(data.table)
    if(!file.exists(path))
        stop('The input file does not exist')
    
    runid = str_replace(basename(path), '_fastqc.zip', '')
    a = try(read.csv(unz(path, paste(runid, "_fastqc/fastqc_data.txt", sep="")), 
                     stringsAsFactors = FALSE))
    
    desc = data.table(do.call(rbind, strsplit(a[3:9,1],'\\t')))
    
    stats = grep('>>', a[,1], value=TRUE)
    stats = grep('END', stats, value=TRUE, invert=TRUE)
    stats = str_replace(stats,'>>','')
    stats = data.table(do.call(rbind, strsplit(stats,'\\t')))
    
    tab = rbind(desc, stats)
    colnames(tab) = c('description',str_replace(desc$V2[1],'.fastq.gz',''))
    
    return(tab)
}

# -------------------------------------------------------- #
parse_FastQC = function(path){
    
    source(file.path(Sys.getenv('MYLIB'),'RFun','ScanLib.R'))
    files = list.files(path, full.names=TRUE, pattern='fastqc.zip', )
    l = lapply(files, parse_FastQCfile)
    ord = data.frame(ord=as.character(l[[1]]$description))
    mqc = MergeDataTable(l, key='description', all=TRUE)
    mqc = mqc[match(ord$ord, mqc$description),]
    return(mqc)
    
}









