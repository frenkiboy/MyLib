### INFO: Parses FastQC output
### DATE: %DATE
### AUTHOR: %AUTH

# -------------------------------------------------------- #
run_BBduk2= function(params.all){
    
    params = paste(with(subset(params.all, !param %in% c('bbduk','log')), paste(param,value,sep='=')), collapse=' ')
    command = paste(params.all$value[params.all$param=='bbduk'],
                    params,
                    paste0('2>',params.all$value[params.all$param == 'log']))
    system(command, intern=TRUE)
}


run_Trimming_BBduk2 = function(inpath, params.tab, nproc=16){
    
    if(basename(inpath) != 'Fastq')
        stop('Something is wrong with the folders')
    
    outpath = file.path(dirname(inpath),'Cleaned','Bbduk2')
        dir.create(outpath, showWarnings=FALSE, recursive=TRUE)
    
    infiles = data.frame(files=list.files(inpath, full.names=TRUE, recursive=TRUE, pattern='fastq'))
    infiles$set = as.numeric(as.factor(sub('.fastq.*','',basename(infiles$files))))
    infiles = split(infiles, infiles$set)
    
    require(doMC)
    registerDoMC(nproc)
    foreach(i = 1:length(infiles))%dopar{
        
        print(i)
        files = infiles[[i]]
        files = files[order(files$files),] 
        name = sub('[(_r\\d+)]?.fastq.*','',basename(files$files)[1])
        outdir = file.path(outpath, name)
            dir.create(outdir, showWarnings=FALSE)
        
        if(nrow(files) > 1){
            ins = data.frame(param=c('in','in2','out','out2'),
                             value=c(files$files[1],
                                     files$files[2],
                                     file.path(outdir, paste0(name,'r1.cl.fastq.gz')),
                                     file.path(outdir, paste0(name,'r2.cl.fastq.gz'))))  
             
        }else{
            ins = data.frame(param=c('in','out'),
                             value=c(files$files[1],
                                     file.path(outdir, paste0(name,'.cl.fastq.gz'))))  
        }
        
        
        statfiles = data.frame(param = c('stats','log'),
                               value = c(file.path(outdir, paste0(name,'.cl.stats.txt')),
                                         file.path(outdir, paste0(name,'.cl.log'))))
        params.all = rbind(params.tab, ins, statfiles)
        
        run_BBduk2(params.all)
    }
    
}








