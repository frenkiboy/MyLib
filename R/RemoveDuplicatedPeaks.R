#!/usr/bin/Rscript
### INFO: Removes duplicates in Macs2 peak calls
### DATE: 151909
### AUTHOR: frenkiboy

# {0} TEST DATA
#/{0} TEST DATA

# {1} LIBRARIES
suppressMessages(library(GenomicRanges))
suppressMessages(library(genomation))
suppressMessages(library(argparser))
suppressMessages(library(stringr))

#/{1} LIBRARIES


# {2} CODE
    # {{1}} FUNCTIONS
    #/{{1}} FUNCTIONS
    
    
    # {{2}} INPUT VARIABLES 
    
        # {{{1}}} ARGUMENTS
        message('Parsing arguments...')
        p <- arg_parser("Removes duplicated peak regions from the peak file")
        p <- add_argument(p, "--input", help="input file",  default=NULL)
        p <- add_argument(p, "--output", help="output file", default=NULL)
        p <- add_argument(p, "--sortby", help="Columnt to sort by", default='qvalue')
        args = parse_args(p, argv = commandArgs(trailingOnly = TRUE))
        
        #/{{{1}}} ARGUMENTS
    
        # {{{2}}} SCRIPT PARAMS
        #/{{{2}}} SCRIPT PARAMS
    
    #/{{2}} INPUT VARIABLES
    
    
    # {{3}} MAIN
        
        message('Loading files...')
        if(is.null(args$i))
            stop('input file is not designated')

        message('Doing stuff...')
        infile = args$i
        peaks = readNarrowPeak(infile)
        peaks = peaks[order(-values(peaks)[[args$s]])]
        peaks.sel = unique(peaks)
        peaks.sel = sort(peaks.sel)		
		
        message(paste('orig.peaks:',length(peaks),'red.peaks:',length(peaks.sel),'perc:',round(length(peaks.sel)/length(peaks),2)))
        
        message('Writing files...')
        outfile= args$o
        if(is.null(outfile))
            outfile=str_replace(infile,'narrowPeak','uniq.narrowPeak')
        

        d = as.data.frame(peaks.sel)
        d = cbind(d[,c('seqnames','start','end','name','score')],'.',d[,c('signalValue','pvalue','qvalue','peak')])
        write.table(d, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
		
		message('Done! Have a nice day!')

    #/{{3}} MAIN
#/{2} CODE
