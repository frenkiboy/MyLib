### INFO: Functions for exporting of various formats
### DATE: 15.03.2010.
### AUTHOR: v+

lib.path=file.path(Sys.getenv('HOME'),'bin/MyLib/RFun')
source(file.path(lib.path, 'ScanLib.R'))

# {1} 
    # --------------------------------------------------- #
	WriteBed = function(tab, path){
        
	   write.table(tab, path, row.names=F, col.names=F, quote=F, sep='\t')
	}

	# --------------------------------------------------- #
    # Exports a GRanges object to a bed file
	WriteRegsToBed = function(r, outpath=NULL, outname=NULL, stranded=FALSE){
	    
        if(is.null(outpath))
            outpath='./'
        
        if(is.null(outname))
	        outname='Regs'
	   
        if(!stranded){
            message('Exporting nonstranded ...')
            outfile = file.path(outpath, DateNamer(paste(outname,'bed',sep='.')))
            write.table(as.data.frame(r)[,1:3], outfile, row.names=F, col.names=F, quote=F, sep='\t')
        
        }else{
    	    for(i in c('+','-')){
    	        message('Exporting',i,'...')
    	        s = ifelse(i=='+','p','m')
    	        outfile = file.path(outpath, DateNamer(paste(outname,s,'bed',sep='.')))
    	        write.table(as.data.frame(r[strand(r) == i])[,1:3], outfile, row.names=F, col.names=F, quote=F, sep='\t')
    	    }
        }
	}
	
    
#/{1} 
