### designates the transcripts as outer, internal or body
#copyright V+


library(rparallel)
Creator <- function(bed, opr, ipr)
{
	ids = unique(bed[,4])
	result = data.frame(matrix(ncol = (ncol(bed2)+1)))
	result = na.omit(result)
	
	if( "rparallel" %in% names( getLoadedDLLs()) )
    {
       runParallel( resultVar="result", resultOp="rbind", nWorkers = 3 )
    }
    else 
	{
		for(i in 1:length(ids))
		{
			cat(i,"\n")
			trans = bed[bed[,4] == ids[i],]
			gene = data.frame(matrix(ncol = ncol(trans)))
			gene = na.omit(gene)
			strand = as.character(unique(trans[,5]))
			if ( strand == "+")
			{
				gene = rbind(gene,trans[trans[,2] == min(trans[,2]),])
				gene = reducer(gene,1)
				gene[1,3] = gene[1,2] + opr
				gene[1,2] = gene[1,2] - opr
				gene = rbind(gene,trans[trans[,3] == max(trans[,3]),])
				gene = reducer(gene,2)
				gene[2,2] = gene[1,2]
				trans = trans[trans[,2] != min(trans[,2]),]
				trans[,3] = trans[,2] + ipr
				trans[,2] = trans[,2] - ipr
				
			}else
			{
				gene = rbind(gene, trans[trans[,3] == max(trans[,3]),])
				gene = reducer(gene,1)
				gene[1,2] = gene[1,3] - opr
				gene[1,3] = gene[1,3] + opr
				gene = rbind(gene,trans[trans[,2] == min(trans[,2]),])
				gene = reducer(gene,2)
				gene[2,3] = gene[1,3]
				trans = trans[trans[,3] != max(trans[,3]),]
				trans[,2] = trans[,3] - ipr
				trans[,3] = trans[,3] + ipr
				
			}
			
			gene = rbind(gene,trans)
			design = rep('I', times = nrow(trans))
			design = c('O','B',design)
			gene = cbind(gene, design)
			names(gene) = names(result)
			result = rbind(result, gene)
		}
	}
	return(result)
}


reducer <- function(frame, rows)
{
	frame = frame[1:rows,]
	return(frame)
}
