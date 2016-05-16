# ---------------------------------------------------------------- #
# converts an ensembl output to gtf
EnsemblExonsTOGtf = function(path, oc=T){

	s = scan(path, what='character', sep='\n')
	a = strsplit(s, split='\t')
	a[[1]] = gsub(' \\(bp\\)','', a[[1]])
	a[[1]] = gsub('\\s+','.', a[[1]])

	d = do.call(rbind, a)
	d = data.frame(d, stringsAsFactors=F)
	names(d) = d[1,]
	d = d[-1,]
	d$Strand = ifelse(d$Strand == '1', '+', '-')
	d = d[,!grepl('UTR',names(d))]
	d = d[nchar(d$Chromosome.Name) <= 2,]
	d$Chromosome.Name[d$Chromosome.Name == 'MT'] = 'M'
	
	gtf = data.frame(
			chr = paste('chr',d$Chromosome.Name,sep=''),
			id = rep('mm9_ensGene', nrow(d)),
			feature = rep('exon', nrow(d)),
			start = d$Exon.Chr.Start,
			end = d$Exon.Chr.End,
			score = rep(0, nrow(d)),
			strand = d$Strand,
			frame = rep('.', nrow(d)),
			gene.id = paste('gene_id', paste('"',d$Ensembl.Gene.ID,'";',sep=''), sep=' '),
			transcript.id = paste('transcript_id', paste('"',d$Ensembl.Transcript.ID,'";',sep=''), sep=' ')
			)
	gtf = unique(gtf)
	attr(gtf, 'cord.sys') = 1		
	return(gtf)
}
# ---------------------------------------------------------------- #

# ---------------------------------------------------------------- #



# ---------------------------------------------------------------- #