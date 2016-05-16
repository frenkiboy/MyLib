ercc.annot.path = '/data/akalin/Base/ERCC/Annotation/ERCC92.txt'
read_ERCC = function(ercc.bamfiles){

		ercc.cnts = lapply(ercc.bamfiles, function(x){
													name = BamName(x)
													d = data.frame(do.call(rbind, strsplit(system(paste('samtools idxstats',x), intern=TRUE),'\\t')))[,c(1,3)]
													colnames(d) = c('id',name)
													data.table(d)
													})
		ercc.cnts = MergeDataTable(ercc.cnts,'id')
		ercc.cnts = ercc.cnts[-1,]
		ercc.annot = read.table(ercc.annot.path, header=TRUE, sep='\t')[,c(2,3,4,5)]
		colnames(ercc.annot) = c('id','subgroup','cont.mix1','cont.mix2')
		ercc.all = merge(ercc.annot, ercc.cnts, by='id')
		return(ercc.all)
}
