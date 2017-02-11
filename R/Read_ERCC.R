read_ERCC = function(ercc.bamfiles, ercc.annot.path){

  library(Rsamtools)
		ercc.cnts = lapply(ercc.bamfiles, function(x){
													name = BamName(x)
                          d = idxstatsBam(x)[,c(1,3)]
													colnames(d) = c('id',name)
													data.table(d)
													})
		ercc.cnts = MergeDataTable(ercc.cnts,'id')
		ercc.cnts = ercc.cnts[-1,]

    stat = Rsamtools::idxstatsBam(ercc.bamfiles[1])[,c('seqnames','seqlength')]
    ercc.cnts = merge(stat, dbam, by.x='seqnames', by.y='id')
    setnames(ercc.cnts, 'seqnames', 'id')

		ercc.annot = read.table(ercc.annot.path, header=TRUE, sep='\t')[,c(2,3,4,5)]
		colnames(ercc.annot) = c('id','subgroup','cont.mix1','cont.mix2')
		ercc.all = merge(ercc.annot, ercc.cnts, by='id')
		return(ercc.all)
}
