### INFO: Functions for working with genes
### DATE: 27.01.2011
### AUTHOR: v+


# {1} LIBRARIES
#/{1} LIBRARIES


# {2} CODE

	# -------------------------------- #
	# {{1}} GenomicRegionMaker
	# takes a gene set and creates a bed table with the coordinates corresponding to the promotor, body, upstream and downstream
	GenomicRegionMaker = function(trans, prom.down=1000, prom.up=1000, reg.up=5000, reg.down=5000){

		if(class(trans) != 'GRangesList')
			stop('genes needs to be a GRanges class')
		
		
		cat('Making the TSS...\n')
        body = unlist(range(trans))
        tss = promoters(body, downstream=prom.down, upstream=prom.up)
		
		cat('Making the upstream regions...\n')
		upstream = resize(body, width=width(body)+reg.up, fix='end')

		cat('Making the downstream regions...\n')
		downstream = resize(body, width=width(body)+reg.down, fix='start')
		
		list(tss = tss,
             exon = unlist(trans),
			 intron = body, 
			 upstream = upstream,
			 downstream = downstream)
	}

	#/{{1}}
	
	# -------------------------------- #
	# {{2}} CpgDesignator
	CpGDesignator = function(genes, cpg, downstream=500, upstream=1000){
	
		genes = BedKey(genes)
		genes.small = (genes$end - genes$start) <= downstream
		genes.p = genes$strand == '+'
		genes.m = genes$strand == '-'
		


		genes.prom = genes
		genes.prom$end[genes.p]   = genes.prom$start[genes.p] + downstream
		genes.prom$start[genes.p] = genes.prom$start[genes.p] - upstream
		genes.prom$start[genes.m] = genes.prom$end[genes.m]   - downstream
		genes.prom$end[genes.m]   = genes.prom$end[genes.m]   + upstream
		# defines the promotor area for genes < downstream
		genes.prom$end[genes.small & genes.p] = genes$end[genes.small & genes.p]
		genes.prom$start[genes.small & genes.m] = genes$start[genes.small & genes.m]

		genes.body = genes
		genes.body$start[genes.p] = genes.body$start[genes.p] + downstream + 1 
		genes.body$end[genes.m]   = genes.body$end[genes.m] - downstream - 1
		# defines the gene body for the small genes
		genes.body$start[genes.small] = genes$start[genes.small]
		genes.body$end[genes.small]   = genes$end[genes.small]

		# finds the overlaps between the genes and cpgs
		genes.cpg.prom = FeatureOverlapBED(genes.prom, cpg)
		genes.cpg.body = FeatureOverlapBED(genes.body, cpg)

		# annotates the genes by cpg overlapp
		genes$cpg.prom = 'no'
		genes$cpg.body = 'no'
		genes$cpg.prom[genes$key %in% genes.cpg.prom$key.a] = 'yes'
		genes$cpg.body[genes$key %in% genes.cpg.body$key.a] = 'yes'
		
		genes = genes[,-grep('key',names(genes))]
		return(genes)
	}
	# genes = genes[,which(names(genes) != 'key')]

	# -------------------------------- #
	### Takes a gtf gene annotation and selects the longest transcript
    GtfSelectTranscript = function(gtf){
    
        library(data.table)
        library(GenomicRanges)
        library(stringr)
        
        if(!is.null(gtf$exon_id)){
            gtf = gtf[gtf$exon_id != 'CCDS']
            gtf = gtf[!str_detect(gtf$exon_id, 'mRNA')]
            gtf = gtf[!str_detect(gtf$exon_id, 'cdna')]
        }
        sl = seqlevels(gtf)
        seqlevels(gtf, pruning.mode='coarse') = sl[!(str_detect(sl,'HG') | str_detect(sl,'GL') | str_detect(sl,'MHC') | str_detect(sl,'HS') | str_detect(sl,'MT'))]
        d = data.table(as.data.frame(values(gtf)))
        d$strand = as.character(strand(gtf))
        d$width = width(gtf)
        d[d$strand == '+' , `:=`( COUNT = .N , IDX = 1:.N ) , by = transcript_id[strand == '+']]
        d[d$strand == '-' , `:=`( COUNT = .N , IDX = .N:.1) , by = transcript_id[strand == '-']]
        d.cnt = d[, c('gene_id','transcript_id','COUNT','width'), with=F]
        d.cnt = d[, list(COUNT=unique(COUNT), width=sum(width)), by=c('gene_id', 'transcript_id')]
        d.cnt = d.cnt[order(-d.cnt$COUNT, -d.cnt$width),]
        d.cnt = d.cnt[!duplicated(d.cnt$gene_id)]
        
        gtf.t = gtf[gtf$transcript_id %in% d.cnt$transcript_id]
        return(gtf.t)
    }
	
    
	
	# -------------------------------- #
	# {{3}} 
	### takes a set of gene transcripts and creates gene models
	MakeGeneModels = function(gffs, return='data.frame'){
	
		if(!class(gffs) == 'GRangesList')
			stop('gffs needs to be a GRangesList object')
			
		if(!return %in% c('data.frame','GRanges'))
			stop('return can only be data.frame or GRanges')
		
		chr =  as.character(unlist(runValue(seqnames(gffs))))
		start =  min(start(gffs))
		end =  max(end(gffs))
		strand = as.character(unlist(runValue(strand(gffs))))
		ens.gene.id = names(gffs)
		ex.width = sapply(width(gffs), sum)
		ex.num = elementLengths(gffs)
		genes = data.frame(chr =  chr,
						   start = start,
						   end = end,
					       strand = strand,
					       ens.gene.id = ens.gene.id, 
					       ex.width=ex.width,
					       ex.num=ex.num)
		
		if(return == 'data.frame'){
			return(genes)
		}else if(return == 'GRanges'){
			return(BedToGRanges(genes, values=T))
		}
	}
	

	# -------------------------------- #
	# {{4}}
	AnnotateRegions = function(region, l.pos, ordering=1:length(l.pos), nclass='Int'){
	
		if(! class(region) == 'GRanges')
			stop('region needs to be GRanges')

		if(! class(l.pos) == 'GRangesList')
			stop('region needs to be GRangesList')
			
		d = do.call(cbind, lapply(l.pos, function(x)ifelse(countOverlaps(region, x) == 0, 0, 1)))
		d = apply(t(t(d) * ordering), 1, function(x)ifelse(sum(x) == 0, 0, min(x[x!=0])))
		n = rep(nclass, length(d))
		n[d != 0] = names(l.pos)[d[d != 0]]
		return(n)
			
	}
	
#/{2} CODE



