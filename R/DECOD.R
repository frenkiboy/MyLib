
	# Function for parsing DECOD discriminative motif discovery

	DECOD_PATH = '/home/vfranke/bin/Software/MotifDiscovery/DECOD-20111024.jar'
	# ------------------------------------------------------------------------ #
	# run DECOD
	
	runDECOD = function(posfile, negffile, outpath='./Decod.out', 
						w=7, nmotif=10, strand='forward', decod.path=DECOD_PATH){

							
		if(!file.exists(posfile))
			stop('posfile needs to exist')
			
		if(!file.exists(negfile))
			stop('negfile needs to exist')
			
		message('Running DECOD...')
		decod = paste('java -jar',decod.path,'-nogui')
		args = paste('-pos',     posfile,
					 '-neg',     negfile,
					 '-o',       outpath,
					 '-w', w,
					 '-nmotif', nmotif,
					 '-strand', strand)
		
		command = paste(decod, args)
		message(command)
		
		system(command, wait=FALSE, intern=FALSE)
		
		message('Returning results...')
		res = readDecod(outpath)
		invisible(return(res))
	}
	
	
	
	# ------------------------------------------------------------------------ #
	# read DECOD
	readDECOD = function(decod.file){

		decod.s = scan(decod.file, what='character', sep='\t')
        motind  = which(str_detect(decod.s, '>Motif'))
		decod.l = lapply(motind,function(x)decod.s[(x+1):(x+4)])
		decod.l = lapply(decod.l, function(x){
								m = str_replace(x,'^..','')
								m = str_replace(m,'\\[','')
								m = str_replace(m,'\\]','')
								m = do.call(rbind, strsplit(m, '\\s'))
								m = apply(m,2,as.numeric)
								colnames(m) = paste('p',1:ncol(m),sep='')
								rownames(m) = c('A','C','G','T')
								m
		})
		names(decod.l) = paste('M', 1:length(decod.l),sep='')
        
        mwhich = which(str_detect(decod.s,'Motif instances'))
        sep = data.frame(start = mwhich[seq(1,length(mwhich),2)] + 1, 
                         end   = mwhich[seq(2,length(mwhich),2)] - 1)
        sen = data.frame(start = mwhich[seq(2,length(mwhich),2)] + 1, 
                         end   = c(tail(motind,-1) - 1,length(decod.s)))
        lhits = list()
        for(i in 1:nrow(sep)){
            name = names(decod.l)[i]
            message(name)
            
            mpos = data.frame(matrix(decod.s[sep[i,1]:sep[i,2]], ncol=4, byrow=TRUE))
            colnames(mpos) = c('id','n1','seq','n2')
            mpos$id = str_replace(mpos$id,'>','')
            
            mneg = data.frame(matrix(decod.s[sen[i,1]:sen[i,2]], ncol=4, byrow=TRUE))
            colnames(mneg) = c('id','n1','seq','n2')
            mneg$id = str_replace(mneg$id,'>','')
            lhits[[name]]$pos = mpos
            lhits[[name]]$neg = mneg
            
        }
        lres = lapply(names(decod.l), function(x)
                                            list(mat=decod.l[[x]], 
                                                 pos=lhits[[x]]$pos,
                                                 neg=lhits[[x]]$neg))
        names(lres) = names(decod.l)
		return(lres)
	}
	
	# ------------------------------------------------------------------------ #
	annotateDECOD = function(decod.l, database.l, nhit=3){
	
        
        if(!all(sapply(decod.l, class) == 'matrix'))
            stop('the input file needs to be a list of matrices')
		require(MotIV)
		
        dbscore = generateDBScores(database.l,cc="PCC",align="SWU",nRand=1000,go=1,ge=0.5)

        matches=motifMatch(decod.l, database.l, DBscores=dbscore, cc="PCC", align="SWU", top=5, go=1, ge=0.5)
        bm = matches@bestMatch
        lbm = list()
        for(i in 1:length(bm)){
            
            name = bm[[i]]@name
            print(name)
            lh = list()
            for(h in 1:nhit){
                obj = bm[[i]]@aligns[[h]]
                lh[[as.character(h)]] = data.frame(motif=name,
                                                   tf = obj@TF@name,
                                                   evalue = obj@evalue,
                                                   sequence = obj@sequence,
                                                   match = obj@match, 
                                                   strand = obj@strand,
                                                   hit=h)
            }
            lbm[[name]] = do.call(rbind, lh)
            
        }
        dbm = do.call(rbind, lbm)    
        return(dbm)
	}
	
	