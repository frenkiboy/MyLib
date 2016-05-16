### gets the rna folds using rnafold algorithm
	GetFolds = function(s, outpath, setname){
	
	
		rna.fold = '/common/USERS/vfranke/bin/ViennaRNA-1.8.5/RNAfold/bin/RNAfold'
		cat('Creating the output directory...\n')
		outdir = file.path(outpath, setname)
			dir.create(outdir, showWarnings=F)
		setwd(outdir)
		
		l.fold = list()
		cat('Folding the sequences...\n')
		for(i in 1:length(s)){
		
			name = names(s)[i]
			cat(name,'\r')
			se = gsub('T','U',as.character(s[i]))
			command = paste(rna.fold)
			fold = system(command = command, intern=T, input = se)
			command = paste('mv', file.path(outdir,'rna.ps') ,file.path(outdir,paste(i, name, setname, 'ps', sep='.')))
			system(command, intern=F)
			
			sp = unlist(strsplit(fold, ' '))
			l.fold[[name]]$seq = sp[1]
			l.fold[[name]]$fold = sp[2]
			l.fold[[name]]$gfe = gsub('\\)|\\(','',sp[3])
		}
		l.fold
	}