### INFO: Functions for drawing positional distributions around regions
### DATE: 29.07.2011.
### AUTHOR: Vedran Franke



# {1} FUNCTIONS

	# {{1}}
	DrawHeatmapClusters = function(l, distance = 'correlation', indices=1:length(l), outpath, set='file', width=150, height=800, colnum=50){
		
		expr.palette = colorRampPalette(c('gray', 'red'),bias=1, interpolate='spline')
		cpg.palette = colorRampPalette(c('darkblue', 'darkorange'))
		histone.palette=colorRampPalette(c('white', 'darkblue'),bias=1, interpolate='spline')

		library(flashClust)
		for(i in indices){
			m=l[[i]]
			name = names(l)[i]
			print(name)
			if(distance == 'correlation'){
				# d = as.dist(1-cor(t(m)+runif(length(m), min=0, max=0.01)))
				d = as.dist(1-cor(t(m + runif(length(m), 0, min(m[m != 0])))))
			}else if(distance == 'euclidean'){
				d = dist(m)
			}
			h = hclust(d)
			o = h$order
				
			CairoPNG(file.path(outpath, paste(set, name, distance, 'image', 'png', sep='.')), width=width*length(l), height=height)
				layout(matrix(c(1:(length(l))), nrow=1))
				par(oma=c(1,1,1,1), mar=c(2,.3,4,.1), cex.axis=1.5, cex.main = width/125, yaxt='n', xaxt='n')
				for(j in seq(along=l)){
					print(j)
					plot.name = names(l)[j]
					m1 = l[[j]]
					if(j == 1){
						m1 = m1+min(m1[m1 != 0])
						image(log10(t(m1[o,])), col=expr.palette(colnum), main='expression', xlab='-1000:1000')
					}else if(j %in% c(2,3)){
						image(t(matrix(m1[o], ncol=1)), col=cpg.palette(2), main='cpg.prom', xlab='-1000:1000')
						
					}else{
						image(t(m1[o,]), col=histone.palette(colnum), mar=c(0,0,0,0), main=plot.name, xlab='-1000:1000')
						abline(v=.5, col='darkorange')
					}
						
				}
			dev.off()
		}
	}
	#/{{1}}
#/{1} FUNCTIONS
	




