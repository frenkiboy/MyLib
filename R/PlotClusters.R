PlotClusters = function(m, k, outpath, name){
	kn = length(unique(k))
	cat('Plotting the parcoord...\n')
	library(reshape2)
	library(ggplot2)
	d = data.frame(cbind(m, cl=k, id.var=1:nrow(m)))
	longDF <- melt(d, id=c("cl", "id.var"))
	CairoPNG(file.path(outpath, paste(name,'k',kn,'kmeans.profiles.png', sep='.')), width=2500, height=2500)
		p <- ggplot(data = longDF, aes(x = variable, y = value, group = id.var)) + geom_line(cex=.1, col='darkred', alpha=.5) + stat_summary(aes(group = 1), fun.y='mean',geom='line', size=1.2, col=gray(.1))  + facet_wrap('cl', ncol=3) + geom_hline(aes(x = 0))
		p <- p + theme(axis.text.x = element_text(size=20, color = rainbow(5)), axis.text.y = element_text(size = 20))
		print(p)
	dev.off()
}