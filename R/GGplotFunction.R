ggplotScatter = function(d, x, y, diff=FALSE, plot.cor=FALSE){
	g=ggplot(d, aes(x=x,y=y)) + 
	 		 geom_point(color='black', cex=2)
	if(diff){
		if(!diff %in% colnames(d))
			stop(paste(diff, 'is not a designated colnames variable'))
	}		 
	
	g = g + xlab(x) + ylab(y) 
	g = g + annotate("text", x =min(5), y = max(15), label = paste('cor:', cort), size=10)
	g = g + labs(title=name)+
	        theme(plot.title = element_text(size = rel(2)), 
		 		  axis.text  = element_text(colour = "black", size=rel(2)),
				  axis.title = element_text(size = rel(2)),
				  axis.line  = element_line(size = rel(3)),
				  panel.background = element_rect(colour = "white", fill = "white"))
	return(g)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}