### calculates the correlation matrix and iteratively removes columns from the original matrix with corr greater than the limit value
corrRemover = function(b, lim)
{
	while(ncol(b) > 0){
		d = cor(b)
		g = table(which(d > lim, arr.ind = T))
		col = as.vector(which(g == max(g)))
		if ((sum(g) == min(g)*length(g))){
			break
		}
		b = b[,-col]
	}
	return(b)
}
	
