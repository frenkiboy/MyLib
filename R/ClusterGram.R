library(amap)
ks.default <- function(rows) seq(2, max(3, rows %/% 4))

many_kmeans <- function(x, ks = ks.default(nrow(x)), nproc=5, ...) {

	require(doMC)
	registerDoMC(nproc)
	l.clust = foreach(i = ks)%dopar%{
		print(i)
		k = Kmeans(x, i, iter.max=1000, method='correlation')
		data.frame(obs = seq_len(nrow(x)), i = i, k = ks[i], cluster = k$cluster)
	}
	do.call(rbind, l.clust)
}

center <- function(x) x - mean(range(x))

#' @param clusters data frame giving cluster assignments as produced by 
#'   many_kmeans or all_hclust
#' @param y value to plot on the y-axis.  Should be length
#'   \code{max(clusters$obs)}
clustergram <- function(clusters, y, line.width = NULL) {
  clusters$y <- y[clusters$obs]
  clusters$center <- ave(clusters$y, clusters$i, clusters$cluster)  

  if (is.null(line.width)) {
    line.width <- 0.5 * diff(range(clusters$center, na.rm = TRUE)) / 
      length(unique(clusters$obs))
  }
  clusters$line.width <- line.width
  
  # Adjust center positions so that they don't overlap  
  clusters <- clusters[with(clusters, order(i, center, y, obs)), ]  
  clusters <- ddply(clusters, c("i", "cluster"), transform, 
    adj = center + (line.width * center(seq_along(y)))
  )
  
  structure(clusters, 
    class = c("clustergram", class(clusters)),
    line.width = line.width)
}

plot.clustergram <- function(x) {
	require(ggplot2)
  i_pos <- !duplicated(x$i)
  
  means <- ddply(x, c("cluster", "i"), summarise, 
    min = min(adj), max = max(adj))
  
  ggplot(x, aes(i)) +
    geom_ribbon(aes(y = adj, group = obs, fill = y, ymin = adj - line.width/2, ymax = adj + line.width/2, colour = y)) + 
    geom_errorbar(aes(ymin = min, ymax = max), data = means, width = 0.1) + 
    scale_x_continuous("cluster", breaks = x$i[i_pos], labels = x$k[i_pos]) +
    labs(y = "Cluster average", colour = "Obs\nvalue", fill = "Obs\nvalue")
    
}

iris_s <- scale(iris[,-5]) 
k_def <- many_kmeans(iris_s)
k_10 <- many_kmeans(iris_s, 2:10)

pr <- princomp(iris_s)
pr1 <- predict(pr)[, 1]
pr2 <- predict(pr)[, 2]

plot(clustergram(k_def, pr1))
plot(clustergram(k_rep, pr1))
plot(clustergram(k_rep, pr2)