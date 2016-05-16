# genomic functions for overlapping features, mostly Galaxy-like functions
# by Altuna Akalin



# interval1 and interval2 are data.frames with chr,start,end and strand(optional) positions
.overlappingIntervals<-function(interval1,interval2)
  {
     if(dim(interval1)[1]==0 || dim(interval2[1])==0 ){return(NULL)}
     s1=as.numeric(as.vector(interval1[,2]))
     s2=as.numeric(interval2[,2])
     e1=as.numeric(as.vector(interval1[,3]))
     e2=as.numeric(interval2[,3])
     p1=outer(s1,e2,function(x,y) x<=y )
     p2=outer(e1,s2,function(x,y) x>=y )
     p=p1+p2
     indeces=which(p==2,arr.ind=T)
     if(dim(indeces)[1]==0){return(NULL)}
     return(cbind(interval1[indeces[,1],] , interval2[indeces[,2],] ))

  }

# interval1 and interval2 are data.frames with chr,start,end and strand(optional) positions
# first dataset shud be contained in the boundaries of the second one
#    S1---------E1
#  S2---------------E2
#  S1>=S2
#  E2>=E1
.containedIntervals<-function(interval1,interval2)
  {
     if(dim(interval1)[1]==0 || dim(interval2[1])==0 ){return(NULL)}
     s1=as.numeric(as.vector(interval1[,2]))
     s2=as.numeric(interval2[,2])
     e1=as.numeric(as.vector(interval1[,3]))
     e2=as.numeric(interval2[,3])
     p1=outer(s1,s2,function(x,y) x>=y )
     p2=outer(e1,e2,function(x,y) y>=x )
     p=p1+p2
     indeces=which(p==2,arr.ind=T)
     if(dim(indeces)[1]==0){return(NULL)}
     return(cbind(interval1[indeces[,1],] , interval2[indeces[,2],] ))

  }

# takes two bed files
# and outputs the rows of BED1 that are contained in the rows of BED2
featureContainBED.strand<-function(bed1,bed2)
{
    #get unique shared chromosomes between 2 files
    chrs=unique(as.vector(bed1[,1]))
    chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
    result=cbind(bed1[1,],bed2[1,],deparse.level=0)
    for (i in 1:length(chrs) )
    {
        overlap.plus=.containedIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="+",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="+",])

        overlap.minus=.containedIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="-",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="-",])
        
        result=rbind(result,overlap.plus,overlap.minus,deparse.level=0)
    }
    if(dim(result)[1]==1){return(NULL)}
    return(result[2:dim(result)[1],])
}




# takes two bed files
# and outputs the overlaping rows between them
featureOverlapBED<-function(bed1,bed2)
{
    #get unique shared chromosomes between 2 files
    chrs=as.vector(unique(bed1[,1]))
    chrs=chrs[chrs %in% unique(bed2[,1])]
    result=cbind(bed1[1,],bed2[1,],deparse.level=0)
    for (i in 1:length(chrs))
    {
        overlap=.overlappingIntervals(
        bed1[bed1[,1]==as.character(chrs[i]),],
        bed2[bed2[,1]==as.character(chrs[i]),])
	 
        result=rbind(result,overlap,deparse.level=0)
    }

    if(dim(result)[1]==1){return(NULL)}
    return(result[2:dim(result)[1],])
}

featureOverlapBED.para<-function(bed1,bed2)
{
    #get unique shared chromosomes between 2 files
    chrs=as.vector(unique(bed1[,1]))
    chrs=chrs[chrs %in% unique(bed2[,1])]
    result=NULL
    if( "rparallel" %in% names( getLoadedDLLs()) )
    {
      runParallel( resultVar="result", resultOp="rbind", nWorkers=3,verbose="info")
    }
    else
    {


      for (i in 1:length(chrs))
      {
		overlap=.overlappingIntervals(
        bed1[bed1[,1]==as.character(chrs[i]),],
        bed2[bed2[,1]==as.character(chrs[i]),])
        
        result=rbind(result,overlap)
      }

      #if(dim(result)[1]==1){return(NULL)}
      return(result)
    }
}

# takes two bed files
# and outputs the overlaping rows between them
featureOverlapBED.strand.para<-function(bed1,bed2)
{
    require(rparallel)
    #get unique shared chromosomes between 2 files
    chrs=unique(as.vector(bed1[,1]))
    chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
    result<-NULL
    if( "rparallel" %in% names( getLoadedDLLs()) )
    {
       runParallel( resultVar="result", resultOp="rbind", nWorkers=3,verbose="info")
    }
    else
    {
      for (i in 1:length(chrs) )
      {
        overlap.plus=.overlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="+",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="+",])

        overlap.minus=.overlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="-",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="-",])
        
        result=rbind(result,overlap.plus,overlap.minus,deparse.level=0)
      }
    }
    
    return(result)
}



# takes two bed files
# and outputs the overlaping rows between them
featureOverlapBED.strand<-function(bed1,bed2)
{
    #get unique shared chromosomes between 2 files
    chrs=unique(as.vector(bed1[,1]))
    chrs=chrs[chrs %in% unique(as.vector(bed2[,1]))]
    result=cbind(bed1[1,],bed2[1,],deparse.level=0)
    for (i in 1:length(chrs) )
    {
        overlap.plus=.overlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="+",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="+",])

        overlap.minus=.overlappingIntervals(
        bed1[bed1[,1]==chrs[i] & bed1[,6]=="-",],
        bed2[bed2[,1]==chrs[i] & bed2[,6]=="-",])
        
        result=rbind(result,overlap.plus,overlap.minus,deparse.level=0)
    }
    if(dim(result)[1]==1){return(NULL)}
    return(result[2:dim(result)[1],])
}


#bed1=data.frame(chr=c("chr1","chr1","chr1","chr1","chr1","chr1"),start=c(1,3,5,7,10,12),end=c(2,4,6,8,11,13))
#clusters the bed file in a very beatifule way

#bed1=bed1[order(bed1[,2]),]

#dist=outer(bed1[,2],bed1[,3],function(x,y) x-y)
#diag(dist)=rep(0,length(diag(dist)))
#z=hclust(as.dist((dist),upper=T,diag=T))

#diff(bed1[,2],bed1[,3],lag=1)





clusterFeaturesBED<-function(bed1,thold=1000,singleton.clusters=T)
{

    chrs=unique(bed[,1])
    result=bed[1,]
    for (i in 1:length(chrs))
    {
      bed1=bed[bed[,1]==chrs[i],]
      bed1=bed1[order(bed1[,2]),]
      dist=bed1[2:nrow(bed1),2]-bed1[1:nrow(bed1)-1,3]

      locations <-dist<=thold
      i <- c(FALSE,locations)
      j <- c(locations,FALSE)
      start <- which(j&!i)
      end <- which(i&!j)
      clusters=data.frame(chr=as.vector(bed1[start,1]),start=as.numeric(bed1[start,2]),end=as.numeric(bed1[end,3]))

      if(singleton.clusters==T){
        sing.ind=which(locations==FALSE)
        sing.ind=sing.ind[!sing.ind %in% start & !sing.ind %in% end]
 

        sings=bed1[sing.ind,]
        rbind(clusters,sings)
      }
    }

    return(result[2:dim(result)[1],])
}


convertToBED<-function(data,chr=1,start=2,end=3,score=NA,name=NA,strand=NA)
{



}


DGE.sig<-function(tag1,tag2,lib.size1,lib.size2)
{
  a=lib.size2/lib.size1
  b=tag1+tag2
  p.val=(a**tag2)*(factorial(b)/(factorial(tag1)*factorial(tag2)))/(((1+a)**(b+1))     )
  return(p.val)
}

DGE.fisher<-function(tag1,tag2,lib.size1,lib.size2)
{
  p.val=fisher.test(matrix(c(tag1,tag2,lib.size1-tag1,lib.size2-tag2),nrow=2))$p.value
  return(p.val)
}


DGE.sig.cum<-function(tag1,tag2,lib.size1,lib.size2)
{
   if(tag1==0 && tag2==0){return(1)}
   p.sum=0
    y=tag2
    while(y>=0)
      {
        p.itr=DGE.sig(tag1,y,lib.size1,lib.size2)
        
         
            p.sum=p.sum+DGE.sig(tag1,y,lib.size1,lib.size2)
         
        y=y-1
      }

   if(tag2/lib.size2 < tag1/lib.size1)
  {
    return(p.sum)
  }
  else
  {

     return(1-p.sum)
  }


}

### takes two very large bed files (eg. more than 15000x15000 rows per chr) and does the overlap
BigFeatureOverlapBED<-function(bed1,bed2)
{
    #get unique shared chromosomes between 2 files
	chrs=as.vector(unique(bed1[,1]))
    chrs=chrs[chrs %in% unique(bed2[,1])]
    result=cbind(bed1[1,],bed2[1,],deparse.level=0)
	#maximum size of dataframe before entering dataframe decomposition = 2*max.size
	max.size = 15000
    for (i in 1:length(chrs))
    {
        bed1.chr = bed1[bed1[,1]==as.character(chrs[i]),]
		bed2.chr = bed2[bed2[,1]==as.character(chrs[i]),]
		
		if ((nrow(bed1) + nrow(bed2)) > 2 * max.size)
		{
			bed1.chr = Designator(bed1.chr, max.size)
			bed2.chr = Designator(bed2.chr, max.size)
			tmp.result = data.frame(matrix(ncol = (ncol(bed1.chr) + ncol(bed2))))
			tmp.result = na.omit(tmp.result)
			for (j in 1:max(bed1.chr$designator))
			{
				for(k in 1:max(bed2.chr$designator))
				{
					overlap=.overlappingIntervals(
					bed1.chr[bed1.chr$designator == j,],
					bed2.chr[bed2.chr$designator == k,]
					)
					tmp.result = rbind(tmp.result,overlap,deparse.level=0)
				}
			}
			
			tmp.result = tmp.result[,-c(ncol(bed1.chr),ncol(tmp.result))]
			if(ncol(tmp.result) == ncol(result))
			{
				names(tmp.result) = names(result)
			}else
			{
				cat("differing number of columns in tmp.result and result\n")
			}
			
			result=rbind(result,tmp.result,deparse.level=0)
		
		}else
		{
			overlap=.overlappingIntervals(
			bed1.chr,
			bed2.chr
			)
			result=rbind(result,overlap,deparse.level=0)
		}
	}

    if(dim(result)[1]==1){return(NULL)}
    return(result[2:dim(result)[1],])
}


### takes two very large bed files (eg. more than 15000x15000 rows per chr) and does the overlap
### makes use of multiple processors
### SUCKS UP MEMMORY!!!
### let the first bed file be the bigger one
BigFeatureOverlapBED.para<-function(bed1,bed2)
{
	library(rparallel)
    #get unique shared chromosomes between 2 files
	chrs=as.vector(unique(bed1[,1]))
    chrs=chrs[chrs %in% unique(bed2[,1])]
    result=cbind(bed1[1,],bed2[1,],deparse.level=0)
	#maximum size of the dataframe before entering te dataframe decomposition = 2*max.size
	max.size = 15000
    for (i in 1:length(chrs))
    {
        bed1.chr = bed1[bed1[,1]==as.character(chrs[i]),]
		bed2.chr = bed2[bed2[,1]==as.character(chrs[i]),]
		
		if ((nrow(bed1) + nrow(bed2)) > 2 * max.size)
		{
			bed1.chr = Designator(bed1.chr, max.size)
			bed2.chr = Designator(bed2.chr, max.size)
			tmp.result = data.frame(matrix(ncol = (ncol(bed1.chr) + ncol(bed2))))
			tmp.result = na.omit(tmp.result)
			if( "rparallel" %in% names( getLoadedDLLs()) )
			{
				runParallel( resultVar="tmp.result", resultOp="rbind", nworkers = 3)
			}
			else 
			{     
				for (j in 1:max(bed1.chr$designator))
				{
					for(k in 1:max(bed2.chr$designator))
					{
						overlap=.overlappingIntervals(
						bed1.chr[bed1.chr$designator == j,],
						bed2.chr[bed2.chr$designator == k,]
						)
						tmp.result = rbind(tmp.result,overlap,deparse.level=0)
					}
				}
			}
			
			tmp.result = tmp.result[,-c(ncol(bed1.chr),ncol(tmp.result))]
			result=rbind(result,tmp.result,deparse.level=0)
		
		}else
		{
			overlap=.overlappingIntervals(
			bed1.chr,
			bed2.chr
			)
			result=rbind(result,overlap,deparse.level=0)
		}
	}

    if(dim(result)[1]==1){return(NULL)}
    return(result[2:dim(result)[1],])
}


### takes one large bed file and designates every max.size elements
Designator <- function(bed, max.size)
{
	num = floor(nrow(bed)/max.size)
	rest = nrow(bed) - (num*max.size)
	designator = rep(1:num, each=max.size)
	designator.rest = rep(num+1, each=rest)
	designator = c(designator, designator.rest)
	bed = cbind(bed, designator)
	
	return(bed)
}