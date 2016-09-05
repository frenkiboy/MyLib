library(rtracklayer)
library(stringr)
library(data.table)
library(dplyr)

files = list.files('./', full.names=TRUE, pattern='bedGraph')

lcov = list()
for(i in 1:length(files)){
    
    file = files[i]
    name = basename(file)
    name = str_replace(name,'\\..+','')
    print(name)
    f = fread(file) %>% 
        filter(V4 > 0.25) %>% 
        mutate(V4 = round(V4,2)) %>% 
        as.data.frame()
    g = with(f, GRanges(V1, IRanges(V2+1,V3), weight=V4))
    cov = coverage(g, weight=as.integer(g$weight*100))
    lcov[[name]] = cov
}
lcov = Reduce(c, lcov)
saveRDS(lcov, 'hg19.phastCons46way.placental.RleList.rds')