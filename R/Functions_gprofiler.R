# ------------------------------------------------------------------------ #
get_Profile = function(setl, 
                       background='', 
                       organism='mmusculus',
                       min_set_size=3,
                       max_set_size=3000,
                       min_insect_size=2,
                       correction_method='gSCS',
                       hier_filtering='moderate',
                       significant=FALSE){
    
    require(gProfileR)
    gpl = list()
    for(i in 1:length(setl)){
        
        sname = names(setl)[i]
        print(sname)
        set = setl[[sname]]
        gp = gprofiler(
            query             = set,
            organism          = organism,
            min_set_size      = min_set_size,
            max_set_size      = max_set_size,
            min_isect_size    = min_insect_size,
            correction_method = correction_method,
            hier_filtering    = hier_filtering,
            custom_bg=background,
            significant=significant
        )
        gpl[[sname]] = gp
        
    }
    return(gpl)
}

# ------------------------------------------------------------------------ #
select_Profile = function(gpl, pval=0.01){
    
    gpd = lapply(names(gpl), function(x){
        message(x)
        d=gpl[[x]][,c('term.name','domain','p.value')];
        if(nrow(d)==0)
            return(data.table(term.name=NA, domain=NA, p.value=0, set=x))
        d[,ncol(d)] = -log10(d[,ncol(d)])
        d$set=x
        data.table(d)
    })
    gdd = rbindlist(gpd)
    gdd$set = factor(gdd$set, levels=names(gpl))
    gdd = dcast(term.name+domain~set, value.var='p.value', fill=0, data=gdd,fun.aggregate=sum)
    gdd = data.table(gdd)
    gdd = gdd[rowSums(gdd[,-c(1:2),with=FALSE] > -log10(pval)) > 0,]
    gdd = na.omit(gdd)
    return(gdd)
}

# ------------------------------------------------------------------------ #
plot_Profile = function(dtab, name, path.out=NULL, width=12, height=12, col='red',
                        cwidth=2, cfont=10, cheight=6){
    
    require(ComplexHeatmap)
    tabl = split(dtab, dtab$domain)
    
    if(!is.null(path.out))
        pdf(file.path(path.out, DateNamer(paste(name, 'pdf',sep='.'))),
            width=width,
            height=height)
    
    for(i in 1:length(tabl)){
        name = names(tabl)[i]
        print(name)
        dat = tabl[[name]]
        if(nrow(dat)>1){
            m = as.matrix(dat[,-c(1:2),with=F])
            rownames(m) = dat$term.name
            h1=Heatmap(m, column_title=name, 
                       col=colorRampPalette(c(gray(.97),col))(50),
                       width=unit(cwidth,'cm'),
                       cluster_columns=FALSE,
                       column_names_gp = gpar(fontsize = cfont),
                       column_names_max_height = unit(cheight,'cm')
            )
            
            draw(h1, heatmap_legend_side='left')
        }
    }
    
    if(!is.null(path.out))
        dev.off()
}
