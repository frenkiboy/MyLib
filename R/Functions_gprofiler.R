# ------------------------------------------------------------------------ #
get_Profile = function(setl, background='', organism='mmusculus'
                       min_set_size=3
                       max_set_size=3000,
                       min_insect_size=2,
                       correction_method='gSCS',
                       hier_filtering='moderate'){
    
    
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
            custom_bg=background
        )
        gpl[[sname]] = gp
        
    }
    return(gpl)
}

# ------------------------------------------------------------------------ #
select_Profile = function(gpl){
    
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
    gdd = dcast(term.name+domain~set, value.var='p.value', fill=0, data=gdd)
    gdd = na.omit(gdd)
    return(gdd)
}

# ------------------------------------------------------------------------ #
plot_Profile = function(dtab, name, path.out, width=12, height=12, col='red',
                        cwidth=2){
    
    require(ComplexHeatmap)
    tabl = split(dtab, dtab$domain)
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
                       cluster_columns=FALSE
            )
            
            draw(h1, heatmap_legend_side='left')
        }
    }
    dev.off()
}
