# ------------------------------------------------------------------------ #
.gene.set.libraries = c(
    'GO_Molecular_Function_2015',
    'GO_Biological_Process_2015',
    'GO_Cellular_Component_2015',
    'KEGG_2016','Panther_2016',
    'Reactome_2016',
    'WikiPathways_2016',
    'KEGG_2016',
    'OMIM_Disease',
    'GeneSigDB',
    'MSigDB_Oncogenic_Signatures',
    'MGI_Mammalian_Phenotype_Level_4',
    'Single_Gene_Perturbations_from_GEO_up',
    'Single_Gene_Perturbations_from_GEO_down',
    'Virus_Perturbations_from_GEO_up',
    'Virus_Perturbations_from_GEO_down',
    'Drug_Perturbations_from_GEO_up',
    'Drug_Perturbations_from_GEO_down',
    'Disease_Perturbations_from_GEO_up',
    'Disease_Perturbations_from_GEO_down',
    'TargetScan_microRNA',
    'Epigenomics_Roadmap_HM_ChIP-seq',
    'ENCODE_TF_ChIP-seq_2015',
    'CORUM',
    'ESCAPE',
    'Tissue_Protein_Expression_from_ProteomicsDB',
    'TF-LOF_Expression_from_GEO',
    'Panther_2016',
    'Genome_Browser_PWMs')

get_Enrichr_list = function(gene.lists){
    
    library(httr)
    library(xml2)
    if(is.character(gene.lists))
        gene.lists = list(gene_list=gene.lists)
    
    lres = list()
    for(name in names(gene.lists)){
        message(name)    
        lres[[name]] = list()
        
        for(i in .gene.set.libraries){
            message(i)
            lres[[name]][[i]] = get_Enrichr(gene.lists[[name]], i)
        }
    }
    return(lres)
}
    
# ------------------------------------------------------------------------ #
combine_Enrichr = function(lres){
    
    library(data.table)
    lmat=list()
    for(i in .gene.set.libraries){
        l = lapply(names(lres), function(x){
            d = lres[[x]][[i]]
            if(nrow(d)> 0){
                d$set = x
                data.table(d)
            }
        })
        d = rbindlist(l)
        if(nrow(d)>0){
            d$log = -log10(d$'Adjusted P-value')
            d$log[is.infinite(d$log)] = 0
            dmat = dcast(d, Term~set, value.var = 'log', fill=0)
            dmat = data.table(domain = i, dmat)
            lmat[[i]] = dmat
        }
    }
    return(rbindlist(lmat))
}

# ------------------------------------------------------------------------ #
select_Enrichr = function(mat, pval=0.05){
    
    if(nrow(mat) > 0){
        rind = rowSums(mat[,-(1:2),with=FALSE] > -log10(pval)) > 0
        if(sum(rind) > 0){
            mat = mat[rind,]
        }else{
            mat = NULL
        }
    }
    
    return(mat)
}

# ------------------------------------------------------------------------ #
get_Enrichr = function(gene.list = NULL,
                       gene.set.library = NULL,
                       description = 'gene_list'
                       ){
   
    if(is.null(gene.list))
        stop('Please specify the gene list')
    
    if(is.null(gene.set.library))
        stop('Please specify the library')
    
    library(httr)
    library(xml2)
    enrichr='http://amp.pharm.mssm.edu/Enrichr/'
    
    # ------------------------------------------------------------------------ #
    # POST
    p = POST(file.path(enrichr, 'addList'), 
             body=list(list=paste(gene.list, collapse='\n'), 
                       description=description))
    user.list.id = sub('\\n.+','',unlist(strsplit(unlist(as_list(content(p))),' ')))[3]
    
    # ------------------------------------------------------------------------ #
    # GET
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'

    g = GET(paste0(file.path(enrichr, 'export'),  
                   sprintf(query_string, user.list.id, 'file.txt', gene.set.library)))
    
    tab = lapply(strsplit(intToUtf8(content(g)),'\n'), strsplit, '\t')
    tab = do.call(rbind, lapply(tab[[1]], function(x)as.data.frame(t(x))))
    colnames(tab) = tab[1,];tab=tab[-1,]
    tab[,c(3,4,5,6,7,8)] = round(as.data.frame(lapply(tab[,c(3,4,5,6,7,8)], as.numeric)),3)
    return(tab)
}

# ------------------------------------------------------------------------ #
plot_Enrichr = function(tabl, outname, path.out=NULL, width=12, height=12, col='red',
                        cwidth=2, cfont=10, cheight=6){
        library(ComplexHeatmap)
        if(!is.null(path.out))
            pdf(file.path(path.out, DateNamer(paste(outname, 'pdf',sep='.'))),
                width=width,
                height=height)
    
        if(any(class(tabl) %in% 'data.table'))
            tabl = split(tabl, tabl$domain)
        
        for(i in 1:length(tabl)){
            name = names(tabl)[i]
            print(name)
            dat = tabl[[name]]
            if(nrow(dat)>1){
                m = as.matrix(dat[,-c(1:2),with=F])
                rownames(m) = dat$Term
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
    

# ------------------------------------------------------------------------ #
#' Enrichr_all
#'
#' @param glist a list of gene names corresponding to different conditions
#' @param pval p value for finding significantly enriched categories
#'
#' @return a list containing the original gene list, enriched categories 
#' per element of the gene list, combined data table, and selected categories
Enrichr_all = function(glist, pval=0.05){
    
    enrichr = get_Enrichr_list(glist)
    combine = combine_Enrichr(enrichr)
    ensel   = select_Enrichr(combine, pval)
    
    lout = list(glist = glist, enrichr=enrichr, combine=combine, ensel=ensel)
    return(lout)
}

# ------------------------------------------------------------------------ #
test_Enrichr = function(){
    gene.list = c("PHF14","RBM3","MSL1","PHF21A","ARL10","INSR","JADE2","P2RX7","LINC00662","CCDC101","PPM1B","KANSL1L","CRYZL1","ANAPC16","TMCC1","CDH8","RBM11","CNPY2","HSPA1L","CUL2","PLBD2","LARP7","TECPR2","ZNF302","CUX1","MOB2","CYTH2","SEC22C","EIF4E3","ROBO2","ADAMTS9-AS2","CXXC1","LINC01314","ATF7","ATP5F1")
    lres = get_Enrichr_list(gene.list) 
    
    l = list(list1=lres[[1]], list2=lres[[1]])
    comb = combine_Enrichr(l)
    sel = select_Enrichr(comb)
}

