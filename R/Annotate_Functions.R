
# ---------------------------------------------------------------------------- # annotates the ranges with the corresponding list
setGeneric("AnnotateRanges",
           function(region, annotation,
                    ignore.strand=FALSE,
                    type='precedence',
                    null.fact='None',
                    collapse.char=':',
                    precedence=NULL,
					id.col=NULL)
               standardGeneric("AnnotateRanges") )

setMethod("AnnotateRanges",signature("GRanges","GRangesList"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':'){

    if(! class(region) == 'GRanges')
        stop('Ranges to be annotated need to be GRanges')

    if(! all(sapply(annotation, class) == 'GRanges'))
        stop('Annotating ranges need to be GRanges')

    if(!type %in% c('precedence','all'))
        stop('type may only be precedence and all')

    require(data.table)
    require(GenomicRanges)
    cat('Overlapping...\n')
    if(any(names(is.null(annotation))))
        stop('All annotations need to have names')

    if(class(annotation) != 'GRangesList')
        annotation = GRangesList(lapply(annotation, function(x){values(x)=NULL;x}))

    a = suppressWarnings(data.table(as.matrix(findOverlaps(region, annotation, ignore.strand=ignore.strand))))
    a$id = names(annotation)[a$subjectHits]
    a$precedence = match(a$id,names(annotation))
    a = a[order(a$precedence)]

    if(type == 'precedence'){
        cat('precedence...\n')
        a = a[!duplicated(a$queryHits)]
        annot = rep(null.fact, length(region))
        annot[a$queryHits] = a$id
    }
    if(type == 'all'){
        cat('all...\n')
        a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
        annot = rep(null.fact, length(region))
        annot[a$queryHits] = a$id

    }
    return(annot)

})

setMethod("AnnotateRanges",signature("GRanges","list"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':'){

		  if(!all(unlist(lapply(annotation, 'class')) == 'GRanges'))
			stop('all elements of annotation need to be GRanges objects')

			annotation = GRangesList(lapply(annotation, function(x){values(x)=NULL;x}))
			AnnotateRanges(region, annotation, ignore.strand, type, null.fact, collapse.char)

		  })



setMethod("AnnotateRanges",signature("GRanges","GRanges"),
          function(region, annotation, ignore.strand=FALSE, type = 'precedence', null.fact = 'None',collapse.char=':', precedence=NULL, id.col=NULL){

				if(is.null(id.col))
					stop('Annotation needs to have a specified id')

				if(is.null(values(annotation)[[id.col]]))
					stop('id.col is not a valid column')


				if(!is.null(precedence)){
                if(!all(precedence %in% annotation$id))
                    stop('all precednce leveles have to be in the annotation id')
				}else{
					type='all'
					message('type set to all when precedence is not defined')
				}

				if(!type %in% c('precedence','all'))
					stop('type may only be precedence and all')

				require(data.table)
				require(GenomicRanges)
				cat('Overlapping...\n')
				if(any(names(is.null(annotation))))
					stop('All annotations need to have names')

				a = suppressWarnings(data.table(as.matrix(findOverlaps(region, annotation, ignore.strand=ignore.strand))))
				a$id = values(annotation)[[id.col]][a$subjectHits]


              if(type == 'precedence'){
                  cat('precedence...\n')
                  a$precedence = match(a$id,precedence)[a$subjectHits]
                  a = a[order(a$precedence)]
                  a = a[!duplicated(a$queryHits)]

              }
              if(type == 'all'){
                  cat('all...\n')
                  a = a[,list(id=paste(unique(id),collapse=collapse.char)),by='queryHits']
              }
			  annot = rep(null.fact, length(region))
              annot[a$queryHits] = a$id
              return(annot)

          })
