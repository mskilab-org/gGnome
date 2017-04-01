## DATA STRUCTURE CHALLENGE:
## Do we need a dynamic "mapping" between indexing of segs and the node ID in
## iGraph such that changes/modifications in either side would be propgrapagated,
## e.g. if you change nodes/edges, changes segs; change segs, change nodes/edges
## Challenge between iGraph and GRanges
## shortest path: easy in iGraph, hard to map back to GRanges
## genomic operations, e.g. merge two nodes: easy in GRanges `reduce()`,
## but how do you merge two nodes and keep their identity in iGraph?

## DATA STRUCTURE CHALLENGE: Do we need a dynamic "mapping" between indexing of segs and the node ID in
## iGraph such that changes/modifications in either side would be propgrapagated, e.g. if you change nodes/edges, changes segs; change segs, change nodes/edges
## Challenge between iGraph and GRanges
## shortest path: easy in iGraph, hard to map back to GRanges
## genomic operations, e.g. merge two nodes: easy in GRanges `reduce()`, but how do you merge two nodes and keep their identity in iGraph?

## preset environmental variables
## TODO: I want to print a message of which reference genome is used and how to change it,
## when the package is loaded
Sys.setenv(REF_GENOME = "data/hg19.broad.BSgenome.rds")   # created from BSgenome.Hsapiens.UCSC.hg19 in `batch.R`, no chr
GENOME = readRDS(Sys.getenv("REF_GENOME"))

## Question: what should we do about mtDNA and unmapped (?) U* regions?

## TODO: put some gTrack pars in environment



#' Junctions: R6 class extended from GRangesList to represent aberrant genomic SVs
#' with respect to a reference genome
#'
#' @import GenomicRanges
#' @import R6
#' @importFrom R6 R6Class
#' @export
#'
#' TODO: overload 'c', '[<-', '%Q%' operators!!!
junctions = R6Class("junctions",
                    public = list(
                        refG = "GENOME",
                        initialize = function(grl=NULL, raFile=NULL){
                            if (!is.null(grl) & length(grl)>0 & is(grl, "GRangesList")){
                                ## if non-empty GRangesList given, test its validity and setup
                                if (any(elementNROWS(grl)!=2)){
                                    stop("Input GRL must be length 2 for every element.")
                                }
                                ## if any width > 2, stop
                                if (any(any(width(grl)))>2){
                                    stop("Ambiguous breakpoint.")
                                }
                                ## ## if any width < 2, resize
                                ## if (any(any(width(grl)))<2){
                                ##     warning("At least 1 breakpoint defined on single nt. Assume it is left.")
                                ##     tmpMcols = mcols(grl)
                                ##     grl = GRangesList( lapply(grl, function(gr) return(resize(gr, width=2))) )
                                ##     mcols(grl) = tmpMcols
                                ## }
                                private$juncGrl = grl
                                return(self)
                            } else if (!is.null(raFile)){
                                ## if grl not given, raFile given, read the file
                                if (!file.exists(raFile)){
                                    stop("Couldn't find nput junction file.")
                                }

                                grl = ra_breaks(raFile, seqlengths=seqlengths(get(self$refG)))
                                ## if any width > 2, stop
                                if (any(any(width(grl)))>2){
                                    stop("Ambiguous breakpoint.")
                                }
                                ## if any width < 2, resize
                                if (any(any(width(grl)))<2){
                                    warning("At least 1 breakpoint defined on single nt. Assume it is left.")
                                    tmpMcols = mcols(grl)
                                    grl = GRangesList( lapply(grl, function(gr) return(resize(gr, width=2))) )
                                    mcols(grl) = tmpMcols
                                }
                                private$juncGrl = grl
                                return(self)
                            } else {
                                return(self)
                            }
                        },

                        ## subset, concatenation
                        subset = function(idx=NULL){
                            ## test if idx is NULL, logical, or numeric
                            if (length(idx)==0){
                                ## idx empty
                                return(junctions$new())
                            } else if (is.numeric(idx)){
                                ## idx all num
                                idx = as.integer(idx)
                                ## out of range if contains 0 or larger than max len
                                if (0 %in% idx | max(abs(idx))>length(private$juncGrl)){
                                    stop("Subscript out of bound.")
                                } else if (any(idx<0)){
                                    ## if any negative, must be all negative
                                    if (any(idx>0)){
                                        stop("Cannot mix negative subscripts with positive.")
                                    }
                                    newIdx = setdiff(seq_along(private$juncGrl), -idx)
                                    return(self$subset(newIdx))
                                } else {
                                    ## idx all num and positive
                                    newGrl = private$juncGrl[idx]
                                    return(junctions$new(newGrl))
                                }
                            } else if (is.logical(idx)){
                                return(self$subset(which(idx)))
                            } else {
                                stop("Invalid input.")
                            }
                        },
                        append = function(newJunc=NULL){
                            if (is(newJunc, "GRangesList")){
                                private$juncGrl = append(private$juncGrl, newJunc)
                                return(self)
                            } else if (is(newJunc, "junctions")){
                                self$append(newJunc$grl)
                                return(self)
                            } else {
                                stop("Invalid input.")
                            }
                        }
                    ),
                    private = list(
                        juncGrl = GRangesList()
                    ),
                    active = list(
                        ## getter for the grl
                        grl = function(){
                            return(private$juncGrl)
                        },

                        ## given a character vector of col names return that subset of mcols
                        values = function(cols="."){
                            ## input must be character
                            if (!is.character(cols)){
                                stop("Colnames must be character.")
                            }

                            if (cols == "."){
                                ## by default return everything
                                return(mcols(private$juncGrl))
                            } else {
                                ## otherwise return whatever we have
                                validCols = intersect(cols, colnames(mcols(private$jAnnotation)))
                                if (length(validCol) < length(cols)) {
                                    ## if we don't have that, warn and give up
                                    warning(paste("Omitted columns:", setdiff(cols, validCols)))
                                }
                                return(mcols(private$juncGrl)[, validCols, drop=F])
                            }
                        },
##                        mcols = values, ## can I do this??

                        length = function(){
                            return(length(private$juncGrl))
                        }
                    ))
## overload some S3 methods
#' `[` -- subsetting a junctions object
#'
#' @param numeric vector of indices
#'
`[.junctions` <- function(x, idx=NULL){
    return(x$subset(idx))
}

#' c -- concate junctions objects
#'
#'
#'
c.junctions <- function(...){
    ## TODO: think about the fastest way to implement `c`
}

#' @import gTrack
#' @import GenomicRanges
#' @import BSgenome
#' @import igraph
#' @import Matrix
#' @import parallel
#' @import data.table
#' @import gUtils
#' @import R6
#' @importFrom R6 R6Class
#' @export
gGraph = R6Class("gGraph",
                 public = list(
                     ## public fields
                     ## name = NULL,
                     refG = "GENOME", ## seqinfo of ref genome

                     ## constructor
                     initialize = function(tile=NULL, junctions=NULL, jabba=NULL, weaver=NULL,
                                           segs=NULL, es=NULL, ploidy=NULL, purity=NULL){
                         ## control how to construct
                         if (!is.null(segs) &
                                    !is.null(es) &
                                    !is.null(junctions) &
                                    !is.null(ploidy) &
                                    !is.null(purity)) {
                             private$gGraphFromScratch(segs, es, junctions, ploidy, purity)
                         } else if (!is.null(tile) | !is.null(junctions)) {
                             message("Initializing with 'tile' and 'junctions'")
                             self$karyograph(tile, junctions)
                         } else if (!is.null(jabba)) {
                             message("only use 'jabba' or 'weaver' field, not both")
                             self$jabba2gGraph(jabba)
                         } else if (!is.null(weaver)) {
                             message("only use 'jabba' or 'weaver' field, not both")
                             self$weaver2gGraph(weaver)
                         } else {
                             ## generate null graph
                             self$nullGGraph()
                         }
                     },

                     ## initialize from global ref genome seqinfo
                     nullGGraph = function(regular=TRUE){
                         tryCatch({
                             tmp = si2gr(get(self$refG)) %Q% (order(strand, seqnames, start))
                             names(tmp) = seq_along(tmp)
                         },
                         error = function(e){
                             cat(conditionMessage(e))
                             cat('\n')
                             cat('Reference genome not set or not Granges or BSgenome object.\n')
                         })
                         if (regular){
                             regularChr = c(as.character(1:22), "X", "Y")
                             tmp = tmp %Q% (seqnames %in% regularChr)
                         }
                         private$segs = c(tmp, gr.flipstrand(tmp)) ## null segs are ref
                         private$segs$cn = private$ploidy ## null cn is ploidy
                         private$segs$loose = FALSE ## all non-loose end

                         private$g = make_empty_graph(n=length(private$segs), directed=T)
                         private$junction = junctions$new()

                        ## close the circle on mitochondrial DNA
                         sinfo = as.data.frame(seqinfo(get(self$refG)))
                         sinfo = data.table(seqnames=rownames(sinfo), sinfo)
                         circChr = sinfo[isCircular==T, seqnames]
                         cat(paste('There is', length(circChr), 'circular contig(s): '))
                         cat(circChr, '\n')
                         ## QUESTION: why Rle doesn't allow match function???
                         circIx = which(as.vector(seqnames(private$segs)) %in% circChr)
                         if ( length(circIx)>0 ){
                             ## constructing edges: 5 required columns
                             ## from, to, cn, type, weight (len o' source node)
                             private$es = rbindlist(
                                 list(private$es, list(
                                                      from=circIx,
                                                      to=circIx,
                                                      cn=private$segs$cn[circIx],
                                                      type=rep("ref", length(circIx)),
                                                      weight=width(private$segs[circIx])
                                                  )
                                      )
                             )
                             ## Q: constructing an igraph easier than modifying?
                             private$g = add_edges(private$g,
                                                   t(as.matrix(private$es[,.(from,to)])),
                                                   ## from-to-from-to-...
                                                   attr = as.list(private$es))
                         }
                         ## assign terminals
                         ## DONE: init terminal field
                         whichTerminal = private$es[, union(
                                                     setdiff(seq_along(private$segs),
                                                             union(from, to)),setxor(from, to))]

                         return(self)
                     },

                     ## initialize from segmenatation AND/OR rearrangement junctions
                     addJuncs = function(junc){
                         ## DONE: populate abEdges while adding new junction!!!!
                         ## NOTE: the bps in junc must be width 2
                         ## TODO: what if junctions come with a CN?
                         ## TODO: bar GR input
                         ## ALERT: convention of junction orientation!!!
                         "Given a GRL of junctions add them to this gGraph."
                         ## 1. every single junction has 2 breakpoints,
                         ## break nodes in graph by these breakpoints
                         ## 2. based on oreintation of the junctions,
                         ## connect those nodes; introduce corresonding edges to graph
                         ## DONE: check if every bp within the ref genome
                         ## if not we need to resolve, maybe by creating new seqnames with warning
                         if (!is(junc, "junctions")){
                             ## NOTE: for a GRL to be junctions class,
                             ## must be 1) each element length 2 and with strand
                             ## 2) width 2, if not, convert
                             junc = junctions$new(junc)
                         }

                         junctions = junc$grl
                         ## save the junctions in the object
                         private$junction$append(junc)

                         if (length(junctions)==0){
                             return(self)
                         }

                         ## from this point only deal with overlapping junctions
                         jIn = grl.in(junctions, private$segs)
                         if (any(jIn==F)){
                             warning(
                                 paste(
                                     sum(jIn==F),
                                     "junctions not overlapping with any segment."))
                         }
                         junctions = junctions[jIn]

                         ## resize to width 2
                         jUl = unlist(junctions)
                         jUl = resize(jUl, 2,
                                      fix=ifelse(strand(jUl)=="+",
                                                 "start", "end"))
                         ## DONE: remember using split()!!
                         mc = mcols(junctions)
                         junctions = split(jUl, rep(seq_along(junctions), each=2))
                         mcols(junctions) = mc

                         ## start processing
                         ## TO DO: write as JaBbA::karyograph() with modifications
                         ## e.g. (30, 2) --> pivot (2, 30)
                         bp.p = grl.pivot(junctions)
                         bp.p = gr.fix(bp.p, get(self$refG))
                         juncTile = c(bp.p[[1]], bp.p[[2]])
                         ## BP 1 and 2, retaining strand-orientation info
                         ## bp1 = gr.start(bp.p[[1]], ignore.strand = T)
                         ## bp2 = gr.start(bp.p[[2]], ignore.strand = T)
                         ## tmpBps = c(bp1, bp2)
                         self$addSegs(juncTile) ## TODO: addSegs is giving dup edges!!!

                         ## now convert bp1 and bp2 to data.table
                         ## bp1 --> every bp associated fwith 4 nodes:
                         ## left +, left -, right +, right -
                         end1 = as.data.table(
                             values(gr.findoverlaps(private$segs, bp.p[[1]]))
                         )
                         ## orientation of the junction at bp1
                         end1[, ":="(ori=as.vector(strand(bp.p[[1]][subject.id])))]
                         end1 = cbind(gr2dt(private$segs[end1$query.id]), end1)
                         ## smaller coor on ref is called "left"
                         end1[, side := ifelse(start==min(start), "left", "right"),
                              by=subject.id]

                         ## bp2, similarly
                         end2 = as.data.table(
                             values(gr.findoverlaps(private$segs, bp.p[[2]]))
                         )
                         end2[, ":="(ori=as.vector(strand(bp.p[[2]][subject.id])))]
                         end2 = cbind(gr2dt(private$segs[end2$query.id]), end2)
                         end2[, side := ifelse(start==min(start), "left", "right"),
                              by=subject.id]

                         ## locate aberrant edges: 1 junction --> 2 edges
                         abEs =
                             Reduce("merge",
                                    list(
                                        end1[which((ori=="-" & side=="left" & strand=="+") |
                                                   (ori=="+" & side=="right" & strand=="-")),
                                             .(from1 = query.id), by=subject.id],
                                        ## node at bp1 sending out the edge
                                        end1[which((ori=="-" & side=="left" & strand=="-") |
                                                   (ori=="+" & side=="right" & strand=="+")),
                                             .(to1 = query.id), by=subject.id],
                                        ## node at bp1 receiving the edge
                                        end2[which((ori=="-" & side=="left" & strand=="+") |
                                                   (ori=="+" & side=="right" & strand=="-")),
                                             .(from2 = query.id), by=subject.id],
                                        ## node at bp2 sending out the edge
                                        end2[which((ori=="-" & side=="left" & strand=="-") |
                                                   (ori=="+" & side=="right" & strand=="+")),
                                             .(to2 = query.id), by=subject.id]))
                         ## node at bp2 receiving the edge
                         ## final edges: from1 --> to2, from2 --> to1
                         ## NOTE: ab.edges set to CN 1 at this point,
                         ## ignoring junction balance!!!

                         mP = as.matrix(abEs[, .(from=from1, to=to2, edge.ix=.I)])
                         mP = rbind(private$abEdges[,,"+"], mP)
                         mN = as.matrix(abEs[, .(from=from2, to=to1, edge.ix=.I)])
                         mN = rbind(private$abEdges[,,"-"], mN)
                         private$abEdges = array(c(mP, mN),
                                                 dim=c(nrow(mP), 3, 2),
                                                 dimnames=list(NULL,
                                                               c("from", "to", "edge.ix"),
                                                               c("+","-")))

                         abEs = abEs[,.(from=c(from1, from2), to=c(to2, to1),
                                        cn=1, type="aberrant"), by=subject.id][
                           , weight := width(private$segs[from])]

                         private$es = rbind(private$es, abEs[, -1])
                         ## connecting the junctions
                         private$g = add_edges(graph = private$g,
                                               edges = as.vector(t(as.matrix(abEs[, .(from, to)]))),
                                               attr = as.list(abEs[,.(junctionId=which(jIn)[subject.id], cn, type)]))
                         return(self)
                     },

                     addSegs = function(tile){
                         ## Given a GRanges obj of a segmentation (complete or not),
                         ## break the gGraph at their ends.
                         ## extract breakpoints
                         bps = reduce(c(gr.start(tile), gr.end(tile)))

                         ## break it
                         private$makeSegs(bps)

                         ## DONE: back tracing old node, connect its incoming edge to
                         ## the first fragment in the new nodes, its outgoing edge to
                         ## the last fragment. Between consecutive new fragments, introduce
                         ## reference edges. (Don't worry about junction balance yet!)
                         tmpDt = gr2dt(private$segs)
                         ## map the 5' end seg in new segs for every old seg,
                         ## they will receive old segs' incoming edges
                         tmpDt[, isHead :=
                                     (start == min(start) & strand=="+") |
                                     (end == max(end) & strand=="-"), by=qid]
                         ## map the 3' end seg in new segs for every old seg,
                         ## they will send out old segs' outgoing edges
                         tmpDt[, isTail :=
                                     (start == min(start) & strand=="-") |
                                     (end == max(end) & strand=="+"), by=qid]
                         ## enumerate all edges in es:
                         private$es[, .(from = tmpDt[, which(qid %in% from & isTail==T)],
                                        to = tmpDt[, which(qid %in% to & isHead==T)],
                                        cn, type, weight)] -> newEs
                         ## introduce ref edges between new breakpoints
                         refEs = tmpDt[, .(from=.I[isTail==F], to=.I[isHead==F]), by=qid]
                         refEs[, ":="(cn=tmpDt[from, cn], type="ref")]
                         refEs = refEs[!is.na(from) & !is.na(to),-1,with=F]
                         refEs[, weight := width(private$segs[from])]
                         newEs = rbindlist(list(newEs, refEs)) ## combine the two parts
                         newEs[!duplicated(newEs)]

                         ## update: es, g
                         private$es = newEs
                         private$g = graph_from_data_frame(
                             newEs, directed = TRUE,
                             vertices=data.frame(name = as.numeric(rownames(tmpDt)),
                                                 cn = tmpDt[, cn]))
                         ## reset tmpSegs
                         private$tmpSegs = NULL

                         return(self)
                     },
                     ## karograph: initialize `nullGGraph()`,
                     ## add junctions to it, then add tiles to it
                     karyograph = function(tile=NULL, junctions=NULL){
                         self$nullGGraph()
                         ## if there is tile, add tile
                         if (!is.null(tile) & length(tile)>0){
                             self$addSegs(tile)
                         }
                         ## if there is junctions, add junctions
                         if (!is.null(junctions) & length(junctions)>0){
                             ## if empty, ignore these GRanges lists
                             self$addJuncs(junctions)
                         }
                         return(self)
                     },
                     ## initialize from JaBbA output
                     jabba2gGraph = function(jabba){
                         ## make sure required mcol is filled
                         private$segs = gr.fix(jabba$segstats, get(self$refG))
                         tmpDt = gr2dt(private$segs)
                         tmpDt[, first:=min(start), by=seqnames]
                         tmpDt[, `:=`(first=min(start),last=max(end)), by=seqnames]
                         ## terminal defined as left/right most seg of a chr
                         ## or loose end decoy segs

                         ## DONE: redefine terminal, node wo both of in/out edges
                         private$segs$terminal = tmpDt[, (loose | start==first | end == last)]

                         private$es = as.data.table(jabba$edges[,1:4])
                         private$es[, weight := width(private$segs[from])]
                         private$g = make_directed_graph(
                             t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))

                         ## DONE: get the union of node ix wo in edge and out edge
                         ## SLOW!!!!
                         whichTerminal = private$es[, setxor(from, to)]


                         private$segs$terminal = seq_along(private$segs) %in% whichTerminal

                         private$junction = junctions$new(jabba$junctions)

                         private$abEdges = jabba$ab.edges
                         private$ploidy = jabba$ploidy
                         private$purity = jabba$purity
                         return(self)
                     },
                     ## initialize from Weaver result
                     weaver2gGraph = function(weaver){

                     },

                     ## public methods
                     ## I/O
                     print = function(){
                         cat('A gGraph object.\n')
                         cat('Based on reference genome: ')
                         cat(self$refG)
                         cat('\n\n')
                         cat('Total non-loose segmentation:')
                         cat(length(private$segs %Q% (loose==F & strand=="+")))
                         cat('\n\n')
                         cat('Junction counts:\n')
                         print(private$es[, table(type)/2])
                     },
                     plot = function(){},
                     summary = function(){
                         summ = list()
                         summ$v = length(private$segs)
                         summ$e = nrow(private$es)
                         return(summ)
                     },
                     length = function(){
                         ## DONE
                         if (is.null(private$partition)){
                             private$partition = self$components()
                         }
                         return(length(private$parition))
                     },
                     ##

                     gGraph2gTrack = function(){
                         "Create gTrack for static genome browser-style viz."
                         ## DONE: replicate classic JaBbA viz
                         ## plotting segments
                         ## if loose, make it white, lift it up
                         ss = private$segs
                         ed = private$es

                         ## set edge apperances
                         ## lwd, lty, col, cex.arrow, v, not.flat, h, dangle.w
                         ed[, ":="(lwd = ifelse(type=="aberrant", log2(0.2*cn+2)+1, 1),
                                   lty = ifelse(type=='loose', 3, 1),
                                   col = ifelse(type=="aberrant",
                                         ifelse(cn>0,
                                                alpha("red", 0.4),
                                                alpha("purple", 0.3)),
                                         ifelse(type=="loose",
                                                alpha("blue",0.6),
                                                alpha("grey",0.2))),
                                   cex.arrow = 0,
                                   not.flat = type=="aberrant",
                                   v = ifelse(type=="aberrant", 2, 1),
                                   h = ifelse(type=="aberrant", 2, 1),
                                   dangle.w = 0.5)]
                         ed = ed[!is.na(from) & !is.na(to)]## TODO: handle so/si-less edges, omit for now

                         ## set segment apperances
                         ## if loose, change its cn to slightly higher than it's incident node
                         if (any(ss$loose==T)){
                             lid = which(ss$loose)
                             ## find partner indices for loose ends
                             pid = sapply(lid,
                                          function(i) ed[from==i | to==i,
                                                         ifelse(from==i, to, from)])
                             ss$cn[lid] = ss$cn[pid]*1.2
                         }
                         ## col, border, ywid
                         ss$col = ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5))
                         ss$border = ifelse(ss$loose, ss$col, alpha("black", 0.5))
                         ss$ywid = ifelse(ss$loose, 0.001, 0.8)

                         gt = gTrack(ss, y.field="cn", edges=ed, name="gGraph", angle=0)
                         return(gt)
                     },
                     gGraph2iGraphViz = function(){
                         "Create igraph layout for static viz of graph structure."
                         ## TODO: channel data to igraph layout
                     },
                     gGraph2json = function(file=NULL,
                                            maxcn=100,
                                            maxweight=100,
                                            trim = TRUE ## trim will only output seqnames that are relevant to the plot
                                            ){
                         system(paste('mkdir -p', file))
                         system(sprintf('cp -r %s %s',
                                        paste0(system.file("extdata", "gTrack.js", package = 'gGnome'), '/*'),
                                        paste0(file, '/')))
                         "Create json file for interactive visualization."
                         qw = function(x) paste0('"', x, '"') ## quote

                         ## range of CN
                         ymin=0
                         ymax=maxcn

                         ## processing nodes
                         ## reduce strand
                         ## remove loose nodes
                         oid = gr2dt(private$segs)[, which(strand == "+" & loose==F & !is.na(cn))]
                         ## ori ind of rev comps
                         rid = gr2dt(private$segs)[, which(strand == "-" & loose==F & !is.na(cn))]
                         nodes = private$segs[oid]
                         ## ori ix of loose nodes
                         loose.id = which(private$segs$loose==T)

                         ## binding into dt
                         node.dt = data.table(
                             ## each row is a non-loose positive segment
                             oid,
                             rid,
                             iid = seq_along(nodes),
                             chromosome = qw(as.character(seqnames(nodes))),
                             startPoint = as.character(start(nodes)), ## smaller coor side
                             strand = "*",
                             endPoint = as.character(end(nodes)),
                                        #                             title = as.character(seq_along(nodes)),
                             title = paste0(seq_along(nodes), ' (', oid, '|', rid, ')'), ## keep track of gGraph node ids
                             type = "interval",
                             y = pmin(maxcn, nodes$cn)
                         )
                         ##setkey(node.dt, "oid")

                         ## processing edges
                         ed = private$es

                         ## TMPFIX: remove NA edges .. not clear where these are coming from
                         ## but likely the result of trimming / hood
                         ed = ed[!is.na(from) & !is.na(to), ]

                         ed[,":="(soStr = as.character(strand(private$segs[from])),
                                  siStr = as.character(strand(private$segs[to])))]
                         edByType = by(ed, ed$type, function(x) x)

                         ## see which of the ab edges are "+"
                         abe = edByType$aberrant
                         if (!is.null(abe)){
                             abe[, key := paste(from, to, sep="_")]
                             setkey(abe, "key")
                             ## info in ab.edges field

### TMPFIX: until private$abEdges gets updated with $hood $trim
                             ##posAbEd = as.data.table(private$abEdges[,1:2,"+"])[!is.na(from+to)]
                             ##abe = abe[posAbEd[, paste(from, to, sep="_")],-c("key")]
                             abe = abe[,-c("key")]
                         }

                         ## put 3 back together
                         ed = rbindlist(list(edByType$reference[soStr=="+"],
                                             edByType$loose[soStr=="+"],
                                             abe))

                         ## if encountered, switch to 0
                         ## mapping from type field to label in json
                         eType = setNames(c("REF", "ALT", "LOOSE"), c("reference", "aberrant", "loose"))
                         ## processing edges, cont.
                         fmap = node.dt[, .(oid, iid)]; setkey(fmap, oid);
                         rmap = node.dt[, .(rid, iid)]; setkey(rmap, rid);

                         if (nrow(ed)>0){
                             ed.dt =
                                 ed[from %in% c(oid, rid) & to %in% c(oid, rid), ## fix to remove junctions linking NA nodes
                                    .(from,
                                      to,
                                      so = ifelse(soStr=="+", fmap[list(from), iid], rmap[list(from), iid]),
                                      si = ifelse(siStr=="+", fmap[list(to), iid], rmap[list(to), iid]),
                                      so.str = ifelse(soStr=="+",1,-1),
                                      si.str = ifelse(siStr=="+",1,-1),
                                      weight=pmin(maxweight, cn), ## diff than defined in es field
                                      title = paste(' ', from, '->', to),
                                      type = eType[type])] ## removed "by"

                             ##TMPFIX: quick hack to remove dup edges
                             ed.dt = ed.dt[
                                 -which(duplicated(paste(apply(cbind(so*so.str, -si*si.str), 1, function(x) paste(sort(x), collapse = ' '))))), ]

                             ## ## was previously (added filter, removed by and added fmap rmap)
                             ## ed.dt = ed[,.(from,
                             ##               to,
                             ##               so = ifelse(soStr=="+", node.dt[oid == from, iid], node.dt[rid == from, iid]),
                             ##               si = ifelse(siStr=="+", node.dt[oid == to, iid], node.dt[rid == to, iid]),
                             ##               so.str = ifelse(soStr=="+",1,-1),
                             ##               si.str = ifelse(siStr=="+",1,-1),
                             ##               weight=pmin(maxweight, cn), ## diff than defined in es field
                             ##               title = "",
                             ##               type = eType[type]),
                             ##            by=1:nrow(ed)]

                             ## ## need to flip the negative segs to positive
                             ## ed.dt[, sig := ifelse(so<si, ## assuming the sorting of segs
                             ##                       paste0(so * so.str, '_', -si*si.str),
                             ##                       paste0(-si * si.str, '_', so*so.str))]
                             ## ed.dt[!duplicated(sig), ][, cid := seq_along(.I)]
                             ed.dt[, cid := 1:length(from)]
                             ed.dt[,":="(so = so*so.str, si = -si*si.str)]

                             connections.json = ed.dt[, paste0(
                                 c("connections: [", paste(
                                                         "\t{",
                                                         "cid: ", cid,
                                                         ifelse(is.na(so), "", ", source: "),
                                                         ifelse(is.na(so), "", so),
                                                         ifelse(is.na(si), "", ", sink: "),
                                                         ifelse(is.na(si), "", si),
                                                         ", title: ", qw(title),
                                                         ", type: ", qw(type),
                                                         ", weight: ", pmin(maxweight, weight),
                                                         "}",
                                                         sep = "",
                                                         collapse = ',\n'),
                                   "]"),
                                 collapse = '\n')]

                         }

                         ## processing intervals
                         intervals.json = node.dt[, paste0(
                             c("intervals: [", paste(
                                                   "\t{",
                                                   "iid: ", iid,
                                                   ", chromosome: ", chromosome,
                                                   ", startPoint: ", startPoint,
                                                   ", endPoint: ", endPoint,
                                                   ", y: ", y,
                                                   ", title: ", qw(title),
                                                   ", type: ", qw(type),
                                                   ", strand: ", qw(strand),
                                                   "}",
                                                   sep = "",
                                                   collapse = ',\n'),
                               "]"),
                             collapse = '\n')
                             ]

                         ## processing meta info
                         ## DONE: seqlengths
                         require(RColorBrewer)

                         chrs = self$getSeqInfo()
                         if (trim)
                             chrs = chrs[seqnames %in% as.character(seqnames(private$segs))]
                         else
                             chrs = chrs[seqnames %in% levels(seqnames(private$segs))]

                         meta.json =
                             paste('\tmetadata: [\n',
                                   chrs[, paste("\t\t{",
                                                " chromosome: ", qw(seqnames),
                                                ", startPoint: ", 1,
                                                ", endPoint: ", seqlengths,
                                                ", color: ", qw(substr(tolower(brewer.master( max(.I), 'BrBG' )), 1, 7)), " }",
                                        #                                                ", color: ", qw(substr(tolower(rainbow( max(.I) )), 1, 7)), " }",
                                                collapse=",\n",
                                                sep="")],
                                   '\n]')

                         ## assembling the JSON
                         out = paste(c("var data = {",
                                       paste(
                                           c(meta.json,
                                             intervals.json,
                                             connections.json),
                                           collapse = ',\n'
                                       ),"}"),
                                     sep = "")

                         ## TODO: remove any NA. Not legal.
                         ## if (!is.null(file)){
                         ##     writeLines(out, file)
                         ## }
                                        #                         return(out)

                         writeLines(out, paste0(file, '/js/data.js'))
                         message(sprintf('Wrote JSON file of gGraph to %s/index.html', file))
                     },

                     ## self-annotating functions

                     ## dicing up the graph
                     components = function(){
                         private$partition = components(private$g)
                         ## merge +/- complements into 1
                         ## 1) if there are pairs of partitions with the same number of V
                         ## 2) test if they are +/- complement
                         ## 3) if so anneal

                         ## TODO: define a compound gGraph class for holding a series of them
                         nComp = private$partition$no ## total N of parts
                         allComponents = setNames(vector("list", nComp), 1:nComp)
                         eqSizeComps = duplicated(private$partition$csize) ## potential dups

                         ## so I had to write my first for loop in here!!!
                         ## TODO: put 2 disconnected strands back to one graph
                         for (i in 1:nComp) {
                             allComponents[[as.character(i)]] =
                                 self$subgraph(which(private$partition$membership==i))
                         }
                         return(allComponents)
                     },

                     subgraph = function(v=numeric(0)){
                         "Given a numeric vector of vertices, return the subgraph consists only thesevertices "
                         if (length(v)==0){
                             ## nothing provided, nothing happens
                             return(self)
                         } else if (is.numeric(v)){
                             ## at least they are num
                             if (!is.integer(v)){
                                 ## if not integer, convert
                                 v = as.integer(v)
                             }
                             if (!all(v %in% seq_along(private$segs))){
                                 v = v[which(v %in% seq_along(private$segs))]
                                 warning("Some v subscripts out of bound! Ignore!")
                             }
                             ## get the subgraph
                             newSegs = private$segs[v]
                             newId = setNames(seq_along(v), v)
                             newEs = private$es[from %in% v | to %in% v,][, ":="(from = newId[as.character(from)],
                                                                                 to = newId[as.character(to)])]

                             jIdx = which(grl.in(private$junction$grl, newSegs, only=T))
                             newJuncs = private$junction[unique(jIdx)]

                             out = gGraph$new(segs=newSegs,
                                              es=newEs,
                                              junctions=newJuncs,
                                              ploidy=private$ploidy,
                                              purity=private$purity)
                             return(out)
                         } else {
                             stop("Invalid input.")
                         }
                     },
                     trim = function(gr=NULL){
                         ## TODO
                         "Given a GRanges, return the trimmed subgraph overlapping it."
                         if (is.null(gr))
                             return(self)

                         gr = gr.fix(gr, get(self$refG))
                         gr = gr.stripstrand(gr)
                         if (!isDisjoint(gr))
                             gr = gr.reduce(gr)

                         v = which(gr.in(private$segs, gr))
                         sg = self$subgraph(v)
                         if (length(v)<=2)
                             return(sg)## TODO: resolve the edge case where gr is contained in single node

                         nss = sg$segstats
                         grS = gr.start(gr)
                         grE = gr.end(gr)
                         ## find if any gr start/end inside segs
                         sInSeg = gr.findoverlaps(grS, nss)
                         eInSeg = gr.findoverlaps(grE, nss)

                         ## for start point inside segs, split node and keep the right part
                         brByS = gr.breaks(nss[sInSeg$subject.id], grS)
                         lastCol = ncol(mcols(brByS))
                         spByS = by(brByS, brByS$qid,
                                    function(gr) {
                                        if (!is(gr, "GRanges"))
                                            gr = GRanges(gr)
                                        gr[2, -lastCol]
                                    })
                         nss[sInSeg$subject.id] = Reduce("c",unlist(spByS))

                         ## for end point inside segs, split node and keep the left part
                         brByE = gr.breaks(nss[eInSeg$subject.id], grE)
                         spByE = by(brByE, brByE$qid,
                                    function(gr) {
                                        if (!is(gr, "GRanges"))
                                            gr = GRanges(gr)
                                        gr[1, -lastCol]
                                    })
                         nss[eInSeg$subject.id] = Reduce("c",unlist(spByE))

                         whichTerminal = sg$edges[, setxor(from, to)]
                         nss$terminal = seq_along(nss) %in% whichTerminal

                         newSg = gGraph$new(segs=nss,
                                            es=sg$edges,
                                            junctions=sg$junctions,
                                            ploidy=private$ploidy,
                                            purity=private$purity)
                         return(newSg)
                     },

                     getSeqInfo = function(){
                         as.data.table(attributes(seqinfo(get(self$refG))))
                     },
                     makeAbEdges = function(){
                         ## TODO: derive abEdges from junction
                         if (length(junctions)==0){
                             return(
                                 array(dim=c(0,3,2),
                                       dimnames=list(NULL,
                                                     c("from", "to", "edge.ix"),
                                                     c("+","-")))
                             )
                         } else {
                             ## based on junctions, get
                             junc = private$junction$grl
                             abe = private$es[type=="aberrant"]
                             abEdges = array(dim=c(length(junc),3,2),
                                   dimnames=list(NULL, c("from", "to", "edge.ix"), c("+","-")))
                             ## find coresponding edge.ix for abe
                             jUl = unlist(junc)
                             jUl$jix = rep(seq_along(junc), each=2)
                             seg = private$segs %Q% (loose==F)
                             sx = jUl %*% seg[,c()]

                             for (i in seq_along(junc)){
                                 segid = (sx %Q% (jix == i))$subject.id
                                 if (length(unique(segid)!=4) | length(segid)==0){
                                     next
                                 } else {
                                     edge.ix = abe[, which(from %in% segid & to %in% segid)]
                                     thisAbe =
                                         t(as.matrix(abe[edge.ix, .(from, to, edge.ix=edge.ix)]))
                                     abEdges[i,,] = thisAbe
                                 }
                             }
                             return(abEdges)
                         }

                         return(private$abEdges)
                     },
                     getAdj = function(flat=FALSE){
                         adjMat = as_adj(private$g)
                         if (flat) {
                             return(adjMat)
                         } else {
                             adjMat[as.matrix(private$es[,.(from, to)])] = private$es$cn
                             return(adjMat)
                         }
                     },
                     ## some query functions
                     hood = function(win, d=NULL, k=NULL, pad=0,
                                     bagel=FALSE, ignore.strand=T, verbose=FALSE){
                         "Get the trimmed subgraph around a given GRanges within a distance on the graph."
                         if (ignore.strand)
                             win = gr.stripstrand(win)

                         ## TODO: what to do when win is larger than segs?????
                         ## ans: return self
                         ## overlapping window and segs, removing loose ends
                         interGr = gr.findoverlaps(private$segs, win, ignore.strand=ignore.strand)
                         lid = which(private$segs$loose==T)
                         interGr = interGr %Q% (!query.id %in% lid)
                         qix = interGr$query.id

                         if (is.null(k)){
                             ## TODO!!!
                             ## no k, use distance
                             if (is.null(d) | d < 0)
                                 stop("Must provide either valid k or d.")

                             ## blend window with segs
                             win = gr.fix(win, get(self$refG))## fix seqinfo
                             ss = tryCatch(c(private$segs[private$segs$loose == F, c()],
                                             win[, c()]), error = function(e) NULL)

                             if (is.null(ss))
                                 ss = grbind(c(jab$segstats[jab$segstats$loose == FALSE, c()], win[, c()]))

                             if (ignore.strand)
                                 ss = gr.stripstrand(ss)

                             ## break it into non-overlapping segs
                             ss = disjoin(ss)
                             ## update win
                             win = gr.findoverlaps(ss, win, ignore.strand = ignore.strand)

                             ## start/end of ss
                             seg.s = suppressWarnings(gr.start(ss, ignore.strand = TRUE))
                             seg.e = suppressWarnings(gr.end(ss, ignore.strand = TRUE))
                             ## distance from all of win to start/end of ss
                             ## TODO: connect with dist
                             D.s = suppressWarnings(self$dist(win, seg.s, verbose = verbose))
                             D.e = suppressWarnings(self$dist(win, seg.e, verbose = verbose))

                             ## shortest path distance
                             min.s = apply(D.s, 2, min, na.rm = TRUE)
                             min.e = apply(D.e, 2, min, na.rm = TRUE)
                             ## which idx bear them?
                             s.close = min.s<=d
                             e.close = min.e<=d

                             ## generate new gGraph, trim the subgraph
                             out = GRanges()
                             if (any(s.close))
                                 out = c(out,
                                         GenomicRanges::flank(seg.s[s.close],
                                                              -(d-min.s[s.close]))
                                         )

                             if (any(e.close))
                                 out = c(out,
                                         GenomicRanges::shift(flank(seg.e[e.close],
                                                                    d-min.e[e.close]),1)
                                         )

                             if (!bagel)
                                 out = streduce(c(win[, c()], out[, c()]))

                             hoodRange = streduce(out, pad)
                             return(self$trim(hoodRange))
                         } else {
                             ## with k, go no more k steps
                             kNeighbors = unique(unlist(ego(private$g, qix, order=k)))
                             return(self$subgraph(kNeighbors)) ## not garanteed size to scale
                         }
                     },

                     dist = function(gr1, gr2,
                                     matrix=T, EPS=1e-9,
                                     include.internal=TRUE, ## consider bp within feature "close"
                                     directed=FALSE, ## if TRUE, only consider gr1-->gr2 paths
                                     verbose=FALSE){
                         "Given two GRanges, return pairwise shortest path distance."
                         ## TODO
                         if (verbose)
                             now = Sys.time()

                         intersect.ix = gr.findoverlaps(gr1, gr2, ignore.strand = FALSE)

                         ngr1 = length(gr1)
                         ngr2 = length(gr2)

                         tiles = private$segs
                         G = private$g

                         ## keep track of original ids when we collapse
                         gr1$id = 1:length(gr1)
                         gr2$id = 1:length(gr2)

                         ## check for double stranded intervals
                         ## add corresponding nodes if present
                         if (any(ix <- strand(gr1)=='*'))
                         {
                             strand(gr1)[ix] = '+'
                             gr1 = c(gr1, gr.flipstrand(gr1[ix]))
                         }

                         if (any(ix <- strand(gr2)=='*'))
                         {
                             strand(gr2)[ix] = '+'
                             gr2 = c(gr2, gr.flipstrand(gr2[ix]))
                         }

                         ## expand nodes by jabba model to get internal connectivity
                         if (include.internal)
                         {
                             gr1 = gr1[, 'id'] %**% private$segs
                             gr2 = gr2[, 'id'] %**% private$segs
                         }

                         if (verbose)
                         {
                             message('Finished making gr objects')
                             print(Sys.time() -now)
                         }

                         tmp = get.edges(G, E(G))
                         E(G)$from = tmp[,1]
                         E(G)$to = tmp[,2]
                         E(G)$weight = width(tiles)[E(G)$to]

                         gr1.e = gr.end(gr1, ignore.strand = FALSE)
                         gr2.s = gr.start(gr2, ignore.strand = FALSE)

                         if (!directed)
                         {
                             gr1.s = gr.start(gr1, ignore.strand = FALSE)
                             gr2.e = gr.end(gr2, ignore.strand = FALSE)
                         }

                         ## graph node corresponding to end of gr1.ew
                         gr1.e$ix = gr.match(gr1.e, tiles, ignore.strand = F)
                         ## graph node corresponding to beginning of gr2
                         gr2.s$ix= gr.match(gr2.s, tiles, ignore.strand = F)

                         if (!directed)
                         {
                             ## graph node corresponding to end of gr1.ew
                             gr1.s$ix = gr.match(gr1.s, tiles, ignore.strand = F)
                             ## graph node corresponding to beginning of gr2
                             gr2.e$ix= gr.match(gr2.e, tiles, ignore.strand = F)
                         }

                         ## 3' offset from 3' end of query intervals to ends of jabba segs
                         ## to add / subtract to distance when query is in middle of a node
                         off1 = ifelse(as.logical(strand(gr1.e)=='+'),
                                       end(tiles)[gr1.e$ix]-end(gr1.e),
                                       start(gr1.e) - start(tiles)[gr1.e$ix])
                         off2 = ifelse(as.logical(strand(gr2.s)=='+'),
                                       end(tiles)[gr2.s$ix]-end(gr2.s),
                                       start(gr2.s) - start(tiles)[gr2.s$ix])

                         ## reverse offset now calculate 3' offset from 5' of intervals
                         if (!directed)
                         {
                             off1r = ifelse(as.logical(strand(gr1.s)=='+'),
                                            end(tiles)[gr1.s$ix]-start(gr1.s),
                                            end(gr1.s) - start(tiles)[gr1.s$ix])
                             off2r = ifelse(as.logical(strand(gr2.e)=='+'),
                                            end(tiles)[gr2.e$ix]-start(gr2.e),
                                            end(gr2.e) - start(tiles)[gr2.e$ix])
                         }

                         ## compute unique indices for forward and reverse analyses
                         uix1 = unique(gr1.e$ix)
                         uix2 = unique(gr2.s$ix)

                         if (!directed)
                         {
                             uix1r = unique(gr1.s$ix)
                             uix2r = unique(gr2.e$ix)
                         }

                         ## and map back to original indices
                         uix1map = match(gr1.e$ix, uix1)
                         uix2map = match(gr2.s$ix, uix2)

                         if (!directed)
                         {
                             uix1mapr = match(gr1.s$ix, uix1r)
                             uix2mapr = match(gr2.e$ix, uix2r)
                         }

                         adj = self$getAdj()
                         self.l = which(diag(adj)>0)

                         if (verbose)
                         {
                             message('Finished mapping gr1 and gr2 objects to jabba graph')
                             print(Sys.time() -now)
                         }

                         ## need to take into account forward and reverse scenarios of "distance" here
                         ## ie upstream and downstream connections between query and target
                         ## edges are annotated with width of target

                         ## so for "downstream distance"  we are getting matrix of shortest paths between from uix1 and uix2 node pair
                         ## and then correcting those distances by (1) adding the 3' offset of uix1 from its node
                         ## and (2) subtracting the 3' offset of uix2
                         Df = sweep(
                             sweep(
                                 shortest.paths(G, uix1, uix2, weights = E(G)$weight,
                                                mode = 'out')[uix1map, uix2map, drop = F],
                                 1, off1, '+'), ## add uix1 3' offset to all distances
                             2, off2, '-') ## subtract uix2 3' offset to all distances


                         if (!directed)
                         {
                             ## now looking upstream - ie essentially flipping edges on our graph - the edge weights
                             ## now represent "source" node widths (ie of the flipped edges)
                                        # need to correct these distances by (1) subtracting 3' offset of uix1 from its node
                             ## and (2) adding the 3' offset of uix2
                             ## and using the reverse indices
                             Dr = sweep(
                                 sweep(
                                     t(shortest.paths(G, uix2r, uix1r, weights = E(G)$weight, mode = 'out'))[uix1mapr, uix2mapr, drop = F],
                                     1, off1r, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                                 2, off2r , '+') ## add uix2 offset to all distances

                             Df2 = sweep(
                                 sweep(
                                     shortest.paths(G, uix1r, uix2, weights = E(G)$weight, mode = 'out')[uix1mapr, uix2map, drop = F],
                                     1, off1r, '+'), ## add uix1 3' offset to all distances
                                 2, off2, '-') ## subtract uix2 3' offset to all distances

                             Dr2 = sweep(
                                 sweep(
                                     t(shortest.paths(G, uix2r, uix1, weights = E(G)$weight, mode = 'out'))[uix1map, uix2mapr, drop = F],
                                     1, off1, '-'), ## substract  uix1 offset to all distances but subtract weight of <first> node
                                 2, off2r , '+') ## add uix2 offset to all distances
                             D = pmin(abs(Df), abs(Dr), abs(Df2), abs(Dr2))
                         }
                         else
                             D = Df

                         if (verbose)
                         {
                             message('Finished computing distances')
                             print(Sys.time() -now)
                         }


                         ## take care of edge cases where ranges land on the same node, since igraph will just give them "0" distance
                         ## ij contains pairs of gr1 and gr2 indices that map to the same node
                         ij = as.matrix(merge(cbind(i = 1:length(gr1.e), nid = gr1.e$ix), cbind(j = 1:length(gr2.s), nid = gr2.s$ix)))

                         ## among ij pairs that land on the same (strand of the same) node
                         ##
                         ## several possibilities:
                         ## (1) if gr1.e[i] < gr2.s[j] then keep original distance (i.e. was correctly calculated)
                         ## (2) if gr1.e[i] > gr2.s[j] then either
                         ##   (a) check if there is a self loop and adjust accordingly (i.e. add back width of current tile)
                         ##   (b) PITA case, compute shortest distance from i's child(ren) to j

                         if (nrow(ij)>0)
                         {
                             ## rix are present
                             rix = as.logical((
                                 (strand(gr1)[ij[,'i']] == '+' & strand(gr2)[ij[,'j']] == '+' & end(gr1)[ij[,'i']] <= start(gr2[ij[,'j']])) |
                                 (strand(gr1)[ij[,'i']] == '-' & strand(gr2)[ij[,'j']] == '-' & start(gr1)[ij[,'i']] >= end(gr2)[ij[,'j']])))

                             ij = ij[!rix, , drop = F] ## NTD with rix == TRUE these since they are calculated correctly

                             if (nrow(ij)>0) ## any remaining will either be self loops or complicated loops
                             {
                                 selfix = (ij[, 'nid'] %in% self.l)

                                 if (any(selfix)) ## correct distance for direct self loops (add back width of current node)
                                     D[ij[selfix, c('i', 'j'), drop = F]]  = D[ij[selfix, c('i', 'j'), drop = F]] + width(tiles)[ij[selfix, 'nid']]

                                 ij = ij[!selfix, , drop = F]

                                 if (nrow(ij)>0) ## remaining are pain in the ass indirect self loops
                                 {
                                     ch = G[[ij[, 'nid']]] ## list of i nodes children for all remaining ij pairs
                                     chu = munlist(ch) ## unlisted children, third column are the child id's, first column is the position of nrix

                                     ## now find paths from children to corresponding j
                                     epaths = suppressWarnings(get.shortest.paths(G, chu[, 3], ij[chu[,'ix'], 'nid'], weights = E(G)$weight, mode = 'out', output = 'epath')$epath)
                                     epathw = sapply(epaths, function(x,w) if (length(x)==0) Inf else sum(w[x]), E(G)$weight) ## calculate the path weights
                                     epathw = epathw + width(tiles)[chu[, 3]] + off1[ij[chu[, 'ix'], 'i']] + off2[ij[chu[,'ix'], 'j']] - width(tiles)[ij[chu[, 'ix'], 'nid']]

                                     ## aggregate (i.e. in case there are multiple children per node) by taking min width
                                     D[ij[, c('i', 'j'), drop = F]] = vaggregate(epathw, by = list(chu[, 'ix']), min)[as.character(1:nrow(ij))]
                                 }
                             }
                         }

                         if (verbose)
                         {
                             message('Finished correcting distances')
                             print(Sys.time() -now)
                         }

                         ## need to collapse matrix ie if there were "*" strand inputs and if we are counting internal
                         ## connections inside our queries ..
                         ## collapsing +/- rows and columns by max value based on their id mapping to their original "*" interval

                         ## melt distance matrix into ij
                         ij = as.matrix(expand.grid(1:nrow(D), 1:ncol(D)))
                         dt = data.table(i = ij[,1], j = ij[,2], value = D[ij])[, id1 := gr1$id[i]][, id2 := gr2$id[j]]

                         tmp = dcast.data.table(dt, id1 ~ id2, fun.aggregate = function(x) min(as.numeric(x)))
                         setkey(tmp, id1)
                         D = as.matrix(tmp[list(1:ngr1), -1, with = FALSE])[, as.character(1:ngr2), drop = FALSE]

                         ## finally zero out any intervals that actually intersect
                         ## (edge case not captured when we just examine ends)
                         if (length(intersect.ix)>0)
                             D[cbind(intersect.ix$query.id, intersect.ix$subject.id)] = 0

                         if (verbose)
                         {
                             message('Finished aggregating distances to original object')
                             print(Sys.time() -now)
                         }

                         return(D)
                     },
                     proximity = function(query, subject,
                                          verbose=F, mc.cores=1,
                                          max.dist=1e6){

                         ## TODO:
                         adj = self$getAdj()
                         ix = which(adj[private$abEdges[,1:2,1]]>0)
                         if (length(ix)>0) {
                             ra1 = gr.flipstrand(gr.end(private$segs[private$abEdges[ix,1,1]], 1, ignore.strand = F))
                             ra2 = gr.start(private$segs[private$abEdges[ix,2,1]], 1, ignore.strand = F)
                             ra1 = GenomicRanges::shift(ra1, ifelse(as.logical(strand(ra1)=='+'), -1, 0))
                             ra2 = GenomicRanges::shift(ra2, ifelse(as.logical(strand(ra2)=='+'), -1, 0))
                             ra = grl.pivot(GRangesList(ra1,ra2))
                         }

                         if (!is(query, "GRanges") & !is(query, "GRanges"))
                             stop("Invalid input")

                         if (length(query)==0 | length(subject)==0)
                             return(list())

                         if (is.null(names(query)))
                             names(query) = 1:length(query)

                         if (is.null(names(subject)))
                             names(subject) = 1:length(subject)

                         query.nm = names(query);
                         subject.nm = names(subject);

                         query = query[, c()]
                         subject = subject[, c()]

                         query$id = 1:length(query)
                         subject$id = 1:length(subject)

                         qix.filt = gr.in(query, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
                         query = query[qix.filt]

                         six.filt = gr.in(subject, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
                         subject = subject[six.filt]

                         if (length(query)==0 | length(subject)==0)
                             return(list())

                         query$type = 'query'
                         subject$type = 'subject'

                         subject = gr.fix(subject, get(self$refG))
                         query = gr.fix(query, get(self$refG))
                         gr = c(query, subject)

                         kg = karyograph(ra, gr)
                         ## TODO: make karyograph output compatible with Marcin's!!!
                         ## kg2 = gGraph$new()$karyograph(gr, ra)

                         ## node.start and node.end delinate the nodes corresponding to the interval start and end
                         ## on both positive and negative tiles of the karyograph
                         gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;

                         ## start and end indices of nodes
                         tip = which(strand(kg$tile)=='+')
                         tin = which(strand(kg$tile)=='-')
                         gr$node.start = tip[gr.match(gr.start(gr,2), gr.start(kg$tile[tip]))]
                         gr$node.end = tip[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tip]))]
                         gr$node.start.n = tin[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tin]))]
                         gr$node.end.n = tin[gr.match(gr.start(gr,2), gr.start(kg$tile[tin]))]

                                        #    gr$node.start = gr.match(gr.start(gr-1,2), gr.start(kg$tile))
                                        #    gr$node.end = suppressWarnings(gr.match(gr.end(gr+1,2), gr.end(kg$tile)))

                         ## so now we build distance matrices from query ends to subject starts
                         ## and subject ends to query starts

                         ## so for each query end we will find the shortest path to all subject starts
                         ## and for each query start we will find the shortest.path from all subject ends
                         ix.query = which(gr$type == 'query')
                         ix.subj = which(gr$type == 'subject')

                         node.start = gr$node.start
                         node.end = gr$node.end
                         node.start.n = gr$node.start.n
                         node.end.n = gr$node.end.n

                         w = width(kg$tile)

                         E(kg$G)$weight = width(kg$tile)[E(kg$G)$to]

                         ## ix.query and ix.subj give the indices of query / subject in gr
                         ## node.start, node.end map gr to graph node ids
                         ##
                         ## these matrices are in dimensions of query and subject, and will hold the pairwise distances between
                         ##
                         D.rel = D.ra = D.ref = D.which = Matrix(data = 0, nrow = length(ix.query), ncol = length(ix.subj))

                         ## "reference" graph (missing aberrant edges)
                         G.ref = subgraph.edges(kg$G, which(E(kg$G)$type == 'reference'), delete.vertices = F)

                         EPS = 1e-9

                         ## for (i in ix.query)
                         tmp = mclapply(ix.query, function(i)
                         {
                             if (verbose)
                                 cat('starting interval', i, 'of', length(ix.query), '\n')

                             ## D1 = shortest query to subject path, D2 = shortest subject to query path, then take shortest of D1 and D2
                             ## for each path, the edge weights correspond to the interval width of the target node, and to compute the path
                             ## length we remove the final node since we are measuring the distance from the end of the first vertex in the path
                             ## to the beginning of the final vertex

                             u.node.start = unique(node.start[ix.subj]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
                             u.node.end = unique(node.end[ix.subj])

                             uix.start = match(node.start[ix.subj], u.node.start)
                             uix.end = match(node.end[ix.subj], u.node.end)

                             tmp.D1 = (shortest.paths(kg$G, node.end[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
                             tmp.D2 = (shortest.paths(kg$G, node.start[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start[i]])[uix.end]
                             tmp.D3 = (shortest.paths(kg$G, node.end.n[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
                             tmp.D4 = (shortest.paths(kg$G, node.start.n[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
                             tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
                             ix = which(tmp.D<max.dist)
                             D.ra[i, ix] = tmp.D[ix]+EPS
                             D.which[i, ix] = apply(cbind(tmp.D1[ix], tmp.D2[ix], tmp.D3[ix], tmp.D4[ix]), 1, which.min)

                             u.node.start = unique(node.start[ix.subj][ix]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
                             u.node.end = unique(node.end[ix.subj][ix])

                             uix.start = match(node.start[ix.subj][ix], u.node.start)
                             uix.end = match(node.end[ix.subj][ix], u.node.end)

                             tmp.D1 = (shortest.paths(G.ref, node.end[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
                             tmp.D2 = (shortest.paths(G.ref, node.start[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start[i]])[uix.end]
                             tmp.D3 = (shortest.paths(G.ref, node.end.n[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
                             tmp.D4 = (shortest.paths(G.ref, node.start.n[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
                             tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
                             D.ref[i, ix] = tmp.D+EPS

                             ## if subject and query intersect (on the reference) then we count both RA and Ref distance as 0
                             ## (easier to do a simple range query here)
                             ix.zero = gr.in(subject[ix], query[i])
                             if (any(ix.zero))
                             {
                                 D.ra[i, ix[ix.zero]] = 0
                                 D.ref[i, ix[ix.zero]] = 0
                             }

                             D.rel[i, ix] = ((D.ra[i, ix]-EPS) / (D.ref[i, ix]-EPS)) + EPS


                             if (verbose)
                                 cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')

                             return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))
                         }, mc.cores = mc.cores)

                         for (i in 1:length(tmp))
                         {
                             if (class(tmp[[i]]) != 'list')
                                 warning(sprintf('Query %s failed', ix.query[i]))
                             else
                             {
                                 D.rel = D.rel + tmp[[i]]$D.rel
                                 D.ra = D.ra + tmp[[i]]$D.ra
                                 D.ref = D.ref + tmp[[i]]$D.ref
                                 D.which = D.which + tmp[[i]]$D.which
                             }
                         }

                         ## "full" size matrix
                         rel = ra = ref = ra.which =
                             Matrix(data = 0, nrow = length(qix.filt), ncol = length(six.filt), dimnames = list(dedup(query.nm), dedup(names(subject.nm))))
                         rel[qix.filt, six.filt] = D.rel
                         ra[qix.filt, six.filt] = D.ra
                         ref[qix.filt, six.filt] = D.ref
                         ra.which[qix.filt, six.filt] = D.which

                         ## summary is data frame that has one row for each query x subject pair, relative distance, ra distance, and absolute distance
                         tmp = which(rel!=0, arr.ind = T)
                         colnames(tmp) = c('i', 'j');
                         sum = as.data.frame(tmp)

                         if (!is.null(query.nm))
                             sum$query.nm = query.nm[sum$i]

                         if (!is.null(subject.nm))
                             sum$subject.nm = subject.nm[sum$j]

                         sum$rel = rel[tmp]
                         sum$ra = ra[tmp]
                         sum$wt = ref[tmp]

                         sum = sum[order(sum$rel), ]
                         sum = sum[sum$rel<1, ] ## exclude those with rel == 1

                         ## reconstruct paths
                         vix.query = matrix(NA, nrow = length(qix.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
                         vix.subject = matrix(NA, nrow = length(six.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
                         vix.query[qix.filt, ] = cbind(values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start.n')], values(gr)[ix.query, c('node.end.n')])
                         vix.subject[six.filt] = cbind(values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start.n')], values(gr)[ix.subj, c('node.end.n')])

                         browser()

                         sum.paths = mapply(function(x, y)
                         {
                             if ((ra.which[x, y]) == 1)
                                 get.shortest.paths(kg$G, vix.query[x, 'end'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
                             else if ((ra.which[x, y]) == 2)
                                 rev(get.shortest.paths(kg$G, vix.query[x, 'start'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
                             else if ((ra.which[x, y]) == 3)
                                 get.shortest.paths(kg$G, vix.query[x, 'end.n'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
                             else if ((ra.which[x, y]) == 4)
                                 rev(get.shortest.paths(kg$G, vix.query[x, 'start.n'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
                         }, sum$i, sum$j, SIMPLIFY = F)

                                        #    sum$paths = lapply(sum.paths, function(x) x[-c(1, length(x))])
                         sum$paths = sum.paths
                         sum$ab.edges = lapply(sum.paths, function(p) setdiff(E(kg$G, path = p)$bp.id, NA))
                         return(list(sum = sum, rel = rel, ra = ra, wt = ref, G = kg$G, G.ref = G.ref, tile = kg$tile, vix.query = vix.query, vix.subject = vix.subject))

                     },
                     jGraph = function(){
                         ##TODO: migrate the jGraph function here

                     },

                     ## property constraints
                     isJunctionBalanced = function(){
                         ## TODO: use adj to calc if every segment is balanced on both sides
                     },
                     isDoubleStrand = function(){
                         ## TODO: test if segs come in +/- pairs
                         identical((ss %Q% (strand=="-"))[, c()],
                                   gr.flipstrand(ss %Q% (strand=="+"))[, c()])
                     }
                 ),

                 private = list(
                     ## private fields
                     ## node/vertex, a GRanges obj of strand-specific ranges
                     segs = NULL,
                     ## temporary segs for backtrace when modifying segs
                     tmpSegs = NULL,
                     ## igraph obj representing the graph structure
                     g = NULL,
                     ## data.table of all edges in g, from, to, cn, type
                     ## type can be ref, aberrant, loose
                     es = data.table(from=integer(0), to=integer(0),
                                     cn=integer(0), type=character(0), weight=numeric(0)),
                     ## putative junctions, junctions
                     junction = junctions$new(),
                     abEdges = array(dim=c(0,3,2), dimnames=list(NULL, c("from", "to", "edge.ix"), c("+","-"))),
                     ## ploidy is set to 2, only to init null graph,
                     ## otherwise inferred from segs
                     ploidy = 2,
                     ## tumor cell proportion
                     purity = 1,
                     partition = NULL,

                     ## private methods
                     ## break the current segments into new segments
                     makeSegs = function(bps){
                         ## DONE: once finished, move to private methods
                         private$tmpSegs = private$segs
                         private$segs = gr.breaks(private$segs, bps)
                         return(self)
                     },
                     ## initialize by directly giving fields values
                     gGraphFromScratch = function(segs, es, junctions, ploidy, purity){

                         private$segs = segs
                         private$es = es
                         private$g = make_directed_graph(
                             t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))
                         private$junction$append(junctions)
                         private$abEdges = self$makeAbEdges()
                         private$ploidy = ploidy
                         private$purity = purity
                     },
                     ## collapse strand info
                     getSs = function(){
                         "Return simple segs, with names, tile.id, is.tel, ab.source, ab.target."

                         ## ## TODO: think about how did he plot loose ends!!!
                         ## ## processing nodes
                         ## ## reduce strand
                         ## ## remove loose nodes
                         ## oid = gr2dt(private$segs)[, which(strand == "+" & loose==F)]
                         ## ## ori ind of rev comps
                         ## rid = gr2dt(private$segs)[, which(strand == "-" & loose==F)]

                         ## ## single strand
                         ## ss = gr.stripstrand(private$segs[oid])
                         ## newMap = match(gr.stripstrand(private$segs), ss)

                         ## ## ori ix of loose nodes
                         ## lid = which(private$segs$loose==T)

                         ## ## processing edges
                         ## ed = private$es
                         ## ed[,":="(soStr = as.character(strand(private$segs[from])),
                         ##          siStr = as.character(strand(private$segs[to])))]
                         ## edByType = by(ed, ed$type, function(x) x)

                         ## ## see which of the ab edges are "+"
                         ## abe = edByType$aberrant
                         ## if (!is.null(abe)){
                         ##     abe[, key := paste(from, to, sep="_")]
                         ##     setkey(abe, "key")
                         ##     ## info in ab.edges field
                         ##     posAbEd = as.data.table(private$abEdges[,1:2,"+"])[!is.na(from+to)]
                         ##     abe = abe[posAbEd[, paste(from, to, sep="_")],-c("key")]
                         ## }

                         ## ## put 3 back together
                         ## ed = rbindlist(list(edByType$reference[soStr=="+"],
                         ##                     edByType$loose[soStr=="+"],
                         ##                     abe))

                         ## ## processing edges, cont.
                         ## if (nrow(ed)>0){
                         ##     ed[, ":="(newFr = newMap[from], newTo = newMap[to])]
                         ## }

                     }
                 ),

                 active = list(
                     edges = function(){
                         return(private$es)
                     }, ## igraph.es
                     segstats = function(){
                         return(private$segs)
                     },
                     td = function(){
                         return(self$gGraph2gTrack())
                     },
                     G = function(){
                         return(private$g)
                     },
                     tmpS = function(){
                         return(private$tmpSegs)
                     },
                     seqInfo = function(){
                         return(self$getSeqInfo())
                     },
                     junctions = function(){
                         return(private$junction$grl)
                     },
                     igPlot = function(){
                         ## TODO: make igraph plot
                         return(self$gGraph2iGraphViz())
                     },
                     parts = function(){
                         ## DONE: use the correct components function
                         return(self$components())
                     },
                     json = function(file='~/public_html/gGraph'){
                         return(self$gGraph2json(file=file))
                     },
                     adj = function(){

                         return(self$getAdj())
                     },
                     ab.edges = function(){
                         return(private$abEdges)
                     },
                     tile = function(){
                         return(self$getSs())
                     }
                 )
                 )

components <- function (x, ...) {
    UseMethod("components", x)
}
components.igraph <- function(iGraph){
    return(igraph::components(iGraph))
}
#' components -- same idea as igraph::components, returns a list of gGraph objects
#'
#' @return a list of gGraph objects representing each partition of the input
#'
components.gGraph <- function(gGraph){
    ## input must be a gGraph!
    if (!is(gGraph, "gGraph")){
        stop("Invalid input.")
    }
    return(gGraph$components())
}

#' Descendant of gGraph class, where junction balance restraint must be met all the time
#'
bGraph = R6Class("bGraph",
                 inherit = "gGraph",
                 public = list(
                     ## overwrite constructor: restrict about junction balance
                     initialize = function(gG=NULL, jabba=NULL){
                         if (is.null(jabba)){
                             self$nullGGraph()
                         } else if (!is.null(gG)) {
                             if (is(gG, "gGraph") & isJunctionBalanced(gG)){
                                 gG
                             } else {
                                 stop("Invalid input gG.")
                             }
                           } else {
                               self$jabba2gGraph(jabba)
                           }
                       },

                       ## accumulate new events
                       ## given a ref range, if doable, do it
                       dsb = function(){},
                       del = function(){},
                       tDup = function(){},
                       invs = function(){},

                       ## decompose graph into all possible haplotypes
                       walk = function(){
                           "Give all the possible multiset of walks that can be represented by this graph."
                           ## TODO: only balanced graph can walk

                       }
                   ),
                   private = list(

                   ),
                   active = list())
#'
#' gWalks: subclass to gGraph
gWalks = R6Class("gWalks",
                 public=list(
                     initialize = function(){

                     }
                 ),
                 private=list(),
                 active=list())

## Utility functions
#' ra_breaks: utility function to read junction data from various common formats
#'
#' @name ra_breaks
#' @import VariantAnnotation
#' @export
#'
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T,
                     snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, skip = NA,
                     get.loose = FALSE){## if TRUE will return a list with fields $junctions and $loose.ends
    if (is.character(rafile))
    {
        if (grepl('(.bedpe$)', rafile))
        {
            ra.path = rafile
            cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

            ln = readLines(ra.path)
            if (is.na(skip))
            {
                nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                if (is.infinite(nh))
                    nh = 1
            }
            else
                nh = skip


            if ((length(ln)-nh)==0)
                if (get.loose)
                    return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                else
                    return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                                        #                          return(GRangesList())


            if (nh ==0)
                rafile = fread(rafile, header = FALSE)
            else
            {

                rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                if (is.null(rafile))
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)

                if (is.null(rafile))
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)

                if (is.null(rafile))
                    stop('Error reading bedpe')
            }
            setnames(rafile, 1:length(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
                                        #                      rafile[, str1 := ifelse(str1=='+', '-', '+')]
                                        #                      rafile[, str2 := ifelse(str2=='+', '-', '+')]

        }
        else if (grepl('(vcf$)|(vcf.gz$)', rafile))
        {
            library(VariantAnnotation)

            vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))
            if (!('SVTYPE' %in% names(info(vcf)))) {
                warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                return(GRangesList());
            }

            ## vgr = rowData(vcf) ## parse BND format
            vgr = read_vcf(rafile, swap.header = swap.header)

            ## no events
            if (length(vgr) == 0)
                return (GRangesList())

            ## fix mateids if not included
            if (!"MATEID"%in%colnames(mcols(vgr))) {
                nm <- vgr$MATEID <- names(vgr)
                ix <- grepl("1$",nm)
                vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                vgr$SVTYPE="BND"
            }

            if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                stop("MATEID or SVTYPE not included. Required")

            vgr$mateid = vgr$MATEID
            if (is.null(vgr$SVTYPE))
                vgr$svtype = vgr$SVTYPE
            else
                vgr$svtype = vgr$SVTYPE

            if (!is.null(info(vcf)$SCTG))
                vgr$SCTG = info(vcf)$SCTG

            if (force.bnd)
                vgr$svtype = "BND"

            if (sum(vgr$svtype == 'BND')==0)
                warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')

            if (!all(vgr$svtype == 'BND'))
                warning(sprintf('%s rows of vcf do not have svtype BND, ignoring these', sum(vgr$svtype != 'BND')))

            bix = which(vgr$svtype == "BND")
            vgr = vgr[bix]
            alt <- sapply(vgr$ALT, function(x) x[1])
            vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
            vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
            vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
            vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
            vgr$mcoord = gsub('chr', '', vgr$mcoord)

            if (all(is.na(vgr$mateid)))
                if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                {
                    warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                    vgr$mateid = paste(gsub('::\\d$', '', names(vgr)), (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
                }
                else if (!is.null(vgr$SCTG))
            {
                warning('MATEID tag missing, guessing BND partner from coordinates and SCTG')
                require(igraph)
                ucoord = unique(c(vgr$coord, vgr$mcoord))
                vgr$mateid = paste(vgr$SCTG, vgr$mcoord, sep = '_')

                if (any(duplicated(vgr$mateid)))
                {
                    warning('DOUBLE WARNING! inferred mateids not unique, check VCF')
                    bix = bix[!duplicated(vgr$mateid)]
                    vgr = vgr[!duplicated(vgr$mateid)]
                }
            }
                else
                    stop('MATEID tag missing')

            vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

            pix = which(!is.na(vgr$mix))

            vgr.pair = vgr[pix]

            if (length(vgr.pair)==0)
                stop('No mates found despite nonzero number of BND rows in VCF')
            vgr.pair$mix = match(vgr.pair$mix, pix)
            vix = which(1:length(vgr.pair)<vgr.pair$mix )
            vgr.pair1 = vgr.pair[vix]
            vgr.pair2 = vgr.pair[vgr.pair1$mix]

            ## now need to reorient pairs so that the breakend strands are pointing away from the breakpoint

            ## if "first" and "right" then we set this entry "-" and the second entry "+"
            tmpix = vgr.pair1$first & vgr.pair1$right
            if (any(tmpix))
            {
                strand(vgr.pair1)[tmpix] = '-'
                strand(vgr.pair2)[tmpix] = '+'
            }

            ## if "first" and "left" then "-", "-"
            tmpix = vgr.pair1$first & !vgr.pair1$right
            if (any(tmpix))
            {
                strand(vgr.pair1)[tmpix] = '-'
                strand(vgr.pair2)[tmpix] = '-'
            }

            ## if "second" and "left" then "+", "-"
            tmpix = !vgr.pair1$first & !vgr.pair1$right
            if (any(tmpix))
            {
                strand(vgr.pair1)[tmpix] = '+'
                strand(vgr.pair2)[tmpix] = '-'
            }

            ## if "second" and "right" then "+", "+"
            tmpix = !vgr.pair1$first & vgr.pair1$right
            if (any(tmpix))
            {
                strand(vgr.pair1)[tmpix] = '+'
                strand(vgr.pair2)[tmpix] = '+'
            }

            pos1 = as.logical(strand(vgr.pair1)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
            if (any(pos1))
            {
                start(vgr.pair1)[pos1] = start(vgr.pair1)[pos1]-1
                end(vgr.pair1)[pos1] = end(vgr.pair1)[pos1]-1
            }

            pos2 = as.logical(strand(vgr.pair2)=='+') ## positive strand junctions shift left by one (i.e. so that they refer to the base preceding the break for these junctions
            if (any(pos2))
            {
                start(vgr.pair2)[pos2] = start(vgr.pair2)[pos2]-1
                end(vgr.pair2)[pos2] = end(vgr.pair2)[pos2]-1
            }
            ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

            this.inf = values(vgr)[bix[pix[vix]], ]

            if (is.null(this.inf$POS))
                this.inf = cbind(data.frame(POS = ''), this.inf)
            if (is.null(this.inf$CHROM))
                this.inf = cbind(data.frame(CHROM = ''), this.inf)

            if (is.null(this.inf$MATL))
                this.inf = cbind(data.frame(MALT = ''), this.inf)

            this.inf$CHROM = seqnames(vgr.pair1)
            this.inf$POS = start(vgr.pair1)
            this.inf$MATECHROM = seqnames(vgr.pair2)
            this.inf$MATEPOS = start(vgr.pair2)
            this.inf$MALT = vgr.pair2$AL

            ## NOT SURE WHY BROKEN
            ## tmp = tryCatch(cbind(values(vgr)[bix[pix[vix]],], this.inf), error = function(e) NULL)
            ## if (!is.null(tmp))
            ##     values(ra) = tmp
            ## else
            ##     values(ra) = cbind(vcf@fixed[bix[pix[vix]],], this.inf)

            values(ra) = this.inf

            if (is.null(values(ra)$TIER))
                values(ra)$tier = ifelse(values(ra)$FILTER == "PASS", 2, 3) ## baseline tiering of PASS vs non PASS variants
            else
                values(ra)$tier = values(ra)$TIER

            if (!get.loose)
                return(ra)
            else
            {
                npix = is.na(vgr$mix)
                vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

                ## NOT SURE WHY BROKEN
                tmp =  tryCatch( values(vgr)[bix[npix], ],
                                error = function(e) NULL)
                if (!is.null(tmp))
                    values(vgr.loose) = tmp
                else
                    values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])

                return(list(junctions = ra, loose.ends = vgr.loose))
            }
        }
        else
            rafile = read.delim(rafile)
    }

    if (is.data.table(rafile))
    {
        require(data.table)
        rafile = as.data.frame(rafile)
    }

    if (nrow(rafile)==0)
    {
        out = GRangesList()
        values(out) = rafile
        return(out)
    }

    if (snowman) ## flip breaks so that they are pointing away from junction
    {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
    }

    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
    {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]
    }


    if (is.null(rafile$str1))
        rafile$str1 = rafile$strand1

    if (is.null(rafile$str2))
        rafile$str2 = rafile$strand2
    if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
    {
        if (breakpointer)
        {
            rafile$pos1 = rafile$T_BPpos1
            rafile$pos2 = rafile$T_BPpos2
        }

        if (!is.numeric(rafile$pos1))
            rafile$pos1 = as.numeric(rafile$pos1)

        if (!is.numeric(rafile$pos2))
            rafile$pos2 = as.numeric(rafile$pos2)

        ## clean the parenthesis from the string

        rafile$str1 <- gsub('[()]', '', rafile$str1)
        rafile$str2 <- gsub('[()]', '', rafile$str2)

        ## goal is to make the ends point <away> from the junction where - is left and + is right
        if (is.character(rafile$str1) | is.factor(rafile$str1))
            rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))

        if (is.character(rafile$str2) | is.factor(rafile$str2))
            rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))


        if (is.numeric(rafile$str1))
            rafile$str1 = ifelse(rafile$str1>0, '+', '-')

        if (is.numeric(rafile$str2))
            rafile$str2 = ifelse(rafile$str2>0, '+', '-')

        rafile$rowid = 1:nrow(rafile)

        bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0

        rafile = rafile[which(!bad.ix), ]

        if (nrow(rafile)==0)
            return(GRanges())

        seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                    data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

        if (chr.convert)
            seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))

        out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
        out = split(out, out$ra.index)
    }
    else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
    {
        ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
        ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
        out = grl.pivot(GRangesList(ra1, ra2))
    }


    if (keep.features)
        values(out) = rafile[, ]

    if (!get.loose)
        return(out)
    else
        return(list(junctions = out, loose.ends = GRanges()))

    return(out)
}

#' @name jab2json
#' @title jab2json
#'
#' @description
#'
#' Dumps JaBbA graph into json
#'
#' @param jab input jab object
#' @param file output json file
#' @author Marcin Imielinski
jab2json = function(jab, file, maxcn = 100, maxweight = 100)
{

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;

    nodes = jab$segstats %Q% (strand == "+")
    id = rep(1:length(nodes), 2)
    id.type = ifelse(nodes$loose, 'loose_end', 'interval')
    str = ifelse(as.character(strand(jab$segstats))=='+', 1, -1)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = "*",
        endPoint = as.character(end(nodes)),
        title = as.character(1:length(nodes)),
        type = ifelse(nodes$loose, "loose_end", "interval"),
        y = pmin(maxcn, nodes$cn))

    aadj = jab$adj*0
    rix = which(rowSums(is.na(jab$ab.edges[, 1:2, '+']))==0)
    aadj[rbind(jab$ab.edges[rix, 1:2, '+'], jab$ab.edges[rix, 1:2, '+'])] = 1
    ed = which(jab$adj!=0, arr.ind = TRUE)

    if (nrow(ed)>0)
        {
            ed.dt = data.table(
                so = id[ed[,1]],
                so.str = str[ed[,1]],
                si = id[ed[,2]],
                weight = jab$adj[ed],
                title = "",
                type = ifelse(aadj[ed], 'ALT', 'REF'),
                si.str = str[ed[,2]])[, sig := ifelse(so<si,
                                                      paste0(so * so.str, '_', -si*si.str),
                                                      paste0(-si * si.str, '_', so*so.str)
                                                      )][!duplicated(sig), ][, cid := 1:length(weight), ][,
                                                                                                          ":="(so = so*so.str, si = -si*si.str)]
            connections.json = ed.dt[, paste0(
                c("connections: [", paste(
                                        "\t{",
                                        "cid: ", cid,
                                        ", source: ", so,
                                        ", sink:", si,
                                        ", title: ", qw(title),
                                        ", type: ", qw(type),
                                        ", weight: ", pmin(maxweight, weight),
                                        "}",
                                        sep = "",
                                        collapse = ',\n'),
                  "]"),
                collapse = '\n')
                ]
        }

    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                              ", type: ", qw(type),
                              ", strand: ", qw(strand),
                              "}",
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var json = {",
                  paste(
                      c(meta.json,
                      intervals.json,
                      connections.json),
                      collapse = ',\n'
                  ),"}"),
                  sep = "")

    writeLines(out, file)
}



#' @name gr2json
#' @title gr2json
#'
#' @description
#'
#' Dumps GRanges into JSON with metadata features as data points in  "intervals"
#'
#'
#'
#' @param GRange input jab object
#' @param file output json file
#' @author Marcin Imielinski
#' @export
#'
gr2json = function(intervals, file, y = rep("null", length(intervals)), labels = '', maxcn = 100, maxweight = 100)
{

    #' ++ = RL
    #' +- = RR
    #' -+ = LL
    qw = function(x) paste0('"', x, '"')

    ymin = 0;
    ymax = maxcn;

    nodes = intervals
    id = rep(1:length(nodes), 2)

    node.dt = data.table(
        iid = 1:length(nodes),
        chromosome = qw(as.character(seqnames(nodes))),
        startPoint = as.character(start(nodes)),
        strand = as.character(strand(nodes)),
        endPoint = as.character(end(nodes)),
        y = y,
        title = labels)

    oth.cols = setdiff(names(values(nodes)), colnames(node.dt))
    node.dt = as.data.table(cbind(node.dt, values(nodes)[, oth.cols]))

    oth.cols = union('type', oth.cols)
    if (is.null(node.dt$type))
        node.dt$type = 'interval'

    intervals.json = node.dt[, paste0(
        c("intervals: [", paste(
                              "\t{",
                              "iid: ", iid,
                              ", chromosome: ", chromosome,
                              ", startPoint: ", startPoint,
                              ", endPoint: ", endPoint,
                              ", y: ", y,
                              ", title: ", qw(title),
                               ", strand: ", qw(strand),
                              eval(parse(text = ## yes R code making R code making JSON .. sorry .. adding additional columns
                                             paste0("paste0(",
                                                    paste0('", ', oth.cols, ':", qw(', oth.cols, ')', collapse = ','),
                                                    ")", collapse = ''))),
                              "}",
                              sep = "",
                              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    meta.json =
        paste('meta: {\n\t',
              paste(
                  c(paste('"ymin:"', ymin),
                  paste('"ymax:"', ymax)),
                  collapse = ',\n\t'),
              '\n}')

    out = paste(c("var data = {",
                  paste(
                      c(meta.json,
                        intervals.json
                        ),
                      collapse = ',\n'
                  ),"}"),
                sep = "")

    writeLines(out, file)
}

###########################
#' proximity
#'
#' Takes a set of n "query" elements (GRanges object, e.g. genes) and determines their proximity to m "subject" elements
#' (GRanges object, e.g. regulatory elements) subject to set of rearrangement adjacencies (GRangesList with width 1 range pairs)
#'
#' This analysis makes the (pretty liberal) assumption that all pairs of adjacencies that can be linked on a karyograph path are in
#' cis (i.e. share a chromosome) in the tumor genome.
#'
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param ra GRangesList of junctions (each a length 2 GRanges, similar to input to karyograph)
#' @param jab existing JaBbA object (overrides ra input)
#' @param verbose logical flag
#' @param mc.cores how many cores (default 1)
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
#' @return
#' list of n x m sparse distance matrices:
#' $ra = subject-query distance in the rearranged genome for all loci < max.dist in tumor genome
#' $wt = subject-query distance in the reference genome for all loci < max.dist in tumor genome
#' $rel = subject-query distance in ra relative to wild type for above loci
#' NOTE: values x_ij in these matrices should be interpreted with a 1e-9 offset to yield the actual value y_ij
#' i.e. y_ij = x_ij-1e-9, x_ij>0, y_ij = NA otherwise (allows for sparse encoding of giant matrices)
#' @export
############################################
proximity = function(query, subject, ra = GRangesList(), jab = NULL, verbose = F, mc.cores = 1,
                     max.dist = 1e6 ## max distance to store / compute in the output matrix.cores
                     )
{
    if (!is.null(jab))
    {
        ix = which(jab$adj[jab$ab.edges[,1:2,1]]>0)
        if (length(ix)>0)
        {
            ra1 = gr.flipstrand(gr.end(jab$segstats[jab$ab.edges[ix,1,1]], 1, ignore.strand = F))
            ra2 = gr.start(jab$segstats[jab$ab.edges[ix,2,1]], 1, ignore.strand = F)
            ra1 = GenomicRanges::shift(ra1, ifelse(as.logical(strand(ra1)=='+'), -1, 0))
            ra2 = GenomicRanges::shift(ra2, ifelse(as.logical(strand(ra2)=='+'), -1, 0))
            ra = grl.pivot(GRangesList(ra1,ra2))
        }
    }

    if (length(ra)==0)
        return(list())

    if (length(query)==0 | length(subject)==0)
        return(list())

    if (is.null(names(query)))
        names(query) = 1:length(query)

    if (is.null(names(subject)))
        names(subject) = 1:length(subject)

    query.nm = names(query);
    subject.nm = names(subject);

    query = query[, c()]
    subject = subject[, c()]

    query$id = 1:length(query)
    subject$id = 1:length(subject)

    qix.filt = gr.in(query, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
    query = query[qix.filt]

    six.filt = gr.in(subject, unlist(ra)+max.dist) ## to save time, filter only query ranges that are "close" to RA's
    subject = subject[six.filt]

    if (length(query)==0 | length(subject)==0)
        return(list())

    query$type = 'query'
    subject$type = 'subject'

    gr = gr.fix(c(query, subject))

    kg = karyograph(ra, gr)

    ## node.start and node.end delinate the nodes corresponding to the interval start and end
    ## on both positive and negative tiles of the karyograph
    gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;

    ## start and end indices of nodes
    tip = which(strand(kg$tile)=='+')
    tin = which(strand(kg$tile)=='-')
    gr$node.start = tip[gr.match(gr.start(gr,2), gr.start(kg$tile[tip]))]
    gr$node.end = tip[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tip]))]
    gr$node.start.n = tin[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tin]))]
    gr$node.end.n = tin[gr.match(gr.start(gr,2), gr.start(kg$tile[tin]))]

                                        #    gr$node.start = gr.match(gr.start(gr-1,2), gr.start(kg$tile))
                                        #    gr$node.end = suppressWarnings(gr.match(gr.end(gr+1,2), gr.end(kg$tile)))

    ## so now we build distance matrices from query ends to subject starts
    ## and subject ends to query starts

    ## so for each query end we will find the shortest path to all subject starts
    ## and for each query start we will find the shortest.path from all subject ends
    ix.query = which(gr$type == 'query')
    ix.subj = which(gr$type == 'subject')

    node.start = gr$node.start
    node.end = gr$node.end
    node.start.n = gr$node.start.n
    node.end.n = gr$node.end.n

    w = width(kg$tile)

    E(kg$G)$weight = width(kg$tile)[E(kg$G)$to]

    ## ix.query and ix.subj give the indices of query / subject in gr
    ## node.start, node.end map gr to graph node ids
    ##
    ## these matrices are in dimensions of query and subject, and will hold the pairwise distances between
    ##
    D.rel = D.ra = D.ref = D.which = Matrix(data = 0, nrow = length(ix.query), ncol = length(ix.subj))

    ## "reference" graph (missing aberrant edges)
    G.ref = subgraph.edges(kg$G, which(E(kg$G)$type == 'reference'), delete.vertices = F)

    EPS = 1e-9

                                        #    for (i in ix.query)
    tmp = mclapply(ix.query, function(i)
    {
        if (verbose)
            cat('starting interval', i, 'of', length(ix.query), '\n')

        ## D1 = shortest query to subject path, D2 = shortest subject to query path, then take shortest of D1 and D2
        ## for each path, the edge weights correspond to the interval width of the target node, and to compute the path
        ## length we remove the final node since we are measuring the distance from the end of the first vertex in the path
        ## to the beginning of the final vertex

        u.node.start = unique(node.start[ix.subj]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj])

        uix.start = match(node.start[ix.subj], u.node.start)
        uix.end = match(node.end[ix.subj], u.node.end)

        tmp.D1 = (shortest.paths(kg$G, node.end[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(kg$G, node.start[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(kg$G, node.end.n[i], u.node.start, weights = E(kg$G)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(kg$G, node.start.n[i], u.node.end, weights = E(kg$G)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        ix = which(tmp.D<max.dist)
        D.ra[i, ix] = tmp.D[ix]+EPS
        D.which[i, ix] = apply(cbind(tmp.D1[ix], tmp.D2[ix], tmp.D3[ix], tmp.D4[ix]), 1, which.min)

        u.node.start = unique(node.start[ix.subj][ix]) ## gets around annoying igraph::shortest.path issue (no dups allowed)
        u.node.end = unique(node.end[ix.subj][ix])

        uix.start = match(node.start[ix.subj][ix], u.node.start)
        uix.end = match(node.end[ix.subj][ix], u.node.end)

        tmp.D1 = (shortest.paths(G.ref, node.end[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D2 = (shortest.paths(G.ref, node.start[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start[i]])[uix.end]
        tmp.D3 = (shortest.paths(G.ref, node.end.n[i], u.node.start, weights = E(G.ref)$weight, mode = 'out') - w[u.node.start])[uix.start]
        tmp.D4 = (shortest.paths(G.ref, node.start.n[i], u.node.end, weights = E(G.ref)$weight, mode = 'in') - w[node.start.n[i]])[uix.end]
        tmp.D = pmin(tmp.D1, tmp.D2, tmp.D3, tmp.D4)
        D.ref[i, ix] = tmp.D+EPS

        ## if subject and query intersect (on the reference) then we count both RA and Ref distance as 0
        ## (easier to do a simple range query here)
        ix.zero = gr.in(subject[ix], query[i])
        if (any(ix.zero))
        {
            D.ra[i, ix[ix.zero]] = 0
            D.ref[i, ix[ix.zero]] = 0
        }

        D.rel[i, ix] = ((D.ra[i, ix]-EPS) / (D.ref[i, ix]-EPS)) + EPS

        if (verbose)
            cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')

        return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))
    }, mc.cores = mc.cores)

    for (i in 1:length(tmp))
    {
        if (class(tmp[[i]]) != 'list')
            warning(sprintf('Query %s failed', ix.query[i]))
        else
        {
            D.rel = D.rel + tmp[[i]]$D.rel
            D.ra = D.ra + tmp[[i]]$D.ra
            D.ref = D.ref + tmp[[i]]$D.ref
            D.which = D.which + tmp[[i]]$D.which
        }
    }

    ## "full" size matrix
    rel = ra = ref = ra.which =
        Matrix(data = 0, nrow = length(qix.filt), ncol = length(six.filt), dimnames = list(dedup(query.nm), dedup(names(subject.nm))))
    rel[qix.filt, six.filt] = D.rel
    ra[qix.filt, six.filt] = D.ra
    ref[qix.filt, six.filt] = D.ref
    ra.which[qix.filt, six.filt] = D.which

    ## summary is data frame that has one row for each query x subject pair, relative distance, ra distance, and absolute distance
    tmp = which(rel!=0, arr.ind = T)
    colnames(tmp) = c('i', 'j');
    sum = as.data.frame(tmp)

    if (!is.null(query.nm))
        sum$query.nm = query.nm[sum$i]

    if (!is.null(subject.nm))
        sum$subject.nm = subject.nm[sum$j]

    sum$rel = rel[tmp]
    sum$ra = ra[tmp]
    sum$wt = ref[tmp]

    sum = sum[order(sum$rel), ]
    sum = sum[sum$rel<1, ] ## exclude those with rel == 1

    ## reconstruct paths
    vix.query = matrix(NA, nrow = length(qix.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.subject = matrix(NA, nrow = length(six.filt), ncol = 4, dimnames = list(NULL, c('start', 'end', 'start.n', 'end.n')))
    vix.query[qix.filt, ] = cbind(values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start')], values(gr)[ix.query, c('node.start.n')], values(gr)[ix.query, c('node.end.n')])
    vix.subject[six.filt] = cbind(values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start')], values(gr)[ix.subj, c('node.start.n')], values(gr)[ix.subj, c('node.end.n')])

    sum.paths = mapply(function(x, y)
    {
        if ((ra.which[x, y]) == 1)
            get.shortest.paths(kg$G, vix.query[x, 'end'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 2)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
        else if ((ra.which[x, y]) == 3)
            get.shortest.paths(kg$G, vix.query[x, 'end.n'], vix.subject[y, 'start'], weights = E(kg$G)$weight, mode = 'out')$vpath[[1]]
        else if ((ra.which[x, y]) == 4)
            rev(get.shortest.paths(kg$G, vix.query[x, 'start.n'], vix.subject[y, 'end'], weights = E(kg$G)$weight, mode = 'in')$vpath[[1]])
    }, sum$i, sum$j, SIMPLIFY = F)

                                        #    sum$paths = lapply(sum.paths, function(x) x[-c(1, length(x))])
    sum$paths = sum.paths
    sum$ab.edges = lapply(sum.paths, function(p) setdiff(E(kg$G, path = p)$bp.id, NA))

    return(list(sum = sum, rel = rel, ra = ra, wt = ref, G = kg$G, G.ref = G.ref, tile = kg$tile, vix.query = vix.query, vix.subject = vix.subject))
}

##########
#' karyograph
#'
#' @details
#' builds graph from rearrangement breakpoints +/- copy number endpoints
#' used for downstream jbaMIP and karyoMIP functions
#'
#' Input bpp is a GRangesList of signed locus pairs describing aberrant adjacencies.
#' The convention is as follows: Each locus in the input breakpoint pair points to the direction that
#' is being joined by the adjacencies i.e.
#' (-) bp points to "left" or preceding segment
#' (+) bp points to the  "right" or the following segment
#'
#' eg imagine a|bp1|b
#'            c|bp2|d
#'  "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)
#'
#' Input "tile" is a set of intervals whose endpoints are also used to partition the genome prior to the building of the
#' karyograph.
#'
#' Output karyograph connects signed genomic intervals (in a signed tiling of the reference genome) with "aberrant" and "reference" edges.
#' Reference edges connect intervals that are adjacent in the reference genome, and aberrant edges are inferred (upstream
#' of this) through cancer genome paired end analysis.
#' Note that every node, edge, and path in this karyograph has a "reciprocal path"
#'
#' @param junctions GRangesList of junctions, where each item is a length GRanges of signed locations
#' @param tile GRanges optional existing tiling of the genome (eg a copy number segmentation) from which additional segments will be created
#' @return
#'  a list with the following fields
#' $tile = GRanges of length 2*n tiling of the genome corresponding to union of rearrangement breakpoints and copy number endpoints
#' $G = igraph object representing karyograph, here are the edge and vertex features
#'      vertex features: $chrom, $start, $end, $width, $strand, $size, $shape, $border.width, $label, $chrom.ord, $y, $col, $weight
#'      edge features:$ $bp.id, $weight, $from, $to, $col, $type, $line.style, $arrow.shape, $width,
#'      important:  $type specifies which edges are "aberrant" and "reference", $bp.id specifies which input rearrangement (item in junctions)
#'      a given aberrant edge came from (and is NA for reference edges)
#' $adj = 2n x 2n adjacency matrix whose nonzero entries ij show the edge.id in $G
#' $ab.adj = 2n x 2n binary matrix specifying aberrant edges
#' $ab.edges = length(junctions) x {'from', 'to'} x {'+', '-'} mapping junction id's (indices into input junctions lists) to source and sink vertices,
#'             in both orientations
#' @export
############################################
karyograph = function(junctions, ## this is a grl of breakpoint pairs (eg output of ra_breaks(dranger.df) where dranger is df of dranger output)
                      tile = NULL, # pre-existing set of intervals on top of which to build a graph (eg endpoints from a copy number based segmentation)
                      label.edges = FALSE
                      )
{
    require(gplots)
    require(igraph)
    require(Matrix)
    require(data.table)
    require(RColorBrewer)

    if (length(junctions)>0)
    {
        bp.p = grl.pivot(junctions)
        bp1 = gr.end(gr.fix(bp.p[[1]]), 1, ignore.strand = F)
        bp2 = gr.start(gr.fix(bp.p[[2]]), 1, ignore.strand = F)

        if (any(as.logical(strand(bp1) == '*') | as.logical(strand(bp2) == '*')))
            stop('bp1 and bp2 must be signed intervals (i.e. either + or -)')

        if (length(bp1) != length(bp2))
            stop('bp1 and bp2 inputs must have identical lengths')

                                        #    if (sum(width(reduce(bp1))) != sum(width(bp1)) | sum(width(reduce(bp2))) != sum(width(bp2)))
                                        #      stop('bp1 or bp2 cannot have duplicates / overlaps (with respect to location AND strand)')

        values(bp1)$bp.id = 1:length(bp1);
        values(bp2)$bp.id = 1:length(bp1)+length(bp1);

        pgrid = sgn1 = c('-'=-1, '+'=1)[as.character(strand(bp1))]
        sgn2 = c('-'=-1, '+'=1)[as.character(strand(bp2))]

### HACK HACK to force seqlengths to play with each other if malformedo
        tmp.sl = seqlengths(grbind(bp1, bp2))
        tmp.sl.og = tmp.sl
                                        #        tmp.sl = gr2dt(grbind(bp1, bp2))[, max(end, na.rm = TRUE), keyby = seqnames][, sl := pmax(V1+2, tmp.sl[as.character(seqnames)], na.rm = TRUE)][, structure(sl, names = as.character(seqnames))]
        tmp.sl = gr2dt(grbind(bp1, bp2))[, max(end+1, na.rm = TRUE), keyby = seqnames][names(tmp.sl.og), structure(pmax(V1, tmp.sl.og, na.rm = TRUE), names = names(tmp.sl.og))]
        bp1 = gr.fix(bp1, tmp.sl)
        bp2 = gr.fix(bp2, tmp.sl)
                                        # first we tile the genome around the combined breakpoints
    }
    else
    {
        if (is.null(tile))
        {
            tile = si2gr(junctions)
            if (length(tile)==0)
            {
                warning('Empty input given, producing empty output')
                return(NULL)
            }
            A = sparseMatrix(1,1, x = 0, dims = rep(length(tile), 2))
            return(
                list(tile = tile, adj = A,
                     G = graph.adjacency(A), ab.adj = A != 0, ab.edges = NULL, junctions = junctions))
        }

        junctions = GRangesList()
        bp1 = bp2 = GRanges()
    }

    if (!is.null(tile))
    {
        ## find disjoint union of tile and join with gaps
        tile = gr.fix(tile)
        tile = gr.fix(tile, bp1)
        bp1 = gr.fix(bp1, tile) ## argh argh argh .. more pain avoiding hacks
        strand(tile) = '+'
        tile = disjoin(tile)
        tile = sort(c(tile, gaps(tile)))

        ## make sure seqlevels / seqinfo are identical
        if (!identical(sort(seqlevels(tile)), seqlevels(junctions)))
        {
            tile = gr.fix(tile, junctions)
            junctions = gr.fix(junctions, tile)
        }

        if(length(junctions)>0)
        {
            tbp = setdiff(gr.stripstrand(gr.trim(tile, 1)), gr.stripstrand(grbind(bp1, bp2)))
            bp1 = gr.fix(bp1, tbp)
            bp2 = gr.fix(bp2, tbp) ## seqlengths pain
            tbp = gr.fix(tbp, bp1)
        }
        else
            tbp = gr.stripstrand(gr.trim(tile, 1))

        tbp = tbp[start(tbp)!=1]

        if (length(tbp)>0)
            tbp$seg.bp = TRUE
    }
    else
        tbp = NULL;

    if (length(junctions)>0)
        if (length(tbp)>0)
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()], tbp[, c()]))))
        else
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()]))))
    else
        g = gaps(gr.stripstrand(sort(tbp)));

    g = g[strand(g)=='*'];
    strand(g) = '+';

    values(g)$bp.id = NA
    values(g)$seg.bp = NA

    ## combine tiles and find disjoint set
    tile.og = tile
    tile = grbind(bp1, bp2, g, tbp);
    tile = disjoin(gr.stripstrand(tile[order(gr.stripstrand(tile))]))
    strand(tile) = '+'
    tile = gr.fix(tile);
    tile$is.tel = start(tile)==1 | end(tile) == seqlengths(tile)[as.character(seqnames(tile))]
    values(tile)$tile.id = 1:length(tile);

                                        # find "breakpoint" i.e. bp associated intervals, i.e. width 1 intervals that end with a bp1 or bp2 location
    junc.bp = grbind(bp1, bp2)
    junc.bpix = numeric()
    if (length(junc.bp)>0)
        junc.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(junc.bp), start(junc.bp)))

    ## make sure all seqlenths are compatible (so effing annoying)
    tile = gr.fix(tile, bp1)
    tile = gr.fix(tile, bp2)
    bp1 = gr.fix(bp1, tile)
    bp2 = gr.fix(bp2, tile)


    ## also keep track of tbp associatd bp.ix
    all.bp = grbind(bp1, bp2, tbp)
    all.bpix = numeric()

    if (length(all.bp)>0)
        all.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(all.bp), start(all.bp)))

    ## now to build the graph, we would like to fuse all the bp associated intervals with their previous interval
    ## UNLESS they are preceded by another bp associated interval
    ##
    if (length(all.bpix>0))
    {
        to.fuse = all.bpix[which(all.bpix>1 & !((all.bpix-1) %in% all.bpix))]
        end(tile)[to.fuse-1] = end(tile)[to.fuse-1]+1
        tile = tile[-to.fuse]
    }

    if (length(junc.bpix)>0)
    {
        ## we have a partition of genomic segments flanked by tile endpoints and/or ra junctions
        ##
        ## Input junction syntax is interpreted as follows:
        ## a- b+ junctions connect seg ending with position a to seg starting with b+1
        ## a- b- junctions connect seg ending with position a to seg ending with position b (on neg strand)
        ## a+ b+ junctions connect seg starting with position a+1 (on negative strand) to seg starting with position b+1
        ## a+ b- junctions connect seg starting with position a+1 (on negative strand) to seg ending with position b (on neg strand)

                                        # collect all pairwise adjacencies implied by breakpoints
                                        # eg imagine a|bp1|b
                                        #            c|bp2|d
                                        # "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)

        ab.pairs = cbind(
            ifelse(as.logical(strand(bp1)=='+'), gr.match(GenomicRanges::shift(gr.start(bp1), 1), gr.start(tile)),
                   gr.match(gr.start(bp1), gr.end(tile))),
            ifelse(as.logical(strand(bp2)=='+'), gr.match(GenomicRanges::shift(gr.start(bp2), 1), gr.start(tile)),
                   gr.match(gr.start(bp2), gr.end(tile)))
        )


        ## ab.pairs = cbind(
        ##   ifelse(as.logical(strand(bp1)=='+'), match(paste(seqnames(bp1), start(bp1)+1), paste(seqnames(tile), start(tile))),
        ##          match(paste(seqnames(bp1), start(bp1)), paste(seqnames(tile), end(tile)))),
        ##   ifelse(as.logical(strand(bp2)=='+'), match(paste(seqnames(bp2), start(bp2)+1), paste(seqnames(tile), start(tile))),
        ##          match(paste(seqnames(bp2), start(bp2)), paste(seqnames(tile), end(tile))))
        ##   )
        ab.pairs.bpid = bp1$bp.id
        pp = (sgn1*sgn2)>0 & sgn1>0;
        mm = (sgn1*sgn2)>0 & sgn1<0;
        mp = sgn1>0 & sgn2<0
        ab.pairs[pp,1] = -ab.pairs[pp,1] # ++ breakpoints --> (-b)d adjacency
        ab.pairs[mm,2] = -ab.pairs[mm,2] # -- breakpoints --> a(-c) adjacency
        ab.pairs[mp, ] = -ab.pairs[mp, ] # +- breakpoints --> (-b)(-c) adjacency

                                        # clean up adj pairs
                                        # remove any that have crossed a chromosome boundary from their breakpoint
                                        # this will occur in cases of badly formed breakpoint input (eg breakpoints that point outward
                                        # from their telomeres)
        edge.id = rep(1:nrow(ab.pairs), 2)
        ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[,2], -ab.pairs[,1]));
        ab.pairs.bpid = c(ab.pairs.bpid, ab.pairs.bpid)

                                        # build "aberrant" adjacency matrix representing directed graph of edges connecting
                                        # <signed> nodes.
                                        # note: indices of matrix represent edge labels
        adj.ab = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
                        dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
        tmp.ix = cbind(match(as.character(ab.pairs[,1]), rownames(adj.ab)),
                       match(as.character(ab.pairs[,2]), colnames(adj.ab)))
        adj.ab[tmp.ix[!duplicated(tmp.ix), , drop = F]] = ab.pairs.bpid[!duplicated(tmp.ix)]
    }
    else
    {
        ab.pairs.bpid = edge.id = c()
        ab.pairs = matrix(nrow = 0, ncol = 2);
        adj.ab = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                        dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
    }

                                        # build reference adjacency matrix (representing consecutive segments on the reference genome)
                                        # note: indices of matrix represent edge labels
    seg.ix = 1:length(tile)
    ref.pairs = cbind(seg.ix[1:(length(seg.ix)-1)], seg.ix[2:(length(seg.ix))])
                                        # ref.pairs = ref.pairs[ref.pairs[,1]>0 & ref.pairs[,2]!=length(tile), ]
    ref.pairs = ref.pairs[which(as.character(seqnames(tile[ref.pairs[,1]])) == as.character(seqnames(tile[ref.pairs[,2]]))), ]

    if (nrow(ref.pairs)>0)
    {
        edge.id = c(edge.id, max(edge.id) + rep(1:nrow(ref.pairs), 2))
        ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[,2], -ref.pairs[,1])) # reverse ref pairs
        adj.ref = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
        adj.ref[cbind(match(as.character(ref.pairs[,1]), rownames(adj.ref)),
                      match(as.character(ref.pairs[,2]), colnames(adj.ref)))] = nrow(ab.pairs)+1:nrow(ref.pairs)
    }
    else
    {
        adj.ref = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
    }

    ## current tile is partition of genome only in positive orientation + dummy intervals for breakpoints
    ## output tile is forward partition and followed by reverse partition
    ## (this is what is currently referenced by adj.ref and adj.ab)
    ## TODO: clean up this part
    tmp.nm = as.character(c(1:length(tile), -(1:length(tile))))
    tile = c(tile, gr.flipstrand(tile))
    names(tile) = tmp.nm

    ## apply ix to adj.ref and adj.ab, and create "adj" which has union of reference and aberrant junctions
    ## and adj.source which remembers whether edge ij was reference (value = 1) or aberrant (value = 2)
    adj.source = sign(adj.ref)+2*sign(adj.ab)
    adj = sign(adj.ref)+sign(adj.ab)
    edges = which(adj!=0, arr.ind=T) ## num edge x 2 matrix of vertex pairs
    adj[edges] = 1:nrow(edges) ## re number edges across edge set
    rownames(adj) = colnames(adj) = 1:nrow(adj)
    G = graph.adjacency(adj ,weighted = 'edge.ix') ## edge.ix will allow us to match up edges in the adj matrix with edges in the igraph
    node.ind = abs(as.numeric(V(G)$name))

    ## add vertex features including formatting to igraph
    V(G)$chrom = as.character(seqnames(tile))[node.ind]
    V(G)$start = start(tile)[node.ind]
    V(G)$end = end(tile)[node.ind]
    V(G)$width = width(tile)[node.ind]
    V(G)$strand = sign(as.numeric(V(G)$name))
    V(G)$size = 5;
    V(G)$shape= c('rectangle', 'crectangle')[1 + as.numeric(V(G)$strand<0)];
    V(G)$border.width = c(1, 2)[1 + as.numeric(V(G)$strand=='-')] ;
    V(G)$label = paste(V(G)$chrom, ':', round(V(G)$start/1e6,0), '-', round(V(G)$end/1e6,0), sep = '')
    V(G)$label[V(G)$strand<0] = paste(V(G)$chrom, ':', round(V(G)$end/1e6,0), '-', round(V(G)$start/1e6,0), sep = '')[V(G)$strand<0]
    col.map = structure(brewer.master(length(seqlevels(tile))), names = seqlevels(tile))
    V(G)$chrom.ord = levapply(as.numeric(V(G)$start), list(V(G)$chrom), 'rank')
    V(G)$y = V(G)$chrom.ord*30
    V(G)$x = chr2num(V(G)$chrom)*300 + 100*rep(c(0,1), each = length(tile)/2)
    V(G)$col = col.map[V(G)$chrom]

    ## add edge features including formatting to igraph
    E(G)$weight = 1
    E(G)$from = edges[E(G)$edge.ix, 1]
    E(G)$to = edges[E(G)$edge.ix, 2]
    E(G)$col = c(col2hex('gray20'), col2hex('red'))[adj.source[edges[E(G)$edge.ix, ]]]
    E(G)$type = c('reference', 'aberrant', 'aberrant')[adj.source[edges[E(G)$edge.ix, ]]]
    E(G)$line.style = 'SEPARATE_ARROW'
    E(G)$arrow.shape = 'ARROW'
    E(G)$width = 1
    ab.ix = E(G)$type=='aberrant'  ## keep track of bp.id leading to edge
    E(G)$bp.id = NA;
    if (length(ab.pairs.bpid)>0)
        E(G)$bp.id[ab.ix] = ab.pairs.bpid[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
    E(G)$eid = NA; ## what is edge ID??? how is different from edge.ix?
    E(G)$eid[ab.ix] = edge.id[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
    E(G)$eid[!ab.ix] = edge.id[adj.ref[cbind(E(G)$from[!ab.ix], E(G)$to[!ab.ix])]]
    values(tile) = values(tile)[, c('tile.id', 'is.tel')]
    tile$ab.source = 1:length(tile) %in% E(G)$from[ab.ix]
    tile$ab.target = 1:length(tile) %in% E(G)$to[ab.ix]

                                        # important: map input ra to aberrant graph edges, i.e. ab.edges matrix with $from $to and $edge.ix columns
                                        # and one row for each aberrant edge
    ab.edges = array(NA, dim = c(length(junctions), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))
    dupped = duplicated(ab.pairs.bpid)
    ab.edges[,1:2,1] = cbind(match(ab.pairs[!dupped,1], names(tile)), match(ab.pairs[!dupped,2], names(tile)))
    ab.edges[,1:2,2] = cbind(match(ab.pairs[dupped,1], names(tile)), match(ab.pairs[dupped,2], names(tile)))
    ab.edges[,3, 1] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to)) ## must be easier way to perform this taks
    ab.edges[,3, 2] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to))

    if (label.edges & nrow(ab.edges)>0)
    {
        ix = c(ab.edges[,1,1], ab.edges[,2,1], ab.edges[,1,2], ab.edges[,2,2])
        tile$edges.out = tile$edges.in = ''
        tile$edges.in[ix]= sapply(ix,
                                  function(x) {ix = which(adj[,x]!=0); paste(ix, '->', sep = '', collapse = ',')})
        tile$edges.out[ix] = sapply(ix,
                                    function(x) {ix = which(adj[x, ]!=0); paste('->', ix,  sep = '', collapse = ',')})
    }

    return(list(tile = tile, adj = adj, G = G, ab.adj = adj.ab != 0, ab.edges = ab.edges, junctions = junctions))
}
