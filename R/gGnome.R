#' @title gGnome
#'
#' @description
#' Reference-based graph representation of structurally altered genome
#' employing GenomicRanges framework.
#'
#'
#' Copyright (C) 2017  Xiaotong Yao, Marcin Imielinski
#'
#'    This program is free software: you can redistribute it and/or modify
#'    it under the terms of the GNU General Public License as published by
#'    the Free Software Foundation, either version 3 of the License, or
#'    (at your option) any later version.
#'
#'    This program is distributed in the hope that it will be useful,
#'    but WITHOUT ANY WARRANTY; without even the implied warranty of
#'    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#'    GNU General Public License for more details.
#'
#'    You should have received a copy of the GNU General Public License
#'    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#'
#'    Github: https://github.com/mskilab/gGnome
#'    For questions: xiaotong.yao23@gmail.com
#'
#' @import R6
#' @import igraph
#' @import Matrix
#' @import S4Vectors
#' @import gUtils
#' @import gTrack
NULL

## TODO: welcome msg when loading

#' Junctions
#'
#' R6 class extended from GRangesList to represent aberrant genomic SVs
#' with respect to a reference genome
#'
#' @import gUtils
#' @import R6
#' @export
junctions = R6Class("junctions",
                    public = list(
                        refG = "GENOME",
                        initialize = function(grl=NULL, raFile=NULL){
                            GENOME = readRDS(system.file("extdata", "hg19.broad.BSgenome.rds", package="gGnome"))
                            regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
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
                        },
                        length = function(){
                            return(length(private$juncGrl))
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
                        ## given a character vector of col names return that subset of mcols
                        mcols = function(cols="."){
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

length.junctions <- function(junc){
    return(junc$length())
}

#' gGraph
#'
#' the central class for rearrangement graphs
#'
#' @import Matrix
#' @import igraph
#' @import gUtils
#' @import gTrack
#' @import R6
#' @export
gGraph = R6Class("gGraph",
                 public = list(
                     ## public fields
                     ## name = NULL,
                     refG = "GENOME", ## seqinfo of ref genome

                     ## constructor
                     initialize = function(tile=NULL, junctions=NULL,
                                           jabba=NULL, weaver=NULL, prego=NULL,
                                           segs=NULL, es=NULL, ploidy=NULL, purity=NULL){
                         ## control how to construct
                         if (!is.null(segs) &
                                    !is.null(es) &
                                    !is.null(ploidy) &
                                    !is.null(purity)) {
                             private$gGraphFromScratch(segs, es, junctions, ploidy, purity)
                         } else if (!is.null(tile) | !is.null(junctions)) {
                             message("Initializing with 'tile' and 'junctions'")
                             self$karyograph(tile, junctions)
                         } else if (!is.null(jabba)) {
                             ## message("only use 'jabba' or 'weaver' field, not both")
                             message("Reading JaBbA output")
                             self$jabba2gGraph(jabba)
                         } else if (!is.null(weaver)) {
                             ## message("only use 'jabba' or 'weaver' field, not both")
                             message("Reading Weaver output")
                             self$weaver2gGraph(weaver)
                         } else if (!is.null(prego)) {
                             message("Reading Prego output")
                             self$prego2gGraph(prego)
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
                             regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
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
                                                      type=rep("reference", length(circIx)),
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

                         ## resize to width 1, left
                         jUl = unlist(junctions)
                         if (!all(width(jUl))==1){
                             jUl = gr.start(jUl)
                         }
                         jUl = resize(jUl, 2,
                                      fix=ifelse(strand(jUl)=="+",
                                                 "start", "end"))
                         ## ## DONE: remember using split()!!
                         ## mc = mcols(junctions)
                         ## junctions = split(jUl, rep(seq_along(junctions), each=2))
                         ## mcols(junctions) = mc

                         ## start processing
                         ## DONE: write as JaBbA::karyograph() with modifications
                         ## e.g. (30, 2) --> pivot (2, 30)
                         bp.p = split(jUl, rep(1:2, length(junctions)))
                         bp.p = gr.fix(bp.p, get(self$refG))
                         juncTile = c(bp.p[[1]], bp.p[[2]])
                         ## BP 1 and 2, retaining strand-orientation info
                         ## bp1 = gr.start(bp.p[[1]], ignore.strand = T)
                         ## bp2 = gr.start(bp.p[[2]], ignore.strand = T)
                         ## tmpBps = c(bp1, bp2)
                         self$addSegs(juncTile) ## DONE: addSegs is giving dup edges!!!

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
                                               attr = as.list(abEs[,.(cn, type)]))

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
                         refEs[, ":="(cn=tmpDt[from, cn], type="reference")]
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
                     jabba2gGraph = function(jabba, regular.only=F){
                         ptm = proc.time()
                         ## make sure required mcol is filled
                         private$segs = gr.fix(jabba$segstats, get(self$refG))

                         tmpDt = gr2dt(private$segs)
                         tmpDt[, first:=min(start), by=seqnames]
                         tmpDt[, `:=`(first=min(start),last=max(end)), by=seqnames]
                         ## terminal defined as left/right most seg of a chr
                         ## or loose end decoy segs

                         ## DONE: redefine terminal, node wo both of in/out edges
                         private$segs$terminal = tmpDt[, (loose | start==first | end == last)]

                         ## DEBUG: hopefully this will deal with empty edges
                         if (!nrow(jabba$edges)==0){
                             private$es = as.data.table(jabba$edges[,1:4])
                         }

                         private$es[, weight := width(private$segs[from])]
                         private$g = make_directed_graph(
                             t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))

                         ## DONE: get the union of node ix wo in edge and out edge
                         ## SLOW!!!!
                         whichTerminal = private$es[, setxor(from, to)]
                         private$segs$terminal = seq_along(private$segs) %in% whichTerminal

                         cat("segstats done")
                         print(proc.time() - ptm)
                         cat("\n")

                         private$junction = junctions$new(jabba$junctions)

                         cat("junction done")
                         print(proc.time() - ptm)
                         cat("\n")

                         private$abEdges = jabba$ab.edges
                         private$ploidy = jabba$ploidy
                         private$purity = jabba$purity

                         if (regular.only==T){
                             regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
                             v = which(as.vector(seqnames(private$segs)) %in% regularChr)
                             self$subgraph(v, na.rm=F)
                         }
                         cat("subgraph done")
                         print(proc.time() - ptm)
                         cat("\n")

                         return(self)
                     },
                     ## initialize from Weaver result
                     weaver2gGraph = function(weaver){
                         ## TODO: get Weaver done!!!! GEt it done@!!!
                     },

                     ## initialize from Prego result
                     prego2gGraph = function(prego){
                         ## ALERT: I don't check file integrity here!
                         ## first part, Marcin's read_prego
                         res.tmp = readLines(fn)
                         res = structure(lapply(split(res.tmp, cumsum(grepl("edges", res.tmp))),
                                                function(x) read.delim(textConnection(x),
                                                                       strings = F, skip = 1,
                                                                       header = F,
                                                                       col.names = c("node1", "chr1",
                                                                                     "pos1", "node2",
                                                                                     "chr2", "pos2", "cn"))), 
                                         names = gsub(":", "", grep("edges", res.tmp, value = T)))
                         res[[1]]$tag = paste(res[[1]]$node1, ":", res[[1]]$node2, 
                                              sep = "")
                         ## turn into our segstats
                         segstats = GRanges(res[[1]]$chr1,
                                            IRanges(res[[1]]$pos1, 
                                                    res[[1]]$pos2),
                                            strand = "+", cn = res[[1]]$cn,
                                            left.tag = res[[1]]$node1, 
                                            right.tag = res[[1]]$node2)
                         segstats = gr.fix(c(segstats, gr.flipstrand(segstats)))
                         neg.ix = which(strand(segstats) == "-")
                         tag1 = segstats$right.tag
                         tag1[neg.ix] = segstats$left.tag[neg.ix]
                         tag2 = segstats$left.tag
                         tag2[neg.ix] = segstats$right.tag[neg.ix]
                         private$segs = segstats

                         ## adjacency in copy number
                         adj.cn = matrix(0, nrow = length(segstats), ncol = length(segstats), 
                                             dimnames = list(tag1, tag2))
                         adj.cn[cbind(res[[2]]$node1, res[[2]]$node2)] = res[[2]]$cn
                         adj.cn[cbind(res[[2]]$node2, res[[2]]$node1)] = res[[2]]$cn
                         adj.cn[cbind(res[[3]]$node1, res[[3]]$node2)] = res[[3]]$cn
                         adj.cn[cbind(res[[3]]$node2, res[[3]]$node1)] = res[[3]]$cn

                         ## adjacency in edge type
                         adj.type = matrix("", nrow = length(segstats), ncol = length(segstats), 
                                               dimnames = list(tag1, tag2))
                         adj.type[cbind(res[[2]]$node1, res[[2]]$node2)] = "reference"
                         adj.type[cbind(res[[2]]$node2, res[[2]]$node1)] = "reference"
                         adj.type[cbind(res[[3]]$node1, res[[3]]$node2)] = "aberrant"
                         adj.type[cbind(res[[3]]$node2, res[[3]]$node1)] = "aberrant"

                         ## create es
                         ed = as.data.table(which(adj.cn>0, arr.ind=T))
                         colnames(ed) = c("from", "to")                         
                         ed[, ":="(cn = adj.cn[as.matrix(ed[, .(from, to)])],
                                   type = adj.type[as.matrix(ed[, .(from, to)])],
                                   weight = width(segstats[from]))]

                         ## create g
                         g = make_directed_graph(
                             t(as.matrix(ed[,.(from,to)])))
                         private$g = g

                         ## junctions, many of them are copy 0
                         ve = data.table(res$`variant edges`)
                         bp1 = dt2gr(ve[, .(seqnames = chr1, start = pos1, end = pos1)])
                         bp2 = dt2gr(ve[, .(seqnames = chr2, start = pos2, end = pos2)])
                         ## vid1
                         ss = gr.stripstrand(segstats %Q% (strand=="+"))
                         strand(bp1) = ifelse(is.na(match(bp1, gr.end(ss))),
                         bp2$vid = match(bp1, gr.stripstrand(gr.end(segstats, ignore.strand=T)))
                         
                         ## create abEdges
                         abE = array(dim=c(length(junc),3,2),
                                             dimnames=list(NULL,
                                                           c("from", "to", "edge.ix"),
                                                           c("+","-")))
                         if (ed[, any(type=="aberrant")]){
                             ## adding edges to abE
                             abE = self$makeAbEdges()
                         }
                         private$abEdges = abE
                     }

                     ## For DEBUG purpose only
                     ## TODO: delete this in official release
                     hardset = function(segs=NULL, es=NULL, juncs=NULL, g=NULL){
                         if (!is.null(segs)) private$segs = segs
                         if (!is.null(es)) private$es = es
                         if (!is.null(juncs)) private$junction = junctions$new(juncs)
                         if (!is.null(g)) private$g = g
                         return(self)
                     },

                     ## public methods
                     ## I/O
                     print = function(){
                         cat('A gGraph object.\n')
                         cat('Based on reference genome: ')
                         cat(self$refG)
                         cat('\n\n')
                         cat('Total segmentation:')
                         cat(length(private$segs %Q% (loose==F & strand=="+")))
                         cat('\n\n')
                         cat('Junction counts:\n')
                         print(private$es[, table(type)/2])
                     },
                     plot = function(pad=1e3, colorful=FALSE){
                         td = self$gGraph2gTrack()
                         if (colorful == TRUE){
                             ## ASSUMPTION: segs are sorted by strand first
                             td@data[[1]]$segment = rep(LETTERS[1:(length(private$segs)/2)], 2)
                             td$gr.colorfield = "segment"
                             td@data[[1]]$lbl = rep(LETTERS[1:(length(private$segs)/2)], 2)
                             td$xaxis.ticklen = 0.25
                             td$xaxis.chronly = T
                             td$xaxis.unit = 1e6
                             td$xaxis.suffix = "Mb"
                             td$xaxis.round = 2
                             td$xaxis.interval = 1e6
                             td$xaxis.cex.tick = 0.75
                             td$xaxis.cex.label = 0
                             td$yaxis.pretty = 2
                             td$yaxis.cex = 0.75
                             td$lwd.border = 2
                             td$gr.labelfield = "lbl"
                             td$gr.cex.label = 1.2
                             td$sep.lwd = 0.5
                         }
                         ## DONE: plot all junctions on top
                         win = streduce(private$segs) + pad
                         ## decide X gap on the fly
                         plot(td, win, links = private$junction$grl)
                     },
                     layout = function(){
                         ## TODO: return the plot value
                         vcolor = ifelse(strand(private$segs)=="+", "salmon", "skyblue")
                         c3 = setNames(skitools::brewer.master(n = 3, palette = "Set1"),
                                       nm = c("aberrant", "loose", "reference"))
                         ed = private$es
                         ed[, ecolor := c3[type]]

                         plot.igraph(private$g,
                                     ## layout
                                     layout = layout_with_gem,
                                     ## vertex pars
                                     vertex.size=log(private$segs$cn,1.4), vertex.color= vcolor,
                                     vertex.shape="circle", vertex.label.cex = 0.75,
                                     vertex.frame.color=NA, vertex.label.color = "black",

                                     ## edge pars
                                     edge.lty=3, edge.arrow.width=0.3, edge.arrow.size=0.25,
                                     edge.width=log(private$es$cn, base = 7)+0.3,
                                     edge.color=ed$ecolor)
                     },
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
                         return(private$parition$no)
                     },
                     ##

                     gGraph2gTrack = function(seg.col){
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
                         ed = ed[!is.na(from) & !is.na(to)]
                         ## TODO: handle so/si-less edges, omit for now

                         ## set segment apperances
                         ## if loose, change its cn to slightly higher than it's incident node
                         if (any(ss$loose==T)){
                             lid = which(ss$loose)
                             ## find partner indices for loose ends
                             pid = sapply(lid,
                                          function(i) ed[from==i | to==i,
                                                         ifelse(from==i, to, from)],
                                          simplify=T)
                             ss$cn[lid] = ss$cn[pid]*1.2
                         }
                         ## col, border, ywid
                         ss$col = ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5))
                         ss$border = ifelse(ss$loose, ss$col, alpha("black", 0.5))
                         ss$ywid = ifelse(ss$loose, 0.001, 0.8)

                         gt = gTrack(ss, y.field="cn", edges=ed, name="CN", angle=0)
                         return(gt)
                     },
                     gGraph2json = function(file=NULL,
                                            maxcn=100,
                                            maxweight=100,
                                            trim = TRUE ## trim will only output seqnames that are relevant to the plot
                                            ){
                         system(paste('mkdir -p', file))
                         system(sprintf('cp -r %s %s',
                                        paste0(system.file("extdata", "gTrack.js/complete-genome-interval-graph", package = 'gGnome'), '/*'),
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
                                 c(paste0(qw("connections"),": ["), paste(
                                                         "\t{",
                                                         qw("cid"), ":", cid,
                                                         ifelse(is.na(so), "", paste0(",",qw("source"),":")),
                                                         ifelse(is.na(so), "", so),
                                                         ifelse(is.na(si), "", paste0(",",qw("sink"),":")),
                                                         ifelse(is.na(si), "", si),
                                                         ",", qw("title"), ":", qw(title),
                                                         ",", qw("type"), ":", qw(type),
                                                         ",", qw("weight"), ": ", pmin(maxweight, weight),
                                                         "}",
                                                         sep = "",
                                                         collapse = ',\n'),
                                   "]"),
                                 collapse = '\n')]

                         }

                         ## processing intervals
                         intervals.json = node.dt[, paste0(
                             c(paste0(qw("intervals"),": ["), paste(
                                                   "\t{",
                                                   qw("iid"), ":", iid,
                                                   ",", qw("chromosome"), ":", chromosome,
                                                   ",", qw("startPoint"), ":", startPoint,
                                                   ",", qw("endPoint"), ":", endPoint,
                                                   ",", qw("y"), ":", y,
                                                   ",", qw("title"), ":", qw(title),
                                                   ",", qw("type"), ":", qw(type),
                                                   ",", qw("strand"), ":", qw(strand),
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
                             paste(paste0('\t',qw("metadata"),': [\n'),
                                   chrs[, paste("\t\t{",
                                                qw("chromosome"),":", qw(seqnames),
                                                ",", qw("startPoint"),":", 1,
                                                ",", qw("endPoint"), ":", seqlengths,
                                                ",", qw("color"),
                                                ":", qw(substr(tolower(brewer.master( max(.I), 'BrBG' )), 1, 7)), " }",
                                                collapse=",\n",
                                                sep="")],
                                   '\n]')

                         ## assembling the JSON
                         out = paste(c("{",
                                       paste(
                                           c(meta.json,
                                             intervals.json,
                                             connections.json),
                                           collapse = ',\n'
                                       ),"}"),
                                     sep = "")

                         ## DONE: remove any NA. Not legal.
                         ## if (!is.null(file)){
                         ##     writeLines(out, file)
                         ## }
                                        #                         return(out)

                         writeLines(out, paste0(file, '/data.json'))
                         message(sprintf('Wrote JSON file of gGraph to %s/index.html', file))
                     },

                     ## self-annotating functions
                     hydrogenBonds = function(){
                         ## collapse +/- strand
                         ss = unique(gr.stripstrand(private$segs))
                         idss = match(gr.stripstrand(private$segs), ss)
                         if (!all(table(idss)==2)){
                             stop("Malformed object. Suggest creation again.")
                         }
                         tmpDt = data.table(ssid = seq_along(ss))
                         tmpDt[, ":="(n1 = which(idss==ssid)[1],
                                      n2 = which(idss==ssid)[2]), by=ssid]
                         hydrogenBs = tmpDt[, .(from = n1, to = n2,
                                                   type="hydrogen", cn=0, weight=0)]
                         return(hydrogenBs)
                     },

                     ## dicing up the graph
                     components = function(mc.cores=1){
                         ## create a sticky graph where pairs of +/- are connected by hydro edges
                         stickyG = private$g
                         hB = self$hydrogenBonds()
                         ## update es and g
                         stickyG = add_edges(stickyG, t(as.matrix(hB[, .(from, to)])))

                         private$partition = components(stickyG)
                         ## merge +/- complements into 1

                         ## DONE!!! faster subgraph construction
                         ## split nodes/edges by membership!!! rather than subgraph!!!
                         ## define a compound gGraph class for holding a series of them
                         nComp = private$partition$no ## total N of parts

                         allComponents = lapply(1:nComp,
                                  function(i){
                                      v = which(private$partition$membership==i)
                                      thisComp = self$subgraph(v, na.rm=F, mod=F)
                                      return(thisComp)
                                  })

                         return(allComponents)
                     },
                     melt = function(){
                         ## TODO: think if I really need this
                         "Return a pair of ssgGraph"
                     },
                     ## DONE:
                     ## if na.rm==F, balanced graph's subgraph should always be balanced!!!!!
                     subgraph = function(v=numeric(0), na.rm=T, mod=T){
                         "Given a numeric vector of vertices, change this gGraph to its subgraph consists only these vertices."
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

                             ## DONE: also recover v's missing reverse complements
                             hB = self$hydrogenBonds()
                             vid = sort(unique(c(v, hB[from %in% v, to], hB[to %in% v, from])))

                             ## get the subgraph
                             newSegs = private$segs[vid]
                             newId = setNames(seq_along(vid), vid)

                             if (na.rm==T){
                                 newEs = private$es[from %in% vid & to %in% vid,]
                             } else {
                                 ## DONE: if na.rm==FALSE, which is the default when calling
                                 ## from bGraph, turn the NA edges to new loose ends
                                 ## except for new "telomere"
                                 newEs = private$es[from %in% vid | to %in% vid,]

                                 newLooseIn = newEs[!from %in% vid]
                                 newLooseOut = newEs[!to %in% vid]
                                 ## only mod newEs when there are extra loose ends
                                 if (nrow(newLooseIn)>0){
                                     ## create new id mapping
                                     newLoose = rbind(newLooseIn, newLooseOut)
                                     looseId = setNames(1:nrow(newLoose)+length(newId), c(newLooseIn$from, newLooseOut$to))
                                     ## create new segs
                                     lin = gr.start(private$segs[newLooseIn$to], ignore.strand=F)
                                     lout = gr.end(private$segs[newLooseOut$from], ignore.strand=F)
                                     ## append new loose end nodes
                                     newL = c(lin, lout)
                                     newL$loose=TRUE
                                     newSegs = c(newSegs, newL)
                                     ## append new loose IDs
                                     newId = c(newId, looseId)
                                 }
                             }

                             newEs[, ":="(from = newId[as.character(from)],
                                          to = newId[as.character(to)])]
                             ## mark extra loose edges
                             newEs[from > length(vid) | to > length(vid), type := "loose"]

                             jIdx = which(grl.in(private$junction$grl, newSegs, only=T))
                             newJuncs = private$junction[unique(jIdx)]

                             if (mod==T){
                                 private$gGraphFromScratch(segs=newSegs,
                                                           es=newEs,
                                                           junc=newJuncs,
                                                           ploidy=private$ploidy,
                                                           purity=private$purity)
                                 return(self)
                             } else {
                                 out = gGraph$new(segs=newSegs,
                                                  es=newEs,
                                                  junctions=newJuncs,
                                                  ploidy=private$ploidy,
                                                  purity=private$purity)
                                 return(out)
                             }
                         } else {
                             stop("Invalid input.")
                         }
                     },
                     trim = function(gr=NULL){
                         ## DONE
                         ## if input gr is super set of private$segs, do nothing!
                         ## Only returning new obj
                         "Given a GRanges, return the trimmed subgraph overlapping it."
                         if (is.null(gr))
                             return(self)

                         gr = gr.fix(gr, get(self$refG))
                         gr = gr.stripstrand(gr)
                         if (!isDisjoint(gr))
                             gr = gr.reduce(gr)

                         ## TODO: causing weird error with seqlengths
                         ## if (length(setdiff(streduce(private$segs), gr))==0)
                         ##     return(self)

                         v = which(gr.in(private$segs, gr))
                         sg = self$subgraph(v, na.rm=F, mod=F)
                         ## if (length(v)<=2)
                         ##     return(sg)
                         ## DONE: resolve the edge case where gr is contained in single node

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
                                        gr[length(gr), -lastCol]
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
                     simplify = function(){
                         ## TODO: rm 0 copy edges, non reg chr
                     },
                     getSeqInfo = function(){
                         as.data.table(attributes(seqinfo(get(self$refG))))
                     },
                     makeAbEdges = function(){
                         ## DONE: derive abEdges from junction
                         if (length(private$junction$grl)==0){
                             return(
                                 array(dim=c(0,3,2),
                                       dimnames=list(NULL,
                                                     c("from", "to", "edge.ix"),
                                                     c("+","-")))
                             )
                         } else {
                             ## based on junctions, get
                             junc = private$junction$grl

                             abEdges = array(dim=c(length(junc),3,2),
                                             dimnames=list(NULL,
                                                           c("from", "to", "edge.ix"),
                                                           c("+","-")))

                             ## find coresponding edge.ix for abe
                             lbp = unlist(junc) ## ASSUMPTION: junctions are width 1, marking the left nt of a bp
                             lbp$jix = rep(seq_along(junc), each=2) ## which junction?
                             lbp$bix = rep(1:2, length(junc)) ## breakpoint 1 or 2?
                             ## tell if bps are width 2 or 1
                             if (all(width(lbp)==1)){
                                 rbp = lbp %+% 1
                             } else {
                                 lbp = gr.start(lbp)
                                 rbp = lbp %+% 1
                             }

                             jUl = c(lbp, rbp) ## left bp, right bp
                             jUl$side = rep(c("left","right"), each=length(lbp))

                             seg = private$segs %Q% (loose==F)## ASSUMPTION: segs are sorted by loose first, loose ends at the tail
                             xStart = gr.stripstrand(jUl) %*% gr.stripstrand(gr.start(seg[,c()]))
                             xEnd = gr.stripstrand(jUl) %*% gr.stripstrand(gr.end(seg[,c()]))

                             ## put together
                             tmpDt = data.table(rbind(as.data.frame(mcols(xStart)), as.data.frame(mcols(xEnd))))
                             if (nrow(tmpDt)==0) return(abEdges)
                             ## only want to retain the jix that have 8/4 rows in this dt
                             ## aka, none of their bps maps to the middle of a segment
                             tmpDt = tmpDt[jix %in% tmpDt[, names(which(table(jix)==8))]]
                             tmpDt[, ori := as.vector(strand(jUl[query.id]))]
                             tmpDt[, str := as.vector(strand(seg[subject.id]))]

                             froms = tmpDt[
                                 side==ifelse(ori=="-", "left", "right"),
                                 .SD,
                                 by=bix][
                                 ori!=str,
                                 .(froms=subject.id, strand=ifelse(bix==1,"+","-")),
                                 by=jix]

                             tos = tmpDt[
                                 side==ifelse(ori=="-", "left", "right"),
                                 .SD,
                                 by=bix][
                                 ori==str,
                                 .(tos=subject.id, strand=ifelse(bix==1,"-","+")),
                                 by=jix]

                             allE = merge(froms, tos, by=c("jix", "strand"))
                             setkey(private$es, from, to)
                             private$es[, edge.ix := .I]
                             allE[, edge.ix := private$es[.(froms, tos), edge.ix], by=.I]

                             allE = unique(allE) ## ALERT: strange edge cases where same bp mapped twice!!!
                             jix = unique(allE[, jix])
                             ## fill in the right blanks
                             abEdges[jix, , "+"] = allE[jix %in% jix & strand=="+", c(froms, tos, edge.ix)]
                             abEdges[jix, , "-"] = allE[jix %in% jix & strand=="-", c(froms, tos, edge.ix)]

                             return(abEdges)
                         }
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

                         ## DONE: what to do when win is larger than segs?????
                         ## ans: return self
                         if (length(setdiff(streduce(private$segs), win))==0)
                             return(self)

                         ## overlapping window and segs, removing loose ends
                         interGr = gr.findoverlaps(private$segs, win, ignore.strand=ignore.strand)
                         lid = which(private$segs$loose==T)
                         interGr = interGr %Q% (!query.id %in% lid)
                         qix = interGr$query.id

                         if (is.null(k)){
                             ## DONE!!!
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
                             ## DONE: connect with dist
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
                             return(self$subgraph(kNeighbors, mod=F)) ## not garanteed size to scale
                         }
                     },

                     dist = function(gr1, gr2,
                                     matrix=T, EPS=1e-9,
                                     include.internal=TRUE, ## consider bp within feature "close"
                                     directed=FALSE, ## if TRUE, only consider gr1-->gr2 paths
                                     verbose=FALSE){
                         "Given two GRanges, return pairwise shortest path distance."
                         ## DONE
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
                     ## NOW TODO:
                     proximity = function(query, subject,
                                          verbose=F, mc.cores=1,
                                          max.dist=1e6){

                         ## TODO:
                         adj = self$getAdj()
                         ix = which(adj[private$abEdges[,1:2,1]]>0)
                         if (length(ix)>0) {
                             ra1 = gr.flipstrand(
                                 gr.end(private$segs[private$abEdges[ix,1,1]],
                                        width=1, ignore.strand = F))
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


                         ## ALERT! TODO! There are NAs in node.start/end etc
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
                         ## DONE: use adj to calc if every segment is balanced on both sides
                         adj = self$getAdj()
                         whichTerminal = which(private$segs$terminal==T)
                         whichNa = which(is.na(private$segs$cn))
                         validTerminal = setdiff(whichTerminal, whichNa)

                         ## balanced on both sides for non-terminal nodes
                         middleTrue = (Matrix::colSums(adj)[-c(whichTerminal, whichNa)] ==
                                       private$segs[-c(whichTerminal, whichNa)]$cn) &
                             (private$segs[-c(whichTerminal, whichNa)]$cn ==
                              Matrix::rowSums(adj)[-c(whichTerminal, whichNa)])
                         ## balanced on either end for terminal nodes
                         tCsum = colSums(adj)[validTerminal]
                         tRsum = rowSums(adj)[validTerminal]
                         terminalConSide = ifelse(tCsum==0, tRsum, tCsum)
                         terminalTrue = terminalConSide == private$segstats[validTerminal]$cn
                         return(all(middleTrue) & all(terminalTrue))
                     },
                     isDoubleStrand = function(){
                         ## DONE: test if segs come in +/- pairs
                         identical((ss %Q% (strand=="-"))[, c()],
                                   gr.flipstrand(ss %Q% (strand=="+"))[, c()])
                     },
                     getLooseEnds = function(){
                         ## TODO: return all loose ends as a GRanges
                     },

                     walk = function(){
                         ## enumerate all paths through the graph even if it is not balanced

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
                     gGraphFromScratch = function(segs, es, junc, ploidy, purity){
                         private$segs = segs
                         private$es = es
                         ## relabel the terminals!
                         ## whichTerminal = private$es[, setxor(from, to)]
                         ## private$segs$terminal = seq_along(private$segs) %in% whichTerminal
                         private$g = make_directed_graph(
                             t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))

                         if (is.null(junc)){
                             warning("Junctions not provided. Inferring from edges.")
                             hB = self$hydrogenBonds()
                             map = hB[, c(setNames(from, to), setNames(to, from))]

                             bp.from =
                                 gr.end(private$segs[private$es$from], ignore.strand=FALSE)[,c()]
                             bp.from[which(strand(bp.from)=="-")] =
                                 bp.from[which(strand(bp.from)=="-")] %-% 1
                             bp.from = gr.flipstrand(bp.from)
                             ## the other end
                             bp.to =
                                 gr.start(private$segs[private$es$to], ignore.strand=FALSE)[,c()]
                             bp.to[which(strand(bp.from)=="+")] =
                                 bp.to[which(strand(bp.from)=="+")] %-% 1

                             junc = grl.pivot(GRangesList(bp.from, bp.to))
                             ## TODO: write grl.duplicated function
                             junc = junc[!grl.duplicated(junc)]
                         }
                         private$junction = junctions$new(junc)

                         private$abEdges = self$makeAbEdges()
                         private$ploidy = ploidy
                         private$purity = purity
                     }## ,
                     ## ## collapse strand info
                     ## getSs = function(){
                     ##     "Return simple segs, with names, tile.id, is.tel, ab.source, ab.target."

                     ##     ## ## TODO: think about how did he plot loose ends!!!
                     ##     ## ## processing nodes
                     ##     ## ## reduce strand
                     ##     ## ## remove loose nodes
                     ##     ## oid = gr2dt(private$segs)[, which(strand == "+" & loose==F)]
                     ##     ## ## ori ind of rev comps
                     ##     ## rid = gr2dt(private$segs)[, which(strand == "-" & loose==F)]

                     ##     ## ## single strand
                     ##     ## ss = gr.stripstrand(private$segs[oid])
                     ##     ## newMap = match(gr.stripstrand(private$segs), ss)

                     ##     ## ## ori ix of loose nodes
                     ##     ## lid = which(private$segs$loose==T)

                     ##     ## ## processing edges
                     ##     ## ed = private$es
                     ##     ## ed[,":="(soStr = as.character(strand(private$segs[from])),
                     ##     ##          siStr = as.character(strand(private$segs[to])))]
                     ##     ## edByType = by(ed, ed$type, function(x) x)

                     ##     ## ## see which of the ab edges are "+"
                     ##     ## abe = edByType$aberrant
                     ##     ## if (!is.null(abe)){
                     ##     ##     abe[, key := paste(from, to, sep="_")]
                     ##     ##     setkey(abe, "key")
                     ##     ##     ## info in ab.edges field
                     ##     ##     posAbEd = as.data.table(private$abEdges[,1:2,"+"])[!is.na(from+to)]
                     ##     ##     abe = abe[posAbEd[, paste(from, to, sep="_")],-c("key")]
                     ##     ## }

                     ##     ## ## put 3 back together
                     ##     ## ed = rbindlist(list(edByType$reference[soStr=="+"],
                     ##     ##                     edByType$loose[soStr=="+"],
                     ##     ##                     abe))

                     ##     ## ## processing edges, cont.
                     ##     ## if (nrow(ed)>0){
                     ##     ##     ed[, ":="(newFr = newMap[from], newTo = newMap[to])]
                     ##     ## }

                     ## }
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
                         ## DONE: make igraph plot
                         return(self$layout())
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
                     },
                     parts = function(){
                         if (is.null(private$partition))
                             tmp = self$components()
                         return(private$partition)
                     }
                 )
                 )
#'
#'
#'
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
#'
components.gGraph <- function(gGraph){
    ## input must be a gGraph!
    if (!is(gGraph, "gGraph")){
        stop("Invalid input.")
    }
    return(gGraph$components())
}

length.gGraph <- function(gGraph){
    ## input must be a gGraph!
    if (!is(gGraph, "gGraph")){
        stop("Invalid input.")
    }
    if (is.null(gGraph$parts)){
        cs = gGraph$components()
    }
    return(gGraph$parts$no)
}

#' Descendant of gGraph class, where junction balance restraint must be met all the time
#'
#' @import R6
#'
bGraph = R6Class("bGraph",
                 inherit = gGraph,
                 public = list(
                     ## overwrite constructor: restrict about junction balance
                     initialize = function(gG=NULL, jabba=NULL){
                         if (!is.null(gG)){
                             if (is(gG, "gGraph")){
                                 private$gGraphFromScratch(gG$segstats, gG$edges, gG$junctions, gG$ploidy, gG$purity)
                                 return(self)
                             } else {
                                 stop("Invalid input gG.")
                             }
                         } else if (!is.null(jabba)) {
                             regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
                             allRegChr = all(
                                 as.vector(seqnames(unlist(jabba$junctions))) %in% regularChr
                             )
                             self$jabba2gGraph(jabba=jabba, allRegChr)
                             if (self$isJunctionBalanced()){
                                 return(self)
                             } else {
                                 stop("Invalid input gG.")
                             }
                         } else {
                             self$nullGGraph()
                             return(self)
                         }
                     },

                     print = function(){
                         cat('A bGraph object.\n')
                         cat('Based on reference genome: ')
                         cat(self$refG)
                         cat('\n\n')
                         cat('Total non-loose segmentation:')
                         cat(length(private$segs %Q% (loose==F & strand=="+")))
                         cat('\n\n')
                         cat('Junction counts:\n')
                         print(private$es[, table(type)/2])
                     },

                     ## DONE: the bGraph created from jab different that from gGraph!!!
                     subgraph = function(v=numeric(0), na.rm=F, mod=F){
                         if (mod == T){
                             super$subgraph(v, na.rm=F, mod=mod)
                             return(self)
                         } else {
                             out = super$subgraph(v, na.rm=F, mod=mod)
                             return(out)
                         }
                     },

                     dsb = function(){},
                     del = function(){},
                     tDup = function(){},
                     invs = function(){},

                     ## decompose graph into all possible haplotypes
                     walk = function(outdir="tmp.walk",
                                     max.iteration = Inf,
                                     mc.cores = 1,
                                     verbose = T,
                                     nsolutions = 100,
                                     tilim = 100){
                         "Give all the possible multiset of walks that can be represented by this graph."
                         ## TODO: do components one by one
                         if (length(self)>1){
                             return("Hasn't implemented this")
                         }

                         segs = private$segs ## ALERT: don't mod the field
                         ss = unique(gr.stripstrand(segs))
                         idss = match(gr.stripstrand(segs), ss)
                         if (!all(table(idss)==2)){
                             stop("Malformed object. Suggest creation again.")
                         }

                         whichSeg = which(segs$loose==F) ## non-loose id in segs
                         A = self$getAdj()[whichSeg, whichSeg]
                         ## get incidence matrix
                         ## vertices x edges
                         ed0 = data.table(which(A!=0, arr.ind=T))
                         colnames(ed0) = c("from", "to")
                         ed0[, ":="(cn = A[ as.matrix(ed0[, .(from, to)]) ],
                                    type="nonslack")]
                         ## uniquely map rev-comp edges
                         ## SHIT, this is wrong
                         ## 2->4 and 4->2 are two different edges!
                         ## ed0[, ":="(fss=idss[from], tss=idss[to])]
                         ## ed0[, ":="(mx=max(fss, tss), mn=min(fss, tss)), by=1:nrow(ed0)]
                         ## ed0[, eclass := as.numeric(as.factor(paste(mn, mx, sep=".")))]
                         ## TODO: find the right way to map edges


                         ifl = Matrix::colSums(A) ## net in flux for every seg
                         ofl = Matrix::rowSums(A) ## net out flux for every seg

                         ## TODO: what to do about self-edge!
                         ## ANS: nothing, but maybe useful to place a micro value
                         ## for the reconstruction of paths

                         ## in flux < cn, need 5' slack
                         slack.in = segs[whichSeg]$cn - ifl
                         ## out flux < cn, need 3' slack
                         slack.out = segs[whichSeg]$cn - ofl

                         ## data.table to keep track of new slacks added
                         slacks = data.table(vid = whichSeg[c(which(slack.in>0), which(slack.out>0))])
                         slacks[, ":="(type = rep(c("slack.in", "slack.out"), each=nrow(slacks)/2),
                                       cn = c(slack.in[slack.in>0], slack.out[slack.out>0]),
                                       ssid = idss[vid],
                                       strand=as.vector(strand(segs[vid])))]
                         slacks[, side := ifelse(
                         (strand=="+" & type=="slack.in") | (strand=="-" & type=="slack.out"),
                         "left", "right")]
                         slacks[, eclass := as.numeric(as.factor(paste(ssid, side, sep=".")))]

                         ## extend ed to contain slacks
                         mx.eclass = ed0[, max(eclass)]
                         ed = rbind(ed0[, .(from, to, cn, etype=type, eclass)],
                                    slacks[type=="slack.in", .(from=NA, to=vid, cn, etype=type, eclass=eclass+mx.eclass)],
                                    slacks[type=="slack.out", .(from=vid, to=NA, cn, etype=type, eclass=eclass+mx.eclass)])

                         ## TODO: assemble h, input to karyoMIP -- e, e.ij, B, eclass, etype
                         ## ASSUMPTION: private$segs is sorted by loose then strand
                         ii = c(ed[etype=="nonslack", c(from, to)],
                                ed[etype=="slack.in", to],
                                ed[etype=="slack.out", from])
                         jj = c(ed[, rep(which(etype=="nonslack"), 2)],
                                ed[, which(etype=="slack.in")],
                                ed[, which(etype=="slack.out")])
                         xx = c(rep(c(-1, 1), each=nrow(ed0)),
                                rep(1, nrow(slacks)/2),
                                rep(-1, nrow(slacks)/2))
                         B = sparseMatrix(i = ii, j = jj, x = xx, dims=c(nrow(A), nrow(ed)))

                         h = list(e = ed[, cn],
                                  e.ij = as.matrix(ed[, .(from, to)]),
                                  B = B,
                                  eclass = ed[, eclass],
                                  etype = ed[, ifelse(grepl("slack", etype),
                                                      "slack", "nonslack")])

                         ## compute convex basis of B
                         K = convex.basis(B)
                         prior = rep(1, ncol(K))

                         browser()
                         ## TODO: convert karyoMIP solution to gWalks object
##                         is.cyc = Matrix::colSums(K[h$etype == 'slack', ])==0 & Matrix::colSums((Bc %*% K)!=0)==0
                         karyo.sol = karyoMIP(K, h$e, h$eclass,
                                              nsolutions = nsolutions,
                                              tilim = tilim,
                                              cpenalty = 1/prior)

                         ## if (saveAll){
                         ##     saveRDS(karyo.sol, "temp.walk/allSol.rds")
                         ## }
                         kag.sol = karyo.sol[[1]]
                         p = karyoMIP.to.path(kag.sol, K, h$e.ij, segs[whichSeg])
                         p$paths = mclapply(p$paths, as.numeric, mc.cores=mc.cores)

                         ## construct gWalks as result
                         gw = gWalks$new(segs=segs[whichSeg], paths=p$paths, isCyc=p$is.cyc, cn = p$cn)
                         browser()
                         return(gw)
                     },

                     ## TODO: hurestic walk decomposition
                     ## new idea: if we assign weight
                     #' @name walk2
                     #' @title walk2
                     #' @description
                     #'
                     #' Computes greedy collection (i.e. assembly) of genome-wide walks (graphs and cycles) by finding shortest paths in JaBbA graph.
                     #'
                     #' @param jab JaBbA object
                     #' #
                     #' @return GRangesList of walks with copy number as field $cn, cyclic walks denoted as field $is.cycle == TRUE, and $wid (width) and $len (segment length) of walks as additional metadata
                     #' @export
                     walk2 = function(verbose = FALSE)
                     {
                         cn.adj = self$getAdj()
                         adj = as.matrix(cn.adj)
                         adj.new = adj*0
                         adj[which(adj!=0, arr.ind = TRUE)] = width(private$segstats)[which(adj!=0, arr.ind = TRUE)[,2]] ## make all edges a large number by default
                         if (verbose)
                             message('Setting edge weights to destination widths for reference edges and 1 for aberrant edges')

                         ab.edges = rbind(private$abEdges[,1:2, 1], private$abEdges[,1:2, 2])
                         ab.edges = ab.edges[rowSums(is.na(ab.edges))==0, ]
                         adj[ab.edges] = sign(cn.adj[ab.edges]) ## make ab.edges = 1
                         adj[is.na(adj)] = 0
                         cn.adj[which(is.na(cn.adj))] = 0

                         G = graph.adjacency(adj, weighted = 'weight')

                         ## define ends not using degree (old method) but using either telomeres or loose ends
                         ## (otherwise lots of fake ends at homozygous deleted segments)
                         ss = gr2dt(private$segstats)[ , vid:= 1:length(seqnames)]
                         ss[loose == TRUE, is.end := TRUE]
                         ss[loose == FALSE, is.end := 1:length(loose) %in% c(which.min(start), which.max(end)), by = list(seqnames, strand)]
                         ends = which(ss$is.end)

                         ## sanity check
                         unb = which(!ss$is.end & rowSums(private$getAdj(), na.rm = TRUE) != colSums(private$getAdj(), na.rm = TRUE))

                         if (length(unb)>0)
                         {
                             message(sprintf('JaBbA model not junction balanced at %s non-ends! Adding these to "ends"', length(unb)))
                             ends = c(ends, unb)         ## shameless HACK ... TOFIX
                         }
## BOOKMARK TODO
                         ## ends = which(degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
                         i = 0
                         D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

                         ## sort shortest paths
                         ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]

                         maxrow = length(ends)*max(cn.adj[ends, ends], na.rm = TRUE)
                         vpaths = rep(list(NA), maxrow)
                         epaths = rep(list(NA), maxrow)
                         cns = rep(NA, maxrow)


                         #' first peel off "simple" paths i.e. zero degree
                         #' ends with >0 copy number
                         psimp =  which(degree(G, mode = 'out')==0 & degree(G, mode = 'in')==0 & jab$segstats$cn>0)
                         i = 0
                         if (length(psimp)>0)
                         {
                             vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
                             epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA)) ## there is no "edge" associated with a zero total degree node
                             cns[1:length(psimp)] = jab$segstats$cn[psimp]
                             i = length(psimp)
                         }

                         ## now iterate from shortest to longest path
                         ## peel that path off and see if it is still there ..
                         ## and see if it is still there
                         ## peel off top path and add to stack, then update cn.adj

                         tile.map = gr2dt(jab$segstats)[, .(id = 1:length(tile.id), tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
                         rtile.map = gr2dt(jab$segstats)[, .(id = 1:length(tile.id), tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
                         setkey(tile.map, id)
                         setkey(rtile.map, tile.id)


                         while (nrow(ij)>0)
                         {
                             if (verbose)
                                 message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' ends to resolve' )
                             i = i+1
                             vpaths[[i]] = p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])
                             epaths[[i]] = cbind(p[-length(p)], p[-1])
                             cns[i] = min(cn.adj[epaths[[i]]])

                             rvpath = rtile.map[list(tile.map[list(vpaths[[i]]), -rev(tile.id)]), id]
                             repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
                             plen = length(rvpath)
                             hplen = floor(length(rvpath)/2)
                             ## ## (awkward) check for palindromicity for odd and even length palindromes
                             ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
                             ##         palindromic = TRUE
                             ## else
                             ## {
                             ##     vpaths[[i+1]] = rvpath
                             ##     epaths[[i+1]] = repath
                             ##     cns[i+1] = cns[i]
                             ##     palindromic = FALSE
                             ## }
                             palindromic = TRUE ## set to true while we "figure things out"

                             #' so now we want to subtract that cn units of that path from the graph
                             #' so we want to update the current adjacency matrix to remove that path
                             #' while keeping track of of the paths on the stack
                             cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]

                             if (!palindromic) ## update reverse complement unless palindromic
                                 cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]

                             if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
                             {
                                 message('backtracking ...') ## maybe we got stuck in a quasi-palindrome and need to backtrack
                                 cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]
                                 if (!palindromic) ## update reverse complement unless palindromic
                                     cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
                                 i = i-1
                                 ij = ij[-1, ]
                             }
                             else ## continue, reduce
                             {
                                 adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
                                 if (!palindromic)
                                     adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

                                 ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                                 ##     browser()

                                 to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,drop = FALSE]
                                 if (!palindromic) ## update reverse complement
                                     to.rm = rbind(to.rm, epaths[[i+1]][which(cn.adj[epaths[[i+1]]]==0), ,drop = FALSE])

                                 if (nrow(to.rm)>0)
                                 {
                                     adj[to.rm] = 0
                                     G = graph.adjacency(adj, weighted = 'weight')
                                     new.ends = setdiff(which(
                                     (degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
                                     & degree(G)>0), ends)

                                     if (length(new.ends)>0)
                                         browser()

                                     D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
                                     ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
                                 }
                                 else
                                     ij = ij[-1, ]

                                 if (!palindromic) ## increase extra counter to account for reverse complement
                                     i = i+1
                             }
                         }
                         if (verbose)
                             message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )

                         vpaths = vpaths[1:i]
                         epaths = epaths[1:i]
                         cns = cns[1:i]


                         #' how to deal with cycles?
                         #' first peel off simple (ie one off) cycles
                         vcycles = rep(list(NA), maxrow)
                         ecycles = rep(list(NA), maxrow)
                         ccns = rep(NA, maxrow)

                         csimp = which(diag(cn.adj)!=0)
                         i = 0
                         if (length(csimp)>0)
                         {
                             vcycles[1:length(csimp)] = split(csimp, 1:length(csimp))
                             ecycles[1:length(csimp)] = lapply(csimp, function(x) cbind(x, x))
                             ccns[1:length(csimp)] = diag(cn.adj)[csimp]
                             cn.adj[cbind(csimp, csimp)] = 0
                             adj[cbind(csimp, csimp)] = 0
                             i = length(csimp)

                             for (j in 1:length(csimp))
                                 adj.new[ecycles[[j]]] = adj.new[ecycles[[j]]] + ccns[j]
                         }

                         ## sort shortest paths and find which connect a node to its ancestor (i.e. is a cycle)
                         .parents = function(adj)
                         {
                             tmp = apply(adj, 2, function(x) which(x!=0))
                             ix = which(sapply(tmp, length)>0)
                             if (length(ix)>0)
                             {
                                 parents = rbindlist(lapply(ix, function(x) data.table(x, tmp[[x]])))
                                 setnames(parents, c('node', 'parent'))
                                 setkey(parents, node)
                             } else {
                                 parents = data.table(node = 0, parent = NA)
                                 setkey(parents, node)
                             }
                         }

                         parents = .parents(adj)

                         #' then find paths that begin at a node and end at (one of its) immediate upstream neighbors
                         #' this will be a path for whom col index is = parent(row) for one of the rows
                         G = graph.adjacency(adj, weighted = 'weight')
                         D = shortest.paths(G, mode = 'out', weight = E(G)$weight)

                         ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ]
                         ij[, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]

                         ## now iterate from shortest to longest path
                         ## peel that path off and see if it is still there ..
                         ## and see if it is still there

                         ## peel off top path and add to stack, then update cn.adj
                         while (nrow(ij)>0)
                         {
                             if (verbose)
                                 message('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
                             i = i+1
                             p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])
                             vcycles[[i]] = p
                             ecycles[[i]] = cbind(p, c(p[-1], p[1]))
                             ccns[i] = min(cn.adj[ecycles[[i]]])

                             rvcycle = rtile.map[list(tile.map[list(vcycles[[i]]), -rev(tile.id)]), id]
                             recycle = cbind(rvcycle[-length(rvcycle)], rvcycle[-1])
                             clen = length(rvcycle)
                             hclen = floor(length(rvcycle)/2)
                             ## ## (awkward) check for palindromicity for odd and even length palindromes
                             ## if (all((vcycles[[i]]==rvcycle)[c(1:hclen,(clen-hclen+1):clen)]))
                             ##         palindromic = TRUE
                             ## else
                             ## {
                             ##     vcycles[[i+1]] = rvcycle
                             ##     ecycles[[i+1]] = recycle
                             ##     ccns[i+1] = ccns[i]
                             ##     palindromic = FALSE
                             ## }
                             palindromic = TRUE ## set to true while we "figure things out"

                             #' so now we want to subtract that cn units of that path from the graph
                             #' so we want to update the current adjacency matrix to remove that path
                             #' while keeping track of of the cycles on the stack
                             cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
                             if (!palindromic) ## update reverse complement unless palindromic
                                 cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]

                             if (!all(cn.adj[ecycles[[i]]]>=0))
                             {
                                 message('backtracking')
                                 cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]+ccns[i]
                                 if (!palindromic) ## update reverse complement unless palindromic
                                     cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]+ccns[i+1]
                                 i = i-1
                                 ij = ij[-1, ]
                             }
                             else
                             {
                                 adj.new[ecycles[[i]]] = adj.new[ecycles[[i]]] + ccns[i]

                                 if (!palindromic)
                                     adj.new[ecycles[[i+1]]] = adj.new[ecycles[[i+1]]] + ccns[i]

                                 ## intermediate cross check
                                 ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                                 ##     browser()

                                 to.rm = ecycles[[i]][which(cn.adj[ecycles[[i]]]==0), ,drop = FALSE]

                                 if (!palindromic) ## update reverse complement
                                     to.rm = rbind(to.rm, ecycles[[i+1]][which(cn.adj[ecycles[[i+1]]]==0), ,drop = FALSE])

                                 if (nrow(to.rm)>0)
                                 {
                                     adj[to.rm] = 0
                                     parents = .parents(adj)
                                     G = graph.adjacency(adj, weighted = 'weight')
                                     D = shortest.paths(G, mode = 'out', weight = E(G)$weight)
                                     ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
                                 }
                                 else
                                     ij = ij[-1, ]

                                 if (!palindromic) ## increase extra counter to account for reverse complement
                                     i = i+1
                             }
                         }
                         if (verbose)
                             message('Cycle peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )


                         if (i>0)
                         {
                             vcycles = vcycles[1:i]
                             ecycles = ecycles[1:i]
                             ccns = ccns[1:i]
                         }
                         else
                         {
                             vcycles = NULL
                             ecycles = NULL
                             ccns = NULL
                         }

                         vall = c(vpaths, vcycles)
                         eall = c(epaths, ecycles)
                         ecn = c(cns, ccns)

                         tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
                         ix = which(rowSums(is.na(tmp[, 1:2]))==0)

                         if (length(ix)>0)
                             adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
                         else
                             adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
                         vix = munlist(vall)
                         paths = split(jab$segstats[vix[,3]], vix[,1])
                         values(paths)$ogid = 1:length(paths)
                         values(paths)$cn = ecn[as.numeric(names(paths))]
                         values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
                         values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
                         values(paths)$numsegs = elementNROWS(paths)
                         values(paths)$wid = sapply(lapply(paths, width), sum)

                         check = which((adj.new - jab$adj) !=0, arr.ind = TRUE)
                         if (length(check)>0)
                             stop('Alleles do not add up to marginal copy number profile!')
                         else if (verbose)
                             message('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')

                         ## match up paths and their reverse complements
                         psig = lapply(paths, function(x) ifelse(as.logical(strand(x)=='+'), 1, -1)*x$tile.id)
                         psig.flip = sapply(psig, function(x) -rev(x))

                         unmix = data.table(
                             ix = 1:length(paths),
                             mix = match(sapply(psig, paste, collapse = ','), sapply(psig.flip, paste, collapse = ',')))[, pos := 1:length(mix)<mix][order(!pos), ]
                         setkey(unmix, ix)
                         unmix[is.na(mix), pos := TRUE] ## if we have paths with no reverse complement i.e. NA mix, then call "+" for now

                         remix = rbind(
                             unmix[pos == TRUE, ][, id := 1:length(ix)],
                             unmix[list(unmix[pos == TRUE, mix]), ][, id := 1:length(ix)][!is.na(ix), ]
                         )


                         paths = paths[remix$ix]
                         names(paths) = paste(remix$id, ifelse(remix$pos, '+', '-'), sep = '')
                         values(paths)$id = remix$id
                         values(paths)$str = ifelse(remix$pos, '+', '-')

                         if (length(setdiff(values(paths)$ogid, 1:length(paths))))
                             message('Warning!!! Some paths missing!')

                         return(paths)
                     }
                     ),
                 private = list(

                 ),
                 active = list())

## Utilities
ul = function(x, n=6){
    n = pmin(pmin(dim(x)), n)
    return(x[1:n, 1:n])
}

#' getPloidy
#'
#' @export
getPloidy = function(segs){
    if (!is(segs, "GRanges")) stop("Not a GRanges!")
    if (!isDisjoint(segs)) stop("Must be disjoint!")
    if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case=T))==0) print("No copy number (cn) column!")

    cn = mcols(segs)[, cnix]
    wd = width(segs)

    pl = weighted.mean(cn, wd)
    return(pl)
}

#' grl.duplicated
#'
#' @export
grl.duplicated = function(x, as.tuple=FALSE, mc.cores=1){
    if (!is(x, "GRangesList")) stop("Not a GRangesList!")

    ## only recurrent
    dt = data.table(ii = seq_along(x), elen = elementNROWS(x), duplicated=FALSE)
    dt[, tlen := nrow(.SD), by=elen]
    ##dt[tlen>1, x[[ii]], by=elen]

    trueId = mclapply(dt[tlen>1, setNames(unique(elen), unique(elen))],
                      function(el){
                          iis = dt[elen==el, ii]
                          ix = combn(iis, 2)
                          thisIdIx = apply(ix, 2, function(iix){
                              if (!as.tuple){
                                  iid = identical(sort(x[iix][1]), sort(x[iix][2]))
                              } else {
                                  iid = identical(x[iix][1], x[iix][2])
                              }
                              if (iid) return(max(iix))
                              else return(NULL)
                          })
                      },
                      mc.cores = mc.cores)

    trueId = unlist(trueId)
    if (length(trueId)>0) set(dt, trueId, 'duplicated', TRUE)
    return(dt[, duplicated])
}

##################################
#' @name vaggregate
#' @title vaggregate
#'
#' @description
#' same as aggregate except returns named vector
#' with names as first column of output and values as second
#'
#' Note: there is no need to ever use aggregate or vaggregate, just switch to data.table
#'
#' @param ... arguments to aggregate
#' @return named vector indexed by levels of "by"
#' @author Marcin Imielinski
#' @export
##################################
vaggregate = function(...)
{
    out = aggregate(...);
    return(structure(out[,ncol(out)], names = do.call(paste, lapply(names(out)[1:(ncol(out)-1)], function(x) out[,x]))))
}

######################################################
#' @name mmatch
#' @title mmatch
#'
#' @description
#' match rows of matrix A to matrix B
#'
#' @param A query matrix k1 x n
#' @param B subject matrix k2 x n
#' @param dir 1
#' @return length k1 vector specifying first row of B matching row i of A
#' @export
#' @author Marcin Imielinski
######################################################
mmatch = function(A, B, dir = 1)
{
    SEP = ' ';
    Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
    Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))

    return(match(Atxt, Btxt))
}

#' @name alpha
#' @title alpha
#' @description
#' Give transparency value to colors
#'
#' Takes provided colors and gives them the specified alpha (ie transparency) value
#'
#' @author Marcin Imielinski
#' @param col RGB color
#' @keywords internal
#' @export
alpha = function(col, alpha)
{
    col.rgb = col2rgb(col)
    out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
}

####################################
#' .e2class
#'
#' edge to contig class conversion
#'
#' given matrix K of k contigs over e edges, each belonging to cardinality 1 or cardinality 2 equivalence classes,
#' assigns id's to equivalent contigs
#'
####################################
.e2class = function(K, eclass)
{
    eclass = factor(as.character(eclass))

    if (length(eclass)!=nrow(K))
        stop('eclass must be of the same length as number of rows in K')

    eclass = factor(as.character(eclass))
    class.count = table(eclass);

    if (any(class.count)>2)
        stop('Edge equivalence classes can have at most 2 members')

    biclasses = names(class.count)[class.count==2];  # classes with two edges

    if (length(biclasses)>0)
    {
                                        # edge class rotation matrix
        R = diag(!(eclass %in% biclasses));  ## edges belonging to classes of cardinality 1 are on the diagonal

        ix = matrix(unlist(split(1:length(eclass), eclass)[biclasses]), ncol = 2, byrow = T); # index pairs corresponding to edges in biclasses
        R[ix[, 1:2]] = 1
        R[ix[, 2:1]] = 1

        Kr = R %*% K
        eix = mmatch(t(Kr), t(K))
        eix[is.na(eix)] = 0
        pairs = t(apply(cbind(1:length(eix), eix), 1, sort))
        pairs = pairs[!duplicated(pairs) & rowSums(pairs==0)==0, , drop = FALSE]

        kclass = rep(NA, ncol(K))
        kclass[pairs[,1]] = 1:nrow(pairs);
        kclass[pairs[,2]] = 1:nrow(pairs);
        kclass[is.na(kclass)] = nrow(pairs) + (1:sum(is.na(kclass)))
    }
    else
        kclass = 1:ncol(K)

    return(kclass)
}
##################################
#' convex.basis
#'
#' Outputs a matrix K of the convex basis of matrix A
#'
#' i.e. each column x = K[,i] is a minimal solution (with respect to sparsity) to
#' Ax = 0, x>=0
#'
#' exclude.basis =  0, 1 matrix of dimension k x ncol(A) specifying k sparsity patterns that we would
#' like to exclude from the convex.basis.  This can speed up computation since any non-negative
#' combinations of vectors that satisfy an exclusion property will also be excludable, and thus
#' we can remove such vectors as soon as we detect them..
#'
#' exclude.range = 9, 1 matrix of dimension k x nrow(A) specifying k sparsity patterns that we would like
#' exclude, but these are specified in terms of the range of abs(A) .. i.e. we want to exclude all
#' basis vectors v such that nz(exclude.ranges[i, ]) C  nz(abs(A)*v)) for some pattern i.  Again
#' any non-neg linear comb of any intermediate-basis vector that satisfies this property will satisfy it,
#' as a result we can exclude these vectors when we see them.
#'
#'
#'
##################################
convex.basis = function(A, interval = 80, chunksize = 100, exclude.basis = NULL, exclude.range = NULL, maxchunks = Inf,
                        verbose = F)
{
    ZERO = 1e-8;
    remaining = 1:nrow(A);
    iter = 0;
    i = 0;
                                        #    order = c()
    numelmos = c()
    K_i = I = as(diag(rep(1, ncol(A))), 'sparseMatrix');
                                        #    A_i = as(A %*% K_i, 'sparseMatrix');
    K_i = I = diag(rep(1, ncol(A)))
    A_i = A %*% K_i

    if (!is.null(exclude.basis))
    {
        exclude.basis = sign(exclude.basis)
        exclude.basis = exclude.basis[rowSums(exclude.basis)>0, ]
        if (nrow(exclude.basis) == 0)
            exclude.basis = NULL
    }

    if (!is.null(exclude.range))
    {
        exclude.range = sign(exclude.range)
        exclude.range = exclude.range[rowSums(exclude.range)>0, ]
        if (nrow(exclude.range) == 0)
            exclude.range = NULL
    }

                                        # vector to help rescale matrix (avoid numerical issues)
    mp  = apply(abs(A), 1, min); # minimum value of each column
    mp[mp[ZERO]] = 1; # columns with zero minimum get scale "1"

    st = Sys.time()
                                        # iterate through rows of A, "canceling" them out
    while (length(remaining)>0)
    {
        if (nrow(K_i)==0 | ncol(K_i)==0) ## TODO figure out why we have to check this so many times
            return(matrix())

        iter = iter+1;
        K_last = K_i;

        if (verbose)
            print(Sys.time() - st)

        if (verbose)
            cat('Iter ', iter, '(of',  nrow(A_i),  ') Num basis vectors: ', nrow(K_i), " Num active components: ", sum(Matrix::rowSums(K_i!=0)), "\n")

        i = remaining[which.min(Matrix::rowSums(A_i[remaining,, drop = FALSE]>=ZERO)*Matrix::rowSums(A_i[remaining,, drop = FALSE]<=(-ZERO)))]  # chose "cheapest" rows

        remaining = setdiff(remaining, i);
                                        #        order = c(order, i);

        zero_elements = which(abs(A_i[i, ]) <= ZERO);
        K_i1 = K_last[zero_elements, , drop = FALSE];  ## K_i1 = rows of K_last that are already orthogonal to row i of A
        K_i2 = NULL; ## K_i1 = will store positive combs of K_last rows that are orthogonal to row i of A (will compute these below)

        pos_elements = which(A_i[i, ]>ZERO)
        neg_elements = which(A_i[i, ]<(-ZERO))

        if (verbose)
            cat('Iter ', iter, " Row ", i, ":", length(zero_elements), " zero elements ", length(pos_elements), " pos elements ", length(neg_elements), " neg elements \n")

        if (length(pos_elements)>0 & length(neg_elements)>0)
            for (m in seq(1, length(pos_elements), interval))
                for (l in seq(1, length(neg_elements), interval))
                {
                    ind_pos = c(m:min(c(m+interval, length(pos_elements))))
                    ind_neg = c(l:min(c(l+interval, length(neg_elements))))

                    indpairs = cbind(rep(pos_elements[ind_pos], length(ind_neg)),
                                     rep(neg_elements[ind_neg], each = length(ind_pos))); # cartesian product of ind_pos and ind_neg
                    pix = rep(1:nrow(indpairs), 2)
                    ix = c(indpairs[,1], indpairs[,2])
                                        #                coeff = c(-A_i[i, indpairs[,2]], A_i[i, indpairs[,1]])  ## dealing with Matrix ghost
                    coeff = c(-A_i[i, ][indpairs[,2]], A_i[i, ][indpairs[,1]])  ##
                    combs = sparseMatrix(pix, ix, x = coeff, dims = c(nrow(indpairs), nrow(K_last)))
                    combs[cbind(pix, ix)] = coeff;

                    H = combs %*% K_last;

                                        # remove duplicated rows in H (with respect to sparsity)
                    H = H[!duplicated(as.matrix(H)>ZERO), ];

                                        # remove rows in H that have subsets in H (with respect to sparsity) ..
                    if ((as.numeric(nrow(H))*as.numeric(nrow(H)))>maxchunks)
                    {
                        print('Exceeding maximum number of chunks in convex.basis computation')
                        stop('Exceeding maximum number of chunks in convex.basis computation')
                    }
                    keep = which(Matrix::colSums(sparse_subset(abs(H)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))<=1) # <=1 since every H is its own subset
                    H = H[keep, , drop = FALSE]

                                        # remove rows in H that have subsets in K_i2
                    if (!is.null(K_i2))
                        if (nrow(K_i2)>0)
                        {
                            if ((as.numeric(nrow(K_i2))*as.numeric(nrow(H)))>maxchunks)
                            {
                                print('Exceeding maximum number of chunks in convex.basis computation')
                                stop('Exceeding maximum number of chunks in convex.basis computation')
                            }
                            keep = which(Matrix::colSums(sparse_subset(abs(K_i2)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                            H = H[keep, , drop = FALSE]
                        }

                                        # remove rows in H that have subsets in K_i1
                    if (!is.null(K_i1))
                        if (nrow(K_i1)>0)
                        {
                            if ((as.numeric(nrow(K_i1))*as.numeric(nrow(H)))>maxchunks)
                            {
                                print('Exceeding maximum number of chunks in convex.basis computation')
                                stop('Exceeding maximum number of chunks in convex.basis computation')
                            }
                            keep = which(Matrix::colSums(sparse_subset(abs(K_i1)>ZERO, abs(H)>ZERO, chunksize = chunksize, quiet = !verbose))==0)
                            H = H[keep, , drop = FALSE]
                        }

                                        # maintain numerical stability
                    if ((iter %% 10)==0)
                        H = diag(1/apply(abs(H), 1, max)) %*% H

                                        #                K_i2 = rBind(K_i2, H)
                    K_i2 = rbind(K_i2, as.matrix(H))
                }

                                        #        K_i = rBind(K_i1, K_i2)
        K_i = rbind(K_i1, K_i2) ## new basis set

        if (nrow(K_i)==0)
            return(matrix())

        if (!is.null(exclude.basis)) ## only keep vectors that fail to intersect all vectors "exclude" in matrix
        {
            if ((as.numeric(nrow(exclude.basis))*as.numeric(nrow(K_i)))>maxchunks)
            {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
            }
            keep = Matrix::colSums(sparse_subset(exclude.basis>0, K_i>ZERO))==0
            if (verbose)
                cat('Applying basis exclusion and removing', sum(keep==0), 'basis vectors\n')
            K_i = K_i[keep, , drop = F]
        }

        if (!is.null(exclude.range)) ## only keep vectors that fail to intersect all vectors "exclude" in matrix
        {
            A_i_abs = abs(A) %*% t(K_i)
            if ((as.numeric(nrow(exclude.range))*as.numeric*ncol(A_i_abs))>maxchunks)
            {
                print('Exceeding maximum number of chunks in convex.basis computation')
                stop('Exceeding maximum number of chunks in convex.basis computation')
            }
            keep = Matrix::colSums(sparse_subset(exclude.range>0, t(A_i_abs), quiet = !verbose))==0
            if (verbose)
                cat('Applying range exclusion and removing', sum(keep==0), 'basis vectors\n')
            K_i = K_i[keep, , drop = F]
        }

        A_i = A %*% t(K_i)
    }

    return(t(K_i))
}

#'
#' hGraph: haplotype gGraph
#' has to be balanced and one unique path through all segs
#'
#' @import gUtils
#' @import gTrack
#' @import igraph
#' @export
hGraph = R6Class("hGraph",
                 inherit=bGraph,
                 public = list(),
                 private = list(),
                 active = list())


#'
#' gWalks: subclass to gGraph
#'
#' @import gUtils
#' @import gTrack
#' @import R6
#'
#' @export
gWalks = R6Class("gWalks",
                 public=list(
                     refG = "GENOME",

                     initialize = function(grl=NULL, segs=NULL, paths=NULL, isCyc=NULL, cn=NULL, str=NULL){
                         if (!is.null(segs)){
                             private$gwFromScratch(segs, paths, isCyc, cn)
                         } else if (!is.null()) {
                             self$grl2gw(grl)
                         } else {
                             self$nullGWalks()
                         }
                     },
                     ## TODO: construct null gWalks
                     nullGWalks = function(){},

                     ## TODO: gw2bg convert to a list of bGraphs
                     gw2gg = function(){
                         ## TODO:
                         ## proceed only if it passes strand pair test
                         ## if (self$isStrandPaired()){}

                         if (is.null(private$es)) private$es = self$path2edges()

                         amp = rep(private$cn, elementNROWS(private$paths))
                         cns = table(rep(unlist(private$paths), amp))

                         private$segs$cn = as.numeric(cns[as.character(seq_along(private$segs))])

                         pl = getPloidy(private$segs)

                         gg = gGraph$new(segs = private$segs,
                                         es = private$es,
                                         ploidy = pl,
                                         purity = 1)
                         return(gg)
                     },

                     gw2grl = function(){
                         segs = private$segs
                         segs$tile.id = rep(LETTERS[1:23][1:(length(private$segs)/2)], 2)
                         grl = lapply(private$paths,
                                      function(pt){
                                          return(segs[pt])
                                      })
                         grl = GRangesList(grl)
                         mcols(grl) = data.table(contig.cn = private$cn, is.cyc = private$isCyc)
                         private$grl = grl
                         return(grl)
                     },
                     grl2gw = function(){

                     },
                     gw2gTrack = function(mc.cores=1, colorful=FALSE){
                         gts = gTrack(private$grl, draw.path=T)
                         if (colorful) {

                         }
                         gts$name = paste0("cycle=", private$isCyc, "\ncn=", private$cn)

                         return(gts)
                     },

                     ## TODO: helper function to turn paths into edges
                     path2edges = function(mc.cores=1){

                         ## whenever this function runs, it will assign result to
                         ## private$es, which will be refreshed to NULL whenever
                         ## a modifying action happens
                         es = do.call('rbind',
                                      mclapply(1:length(private$paths),
                                      function(i){
                                          thisPath = private$paths[[i]]
                                          ll = length(thisPath)
                                          thisFrom = thisPath[1:(ll-1)]
                                          thisTo = thisPath[2:ll]

                                          if (private$isCyc[[i]]){
                                              thisFrom[ll] = thisPath[ll]
                                              thisTo[ll] = thisPath[1]
                                          }

                                          thisWeight = width(private$segs)[thisFrom]
                                          thisEs = data.table(from = thisFrom,
                                                              to = thisTo,
                                                              cn = private$cn[i],
                                                              type = "unknown",
                                                              weight = thisWeight)

                                      },
                                      mc.cores=mc.cores))

                         ## if same edges shows up more than once, dedup and populate cn
                         es[, tmp.id := paste(from, to, sep="-")]
                         setkey(es, "tmp.id")
                         es[, cn := sum(cn), by=tmp.id]
                         es = es[!duplicated(tmp.id), .(from, to, cn, type, weight)]

                         ## finally, identify type of edge by segs
                         es[, ":="(fromChr = as.vector(seqnames(private$segs[from])),
                                       fromStr = as.vector(strand(private$segs[from])),
                                       fromStart = start(private$segs[from]),
                                       fromEnd = end(private$segs[from]),
                                       toChr = as.vector(seqnames(private$segs[to])),
                                       toStr = as.vector(strand(private$segs[to])),
                                       toStart = start(private$segs[to]),
                                       toEnd = end(private$segs[to]))]

                         es[(fromChr!=toChr | fromStr!=toStr), type := "aberrant"]

                         es[(fromChr==toChr & fromStr==toStr &
                                 fromStr=="+" & (fromEnd - toStart)==-1), type := "reference"]
                         es[(fromChr==toChr & fromStr==toStr &
                                 fromStr=="-" & (fromStart - toEnd)==1), type := "reference"]
                         es[type=="unknown", type := "aberrant"]

                         es = es[, .(from, to, cn, type, weight)]
                         return(es)
                     },
                     subset = function(ix){
                         ## TODO:subset and overwrite `[`
                     },

                     print = function(){
                         str = paste0("gWalks:\n",
                                      "\t", length(private$paths), " contigs\n",
                                      "\t", length(private$segs), " ranges\n")
                         cat(str)
                         return(str)
                     },
                     len = function(){
                         return(length(private$paths))
                     },
                     metaCols = function(){
                         mdt = data.table(isCyc = private$isCyc,
                                          cn = private$cn,
                                          str = private$str)
                         return(mdt)
##                         paths = sapply(private$paths, paste, collapse=",")
                     },
                     window = function(ix = NULL, pad=1e3){
                         if (is.null(ix))
                             return(trim(streduce(private$segs) + pad))

                         ix = unlist(private$paths[ix])
                         return(trim(streduce(private$segs[ix]) + pad))
                     },
                     plot = function(){
                         td = self$gw2gTrack()
                         win = self$window()
                         plot(td, win)
                     },

                     ## tests
                     isStrandPaired = function(){
                         ## check point 1
                         if (!all(table(gr.match(segs3, segs3))==2)) return(FALSE)

                         ## check point 2


                         return(TRUE)
                     }
                 ),
                 private=list(
                     segs = NULL,
                     paths = NULL,
                     isCyc = NULL,
                     cn = NULL,
                     str = NULL,
                     es = NULL,
                     grl = NULL,

                     ## private constructor
                     gwFromScratch = function(segs, paths=NULL, isCyc=NULL, cn=NULL, str=NULL){
                         ## segs must be a GRanges
                         if (!is(segs, "GRanges")) stop("segs needs to be a GRanges.")

                         ## ... and disjoint
                         if (!isDisjoint(segs)) stop("segs must be disjoint.")

                         private$segs = segs

                         ## paths must be a list of numeric vectors
                         if (!is.list(paths) |
                             !all(sapply(paths, is.numeric)) |
                             any(sapply(paths, max)>length(segs)) |
                             any(sapply(paths, min)<1))
                             stop("paths must be list of indices, within 1:length(segs)")

                         private$paths = paths

                         if (!is.null(isCyc) & is.logical(isCyc) & length(isCyc)==length(paths)) private$isCyc = isCyc
                         else private$isCyc = rep(FALSE, length(private$paths))

                         if (!is.null(cn) & is.numeric(cn) & length(cn)==length(paths)) private$cn = cn
                         else private$cn = rep(1, length(private$paths))

                         if (!is.null(str) & is.numeric(str) & length(str)==length(paths)) private$str = str
                         else private$str = rep("+", length(private$paths))

                         return(self)
                     }
                 ),
                 active=list(
                     segstats = function(){
                         return(private$segs)
                     },

                     td = function(){
                         return(self$gw2gTrack())
                     },
                     path = function(){
                         return(private$paths)
                     },
                     as.grl = function(){
                         if (!is.null(private$grl)) return(private$grl)
                         return(self$gw2grl())
                     },
                     values = function(){
                         return(self$metaCols())
                     }
                 ))

## Utility functions
setxor = function (A, B)
{
    return(setdiff(union(A, B), intersect(A, B)))
}

## ## Euler characteristics of a graph
## euler.chi = function(G){
##     if (!is(G, "igraph")) stop("Invalid input.")

##     ne = ecount(G)
##     nv = vcount(G)
##     nf =
## }

#############################################################
#' @name munlist
#' @title munlist
#'
#' @description
#' unlists a list of vectors, matrices, data frames into a n x k matrix
#' whose first column specifies the list item index of the entry
#' and second column specifies the sublist item index of the entry
#' and the remaining columns specifies the value(s) of the vector
#' or matrices.
#'
#' force.cbind = T will force concatenation via 'cbind'
#' force.rbind = T will force concatenation via 'rbind'
#'
#' @param x list of vectors, matrices, or data frames
#' @param force.rbind logical flag to force concatenation via rbind (=FALSE), otherwise will guess
#' @param force.cbind logical flag to force concatenation via cbind (=FALSE), otherwise will guess
#' @param force.list logical flag to force concatenation via unlist (=FALSE), otherwise will guess
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski9
#' @export
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
{
    if (!any(c(force.list, force.cbind, force.rbind)))
    {
        if (any(sapply(x, function(y) is.null(dim(y)))))
            force.list = T
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1)
            force.rbind = T
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1))
            force.cbind = T
    }
    else
        force.list = T

    if (force.list)
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                     unlist(x)))
    else if (force.rbind)
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                     do.call('rbind', x)))
    else if (force.cbind)
        return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                       iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                       do.call('cbind', x))))
}


#' ra_breaks: utility function to read junction data from various common formats
#'
#' @name ra_breaks
#' @import VariantAnnotation
#' @export
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
gr2json = function(intervals, file, y = rep("null", length(intervals)), labels = '', maxcn = 100, maxweight = 100)
{

                                        # ++ = RL
                                        # +- = RR
                                        # -+ = LL
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

###############################################################
#' karyoMIP
#'
#' MIP to locally compute walks in an existing JaBbA reconstruction, note: usually many optimal solutions to a given run.
#' Used by jabba.walk.
#'
#' TODO: Make user friendly, still pretty raw
#'
#' takes |E| x k matrix of k extreme paths (i.e. contigs) across e edges of the karyograph
#' and length |E| vector of edge copy numbers (eclass), length |E| vector of edge equivalence classes (both outputs of jbaMIP.process)
#' and computes most likely karyotypes that fit the edge copy number profile subject to some prior likelihood
#' over the k extreme paths
#'
#' @param K  |E| x k binary matrix of k "extreme" contigs across |E| edges
#' @param e  edge copy numbers across |E| edges
#' @param eclass  edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together, each class can have at most 2 members
#' @param prior  prior log likelihood of a given contig being in the karyotype
#' @param cpenalty karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype, should be calibrated to prior, i.e. higher than the contig-contig variance in the prior, otherwise complex karyotypes may be favored
#' @param tilim time limit to optimizatoin
#' @param nsolutions how many equivalent solutions to report
#' @return
#' Rcplex solution list object with additional field $kcn for path copy number, $kclass for k class id, $mval for mval
#' @export
###############################################################
karyoMIP = function(K, # |E| x k binary matrix of k "extreme" contigs across |E| edges
                    e, # edge copy numbers across |E| edges
                    eclass = 1:length(e), # edge equivalence classes, used to constrain strand flipped contigs to appear in solutions together,
                                        # each class can have at most 2 members
                    kclass = NULL,
                    prior = rep(0, ncol(K)), # prior log likelihood of a given contig being in the karyotype
                    cpenalty = 1, # karyotype complexity penalty - log likelihood penalty given to having a novel contig in the karyotype,
                                        # should be calibrated to prior, i.e. higher than the contig-contig variance in the prior,
                                        # otherwise complex karyotypes may be favored
                    tilim = 100, epgap = 1, nsolutions = 50, objsense = 'max', ...)
{
    require(Rcplex)

    M = 1e7;
    K = as(K, 'sparseMatrix')

    if (length(prior)!=ncol(K))
        stop('prior must be of the same length as number of columns in K')

                                        # variable indices
    v.ix = 1:ncol(K)
    M.ix = max(v.ix) + (1:ncol(K))
    n = max(M.ix);

                                        # add big M constraints
    Zero = sparseMatrix(1, 1, x = 0, dims = c(n, n)) # upper bound is infinity if indicator is positive
    Amub = Zero[1:length(M.ix), ]
    Amub[cbind(1:length(M.ix), v.ix)] = 1
    Amub[cbind(1:length(M.ix), M.ix)] = -M

    Amlb = Zero[1:length(M.ix), ] # lower bound > 0 if indicator is positive
    Amlb[cbind(1:length(M.ix), v.ix)] = 1
    Amlb[cbind(1:length(M.ix), M.ix)] = -0.1

    if (is.null(kclass))
        kclass = .e2class(K, eclass)

    kclass.counts = table(kclass)
    if (any(kclass.counts>1)) ## any equiv i.e. strand flipped contig pairs? then make sure they appear in solutions togethrer
    {
        bikclass = which(kclass.counts>1)
        Ac = Zero[1:length(bikclass), , drop=F]
        pairs = matrix(unlist(split(1:length(kclass), kclass)[as.character(bikclass)]), ncol = 2, byrow = T)
        Ac[cbind(1:nrow(pairs), pairs[,1])] = 1
        Ac[cbind(1:nrow(pairs), pairs[,2])] = -1
    }
    else
        Ac = Zero[1,,drop = FALSE]

                                        # combine constraints
    A = rBind(cBind(K, Zero[rep(1, nrow(K)), M.ix]), Amub, Amlb, Ac);
    b = c(e, rep(0, nrow(Amlb)*2), rep(0, nrow(Ac)));
    sense = c(rep('E', nrow(K)), rep('L', nrow(Amlb)), rep('G', nrow(Amlb)), rep('E', nrow(Ac)))
    vtype = c(rep('I', length(v.ix)), rep('B', length(M.ix)))

    cvec = c(rep(0, length(v.ix)), prior-cpenalty*rep(1, length(M.ix)))


    sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))

    if (!is.null(sol$xopt))
        sol = list(sol)

    sol = lapply(sol, function(x)
    {
        x$kcn = round(x$xopt[v.ix])
        x$kclass = kclass
        x$mval= round(x$xopt[M.ix])
        return(x)
    })

    return(sol)
}

##############################################################
#' karyoMIP.to.path
#'
#' for a karyoMIP solution and associated K matrix of n x e elementary paths  (input to karyoMIP), and v x e edge signed incidence matrix
#'
#'
#' @param sol solution to karyoMIP
#' @param K matrix of elementary paths (input to karyoMIP)
#' @param e nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
#' @param gr optional GRanges whose names are indexed by rownames of B
#' @param mc.cores integer number of cores
#' @param verbose flag
#' @return
#' A list with following items:
#' $path length k list of paths, cycles (each item i is vector denoting sequence of vertices in G )
#' $is.cycle length k logical vector whose component i denotes whether path i is cyclic
#' $cn  length k integer vector whose component i denotes copy number of contig i
#' $path.grl if path.grl == T
#' @export
##############################################################
karyoMIP.to.path = function(sol, ## karyoMIP solutions, i.e. list with $kcn, $kclass (edges vectors)
                            K, ## K matrix input to karyomip (edges x paths)
                            e, ## nrow(K) x 2 edge matrix representing vertex pairs (i.e. edges to which K is referring to)
                            gr = NULL, ## optional GRanges who names are indexed by <<rownames>> of B
                            mc.cores = 1,
                            verbose = T
                            )
{
    contigs = which(sol$kcn!=0)
    c1 =  contigs[!duplicated(sol$kclass[contigs])]
    c2 = setdiff(contigs, c1)
    c2 = c2[match(sol$kclass[c2], sol$kclass[c1])]
    contigs = c1
    contigs2 = c2

    nm.gr = names(gr)
    names(gr) = NULL

    if (is.null(nm.gr))
        nm.gr  = 1:length(gr)

    if (any(duplicated(nm.gr)))
        nm.gr = 1:length(gr)

    if (!is.character(e))
        e = matrix(as.character(e), ncol = 2)

    out = list();

    i1 = which(!is.na(e[,1]))
    i2 = which(!is.na(e[,2]))
    B = sparseMatrix(as.numeric(c(e[i1,1], e[i2,2])),  c(i1,i2), x = c(rep(-1, length(i1)), rep(1, length(i1))))
    rownames(B) = 1:nrow(B)

    ## tells us whether the given contig is a cycle .. cycles represent any path lacking net flow in a
    ## non-slack vertex

    is.slack = rowSums(is.na(e))!=0

    out$is.cyc = Matrix::colSums(K[is.slack, contigs, drop = F])==0 & Matrix::colSums((B %*% K[, contigs, drop = F])!=0)==0
    out$cn = sol$kcn[contigs]
    out$kix = contigs;
    out$kix2 = contigs2;

    K = K[, contigs, drop = F]
    out$paths = mclapply(1:length(contigs),
                         function(i)
                         {
                             if (verbose)
                                 cat('contig', i, 'of', length(contigs), '\n')

                             k = K[, i]
                             v.all = setdiff(as.vector(e[k!=0,]), NA)
                             ##      v.all = rownames(B)[which(rowSums(abs(B) %*% k)>0)]  ## vertices associated with edges in path / cycle  k

                             if (length(v.all)==1) ## this is a slack to slack path involving 1 node
                                 return(v.all)

                             ## make subgraph corresponding to edges in this path / cycle
                             ##       B.tmp = B[, which(!is.slack)[k[!is.slack]!=0], drop = F] ##
                             ##       so = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x<0))]
                             ##       si = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x>0))]
                             ##       sG = graph(rbind(so, si))
                             ##       sG = graph(rbind(so, si))

                             tmp.e = e[k!=0, ,drop = F]
                             tmp.e = tmp.e[rowSums(is.na(tmp.e))==0,,drop = F]
                             sG = graph(t(tmp.e))

                             if (out$is.cyc[i])
                             {
                                 p.fwd = names(get.shortest.paths(sG, v.all[1], v.all[pmin(length(v.all), 2)])$vpath[[1]])
                                 p.bwd = names(get.shortest.paths(sG, v.all[pmin(length(v.all), 2)], v.all[1])$vpath[[1]])
                                 return(unique(unlist(c(p.fwd, p.bwd))))
                             }
                             else
                             {
                                 io = as.numeric(B[, !is.slack, drop = F] %*% k[!is.slack])
                                 v.in = rownames(B)[io<0][1]
                                 v.out = rownames(B)[io>0][1]
                                 return(names(get.shortest.paths(sG, v.in, v.out)$vpath[[1]]))
                             }
                         }, mc.cores = mc.cores)

    if (!is.null(gr))
    {
        if (is.null(nm.gr))
            nm.gr = names(B)
        names(gr) = NULL
        out$grl = do.call('GRangesList', lapply(out$paths, function(x) gr[match(x, nm.gr), c()]))  ## match non-slack vertices
        names(out$grl) = paste('Contig ', out$kix, ' (CN = ', out$cn, ')', sep = '')
        values(out$grl)$is.cycle = out$is.cyc
    }

    return(out)
}


####################################################
#' jabba.walk
#'
#' Computes walks around all aberrant edges in JABbA object
#'
#' Takes in JaBbA solution and computes local
#' reconstructions around all aberrant edges (default).  Reconstructions (i.e. Huts) consists
#' of collections of walks, each walk associated with a copy number, and a given
#' region (collection of genomic windows).  The interval sum of walks in a given region, weighted
#' by copy numbers will recapitulate the marginal copy profile (as estimated by JaBbA).
#' The reconstruction is chosen to maximize parsimony.
#'
#' Optional flags allow making huts around specific junctions or specified loci (GRangesList)
#'
#' Walks are reconstructed locally within "clustersize" nodes of each aberrant edge, where
#' clustersize is measured by the number of total edges.  Larger cluster sizes may fail to be
#' computationally tractable, i.e. with a highly rearranged genome in an area of dense interconnectivity.
#'
#' @param sol JaBbA object
#' @param outdir output directory
#' @param junction.ix junction indices around which to build walks (default is all junctions)
#' @param loci  loci around which to build walks (over-rides junction.ix), alternatively can be a list of  "all.paths" objects (i.e. each a list utput of initial all.paths = TRUE run  +/- field $prior for walk to re-eval a given all.paths combo
#' @param clustersize size of the cluster to output around the locus or junction of interest
#' @param trim logical flag whether trim in neighborhood of junction (only applicable if loci = NULL, default = TRUE)
#' @param trim.w integer width to which to trim to
#' @param prune flag whether to prune trivial walks for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
#' @param prune.d1 local distance threshold for walk pruning
#' @param prune.d2 referenc distance threshold for walk pruning
#' @param mc.cores number of cores to use, default 1
#' @param genes character vector of gene symbols with which to annotate walk (eg cancer genes)
#' @param verbose logical flag
#' @return list of walk set around each locus or junction that is inputted to analysis, each list item is a list with the following fields
#' $win = input locus of interest, $grl = GRangesList of walks, $grs is a collapsed footprint of all walks in the walk list for this locu
#' $td gTrack of of the output, additional outputs for debugging: $sol, $K, $Bc, $eix, $vix, $h
#' @export
####################################################
jabba.walk = function(sol, kag = NULL, digested = T, outdir = 'temp.walk', junction.ix = NULL, loci = NULL, clustersize = 100,
                      trim = FALSE, ## whether to trim around junction (only applicable when loci = NULL)
                      trim.w = 1e6, ## how far to trim in neighborhood of junction (only applicable when loci = NULL
                      prune = FALSE, ## whether to prune trivial walks i.e. those for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
                      prune.d1 = 1e5, ## local distance threshold for walk pruning
                      prune.d2 = 1e5, ## reference distance threshold for walk pruning
                      maxiterations = Inf, mc.cores = 1, genes = read.delim('~/DB/COSMIC/cancer_gene_census.tsv', strings = F)$Symbol, verbose = T, max.threads = 4, customparams = T, mem = 6, all.paths = FALSE, nomip = F, tilim = 100, nsolutions = 100, cb.interval = 1e4, cb.chunksize = 1e4, cb.maxchunks = 1e10)
{
    system(paste('mkdir -p', outdir))
    ## awkward workaround to limit the number of processors Cplex will gobble up
    ##

    if (customparams)
    {
        out.file = paste(outdir, 'tmp.prm', sep = '/')
        max.threads = Sys.getenv("LSB_DJOB_NUMPROC")
        if (nchar(max.threads) == 0)
            max.threads = Inf
        else
            max.threads = as.numeric(max.threads)
        max.threads = min(max.threads, mc.cores)
        if (is.infinite(max.threads))
            max.threads = 0

        param.file = paste(out.file, '.prm', sep = '')
        .cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)

        Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
        print(Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
    }


    if (is.null(sol))
        sol = kag

    if (is.null(sol$segstats))
    {
        sol$segstats = sol$tile
        sol$segstats$cn = 2
        sol$segstats$eslack.out = 0
        sol$segstats$eslack.in = 0
    }

    if (is.null(kag))
        kag = sol


    out = list()
    tmp.adj = sol$adj
    if (digested)  ## if input is already "digested", then don't need to bother with slacks
    {
        sol$segstats$eslack.in = 0
        sol$segstats$eslack.out = 0
        G = sol$G
    }
    else ## soon to be deprecated
    {
        ix = which(sol$segstats$eslack.in!=0 | sol$segstats$eslack.out!=0)
        tmp.adj[ix, ix] = 0
        pos.ix = which(strand(sol$segstats)=='+')
        sol$segstats$tile.id = match(gr.stripstrand(sol$segstats), gr.stripstrand(sol$segstats[pos.ix]))
        G = graph.adjacency(tmp.adj!=0)
    }

    h = jbaMIP.process(sol)

    if (verbose)
        cat(paste('Finished processing JaBbA, getting ready to construct walks\n'))

    if (is.null(junction.ix) & is.null(loci))
        junction.ix = 1:nrow(kag$ab.edges)

    if (!is.null(junction.ix))
        if (is.null(names(junction.ix)))
            names(junction.ix) = 1:length(junction.ix)


    if (is.null(loci)) ## junction.ix should be not null here
    {
        loci = do.call('GRangesList', mclapply(junction.ix, function(i)
        {
            if (verbose)
                cat(paste('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nDefining subgraph around junction', i, '\n'))
            vix = vix.i = setdiff(kag$ab.edges[i, 1:2, ], NA)
            if (length(vix)==0)
                return(GRanges())
            k = 0
            last.clustersize = 0
            while (length(vix)<clustersize & k < maxiterations & length(vix)>last.clustersize)
            {
                k = k + 1
                last.clustersize = length(vix)
                vix = unique(unlist(neighborhood(G, vix.i, order = k)))
            }
            if (verbose)
                cat(paste('Outputting', length(vix), 'vertices around junction', i, '\n'))

            return(kag$segstats[vix])
        }
      , mc.cores = mc.cores))

        names(loci) = names(junction.ix)
        loci = loci[sapply(loci, length)>0]
    }
    else ## if loci are provided (i.e. not junction centric) then we will not trim or prune
    {
        trim = F
        prune = F
    }

    if (verbose)
        cat(paste('Finished defining subgraphs\n'))

    starts = gr.start(sol$segstats, ignore.strand = F)
    ends = gr.end(sol$segstats, ignore.strand = F)

    names(sol$segstats) = 1:length(sol$segstats)

    if (is.null(names(loci)))
        lnames =  paste('locus', 1:length(loci), sep = '')
    else
        lnames = names(loci)

    all.junc.pair = c(paste(sol$ab.edges[, 1, 1], sol$ab.edges[, 2, 1], sep = ','), paste(sol$ab.edges[, 1, 2], sol$ab.edges[, 2, 2], sep = ','))
    names(all.junc.pair) = c(1:nrow(sol$ab.edges), -c(1:nrow(sol$ab.edges)))

    if (length(loci)>0)
    {
        out = mclapply(1:length(loci), function(i)
        {

            label = lnames[i]
            outfile.rds = sprintf('%s/%s.rds', outdir, label)
            outfile.pdf = sprintf('%s/%s.pdf', outdir, label)
            outfile.txt = sprintf('%s/%s.txt', outdir, label)
            outfile.allpaths.txt = sprintf('%s/%s.allpaths.txt', outdir, label)
            if (is(loci[[i]], 'GRanges'))
            {
                vix = which(gr.in(kag$segstats, loci[[i]]))
                cat('Number of vertices:', length(vix), '\n')
                eix = which((h$e.ij[,1] %in% vix | h$e.ij[,2] %in% vix) & h$e>0)
                Bc = as.matrix(h$B)[vix, eix]
                K = tryCatch(convex.basis(Bc, interval = cb.interval, chunksize = cb.chunksize, verbose = T, maxchunks = cb.maxchunks), error = function(e) as.character(e))
                if (is.character(K))
                          return(list(README = K))
                      prior = rep(1, ncol(K))
                  }
              else ## assume we are re-heating a previous all.paths = TRUE output (and presumably adding a prior)
                  {
                      K = loci[[i]]$K
                      h = loci[[i]]$h
                      eix = loci[[i]]$eix
                      Bc = loci[[i]]$Bc
                      vix = loci[[i]]$vix
                      prior = rep(1, ncol(K))
                      if (!is.null(loci[[i]]$prior))
                          prior[c(values(loci[[i]]$allpaths.og)$kix,values(loci[[i]]$allpaths.og)$kix2)]  = loci[[i]]$prior
                      loci[[i]] = loci[[i]]$win
                  }

          is.cyc = Matrix::colSums(K[h$etype[eix] == 'slack', ])==0 & Matrix::colSums((h$B[, eix, drop = F] %*% K)!=0)==0
          karyo.sol = karyoMIP(K, h$e[eix], h$eclass[eix], nsolutions = nsolutions, tilim = tilim, cpenalty = 1/prior)
          kag.sol = karyo.sol[[1]]
          p = karyoMIP.to.path(kag.sol, K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores))
          values(p$grl)$cn = p$cn
          values(p$grl)$is.cyc = p$is.cyc
##          td.rg$stack.gap = 5e6

          if (!is.null(kag$junctions))
            {
              values(kag$junctions)$lwd = sol$adj[kag$ab.edges[,1:2, 1]]
              values(kag$junctions)$lty = 1
              values(kag$junctions)$label = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, sol$adj[kag$ab.edges[,1:2, 1]], '')
              values(kag$junctions)$col = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, alpha('red', 0.3), alpha('white', 0))
            }
          win = streduce(sol$segstats[vix], 1e4)

          y1 = max(sol$segstats$cn[gr.in(sol$segstats, win)], na.rm = T)*1.1
          pdf(outfile.pdf, height = 30, width = 24)
          grs = gr.simplify(grl.unlist(p$grl), 'grl.ix', split = T)
          values(grs) = values(p$grl)
          names(grs) = names(p$grl)

          if (!is.null(sol$td))
            {
              td.seg = sol$td
              td.seg$y1 = y1
              td = td.seg
              ## td = c(td.seg, td.rg)
            }
          else
              {
                  td.seg = gTrack(sol$segstats, y.field = 'cn', angle = 0, col ='black', height = 6, labels.suppress = T, y1 = y1)

                                        #          td = c(gTrack(grs, draw.paths = T, path.cex.arrow = 0, border = NA, angle = 0, ywid = 0.5, path.stack.x.gap = 1e6, height = 20, labels.suppress.gr = T),

                  gt.walk = gTrack(grs, draw.paths = T, border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
                  gt.walk$path.cex.arrow = 0
                  gt.walk$path.stack.x.gap = 1e6
                  td = c(
                      gt.walk,
                      td.seg)

                  plot(td,
                       windows = win, links = kag$junctions)
                  dev.off()
              }

          df = data.frame(label = label, cn = p$cn, walk = sapply(grs, function(x) paste(gr.string(x, mb = F), collapse = ',')), widths = sapply(grs, function(x) paste(width(x), collapse = ',')), width = sapply(grs, function(x) sum(width(x))), numpieces = sapply(grs, length), type = 'walk')
          df = rbind(data.frame(label = label, cn = NA, walk = paste(gr.string(win, mb = F), collapse = ','), widths = paste(width(win), collapse = ','), width = sum(width(win)), type = 'window', numpieces = length(win)), df)
          write.tab(df, outfile.txt)
          out = list(
            win = win, grl = p$grl, grls = grs, td = td, sol = karyo.sol,
            K = K, Bc = Bc, eix = eix, vix = vix, h = h,
            README = 'win=windows, grl = raw granges list corresponding to paths, grls = simplified granges list corresponding to paths, td = gTrack object plotting walks, sol = solution object from karyoMIP of local walks, K = incidence matrix input to karyomip, Bc = input to convex.basis, eix = eix input to karyomip, vix = vix input corresponding to rows of Bc, h = h input to karyomip')

          if (all.paths)
            {
              outfile.allpaths.pdf = sprintf('%s/%s.allpaths.pdf', outdir, label)

              if (verbose)
                cat('Generating all walks\n')

              ## repurpose karyoMIP.to.path to generate all paths using "fake solution" i.e. all 1 weights,  to karyoMIP as input
              pallp = karyoMIP.to.path(list(kcn = kag.sol$kcn*0 + 1, kclass = kag.sol$kclass), K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores), verbose = verbose)
              allp = pallp$grl

              allps = gr.simplify(grl.unlist(allp), 'grl.ix', split = T)
              allps[values(allp)$is.cycle] = do.call('GRangesList', lapply(which(values(allp)$is.cycle), function(x) c(allps[[x]], allps[[x]])))
              allps.og = allps; ## save for later
              values(allps.og)$kix = pallp$kix
              values(allps.og)$kix2 = pallp$kix2

              ## text encoding of junctions
              if (!is.null(junction.ix))
                junc.pair = paste(sol$ab.edges[junction.ix[i], 1, ], sol$ab.edges[junction.ix[i], 2, ], sep = ',')

              if (trim | prune) ## junction.ix should be not null here (i.e. they were provided as input or loci = NULL)
                {
                  allps.u = grl.unlist(allps)
                  allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                  allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
                  allps = split(allps.u, allps.u$grl.ix)
                  allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
                  allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
                  allps.w = split(width(allps.u), allps.u$grl.ix)
                  allps.endc = split(levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum), allps.u$grl.ix)

                  if (trim) ## only include windows around the junction of interest
                    {
                      ## allps.ix.pairs tells us what junction indices are present in a walk collection
                      allps.ix.pairs = mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F)
                      ## first, which windows contain the junction

                      wix = which(sapply(allps.ix.pairs, length)>0)
                      allps = allps[wix]

                      if (length(allps)>0)
                        {
                          allps.ixs = allps.ixs[wix] ## start interval id of kth interval in ith walk
                          allps.ixe = allps.ixe[wix] ## end interval id of kth interval in ith walk
                          allps.endc = allps.endc[wix] ## end walk coordinate of kth interval in ith walk
                          allps.w = allps.w[wix]
                          allps.ix.pairs = allps.ix.pairs[wix]

                          ## start window for trimming
                          values(allps)$allps.junc.first =
                            pmax(0, mapply(function(x, y) y[x[1]], allps.ix.pairs, allps.endc)) ## walk position of first junction
                          values(allps)$allps.junc.last =
                            pmax(0, mapply(function(x, y) y[x[length(x)]], allps.ix.pairs, allps.endc)) ## walk position of last junction

                          ## check for any quasi-palindromic walks that contain both orientations of a junction
                          ## split each of these into two so we can maintain the width limit
                          pal.wix = which(values(allps)$allps.win.firstix != values(allps)$allps.win.lastix)
                          if (length(pal.wix)>0)
                            {
                              allps.dup = allps[pal.wix]
                              values(allps.dup)$allps.junc.first = values(allps)$allps.junc.last
                              allps = c(allps, allps.dup)
                              allps.endc = c(allps.endc, allps.endc[pal.wix])
                              allps.w = c(allps.w, allps.w[pal.wix])
                            }

                          values(allps)$allps.win.first =
                            pmax(0, values(allps)$allps.junc.first - trim.w) ## walk coordinate of new window start
                          values(allps)$allps.win.last =
                            pmin(sapply(allps.endc, function(x) x[length(x)]), values(allps)$allps.junc.first + trim.w) ## walk coordinate of new window end
                          values(allps)$allps.win.firstix = ## first walk interval to trim to
                            mapply(function(x, y) setdiff(c(which(x>y)[1], 1), NA)[1], allps.endc, values(allps)$allps.win.first)
                          values(allps)$allps.win.lastix = ## last walk interval to trim to
                            mapply(function(x, y) setdiff(c(which(x>y)[1], length(x)), NA)[1], allps.endc, values(allps)$allps.win.last)
                          values(allps)$allps.win.first.keep =
                            mapply(function(p,e,i) e[i] - p, values(allps)$allps.win.first, allps.endc, values(allps)$allps.win.firstix)
                          values(allps)$allps.win.last.keep =
                            mapply(function(p,e,i,w) w[i] - (e[i] - p), values(allps)$allps.win.last, allps.endc, values(allps)$allps.win.lastix, allps.w)
                          ## apply trimming
                          ## we are trimming walks so that they are within trim.w bases of junction
                          allps.u = grl.unlist(allps)
                          iix = mapply(function(x,y) y %in% values(allps)$allps.win.firstix[x]:values(allps)$allps.win.lastix[x], allps.u$grl.ix, allps.u$grl.iix)
                          allps.u = allps.u[iix]
                          allps.u$keep.end = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.firstix[x], values(allps)$allps.win.first.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)
                          allps.u$keep.start = mapply(function(x, y)
                            ifelse(y == values(allps)$allps.win.lastix[x], values(allps)$allps.win.last.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)

                          if (any(tmp.ix <- !is.na(allps.u$keep.start))) ## we keep the end of the first segment
                            allps.u[tmp.ix] = gr.start(allps.u[tmp.ix], allps.u$keep.start[tmp.ix], ignore.strand = F)

                          if (any(tmp.ix <- !is.na(allps.u$keep.end))) ## we keep the beginning of the last segment
                            allps.u[tmp.ix] = gr.end(allps.u[tmp.ix], allps.u$keep.end[tmp.ix], ignore.strand = F)

                          ## if there are multiple walks with the same aberrant junction set, then pick the longest of these

                          ## first need to find the aberrant walks in each set
                          ij = paste(allps.u$ix.e[-length(allps.u)], allps.u$ix.s[-1], sep = ',') ## indices of all walk adjacent interval pairs
                          names(ij) = 1:length(ij)
                          ij = ij[diff(allps.u$grl.ix)==0] ## only pick intra-walk interval pairs
                          ij.ix = names(all.junc.pair)[match(ij, all.junc.pair)]
                          ## then compute the width of each walk

                          allps = split(allps.u, allps.u$grl.ix)
                          ij.ix.l = split(ij.ix, allps.u$grl.ix[as.numeric(names(ij))])[names(allps)]
                          values(allps)$ab.junc = lapply(ij.ix.l, paste, collapse = ',')
                          values(allps)$wid = vaggregate(width(allps.u), by = list(allps.u$grl.ix), FUN = sum)[names(allps)]
                          ix.w = order(-values(allps)$wid)
                          allps = allps[ix.w[which(!duplicated(values(allps)$ab.junc[ix.w]))]] ## only keep the longest non-duplicate walks
                        }
                    }

                  ## now dedup and trim contigs to locus (mainly useful if loci was provided as argument)
                  if (length(allps)>0)
                    {
                      win = reduce(gr.stripstrand(loci[[i]]))
                      allps.u = grl.unlist(allps)

                      ## trim to locus
                      ix = gr.match(allps.u, win)
                      allps.u = allps.u[!is.na(ix)]
                      ix = ix[!is.na(ix)]
                      start(allps.u) = pmax(start(allps.u), start(win)[ix])
                      end(allps.u) = pmin(end(allps.u), end(win)[ix])

                      allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
                      allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)

                      ## remove dups
                      allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of intervals
                      allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of intervals

                      allps.u = allps.u[allps.u$grl.ix %in% which(!duplicated(paste(sapply(allps.ixs, paste, collapse = ','), sapply(allps.ixe, paste, collapse = ','))))]
                      allps = split(allps.u, allps.u$grl.ix)
                    }


                  if (prune & length(allps)>0)
                    ## this is to prune pseudo-aberrant walks that basically consist of short insertions of non-reference
                    ## sequences in a big reference chunk
                    {
                      ## for each walk create graph of intervals by determining whether pair ij is BOTH near on the walk (<= d1)
                      ## and near on the refernce (<= d2)
                      allps.u = grl.unlist(allps)

                      ## what are the ij pairs we want to test from this collapsed list
                      ij = merge(cbind(i = 1:length(allps.u), ix = allps.u$grl.ix), cbind(j = 1:length(allps.u), ix = allps.u$grl.ix))[, c('i', 'j')]

                      tmp = levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum)
                      allps.u.ir = IRanges(tmp - width(allps.u) + 1, tmp)

                      ## distance on the walk
                      D1 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(IRanges(start = end(allps.u.ir[ij[,'i']]), width = 1),
                                   IRanges(start(allps.u.ir[ij[,'j']]), width = 1))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      ## distance on the reference
                      D2 = sparseMatrix(ij[, 'i'],ij[, 'j'],
                        x = suppressWarnings(
                          distance(gr.end(allps.u[ij[,'i']], ignore.strand = F),
                                   gr.start(allps.u[ij[,'j']], ignore.strand = F))) + 1e-5, dims = rep(length(allps.u.ir), 2))

                      D1 = pmin(as.matrix(D1), as.matrix(t(D1)))
                      D2 = pmin(as.matrix(D2), as.matrix(t(D2)))

                      tmp = D1>0 & D1<prune.d1 & D2>0 & D2<prune.d2
                      tmp[which(is.na(tmp))] = FALSE
                      G = graph.adjacency(tmp)
                      cl = clusters(G, 'weak')$membership ## clusters based on this adjacency relationship
                      cls = split(1:length(cl), cl)
                      lens = sapply(allps, length)

                      ## check if there any clusters that contain both the first and last member  of a walk
                      cls.fl = cls[mapply(function(x) all(c(1,lens[allps.u$grl.ix[x[1]]]) %in% allps.u$grl.iix[x]), cls)]

                      if (length(cls.fl)>0)
                        {
                          toprune = allps.u$grl.ix[sapply(cls.fl, function(x) x[1])]
                          if (length(toprune)>0)
                            cat('Pruning', length(toprune), 'walks\n')
                          allps = allps[-toprune]
                        }
                    }
                }

              if (length(allps)>0)
                win = streduce(unlist(allps), 0)
#                win = streduce(unlist(allps), sum(width(unlist(allps)))*0)

              values(allps) = NULL
              out$allpaths = allps
              out$allpaths.og = allps.og ## untouched all.paths if we want to reheat eg after computing 10X support
              gt.walk = gTrack(out$allpaths, draw.paths = T,border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
              gt.walk$path.cex.arrow = 0
              gt.walk$path.stack.x.gap = 1e6
              out$td.allpaths = c(
                  gt.walk,
                  td.seg)
              pdf(outfile.allpaths.pdf, height = 30, width = 24)
              plot(out$td.allpaths,
                      windows = win, links = kag$junctions)
              dev.off()
              out$README = paste(out$README, 'allpaths= all paths through windows (not just optimal ones), td.allpaths = gTrack object of plot of all paths')
            }

          ## if junction.ix was specified then label which positions in the walks represent the rearrangement junction
          if (!is.null(junction.ix) & length(out$allpaths)>0)
            {
              allps = out$allpaths
              allps.u = grl.unlist(allps)
              allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
              allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
              allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
              allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
              allps.ix.pairs = sapply(mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F), paste, collapse = ',')
              values(allps)$junction.id = names(junction.ix)[i]
              values(allps)$junction.ix = allps.ix.pairs
              out$allpaths = allps
            }

          if (length(out$allpaths)>0)
            {
              values(out$allpaths)$string = grl.string(out$allpaths)
              values(out$allpaths)$wid = sapply(out$allpaths, function(x) sum(width(x)))
              values(out$allpaths)$wids = sapply(out$allpaths, function(x) paste(width(x), collapse = ','))
              write.tab(as.data.frame(values(out$allpaths)), outfile.allpaths.txt)
            }

          saveRDS(out, outfile.rds)
          return(out)
        }, mc.cores = mc.cores)
    }

  ## awkward workaround to limit the number of processors Cplex will gobble up
 if (customparams)
    {
      system(paste('rm', param.file))
      Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')
      cat('Finished\n')
    }

  return(out)
}

## accessory function for walk
## cplex set max threads (warning can only do once globally per machine, so be wary of multiple hosts running on same machine)
##
.cplex_customparams = function(out.file, numthreads = 0, nodefileind = NA, treememlim = NA)
{
  param_lines = "CPLEX Parameter File Version 12.6.0.0"

  param_lines = c(param_lines, paste("CPX_PARAM_THREADS", numthreads, sep = '\t'))

  if (!is.na(nodefileind))
    param_lines = c(param_lines, paste("CPX_PARAM_NODEFILEIND", nodefileind, sep = '\t'))

  if (!is.na(treememlim))
    {
#'      param_lines = c(param_lines, paste("CPX_PARAM_WORKDIR", getwd(), sep = '\t'))
      param_lines = c(param_lines, paste("CPX_PARAM_TRELIM", treememlim, sep = '\t'))
    }

  writeLines(param_lines, out.file)
  Sys.setenv(ILOG_CPLEX_PARAMETER_FILE=out.file)
}

sparse_subset = function (A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
{
    nz = Matrix::colSums(A != 0, 1) > 0
    if (is.null(dim(A)) | is.null(dim(B)))
        return(NULL)
    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))
    for (i in seq(1, nrow(A), chunksize)) {
        ixA = i:min(nrow(A), i + chunksize - 1)
        for (j in seq(1, nrow(B), chunksize)) {
            ixB = j:min(nrow(B), j + chunksize - 1)
            if (length(ixA) > 0 & length(ixB) > 0 & !quiet)
                cat(sprintf("\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n",
                  ixA[1], ixA[length(ixA)], nrow(A), ixB[1],
                  ixB[length(ixB)], nrow(B)))
            if (strict)
                C[ixA, ixB] = (sign((A[ixA, , drop = FALSE] !=
                  0)) %*% sign(t(B[ixB, , drop = FALSE] != 0))) *
                  (sign((A[ixA, , drop = FALSE] == 0)) %*% sign(t(B[ixB,
                    , drop = FALSE] != 0)) > 0)
            else C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE] !=
                0) %*% sign(t(B[ixB, nz, drop = FALSE] == 0))) ==
                0
        }
    }
    return(C)
}

## TODO:
