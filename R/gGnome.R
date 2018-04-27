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
#' @import methods
#' @import R6
#' @import data.table
#' @import Matrix
#' @import jsonlite
#' @import GenomicRanges
#' @import igraph
#' @import gUtils
#' @import gTrack
"_PACKAGE"

## rules
## the only required field for nodes is the GRanges itself.
## the nodes must be +/- paired up, eg if a (1[x, y]+) exists, (1[x, y]-) must exist
## there can be duplicated ranges
## cn is optional, loose is optional
## the only required fields for edges is the "from" and "to"
## type is optional, cn is optional
## you can store cn==0 junctions, just their corresponding edges must be cn==0 too
## gg$junctions will include them

## ================ symsegs ============== ##
#' @name symsegs
#' @title nodes of a skew-symmetric genomic graph
#' @docType class
#'
#' @description
#' This is a special class of \code{GRanges} in which there defined a mapping function \code{hb}
#' that maps the pairs of two ssDNA that constitutes the dsDNA molecules such that at any time
#' ssDNA \code{s1} is present in an \code{symsegs} object, \code{s2=hb(s1)} must also be present.
#'
#' The mapping funciton \code{hb} can be think of as getting the reverse complement of a ssDNA.
#' Nevertheless, the properties that sufficiently and necessarily define a unique ssDNA is not
#' merely its genomic ranges, but flexible and extendible by defining the property fields of
#' it. For example, we can have a type-\code{factor} field called "alleles" to represent each
#' different alleles even of the same genomic ranges. Obviously one can define an arbitrary number
#' of property fields on a \code{symsegs} object and the set of the elements is the Cartisan
#' product of all of them. To extend your imagination, we can also have "loose" for loose ends,
#' or "sample.id" for storing multiple samples in the same object. A property field must be
#' boolean, categorical, or nominal.
#'
#' @import methods
#' @import GenomicRanges
#' @importClassesFrom GenomicRanges GRanges
#'
#' @export symsegs
#' @exportClass symsegs
symsegs = setClass("symsegs",
                   contains = "GRanges",
                   representation = list(
                       prop = "character"
                   ))

setValidity("symsegs",
            function(object){
                verbose = getOption("gGnome.verbose")
                if (length(object)==0){
                    if (verbose){message("Empty ranges.")}
                    return(TRUE)
                } else if (any(as.logical(strand(object)=="*"))){
                    if (verbose){message("Every range must be strand specific.")}
                    return(FALSE)
                } else {
                    hb = hbonds(object)
                    if (any(hb[, is.na(from)| is.na(to)])){
                        return(FALSE)
                    } else {
                        return(TRUE)
                    }
                }
            })

## ================ junctions ============== ##
############################################
#' @name junctions-class
#' @title junctions: GRangesList to store junction info
#' @docType class
#'
#' @description
#' S4 class representing the geomic structural variations based on a reference genome. A
#' stuctural variation or junction is a pair of breakpoints represented by strand-specific
#' width 1 ranges.
#'
#' @import methods
#' @import GenomicRanges
#' @importClassesFrom GenomicRanges GRangesList
#'
#' @export junctions
#' @export
junctions = setClass("junctions",
                     contains="GRangesList")

## validity test when intializing
setValidity("junctions",
            function(object){
                verbose = getOption("gGnome.verbose")
                if (length(object)==0){
                    if (verbose){message("Empty junction set.")}
                    return(TRUE)
                }
                else if (!all(IRanges::elementNROWS(object)==2)){
                    if (verbose){message("Each element must be length 2.")}
                    return(FALSE)
                }
                else if (is.element("*", as.character(strand(unlist(object))))){
                    if (verbose){message("All strand info must be present.")}
                    return(FALSE)
                }
                else{
                    return(TRUE)
                }
            })
## ra.duplicated
## ra.bedpe
## ra.dedup
## ra.match

## ----------- functions for junctions --------- ##
#' @name ra.duplicated
#' @title marking duplicated junctions as \code{TRUE}
#' @description
#' Similar to \code{duplicated}, will return \code{logical} of the same length as the input,
#' where the \code{TRUE} denotes junctions that have been found in previous part of the input.
#'
#' @importFrom gUtils ra.overlaps
#'
#' @export
ra.duplicated = function(junc){
    if (!inherits(junc, "junctions")){
        stop("Only works for a junctions object.")
    }

    if (length(junc)==0){
        return(logical(0))
    } else {
        ov = suppressWarnings(gUtils::ra.overlaps(junc, junc))
        ov = data.table::as.data.table(ov)[ra1.ix != ra2.ix]
        if (nrow(ov)==0){
            return(rep(FALSE, length(junc)))
        } else {
            return(seq_along(junc) %in% ov[, unique(pmax(ra1.ix, ra2.ix))])
        }
    }
}

## ================== gGraph class definition =========== ##
#' @export
gGraph = setClass("gGraph")

############################################################
#' @name gGraph-class
#' @title genomic rearrangement graph
#' @docType class
#' @description
#' The main work horse of this package. Rearrangement graph G=(V, E), where V is a set of
#' strand-specific \code{GRanges} that both strand of any range must be present, and E is
#' a set of directed edges connecting adjacent nodes stored in the form of \code{data.table}
#' with two required columns \code{from} and \code{to} that matches the node's index in V.
#'
#' Every gGraph must be defined on a reference genome, and that is stored in the \code{seqinfo}
#' of V. Optional metadata is allowed and appended as extra columns in V or E, some of which
#' are required by the descendant classes like \code{bGraph}.
#'
#' Nodes and edges are necessary and sufficient to define a \code{gGraph} instance, while
#' optional metadata fields like copy numbers, edge attributes can be extended.
#'
#' In the following examples \code{gg} is a gGraph object.
#'
#' @usage
#' Constructor:
#'   gGraph$new(tile=NULL, juncs=NULL, cn = FALSE,
#'              jabba=NULL,
#'              weaver=NULL,
#'              prego=NULL,
#'              segs=NULL, es=NULL,
#'              ploidy=NULL, purity=NULL)
#'
#'   gread(filename)
#'
#' Public fields:
#'   gg$segstats
#'
#'   gg$edges
#'
#'   gg$junctions
#'
#'   gg$G
#'
#'   gg$adj
#'
#'   gg$A
#'
#'   gg$parts
#'
#' Public methods:
#'
#'   gg$seqinfo()
#'
#'   seqinfo(gg)
#'
#'   gg$nullGGraph()
#'
#'   gg$simpleGraph(genome = NULL, chr = FALSE, include.junk = FALSE)
#'
#'   gg$dipGraph(genome = NULL, chr = FALSE, include.junk = FALSE)
#'
#'   gg$karyograph(tile = NULL, juncs = NULL)
#'
#'   gg$addJuncs(juncs)
#'
#'   gg$addSegs(bps)
#'
#'   gg$jab2gg(jabba, regular=NULL)
#'
#'   gg$wv2gg(weaver)
#'
#'   gg$pr2gg(prego)
#'
#'   gg$print()
#'
#'   print(gg)
#'
#'   gg$plot(pad = 1000) ## TODO: add node and edge viz configurations
#'
#'   plot(gg)
#'
#'   gg$layout()
#'   ## TODO: rewrite layout(gg); add layout method option
#'
#'   gg$summary()
#'   ## TODO: rewrite summary(gg)
#'
#'   gg$simplify(mod=TRUE)
#'
#'   gg$decouple(mod=TRUE)
#'
#'   ## TODO: gg$add(); gg$subtract()
#'
#'   gg$length()
#'
#'   length(gg)
#'
#'   gg$gg2td()
#'   ## TODO: add node and edge appearences
#'
#'   gg$json(filename = ".", maxcn = 100, maxweight = 100)
#'
#'   gg$gg2js(filename = ".", maxcn = 100, maxweight = 100, save = TRUE)
#'
#'   gg$html(filename = ".", gGnome.js = Sys.getenv("DEFAULT_GGENOMEJS"), maxcn = 100, maxweight=100)
#'
#'   gg$components(mc.cores = 1)
#'
#'   components(gg)
#'
#'   gg$subgraph(v = numeric(0), na.rm = T, mod = T)
#'
#'   gg$fillin()
#'
#'   gg$fillup()
#'
#'   gg$trim(gr = NULL)
#'
#'   gg$get.g(force = FALSE)
#'
#'   gg$get.adj()
#'
#'   gg$hood()
#'
#'   gg$dist()
#'
#'   gg$isBalance()
#'
#'   gg$get.loose()
#'
#'   gg$get.walk(v=numeric(0), e=numeric(0), peel=FALSE, cn=NULL)
#'
#'   random.walk(start=numeric(0), steps, mode=c("out","in","all"))
#'
#'   chromoplexy()
#'
#'   chromothripsis()
#'
#'   kid.frag()
#'
#'   bfb()
#'
#' @param tile the \code{GRanges} genome segmentation
#' @param juncs the \code{GRangesList} of SV junctions
#' @param cn \code{logical} if \code{TRUE} honor the copy number annotation in the input
#' @param jabba the path to or actual \code{list} of
#' \href{http://github.com/mskilab/JaBbA}{JaBbA} output
#' @param weaver the directory containing the output files from Weaver
#' \href{https://github.com/ma-compbio/Weaver}{Weaver}
#' @param prego the interval.results output file from PREGO
#' \href{http://compbio.cs.brown.edu/projects/prego/}{PREGO}
#' @param segs \code{GRanges} object of the nodes
#' @param es \code{data.table} object of the edges
#' @param ploidy defined as the width weighted copy number of the nodes
#' @param purity the proportion of cells that has rearranged genome described by the graph
#' in the biological sample, the rest is assumed diploid reference
#' @param filename a path to the input or output file
#' @param genome \code{seqinfo}, \code{seqlengths}, or \code{BSgenome} objects of the reference genome
#' @param chr \code{logical} scalar, whether chromosome names should have the "chr" prefixes
#' @param include.junk \code{logical} scalar, whether to keep the unassembled gaps in reference genome
#' @param ploidy \code{numeric} scalar, the ploidy with which to initialize the simple graph
#' @param bps \code{GRanges} of genomic breakpoints or segments
#' @param regular \code{logical} or \code{character}, defining the \code{seqlevels} of the regular chromosomes,
#' or non-unassembled gaps in a reference genome, if \code{TRUE} will read this info from DEFAULT_REGULAR_CHR
#' @param pad the amount of extension up and downstream of the ranges when defining genomic ranges
#' @param mod \code{logical} scalar, whether to modify the instance by reference (i.e. \code{self}) or produce
#' a copy of the output
#' @param maxcn \code{numeric}, any node copy number exceeding will be capped
#' @param maxweight \code{numeric}, similar to \code{\link{maxcn}}, but for edges
#' @param gGnome.js \code{character}, the path to the repository of gGnome.js
#' @param invoke \code{logical} scalar, whether to start gGnome.js server right away
#'
#'
#'
#'
#' @import R6
#' @import data.table
#' @import igraph
#' @import gUtils
#' @import gTrack
#'
#' @section Details:
#' \subsection{Consructors}{
#' To parse the output from genome graph callers like JaBbA and load it as a gGraph/bGraph object,
#' use \code{gread()}.
#' }
#'
#'
#' @exportClass gGraph
#' @export gGraph
#' @export
############################################################
gGraph = R6::R6Class("gGraph",
                 public = list(
                     ## public fields
                     ## name = NULL,
                     ## refG = "GENOME", ## seqinfo of ref genome

                     ## constructor
                     initialize = function(tile=NULL, juncs=NULL, cn = FALSE,
                                           jabba=NULL,
                                           weaver=NULL,
                                           prego=NULL,
                                           segs=NULL, es=NULL,
                                           ploidy=NULL, purity=NULL,
                                           regular=TRUE){
                         ## control how to construct
                         verbose = getOption("gGnome.verbose")
                         if (!is.null(segs) & !is.null(es)){
                             private$gGraphFromScratch(segs, es,
                                                       junctions,
                                                       ploidy, purity)
                         } else if (!is.null(tile) | !is.null(juncs)) {
                             if (verbose){
                                 message("Initializing with 'tile' and 'junctions'")
                             }
                             self$karyograph(tile, juncs, cn = cn)
                         } else if (!is.null(jabba)) {
                             if (verbose) {
                                 message("Reading JaBbA output")
                             }
                             self$jab2gg(jabba)
                         } else if (!is.null(weaver)) {
                             if (verbose) {
                                 message("Reading Weaver output")
                             }
                             self$wv2gg(weaver)
                         } else if (!is.null(prego)) {
                             if (verbose){
                                 message("Reading Prego output")
                             }
                             self$pr2gg(prego)
                         } else {
                             self$nullGGraph(regular)
                         }
                     },

                     set.seqinfo = function(genome=NULL, gname=NULL, drop=FALSE){
                         "Set the seqinfo of the nodes."
                         if (inherits(genome, "gGraph")){
                             genome = seqinfo(genome)
                         }
                         private$segs = gUtils::gr.fix(private$segs,
                                                       genome=genome,
                                                       gname=gname,
                                                       drop=drop)
                         return(self)
                     },

                     ## initialize from global ref genome seqinfo
                     ## This is actually a diploid graph
                     ## null graph should be all zero
                     nullGGraph = function(regular=TRUE, genome=NULL){
                         "Create empty gGraph."
                         if (!is.null(private$segs)){
                             old.si = seqinfo(private$segs)
                             if (length(old.si@seqlengths)>0){
                                 if (getOption("gGnome.verbose")){
                                     message("Adopting reference instance's seqinfo. Ignoring input.")
                                 }
                                 genome = old.si
                             }
                         }
                         private$segs = GRanges(seqinfo = genome)
                         private$g = igraph::make_empty_graph()
                         private$junction = new("junctions", GRangesList())
                         private$es = data.table(from=integer(0),
                                                 to=integer(0),
                                                 type=character(0))
                         return(self)
                     },

                     simpleGraph = function(genome = NULL, chr=FALSE, include.junk=FALSE, ploidy = NULL){
                         if (!is.null(genome)){
                             tmp = tryCatch(si2gr(genome), error=function(e) NULL)
                             if (is.null(tmp)){
                                 warning("Input 'genome' cannot be converted to a seqinfo. Return default.")
                                 genome = NULL
                             }
                         }

                         if (is.null(genome)){
                             sl = gUtils::hg_seqlengths(genome=genome,
                                                        chr=chr,
                                                        include.junk=include.junk)
                             tmp = si2gr(sl)
                             names(tmp) = NULL
                         }

                         private$segs = c(tmp, gUtils::gr.flipstrand(tmp)) ## null segs are ref
                         if (!is.null(ploidy)){
                             private$segs$cn = ploidy
                         }
                         private$segs$loose = FALSE ## all non-loose end

                         private$es = data.table(from = integer(0),
                                                 to = integer(0),
                                                 type = character(0))

                         sinfo = seqinfo(private$segs)
                         if (any(sinfo@is_circular, na.rm=T)){
                             ## close the circle on mitochondrial DNA
                             circChr = sinfo@seqnames[which(sinfo@is_circular==T)]
                             if (getOption("gGnome.verbose")){
                                 cat(paste('There is', length(circChr), 'circular contig(s): '))
                                 cat(circChr, '\n')
                             }

                             ## QUESTION: why Rle doesn't allow match function???
                             circIx = which(as.vector(seqnames(private$segs)) %in% circChr)
                             if ( length(circIx)>0 ){
                                 ## constructing edges: 5 required columns
                                 ## from, to, cn, type, weight (len o' source node)
                                 private$es = rbindlist(
                                     list(private$es, list(
                                                          from=circIx,
                                                          to=circIx,
                                                          type=rep("reference", length(circIx)))
                                          )
                                 )
                             }
                         }

                         private$reset()
                         private$.ploidy = ploidy
                         return(self)
                     },

                     ## start as diploid chromosomes
                     ## NOTE: depends it on the input's
                     dipGraph = function(genome = NULL,
                                         chr=FALSE,
                                         include.junk=FALSE){
                         "Create diploid reference genome where each chr is a node pair without edges."
                         ## default seqlengths
                         self$simpleGraph(genome = NULL,
                                          chr = FALSE,
                                          include.junk = FALSE,
                                          ploidy = 2)
                         return(self)
                     },

                     ## initialize from segmenatation AND/OR rearrangement junctions
                     addJuncs = function(juncs, cn=TRUE){
                         ## DONE: populate abEdges while adding new junction!!!!
                         ## NOTE: the bps in junc must be width 2

                         ## ALERT: convention of junction orientation!!!
                         if (getOption("gGnome.verbose")){
                             message("Given a GRL of junctions add them plainly to this gGraph.")
                         }

                         if (is.null(juncs)){
                             stop("There has to be some input.")
                         }

                         ## 1. every single junction has 2 breakpoints,
                         ## break nodes in graph by these breakpoints
                         ## 2. based on oreintation of the junctions,
                         ## connect those nodes; introduce corresonding edges to graph
                         ## DONE: check if every bp within the ref genome
                         ## if not we need to resolve, maybe by creating new seqnames with warning
                         if (!inherits(juncs, "junctions")){
                             ## NOTE: for a GRL to be junctions class,
                             ## must be 1) each element length 2 and with strand
                             ## 2) width 2, if not, convert
                             juncs = tryCatch(junctions(juncs),
                                             error = function(e){
                                                 NULL
                                             })
                             if (is.null(juncs)){
                                 stop("Input is not a valid junctions set.")
                             }
                         }

                         if (length(juncs)==0){
                             return(self)
                         }
                         ##################################################
                         ## start processing
                         ## DONE: write as JaBbA::karyograph() with modifications
                         ## e.g. (30, 2) --> pivot (2, 30)
                         ## jadd = jadd[j.in==TRUE & cn>0, jix]
                         if (!"cn" %in% colnames(values(juncs))){
                             values(juncs)$cn = 1
                         }

                         ## it has to be cn>0 if required
                         if (cn){
                             j.non.empty = values(juncs)$cn>0
                         } else {
                             j.non.empty = rep(TRUE, length(juncs))
                         }

                         ## both breakpoints must be in scope
                         j.in.scope = grl.in(juncs, private$segs, only=TRUE)

                         jadd = which(j.in.scope & j.non.empty)

                         ## if nothing to add
                         if (length(jadd)==0){
                             return(self)
                         }

                         ## resize to width 1, left
                         jUl = grl.unlist(juncs)
                         if (!all(width(jUl))==1){
                             jUl = gr.start(jUl)
                         }
                         names(jUl) = NULL

                         bp.p = split(jUl %Q% (grl.ix %in% jadd),
                                      rep(1:2, length(jadd)))
                         bp.p = gr.fix(bp.p, private$segs)
                         juncTile = c(bp.p[[1]], bp.p[[2]])

                         ## BP 1 and 2, retaining strand-orientation info
                         self$addSegs(juncTile, cn=FALSE) ## DONE: addSegs is giving dup edges!!!

                         ## sanity check: at this point, all bps should map to only right bound
                         ## ALERT: THERE MAY BE JUNCTIONS OUT OF SCOPE
                         if (!all(unlist(bp.p) %^% gr.end(private$segs))){
                             ## browser()
                             stop("Error: Something went wrong when breaking up the segs!")
                         }

                         ## ##################################################
                         ## TODO: design change:
                         ## don't assume one breakpoint is only on one node!
                         ## there may be any number of alleles per segment
                         hb = hydrogenBonds(private$segs)
                         hb = hb[, setNames(from, to)]
                         ## now convert bp1 and bp2 to data.table
                         ## bp1 --> every bp associated fwith 4 nodes:
                         ## left +, left -, right +, right -
                         end1 = bp.p[[1]] %+% as.numeric(as.logical(strand(bp.p[[1]])=="+"))
                         anc1 = gr.match(end1, private$segs, ignore.strand=FALSE)
                         to1 = anc1
                         from1 = hb[as.character(anc1)]

                         end2 = bp.p[[2]] %+% as.numeric(as.logical(strand(bp.p[[2]])=="+"))
                         anc2 = gr.match(end2, private$segs, ignore.strand=FALSE)
                         to2 = anc2
                         from2 = hb[as.character(anc2)]
                         ## for the time being we will leave it as is
                         ## Friday, Feb 16, 2018 06:56:30 PM
                         ## #################################################

                         if ("cn" %in% colnames(values(bp.p[[1]]))){
                             abEs = rbind(data.table(from = from1, to = to2, edge.ix = seq_along(anc1), cn = bp.p[[1]]$cn),
                                          data.table(from = from2, to = to1, edge.ix = seq_along(anc2), cn = bp.p[[2]]$cn))
                         } else {
                             abEs = rbind(data.table(from = from1, to = to2, edge.ix = seq_along(anc1)),
                                          data.table(from = from2, to = to1, edge.ix = seq_along(anc2)))
                         }

                         abEs = abEs[,.(from, to, type="unknown", cn)]
                         abEs[, eid := paste(from, to)]
                         abEs[, reid := paste(hb[as.character(to)], hb[as.character(from)])]

                         ## collapse CN!!!
                         abEs[, ":="(cn = sum(cn)), keyby=eid]
                         abEs = abEs[!duplicated(eid)]

                         ematch = abEs[, match(eid, reid)]
                         ## if (any(ub <- abEs[, cn] != abEs[ematch, cn])){
                         ##     message("BOOM")
                         ##     browser()
                         ## }

                         ## make edges
                         if (!"type" %in% colnames(private$es)){
                             private$es$type = "unknown"
                         }

                         if (cn==TRUE){
                             ## haven't decided yet, but a few rules to keep
                             ## 1) if abEs has "cn", we must honor it
                             ## 2) if cn==TRUE, private$es must have "cn" after this
                             ## 3) if segs has "cn", it cannot be altered
                             ## 4) no duplicated edges, if there is, add up the cn
                             if ("cn" %in% colnames(private$es)){
                                 new.es = rbind(abEs[, .(from, to, type, cn)],
                                                private$es[, .(from, to, type, cn)])
                             } else {
                                 new.es = rbind(abEs[, .(from, to, type, cn)],
                                                private$es[, .(from, to, type, cn=0)])
                             }

                             ## sum up the edge cn
                             new.es[, eid := paste(from, to)]
                             ## DONE: this dedups rows
                             new.es[, .(from, to, type, cn=sum(cn)), keyby=eid]
                             new.es[, reid := paste(hb[as.character(to)], hb[as.character(from)])]
                             ematch = new.es[, match(eid, reid)]
                             ## if (any(ub <- new.es[, cn] != new.es[ematch, cn])){
                             ##     message("BAM")
                             ##     browser()
                             ## }
                             new.es = new.es[!duplicated(eid)]
                         } else {
                             new.es = rbind(private$es[, .(from, to, type)],
                                                abEs[, .(from, to, type)])
                             new.es = new.es[!duplicated(paste(from, to))]
                         }

                         et = etype(private$segs, new.es, force=T, both=T)
                         private$segs = et$segs
                         private$es = et$es
                         return(self)
                     },

                     addSegs = function(tile, cn=TRUE){
                         ## How do we consider "cn" field?
                         ## if current segs has no "cn", tile has no "cn", then the output has
                         ## no "cn"; if current doesn't but tile does, we'll assign the "cn"
                         ## of tile to the output; if both have "cn", we'll assign the sum of
                         ## the two to the output.
                         ## Given a GRanges obj of a segmentation (complete or not),
                         ## break the gGraph at their ends.
                         ## extract breakpoints
                         ## bps = reduce(c(gr.start(tile), gr.end(tile)))
                         if (is.null(tile)){
                             stop("There has to be some input.")
                         }

                         if (!inherits(tile, "GRanges")){
                             tile = tryCatch(GRanges(tile),
                                            error=function(e){
                                                NULL
                                            })
                             if (is.null(tile)){
                                 stop("Input cannot be converted into a GRanges object.")
                             }
                         }

                         ## break it
                         private$makeSegs(disjoin(tile))

                         if (cn){
                             if ("cn" %in% colnames(values(tile))){
                                 if ("loose" %in% colnames(values(tile))){

                                     cn.tile = as(coverage(tile %Q% (strand=="+" & loose==FALSE),
                                                           weight="cn"),
                                                  "GRanges")
                                 } else {
                                     cn.tile = as(coverage(tile %Q% (strand=="+"),
                                                           weight="cn"),
                                                  "GRanges")
                                 }

                                 ## private$segs$cn = gr.val(private$segs, tile[, "cn"])$value
                                 ## MOMENT
                                 ## how to make sure that CN is passed????
                                 private$segs = private$segs %$% cn.tile
                                 private$segs$cn = private$segs$score
                                 private$segs$score = NULL
                             }
                         }

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
                         if (is.null(private$es)){
                             ## NOTE: don't understand why es is NULL sometimes
                             private$es = data.table(from = integer(0),
                                                     to = integer(0),
                                                     type = character(0),
                                                     cn = numeric(0))
                         }

                         if (!"cn" %in% colnames(private$es)){
                             private$es[, cn := as.numeric(NA)]
                         }

                         newEs = copy(private$es)
                         newEs[, ":="(from = tmpDt[, which((qid %in% from) & isTail==T)],
                                      to = tmpDt[, which((qid %in% to) & isHead==T)])]

                         ## introduce ref edges between new breakpoints
                         refEs = tmpDt[, .(from=.I[isTail==F],
                                           to=.I[isHead==F],
                                           type = "reference"),
                                       by=qid]
                         refEs = refEs[!is.na(from) & !is.na(to), -c("qid"), with=F]
                         if ("cn" %in% colnames(tmpDt)){
                             cns = tmpDt[, cn]
                             refEs[, ":="(from.cn = cns[from], to.cn = cns[to])]
                             ## if (refEs[, any(from.cn != to.cn, na.rm=T)]){
                             ##     browser()
                             ## }
                             refEs[, cn := pmin(from.cn, to.cn)]
                             if (verbose <- getOption("gGnome.verbose")){
                                 message("We don't keep any edge incident to NA copy nodes.")
                             }
                             refEs = refEs[!is.na(cn)]
                         }

                         if (nrow(newEs)>0){
                             if (all(c("type", "cn") %in% colnames(newEs)) &
                                 all(c("type", "cn") %in% colnames(refEs))){
                                 ## combine the two parts
                                 newEs = rbindlist(list(newEs[, .(from, to, type, cn)],
                                                        refEs[, .(from, to, type, cn)]))
                             } else if ("cn" %in% colnames(refEs)){
                                 ## combine the two parts
                                 newEs = rbind(newEs[, .(from, to,
                                                         type="unknown",
                                                         cn = 0)],
                                               refEs[, .(from, to, type, cn)])
                             } else {
                                 newEs = rbind(newEs[, .(from, to,
                                                         type="unknown",
                                                         cn = 0)],
                                               refEs[, .(from, to, type, cn=0)])
                             }
                         } else {
                             newEs = refEs
                         }
                         newEs = etype(private$segs, newEs, force=TRUE)

                         ## update: es, g
                         private$es = newEs
                         ## reset
                         private$reset()

                         return(self)
                     },

                     ## karograph: initialize `dipGraph()`,
                     ## add junctions to it, then add tiles to it
                     karyograph = function(tile=NULL,
                                           juncs=NULL,
                                           cn=FALSE,
                                           regular=FALSE,
                                           genome = NULL){
                         if (is.null(tile) & is.null(juncs)){
                             return(self)
                         }

                         ## TODO: make this compatible with JaBbA!!?
                         self$simpleGraph(genome = genome)
                         ## no tile, no cn
                         if (is.null(tile)){
                             cn = FALSE
                         }
                         else if (!"cn" %in% colnames(values(tile))){
                             cn = FALSE
                         }

                         if (!is.null(juncs) & length(juncs)>0){
                             if ("cn" %in% colnames(values(juncs))){
                                 jadd = which(values(juncs)$cn > 0)
                             } else {
                                 jadd = seq_along(juncs)
                             }
                         }

                         ## if there is tile, add tile
                         if (!is.null(tile) & length(tile)>0 & !is.null(juncs) & length(juncs)>0){
                             ## self$addSegs(c(gr.stripstrand(tile[,c()]),
                             ##                gr.stripstrand(unlist(juncs[jadd])[,c()])))
                             ## self$addJuncs(juncs)
                             self$addSegs(tile, cn=cn)$addJuncs(juncs, cn=cn)
                         } else if (!is.null(tile) & length(tile)>0){
                             self$addSegs(tile, cn=cn)
                         } else if (!is.null(juncs) & length(jadd)>0){
                             ## if empty, ignore these GRanges lists
                             self$addJuncs(juncs[jadd], cn=cn)
                         }
                         return(self)
                     },

                     simplify = function(mod=TRUE){
                         ## if two or more segment are only connected by ref edges
                         ## and they have the same copy number
                         ## merge them into one node
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Merge all pairs of noded only connected by reference edge.")
                         }

                         if (length(private$segs)==0 | is.null(private$es)){
                             if (verbose){
                                 warning("Empty graph. Cannot be simpler.")
                             }
                             return(self)
                         }

                         ## get the part of the graph where the nodes are
                         ## those at least one side is connecting a single reference edge
                         ## find the ones with indgree==1 & outdgree==1 & aberrant degree==0
                         node.dt = gr2dt(private$segs)
                         node.dt[, id := seq_along(private$segs)]
                         node.dt[, ":="(in.d = private$es[, table(to)[as.character(1:nrow(node.dt))]],
                                        out.d = private$es[, table(from)[as.character(1:nrow(node.dt))]])]
                         node.dt[is.na(in.d), in.d:=0]
                         node.dt[is.na(out.d), out.d := 0]

                         ref.es = private$es[type=="reference"]
                         ref.es[, ":="(to.in.d = node.dt[to, in.d],
                                       from.out.d = node.dt[from, out.d])]
                         ref.es[to.in.d==1 & from.out.d==1]

                         if ("cn" %in% colnames(node.dt)){
                             ref.es[, ":="(to.cn = node.dt[to, cn],
                                           from.cn = node.dt[from, cn])]
                             ref.es[, leveled := from.cn==to.cn]
                         } else {
                             ref.es[, leveled := TRUE]
                         }

                         node.dt[, ":="(ref.in.d = ref.es[which(leveled),
                                                          table(to)[as.character(1:nrow(node.dt))]],
                                        ref.out.d = ref.es[which(leveled),
                                                           table(from)[as.character(1:nrow(node.dt))]])]
                         node.dt[is.na(ref.in.d), ref.in.d := 0]
                         node.dt[is.na(ref.out.d), ref.out.d := 0]

                         ## now with candidate ref edges, merge the two nodes one by one
                         nix.to.simplify = ref.es[which(leveled), sort(unique(c(from, to)))]
                         ref.A = Matrix::sparseMatrix(i = ref.es[which(leveled), from],
                                                      j = ref.es[which(leveled), to],
                                                      x = 1,
                                                      dims=dim(self$get.adj()))
                         ref.G = igraph::graph_from_adjacency_matrix(ref.A)

                         ## find all weak components here
                         ref.comps = igraph::components(ref.G, "weak")

                         ## within each cluster, identify the start and end of the Hamilton path
                         node.dt[, map := ref.comps$membership]
                         node.dt[, north.bound := ref.in.d == 0]
                         node.dt[, south.bound := ref.out.d == 0]

                         ## check point: every weakly connected component should have one head and one tail!!!!!
                         if (!node.dt[, .(sum(north.bound)==1, sum(south.bound)==1), by = map][, all(V1) & all(V2)]){
                             stop("You found a bug in gg$simplify!")
                         }

                         ## start merging
                         segs = private$segs
                         segs$map = node.dt[, map]
                         new.segs = gr.simplify(segs, field="map")
                         new.head = node.dt[north.bound==TRUE, setNames(map, id)]
                         new.tail = node.dt[south.bound==TRUE, setNames(map, id)]
                         if ("cn" %in% colnames(private$es)){
                             new.es = private$es[, .(from = new.tail[as.character(from)],
                                                     to = new.head[as.character(to)],
                                                     type, cn)] ## assume there is always "type" in es
                         } else {
                             new.es = private$es[, .(from = new.tail[as.character(from)],
                                                     to = new.head[as.character(to)],
                                                     type)] ## assume there is always "type" in es
                         }
                         new.es = new.es[!is.na(from) & !is.na(to)]

                         ## new.et = etype(new.segs, new.es, force=T, both=T)
                         ## new.es = new.et$es
                         ## new.segs = new.et$segs

                         if (mod){
                             private$reset()
                             private$gGraphFromScratch(new.segs, new.es)
                             return(self)
                         } else {
                             out = gGraph$new(segs = new.segs, es = new.es)
                             if (class(self)[1]=="bGraph"){
                                 out = bGraph$new(gG = out)
                             }
                             return(out)
                         }
                     },

                     ## TODO: why do we miss GL contigs in here?????
                     decouple = function(mod=TRUE){
                         ## DEBUG
                         ## NOTE: the problem is here, we added extra copies of certain nodes!!
                         if (verbose <- getOption("gGnome.verbose")){
                             message("When there's overlapping nodes, break them down and reconnect.")
                         }

                         if (!"loose" %in% colnames(values(private$segs))){
                             private$segs = etype(private$segs, private$es, both=TRUE)$segs
                         }

                         ## ASSUMPTION: nodes of a gGraph are always skew-symmetric
                         if (isDisjoint(private$segs %Q% (strand=="+" & loose==FALSE))){
                             return(self)
                         }

                         ## start decoupling
                         ## temporary fix
                         ## This is wrong because same cn wouldn't break a segment
                         ## pos.cov = coverage(private$segs %Q% (loose==FALSE & strand=="+"),
                         ##                    weight="cn")
                         ## neg.cov = coverage(private$segs %Q% (loose==FALSE & strand=="-"),
                         ##                    weight="cn")

                         ## pos.gr = tryCatch(as(pos.cov, "GRanges"),
                         ##                   error = function(e){
                         ##                       reg.chr = Sys.getenv("DEFAULT_REGULAR_CHR")
                         ##                       if (file.exists(reg.chr)){
                         ##                           reg.chr = fread(reg.chr)[, setNames(V2, V1)]
                         ##                       }
                         ##                       return(as(pos.cov[names(reg.chr)], "GRanges"))
                         ##                   })
                         ## neg.gr = tryCatch(as(neg.cov, "GRanges"),
                         ##                   error = function(e){
                         ##                       reg.chr = Sys.getenv("DEFAULT_REGULAR_CHR")
                         ##                       if (file.exists(reg.chr)){
                         ##                           reg.chr = fread(reg.chr)[, setNames(V2, V1)]
                         ##                       }
                         ##                       return(as(neg.cov[names(reg.chr)], "GRanges"))
                         ##                   })

                         ## strand(pos.gr) = "+"
                         ## strand(neg.gr) = "-"
                         ## segs = c(pos.gr, neg.gr)

                         ## segs$cn = segs$score; segs$score = NULL
                         ## segs = segs %Q% (!is.na(cn))


                         ## this messed CN!!!!!!!!!!!!!!!!!!!!!!
                         ## full.segs = gUtils::streduce(private$segs)
                         ## segs = gUtils::gr.breaks(full.segs, private$segs %Q% (loose == FALSE))
                         segs = disjoin(private$segs %Q% (loose==FALSE))
                         segs = gUtils::gr.val(segs,
                         (private$segs %Q% (loose==FALSE & strand=="+"))[,c("cn")],
                         val="cn", FUN=sum,
                         weighted=FALSE,
                         ignore.strand=FALSE)

                         ## NOTE: once breaking the segs, there will be a lot of ref edges missing
                         ## NOTE: this is because the 0 CN nodes are discarded!!
                         all.j = e2j(private$segs, private$es, etype = "aberrant")

                         if (mod==T){
                             self$karyograph(tile = segs, juncs = all.j, cn = TRUE, genome = seqinfo(segs))
                             return(self)
                         } else {
                             out = gGraph$new(tile = segs, juncs = all.j, cn = TRUE)
                             return(out)
                         }
                     },

                     add = function(gg,
                                    mod = FALSE,
                                    decouple = TRUE){
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Simply put two gGraphs together.")
                         }

                         if (!inherits(gg, "gGraph"))
                             stop("Error: Can only deal with addition of two gGraph objects.")

                         ## bare GRanges
                         new.segs = c(private$segs[,c()], gg$segstats[,c()])
                         common.segs.mc = intersect(colnames(values(private$segs)),
                                                    colnames(values(gg$segstats)))
                         if (length(common.segs.mc)>0){
                             values(new.segs) = rbind(values(private$segs)[, common.segs.mc],
                                                      values(gg$segstats)[, common.segs.mc])
                         }

                         ## only from and to are required
                         new.es = rbind(private$es[,.(from, to)],
                                        gg$edges[,.(from = from + length(private$segs),
                                                    to = to + length(private$segs))])
                         common.es.mc = setdiff(intersect(colnames(private$es),
                                                          colnames(gg$edges)),
                                                c("from", "to"))

                         if (mod){
                             private$gGraphFromScratch(segs = new.segs,
                                                       es = new.es)
                             if (decouple){
                                 self$decouple()
                             }
                             return(self)
                         } else {
                             gg = gGraph$new(segs = new.segs, es = new.es)
                             if (decouple){
                                 gg$decouple()
                             }
                             return(gg)
                         }
                     },

                     ## initialize from JaBbA output
                     jab2gg = function(jabba, regular=FALSE){
                         if (is.list(jabba)) {
                             if (!all(is.element(c("segstats", "adj", "ab.edges",
                                                   "purity", "ploidy", "junctions"),
                                                 names(jabba))))
                                 stop("The input is not a JaBbA output.")
                         } else if (is.character(jabba) & grepl(".rds$", jabba)){
                             if (file.exists(jabba)){
                                 jabba = readRDS(jabba)
                             }
                         } else {
                             stop("Error: Input must be either JaBbA list output or the RDS file name that contains it!")
                         }

                         ## make sure required mcol is filled
                         private$segs = jabba$segstats
                         if (all(is.na(s.sl <- seqlengths(jabba$segstats)))){
                             ## JaBbA's output segstats is without seqlengths!!
                             ## resort to junctions
                             if (!any(is.na(j.sl <- seqlengths(jabba$junctions)))){
                                 if (getOption("gGnome.verbose")){
                                     warning("No valid seqlengths found in $segstats, force to the same as $junctions.")
                                 }
                                 private$segs = gUtils::gr.fix(private$segs, jabba$junctions)
                             } else {
                                 if (getOption("gGnome.verbose")){
                                     warning("No valid seqlengths found anywhere in input, force to DEFAULT.")
                                 }
                                 default.sl = data.table::fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]
                                 private$segs = gUtils::gr.fix(private$segs, default.sl)
                             }
                         }

                         if ("edges" %in% names(jabba)){
                             private$es = as.data.table(jabba$edges)
                         } else {
                             ## DEBUG: hopefully this will deal with empty edges
                             private$es = as.data.table(which(jabba$adj>0, arr.ind=T))
                             colnames(private$es) = c("from", "to")
                         }
                         private$es = etype(private$segs, private$es)
                         ## private$g = igraph::make_directed_graph(t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))

                         ## MARCIN DEBUG: THIS MISSES WHOLE CHROMOSOMES WHICH
                         ## HAVE BOTH ZERO IN AND ZERO OUT DEGREE!!
                         private$segs$terminal  = !(1:length(private$segs) %in% private$es$from) | !(1:length(private$segs) %in% private$es$to)

                         ## private$abEdges = jabba$ab.edges
                         private$.ploidy = jabba$ploidy
                         private$.purity = jabba$purity

                         if (inherits(regular, "character")){
                             regularChr = seqinfo()
                         } else if (regular==T){
                             if (getOption("gGnome.verbose")){
                                 warning("Forcing regular chromosomes. Will try default. See `Sys.getenv('DEFAULT_REGULAR_CHR')`.")
                             }
                             regularChr = si2gr(data.table::fread(Sys.getenv('DEFAULT_REGULAR_CHR'))[, setNames(V2, V1)])
                             self$trim(regularChr)
                         }
                         return(self)
                     },

                     ## initialize from Weaver result
                     wv2gg = function(weaver){
                         ## DONE: get Weaver done!!!! GEt it done@!!!
                         ## input weaver: directory that contains three and only three files
                         ## named: REGION_CN_PHASE, SV_CN_PHASE, and SNP_CN_PHASE
                         if (!dir.exists(weaver)){
                             stop("Error: Invalid input weaver directory!")
                         }

                         if (!all(is.element(c("SV_CN_PHASE", "REGION_CN_PHASE"), dir(weaver))) ){
                             stop('Error: Need "SV_CN_PHASE" and "REGION_CN_PHASE".')
                         }

                         sl = fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]

                         region = data.table(read.delim(
                             paste(weaver, "REGION_CN_PHASE", sep="/"),
                             header = FALSE, sep = "\t"))

                         sv.fn = paste(weaver, "SV_CN_PHASE", sep="/")
                         if (file.size(sv.fn)>0){
                             sv = data.table(read.delim(sv.fn, header = FALSE, sep = "\t"))
                             names(sv) = c("chr1", "pos1", "side1", "allele1",
                                           "chr2", "pos2", "side2", "allele2",
                                           "cn", "unknown1", "unknown2", "timing", "class")[1:ncol(sv)]
                         }
                         else {
                             sv = NULL
                         }

                         ## define the columns
                         names(region) = c("seqnames", "start", "end", "acn", "bcn")
                         region[, cn := acn + bcn]
                         ## names(snp) = c("seqnames", "pos", "ref", "alt", "acn", "bcn")

                         ss = dt2gr(region)
                         ss = gr.fix(ss, sl)

                         ## get junctions
                         ## ALERT: in the file, +/- means right/left end of a segment
                         ## exactly reverse of what we define a junction
                         strmap = setNames(c("+", "-"), c("-", "+"))
                         ## sv.select = sv[!is.na(allele1) & !is.na(allele2)]
                         if (!is.null(sv)){
                             sv.select = sv[, which(cn>0)] ## makes more sense?
                             bps = c(
                                 dt2gr(
                                     sv[, .(seqnames = chr1,
                                            start = ifelse(side1=="-", pos1-1, pos1),
                                            end = ifelse(side1=="-", pos1-1, pos1),
                                            jix=.I, ii = 1,
                                            strand = strmap[side1])]),
                                 dt2gr(
                                     sv[, .(seqnames = chr2,
                                            start = ifelse(side2=="-", pos2-1, pos2),
                                            end = ifelse(side2=="-", pos2-1, pos2),
                                            jix=.I, ii = 2,
                                            strand = strmap[side2])]))
                             ## ALERT: nudge 1bp offset for only the "-" bp

                             ## sanity check, all raw.bp at this point should
                             ## locate at left/right boundary of segements
                             ss.ends = c(gr.start(ss), gr.end(ss))
                             if (any(!bps %^% ss.ends)){
                                 warning("Eligible SVs not matching segment ends!")
                             }

                             ## create junctions
                             junc = grl.pivot(split(bps, bps$ii))
                             toget = intersect(c("allele1", "allele2", "cn", "unknown1", "unknown2", "timing", "class"), colnames(sv))
                             values(junc) = sv[, toget, with=F]
                         }
                         else {
                             junc = NULL
                         }

                         ## edges and graph
                         ## ALERT!! ALERT!!
                         ## TODO: still can't use addxxx() functions in a chain
                         ## doing these two steps apart will result in breakpoint missing from
                         ## self$nullGGraph()$addSegs(ss)$addJuncs(junc)
                         ## private$abEdges = self$makeAbEdges()
                         self$karyograph(tile = ss, juncs = junc, cn = TRUE)

                         ## ALERT: Weaver out put is not balanced!
                         ## ## balancing it out
                         ## adj = self$get.adj()
                         ## ifl = Matrix::colSums(adj)
                         ## ofl = Matrix::rowSums(adj)
                         ## cns = private$segs$cn
                         ## ## NA segs conforms to larger of (ifl, ofl)
                         ## private$segs$cn[which(is.na(cns))] = pmax(ifl, ofl)[which(is.na(cns))]
                         ## cns = private$segs$cn
                         ## private$segs[which(cns<ifl)]
                         return(self)
                     },

                     ## initialize from Prego result
                     pr2gg = function(fn){
                         sl = fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]
                         ## ALERT: I don't check file integrity here!
                         ## first part, Marcin's read_prego
                         res.tmp = readLines(fn)
                         chrm.map.fn = gsub(basename(fn), "chrm.map.tsv", fn)

                         if (file.exists(chrm.map.fn)){
                             message(chrm.map.fn)
                             message("Seqnames mapping found.")
                             chrm.map = fread(chrm.map.fn)[,setNames(V1, V2)]
                         } else {
                             warning("Warning: No mapping seqnames info, will throw out all non 1:24 values.")
                         }

                         res = structure(lapply(split(res.tmp, cumsum(grepl("edges", res.tmp))),
                                                function(x) {
                                                    rd = read.delim(textConnection(x),
                                                                    strings = F,
                                                                    skip = 1,
                                                                    header = F,
                                                                    col.names = c("node1", "chr1",
                                                                                  "pos1", "node2",
                                                                                  "chr2", "pos2", "cn"))
                                                    if (exists("chrm.map")){
                                                        rd$chr1 = chrm.map[rd$chr1]
                                                        rd$chr2 = chrm.map[rd$chr2]
                                                    }
                                                    else {
                                                        rd = rd[which(rd$chr1 %in% as.character(1:24) &
                                                                      rd$chr2 %in% as.character(1:24)),]
                                                        rd$chr1 = gsub("24", "Y", gsub("23","X",rd$chr1))
                                                        rd$chr2 = gsub("24", "Y", gsub("23","X",rd$chr2))
                                                    }

                                                    return(rd)
                                                }),
                                         names = gsub(":", "", grep("edges", res.tmp, value = T)))
                         res[[1]]$tag = paste0(res[[1]]$node1, ":", res[[1]]$node2)

                         ## turn into our segstats
                         segstats = GRanges(res[[1]]$chr1,
                                            IRanges(res[[1]]$pos1,
                                                    res[[1]]$pos2),
                                            strand = "+",
                                            cn = res[[1]]$cn,
                                            left.tag = res[[1]]$node1,
                                            right.tag = res[[1]]$node2,
                                            loose=FALSE)
                         segstats = gr.fix(c(segstats, gr.flipstrand(segstats)), sl)
                         neg.ix = which(as.logical(strand(segstats) == "-"))
                         ## tag1 is the 3' end
                         tag1 = segstats$right.tag
                         tag1[neg.ix] = segstats$left.tag[neg.ix]
                         ## tag2 is the 5' end
                         tag2 = segstats$left.tag
                         tag2[neg.ix] = segstats$right.tag[neg.ix]

                         hb = hydrogenBonds(segstats)
                         hb.map = hb[, setNames(from, to)]

                         ## adjacency in copy number
                         adj.cn = matrix(0, nrow = length(segstats), ncol = length(segstats),
                                         dimnames = list(tag1, tag2))
                         adj.cn[cbind(res[[2]]$node1, res[[2]]$node2)] = res[[2]]$cn
                         adj.cn[cbind(res[[2]]$node2, res[[2]]$node1)] = res[[2]]$cn
                         adj.cn[cbind(res[[3]]$node1, res[[3]]$node2)] = res[[3]]$cn
                         adj.cn[cbind(res[[3]]$node2, res[[3]]$node1)] = res[[3]]$cn

                         ## ## adjacency in edge type
                         ## adj.type = matrix("", nrow = length(segstats), ncol = length(segstats),
                         ##                   dimnames = list(tag1, tag2))
                         ## adj.type[cbind(res[[2]]$node1, res[[2]]$node2)] = "reference"
                         ## adj.type[cbind(res[[2]]$node2, res[[2]]$node1)] = "reference"
                         ## adj.type[cbind(res[[3]]$node1, res[[3]]$node2)] = "aberrant"
                         ## adj.type[cbind(res[[3]]$node2, res[[3]]$node1)] = "aberrant"

                         ## create es
                         ed = as.data.table(which(adj.cn>0, arr.ind=T))
                         colnames(ed) = c("from", "to")
                         ed[, ":="(cn = adj.cn[cbind(from, to)])]
                         ed = etype(segstats, ed)
                         private$gGraphFromScratch(segs = segstats,
                                                   es = ed,
                                                   purity = 1)
                         return(self)
                     },

                     ## public methods
                     ## I/O
                     print = function(){
                         cat('A gGraph object.\n')
                         ## cat('Based on reference genome: ')
                         ## cat(genome(private$segs))
                         cat('\n\n')
                         cat('Total segmentation:')
                         if ("loose" %in% colnames(values(private$segs))){
                             cat(length(private$segs %Q% (loose==F & strand=="+")))
                         } else {
                             ## ALERT!!! TODO!!! This means we have to make sure if there is
                             ## loose end, it must be labeled in the nodes.
                             cat(length(private$segs %Q% (strand=="+")))
                         }
                         cat('\n\n')
                         cat('Edge counts:\n')
                         if (is.null(private$es)){
                             cat('None')
                         } else if (nrow(private$es)==0){
                             cat('None')
                         } else {
                             if (!"type" %in% colnames(private$es)){
                                 private$es = etype(private$segs, private$es, force=T)
                             }
                             print(private$es[, table(type)/2])
                         }
                     },

                     ## TODO: find better default settings
                     plot = function(pad=1e3, colorful=FALSE, ...){
                         td = self$gg2td()
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
                         win = self$window(pad)
                         ## decide X gap on the fly
                         plot(td, win, links = private$junction)
                     },

                     window = function(pad=0){
                         if ("cn" %in% colnames(values(private$segs))){
                             win = gUtils::streduce((private$segs %Q% (!is.na(cn))) + pad)
                         } else {
                             win = gUtils::streduce((private$segs) + pad)
                         }
                         return(win)
                     },

                     ## TODO: find better default settings
                     ## we want to layout anything!
                     layout = function(){
                         if (length(private$segs)==0 | is.null(private$es)){
                             return(NULL)
                         }
                         if (!inherits(private$g, "igraph")){
                             self$get.g(force=TRUE)
                         }
                         ## TODO: return the plot value
                         ## TODO: decide best visual parameters depend on the size of the graph!!!
                         vcolor = ifelse(strand(private$segs)=="+", "salmon", "skyblue")
                         c3 = setNames(skitools::brewer.master(n = 3, palette = "Set1"),
                                       nm = c("aberrant", "loose", "reference"))
                         ed = private$es
                         if (!is.null(ed)){
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
                         } else {
                             plot.igraph(private$g,
                                         ## layout
                                         layout = layout_with_gem,
                                         ## vertex pars
                                         vertex.size=log(private$segs$cn,1.4), vertex.color= vcolor,
                                         vertex.shape="circle", vertex.label.cex = 0.75,
                                         vertex.frame.color=NA, vertex.label.color = "black",
                                         )
                         }
                         return(NULL)
                     },

                     ## TODO: make it informative
                     summary = function(){
                         summ = "This is a gGraph object."
                         return(summ)
                     },

                     length = function(){
                         ## DONE
                         if (length(private$segs)==0){
                             return(0L)
                         }
                         if (is.null(private$partition)){
                             private$partition = self$components()
                         }
                         return(private$parition$no)
                     },
                     ##

                     gg2td = function(seg.col, ...){
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Create gTrack for static genome browser-style viz.")
                         }
                         if (length(private$segs)==0 | length(private$es)==0){
                             if (verbose){
                                 warning("Nothing to plot!")
                             }
                             return(NULL)
                         }
                         if (!"loose" %in% colnames(values(private$segs)) |
                             !"type" %in% colnames(private$es)){
                             et = etype(private$segs, private$es, force=T, both=T)
                             private$segs = et$segs
                             private$es = et$es
                         }

                         ## DONE: allow users to define extra fields to annotate segs or edges!!!
                         ## DONE: replicate classic JaBbA viz
                         ## plotting segments
                         ## if loose, make it white, lift it up
                         ss = private$segs
                         ed = private$es

                         if (!is.null(ed))
                         {
                             ## set edge apperances
                             ## lwd, lty, col, cex.arrow, v, not.flat, h, dangle.w
                             if (!is.element("cn", colnames(ed))) {
                                 if (verbose){
                                     warning("Edges has no copy number yet.")
                                 }
                                 ed[, cn := 1]
                             }
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
                         }

                         ## DONE: handle so/si-less edges, omit for now

                         ## set segment apperances
                         ## if loose, change its cn to slightly higher than it's incident node
                         if (any(ss$loose==T)){
                             lid = which(ss$loose)
                             ## find partner indices for loose ends
                             pid = sapply(lid,
                                          function(i) ed[from==i | to==i,
                                                         ifelse(from==i, to, from)],
                                          simplify=T)
                             if (is.list(pid)){
                                 pid = unlist(pid)
                             }
                             ss$cn[lid] = ss$cn[pid]*1.2
                         }

                         ## col, border, ywid
                         ss$col = ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5))
                         ss$border = ifelse(ss$loose, ss$col, alpha("black", 0.5))
                         ss$ywid = ifelse(ss$loose, 0.001, 0.8)

                         if ("cn" %in% colnames(values(ss))){
                             gt = gTrack(ss, y.field="cn", edges=ed, name="CN", angle=0, ...)
                         } else {
                             gt = gTrack(ss, edges=ed, name="CN", angle=0, ...)
                         }
                         return(gt)
                     },

                     json = function(filename='.',
                                     maxcn=100,
                                     maxweight=100){
                         self$gg2js(filename, maxcn, maxweight, save=TRUE)
                     },

                     html = function(filename='.',
                                     gGnome.js=Sys.getenv("DEFAULT_GGNOMEJS"),
                                     maxcn=100,
                                     maxweight=100,
                                     invoke=FALSE){
                         "Dump JSON into a copy of gGnome.js for quick viz"
                         if (!dir.exists(gGnome.js)){
                             message("No gGnome.js repository found on your system.")
                             stop("Get from https://github.com/mskilab/gGnome.js")
                         }

                         if (dir.exists(filename)){
                             basedir = filename
                             filename = paste(basedir, "gGnome.js", "json",
                                              "data.json", sep="/")
                         } else if (grepl(".json$", filename)){
                             basedir = dirname(filename)
                             filename = paste(basedir, "gGnome.js", "json",
                                              basename(filename), sep="/")
                         } else {
                             basedir = filename
                             system(paste("mkdir -p", filename))
                             filename = paste(basedir, "gGnome.js", "json",
                                              "data.json", sep="/")
                         }

                         ## copy the whole directory
                         system(paste("cp",
                                      system.file("extdata", "gGnome.js", package="gGnome"),
                                      basedir
                                      ))

                         ## generating the JSON
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Writing your JSON file to:", filename)
                         }
                         self$gg2js(filename, maxcn, maxweight, trim, all.js=TRUE)

                         if (invoke){
                             system(paste0(paste0(basedir, "/gGnome.js"), "/start.sh"))
                         }
                         return(normalizePath(filename))
                     },

                     gg2js = function(filename='.',
                                      maxcn=100,
                                      maxweight=100,
                                      save = TRUE,
                                      settings = NULL,
                                      no.y = FALSE){
                         if (save){
                             if (grepl('\\.js(on)*$', filename)){
                                 ## if json path was provided
                                 basedir = dirname(filename)
                             }
                             else if (filename==".") {
                                 ## default path was provided
                                 basedir = './'
                                 filename = "data.js"
                             }
                             else {
                                 ## a directory was provided
                                 basedir = filename
                                 filename = paste(filename, 'data.json', sep = '/')
                             }

                             if (!file.exists(basedir)) {
                                 message('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }
                         }

                         if (verbose <- getOption("gGnome.verbose")){
                             message("Create json file for interactive visualization.")
                         }

                         if (is.null(settings)){
                             settings = list(y_axis = list(name = "copy number"))
                         }
                         qw = function(x) paste0('"', x, '"') ## quote

                         ## range of CN
                         ymin=0
                         ymax=maxcn

                         ## ALERT: for a clean viz for now, only contain regular chromosomes
                         ## ALERT: for a clean viz for now, only contain regular chromosomes
                         ## ADDED BY MARCIN: define regularChr
                         ## EDIT BY XT: now we define env default values
                         regular.sl =
                             fread(Sys.getenv("DEFAULT_REGULAR_CHR"))[, setNames(V2, V1)]
                         regsegs.ix = which(as.character(seqnames(private$segs))
                                            %in% names(regular.sl))

                         loose.ix = which(private$segs$loose==TRUE)

                         ed = copy(private$es) ## otherwise change by reference!
                         ## construct intervals
                         node.dt = data.table(oid = which(as.logical(strand(private$segs)=="+")))
                         ## node.dt[, rid := seq_along(private$segs)[-oid][match(private$segs[-oid],
                         ##                                                      gUtils::gr.flipstrand(
                         ##                                                                  private$segs[oid]
                         ##                                                              ))]]
                         hb = hydrogenBonds(private$segs)
                         hb.map = hb[, setNames(from, to)]
                         ## MOMENT
                         node.dt[, rid := hb.map[as.character(oid)]]

                         node.dt = node.dt[oid %in%
                                           which(private$segs$loose==FALSE &
                                                 as.character(seqnames(private$segs))
                                                 %in% names(regular.sl))]

                         node.dt[, iid := 1:.N]
                         setkey(node.dt, "iid")
                         node.dt[, ":="(chr = as.character(seqnames(private$segs[oid])),
                                        start = start(private$segs[oid]),
                                        end = end(private$segs[oid]))]
                         node.map = node.dt[, c(setNames(iid, oid),
                                                setNames(iid, rid))]


                         node.dt[, y := private$segs$cn[oid]]

                         ## TODO: do not assume things are paired up
                         ## do not assume the cn field in the segs is correct
                         node.dt[, title := paste(iid, paste0("(",oid,"|",rid,")"))]
                         node.dt[, type := "interval"]
                         node.dt[, strand := "*"]

                         node.dt.both = rbind(node.dt[, .(nid = oid, iid,
                                                          chr, start, end, y,
                                                          title, type, strand="+")],
                                              node.dt[, .(nid = rid, iid,
                                                          chr, start, end, y,
                                                          title, type, strand="-")])
                         setkey(node.dt.both, "nid")

                         ## NODE.JSON
                         node.json = node.dt[, .(iid,
                                                 chromosome = chr,
                                                 startPoint = start,
                                                 endPoint = end,
                                                 y,
                                                 title,
                                                 type,
                                                 strand)]

                         ## TMPFIX: remove NA edges .. not clear where these are coming from
                         ## but likely the result of trimming / hood, but then it's not balanced
                         ## mapping from type field to label in json
                         eType = setNames(c("REF", "ALT", "LOOSE"),
                                          c("reference", "aberrant", "loose"))

                         ## some edges are out of the scope of regular chrs
                         e.na.ix = ed[, which(is.na(from) |
                                              is.na(to) |
                                              !(from %in% regsegs.ix) |
                                              !(to %in% regsegs.ix))]
                         ed.na = ed[e.na.ix, ]

                         ## if any edge left, process
                         if (nrow(ed)-length(e.na.ix)>0){
                             if (any(e.na.ix)){
                                 ed = ed[-e.na.ix, ]
                             }

                             ## edge's unique identifier
                             ed[, ":="(eid = paste(from, to, sep="-"),
                                       reid = paste(hb.map[as.character(to)],
                                                    hb.map[as.character(from)],
                                                    sep="-"))]

                             ## ALERT: bc strandlessness, I only retained half of the edges
                             ## to map edges in gwalks, we will need strandedness,
                             ## so will retain everything
                             ed[,":="(soStr = as.character(strand(private$segs[from])),
                                      siStr = as.character(strand(private$segs[to])))]

                             ## compute eclass
                             ed[, ":="(ix = 1:.N,
                                       rix = match(reid, eid))]
                             ed[, unique.ix := ifelse(rix>=ix, paste(ix, rix), paste(rix, ix))]
                             ed[, eclass := as.numeric(as.factor(unique.ix))]
                             ed[, iix := 1:.N, by=eclass]

                             ## metadata of edges
                             ed[, ":="(so = node.map[as.character(from)],
                                       si = node.map[as.character(to)],
                                       so.str = ifelse(soStr=="+",1,-1),
                                       si.str = ifelse(siStr=="+",1,-1),
                                       title = "",
                                       type = eType[type],
                                       weight = cn,
                                       cid = eclass)]

                             ed[, ":="(source = so*so.str,
                                       sink = -si*si.str)]

                             ## NOTE: I will never ever manually create/parse a JSON from string myself in my lift
                             ## ppl wrote JSON format to make things standardized and pain-free to use
                             ## Let's trust ppl

                             ## EDGE.JSON
                             ed.json = ed[iix==1, ## only need half of edges
                                          .(cid,
                                            source,
                                            sink,
                                            title,
                                            type,
                                            weight)]

                         } else {
                             ed.json = data.table(cid = numeric(0),
                                                  source = numeric(0),
                                                  sink = numeric(0),
                                                  title = character(0),
                                                  type = character(0),
                                                  weight = numeric(0))
                         }

                         ed.json = ed.json[!is.na(cid)]

                         gg.js = list(intervals = node.json, connections = ed.json)

                         if (no.y){
                             settings$y_axis = list(visible=FALSE)
                             gg.js$intervals[, y := NULL]
                         }

                         if (!is.null(settings)){
                             gg.js = c(list(settings = settings), gg.js)
                         }

                         if (save){
                             if (verbose <- getOption("gGnome.verbose")){
                                 message("Saving JSON to: ", filename)
                             }
                             jsonlite::write_json(gg.js, filename,
                                                  pretty=TRUE, auto_unbox=TRUE, digits=4)
                             return(normalizePath(filename))
                         } else {
                             return(gg.js)
                         }
                     },

                     ## dicing up the graph
                     components = function(mc.cores=1){
                         ## create a sticky graph where pairs of +/- are connected by hydro edges
                         stickyG = private$g
                         hB = hydrogenBonds(private$segs)
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

                     ## DONE:
                     ## if na.rm==F, balanced graph's subgraph should always be balanced!!!!!
                     ## TODO:
                     ## record old segs labels
                     subgraph = function(v=numeric(0),
                                         na.rm=T,
                                         mod=FALSE){
                         "Given a numeric vector of vertices, \
                         change this gGraph to its subgraph consists of only these vertices."
                         if (length(v)==0){
                             ## nothing provided, nothing happens
                             return(self)
                         }
                         else if (is.numeric(v)){
                             ## at least they are num
                             if (!is.integer(v)){
                                 ## if not integer, convert
                                 v = as.integer(v)
                             }

                             if (!all(v %in% seq_along(private$segs))){
                                 v = v[which(v %in% seq_along(private$segs))]
                                 warning("Some v subscripts out of bound! Ignore!")
                             }

                             if (length(loose.v <- intersect(v, which(private$segs$loose==T)))>0){
                                 warning("Some v is loose end. Ignore!")
                                 v = setdiff(v, loose.v)
                             }

                             ## DONE: also recover v's missing reverse complements
                             hB = hydrogenBonds(private$segs)
                             vid = sort(unique(c(v, hB[from %in% v, to], hB[to %in% v, from])))

                             ## get the subgraph
                             newSegs = private$segs[vid]
                             newSegs$last.id = vid

                             newId = setNames(seq_along(vid), vid)
                             newEs = private$es[from %in% vid & to %in% vid]
                             newEs[, ":="(last.from = from,
                                          last.to = to,
                                          from=newId[as.character(from)],
                                          to=newId[as.character(to)])]

                             ## ## DONE: use "fillin" function on the graph if na.rm=F
                             ## jIdx = which(grl.in(private$junction, newSegs, only=T))
                             ## newJuncs = private$junction[unique(jIdx)]

                             if (mod==T){
                                 private$gGraphFromScratch(segs=newSegs,
                                                           es=newEs,
                                                           ## junc=newJuncs,
                                                           purity=private$.purity)
                                 if (na.rm==F){
                                     self$fillin()
                                 }
                                 return(self)
                             }
                             else {
                                 out = gGraph$new(segs=newSegs,
                                                  es=newEs,
                                                  ## junctions=newJuncs,
                                                  purity=private$.purity)
                                 if (na.rm==F){
                                     out$fillin()
                                 }
                                 return(out)
                             }
                         }
                         else {
                             stop("Error: Invalid input.")
                         }
                     },

                     clusters = function(v = numeric(0), ## the vertices list we are looking at
                                         mode = c("weak", "strong"),
                                         use.hb=TRUE){
                         if (length(v)>0){
                             gg = self$subgraph(v, na.rm=FALSE)
                         } else {
                             v = seq_along(private$segs)
                             gg = self
                         }

                         stickyG = gg$G
                         hB = hydrogenBonds(gg$segstats)
                         ## update es and g
                         stickyG = add_edges(stickyG, t(as.matrix(hB[, .(from, to)])))
                         ## compute clusters
                         cl = igraph::clusters(stickyG, mode = mode)
                         if (!setequal(v, seq_along(private$segs))){
                             memb = rep(as.numeric(NA), length(private$segs))
                             memb[v] = cl$membership
                             cl$membership = memb
                         }
                         private$partition = cl
                         ## MOMENT
                         out = lapply(seq_len(cl$no),
                                      function(no){
                                          this.v = which(cl$membership==no)
                                          this.sg = self$subgraph(this.v)
                                          return(this.sg)
                                      })
                         names(out) = seq_len(cl$no)
                         return(out)
                     },

                     ## DONE!!!!!!
                     ## the idea of loose end: accesorries, only exist to BALANCE the graph
                     ## make them transient
                     fillin = function(mod=TRUE){
                         "Fill in the missing copies of edges to make the graph balanced."
                         if (self$isBalance()){
                             if (verbose <- getOption("gGnome.verbose")){
                                 message("Already a balanced graph.")
                             }
                             return(self)
                         }

                         ## GOAL: make loose ends a very free thing, add it, remove it, fuse a
                         ## pair of them or convert to a terminal feature.
                         adj = self$get.adj()
                         ifl = Matrix::colSums(adj)
                         ofl = Matrix::rowSums(adj)
                         cns = private$segs$cn

                         ## sanity check: edge copy number sum should eq adj
                         inE = data.table(toid=seq_along(private$segs))
                         tmp.inE = private$es[, .(cn = sum(cn)), by=to]
                         setkey(tmp.inE, "to")
                         inE[, cn := tmp.inE[.(toid), cn]]
                         inE[, cn := ifelse(is.na(cn), 0, cn)]

                         outE = data.table(fromid=seq_along(private$segs))
                         tmp.outE = private$es[, .(cn = sum(cn)), by=from]
                         setkey(tmp.outE, "from")
                         outE[, cn := tmp.outE[.(fromid), cn]]
                         outE[, cn := ifelse(is.na(cn), 0, cn)]

                         if (!all(inE[, setNames(cn, toid)] == ifl) | !all(outE[, setNames(cn, fromid)] == ofl)){
                             stop("Error: Adjacency not matching edges table!")
                         }

                         if (!"loose" %in% colnames(values(private$segs))){
                             private$segs = etype(private$segs, private$es, both=TRUE)$segs
                         }

                         ## next we determine if it is feasible to fill the slacks
                         ## test if inSum>cns | outSum>cns
                         if (any(ifl>cns | ofl>cns, na.rm = TRUE)){
                             warning("Infeasible graph!!")
                             return(self)
                         }
                         else {
                             ## node.cn = data.table(id = seq_along(private$segs))
                             node.cn = data.table(id = which(private$segs$loose==FALSE))
                             node.cn[, cn := cns[id]]
                             node.cn[, loose.out := (cns - ofl)[id]]
                             node.cn[, loose.in := (cns - ifl)[id]]

                             sl = seqlengths(private$segs)
                             ## right end of terminal
                             node.cn[, str := as.character(strand(private$segs[id]))]
                             node.cn[, east.bound := end(private$segs[id]) ==
                                           sl[as.character(seqnames(private$segs[id]))]]
                             node.cn[, west.bound := start(private$segs[id]) == 1]

                             ## telomeres don't have to have loose ends
                             node.cn[str=="+" & west.bound==TRUE, loose.in := 0]
                             node.cn[str=="+" & east.bound==TRUE, loose.out := 0]
                             node.cn[str=="-" & west.bound==TRUE, loose.out := 0]
                             node.cn[str=="-" & east.bound==TRUE, loose.in := 0]

                             ## before any action, if nothing to be filled in, then stop
                             if (node.cn[, !any(loose.in>0, na.rm=T)] |
                                 node.cn[, !any(loose.out>0, na.rm=T)]){
                                 return(self)
                             }

                             loose.ix = which(private$segs$loose)
                             seg.ix = which(!private$segs$loose)

                             new.loose.in = node.cn[
                                 loose.in>0,
                                 gr.start(private$segs[id], ignore.strand=FALSE)]
                             values(new.loose.in) = NULL
                             values(new.loose.in)$cn = node.cn[loose.in>0, loose.in]
                             values(new.loose.in)$loose = TRUE
                             values(new.loose.in)$terminal = TRUE
                             new.loose.in.es = node.cn[
                                 loose.in>0, .(from = .I,
                                               to=id,
                                               type="loose",
                                               cn = loose.in)]

                             ## construct GR for new loose ends required
                             new.loose.out = node.cn[
                                 loose.out>0,
                                 gr.end(private$segs[id], ignore.strand=FALSE)]
                             values(new.loose.out) = NULL
                             values(new.loose.out)$cn = node.cn[loose.out>0, loose.out]
                             values(new.loose.out)$loose = TRUE
                             values(new.loose.out)$terminal = TRUE
                             new.loose.out.es = node.cn[
                                 loose.out>0, .(from=id,
                                                to=.I,
                                                type="loose",
                                                cn = loose.out)]

                             new.loose = c(new.loose.in, new.loose.out)
                             new.loose.in.ix = setNames(length(private$segs) +
                                                        seq_along(new.loose.in),
                                                        seq_along(new.loose.in))
                             new.loose.out.ix = setNames(length(private$segs) +
                                                         length(new.loose.in) +
                                                         seq_along(new.loose.out),
                                                         seq_along(new.loose.out))

                             ## append new loose ends
                             seg.dt = rbind(gr2dt(private$segs),
                                            gr2dt(new.loose),
                                            fill=TRUE)
                             new.segs = dt2gr(seg.dt)
                             ## new.segs = c(private$segs[, c("cn", "loose")],
                             ##              new.loose[, c("cn", "loose")])

                             new.es = rbind(new.loose.out.es[
                               , .(from,
                                   to = new.loose.out.ix[as.character(to)],
                                   type, cn)],
                                 new.loose.in.es[
                                   , .(from = new.loose.in.ix[as.character(from)],
                                       to,
                                       type, cn)])
                             new.es = rbind(private$es[, .(from, to, type, cn)], new.es)

                             if (mod==FALSE){
                                 new.gg = gGraph$new(segs = new.segs,
                                                     es = new.es)
                                 if (new.gg$isBalance()){
                                     return(bGraph$new(new.gg))
                                 } else {
                                     stop("Fillin does not produce balanced graph!")
                                 }
                             } else {
                                 private$segs = new.segs
                                 private$es = new.es
                                 private$reset()
                                 return(self)
                             }
                         }
                     },

                     trim = function(gr=NULL, mod=FALSE){
                         "Return a trimmed subgraph that all the nodes are wihtin the range defined by gr."
                         ## DONE
                         ## if input gr is super set of private$segs, do nothing!
                         ## Only returning new obj
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Given a GRanges, return the trimmed subgraph overlapping it.")
                         }

                         if (is.null(gr)){
                             return(self)
                         }

                         gr = gr.fix(gr, private$segs)
                         gr = streduce(gr)

                         segs = private$segs
                         ov = gr.findoverlaps(segs, gr)

                         nss = ov

                         ## old segments in corresponding order
                         oss = segs[nss$query.id]

                         strand(nss) = strand(segs)[nss$query.id]
                         nss$eq = nss == oss
                         nss$left = start(ov)==start(oss)
                         nss$right = end(ov)==end(oss)
                         nss$internal = !nss$left & !nss$right

                         mcols(nss) = cbind(mcols(nss), mcols(oss)) ## carry over the metadata

                         ## map the edges
                         if (nrow(private$es)==0){
                             nes = private$es
                         }
                         else {
                             nss.dt = gr2dt(nss)[, nid := 1:.N]
                             nes = private$es[from %in% nss$query.id | to %in% nss$query.id,
                                              .(from, to, cn, type)]

                             ## left side of a + node receives its incoming edges
                             e.in = rbind(nss.dt[eq==TRUE,
                                                 .(oid=query.id, receive = nid)],
                                          nss.dt[eq==FALSE & left==TRUE & strand=="+",
                                                 .(oid=query.id, receive = nid)],
                                          nss.dt[eq==FALSE & right==TRUE & strand=="-",
                                                 .(oid=query.id, receive = nid)])
                             e.out = rbind(nss.dt[eq==TRUE,
                                                  .(oid=query.id, send = nid)],
                                           nss.dt[eq==FALSE & right==TRUE & strand=="+",
                                                  .(oid=query.id, send = nid)],
                                           nss.dt[eq==FALSE & left==TRUE & strand=="-",
                                                  .(oid=query.id, send = nid)])

                             setkey(e.in, "oid")
                             setkey(e.out, "oid")

                             ## if an old node lost its end, it will be NA after mapping
                             ## NOTE: don't forget the keyed query of data.table
                             new.es = cbind(nes[, .(from = e.out[.(from), send],
                                                    to = e.in[.(to), receive])],
                                            nes[, !c("from", "to")])
                         }

                         ## ov should be the only ranges of the returned graph
                         ## ov might be duplicated since we allow overlapping nodes in gGraph now
                         if (any(duplicated(nss[, c()]))){
                             nr.nss = unique(nss[, c()])
                             nmatch = data.table(nid = seq_along(nss),
                                                 nr.nid = match(nss[,c()], nr.nss))
                             if ("cn" %in% colnames(values(nss))){
                                 if (verbose){
                                     warning("Only 'cn' field is carried over.")
                                 }
                                 nmatch[, cn := nss$cn]
                                 nr.nss$cn = nmatch[, .(cn=sum(cn)), by=nr.nid][seq_along(nr.nss), cn]
                             }
                         }

                         ## reorder so easier for human reading
                         ord.nss = nss %Q% (order(loose, strand, seqnames, start))
                         nmatch = data.table(nid = seq_along(nss),
                                             ord.nid = match(nss, ord.nss))
                         setkey(nmatch, "nid")

                         new.es[, ":="(from = nmatch[.(from), ord.nid],
                                       to = nmatch[.(to), ord.nid])]

                         ## finally, recreate the trimmed graph
                         if (mod){
                             private$gGraphFromScratch(segs = ord.nss,
                                                       es = new.es)
                         } else {
                             newSg = gGraph$new(segs=ord.nss,
                                                es=new.es)
                             if (!newSg$isBalance()){
                                 ## NOTE: why do I have to reassign it here??
                                 newSg = newSg$fillin(mod=TRUE)
                             }
                             if (newSg$isBalance()){
                                 if (inherits(self, "bGraph")){
                                     newSg = as(newSg, "bGraph")
                                 }
                             }
                             return(newSg)
                         }
                     },

                     get.g = function(force=FALSE){
                         if (!is.null(private$g) & !force){
                             return(private$g)
                         } else {
                             private$g = igraph::make_directed_graph(
                                                     t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))
                             return(private$g)
                         }
                     },

                     get.adj = function(flat=FALSE){
                         if (is.null(private$g)){
                             self$get.g()
                         }
                         adjMat = igraph::as_adj(private$g)
                         if (flat) {
                             return(adjMat)
                         } else {
                             if (is.element("cn", colnames(private$es))){
                                 adjMat[private$es[,cbind(from, to)]] = private$es$cn
                             }
                             return(adjMat)
                         }
                     },
                     ## some query functions
                     hood = function(win,
                                     d=NULL,
                                     k=NULL,
                                     pad=0,
                                     bagel=FALSE,
                                     mod = FALSE,
                                     ignore.strand=T,
                                     verbose=FALSE){
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Get the trimmed subgraph around a given GRanges within a distance on the graph.")
                         }

                         if (ignore.strand){
                             win = gr.stripstrand(win)
                         }

                         ## DONE: what to do when win is larger than segs?????
                         ## ans: return self
                         if (length(setdiff(streduce(private$segs), win))==0){
                             return(self)
                         }

                         ## overlapping window and segs, removing loose ends
                         interGr = gr.findoverlaps(private$segs, win, ignore.strand=ignore.strand)
                         lid = which(private$segs$loose==T)
                         interGr = interGr %Q% (!query.id %in% lid)
                         qix = interGr$query.id

                         if (is.null(k)){
                             ## DONE!!!
                             ## no k, use distance
                             if (is.null(d) | d < 0){
                                 stop("Must provide either valid k or d.")
                             }

                             ## blend window with segs
                             win = gr.fix(win, private$segs)## fix seqinfo
                             ss = tryCatch(c(private$segs[private$segs$loose == F, c()],
                                             win[, c()]), error = function(e) NULL)

                             if (is.null(ss)){
                                 ss = grbind(c(private$segs[private$segs$loose == FALSE, c()],
                                               win[, c()]))
                             }

                             if (ignore.strand){
                                 ss = gr.stripstrand(ss)
                             }

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
                             if (any(s.close)){
                                 out = c(out,
                                         GenomicRanges::flank(seg.s[s.close],
                                                              -(d-min.s[s.close]))
                                         )
                             }

                             if (any(e.close)){
                                 out = c(out,
                                         GenomicRanges::shift(flank(seg.e[e.close],
                                                                    d-min.e[e.close]),1)
                                         )
                             }

                             if (!bagel){
                                 out = streduce(c(win[, c()], out[, c()]))
                             }

                             hoodRange = streduce(out, pad)

                             return(self$trim(hoodRange, mod=mod))
                         }
                         else {
                             ## with k, go no more k steps
                             kNeighbors = unique(unlist(ego(private$g, qix, order=k)))
                             return(self$subgraph(kNeighbors, mod=mod)) ## not garanteed size to scale
                         }
                     },

                     dist = function(gr1, gr2 = NULL,
                                     matrix=T,
                                     EPS=1e-9,
                                     include.internal=TRUE, ## consider bp within feature "close"
                                     directed=FALSE, ## if TRUE, only consider gr1-->gr2 paths
                                     verbose=FALSE){
                         if (verbose <- getOption("gGnome.verbose")){
                             message("Given two GRanges, return pairwise shortest path distance.")
                         }

                         if (is.null(gr2)){
                             gr2 = gr1 ## self distance
                         }
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
                         if (any(ix <- strand(gr1)=='*')){
                             strand(gr1)[ix] = '+'
                             gr1 = c(gr1, gr.flipstrand(gr1[ix]))
                         }

                         if (any(ix <- strand(gr2)=='*')){
                             strand(gr2)[ix] = '+'
                             gr2 = c(gr2, gr.flipstrand(gr2[ix]))
                         }

                         ## expand nodes by jabba model to get internal connectivity
                         if (include.internal)
                         {
                             gr1 = gr1[, 'id'] %**% private$segs
                             ## BUG here, why does this overlap command give NA seqnames output?
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

                         adj = self$get.adj()
                         self.l = which(Matrix::diag(adj)>0)

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

                     ## e2j = function(etype="aberrant"){
                     ##     if (verbose <- getOption("gGnome.verbose")){
                     ##         message("Return the junctions based on edges in this graph.")
                     ##     }

                     ##     strmap = setNames(c("+", "-"), c("-", "+"))

                     ##     if (!is.element("type", colnames(private$es)) |
                     ##         !is.element("loose", colnames(private$segs))){
                     ##         tmp = etype(private$segs, private$es, force=T, both=TRUE)
                     ##         private$es = tmp$es
                     ##         private$segs = tmp$segs
                     ##     }

                     ##     es = private$es
                     ##     es[, eix := 1:.N]

                     ##     if (etype=="all"){
                     ##         etype = c("aberrant", "reference", "loose")
                     ##     }

                     ##     abe = es[type %in% etype]
                     ##     if (nrow(abe)==0){
                     ##         empty.out = private$junction = junctions()
                     ##         return(empty.out)
                     ##     }

                     ##     if (any(! c("fromChr", "fromStr", "fromStart", "fromEnd",
                     ##                 "toChr", "toStr", "toStart", "toEnd") %in%
                     ##             colnames(abe))){
                     ##         if (verbose){
                     ##             message("Redo the important metadata gathering.")
                     ##         }

                     ##         abe[, fromStr := ":="(fromChr = as.vector(seqnames(private$segs[from])),
                     ##                               fromStr = as.vector(strand(private$segs[from])),
                     ##                               fromStart = start(private$segs[from]),
                     ##                               fromEnd = end(private$segs[from]),
                     ##                               toChr = as.vector(seqnames(private$segs[to])),
                     ##                               toStr = as.vector(strand(private$segs[to])),
                     ##                               toStart = start(private$segs[to]),
                     ##                               toEnd = end(private$segs[to]))]
                     ##     }

                     ##     if (any(!c("eid", "reid") %in% colnames(abe))){
                     ##         hb = hydrogenBonds(private$segs)
                     ##         hb.map = hb[, c(setNames(from, to),
                     ##                         setNames(to, from))]
                     ##         abe[, ":="(eid = paste(from, to),
                     ##                    reid = paste(hb.map[as.character(to)],
                     ##                                 hb.map[as.character(from)]))]
                     ##     }

                     ##     abe[, ":="(ix = 1:.N,
                     ##                rix = match(reid, eid))]
                     ##     abe[, unique.ix := ifelse(rix>=ix,
                     ##                               paste(ix, rix),
                     ##                               paste(rix, ix))]
                     ##     abe[, eclass := as.numeric(as.factor(unique.ix))]
                     ##     abe[, iix := 1:.N, by=eclass]
                     ##     setkeyv(abe, c("eclass", "iix"))

                     ##     jdt = abe[iix==1, .(eclass, from, to,
                     ##                         fromChr, fromStr, fromStart, fromEnd,
                     ##                         toChr, toStr, toStart, toEnd)]

                     ##     bp1 = dt2gr(jdt[, .(seqnames = fromChr,
                     ##                         strand = strmap[fromStr],
                     ##                         start = ifelse(fromStr=="+", fromEnd, fromStart-1),
                     ##                         end = ifelse(fromStr=="+", fromEnd, fromStart-1),
                     ##                         eclass)])

                     ##     bp2 = dt2gr(jdt[, .(seqnames = toChr,
                     ##                         strand = toStr,
                     ##                         start = ifelse(toStr=="+", toStart-1, toEnd),
                     ##                         end = ifelse(toStr=="+", toStart-1, toEnd),
                     ##                         eclass)])

                     ##     bps = gr.fix(GRangesList(bp1, bp2), private$segs)
                     ##     junc = junctions(grl.pivot(bps))
                     ##     values(junc)$eclass = bp1$eclass
                     ##     values(junc)$type = abe[.(values(junc)$eclass, 1), type]
                     ##     values(junc)$from1 = abe[.(values(junc)$eclass, 1), from]
                     ##     values(junc)$to1 = abe[.(values(junc)$eclass, 1), to]
                     ##     values(junc)$from2 = abe[.(values(junc)$eclass, 2), from]
                     ##     values(junc)$to2 = abe[.(values(junc)$eclass, 2), to]

                     ##     private$junction = junc
                     ##     return(junc)
                     ## },

                     jGraph = function(){
                         ##TODO: migrate the jGraph function here

                     },

                     fillup = function(){
                         "Increase node cn to accomodate edges."
                         A = self$get.adj()
                         ifl = Matrix::colSums(A)
                         ofl = Matrix::rowSums(A)
                         if (!is.element("cn", colnames(values(private$segs)))){
                             private$segs$cn = pmax(ifl, ofl)
                             self$fillin()
                         } else {
                             cns = private$segs$cn
                             private$segs$cn = ifelse(cns<pmax(ifl, ofl),
                                                      pmax(ifl, ofl),
                                                      cns)
                             self$fillin()
                         }
                         if (self$isBalance()){
                             return(self)
                         } else {
                             browser()
                             stop("What happened? We should be always able to fill up the segment cn.")
                         }
                     },

                     make.balance = function(mod=TRUE){
                         "Forcing a gGraph to be junction balanced. Will decrease reference edge cn to honor node cn and aberrant edge cn."
                         if (self$isBalance()){
                             return(self)
                         }

                         adj = self$get.adj()
                         cns = private$segs$cn
                         node.cn = data.table(id = seq_along(private$segs),
                                              cn = private$segs$cn,
                                              ifl = Matrix::colSums(adj),
                                              ofl = Matrix::rowSums(adj),
                                              str = as.character(strand(private$segs)))

                         node.cn[, ":="(ifl.ref =
                                            private$es[type=="reference",
                                                       .(ecn = sum(cn)),
                                                       by = to][
                                                       ,setNames(ecn, to)][as.character(id)],
                                        ifl.abe =
                                            private$es[type=="aberrant",
                                                       .(ecn = sum(cn)),
                                                       by = to][
                                                      , setNames(ecn, to)][as.character(id)],
                                        ofl.ref =
                                            private$es[type=="reference",
                                                       .(ecn = sum(cn)),
                                                       by = from][
                                                      , setNames(ecn, from)][as.character(id)],
                                        ofl.abe =
                                            private$es[type=="aberrant",
                                                       .(ecn = sum(cn)),
                                                       by = from][
                                                      , setNames(ecn, from)][as.character(id)])]


                         node.cn[, ":="(ifl.ref = ifelse(is.na(ifl.ref), 0, ifl.ref),
                                        ifl.abe = ifelse(is.na(ifl.abe), 0, ifl.abe),
                                        ofl.ref = ifelse(is.na(ofl.ref), 0, ofl.ref),
                                        ofl.abe = ifelse(is.na(ofl.abe), 0, ofl.abe))]

                         if (node.cn[, any(ifl.abe>cn | ofl.abe>cn, na.rm=T)]){
                             warning("Some aberrant edges oversubscribe the node copies! Try using JaBbA.")
                             return(self)
                         }

                         ## start reducing ref edge CN to squeeze the aberrant edges in
                         lower.ifl.cn = node.cn[ifl.ref+ifl.abe>cn, .(id, target.ref.cn = cn - ifl.abe - ifl.ref)][order(id)]
                         lower.ofl.cn = node.cn[ofl.ref+ofl.abe>cn, .(id, target.ref.cn = cn - ofl.abe - ofl.ref)][order(id)]

                         private$es[type=="reference" & to %in% lower.ifl.cn[, id] & order(to), cn := cn + lower.ifl.cn[, target.ref.cn]]
                         private$es[type=="reference" & from %in% lower.ofl.cn[, id] & order(from), cn := cn + lower.ofl.cn[, target.ref.cn]]

                         ## fillin the loose ends
                         self$fillin()

                         return(self)
                     },

                     ## property constraints
                     isBalance = function(){
                         "Testing if junction balanced."
                         ## ALERT: this is too loose!!!
                         ## TODO: redo this function!!!
                         if (is.null(private$es)){
                             return(NULL)
                         } else if (!"cn" %in% c(colnames(values(private$segs)), colnames(private$es))) {
                             return(FALSE)
                         } else if (length(private$segs)==0){
                             return(TRUE)
                         } else if (nrow(private$es)==0){
                             return(TRUE)
                         }

                         if (!is.element("cn", colnames(values(private$segs))) |
                             !is.element("cn", colnames(private$es))){
                             return(FALSE)
                         }
                         ## DONE: use adj to calc if every segment is balanced on both sides
                         adj = self$get.adj()
                         cns = private$segs$cn

                         node.cn = data.table(id = seq_along(private$segs),
                                              cn = private$segs$cn,
                                              ifl = Matrix::colSums(adj),
                                              ofl = Matrix::rowSums(adj),
                                              str = as.character(strand(private$segs)))

                         whichTerminal = node.cn[,which(ifl==0 | ofl==0)]
                         whichNa = which(is.na(private$segs$cn))
                         validTerminal = setdiff(whichTerminal, whichNa)

                         bal = node.cn[
                           , ifelse((1:.N) %in% validTerminal,
                                    cns == pmax(ifl, ofl) | ((ifl==0) & (ofl==0)),
                                    (cns == ifl) & (cns==ofl))]

                         bal = bal[which(!is.na(cns))]
                         return(all(bal))
                     },

                     get.loose = function(){
                         "Return all loose ends as a GRanges."
                         if (!is.element("loose", colnames(values(private$segs)))){
                             private$segs = etype(segs = private$segs,
                                                  es = private$es,
                                                  force=TRUE, both=TRUE)$segs
                         }
                         return(private$segs %Q% (loose==TRUE))
                     },

                     get.walk = function(v = numeric(0),
                                         e = NULL,
                                         peel = FALSE,
                                         cn = NULL){
                         "Generate gWalks object given by node or edge sequences."
                         ## TODO: if given j, override v
                         ## MOMENT
                         ## test validity (existence and CN) of path defined by e
                         ## if passed, convert to v and recurse
                         ## How to test if a vetices sequence is a valid path on a graph???
                         if (length(v)>0 & is.null(e)){
                             e = as.data.table(cbind(shift(v), v))[-1, , drop=FALSE]
                             colnames(e) = c("from", "to")
                         } else if (!is.null(e)){
                             if (inherits(e, "data.frame")){
                                 if (!all(c("from", "to") %in% colnames(e))){
                                     message("e is provided but no 'from' and 'to' columns.")
                                     return(NULL)
                                 }
                                 if (nrow(e)==0){
                                     message("e is empty")
                                 }
                                 v = c(e$from, e$to[nrow(e)])
                                 e = data.table(e)
                             } else {
                                 message("e needs to be a data.frame-like object.")
                                 return(NULL)
                             }
                         } else {
                             message("No input path given.")
                             return(NULL)
                         }

                         ## the edge list of this path
                         es = private$es
                         es[, eid := paste(from, to)]
                         e[, eid := paste(from, to)]
                         if (!all(e[, eid] %in% es[, eid])){
                             message("Given path is not a valid path in the graph!!")
                             return(NULL)
                         }

                         ## MOMENT
                         ## TODO: if peel==TRUE, deduce this walk from the object
                         wk = gWalks$new(grl = private$segs[v], cn = cn)
                         return(wk)
                     },

                     random.walk = function(start,
                                            steps,
                                            mode = c("out", "in", "all"),
                                            stuck = c("return", "error")){
                         "Generate large numbers of gWalks using the algorithm DeepWalk."
                         ig = self$get.g()
                         v = tryCatch(igraph::random_walk(ig,
                                                 start = start,
                                                 steps = steps,
                                                 mode = mode,
                                                 stuck = stuck),
                                      error = function(e){
                                          warning("Can't generate the desired walks")
                                          return(NULL)
                                      })
                         if (is.null(v)){
                             return(NULL)
                         }
                         return(self$get.walk(v = as.numeric(v)))
                     },


                     ## all t
                     chromoplexy = function(pad = 1e3){
                         "Identifying parts of the graph that are probably produced from chromoplexy events. In gGraph class the method ignores information from CN."
                         if (is.null(private$g)){
                             G = self$get.g()
                         } else {
                             G = private$g
                         }

                         if (is.null(private$junction)){
                             juncs = e2j(private$segs, private$es)
                         } else {
                             juncs = private$junction
                         }
                         ## browser()
                         ## MOMENT
                         ab.edges = data.table(data.frame(values(juncs)))
                         ## Marcin's take
                         if (!is.null(jab))
                         {
                             if (filt.jab)
                             {
                                 nnab = which(rowSums(is.na(rbind(jab$ab.edges[, 1:2, 1])))==0)
                                 edge.ix = which(jab$adj[rbind(jab$ab.edges[nnab, 1:2, 1])]>0)
                                 jab$ab.edges = jab$ab.edges[edge.ix, ,,drop = F]
                             }
                             else
                                 edge.ix = 1:nrow(kag$ab.edges)
                             kag = jab
                             sol = jab
                         }
                         else
                             edge.ix = 1:nrow(kag$ab.edges)

                         G = kag$G

                         if (is.null(kag$tile))
                             kag$tile = kag$segstats

                         nnab = which(rowSums(is.na(rbind(kag$ab.edges[, 1:2, 1])))==0)
                         if (ref.only)
                         {
                             adj2 = kag$adj
                             adj2[kag$ab.edges[nnab, 1:2, 1]] = 0
                             adj2[kag$ab.edges[nnab, 1:2, 2]] = 0
                             G = graph(as.numeric(t(Matrix::which(adj2!=0, arr.ind = T))), n = length(kag$segstats), directed = T)
                         }

                         ## define edge source to edge sink distance
                         ## this is minimum between (1) sum of vertex width of path from e2 source to e1 sink (excluding source and sink)
                         ## and (2) sum of vertex width of path from e1 sink to e2 source (including source and sink)

                         tmp = igraph::get.edges(G, E(G))
                         E(G)$from = tmp[,1]
                         E(G)$to = tmp[,2]
                         E(G)$weights.source = width(kag$tile[E(G)$from])

                         ab.edges = cbind(rbind(kag$ab.edges[nnab, c('from', 'to'), '+'], kag$ab.edges[nnab, c('from', 'to'), '-']), junc.id = rep(nnab,2))

                         if (nrow(ab.edges)==0)
                             return(list(cycles = NULL, paths = NULL))

                                        #    emap = c(1:nrow(kag$ab.edges), -(1:nrow(kag$ab.edges)))
                         emap = c(nnab, -nnab)

                         ## basically constructing the jGraph
                         jg = self$jgraph()
                         D1 = D2 = array(Inf, dim = rep(nrow(ab.edges),2))

                         uix = unique(c(ab.edges[,1], ab.edges[,2]))
                         uixmap1 = match(ab.edges[,1], uix)
                         uixmap2 = match(ab.edges[,2], uix)

                         tmp = shortest.paths(G, uix, uix, weights = E(G)$weights.source, mode = 'out')


                         ## deletion bridge, or reciprocal
                         if (reciprocal)
                         {
                             ## "deletion bridge", i.e. source to sink bridge
                             D1 = t(sweep(tmp[uixmap1, uixmap2],
                                          1, width(kag$tile[ab.edges[,1]]))) ## subtract width of first vertex from path length (second vertex already excluded)
                             D1[do.call('rbind', lapply(ab.edges[,2], function(x) ab.edges[,1] %in% x))] = NA ## edge case where e1 sink = e2 source
                         }

                         ## "amplification bridge", i.e. sink to source bridge
                         if (hijacked)
                         {
                             D2 = sweep(tmp[uixmap2, uixmap1],
                                        2, -width(kag$tile[ab.edges[,1]])) ## add width of last vertex to path (first vertex already included)
                         }


                         D = matrix(pmin(D1, D2, na.rm = T), nrow = nrow(D1), ncol = nrow(D2))
                         D.which = matrix(ifelse(D1<D2, 1, 2), nrow = nrow(D1), ncol = nrow(D2))
                         D.which[is.na(D.which)] = 2
                         D[is.infinite(D)] = NA
                         D[cbind(1:nrow(D), 1:nrow(D))] = NA

                         ## quasi pairvvs are ab edge pairs within a certain distance of each other on the graph
                         quasi.pairs = which(D<dist, arr.ind = T)
                         quasi.pairs.which = D.which[quasi.pairs]

                         ## now need to check .. depending on whether edge pair is deletion bridge or amp bridge or fully reciprocal
                         ## whether associated vertices show a copy change "in the right direction"

                         ## to do this, we need to examine the vertices "in between" for a deletion bridge and the source / sink vertices
                         ## in an amplification bridge, and see if they show a copy change with respect their reference parents

                         ## for reciprocal pairs, the source and sink will be the same

                         adj.ref = kag$adj; adj.ref[ab.edges[, 1:2]] = 0

                         del.bridge.candidate = which(quasi.pairs.which == 1)
                         v1 = ab.edges[quasi.pairs[del.bridge.candidate, 1], 2]
                         v1.parent = apply(adj.ref[, v1, drop = FALSE], 2, function(x) which(x != 0)[1])
                         v1.child = apply(adj.ref[v1, , drop = FALSE], 1, function(x) which(x != 0)[1])
                         v2 = ab.edges[quasi.pairs[del.bridge.candidate, 2], 1]
                         v2.child = apply(adj.ref[v2, , drop = FALSE], 1, function(x) which(x != 0)[1])
                         v2.parent = apply(adj.ref[, v2, drop = FALSE], 2, function(x) which(x != 0)[1])

                         recip = del.bridge.candidate[which(v2.child == v1)]
                         nonrecip = which(v2.child != v1) ## these need to fulfill the "deletion bridge criterion"

                                        # test for deletion bridge criterion, i.e. does v1 have greater copy number than its
                                        # child, and does v2 have greater copy number than its parent?
                         if (!is.null(sol))
                             del.bridge = del.bridge.candidate[nonrecip[(sol$segstats$cn[v2.child[nonrecip]] < sol$segstats$cn[v2[nonrecip]] & sol$segstats$cn[v1.parent[nonrecip]] < sol$segstats$cn[v1[nonrecip]]) | D[quasi.pairs][del.bridge.candidate[nonrecip]] < cn.dist]]
                         else
                             del.bridge = del.bridge.candidate

                         amp.bridge.candidate = which(quasi.pairs.which == 2)
                         v1 = ab.edges[quasi.pairs[amp.bridge.candidate, 1], 2]
                         v1.parent = apply(adj.ref[, v1, drop = FALSE], 2, function(x) which(x != 0)[1])
                         v2 = ab.edges[quasi.pairs[amp.bridge.candidate, 2], 1]
                         v2.child = apply(adj.ref[v2,,drop = FALSE], 1, function(x) which(x != 0)[1])

                                        # test for amp bridge criterion, does v1 have higher copy number than its parent, does v2 have higher copy number than its child?

                         if (!is.null(sol))
                             amp.bridge = amp.bridge.candidate[(sol$segstats$cn[v1.parent] < sol$segstats$cn[v1] & sol$segstats$cn[v2.child] < sol$segstats$cn[v2])
                                                               | D[quasi.pairs][amp.bridge.candidate] < cn.dist]
                         else
                             amp.bridge = amp.bridge.candidate

                         ## now put together all surviving edges into a graph and try to find cycles

                         ## store data frame of edge pairs for bp graph
                         ## NOTE: every node in bp graph is an edge in the original karyograph, and thus edges in the bp graph represent ordered <edge pairs>
                         bp.df = data.frame(
                             e1 = quasi.pairs[c(recip, del.bridge, amp.bridge), 1], e2 = quasi.pairs[c(recip, del.bridge, amp.bridge), 2],
                             from = ab.edges[quasi.pairs[c(recip, del.bridge, amp.bridge), 1], 1],
                             to = ab.edges[quasi.pairs[c(recip, del.bridge, amp.bridge), 2], 2],
                             type = c(rep('recip', length(recip)), rep('del', length(del.bridge)), rep('amp', length(amp.bridge))), stringsAsFactors = F)

                         bp.df = bp.df[!is.na(bp.df$e1) & !is.na(bp.df$e2), ]

                         ## make adj matrix of breakpoints, basically by matching bp1 and bp2 if "to" field of bp1 = "from" field of bp2
                         ## here we are looking for <exact> matches because we are now going to join an edge to another edge if the target
                         ## of one edge is the source of the next
                                        #    adj.bp = matrix(0, nrow = nrow(bp.df), ncol = nrow(bp.df))
                                        #    for (i in 1:ncol(adj.bp))
                                        #      adj.bp[i,] = bp.df$from %in% bp.df$to[i] & !is.na(bp.df$to[i])

                         ## breakpoint graph links every edge to every other edge via "quasi pair" connection
                         ## we find cycles and paths in this graph
                         adj.bp = sparseMatrix(i = bp.df$e1, j = bp.df$e2, x = 1, dims = rep(nrow(ab.edges), 2))

                         if (junc.only){
                             G.bp = igraph::graph_from_adjacency_matrix(adj.bp)
                             comp = components(G.bp, "strong")
                             good.comp = which(comp$csize>1)
                             good.ix = which(comp$membership %in% good.comp)
                             return(unique(ab.edges[good.ix, 3]))
                         } else {
                             if (verbose)
                                 if (paths)
                                     cat(sprintf('Running with paths on breakpoint graph with dim %s vertices and %s edges\n', nrow(adj.bp), sum(adj.bp)))
                                 else
                                     cat(sprintf('Running without paths on breakpoint graph with dim %s vertices and %s edges\n', nrow(adj.bp), sum(adj.bp)))

                             if (prod(dim(adj.bp))>0)
                             {
                                 ## want to exclude any paths involving breaks and their pairs
                                        #        tmp = split(1:nrow(ab.edges), ab.edges[,'junc.id'])
                                        #        exclude.ij = cbind(ab.edges[,3], unlist(tmp))
                                        #        exclude = sparseMatrix(exclude.ij[,1], exclude.ij[,2], x = 1)
                                 exclude = NULL
                                 dt = data.table(i = 1:nrow(adj.bp), j = mmatch(adj.bp, adj.bp[!duplicated(as.matrix(adj.bp)), , drop = F]), key = 'j')
                                 dtu = dt[!duplicated(j), ]
                                 pc = all.paths(adj.bp[dtu$i, dtu$i, drop = FALSE], all = paths, verbose = verbose, interval = interval, chunksize = chunksize, exclude = exclude)
                                 if (length(pc$cycles)>0)
                                     pc$cycles = lapply(pc$cycles, function(x) dtu[x, ]$i)
                                 if (length(pc$paths)>0)
                                     pc$paths = lapply(pc$paths, function(x) dtu[x, ]$i)
                             }
                             else
                                 return(list(paths = c(), cycles = c()))

                             ## if there are other possible "bridge links" between members of a cycle that do not involve
                             ## members of the cycle.  Fix: Best way to fix this would be actually recompute shortest paths after removing
                             ## edges cresponding to edges in the path.
                             .check.pc = function(x, is.cycle = F)
                             {
                                 if (is.cycle)
                                     tmp.edges = cbind(x, c(x[-1], x[1]))
                                 else
                                     tmp.edges = cbind(x[-length(x)], x[-1])
                                 tmp.D.which = D.which[tmp.edges]  ## D.which keeps track of whether we linked these edges via D1 or D2
                                 if (any(ix <- tmp.D.which==1)) ## if 1 then we are looking for path from col 2 to col 1, so flip
                                     tmp.edges[ix,] = tmp.edges[ix, c(2:1)]
                                 tmp.ab.edges = cbind(ab.edges[cbind(tmp.edges[,1], tmp.D.which)], ab.edges[cbind(tmp.edges[,2], ifelse(tmp.D.which == 1, 2, 1))])
                                 tmp.sp = lapply(1:nrow(tmp.ab.edges), function(i)
                                     get.shortest.paths(G, tmp.ab.edges[i,1], tmp.ab.edges[i,2], weights = E(G)$weights.source, mode = 'out')$vpath[[1]])
                                 if (any(ix <- tmp.D.which==1))
                                     tmp.sp[ix] = lapply(tmp.sp[ix], function(x) x[-c(1, length(x))])
                                 bp.id = unique(unlist(lapply(tmp.sp, function(x) E(G, path = x)$bp.id)))
                                 return(any(x %in% bp.id))
                             }


                             ## xtYao modified: mclapply to replace lapply and sapply
                             if (length(pc$cycles)>0)
                             {
                                 pc$cycles = pc$cycles[!unlist(mclapply(pc$cycles, .check.pc, is.cycle = T, mc.cores=mc.cores))]
                                 pc$cycles = mclapply(pc$cycles, function(x) sign(emap[x])*edge.ix[abs(emap[x])], mc.cores=mc.cores)
                                 pc$cycles = pc$cycles[!duplicated(unlist(mclapply(pc$cycles, function(x) paste(unique(sort(x)), collapse = ' '), mc.cores=mc.cores)))]
                                 pc$cycles = pc$cycles[order(-unlist(mclapply(pc$cycles, length, mc.cores = mc.cores)))]
                             }

                             if (length(pc$paths)>0)
                             {
                                 pc$paths = pc$paths[!unlist(mclapply(pc$paths, .check.pc, is.cycle = F, mc.cores = mc.cores))]
                                 pc$paths = mclapply(pc$paths, function(x) sign(emap[x])*edge.ix[abs(emap[x])], mc.cores = mc.cores)
                                 pc$paths = pc$paths[!duplicated(unlist(mclapply(pc$paths, function(x) paste(unique(sort(x)), collapse = ' '), mc.cores = mc.cores)))]
                                 pc$paths = pc$paths[order(-unlist(mclapply(pc$paths, length, mc.cores = mc.cores)))]
                             }

                             return(pc)
                      }
                     },

                     chromothripsis = function(){
                         "Identifying parts of the graph that are probably produced from "
                     },
                     kid.frag = function(){
                         "Putative TIC events."
                     },
                     bfb = function(){
                         "Classic breakage-fusion-bridge cycles."

                     }
                 ),

                 private = list(
                     ## ----- private fields
                     ## ===== required
                     ## node/vertex, a GRanges obj of strand-specific ranges
                     segs = NULL,
                     ## data.table of all edges in g, from, to, cn, type
                     ## type can be ref, aberrant, loose
                     es = NULL,

                     ## ===== optional slots
                     ## ALERT: whenever segs or es changes, all of these need to be reset!
                     ## igraph obj representing the graph structure
                     g = NULL,
                     ## temporary segs for backtrace when modifying segs
                     tmpSegs = NULL,
                     ## putative junctions, junctions
                     junction = NULL,
                     ## ploidy is set to 2, only to init null graph,
                     ## otherwise inferred from segs
                     .ploidy = NULL,
                     ## tumor cell proportion
                     .purity = NULL,
                     ## the partition result of 'g'
                     partition = NULL,

                     ## ----- private methods
                     ## reset optional fields
                     reset = function(){
                         private$g = NULL
                         private$tmpSegs = NULL
                         private$junction = NULL
                         private$.ploidy = NULL
                         private$.purity = NULL
                         private$partition = NULL
                     },

                     ## break the current segments into new segments
                     makeSegs = function(bps){
                         ## DONE: once finished, move to private methods
                         private$tmpSegs = private$segs
                         names(bps) = NULL
                         private$segs = gUtils::gr.breaks(bps, private$segs)
                         return(self)
                     },

                     ## initialize by directly giving fields values
                     gGraphFromScratch = function(segs,
                                                  es,
                                                  junc=NULL,
                                                  ploidy=NULL,
                                                  purity=NULL){
                         if (getOption("gGnome.verbose")){
                             message("Nodes as GRanges, edges as data.frame or adj matrix.")
                         }

                         if (!is.null(names(segs))){
                             if (any(duplicated(names(segs)))){
                                 names(segs) = NULL
                             }
                         }

                         if (length(segs)==0){
                             private$segs = GRanges()
                             return(self)
                         }

                         if (is.null(es)){
                             es = data.table(from = numeric(0),
                                             to = numeric(0),
                                             type = character(0))
                         }
                         if (!inherits(es, "data.frame")){
                             es = data.table(from = numeric(0),
                                             to = numeric(0),
                                             type = character(0))
                         } else {
                             es = data.table(es)
                         }

                         if (any(dim(es))==0){
                             private$es = es
                         }

                         ## check es input
                         if (inherits(es, "data.frame")){
                             if (!all(c("from", "to") %in% colnames(es))){
                                 stop("Error: Given edge data must have 'from' and 'to' fields.")
                             }

                             es = es[from %in% seq_along(segs) & to %in% seq_along(segs)]
                         }
                         else if (inherits(es, "matrix") | inherits(es, "Matrix")){
                             A = es
                             if (!all(dim(A)==length(segs))){
                                 stop("Given adjacency matrix in wrong dimension.")
                             }

                             ## we allow numeric or logical values
                             if (inherits(A[1,1], "numeric")){
                                 ## when it's numeric & non-negative, the value is copy number
                                 if (any(A<0)){
                                     A = A>0
                                 }
                                 else {
                                     es = which(A>0, arr.ind=T)
                                     colnames(es) = c("from", "to")
                                     cn = A[es]
                                     es = as.data.table(es)
                                     es[, cn := cn]
                                 }
                             }

                             if (inherits(A[1,1], "logical")){
                                 es = which(A==TRUE, arr.ind=T)
                                 colnames(es) = c("from", "to")
                                 es = as.data.table(es)
                             }
                         }

                         ## ALERT: if no "loose" col in segs, deduct from edges
                         ## OR if no "type" col in es, do it too
                         if (!is.element("loose",colnames(values(segs))) |
                             !is.element("type", colnames(es))){
                             tmp = etype(segs, es, both=TRUE, force=TRUE)
                             segs = tmp$segs
                             es = tmp$es
                         }

                         ## ALERT: from now on the reference genome
                         ## is completely defined by the "seqinfo" of nodes
                         ## here is to overwrite
                         private$segs = segs
                         ## MOMENT
                         private$segs$tile.id = get.tile.id(private$segs)


                         hb = hydrogenBonds(private$segs)
                         hb.map = hb[, setNames(to, from)]
                         ## finished converting adj matrix to data.table
                         ## now check if it is skew-symmetric
                         es[, eid := paste(from, to)]
                         es[, reid := paste(hb.map[as.character(to)],
                                            hb.map[as.character(from)])]
                         ematch = es[, match(eid, reid)]

                         es[, ":="(ix = 1:.N,
                                   rix = match(reid, eid))]
                         es[, unique.ix := ifelse(rix>=ix,
                                                  paste(ix, rix),
                                                  paste(rix, ix))]
                         es[, eclass := as.numeric(as.factor(unique.ix))]
                         es[, iix := 1:.N, by=eclass]

                         ## ALERT: sometimes there are NAs in ematch!!!
                         ## TODO: how to deal with NA in ematch???
                         if ("cn" %in% colnames(es)){
                             if (all(es[ematch, cn]==es[, cn], na.rm=T)){
                                 if (as.logical(getOption("gGnome.verbose"))){
                                     message("Edge copies balanced!")
                                 }
                             } else {
                                 ## TODO: maybe don't try to do too much????
                                 ## or help the user with this????
                                 stop("Error: Given edge data is not skew-symmetric!!!")
                             }
                         }

                         private$es = es
                         ## relabel the terminals!
                         whichTerminal = private$es[, setdiff(seq_along(private$segs), intersect(from, to))]
                         mcols(private$segs)$terminal = seq_along(private$segs) %in% whichTerminal

                         ## private$segs$terminal = seq_along(private$segs) %in% whichTerminal
                         private$g = igraph::make_directed_graph(t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))
                         private$.ploidy = ploidy
                         private$.purity = purity
                         return(self)
                     }
                 ),

                 active = list(
                     ## ======= getters
                     segstats = function(){
                         return(private$segs)
                     },
                     edges = function(){
                         if (!"type" %in% colnames(private$es)){
                             private$es = etype(private$segs, private$es, force=T)
                         }
                         return(private$es)
                     },
                     junctions = function(){
                         if (is.null(private$junction)){
                             private$junction = e2j(private$segs, private$es)
                         }
                         return(copy(private$junction))
                     },
                     G = function(){
                         if (is.null(private$g)){
                             self$get.g()
                         }
                         return(private$g)
                     },
                     adj = function(){
                         return(self$get.adj())
                     },
                     A = function(){
                         return(self$get.adj())
                     },
                     parts = function(){
                         if (length(private$segs)==0 | is.null(private$es)){
                             return(NULL)
                         } else if (is.null(private$partition)){
                             tmp = self$components()
                         }
                         return(private$partition)
                     },
                     seqinfo = function(){
                         return(seqinfo(private$segs))
                     },
                     purity = function(){
                         return(private$.purity)
                     },
                     ploidy = function(){
                         if (is.null(private$.ploidy)){
                             private$.ploidy = get.ploidy(private$segs)
                         }
                         return(private$.ploidy)
                     },

                     ## ========== viz
                     td = function(){
                         return(self$gg2td())
                     },
                     win = function(){
                         self$window()
                     },
                     ig = function(){
                         ## DONE: make igraph plot
                         return(self$layout())
                     }
                 )
                 )

#' @export
#' @export bGraph
bGraph = setClass("bGraph")

##############################
## bGraph
##############################
#' @name bGraph-class
#' @title junction-balanced graph
#' @docType class
#' @description
#' Descendant of gGraph class, where junction balance restraint must be met at all times.
#'
#' @import R6
#' @import Matrix
#'
#'
#' @exportClass bGraph
#' @export bGraph
##############################
bGraph = R6::R6Class("bGraph",
                     inherit = gGraph,
                     public = list(
                         ## overwrite constructor: restrict about junction balance
                         initialize = function(gG=NULL, jabba=NULL, prego=NULL){
                             if (!is.null(gG)){
                                 if (inherits(gG, "gGraph")){
                                     ## TODO
                                     private$gGraphFromScratch(gG$segstats, gG$edges, gG$junctions, gG$ploidy, gG$purity)
                                     return(self)
                                 } else {
                                     stop("Invalid input gG.")
                                 }
                             } else if (!is.null(jabba)) {
                                 if (is.character(jabba))
                                 {
                                     if (file.exists(jabba))
                                         jabba = readRDS(jabba)
                                     else
                                         stop(paste('file', jabba, 'not found'))
                                 }

                                 self$jab2gg(jabba=jabba)
                                 if (self$isBalance()){
                                     return(self)
                                 } else {
                                     stop("Invalid input gG.")
                                 }
                             } else if (!is.null(prego)){
                                 if (is.character(prego))
                                 {
                                     if (file.exists(prego))
                                         self$pr2gg(fn=prego)
                                     else
                                         stop(paste('file', jabba, 'not found'))
                                 }

                                 if (self$isBalance()){
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
                             ## cat('Based on reference genome: ')
                             ## cat(private$segs)
                             cat('\n\n')
                             cat('Total non-loose segmentation:')
                             if ("loose" %in% colnames(values(private$segs))){
                                 cat(length(private$segs %Q% (loose==F & strand=="+")))
                             } else {
                                 ## ALERT!!! TODO!!! This means we have to make sure if there is
                                 ## loose end, it must be labeled in the nodes.
                                 cat(length(private$segs %Q% (strand=="+")))
                             }
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
                                 out = as(out, "bGraph")
                                 return(out)
                             }
                         },

                         ## decompose graph into all possible haplotypes
                         ## MOMENT
                         walk = function(outdir="tmp.walk",
                                         max.iteration = Inf,
                                         mc.cores = 1,
                                         verbose = T,
                                         all.paths = FALSE,
                                         nsolutions = 100,
                                         tilim = 100,
                                         gurobi = FALSE,
                                         cplex = !gurobi){
                             "Enumerate all the possible multiset of walks or give the most parsimonious ones that can be represented by this graph."
                             ## ASSUMPTION: no duplicated rows in $segs
                             ## TODO: something's wrong here, need redo
                             if (verbose <- getOption("gGnome.verbose")){
                                 message("Enumerating all possible karyotypes with minimum cardinal numbers.")
                             }

                             if (length(private$segs)==0 | is.null(private$es)){
                                 if (verbose){warning("Empty graph. Nothing to decompose.")}
                                 return(NULL)
                             }

                             if (!"type" %in% colnames(private$es) |
                                 !"loose" %in% colnames(values(private$segs))){
                                 et = etype(private$segs, private$es, T, T)
                                 private$segs = et$segs
                                 private$es = et$es
                             }

                             which.loose = which(private$segs$loose==TRUE)
                             if (length(which.loose)>0){
                                 segs = private$segs[-which.loose]
                             } else {
                                 segs = private$segs
                             }
                             new.ix = setNames(seq_along(segs),
                                               setdiff(seq_along(private$segs),
                                                       which.loose))
                             ed0 =
                                 private$es[type!="loose"][, ":="(from = new.ix[as.character(from)],
                                                                  to = new.ix[as.character(to)])]

                             A = self$get.adj()[setdiff(seq_along(segs), which.loose),
                                                setdiff(seq_along(segs), which.loose)]

                             ## node mapping
                             hb = hydrogenBonds(segs)
                             hb.map = hb[, setNames(from, to)]

                             ## keep the network flowing: adding sources and sinks
                             ifl = Matrix::colSums(A)
                             ofl = Matrix::rowSums(A)
                             avail = ifl - ofl

                             if (sum(avail)!=0){
                                 stop("The edge copy flow is not balanced!")
                             }

                             slack = rbind(data.table(from = NA,
                                                      to = which(avail<0),
                                                      cn = abs(avail[which(avail<0)]),
                                                      type = "slack.in"),
                                           data.table(from = which(avail>0),
                                                      to = NA,
                                                      cn = abs(avail[which(avail>0)]),
                                                      type = "slack.out"))

                             ed0 = rbind(ed0[,.(from, to, cn, type)],
                                         slack)
                             ed0[, ":="(eid = paste(from, to),
                                        reid = paste(hb.map[as.character(to)],
                                                     hb.map[as.character(from)]))]

                             ## get eclass
                             ed0[, ":="(ix = 1:.N,
                                        rix = match(reid, eid))]
                             ed0[, unique.ix := ifelse(rix>=ix,
                                                       paste(ix, rix),
                                                       paste(rix, ix))]
                             ed0[, eclass := as.numeric(as.factor(unique.ix))]
                             ed0[, iix := 1:.N, by=eclass]
                             ## rename non-slack edge types
                             ed0[!grepl("slack", type), type := "nonslack"]

                             ## get incidence matrix
                             ## vertices x edges
                             ## TODO: assemble h, input to karyoMIP -- e, e.ij, B, eclass, etype
                             ## ASSUMPTION: private$segs is sorted by loose then strand
                             ## copies going away
                             ii1 = c(ed0[type=="nonslack", from],
                                     ed0[type=="slack.out", from])
                             jj1 = c(ed0[, which(type=="nonslack")],
                                     ed0[, which(type=="slack.out")])
                             xx1 = c(rep(-1, ed0[,sum(type=="nonslack")]),
                                     rep(-1, ed0[,sum(type=="slack.in")]))
                             ## copies coming in
                             ii2 = c(ed0[type=="nonslack", to],
                                     ed0[type=="slack.in", to])
                             jj2 = c(ed0[, which(type=="nonslack")],
                                     ed0[, which(type=="slack.in")])
                             xx2 = c(rep(1, ed0[,sum(type=="nonslack")]),
                                     rep(1, ed0[,sum(type=="slack.in")]))

                             ##
                             B = sparseMatrix(i = ii1,
                                              j = jj1,
                                              x = xx1,
                                              dims=c(length(segs), nrow(ed0))) +
                                 sparseMatrix(i = ii2,
                                              j = jj2,
                                              x = xx2,
                                              dims=c(length(segs), nrow(ed0)))

                             ## form the hypothesis list
                             h = list(e = ed0[, cn],
                                      e.ij = as.matrix(ed0[, .(from, to)]),
                                      B = B,
                                      eclass = ed0[, eclass],
                                      etype = ed0[, ifelse(grepl("nonslack", type),
                                                           "nonslack", "slack")])

                             ## compute convex basis of B
                             K = convex.basis(B)
                             prior = rep(1, ncol(K))

                             if (getOption("gGnome.debug")){
                                 saveRDS(h, "h.rds")
                                 saveRDS(K, "K.rds")
                                 saveRDS(prior, "prior.rds")
                             }

                             ## MOMENT
                             ## browser()
                             if (all.paths)
                             {
                                 ## outfile.allpaths.pdf = sprintf('%s/%s.allpaths.pdf', outdir, label)

                                 ## if (verbose){
                                 ##     cat('Generating all walks\n')
                                 ## }

                                 ## repurpose karyoMIP.to.path to generate all paths
                                 ## using "fake solution" i.e. all 1 weights, to karyoMIP as input
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
                                             if (length(pal.wix)>0)C
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
                                         cl = igraph::clusters(G, 'weak')$membership ## clusters based on this adjacency relationship
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
                                 out$gtrack.allpaths = c(
                                     gt.walk,
                                     td.seg,
                                     td.rg)
                                 pdf(outfile.allpaths.pdf, height = 30, width = 24)
                                 gTrack::plot(out$gtrack.allpaths,
                                              windows = win, links = kag$junctions)
                                 dev.off()
                                 out$README = paste(out$README, 'allpaths= all paths through windows (not just optimal ones), td.allpaths = gTrack object of plot of all paths')
                             }

                             ## MEAT
                             ## TODO: convert karyoMIP solution to gWalks object,
                             ## with the new gw definition
                             ## is.cyc = Matrix::colSums(K[h$etype == 'slack', ])==0 &
                             ## Matrix::colSums((Bc %*% K)!=0)==0
                             karyo.sol = karyoMIP(K, h$e, h$eclass,
                                                  nsolutions = nsolutions,
                                                  tilim = tilim,
                                                  cpenalty = 1/prior,
                                                  gurobi = gurobi)

                             ## if (saveAll){
                             ##     saveRDS(karyo.sol, "temp.walk/allSol.rds")
                             ## }
                             if (cplex){
                                 kag.sol = karyo.sol[[1]]
                             } else {
                                 kag.sol = karyo.sol
                             }

                             p = karyoMIP.to.path(kag.sol, K, h$e.ij, segs)
                             p$paths = mclapply(p$paths, as.numeric, mc.cores=mc.cores)

                             ## construct gWalks as result
                             gw = gWalks$new(segs=segs,
                                             paths=p$paths,
                                             is.cycle=p$is.cyc,
                                             cn = p$cn)
                             return(gw)
                         },

                         ## TODO: hurestic walk decomposition
                         ## new idea: if we assign weight
                         ## BUG: suddenly it doesn't work!!!!!!!!!?????????????
                         walk2 = function(verbose = FALSE,
                                          grl = TRUE,
                                          e.weight = NULL,
                                          gurobi = FALSE,
                                          cplex = !gurobi){
                             "Heuristic for decomposing a junction-balanced graph into a multiset of walks."
                             if (length(private$segs)==0 | is.null(private$es)){
                                 if(verbose){warning("Empty graph. Nothing to decompose.")}
                                 return(NULL)
                             }
                             segs = private$segs
                             sl = seqlengths(segs)
                             hb = hydrogenBonds(segs)
                             hb.map = hb[, setNames(from, to)]
                             cn.adj = self$get.adj()
                             adj = as.matrix(cn.adj)
                             adj.new = adj*0
                             ## ALERT!!! see below
                             adj[which(adj!=0, arr.ind = TRUE)] =
                                 width(segs)[which(adj!=0, arr.ind = TRUE)[,2]]
                             ## make all edges a large number by default

                             if (verbose){
                                 message('Setting edge weights to destination widths for reference edges and 1 for aberrant edges')
                             }

                             ab.edges = private$es[type=="aberrant", cbind(from, to)]
                             ## ALERT!!!
                             if (nrow(ab.edges)>0)
                             {
                                 adj[ab.edges] = sign(cn.adj[ab.edges]) ## make ab.edges = 1
                             }
                             adj[is.na(adj)] = 0
                             cn.adj[which(is.na(cn.adj))] = 0

                             ## ALERT!!! major change
                             G = graph.adjacency(adj, weighted = 'weight')

                             ## define ends not using degree (old method) but using
                             ## either telomeres or loose ends
                             ## (otherwise lots of fake ends at homozygous deleted segments)
                             ss = gr2dt(segs)[ , vid:= 1:length(seqnames)]
                             ss[loose == TRUE, is.end := TRUE]
                             ss[loose == FALSE,
                                is.end := 1:length(loose) %in%
                                    c(which.min(start), which.max(end)),
                                by = list(seqnames, strand)]
                             ends = which(ss$is.end)
                             src = (Matrix::colSums(adj)[ends]==0) ## indicate which are sources

                             ## sanity check
                             unb = which(!ss$is.end & Matrix::rowSums(self$get.adj(), na.rm = TRUE) !=
                                         Matrix::colSums(self$get.adj(), na.rm = TRUE))

                             if (length(unb)>0)
                             {
                                 message(sprintf('JaBbA model not junction balanced at %s non-ends! Adding these to "ends"', length(unb)))
                                 ends = c(ends, unb)         ## shameless HACK ... TOFIX
                             }

                             i = 0
                             ## adjust weight just before creating D
                             ## assign lighter weight to higher copy
                             ## D records distance from ends to every node
                             D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

                             ## sort shortest paths
                             ij = as.data.table(
                                 which(!is.infinite(D), arr.ind = TRUE)
                             )[, dist := D[cbind(row, col)]][
                                 row != col, ][order(dist), ][
                               , row := ends[row]][
                               , col := ends[col]]

                             maxrow = length(ends)*max(cn.adj[ends, ends], na.rm = TRUE)
                             vpaths = rep(list(NA), maxrow)
                             epaths = rep(list(NA), maxrow)
                             cns = rep(NA, maxrow)
                             palindromic.path = rep(FALSE, maxrow)
                             palindromic.cycle = rep(FALSE, maxrow)

                             nb.all = which(Matrix::rowSums(cn.adj) != Matrix::colSums(cn.adj))
                             cn.adj0 = cn.adj
                             G0 = G
                             D0 = D

                             #' first peel off "simple" paths i.e. zero degree
                             #' ends with >0 copy number
                             psimp = which(degree(G, mode = 'out')==0 &
                                           degree(G, mode = 'in')==0 &
                                           segs$cn>0)
                             i = 0
                             if (length(psimp)>0)
                             {
                                 vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
                                 epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA))
                                 ## there is no "edge" associated with a zero total degree node
                                 cns[1:length(psimp)] = segs$cn[psimp]
                                 i = length(psimp)
                             }

                             ## now iterate from shortest to longest path
                             ## peel that path off and see if it is still there ..
                             ## and see if it is still there
                             ## peel off top path and add to stack, then update cn.adj

                             segs$tile.id = get.tile.id(segs)

                             tile.map = data.table(id = seq_along(segs), tile.id = segs$tile.id)
                             rtile.map = data.table(id = seq_along(segs), tile.id = segs$tile.id)
                             setkey(tile.map, id)
                             setkey(rtile.map, tile.id)

                             ## unique pair of edge ids: rev comp of a foldback edge will be identical to itself!!!
                             ed = data.table(private$es)[cn>0, .(from, to , cn)]

                             if (nrow(ed)==0) ## make trivial gwalk of reference segments
                             {
                                 ss = segs[segs$loose == FALSE, ]
                                 paths = split(ss[, 'tile.id'], 1:length(ss))
                                 values(paths)$ogid = 1:length(paths)
                                 values(paths)$cn = ss$cn
                                 values(paths)$label = paste('CN=', ss$cn, sep = '')
                                 values(paths)$is.cycle = FALSE
                                 values(paths)$numsegs = elementNROWS(paths)
                                 values(paths)$num.ab = 0
                                 values(paths)$wid = width(ss)
                             }
                             else {

                                 ed[, ":="(fromss = tile.map[ .(from), tile.id],
                                           toss = tile.map[ .(to), tile.id]),
                                    by = 1:nrow(ed)]
                                 ed[, weight :=  adj[cbind(from, to)]]
                                 print(ed)
                                 ed[fromss*toss > 0,
                                    eclass := ifelse(fromss>0,
                                                     paste(fromss, toss),
                                                     paste(-toss, -fromss))]
                                 ed[fromss*toss < 0,
                                    eclass := ifelse(abs(fromss)<=abs(toss),
                                                     paste(fromss, toss),
                                                     paste(-toss, -fromss))]
                                 ed[, eclass := as.numeric(as.factor(eclass))]
                                 ed[, eid := paste(from, to)]
                                 setkey(ed, "eid")
                                 eclass.cn = ed[!duplicated(eclass), setNames(cn, eclass)]

                                 cleanup_mode = FALSE


                                 while (nrow(ij)>0)
                                 {
                                     if (verbose)
                                         message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' end-pairs to resolve' )
                                     i = i+1
                                     p = get.constrained.shortest.path(cn.adj,
                                                                       G,
                                                                       v = ij[1, 1],
                                                                       to = ij[1, 2],
                                                                       weight = E(G)$weight,
                                                                       edges = ed,
                                                                       verbose = TRUE,
                                                                       mip = cleanup_mode,
                                                                       gurobi = gurobi)

                                     if (is.null(p)){
                                         message('Came up empty!')
                                         i = i -1
                                         ij = ij[-1, , drop = FALSE]
                                     }
                                     else
                                     {
                                         ## Don't forget to update ed here
                                         ed$cn = cn.adj[cbind(ed$from, ed$to)]

                                         vpaths[[i]] = p
                                         epaths[[i]] = cbind(p[-length(p)], p[-1])
                                         eids = paste(epaths[[i]][,1], epaths[[i]][,2])
                                         cns[i] = ed[.(eids),
                                                     if (length(cn)>1) cn/2
                                                     else cn, by = eclass][, floor(min(V1))]
                                         ## update cn correctly
                                         ## adjust constraints for palinrdom edges by 1/2

                                         rvpath = rtile.map[list(tile.map[list(vpaths[[i]]),
                                                                          -rev(tile.id)]), id]
                                         repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
                                         plen = length(rvpath)
                                         hplen = floor(length(rvpath)/2)

                                         ## check for palindromicity for odd and even length palindromes
                                         ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
                                         if (ed[eids, any(table(eclass)>1)])
                                             palindromic.path[i] = TRUE

                                         vpaths[[i+1]] = rvpath
                                         epaths[[i+1]] = repath
                                         cns[i+1] = cns[i]
                                         palindromic.path[i+1] = TRUE

                                         ## so we want to update the current adjacency matrix
                                         ## to remove that path
                                         ## while keeping track of of the paths on the stack
                                         cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]
                                         cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]
                                         ## if (any(is.na(cn.adj))){
                                         ##     browser()
                                         ## }
                                         if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
                                         {
                                             ## maybe we got stuck in a quasi-palindrome and backtrack
                                             message('backtracking ...')

                                             cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]
                                             cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
                                             i = i-1
                                             ij = ij[-1, , drop = FALSE]
                                         }
                                         else ## continue, reduce
                                         {
                                             adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
                                             ## if (!palindromic)
                                             adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

                                             to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,
                                                                 drop = FALSE]
                                             ## if (!palindromic) ## update reverse complement
                                             to.rm = rbind(to.rm, epaths[[i+1]][which(cn.adj[epaths[[i+1]]]==0), ,drop = FALSE])

                                             if (nrow(to.rm)>0)
                                             {
                                                 adj[to.rm] = 0
                                                 ## ALERT!!! major change
                                                 ## adjj = adj/as.matrix(cn.adj)
                                                 ## adjj[which(is.nan(adjj))] = 0
                                                 ## adjj[which(adjj<0)] = 0
                                                 G = graph.adjacency(adj, weighted = 'weight')
                                                 ## G = graph.adjacency(adjj, weighted = 'weight')
                                                 new.ends = setdiff(which(
                                                 (degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
                                                 & degree(G)>0), ends)

                                                 ## ## check if cn.adj out of balance
                                                 ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
                                                 ##     print("Junction OUT OF BALANCE!")
                                                 ##     browser()
                                                 ## }

                                                 ## ## should be no new ends
                                                 ## if (length(new.ends)>0){
                                                 ##     print("Please, no new ends!")
                                                 ##     browser()
                                                 ## }

                                                 ## remain = as.matrix(jab$adj) - adj.new
                                                 ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
                                                 ## if (any(!is.element(nb, nb.all)))
                                                 ##     browser()

                                                 D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
                                                 ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
                                             }
                                             else
                                                 ij = ij[-1, , drop = FALSE]

                                             ## if (!palindromic) ## increase extra counter to account for reverse complement
                                             ## TOFIX: just update counter by 2 above, since we are just doing every path and its rc
                                             i = i+1
                                         }
                                     }


                                     ## DEBUG DEBUG DEBUG
                                     seg.ix = which(as.character(strand(segs))=='+'); seg.rix = which(as.character(strand(segs))=='-');


                                     if (nrow(ij)==0 & cleanup_mode == FALSE)
                                     {
                                         message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR PATHS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                                         ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
                                         cleanup_mode = TRUE
                                     }
                                 }
                                 if (verbose)
                                     message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )

                                 ## ## record G, D, remaining edges at the end of path peeling
                                 ## G1 = G
                                 ## D1 = D
                                 ## remain1 = remain

                                 vpaths = vpaths[1:i]
                                 epaths = epaths[1:i]
                                 cns = cns[1:i]
                                 palindromic.path = palindromic.path[1:i]

                                 vcycles = rep(list(NA), maxrow)
                                 ecycles = rep(list(NA), maxrow)
                                 ccns = rep(NA, maxrow)

                                 csimp = which(diag(cn.adj)!=0)
                                 ipath = i
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
                                 ## ALERT!!! major change
                                 ## adjj = adj/as.matrix(cn.adj)
                                 ## adjj[which(is.nan(adjj))] = 0
                                 ## adjj[which(adjj<0)] = 0
                                 G = graph.adjacency(adj, weighted = 'weight')
                                 ## G = graph.adjacency(adjj, weighted = 'weight')
                                 D = shortest.paths(G, mode = 'out', weight = E(G)$weight)

                                 ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]


                                 ## now iterate from shortest to longest path
                                 ## peel that path off and see if it is still there ..
                                 ## and see if it is still there

                                 ## peel off top path and add to stack, then update cn.adj

                                 cleanup_mode = FALSE
                                 while (nrow(ij)>0)
                                 {
                                     if (verbose)
                                         message('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
                                     i = i+1
                                        #        p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])

                                     p = get.constrained.shortest.path(cn.adj, G, allD = D, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode, gurobi = gurobi)

                                     if (is.null(p)){
                                         message('Came up empty!')
                                         i = i -1
                                         ij = ij[-1, , drop = FALSE]
                                     } else
                                     {

                                         ed$cn = cn.adj[cbind(ed$from, ed$to)]
                                         vcycles[[i]] = p
                                         ecycles[[i]] = cbind(p, c(p[-1], p[1]))
                                         eids = paste(ecycles[[i]][,1], ecycles[[i]][,2])
                                         ccns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palindromic edges by 1/2

                                         rvcycle = rtile.map[list(tile.map[list(vcycles[[i]]), -rev(tile.id)]), id]
                                         recycle = cbind(rvcycle, c(rvcycle[-1], rvcycle[1]))
                                         clen = length(rvcycle)
                                         hclen = floor(length(rvcycle)/2)
                                         ## (awkward) check for palindromicity for odd and even length palindromes

                                         ## if (all((vcycles[[i]]==rvcycle)[c(1:hclen,(clen-hclen+1):clen)]))
                                         if (ed[eids, any(table(eclass)>1)])
                                             palindromic.cycle[i] = TRUE
                                         ## else
                                         ## {
                                         vcycles[[i+1]] = rvcycle
                                         ecycles[[i+1]] = recycle
                                         ccns[i+1] = ccns[i]
                                         palindromic.cycle[i+1] = TRUE
                                         ##     palindromic = FALSE
                                         ## }
                                         ##        palindromic = TRUE ## set to true while we "figure things out"

                                         #' so now we want to subtract that cn units of that path from the graph
                                         #' so we want to update the current adjacency matrix to remove that path
                                         #' while keeping track of of the cycles on the stack
                                         cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
                                         ## if (!palindromic) ## update reverse complement unless palindromic
                                         cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]
                                         ## if (any(is.na(cn.adj))){
                                         ##     browser()
                                         ## }
                                         if (!all(cn.adj[ecycles[[i]]]>=0))
                                         {
                                             message('backtracking')
                                             ## browser()
                                             cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]+ccns[i]
                                             ## if (!palindromic) ## update reverse complement unless palindromic
                                             cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]+ccns[i+1]
                                             i = i-1
                                             ij = ij[-1, , drop = FALSE]
                                         }
                                         else
                                         {
                                             adj.new[ecycles[[i]]] = adj.new[ecycles[[i]]] + ccns[i]

                                             ## ## if (!palindromic)
                                             ##     adj.new[ecycles[[i+1]]] = adj.new[ecycles[[i+1]]] + ccns[i]

                                             ## ## ## make sure I didn't overuse any edge
                                             ## ## if (length(overdue <- which((as.matrix(jab$adj)-adj.new)<0))) {
                                             ## ##     print("Edge copy deficit!")
                                             ## ##     browser()
                                             ## ## }

                                             ## ## ## intermediate cross check
                                             ## ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                                             ## ##     browser()

                                             to.rm = ecycles[[i]][which(cn.adj[ecycles[[i]]]==0), ,drop = FALSE]

                                             ## if (!palindromic) ## update reverse complement
                                             to.rm = rbind(to.rm, ecycles[[i+1]][which(cn.adj[ecycles[[i+1]]]==0), ,drop = FALSE])

                                             if (nrow(to.rm)>0)
                                             {
                                                 adj[to.rm] = 0
                                                 parents = .parents(adj)
                                                 ## G = graph.adjacency(adj, weighted = 'weight')

                                                 ## ALERT!!! major change
                                                 ## adjj = adj/as.matrix(cn.adj)
                                                 ## adjj[which(is.nan(adjj))] = 0
                                                 ## adjj[which(adjj<0)] = 0
                                                 G = graph.adjacency(adj, weighted = 'weight')
                                                 ## G = graph.adjacency(adjj, weighted = 'weight')

                                                 ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
                                                 ##     print("Junction OUT OF BALANCE!")
                                                 ##     browser()
                                                 ## }

                                                 ## remain = as.matrix(jab$adj) - adj.new
                                                 ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
                                                 ## if (any(!is.element(nb, nb.all)))
                                                 ##     browser()

                                                 D = shortest.paths(G, mode = 'out', weight = E(G)$weight)
                                                 ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
                                             }
                                             else
                                                 ij = ij[-1, ,drop = FALSE]

                                             ## if (!palindromic) ## increase extra counter to account for reverse complement
                                             i = i+1
                                         }
                                     }

                                     if (nrow(ij)==0 & cleanup_mode == FALSE)
                                     {
                                         message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR CYCLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                                         ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]

                                         cleanup_mode = TRUE
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

                                 ## ## record G, D, remaining edges at the end of cycle peeling
                                 ## G2 = G
                                 ## D2 = D
                                 ## remain2 = remain
                                 remain = as.matrix(self$get.adj()) - adj.new
                                 remain.ends = which(Matrix::colSums(remain)*Matrix::rowSums(remain)==0 & Matrix::colSums(remain)-Matrix::rowSums(remain)!=0)
                                 if (length(remain.ends)>0){
                                     if (verbose)
                                         message(length(remain.ends), "ends were not properly assigned a path. Do them.")
                                 }

                                 tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
                                 ix = which(Matrix::rowSums(is.na(tmp[, 1:2]))==0)

                                 if (length(ix)>0)
                                     adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
                                 else
                                     adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
                                 vix = munlist(vall)

                                 segs$node.id = 1:length(segs)
                                 pathsegs = segs[vix[,3]]
                                 pathsegs$grl.ix = vix[,1]

                                 ## abjuncs =  as.data.table(ab.edges)[,id := rep(1:(nrow(ab.edges)/2),2)*
                                 ##                                         rep(c(1, -1), each = nrow(ab.edges)/2)][
                                 ##     !is.na(from), ]
                                 abjuncs = as.data.table(ab.edges)
                                 abjuncs[, ":="(eid = paste(from, to),
                                                reid = paste(hb.map[as.character(to)],
                                                             hb.map[as.character(from)]))]
                                 abjuncs[, ":="(ix = 1:.N,
                                                rix = match(reid, eid))]
                                 abjuncs[, unique.ix := ifelse(rix>=ix, paste(ix, rix), paste(rix, ix))]
                                 abjuncs[, eclass := as.numeric(as.factor(unique.ix))]
                                 abjuncs[, iix := 1:.N, by=eclass]
                                 abjuncs[, id := -(iix-1.5)/0.5*eclass]
                                 abjuncs = abjuncs[, tag := structure(paste(from, to), names = id)]
                                 setkey(abjuncs, tag)

                                 ## annotate ab.id (if any) following each segment in each path
                                 pathsegs$ab.id = gr2dt(pathsegs)[ , ab.id := c(abjuncs[paste(node.id[-length(node.id)], node.id[-1]), id], NA), by = grl.ix][, ab.id]

                                 paths = split(pathsegs, vix[,1] )
                                 values(paths)$ogid = 1:length(paths)
                                 values(paths)$cn = ecn[as.numeric(names(paths))]
                                 values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
                                 values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
                                 values(paths)$numsegs = elementNROWS(paths)
                                 values(paths)$num.ab = sapply(paths, function(x) sum(!is.na(x$ab.id)))
                                 values(paths)$wid = sapply(lapply(paths, width), sum)

                                 check = which((adj.new - self$get.adj()) !=0, arr.ind = TRUE)

                                 if (length(check)>0)
                                     stop('Alleles do not add up to marginal copy number profile!')
                                 else if (verbose)
                                     message('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')
                             }

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

                             ## for gGnome compatibiliity
                             if (!grl)
                             {
                                 tmp.dt = as.data.table(copy(paths))[, pid := group_name][, nix := 1:.N, by =pid]
                                 setkeyv(tmp.dt, c('pid', 'nix'))

                                 ## mark nodes that precede a reference junction
                                 tmp.dt[, d.to.next := c((start-data.table::shift(end))[-1], NA), by = pid]
                                 tmp.dt[, d.to.next.neg := c((data.table::shift(start)-end)[-1], NA), by = pid]
                                 tmp.dt[, same.strand := c((strand==data.table::shift(strand))[-1], NA), by = pid]
                                 tmp.dt[, same.chrom := c((as.character(seqnames)==data.table::shift(as.character(seqnames)))[-1], NA), by = pid]
                                 tmp.dt[, last.node := 1:.N == .N, by = pid]
                                 tmp.dt[, before.ref :=
                                              (((d.to.next<=1 & d.to.next>=0 & strand == '+') |
                                                (d.to.next.neg<=1 & d.to.next.neg>=0 & strand == '-')
                                              ) & same.strand & same.chrom)]
                                 tmp.dt[is.na(before.ref), before.ref := FALSE]

                                 ## label reference runs of nodes then collapse
                                 .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), as.integer(NA))
                                 tmp.dt[, ref.run := .labrun(before.ref), by = pid]
                                 tmp.dt[, ref.run.last := data.table::shift(ref.run), by = pid]
                                 tmp.dt[is.na(ref.run) & !is.na(ref.run.last), ref.run := ref.run.last]
                                 tmp.dt[!is.na(ref.run), ref.run.id := paste(pid, ref.run)]

### TODO: store ab.ids in walks
                                        #tmp.dt[loose == TRUE, ref.run.id := NA] ## make sure loose ends stay ungrouped
                                 if (any(!is.na(tmp.dt$ref.run.id)))
                                 {
                                     collapsed.dt = tmp.dt[!is.na(ref.run.id), .(
                                                                                   nix = nix[1],
                                                                                   pid = pid[1],
                                                                                   seqnames = seqnames[1],
                                                                                   start = min(start),
                                                                                   end = max(end),
                                                                                   loose = FALSE,
                                                                                   strand = strand[1]
                                                                               ), by = ref.run.id]
                                     tmp.dt = rbind(
                                         tmp.dt[is.na(ref.run.id),
                                                .(pid, nix, seqnames, start, end, strand, loose)],
                                         collapsed.dt[, .(pid, nix, seqnames, start, end, strand, loose)])

                                 }

                                 ## concatenate back with nodes that precede a non reference junctiono
                                 setkeyv(tmp.dt, c('pid', 'nix'))

                                 tmp.gr = dt2gr(tmp.dt)
                                 tmp.segs = unique(tmp.gr)
                                 tmp.gr$seg.id = match(tmp.gr, tmp.segs)
                                 tmp.paths = split(tmp.gr$seg.id, tmp.gr$pid)
                                 tmp.vals = as.data.frame(values(paths[names(tmp.paths)]))

                                 names(tmp.paths) = ifelse(grepl('\\-', names(tmp.paths)), -1, 1)*as.numeric(gsub('\\D', '', names(tmp.paths)))
                                 ## gw = as(paths, "gWalks")
                                 ## gw$simplify()
                                 ## return(gw)
                                 tmp.segs = gr.fix(tmp.segs, sl)
                                 gw = gWalks$new(segs=tmp.segs,
                                                 paths=tmp.paths,
                                                 metacols=tmp.vals)
                                 return(gw)
                             }
                             return(paths)
                         }),
                     private = list(

                     ),
                     active = list())

## ============= generics of gGraph ============= ##
## ============= generics of gGraph ============= ##
## components
##
## Big problem: how to define extra arguments unique for each dispatch?
##
## setClass("igraph")
## setGeneric("components", function(x, ...) {
##     standardGeneric("components")
## })
## setMethod("components",
##           signature(x = "igraph", mode="character"),
##           function(x) {
##               igraph::components(x, mode="weak")
##           }
##           )
## setMethod("components",
##           c(x = "gGraph"),
##           function(x) {
##               x$components()
##           }
##           )


`%+%.gGraph` <- function(gg1, gg2){
    return(gg1$add(gg2))
}

setMethod("seqinfo",
          c(x = "gGraph"),
          function(x) {
              x$seqinfo
          }
          )

setMethod("seqinfo",
          c(x = "bGraph"),
          function(x) {
              x$seqinfo
          }
          )

#' @name length
#' The number of strongly connected components of the graph
#'
#' @param gGraph a \code{gGraph} object
#'
#' @return the number of strongly connected components in this graph
#'
#' @export
length.gGraph <- function(gGraph){
    ## input must be a gGraph!
    if (!inherits(gGraph, "gGraph")){
        stop("Error: Invalid input.")
    }
    if (is.null(gGraph$parts)){
        cs = gGraph$components()
    }
    return(gGraph$parts$no)
}

#' @name %+%.gGraph
#' Adding two \code{gGraph} instances
#'
#' @param gg1, gg2 instances of the \code{gGraph} class
#'
#' @return a copy of a new \code{gGraph} object that is the simple sum of the two inputs
#'
#' @export
`%+%.gGraph` <- function(gg1, gg2){
    return(gg1$add(gg2))
}

#' @inheritParams %+%.gGraph
#' @export
`%+%.bGraph` <- function(bg1, bg2){
    return(bg1$add(bg2))
}

#' get.constrained.shortest.path
get.constrained.shortest.path = function(cn.adj, ## copy number matrix
                                         G, ## graph with distances as weights
                                         allD=NULL, ## shortest path between all nodes in graph
                                         v,
                                         to,
                                         weight,
                                         edges,
                                         verbose = TRUE,
                                         mip = TRUE,
                                         gurobi = TRUE,
                                         cplex = !gurobi){
    if (is.null(allD)){
        allD = shortest.paths(G, mode="out", weights = weight)
    }

    v = as.numeric(v)
    to = as.numeric(to)

    if (is.infinite(allD[v, to]) | allD[v, to]==0){
        return(NULL)
    }

    edges$cn = cn.adj[cbind(edges$from, edges$to)]

    ## ASSUME: from, to are scalars, within node range, to is reachable from from
    ## ASSUME edges contains eid key and eclass mapping
    tmp.p = as.numeric(get.shortest.paths(G, from=v, to=to, "out", weights=weight)$vpath[[1]])
    tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
    tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
    tmp.eclass = edges[.(tmp.eid), eclass]


    ## the cn of this path is the max number of copies that the network will allow
    ## here we have to group by eclass, i.e. so if there are two edges from an eclass
    ## in a given path then we need to halve the "remaining copies" constraint
    tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))]


    edges[, rationed := cn<(tmp.pcn*2)]

    D.totarget = allD[, as.numeric(to)]
    edges[, distance_to_target :=  D.totarget[to]]
    edges = edges[!is.infinite(distance_to_target) & cn>0, ]

    rationed.edges = edges[rationed == TRUE, ]

    ## find overdrafted eclasses - meaning two instances in this path but only one remaining copy
    overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)

    first.overdraft = which(tmp.eclass %in% overdrafts.eclass & duplicated(tmp.eclass))[1]

    ## no overdrafts?, then return
    ## ALERT!!! tmp.pcn could be NA!!!!
    ## TODO: how to fix it!!!???
    ## if (is.na(first.overdraft) & is.na(tmp.pcn)){
    ##     if (verbose)
    ##         message('Shortest path is good enough!')
    ##     return(tmp.p)
    ## }
    if (is.na(first.overdraft) & tmp.pcn>0)
    {
        if (verbose){
            message('Shortest path is good enough!')
        }
        return(tmp.p)
    }

    if (!mip){
        return(NULL)
    }


    ## use MIP to find constrained path
    edges[, enum := 1:length(eid)]

    ## incidence matrix constraints + 1 for tmp.pcn
    A = sparseMatrix(edges$to, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges))) -
        sparseMatrix(edges$from, edges$enum, x = 1, dims = c(nrow(cn.adj), nrow(edges)))
    b = rep(0, nrow(A))
    b[v] = -1
    b[to] = 1

    ix = which(Matrix::rowSums(A!=0)!=0) ## remove zero constraints

    ## "ration" or reverse complementarity constraints
    tmp.constraints = edges[, list(e1 = enum[1], e2 = enum[2], ub = cn[1]), by = eclass]
    tmp.constraints = tmp.constraints[!is.na(e1) & !is.na(e2), ]

    R = sparseMatrix(rep(1:nrow(tmp.constraints), 2),
                     c(tmp.constraints$e1, tmp.constraints$e2),
                     x = 1, dims = c(nrow(tmp.constraints), nrow(edges)))
    Rb = tmp.constraints$ub

    ## minimize weight of path
    cvec = edges$weight

    if (cplex){
        res = Rcplex::Rcplex(cvec = cvec,
                             Amat = rbind(A[ix,], R),
                             bvec = c(b[ix], Rb),
                             sense = c(rep('E', length(ix)), rep('L', length(Rb))),
                             lb = 0,
                             vtype = "B",
                             objsense = 'min')
    } else {
        model = list(A = rbind(A[ix,], R),
                     obj = cvec,
                     sense = c(rep('=', length(ix)), rep('>=', length(Rb))),
                     rhs = c(b[ix], Rb),
                     vtype = "B",
                     modelsense = "min")
        params = list()
        res = gurobi::gurobi(model, params)
    }

    if (verbose){
        message('YES WE ARE DOING PROPER MIP!!!!')
    }

    if (res$status!=101)
    {
        if (verbose){
            message('No solution to MIP!')
        }

        return(NULL)
    }

    ## use igraph to sort these edges into a path, i.e. make simple graph with one path and extract it using igraph (lazy :)
    tmp.p = as.numeric(get.shortest.paths(graph_from_edgelist(edges[res$xopt!=0, cbind(from, to)]), v, to)$vpath[[1]])

    ## check if overdrafted
    if (verbose)
    {
        tmp.e = cbind(tmp.p[-length(tmp.p)], tmp.p[-1])
        tmp.eid = paste(tmp.e[, 1], tmp.e[, 2])
        tmp.eclass = edges[.(tmp.eid), eclass]
        tmp.pcn = edges[.(tmp.eid), if (length(cn)>1) cn/2 else cn, by = eclass][, min(V1)]
        overdrafts.eclass = intersect(names(which(table(tmp.eclass)==2)), rationed.edges$eclass)
        if (length(overdrafts.eclass)==0){
            message('No overdrafts after MIP')
        }
        else
        {
            message('Still overdraft!')
            ## browser()
        }
    }

    return(tmp.p)
}



################################
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
convex.basis = function(A,
                        interval = 80,
                        chunksize = 100,
                        exclude.basis = NULL,
                        exclude.range = NULL,
                        maxchunks = Inf,
                        verbose = F){
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

############################################
#' hGraph: haplotype gGraph (underdeveloped)
#' has to be balanced and one unique path through all segs
#'
#' @import gUtils
#' @import gTrack
#' @import igraph
#############################################
hGraph = R6::R6Class("hGraph",
                     inherit=bGraph,
                     public = list(),
                     private = list(),
                     active = list())




#' @name grl.match
#' Matching the GRanges elements in a GRangesList
#'
#' @param query the \code{GRanges} attempting to match
#' @param subject the \code{GRanges} being matched to
#' @param ordered if TRUE the order of elements are considered in comparison
#' @param ignore.strand if TRUE the strand info is ignored
#'
#' @return vector same length as query, containing the indices of the first hit in subject
grl.match = function(query, subject,
                     ordered=FALSE,
                     ignore.strand=FALSE,
                     verbose=FALSE){
    if (ignore.strand){
        if (verbose) message("Ignoring the strand info.")
        grl1 = gr.stripstrand(grl1)
        grl2 = gr.stripstrand(grl2)
        return(grl.match(grl1, grl2, ordered=ordered, ignore.strand=FALSE))
    }

    ## TODO: find out the speed optimization of 'match'
    dt1 = copy(gr2dt(grl.unlist(grl1)))
    dt2 = copy(gr2dt(grl.unlist(grl2)))

    ## LATER
}

#' rev.comp
#' Return the reverse complement of a given GRanges
#'
#' @param gr GRanges with strand specified
#'
#' @return reverse complement of the input
rev.comp = function(gr){
    strmap = setNames(c("+", "-"), c("-", "+"))
    if (!inherits(gr, "GRanges")){
        stop("Input must be GRanges.")
    } else if (!all(strand(gr) %in% strmap)) {
        stop("Input must be all strand specific.")
    }
    return(rev(gr.flipstrand(gr)))
}

## 2) converting a gwalks to
## ============= R6 gWalks class definition ============= ##
#' @export
setClass("gWalks")

################################################
#' @name gWalks-class
#' @title genomic contigs
#' @docType class
#' @description
#' This class represents a multiset of walks consist of stranded genomic ranges from 5' to 3'.
#' The component of a contig is represented by a numeric vector of the \code{GRanges} indices
#' in the collection of genomic ranges. A \code{data.table} is associated with each
#' \code{gWalks} instance requiring 3 columns: \code{cn}, \code{is.cycle}, and \code{str}.
#' \code{cn} represents the copy numbers of the contigs. \code{is.cycle} stands for whether the
#' contig is circular (if so it means there is an implicit connection from the last range to
#' the first range). \code{str} is an arbitrary assignment of either of \code{+} or \code{-}
#' indicating the which of the pair of reverse complement walks.
#'
#' It is able to represents any state of a set of DNA sequences based on a reference genome.
#'
#' In the following examples, \code{gw} is a gWalks object.
#'
#' @usage
#' \strong{Constructors:}
#' gWalks$new(grl = NULL,
#'            segs = NULL, paths = NULL,
#'            is.cycle = NULL, cn = NULL, str = NULL, metacols = NULL)
#'
#' \strong{Public fields:}
#' gw$segstats
#'
#' gw$path
#'
#' gw$values
#'
#' gw$edges
#'
#' \strong{Public methods:}
#'
#' @import R6
#' @import gUtils
#' @import gTrack
#'
#' @export
################################################
gWalks = R6::R6Class("gWalks",
                     public=list(
                         ## refG = "GENOME",
                         initialize = function(grl=NULL, segs=NULL, paths=NULL,
                                               is.cycle=NULL, cn=NULL, str=NULL,
                                               metacols = NULL){
                             if (!is.null(segs)){
                                 private$gwFromScratch(segs,
                                                       paths,
                                                       is.cycle,
                                                       cn,
                                                       str,
                                                       metacols=metacols)
                             }
                             else if (!is.null(grl)) {
                                 self$grl2gw(grl)
                             }
                             else {
                                 self$nullGWalks()
                             }
                         },

                         ## TODO: construct null gWalks
                         nullGWalks = function(){
                             return(self)
                         },

                         pairup = function(){
                             verbose = getOption("gGnome.verbose")
                             strmap = setNames(c("+", "-"), c("-", "+"))

                             if (self$isStrandPaired()){
                                 return(self)
                             }

                             try.segs = tryCatch(seg.fill(private$segs),
                                                 error = function(e) NULL) ## TODO

                             if (!is.null(try.segs)){
                                 new.ix = match(private$segs, try.segs)
                                 private$segs = try.segs
                                 private$paths =
                                     relist(new.ix[unlist(private$paths)],
                                            private$paths)
                             } else {
                                 warning("Failed to pair up segs. Return the input.")
                                 return(self)
                             }

                             ## the other side
                             rpaths = rpaths(private$segs, private$paths) ## MOMENT
                             opaths = private$paths
                             plen = length(private$paths)

                             ## register the rev.comp paths that were previously not
                             plen = length(private$paths)
                             if (any(rp.add <- !rpaths %in% private$paths)){
                                 if (verbose){
                                     warning(paste("Appending", plen, "missing rev comp paths."))
                                 }

                                 length(private$paths) = plen + sum(rp.add)
                                 private$paths[(plen+1):length(private$paths)] = rpaths[which(rp.add)]

                                 rp.mdata = copy(private$metacols[which(rp.add)])
                                 rp.mdata[, str := strmap[str]]
                                 ## TODO! What is the way to append rows to data.table??
                                 private$metacols = rbind(private$metacols, rp.mdata)
                             }

                             rp.map = c(ifelse(rp.add,
                                               ## if not registered, rc is appended
                                               plen+cumsum(rp.add),
                                               ## if registered, rc is in opaths
                                               match(rpaths, opaths)),
                                        which(rp.add))
                             private$rp.map = rp.map

                             ## rc path must have same CN, if not, reconcile to higher
                             if (any(private$metacols[rp.map, cn] != private$metacols[, cn])){
                                 if (verbose)
                                     warning(paste("Increasing minor walk copies so there is same dosage on both strand."))

                                 tmp.mc = copy(private$metacols)
                                 tmp.mc[, rp.cn := tmp.mc[rp.map, cn]]
                                 tmp.mc[cn != rp.cn, cn := pmax(cn, rp.cn)]
                                 private$metacols = tmp.mc[, -c("rp.cn")]
                             }
                             return(self)
                         },

                         set.seg.cn = function(){
                             mc.cn = private$metacols$cn
                             mc.cn[which(is.na(mc.cn))] = 0
                             ## amplitude of each walk
                             amp = rep(mc.cn,
                                       IRanges::elementNROWS(private$paths))
                             cns = table(rep(unlist(private$paths), amp))

                             private$segs$cn = ifelse(
                                 as.character(seq_along(private$segs)) %in% names(cns),
                                 as.numeric(cns[as.character(seq_along(private$segs))]),
                                 0)
                             return(self)
                         },

                         gw2gg = function(){
                             verbose = getOption("gGnome.verbose")
                             strmap = setNames(c("+", "-"), c("-", "+"))
                             gw = self$clone()$reduce()

                             if (!gw$isStrandPaired()){
                                 ## MARCIN EDIT: NOT SURE WHY THIS FAILS SOMETIMES
                                 ## first check the segs
                                 ## in case two strands are not both present: fill it in
                                 gw$pairup()
                             }

                             ## summarize the edges and nodes
                             gw$set.seg.cn()
                             es = gw$p2e()
                             if (is.null(es)){
                                 es = data.table(from = numeric(0),
                                                 to = numeric(0),
                                                 type = character(0))
                             }

                             pl = get.ploidy(gw$segstats)

                             ## NOTE: rest assured no seg info is lost!
                             gg = gGraph$new(segs = gw$segstats,
                                             es = es,
                                             ploidy = pl,
                                             purity = 1)

                             ## it must be junction balanced
                             return(as(gg, "bGraph"))
                         },

                         gw2grl = function(ix=NULL, mod=FALSE){
                             if (is.null(ix))
                                 ix = seq_along(private$paths)

                             segs = private$segs

                             if (!"tile.id" %in% colnames(values(segs))){
                                 ss = gr.stripstrand(segs)
                                 names(ss) = NULL
                                 ss = ss[!duplicated(ss)] %Q% (order(seqnames, start))
                                 ss$tile.id = seq_along(ss)

                                 mix = match(gr.stripstrand(segs[,c()]), ss[,c()])
                                 private$segs$tile.id = segs$tile.id = ss$tile.id[mix]
                             }

                             ## segs$tile.id = rep(LETTERS[1:23][1:(length(private$segs)/2)], 2)
                             grl = lapply(private$paths[ix],
                                          function(pt){
                                              return(segs[pt])
                                          })
                             grl = GRangesList(grl)
                             mcols(grl) = private$metacols[ix]
                             if (mod){
                                 private$.grl = grl
                             }
                             return(grl)
                         },

                         grl2gw = function(grl, mc.cores=1){
                             if (length(grl)==0) return(self)

                             ## TODO:
                             ## why would there be cn==0 cycles out of gwalks function?
                             if ("cn" %in% colnames(values(grl))){
                                 grl = grl[which(values(grl)$cn>0)]
                             }

                             ## first of all, both strand of one range must be present
                             ## if not double them
                             grs = grl.unlist(grl)
                             private$segs = grs

                             already.ix = cumsum(table(grs$grl.ix))[as.character(grs$grl.ix-1)]
                             already.ix[is.na(already.ix)] = 0
                             names(already.ix) = NULL
                             private$paths = split(grs$grl.iix + already.ix, grs$grl.ix)

                             mc = as.data.frame(values(grl))
                             if (ncol(mc)==0) {
                                 ## if no meta data columns, default
                                 mc = data.table(is.cycle = rep(FALSE, length(grl)),
                                                 cn = rep(1,length(grl)),
                                                 str=rep("+", length(grl)))
                             } else {
                                 mc = as.data.table(mc)
                                 if (!"cn" %in% colnames(mc)){
                                     mc[,cn:=rep(1,length(grl))]
                                 }
                                 if (!"is.cycle" %in% colnames(mc)){
                                     mc[,is.cycle := rep(FALSE, length(grl))]
                                 }
                                 if (!"str" %in% colnames(mc)){
                                     mc[,str:=rep("+",length(grl))]
                                 }
                             }
                             private$metacols = mc
                             private$.grl = grl
                             return(self)
                         },

                         gw2td = function(ix=NULL,
                                          colorful=FALSE,
                                          mc.cores=1,
                                          ...){
                             if (length(private$paths)==0){
                                 if (verbose <- getOption("gGnome.verbose")){
                                     warning("Nothing to plot!")
                                 }
                                 return(NULL)
                             }
                             grl = self$gw2grl(ix)
                             gts = gTrack(grl, draw.path=T)
                             if (colorful) {
                                 ## TODO
                             }
                             return(gts)
                         },

                         json = function(fn = ".", settings = NULL){
                             return(self$gw2js(fn, settings = settings))
                         },

                         ## MOMENT
                         ## TODO: dedup the cid
                         gw2js = function(filename = ".",
                                          simplify=TRUE,
                                          trim=TRUE,
                                          mc.cores=1,
                                          debug = numeric(0),
                                          settings = NULL){
                             "Convert a gWalks object to JSON format for viz."
                             verbose = getOption("gGnome.verbose")

                             ## BUG: why doesn't the default value for settings work??
                             if (is.null(settings)){
                                 settings = list(y_axis = list(name = "copy number"))
                             }

                             if (length(private$segs)==0){
                                 if (verbose){
                                     warning("Empty walk. Nothing to plot.")
                                 }
                                 return(NULL)
                             }

                             require(jsonlite)

                             ## TODO: match up the cids with the gGraph
                             if (grepl('\\.js(on)*$', filename)){
                                 ## if json path was provided
                                 basedir = dirname(filename)
                                 filename = basename(filename)
                             }
                             else if (filename==".") {
                                 ## default path was provided
                                 basedir = './'
                                 filename = "data.json"
                             } else {
                                 ## a directory was provided
                                 basedir = filename
                                 filename = paste(filename, 'data.json', sep = '/')
                             }

                             if (!file.exists(basedir)) {
                                 message('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }

                             ## no walk, just graph
                             if (length(private$paths)==0){
                                 if (length(private$segs)==0){
                                     stop("Empty walks, empty graph.")
                                 }
                                 gg = self$gw2gg()
                                 gg.js = gg$gg2js(filename = filename,
                                                  trim = trim,
                                                  mc.cores = mc.cores,
                                                  settings = settings)
                                 return(normalizePath(filename))
                             }

                             ## start processing
                             ## first, we prepare the state of the gw that we want to plot
                             gw = self$clone() ## a temporary clone
                             ## need simple walks??
                             if (simplify) {
                                 gw$simplify()
                             }

                             ## are the paths paired??
                             if (!gw$isStrandPaired()){
                                 gw$pairup()
                             }

                             ## Let's reduce dup and unused nodes
                             gw$reduce()

                             ## get the y values for nodes
                             grl = gw$gw2grl()
                             ys = draw.paths.y(grl)

                             ## force correct segs CN
                             gw$set.seg.cn()

                             ## if no es, compute it, we need it here to define loose nodes
                             if (is.null(gw$edges)){
                                 gw$p2e()
                             }

                             gg = gw$gw2gg()$decouple(mod=TRUE)$make.balance()
                             gg.js = gg$gg2js(save=FALSE, settings = settings)

                             segs = copy(gw$segstats)
                             ed = copy(gw$edges)

                             if (!"loose" %in% colnames(values(segs)) | is.null(ed)){
                                 tmp = etype(segs, ed, force=TRUE, both=TRUE)
                                 segs = tmp$segs
                                 ed = tmp$es
                             }

                             hb = hydrogenBonds(segs)
                             if (hb[, any(is.na(from) | is.na(to))]){

                             }
                             hb.map = hb[, setNames(from, to)]

                             ## NOTE: prepare the data.tables anyway. They are needed in path too.
                             ## filter out junk ref contigs
                             regular.chrs = tryCatch(data.table::fread(Sys.getenv("DEFAULT_REGULAR_CHR")),
                                                     error = function(e) NULL)
                             if (!inherits(regular.chrs, "data.table")){
                                 if (verbose) {
                                     warning("The DEFAULT_REGULAR_CHR is not in correct format. See https://software.broadinstitute.org/software/igv/chromSizes for more.")
                                 }
                                 return(NULL)
                             }

                             regular.sl = regular.chrs[, setNames(V2, V1)]
                             regsegs.ix = which(as.character(seqnames(segs))
                                                %in% names(regular.sl))

                             ## know who are loose
                             loose.ix = which(segs$loose==TRUE)

                             ## construct intervals
                             node.dt = data.table(oid = which(as.logical(strand(segs)=="+")))
                             node.dt[, rid := hb.map[as.character(oid)]]

                             ## only keep the segment nodes that are inside regular chromosomes
                             node.dt = node.dt[oid %in%
                                               which(segs$loose==FALSE &
                                                     as.character(seqnames(segs))
                                                     %in% names(regular.sl))]

                             ## interval IDs
                             node.dt[, iid := 1:.N]
                             setkey(node.dt, "iid")
                             ## interval properties
                             node.dt[, ":="(chr = as.character(seqnames(segs[oid])),
                                            start = start(segs[oid]),
                                            end = end(segs[oid]))]
                             ## mapping from node id to interval ids
                             node.map = node.dt[, c(setNames(iid, oid),
                                                    setNames(iid, rid))]

                             ## vertical is the copy number
                             node.dt[, y := segs$cn[oid]]
                             node.dt[, title := paste(iid, paste0("(",oid,"|",rid,")"))]
                             node.dt[, type := "interval"]
                             node.dt[, strand := "*"]

                             ## This is needed later when constructing path.json
                             node.dt.both = rbind(node.dt[, .(nid = oid, iid,
                                                              chr, start, end, y,
                                                              title, type, strand="+")],
                                                  node.dt[, .(nid = rid, iid,
                                                              chr, start, end, y,
                                                              title, type, strand="-")])
                             setkey(node.dt.both, "nid")

                             ## keep the node/edge naming consistent with walks
                             ## node.json = node.dt[, .(iid,
                             ##                         chromosome = chr,
                             ##                         startPoint = start,
                             ##                         endPoint = end,
                             ##                         y,
                             ##                         title,
                             ##                         type,
                             ##                         strand)]
                             node.json = gg.js$intervals
                             ## finished processing nodes

                             ## TMPFIX: remove NA edges .. not clear where these are coming from
                             ## but likely the result of trimming / hood, but then it's not balanced
                             ## mapping from type field to label in json
                             eType = setNames(c("REF", "ALT", "LOOSE"),
                                              c("reference", "aberrant", "loose"))

                             ## some edges are out of the scope of regular chrs
                             e.na.ix = ed[, which(is.na(from) |
                                                  is.na(to) |
                                                  !(from %in% regsegs.ix) |
                                                  !(to %in% regsegs.ix))]
                             ed.na = ed[e.na.ix, ]

                             ## if any edge left, process
                             if (nrow(ed)-length(e.na.ix)>0){
                                 if (any(e.na.ix)){
                                     ed = ed[-e.na.ix, ]
                                 }

                                 ## edge's unique identifier
                                 ed[, ":="(eid = paste(from, to, sep="-"),
                                           reid = paste(hb.map[as.character(to)],
                                                        hb.map[as.character(from)],
                                                        sep="-"))]

                                 ## ALERT: bc strandlessness, I only retained half of the edges
                                 ## to map edges in gwalks, we will need strandedness,
                                 ## so will retain everything
                                 ed[,":="(soStr = as.character(strand(segs[from])),
                                          siStr = as.character(strand(segs[to])))]

                                 ## compute eclass
                                 ed[, ":="(ix = 1:.N,
                                           rix = match(reid, eid))]
                                 ed[, unique.ix := ifelse(rix>=ix, paste(ix, rix), paste(rix, ix))]
                                 ed[, eclass := as.numeric(as.factor(unique.ix))]
                                 ed[, iix := 1:.N, by=eclass]

                                 ## metadata of edges
                                 ed[, ":="(so = node.map[as.character(from)],
                                           si = node.map[as.character(to)],
                                           so.str = ifelse(soStr=="+",1,-1),
                                           si.str = ifelse(siStr=="+",1,-1),
                                           title = "",
                                           type = eType[type],
                                           weight = cn,
                                           cid = eclass)]

                                 ed[, ":="(source = so*so.str,
                                           sink = -si*si.str)]
                                 ## hard set NA on loose
                                 ed[fromLoose==T, source := NA]
                                 ed[toLoose==T, sink := NA]

                                 ## NOTE: I will never ever manually create/parse a JSON from string myself in my lift
                                 ## ppl wrote JSON format to make things standardized and pain-free to use


                                 ## ## EDGE.JSON
                                 ## ed.json = ed[iix==1, ## only need half of edges
                                 ##              .(cid,
                                 ##                source,
                                 ##                sink,
                                 ##                title,
                                 ##                type,
                                 ## weight)]
                             } else {
                                 ## ed.json = data.table(cid = numeric(0),
                                 ##                        source = numeric(0),
                                 ##                        sink = numeric(0),
                                 ##                        title = character(0),
                                 ##                        type = character(0),
                                 ##                        weight = numeric(0))
                             }

                             ## ALERT: we are having unpaired edges!!!!!!!!! WHY???
                             ## TODO: investigate what step created unpaired edge
                             ## ed.json = ed.json[!is.na(cid)]
                             ed.json = gg.js$connections

                             ## PATH.JSON, must be a list
                             path.json =
                                 mclapply(which(gw$metaCols()$str=="+"),
                                          function(pti){
                                              if (is.null(names(gw$path))){
                                                  this.pname = pti
                                              } else if (any(is.na(names(gw$path)))) {
                                                  this.pname = pti
                                              } else {
                                                  this.pname = names(gw$path)[pti]
                                              }

                                              this.npath = gw$path[[pti]]
                                              mc = gw$metaCols()
                                              this.cyc = mc[pti, is.cycle]

                                              ## test loose location
                                              loose.n = which(this.npath %in% loose.ix)
                                              if (length(loose.n)>0){
                                                  if (this.cyc==TRUE){
                                                      ## browser()
                                                      return(NULL)
                                                  } else {
                                                      if (any(loose.n %in%
                                                              setdiff(seq_along(this.npath),
                                                                      c(1, length(this.npath))))){
                                                          ## browser()
                                                          return(NULL)
                                                      }
                                                  }
                                              }

                                              this.non.loose = which(!this.npath %in% loose.ix)

                                              ## MARCIN EDIT: FIX TO DEAL WITH LENGTH 1 CYCLES
                                              this.epath.eid = c()
                                              if (length(this.npath)>1)
                                              {
                                                  this.epath.eid =
                                                      paste(this.npath[1:(length(this.npath)-1)],
                                                            this.npath[2:length(this.npath)],
                                                            sep="-")
                                              }

                                              this.ys = ys[[pti]][this.non.loose]

                                              if (this.cyc){
                                                  ##if (this.cyc & length(this.npath)>1){
                                                  this.epath.eid = c(this.epath.eid,
                                                                     paste(this.npath[length(this.npath)],
                                                                           this.npath[1], sep="-"))
                                              }

                                              this.npath = this.npath[this.non.loose]
                                              this.ndt = data.table(nid = this.npath, y = this.ys)
                                              this.ndt = cbind(this.ndt,
                                                               node.dt.both[.(this.ndt$nid),
                                                                            .(iid,
                                                                              chromosome = chr,
                                                                              startPoint = start,
                                                                              endPoint = end,
                                                                              strand = as.character(strand(segs[this.npath])),
                                                                              title,
                                                                              type)])

                                              ## MARCIN EDIT: RETURN NULL IF NO NODES AFTER regularChr
                                              ## TRIMMING
                                              this.ndt = this.ndt[!is.na(iid) & !is.na(chromosome), ]
                                              if (nrow(this.ndt)==0)
                                                  return(NULL)

                                              this.nids.json = this.ndt[,.(iid = iid,
                                                                           chromosome,
                                                                           startPoint,
                                                                           endPoint,
                                                                           y,
                                                                           title,
                                                                           type,
                                                                           strand)]
                                              if (this.nids.json[, any(duplicated(iid))]){
                                                  ## dup.iid = this.nids.json[duplicated(iid), 1:.N] + this.nids.json[, max(iid)]
                                                  ## this.nids.json[duplicated(iid),
                                                  ##                iid := dup.iid]
                                                  dup.ix =
                                                      this.nids.json[, which(duplicated(iid))]
                                                  dedup.ix = dup.ix +
                                                      this.nids.json[, max(iid)]
                                                  this.nids.json[dup.ix, iid := dedup.ix]
                                              }

                                              ## ALERT: throwing away good edges
                                              ## just bc they are not in ed.dt
                                              ## EDIT BY MARCIN:
                                              ## SOME VALID PATHS WILL HAVE >=1 nodes and NO EDGES
                                              ## if (any(!this.epath.eid %in% ed.dt[, eid])) return(NULL)
                                              ## just remove any edges that are off the grid, should be good enough
                                              if (nrow(ed)>0){
                                                  this.epath.eid = this.epath.eid[this.epath.eid %in% ed[type!="LOOSE", eid]]
                                              }


                                              this.cids.json = data.table(cid = numeric(0),
                                                                          source = numeric(0),
                                                                          sink = numeric(0),
                                                                          title = character(0),
                                                                          type = character(0),
                                                                          weight = numeric(0))


                                              if (length(this.epath.eid)>0)
                                              {
                                                  this.pdt = data.table(eid = this.epath.eid)
                                                  this.pdt = merge(this.pdt, ed[type!="LOOSE"], by="eid")

                                                  if (nrow(this.pdt)>0){
                                                      this.cids.json = this.pdt[this.epath.eid,
                                                                                .(cid = cid,
                                                                                  source,
                                                                                  sink,
                                                                                  title,
                                                                                  type,
                                                                                  weight)]
                                                      if (this.cids.json[, any(duplicated(cid))]){
                                                          ## dup.cid = this.cids.json[duplicated(cid), 1:.N] + this.cids.json[, max(cid)]
                                                          ## this.cids.json[duplicated(cid),
                                                          ##                cid := dup.cid]
                                                          dup.ix =
                                                              this.cids.json[
                                                                , which(duplicated(cid))]
                                                          dedup.ix =
                                                              dup.ix +
                                                              this.cids.json[, max(cid)]
                                                          this.cids.json[dup.ix, cid := dedup.ix]
                                                      }
                                                  }
                                              }

                                              this.mc = as.list(mc[pti,
                                                                   .(pid=pti,
                                                                     cn,
                                                                     type=ifelse(is.cycle,"cycle","path"),
                                                                     strand=str)])

                                              ## dedup the ids


                                              this.mc$cids = this.cids.json
                                              this.mc$iids = this.nids.json

                                              ## if (this.cids.json[, any(duplicated(cid))] |
                                              ##     this.nids.json[, any(duplicated(iid))])
                                              ## {
                                              ##     browser()
                                              ## }
                                              return(this.mc)
                                          },
                                          mc.cores=mc.cores)

                             path.json = path.json[which(
                                 sapply(path.json,
                                        function(x){
                                            out = valid = !is.null(x)
                                            if (valid){
                                                pos = x$strand=="+"
                                                out = valid & pos ## only plotting half of the walks, "+"
                                                return(out)
                                            } else {
                                                return(valid)
                                            }
                                        })
                             )]

                             out.json = list(settings = settings,
                                             intervals = node.json,
                                             connections = ed.json,
                                             walks = path.json)

                             if (verbose <- getOption("gGnome.verbose")){
                                 message("Writing JSON to ",
                                         paste(normalizePath(basedir),filename, sep="/"))
                             }

                             filename = paste(normalizePath(basedir),filename, sep="/")
                             jsonlite::write_json(out.json,
                                                  filename,
                                                  auto_unbox=TRUE, digits=4, pretty=TRUE)

                             return(normalizePath(filename))
                         },

                         v2e = function(mc.cores=1){
                             ## converting default node path into edges paths
                             if (is.null(private$gg)){
                                 gg = self$gw2gg()
                             } else {
                                 gg = private$gg
                             }
                         },

                         ## DONE: helper function to turn paths into edges
                         p2e = function(mc.cores=1){
                             ## whenever this function runs, it will assign result to
                             ## private$es, which will be refreshed to NULL whenever
                             ## a modifying action happens
                             es = do.call(
                                 'rbind',
                                 mclapply(1:length(private$paths),
                                          function(i){
                                              ## There might be NA in path cn
                                              if (private$metacols[i, is.na(cn)]){
                                                  return(NULL)
                                              }

                                              if (private$metacols[i, cn==0]){
                                                  return(NULL)
                                              }

                                              thisPath = private$paths[[i]]
                                              if (length(thisPath)>1){
                                                  ll = length(thisPath)
                                                  thisFrom = thisPath[1:(ll-1)]
                                                  thisTo = thisPath[2:ll]

                                                  if (private$metacols$is.cycle[i]){
                                                      thisFrom[ll] = thisPath[ll]
                                                      thisTo[ll] = thisPath[1]
                                                  }
                                              } else {
                                                  if (private$metacols$is.cycle[i]==TRUE){
                                                      thisFrom = thisTo = thisPath
                                                  } else {
                                                      return(NULL)
                                                  }
                                              }

                                              thisEs = data.table(from = thisFrom,
                                                                  to = thisTo,
                                                                  cn = private$metacols$cn[i],
                                                                  type = "unknown",
                                                                  path.ix = i)
                                              if (thisEs[, any(is.na(to) | is.na(from))]){
                                                  ## browser()
                                                  warning("Path ", i, " invalid, discard.")
                                                  return(NULL)
                                              }

                                              return(thisEs)
                                          },
                                          mc.cores=mc.cores)
                             )

                             ## up until this point, the edge cns are still identical
                             ## if same edges shows up more than once, dedup and populate cn
                             if (!is.null(es)){
                                 es[, eid := paste(from, to)]
                                 es[, cn := sum(cn), by=eid]
                                 path.ix.all =
                                     es[
                                       , .(path.all = paste(path.ix, collapse=","))
                                       , keyby=eid]
                                 es = merge(es, path.ix.all, by="eid")
                                 es = es[
                                     !duplicated(eid), .(from, to, cn, type, path.ix = path.all)]

                                 tmp = etype(private$segs, es, force = TRUE, both=TRUE)
                                 es = tmp$es
                                 private$segs = tmp$segs
                             }

                             private$es = es
                             return(es)
                         },

                         epath = function(mc.cores=1){
                             epaths = mclapply(seq_along(private$paths),
                                               function(i){
                                                   this.npath = private$paths[[i]]
                                                   nl = length(this.npath)
                                                   if (length(this.npath)>1){
                                                       this.epath.eids = paste(this.npath[1:(nl-1)])
                                                   } else {
                                                       this.epath.eids = character(0)
                                                   }

                                                   this.cyc = private$metacols[i, is.cycle]
                                                   if (this.cyc){
                                                       this.epath.eids[nl] = paste(this.npath[nl], this.npath[1])
                                                   }

                                               },
                                               mc.cores=mc.cores)
                             return(epaths)
                         },

                         simplify = function(mod = TRUE,
                                             reorder = FALSE,
                                             mc.cores=1,
                                             reduce = FALSE,
                                             debug=numeric(0)){
                             verbose = getOption("gGnome.verbose")

                             ## TODO: merge ref connected segments into one big
                             ## MOMENT
                             if (is.null(private$es))
                                 es = self$p2e()
                             else
                                 es = private$es

                             es = etype(private$segs, es, force=T)

                             if (es[, !any(type=="reference")]){
                                 return(self)
                             }

                             new.segs = private$segs[, c()]
                             setkey(es, eid)

                             if (verbose){
                                 warning("Updating segs, metadata fields are modified. Proceed with care.")
                             }

                             new.paths =
                                 mclapply(seq_along(private$paths),
                                          function(e){
                                              pth = private$paths[[e]]
                                              ## find out runs of at least one reference edge
                                              ep = data.table(from = shift(pth),
                                                              to = pth)[-1, ]
                                              ep[, eid := paste(from, to)]
                                              ep[, type := es[.(ep$eid), type]]

                                              n.loose = ep[, sum(type=="loose")]

                                              ## why NA?
                                              if (!any(ep$type=="reference")){
                                                  return(pth)
                                                  next
                                              }

                                              ref.eix = ep[, rle(type)]
                                              ref.run = which(ref.eix$values=="reference")

                                              until = cumsum(ref.eix$lengths)[ref.run]
                                              if (identical(ref.run, 1L)){
                                                  since = 1
                                              } else if (ref.run[1]==1) {
                                                  since = c(1, cumsum(ref.eix$lengths)[ref.run-1]+1)
                                              } else {
                                                  since = cumsum(ref.eix$lengths)[ref.run-1]+1
                                              }

                                              ep[, run.ix := 0]
                                              tmp = c(0, seq_along(since))
                                              for (i in seq_along(since)){
                                                  ## TODO: a better way than for loop???
                                                  ep[since[i]:until[i], run.ix := i]
                                              }

                                              to.merge = ep[run.ix>0,
                                                            .(n.ix = c(head(from, 1), to)),
                                                            by=run.ix]

                                              ## create new node
                                              new.node =
                                                  do.call(`c`,
                                                          sapply(unique(to.merge$run.ix),
                                                                 function(ix){
                                                                     ns = reduce(
                                                                         private$segs[to.merge[run.ix==ix, n.ix]])
                                                                     ns$run.ix = ix
                                                                     return(ns)
                                                                 })
                                                          )

                                              ## ALERT: modifying segs
                                              new.node$new.ix = new.node$run.ix + length(new.segs)
                                              if (e %in% debug){
                                                  browser()
                                              }
                                              new.segs <<- c(new.segs, new.node[,c()])
                                              ## was this evaluated?
                                              new.nix = new.node$new.ix[new.node$run.ix]

                                              ## modifying paths
                                              ep[, ":="(last.run.ix = shift(run.ix),
                                                        next.run.ix = c(tail(run.ix, -1), NA))]

                                              ep[run.ix==0 & last.run.ix!=0,
                                                 from := new.nix[last.run.ix]]
                                              ep[run.ix==0 & next.run.ix!=0,
                                                 to := new.nix[next.run.ix]]
                                              new.ep = ep[run.ix==0]

                                              new.n.loose = new.ep[, sum(type=="loose")]
                                              if (new.n.loose != n.loose){
                                                  ## browser()
                                              }

                                              if (nrow(new.ep)==0){
                                                  new.path = new.nix
                                              } else {
                                                  new.path = new.ep[, c(head(from, 1), to)]
                                              }

                                              if (any(is.na(new.path))){
                                                  ## browser()
                                              }
                                              return(new.path)
                                          },
                                          mc.cores = mc.cores)

                             if (reduce){
                                 self$reduce(mod=TRUE)
                             }

                             if (mod==TRUE){
                                 ## modify required fields THIS instance
                                 private$segs = new.segs
                                 private$paths = new.paths
                                 ## reset all optional fields
                                 private$es = NULL
                                 private$.grl = NULL
                                 private$gg = NULL
                                 private$rp.map = NULL
                                 return(self)
                             } else {
                                 out = gWalks$new(segs = new.segs,
                                                  paths = new.paths,
                                                  metacols = private$metacols)
                                 return(out)
                             }
                         },

                         reduce = function(mod=TRUE){
                             ## final step, dedup the segments, relabel the paths
                             ## this deduping is based on ranges only
                             new.segs = private$segs
                             new.paths = private$paths

                             ## loose is just a property field of the GRanges
                             ## IDEA: we can allow more property fields in other classes
                             ## like hGraph can have "alleles"
                             ## TODO: extend GRanges to be skew-symmetric
                             ## if ("loose" %in% colnames(values(private$segs))){
                             ##     segments.ix = which(new.segs$loose==FALSE)
                             ##     loose.ix = which(new.segs$loose==TRUE)

                             ##     segments = new.segs[segments.ix]
                             ##     loose.ends = new.segs[loose.ix]

                             ##     seg.dup = which(duplicated(segments))
                             ##     loose.dup = which(duplicated(loose.ends))

                             ##     dup <- c(segments.ix[seg.dup],
                             ##              loose.ix[loose.dup])
                             ## } else {
                             dup = which(duplicated(new.segs))
                             ## }

                             if (length(dup)>0){
                                 ## if ("loose" %in% colnames(values(private$segs))){
                                 ##     nr.segments = segments[setdiff(seq_along(segments),
                                 ##                                    seg.dup)]
                                 ##     nr.loose.ends = loose.end[setdiff(seq_along(loose.ends),
                                 ##                                       loose.dup)]
                                 ##     match(segments, nr.segments)
                                 ##     nr.map = match(segments, )
                                 ##     ## MOMENT

                                 ## } else {
                                 nr.segs = new.segs[-dup]
                                 nr.map = match(new.segs, nr.segs)
                                 ## }

                                 new.paths = lapply(new.paths, function(x) nr.map[x])
                                 new.segs = nr.segs
                             }

                             ## discard any old seg that's not used anymore
                             if (length(unused <- setdiff(seq_along(new.segs),
                                                          do.call(`c`, new.paths)))>0){
                                 used.segs = new.segs[-unused]
                                 used.map = match(new.segs, used.segs)
                                 new.paths = lapply(new.paths, function(x) used.map[x])
                                 new.segs = used.segs
                             }

                             if (mod){
                                 private$segs = new.segs
                                 private$paths = new.paths
                                 private$reset()
                                 return(self)
                             } else {
                                 out = gWalks$new(segs = new.segs,
                                                  paths = new.paths,
                                                  metacols = private$metacols)
                                 return(out)
                             }
                         },

                         subset = function(ix){
                             new.paths = private$paths[ix]
                             new.mc = private$metaCols[ix,]
                             return(gWalks$new(segs = private$segs,
                                               paths = new.paths,
                                               metacols = new.mc))
                         },

                         print = function(){
                             str = paste0("gWalks:\n",
                                          "\t", length(private$paths), " contigs\n",
                                          "\t", length(private$segs), " ranges\n")
                             cat(str)
                         },

                         len = function(){
                             return(length(private$paths))
                         },

                         metaCols = function(){
                             return(private$metacols)
                         },

                         window = function(ix = NULL, pad=1e3){
                             if (is.null(ix)){
                                 ix = seq_along(private$paths)
                             }
                             ix = unlist(private$paths[ix])
                             return(trim(streduce(private$segs[ix]) + pad))
                         },

                         plot = function(){
                             td = self$gw2td()
                             win = self$window()
                             plot(td, win)
                         },

                         ## tests
                         isStrandPaired = function(){
                             ## TODO: need to update the matching between reverse comps
                             ## check point 1
                             segs = private$segs[, c()]
                             names(segs) = NULL
                             if (any(is.na(
                                 match(gr.stripstrand(segs %Q% (strand=="+")),
                                       gr.stripstrand(segs %Q% (strand=="-")))
                             ))){
                                 return(FALSE)
                             } ## else if (any(is.na(match(private$paths,
                             ##                          self$rpaths())))){
                             ##   return(FALSE)
                             ## }
                             return(TRUE)
                         },

                         label = function(new.ord=NULL, mod=TRUE){
                             verbose = getOption("gGnome.verbose")
                             "Relabel the nodes."
                             seg.dt = gr2dt(private$segs)
                             if (is.null(new.ord)){
                                 if (verbose)
                                     message("Will sort the nodes by 'loose', 'strand', 'seqnames', 'start'.")

                                 if ("loose" %in% colnames(seg.dt))
                                     new.ord = seg.dt[, order(loose, strand, seqnames, start)]
                                 else
                                     new.ord = seg.dt[, order(strand, seqnames, start)]
                             }

                             if (length(new.ord) != length(private$segs))
                                 stop("Given new order is not the same length as the nodes.")

                             if (identical(new.ord, seq_along(private$segs)))
                                 return(self)

                             ## sort segs
                             new.segs = private$segs[new.ord]
                             ## relabel paths
                             new.paths = lapply(private$paths,
                                                function(x) new.ord[as.character(x)])

                             if (mod==TRUE){
                                 ## modify required fields THIS instance
                                 private$segs = new.segs
                                 private$paths = new.paths
                                 ## reset all optional fields
                                 private$es = NULL
                                 private$.grl = NULL
                                 private$gg = NULL
                                 private$rp.map = NULL
                                 return(self)
                             } else {
                                 return(gWalks$new(segs = new.segs, paths = new.paths,
                                                   metacols = private$metacols))
                             }
                         }
                     ),

                     private = list(
                         ## ----- private fields
                         ## ===== required
                         segs = NULL,
                         paths = NULL,
                         metacols = NULL,

                         ## ===== optional, dynamically created, turn back to NULL when being updated
                         es = NULL,
                         .grl = NULL,
                         gg = NULL,
                         rp.map = NULL,

                         ## ----- private methods
                         reset = function(){
                             private$es = NULL
                             private$.grl = NULL
                             private$gg = NULL
                             private$rp.map = NULL
                             return(self)
                         },

                         ## you can give me columns separately as vectors
                         ## or you can give me data.frame as a whole
                         gwFromScratch = function(segs, paths=NULL, is.cycle=NULL,
                                                  cn=NULL, str=NULL, metacols=NULL){
                             ## segs must be a GRanges
                             if (!inherits(segs, "GRanges")) stop("segs needs to be a GRanges.")

                             ## ALERT: sometimes you need to end up at a loose end
                             ## ## ... and disjoint

                             ## if (!isDisjoint(segs)) stop("segs must be disjoint.")

                             ## MARCIN COMMENT:
                             ## if loose ends missing from input segs then $loose flag never set
                             ## to default value, breaks downstream code
                             private$segs = segs

                             ## paths must be a list of numeric vectors
                             if (!is.null(paths)){
                                 if (!is.list(paths) |
                                     !all(sapply(paths, is.numeric)) |
                                     any(sapply(paths, max)>length(segs)) |
                                     any(sapply(paths, min)<1))
                                     stop("paths must be list of indices, within 1:length(segs)")
                                 private$paths = paths
                             } else {
                                 private$paths = list()
                             }

                             ## TODO: construct gw with meta cols or vectors
                             ## metacols always overwrite other arguments
                             ## if those particular columns exist
                             if (length(private$paths)>0){
                                 ## given metacols???
                                 if (!is.null(metacols)){
                                     if (is.data.frame(metacols) &
                                         nrow(metacols)==length(private$paths))
                                         private$metacols = as.data.table(metacols)
                                 }

                                 ## only populate metacols when paths is not empty
                                 ## given is.cycle??
                                 if (!is.null(is.cycle) & is.logical(is.cycle) &
                                     length(is.cycle)==length(paths))
                                     private$metacols$is.cycle = is.cycle
                                 else if (!"is.cycle" %in% colnames(private$metacols))
                                     private$metacols$is.cycle = rep(FALSE, length(private$paths))

                                 ## given cn??
                                 if (!is.null(cn) & is.numeric(cn) & length(cn)==length(paths))
                                     private$metacols$cn = cn
                                 else if (!"cn" %in% colnames(private$metacols))
                                     private$metacols$cn = rep(1, length(private$paths))

                                 ## given strand??
                                 if (!is.null(str) & is.numeric(str) & length(str)==length(paths))
                                     private$metacols$str = str
                                 else if (!"str" %in% colnames(private$metacols))
                                     private$metacols$str = rep("+", length(private$paths))

                                 private$metacols = as.data.table(private$metacols)
                             } else {
                                 private$metacols = data.table(is.cycle=logical(0),
                                                               cn=numeric(0),
                                                               str=character(0))
                             }
                             return(self)
                         }
                     ),
                     active=list(
                         ## notice the returned value is by reference
                         segstats = function(){
                             return(private$segs)
                         },
                         edges = function(){
                             if (is.null(private$es)){
                                 self$p2e()
                             }
                             return(copy(private$es))
                         },
                         grl = function(){
                             if (is.null(private$.grl)){
                                 private$.grl = self$gw2grl()
                             }
                             return(private$.grl)
                         },
                         td = function(){
                             ## default viz, will not show CN==0 or str=="-"
                             ix = private$metacols[, which(cn>0 & str=="+")]
                             return(self$gw2td(ix))
                         },
                         path = function(){
                             return(private$paths)
                         },
                         values = function(){
                             return(self$metaCols())
                         }
                     ))

length.gWalks <- function(x){
    return(length(x$path))
}

`[.gWalks` <- function(x, idx=NULL){
    if (is.null(idx)){
        idx = seq_len(length(x))
    }
    return(x$subset(idx))
}

`[<-.gWalks` <- function(x, idx, subs){

}

## ============= exported functions ============= ##
## explicit coercion and that's it!
setAs("GRangesList", "junctions",
      function(from){
          new("junctions", from)
      })

setAs("gGraph", "bGraph",
      function(from){
          return(bGraph$new(from))
      })

setAs("gWalks", "GRangesList",
      function(from){
          return(from$grl)
      })

setAs("gWalks", "bGraph",
      function(from){
          return(from$gw2gg())
      })

setAs("GRangesList", "gWalks",
      function(from){
          return(gWalks$new(grl = from))
      })

setAs("list", "gWalks",
      function(from){
          pre.grl = tryCatch(GRangesList(from),
                             error = function(e) NULL)
          if (is.null(pre.grl))
              stop("This list must be convertible to GRangesList.")
          else
              return(as(pre.grl, "gWalks"))
      })

############################################
#' @name clusters
#' Clustering the vertices in a grpah by connectivity
#'
#' @param gg a gGraph object
#' @param v indices of the subset of vertices to cluster
#' @param mode either "weak"-ly or "strong"-ly connected components
#' @param use.hb logical if TRUE attach hydrogen bonds before applying clusters
#'
#' @return list of subgraphs
#' @export
############################################
clusters = function(gg,
                    v = numeric(0),
                    mode = c("weak", "strong"),
                    use.hb = TRUE){
    return(gg$clusters(v, mode, use.hb))
}



############################################
#' @name isBalance
#' Test if a gGraph object is junction balanced
#'
#' @param gg a gGraph object
#'
#' @details
#' This function will check if all node copy numbers equal to the in-flow and out-flow
#' of edge copy numbers, except for terminal nodes (e.g. loose ends, telomeres), where
#' one side satisfies the equality and the other is zero. Singleton nodes, NA copy nodes
#' are ignored.
#'
#' @return a logical scalar
#' @export
###########################################
isBalance = function(gg){
    if (!inherits(gg, "gGraph")){
        stop("Input is not a gGraph.")
    }

    if (!"cn" %in% colnames(values(gg$segstats)) |
        !"cn" %in% colnames(gg$edges)){
        return(FALSE)
    }

    return(gg$isBalance())
}

############################################
#' @name gread
#' Parse the outputs from rearrangement graph callers.
#'
#' @param filename filename to JaBbA's rds, PREGO's intervalFile, or Weaver's output directory
#'
#' @details
#' This function will interpret ".rds" files as JaBbA output.
#'
#' @return a proper gGraph family instance
#' @export
###########################################
gread = function(filename){
    verbose = getOption("gGnome.verbose")

    if (is.list(filename)){
        bg = tryCatch(bGraph$new(jabba = filename),
                      error = function(e) NULL)
        if (!is.null(bg)){
            return(bg)
        } else {
            return(gGraph$new(jabba = filename))
        }
    }

    ## decide what output this is
    if (!file.exists(filename)){
        stop("Error: No such file or directory!")
    }

    if (dir.exists(filename)){
        if (verbose){
            message("Given a directory, assume it's Weaver.")
        }
        return(gGraph$new(weaver=filename))
    } else if (grepl(".rds$", filename, ignore.case=TRUE)){
        rds = tryCatch(readRDS(filename),
                       error=function(e)
                           stop("Given file can't be read as RDS."))

        if (inherits(rds, "gGraph")) {
            return(rds)
        }
        else if (inherits(rds, "list")){
            bg = tryCatch(bGraph$new(jabba = rds),
                          error = function(e) NULL)
            if (!is.null(bg)){
                return(bg)
            } else {
                return(gGraph$new(jabba = rds))
            }
        }
    } else if (grepl(".js[on]*$", filename)){
        ## TODO: what's the re for matching 0 or 1 time???

    }
    else {
        prego = gGraph$new(prego = filename)
        prego = as(prego, "bGraph")
        return(prego)
    }
}

############################################
#' @name refresh
#' Refreshing your object with the latest code.
#'
#' @param object any instance of the R6 classes in gGnome
#'
#' @return a proper gGraph family instance with the old data
#' @export
###########################################
refresh = function(object){
    if (inherits(object, "bGraph")){
        return(bGraph$new(object))
    } else if (inherits(object, "gGraph")){
        return(gGraph$new(segs = object$segstats,
                          es = object$edges,
                          purity = object$purity,
                          ploidy = object$ploidy))
    } else if (inherits(object, "gWalks")){
        return(as(object$grl, "gWalks"))
    } else {
        warning("Not a gGnome object.")
        return(NULL)
    }
}

#####################################
#' fusions
#'
#' annotates all gene fusions given an n x n adjacency matrix A of n genomic segments seg and grl of transcripts (eg output of read_RefGene)
#' seg must be (1) a tiling of the genome and (2) have copies of both + and - intervals for each genomic range (eg output of karyograph)
#'
#' alternate input is a pile of junctions of ranges with strand orientation pointing AWAY from breakend
#'
#' cds = gencode cds GRanges gff3 / gtf input
#'
#' "gene_name" GRangesList meta data field is used in annotation and in not creating "splice" fusions that arise from different transcripts of the same gene.
#'
#' @param gg a gGraph object (for example from JaBbA; overrides junctions input)
#' @param junc GRangesList of junctions (each a length 2 GRanges)
#' @param cds CDS annotations (GrangesList of transcript composed of coordinates coding regions of exons)
#' @param promoters GRanges of promoters (same length as transcript)
#' @param query optional query limiting walks to specific regions of interest
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @return GRangesList of walks corresponding to transcript boundaires
#'
#' @export
################################
fusions = function(gg = NULL,
                   junc = NULL,
                   cds = NULL,
                   promoters = NULL,
                   query = NULL,
                   prom.window = 1e3,
                   max.chunk = 1e10,
                   cb.interval = 1e4,
                   cb.chunksize = 1e4,
                   cb.maxchunks = 1e10,
                   exhaustive = FALSE,
                   debug = NULL,
                   mc.cores = 1,
                   verbose = NULL){
    if (is.null(verbose)){
        verbose <- getOption("gGnome.verbose")
    }
    if (!is.logical(verbose)){
        verbose <- FALSE
    }

    ## QC input graph or junctions
    if (!is.null(gg)){

        if (!inherits(gg, "gGraph")){

            if (is.list(gg) &
                !all(is.element(
                     c("segstats", "adj", "ab.edges","purity", "ploidy", "junctions"),
                     names(gg)))){
                stop("Invalid input graph.")
            }

            seg = gg$segstats
            A = gg$adj
        }

    } else if (!is.null(junctions)){
        if (verbose){
            message("Generating graph from junctions.")
        }
        gg = gGraph$new(junc = junc)
        A = gg$adj
        seg = gg$segstats
    } else{
        stop('either gg or junctions input must be non NULL')
    }

    if (is.null(A) | is.null(seg))
        stop('Some essential args are NULL')

    ## QC input cds
    if (!inherits(cds, 'GRangesList')){
        cds = NULL
    } else if (!all(c('Transcript_id', 'Gene_name') %in%
                    names(values(cds)))){
        cds = NULL
    }

    ## loading default cds
    if (is.null(cds))
    {
        if (verbose)
            message('CDS object missing or malformed (e.g. does not contain Transcript_id and Gene_name GRangesList metadata fields\nReading in from gencode CDS via Sys.getenv("DEFAULT_GENE_ANNOTATION")')
        cds = read_gencode(type="cds")
    }

    ## convert back to a GR?
    tx.span = gUtils::seg2gr(values(cds))

    names(values(tx.span))[match(c('transcript_id', 'gene_name'), tolower(names(values(tx.span))))] = c('transcript_id', 'gene_name')

    if (verbose){
        message('got transcript boundaries\n')
    }

    ## determine set of transcript fragments
    ## these correspond to transcripts that intersect a segment boundary
    cds.frag.left = gUtils::gr.findoverlaps(tx.span, gr.start(seg),
                                            qcol = c('gene_name', 'transcript_id'),
                                            ignore.strand = F, max.chunk = max.chunk)
    strand(cds.frag.left) = strand(tx.span)[cds.frag.left$query.id]
    cds.frag.right = gUtils::gr.findoverlaps(tx.span, gr.end(seg),
                                             qcol = c('gene_name', 'transcript_id'),
                                             ignore.strand = F, max.chunk = max.chunk)
    strand(cds.frag.right) = strand(tx.span)[cds.frag.right$query.id]

    ## I want to find all unique walks that involve tx fragments
    if (length(cds.frag.left)>0 & length(cds.frag.right) > 0 )
        tmp = merge(data.frame(i = 1:length(cds.frag.left), key1 = cds.frag.left$query.id, key2 = cds.frag.left$subject.id),
                    data.frame(j = 1:length(cds.frag.right), key1 = cds.frag.right$query.id, key2 = cds.frag.right$subject.id), all = T)
    else
        return(GRangesList())

    pos.right = which(as.logical( strand(cds.frag.right)=='+'))
    pos.left = which(as.logical(strand(cds.frag.left)=='+'))
    neg.right = which(as.logical(strand(cds.frag.right)=='-'))
    neg.left = which(as.logical( strand(cds.frag.left)=='-'))

    ## positive start fragments will be "right" fragments
    cds.start.frag.pos = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% pos.right, ]$j]
    start(cds.start.frag.pos) = start(tx.span)[cds.start.frag.pos$query.id]
    if (length(cds.start.frag.pos)>0)
        cds.start.frag.pos$type = 'start'

    ## positive end fragments will be "left" fragments
    cds.end.frag.pos = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% pos.left, ]$i]
    end(cds.end.frag.pos) = end(tx.span)[cds.end.frag.pos$query.id]
    if (length(cds.end.frag.pos)>0)
        cds.end.frag.pos$type = 'end'

    ## negative start fragments will be "right" fragments
    cds.start.frag.neg = cds.frag.left[tmp[is.na(tmp$j) & tmp$i %in% neg.left, ]$i]
    end(cds.start.frag.neg) = end(tx.span)[cds.start.frag.neg$query.id]
    if (length(cds.start.frag.neg)>0)
        cds.start.frag.neg$type = 'start'

    ## negative end fragments will be "left" fragments
    cds.end.frag.neg = cds.frag.right[tmp[is.na(tmp$i) & tmp$j %in% neg.right, ]$j]
    start(cds.end.frag.neg) = start(tx.span)[cds.end.frag.neg$query.id]
    if (length(cds.end.frag.neg)>0)
        cds.end.frag.neg$type = 'end'

    ## remaining will be "middle" fragments
    middle.frag = cds.frag.left[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$i]
    end(middle.frag) = end(cds.frag.right[tmp[!is.na(tmp$i) & !is.na(tmp$j),]$j])
    if (length(middle.frag)>0)
        middle.frag$type = 'middle'

    ## concatenate fragments
    ## subject.id of frags is the id of the node on the graph

    all.frags = c(cds.start.frag.pos, cds.end.frag.pos, cds.start.frag.neg, cds.end.frag.neg, middle.frag)

    ##
    ## now connect all.frags according to A
    ## i.e. apply A connections to our fragments, so draw an edge between fragments
    ## if
    ## (1) there exists an edge connecting segment and
    ## (2) only allowable connections are 'start' --> 'middle' --> 'middle' --> 'end'
    seg.edges = as.data.frame(which(A!=0, arr.ind = T))
    colnames(seg.edges) = c('from.seg', 'to.seg')
    edges = merge(merge(data.frame(i = 1:length(all.frags), from.seg = all.frags$subject.id),
                        seg.edges), data.frame(j = 1:length(all.frags), to.seg = all.frags$subject.id))

    edges = edges[all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'middle' |
                  all.frags$type[edges$i] == 'start' & all.frags$type[edges$j] == 'end' |
                  all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'middle' |
                  all.frags$type[edges$i] == 'middle' & all.frags$type[edges$j] == 'end', ]

    ## this removes splice variants .. keeping only links that fuse different genes or same transcripts
    edges = edges[tx.span$gene_name[all.frags$query.id[edges$i]] != tx.span$gene_name[all.frags$query.id[edges$j]] |
                  all.frags$query.id[edges$i] == all.frags$query.id[edges$j],]

    if (nrow(edges)==0){
        return(GRangesList())
    }

    if (verbose){
        cat('computed subgraph\n')
    }

    A.frag = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(length(all.frags),2))
    keep.nodes = which(Matrix::rowSums(A.frag)>0 | Matrix::colSums(A.frag)>0)
    A.frag = A.frag[keep.nodes, keep.nodes]
    all.frags = all.frags[keep.nodes]

    sources = which(all.frags$type == 'start')
    sinks = which(all.frags$type == 'end')

    G = graph.adjacency(A.frag)
    C = igraph::clusters(G, 'weak')
    vL = split(1:nrow(A.frag), C$membership)
    paths = do.call('c', mclapply(1:length(vL), function(i) {
        if (verbose & (i %% 10)==0){
            cat(i, ' of ', length(vL), '\n')
        }
        x = vL[[i]]
        if (!is.null(debug)){
            if (i %in% debug){
                browser()
            }
        }

        tmp.source = setdiff(match(sources, x), NA)
        tmp.sink = setdiff(match(sinks, x), NA)
        tmp.mat = A.frag[x, x, drop = FALSE]!=0

        if (length(x)<=1){
            return(NULL)
        }

        if (length(x)==2){
            list(x[c(tmp.source, tmp.sink)])
        }else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1)){
            get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
        } else {
            if (exhaustive){
                lapply(all.paths(A.frag[x,x, drop = FALSE],
                                 source.vertices = tmp.source,
                                 sink.vertices = tmp.sink,
                                 verbose = verbose)$paths,
                       function(y) x[y])
            } else {
                ## ALERT: possible wrong syntax!!!!!!!
                out = do.call('c',
                              lapply(intersect(x, sources),
                                     function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks))
                              )
                out = out[sapply(out, length)!=0]
                if (length(out)>0)
                    out = out[!duplicated(sapply(out, paste, collapse = ','))]
                return(out)
            }
        }
    }, mc.cores = mc.cores))

    if (verbose){
        cat('computed paths\n')
    }

    paths.u = unlist(paths)
    paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
    walks = split(seg[all.frags$subject.id[paths.u]], paths.i)
    values(walks)$seg.id = split(all.frags$subject.id[paths.u], paths.i)

    ## for now just want to pick the non duplicated paths on the original graph and send these to the walk annotation module
    walks = walks[!duplicated(sapply(values(walks)$seg.id, function(x) paste(x, collapse = ',')))]

    ## note aberrant junction ids and filter out trivial walks that don't employ any ab junctions
    A.ab = sparseMatrix(1, 1, x = as.numeric(NA), dims = dim(gg$adj))
    ## ab.ix = !is.na(rowSums(rbind(jab$ab.edges[, 1:2,1])))
    ab.ix = values(gg$junctions)$eclass

    A.ref = sign(gg$adj)-sign(A.ab)

    jdt = data.table(as.data.frame(values(gg$junctions)))
    if (any(ab.ix)){
        A.ab[jdt[,cbind(from1, to1)]] = A.ab[jdt[,cbind(from2, to2)]] = ab.ix
    }

    values(walks)$ab.id = lapply(values(walks)$seg.id, function(x)
        if (length(x)==1) c() else setdiff(A.ab[cbind(x[-length(x)], x[-1])], NA))

    values(walks)$junc.type = lapply(values(walks)$seg.id, function(x)
        if (length(x)==1) c() else sign(A.ref[cbind(x[-length(x)], x[-1])]) + 2*sign(A.ab[cbind(x[-length(x)], x[-1])]))

    walks = walks[!sapply(values(walks)$junc.type, function(x) all(x == 1))]

    values(walks)$seg.id = sapply(values(walks)$seg.id, paste, collapse = ',')
    values(walks)$ab.id = sapply(values(walks)$ab.id, paste, collapse = ',')
    values(walks)$junc.type = NULL

    if (verbose){
        cat(sprintf('Annotating %s walks\n', length(walks)))
    }

    if (length(walks)==0){
        return(walks)
    } else {
        names(walks) = 1:length(walks)
        if (!is.null(query))
            walks = walks[grl.in(walks, query, some = TRUE)]
        ## TODO: return gWalks in the future
        return(annotate.walks(walks,
                              cds,
                              promoters,
                              verbose = verbose,
                              exhaustive = FALSE,
                              mc.cores = mc.cores))
    }
}

## MOMENT
###########################
#' @name proximity
#' @title find new linearly proximal genomic features in light of rearrangements
#' @description
#'
#' Takes a set of n \code{query} elements (\code{GRanges} object, e.g. genes) and determines
#' their proximity to m \code{subject} elements (\code{GRanges} object, e.g. regulatory elements)
#' subject to set of rearrangement adjacencies (\code{junctions}). This redefined genomic
#' distance is computed by \code{gGraph$dist} and is a lower bound when all junctions utilized
#' on the shortest paths are in cis (on the same molecule).
#'
#' @param query GRanges of "intervals of interest" eg regulatory elements
#' @param subject GRanges of "intervals of interest" eg genes
#' @param ra GRangesList of junctions (each a length 2 GRanges, similar to input to karyograph)
#' @param jab existing JaBbA object (overrides ra input)
#' @param verbose logical flag
#' @param mc.cores how many cores (default 1)
#' @param max.dist maximum genomic distance to store and compute (1MB by default) should roughly
#' be the maximum distance at which biological interactions may occur
#'
#' @return
#' list of n x m sparse distance matrices:
#' $ra = subject-query distance in the rearranged genome for all loci < max.dist in tumor genome
#' $wt = subject-query distance in the reference genome for all loci < max.dist in tumor genome
#' $rel = subject-query distance in ra relative to wild type for above loci
#'
#' @details Values x_ij in these matrices should be interpreted with a 1e-9 offset to yield the
#' actual value y_ij, i.e. y_ij = x_ij-1e-9, x_ij>0, y_ij = NA otherwise (allows for sparse
#' encoding of giant matrices)
#'
#' @export
############################################
proximity = function(query,
                     subject,
                     verbose=F,
                     mc.cores=1,
                     max.dist=1e6){
    adj = self$get.adj()
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

    if (!inherits(query, "GRanges") & !inherits(query, "GRanges"))
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

    subject = gr.fix(subject, get(private$segs))
    query = gr.fix(query, get(private$segs))
    gr = c(query, subject)

    ## browser()
    kg = karyograph(ra, gr)
    ## kg2 = gGraph$new()$karyograph(gr, ra)

    ## node.start and node.end delinate the nodes
    ## corresponding to the interval start and end
    ## on both positive and negative tiles of the karyograph
    gr$node.start = gr$node.end = gr$node.start.n = gr$node.end.n = NA;

    ## start and end indices of nodes
    tip = which(as.logical(strand(kg$tile)=='+'))
    tin = which(as.logical(strand(kg$tile)=='-'))
    gr$node.start = tip[gr.match(gr.start(gr,2), gr.start(kg$tile[tip]))]
    gr$node.end = tip[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tip]))]
    gr$node.start.n = tin[gr.match(GenomicRanges::shift(gr.end(gr,2),1), gr.end(kg$tile[tin]))]
    gr$node.end.n = tin[gr.match(gr.start(gr,2), gr.start(kg$tile[tin]))]

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
        if ((ra.which[x, y]) == 1){
            get.shortest.paths(kg$G,
                               vix.query[x, 'end'],
                               vix.subject[y, 'start'],
                               weights = E(kg$G)$weight,
                               mode = 'out')$vpath[[1]]
        }
        else if ((ra.which[x, y]) == 2){
            rev(get.shortest.paths(kg$G,
                                   vix.query[x, 'start'],
                                   vix.subject[y, 'end'],
                                   weights = E(kg$G)$weight,
                                   mode = 'in')$vpath[[1]])
        }
        else if ((ra.which[x, y]) == 3){
            get.shortest.paths(kg$G,
                               vix.query[x, 'end.n'],
                               vix.subject[y, 'start'],
                               weights = E(kg$G)$weight,
                               mode = 'out')$vpath[[1]]
        }
        else if ((ra.which[x, y]) == 4){
            rev(get.shortest.paths(kg$G,
                                   vix.query[x, 'start.n'],
                                   vix.subject[y, 'end'],
                                   weights = E(kg$G)$weight,
                                   mode = 'in')$vpath[[1]])
        }
    }, sum$i, sum$j, SIMPLIFY = F)

                                        #    sum$paths = lapply(sum.paths, function(x) x[-c(1, length(x))])
    sum$paths = sum.paths
    sum$ab.edges = lapply(sum.paths,
                          function(p) setdiff(E(kg$G, path = p)$bp.id, NA))
    return(list(sum = sum,
                rel = rel,
                ra = ra,
                wt = ref,
                G = kg$G,
                G.ref = G.ref,
                tile = kg$tile,
                vix.query = vix.query,
                vix.subject = vix.subject))

}

## ============= Utility functions ============= ##
#' @name rpaths
#' @description
#' Given segments and list of numeric vectors of the indices,
#' return the  reverse complement paths
rpaths = function(segs, paths){
    ## then check the paths
    hb = hydrogenBonds(segs)
    if (hb[, any(is.na(from) | is.na(to))]){
        warning("Some vertices in the paths do not have a reverse complement.")
    }
    hbmap = hb[, setNames(from, to)]

    ## cache the original and rev comp paths
    rpaths = lapply(paths, ## some index doesn't exist anymore!!
                    function(x) {
                        out = rev(hbmap[as.character(x)])
                        names(out) = NULL
                        return(out)
                    })
    return(rpaths)
}

#' @name .gencode_transcript_split
#' @rdname gencode_transcript_split
#' @title .gencode_transcript_split
#'
#' @description
#'
#' splits gencode gr into transcript, taking care of some junky issues in the meantime
#'
#' @author Marcin Imielinski
.gencode_transcript_split = function(gr, tx){
    gr$start.local = gr2dt(gr)[
      , id := 1:length(gr)][
      , tmp.st := 1+c(0, cumsum(width)[-length(width)]), by = transcript_id][
      , keyby = id][
      , tmp.st]

    gr$end.local = gr$start.local + width(gr) -1
    grl = split(gr, gr$transcript_id)
    tmp.val = as.data.frame(tx)[match(names(grl), tx$transcript_id), ]
    rownames(tmp.val) = tmp.val$transcript_id
    names(tmp.val) = capitalize(names(tmp.val))
    values(grl) = tmp.val
    return(grl)
}

###############################################
#' collapse.paths
#'
#' collapse simple paths in a graph G (adjacency matrix or igraph object)
#' returns m x m new adjacency matrix and map of old vertex id's to new ones
#' $adj = m x m matrix
#' #map = length n with indices 1 .. m
#'
###############################################
collapse.paths = function(G, verbose = T){
    if (inherits(G, 'igraph'))
        G = G[,]

    out = G!=0

    if (verbose)
        cat('graph size:', nrow(out), 'nodes\n')

    ## first identify all nodes with exactly one parent and child to do initial collapsing of graph
    singletons = which(Matrix::rowSums(out)==1 & Matrix::colSums(out)==1)

    if (verbose)
        cat('Collapsing simple paths..\n')

    sets = split(1:nrow(G), 1:nrow(G))
    if (length(singletons)>0){
        tmp = out[singletons, singletons]
        cl = igraph::clusters(
            graph(as.numeric(t(Matrix::which(tmp, arr.ind = TRUE))),
                  n = nrow(tmp)), 'weak')$membership
        dix = unique(cl)
        if (length(dix)>0)
        {
            for (j in dix)
            {
                if (verbose)
                    cat('.')

                ## grab nodes in this cluster
                setj = singletons[which(cl == j)]

                ## move all members into a single set
                sets[setj[1]] = list(setj)
                sets[setj[-1]] = list(NULL)

                ## connect this node to the parent and child of the set
                parent = setdiff(which(Matrix::rowSums(out[, setj, drop = FALSE])>0), setj)
                child = setdiff(which(Matrix::colSums(out[setj, , drop = FALSE])>0), setj)
                out[setj, c(setj, child)] = FALSE
                out[c(setj, parent), setj] = FALSE
                out[parent, setj[1]] = TRUE
                out[setj[1], child] = TRUE
            }
        }
    }

    if (verbose){
        cat('done\nnow fixing branches\n')
    }

    todo = rep(FALSE, nrow(G))
    todo[Matrix::rowSums(out)==1 | Matrix::colSums(out)==1] = TRUE

    while (sum(todo)>0){
        sets.last = sets
        out.last = out

        if (verbose)
            if ((sum(todo) %% 200)==0)
                cat('todo:', sum(todo), 'num sets:', sum(!sapply(sets, is.null)), '\n')

        i = which(todo)[1]

        todo[i] = F

        child = which(out[i, ])
        parent = which(out[,i])

        if (length(child)==1 & length(parent)==1){
            ## if there is exactly one child and one parent then we want to merge with one or both
            ## if i-child has no other parents and i-parent has no other child
            ## then merge i, i-parent and i-child
            if (sum(out[,  child])==1 & sum(out[parent, ])==1)
            {
                grandch = which(out[child, ])
                if (length(grandch)>0)
                {
                    out[parent, grandch] = TRUE  ## parent inherits grandchildren of i
                    out[child, grandch] = FALSE
                }
                out[parent, i] = FALSE ## remove node i's edges
                out[i, child] = FALSE
                sets[[parent]] = c(sets[[parent]], sets[[child]], sets[[i]])
                sets[c(i, child)] = list(NULL)
                todo[child] = F ## no longer have to do i-child, since they have already been merged with parent
            }
            ## otherwise if either i-child has no other parent or i-parent has no other children (but not both)
            ## then connect i-parent to i-child, but do not merge them (but merge ONE of them with i)
            else if (sum(out[,  child])==1 | sum(out[parent, ])==1)
            {
                ## if parent has no other children then merge with him
                if (sum(out[parent, ])==1)
                    sets[[parent]] = c(sets[[parent]], sets[[i]])
                else
                    sets[[child]] = c(sets[[child]], sets[[i]])

                out[parent, child] = TRUE
                out[parent, i] = FALSE ## remove node i's edges
                out[i, child] = FALSE
                sets[i] = list(NULL)
            }
        } else if (length(child)==1 & length(parent)>1){
            ## if i has more than one parent but one child, we merge with child if child has no other parents
            if (sum(out[, child])==1)
            {
                sets[[child]] = c(sets[[child]], sets[[i]])
                out[parent, child] = TRUE
                out[parent, i] = FALSE ## remove node i's edges
                out[i, child] = FALSE ## remove node i's edges
                sets[i] = list(NULL)
            }


        } else if (length(child)>1 & length(parent)==1){
            ## if i has more than one child but one parent, then merge with parent if parent has no other children
            if (sum(out[parent, ])==1)
            {
                sets[[parent]] = c(sets[[parent]], sets[[i]])
                out[parent, child] = TRUE
                out[parent, i] = FALSE ## remove node i's edges
                out[i, child] = FALSE ## remove node i's edges
                sets[i] = list(NULL)
            }
        }
    }

    slen = sapply(sets, length)
    ix = which(slen>0)
    map = rep(NA, nrow(G))
    map[unlist(sets)] = match(rep(1:length(sets), slen), ix)
    out = out[ix, ix]
    colnames(out) = rownames(out) = NULL

    return(list(adj = out, map = map, sets = split(1:length(map), map)))
}

#' @name read_gencode
#' @title read_gencode
#'
#' @description
#'
#' reads GENCODE file dump from .rds of GRanges.
#'
#' If that file doesn't exist
#'
#' @param con a path, URL, or \code{GFFFile} object pointing to gene annotation data
#' @param type the value
#'
#' @importFrom rtracklayer import.gff import.gff3
#' @author Marcin Imielinski
#' @export
read_gencode = function(con = Sys.getenv("DEFAULT_GENE_ANNOTATION"),
                        type = NULL,
                        by = NULL){
    ## con is given is some form
    ## path to GTF?
    ## connection to UCSC?
    if (is.character(con) & nchar(con)>0){
        if (file.exists(con)){
            ge = tryCatch(readRDS(con), error=function(e) NULL)
            if (is.null(ge)){
                ge = tryCatch(rtracklayer::import.gff(con), error=function(e) NULL)
                if (is.null(ge)){
                    stop("Given file can't be parsed as RDS or GFF.")
                }
            }
        } else {
            ge = tryCatch(rtracklayer::import.gff(con), error=function(e) NULL)
            if (is.null(ge)){
                stop("Input 'con' is not valid.")
            }
        }
    } else if (inherits(con, "GFFFile")){
        ge = tryCatch(rtracklayer::import.gff(con), error=function(e) NULL)
        if (is.null(ge)){
            stop("Input 'con' is neither string of address or GFFFile object.")
        }
    } else {
        stop("Must provide a valid location of gene annotation data.")
    }

    ## my territory, you can't have "chr"/"Chr" when the default genome doesn't!
    sl = fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]
    if (!any(grepl("chr", names(sl), ignore.case=TRUE))){
        seqlevels(ge) = gsub("chr", "", seqlevels(ge), ignore.case=TRUE)
    }

    TYPES = c('exon', 'gene', 'transcript', 'CDS')
    BY = c('transcript_id', 'gene_id')

    if (!is.null(by))
        by = toupper(by)

    if (!is.null(type))
        type = toupper(type)

    if (!is.null(type))
    {
        type = grep(type, TYPES, value = TRUE, ignore.case = TRUE)[1]
        if (!all(type %in% TYPES))
            stop(sprintf('Type should be in %s', paste(TYPES, collapse = ',')))
        tx = ge[ge$type %in% 'transcript']
        ge = ge[ge$type %in% type]

        if (type == 'CDS' & is.null(by))
            return(.gencode_transcript_split(ge, tx))
    }

    if (!is.null(by))
    {
        by = grep(by, BY, value = TRUE, ignore.case = TRUE)[1]

        if (!(by %in% BY))
            stop(sprintf('Type should be in %s', paste(TYPES, collapse = ',')))

        if (by == 'transcript_id')
            return(.gencode_transcript_split(ge, tx))
        else
            return(split(ge, values(gene)[, by]))
    }
    return(ge)
}

############################
#' @name capitalize
#' @title capitalize
#' @description
#' Capitalize first letter of each character element of vector "string"
#'
#' @param string character vector to capitalize
#' @param un logical flag whether to uncapitalize (=FALSE)
#' @return character vector of strings with capitalized values
##############################
capitalize = function(string, un = FALSE)
{
    if (!un)
    {
        capped <- grep("^[^A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- toupper(substr(string[capped],1, 1))
    }
    else
    {
        capped <- grep("^[A-Z].*$", string, perl = TRUE)
        substr(string[capped], 1, 1) <- tolower(substr(string[capped],1, 1))
    }

    return(string)
}

###############################################
#' annotate.walks
#'
#' Low level function to annotate walks (GRanges list) with cds / promoter annotations
#'
#' given:
#' walks: input grl of walks on the genome
#' tx:  gr annotating transcript boundaries
#' or grl annotating exon level transcripts e.g. refgene or grl with
#' grl-level meta data fields $s1, $s2, $e1, $e2 annotating start and end positions of transcipt and cds
#' respectively, $gene_sym representing gene label, $chr - chromosome, $str strand
#' and gr level features (exon_frame) annotating the frame (0,1,2) of the first exon position in a + transcript
#' and last exon position in a - transcript.
#' Assumes that exons are ordered in each grl item in the order of transcription (i.e. right most exon for negative strand transcripts)
#' (e.g. output of read_refGene(grl = T))
#'
#' @param walks GRangesList of walks to query (eg from traversal of JaBbA object)
#' @param cds  GRangesList of CDS annotation, one per transcript (each a coding region of an exon)
#' @param promoters GRanges of promoters (same length as transcript)
#' @param filter.splice flag whether to filter out splice variants of a given gene
#' @param verbose  flag
#' @param prom.window window to use around each transcript to identify putative promoter if promoter is NULL
#' @param mc.cores number of cores to use
#' @return
#' a grl of putative fusions with annotations in values field:
#' $label
#' $type  e.g. promoter-fusion, 5-UTR fusion, In-frame fusion, 3' truncated fusion, 5' truncated fusion, in-frame poly-fusion, 3' truncated poly-fusion,
#'             5' truncated poly-fusion
#' $genes genes involved in fusion
#' $transcripts transcripts involved in fusion
#'
#' annotates every possible altered transcript in region including
#'   - transcripts with truncated
############################################
annotate.walks = function(walks, cds, promoters = NULL, filter.splice = T, verbose = F, prom.window = 1e3, max.chunk = 1e9, mc.cores = 1, exhaustive = FALSE)
{
    require(igraph)

    if (inherits(walks, 'GRanges'))
        walks = GRangesList(walks)

    if (is(walks, 'list'))
        walks = do.call(GRangesList, walks)

    if (!is(cds, 'GRangesList'))
    {
        if (verbose)
            cat('splitting cds\n')
        cds = .gencode_split(cds, by = 'transcript_id')
    }

    tx.span = seg2gr(values(cds)) ## assumed that transcript span is encoded in the cds metadata (i.e. beginning end including UTR)

    cdsu = gr2dt(grl.unlist(cds)[, c('grl.ix')])
    setkey(cdsu, grl.ix)

    cds.span = cdsu[, list(start = min(start), end = max(end)), keyby = grl.ix][list(1:length(tx.span)), ]

    ## There are negative width CDS in the new GENCODE v27!!!!
    utr.left.dt = gr2dt(tx.span)[
      , list(seqnames = seqnames,
             start = start,
             strand = strand,
             end = cds.span[list(1:length(start)), start],
             transcript_id = Transcript_id,
             transcript_name = Transcript_name,
             gene_name = Gene_name)]

    utr.right.dt = gr2dt(tx.span)[
      , list(seqnames = seqnames,
             start = cds.span[list(1:length(start)), end],
             strand = strand,
             end = end,
             transcript_id = Transcript_id,
             transcript_name = Transcript_name,
             gene_name = Gene_name)]

    ## DO a check of eligibility before converting to GRanges
    ## MOMENT
    trash.ix = which(utr.left.dt[, start>=end] |
                     utr.right.dt[, start>=end])

    if (length(trash.ix)>0){
        message("Throwing out fxxxking trash annotations! Why could there be CDS whose end is smaller than start?")
    }
    utr.left = seg2gr(utr.left.dt[-trash.ix])
    utr.right = seg2gr(utr.right.dt[-trash.ix])

    utr = c(utr.left, utr.right)

    names(values(tx.span))[match(c('Transcript_id', 'Transcript_name', 'Gene_name'), names(values(tx.span)))] = c('transcript_id', 'transcript_name', 'gene_name')

    if (is.null(promoters))
    {
        promoters = flank(tx.span, prom.window)
        values(tx.span) = values(promoters)
    }

    tx.span$type = 'cds'
    promoters$type = 'gene'
    utr.left$type = 'gene'
    utr.right$type = 'gene'
    tx.span$cds.id = 1:length(tx.span)
    tx.span = seg2gr(as.data.table(rrbind(as.data.frame(tx.span), as.data.frame(promoters)))[, list(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1], gene_name = gene_name[1], transcript_id = transcript_id[1], transcript_name = transcript_name[1], cds.id = cds.id), keyby = cds.id][!is.na(cds.id), ], seqlengths = seqlengths(tx.span))

                                        # match up tx.span to walks
    walks.u = grl.unlist(walks)

    ## these are fragments of transcripts that overlap walks
    this.tx.span = gr.findoverlaps(tx.span, walks.u, qcol = c('transcript_id', 'transcript_name', 'gene_name'), verbose = verbose, max.chunk = max.chunk)
    this.tx.span$tx.id = this.tx.span$query.id

    strand(this.tx.span) = strand(tx.span)[this.tx.span$query.id]
    this.tx.span$left.broken = start(this.tx.span) != start(tx.span)[this.tx.span$query.id]
    this.tx.span$right.broken = end(this.tx.span) != end(tx.span)[this.tx.span$query.id]

                                        # remove elements that are "unbroken" by the window boundaries
    this.tx.span = this.tx.span[this.tx.span$left.broken | this.tx.span$right.broken]
    this.tx.span$cds.sign = c('-'= -1, '+' = 1)[as.character(strand(this.tx.span))]
    this.tx.span$window.sign = c('-'= -1, '+' = 1)[as.character(strand(walks.u)[this.tx.span$subject.id])]

                                        # annotate left and right ends (if information available)
                                        # i.e. UTR, CDS, Promoter

    ## we trim cds by 1 nucleotide at both ends so that 'cds' annotation only refers to a cds fragment (not including start and stop)
    ## annotate ends with feature types so that positions internal to cds bases will be 'cds',
    ## external to cds but in cds will 'utr'
    ## and 5' flank (with respect to cds orientation) will be promoter
    this.tx.span$left.feat = NA
    this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
    this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
    this.tx.span$left.feat[gr.findoverlaps(gr.start(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'

    this.tx.span$right.feat = NA
    this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), tx.span, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'cds'
    this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), utr, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'utr'
    this.tx.span$right.feat[gr.findoverlaps(gr.end(this.tx.span), promoters, by = 'transcript_id', max.chunk = max.chunk)$query.id] = 'promoter'

    ## if lands in CDS annotate this.tx.span ends with right and/or left exon frame
    ## (if lands in intron then annotate left end with frame of next exon on right and right end
    ## with frame of next exon on left)
    ## otherwise annotate as NA
    ## we will eventually integrate frames across walks and call a transition "in frame" if the
    ## frame of the right (left) end of the previous + (-) interval

    ## now we want to find the first exon to the right of the left boundary and
    ## the first exon to the left of the right boundary for each fragment
    ##tix = match(this.tx.span$transcript_id, tx.span$transcript_id)

    cds.u = grl.unlist(cds[this.tx.span$tx.id])
                                        #        ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]], resolve.empty = 'start.x'))
    ranges(cds.u) =  ranges(pintersect(cds.u, tx.span[this.tx.span$tx.id[cds.u$grl.ix]]))

    tmp = gr.findoverlaps(this.tx.span, cds.u, scol = c('start.local', 'end.local', 'exon_number'), by = 'transcript_id', verbose = verbose, max.chunk = max.chunk)

    leftmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, start = start(tmp))[, id[which.min(start)], by = qid][, V1]
    rightmost.cds.exon = data.table(id = 1:length(tmp), qid = tmp$query.id, end = end(tmp))[, id[which.max(end)], by = qid][, V1]

    ## now we want to get the frame of the left and right base
    ## of each leftmost and rightmost exon (etc.
    ## this will depend on orientation of the exon and side that we are querying
    ## for left side of - exon, (exonFrame + width) %% 3
    ## for right side of - exon, (exonFrame + width(og.exon) - width %% 3
    ## for left side of + exon, (exonFrame + width(og.exon) - width) %%3
    ## for right side of - exon, (exonFrame + int.exon) %% 3

    leftmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
    (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]),
    (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$start.local[tmp$subject.id[leftmost.cds.exon]]))

    rightmost.coord = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
    (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]),
    (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$start.local[tmp$subject.id[rightmost.cds.exon]]))

    leftmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[leftmost.cds.exon]])=='+'),
    (start(tmp)[leftmost.cds.exon] - start(cds.u)[tmp$subject.id[leftmost.cds.exon]] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3,
    (end(cds.u)[tmp$subject.id[leftmost.cds.exon]] - start(tmp)[leftmost.cds.exon] + cds.u$phase[tmp$subject.id[leftmost.cds.exon]]) %% 3)

    rightmost.frame = ifelse(as.logical(strand(cds.u[tmp$subject.id[rightmost.cds.exon]])=='+'),
    (end(tmp)[rightmost.cds.exon] - start(cds.u)[tmp$subject.id[rightmost.cds.exon]] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3,
    (end(cds.u)[tmp$subject.id[rightmost.cds.exon]] - end(tmp)[rightmost.cds.exon] + cds.u$phase[tmp$subject.id[rightmost.cds.exon]]) %% 3)

    leftmost.frame = leftmost.coord %% 3
    rightmost.frame = rightmost.coord %% 3

    this.tx.span$left.coord = this.tx.span$left.boundary = this.tx.span$right.coord = this.tx.span$right.boundary = NA

    this.tx.span$left.coord[tmp$query.id[leftmost.cds.exon]] = leftmost.coord
    this.tx.span$right.coord[tmp$query.id[rightmost.cds.exon]] = rightmost.coord

    this.tx.span$left.boundary[tmp$query.id[leftmost.cds.exon]] =
        ifelse(strand(this.tx.span)[tmp$query.id[leftmost.cds.exon]]=="+",
               tmp$start.local[leftmost.cds.exon], tmp$end.local[leftmost.cds.exon])

    this.tx.span$right.boundary[tmp$query.id[rightmost.cds.exon]] =
        ifelse(strand(this.tx.span)[tmp$query.id[rightmost.cds.exon]]=="+",
               tmp$end.local[rightmost.cds.exon], tmp$start.local[rightmost.cds.exon])

    this.tx.span$right.exon_del = this.tx.span$left.exon_del = NA;
    this.tx.span$right.frame = this.tx.span$left.frame = NA;
    this.tx.span$right.exon_id= this.tx.span$left.exon_id = NA;

    ## keep track of exon frames to left and right
    this.tx.span$left.frame[tmp$query.id[leftmost.cds.exon]] = leftmost.frame
    this.tx.span$right.frame[tmp$query.id[rightmost.cds.exon]] = rightmost.frame

    this.tx.span$left.exon_id[tmp$query.id[leftmost.cds.exon]] = cds.u[tmp$subject.id[leftmost.cds.exon]]$exon_number
    this.tx.span$right.exon_id[tmp$query.id[rightmost.cds.exon]] = cds.u[tmp$subject.id[rightmost.cds.exon]]$exon_number

    ## keep track of which exons have any sort of sequence deletion
    this.tx.span$left.exon_del[tmp$query.id[leftmost.cds.exon]] =
        start(this.tx.span)[tmp$query.id[leftmost.cds.exon]] > start(cds.u)[tmp$subject.id[leftmost.cds.exon]]
    this.tx.span$right.exon_del[tmp$query.id[rightmost.cds.exon]] =
        end(this.tx.span)[tmp$query.id[rightmost.cds.exon]] <  end(cds.u)[tmp$subject.id[rightmost.cds.exon]]


    ## now traverse each walk and connect "broken transcripts"
    ## each walk is a list of windows, connections will be made with respect to the window walk
    ## paying attention to (1) side of breakage, (2) orientation of adjacent windows in walk (3) orientation of transcript
    ##
    ## right broken + cds upstream of ++ junction attaches to left  broken + cds in next window
    ##              - cds upstream of ++ junction attached to left  broken - cds in next window
    ##              + cds upstream of +- junction attaches to right broken - cds in next window
    ##              - cds upstream of +- junction attaches to right broken + cds in next window
    ##
    ## left  broken + cds upstream of -- junction attaches to right broken + cds in next window
    ##              - cds upstream of -- junction attached to right broken - cds in next window
    ##              + cds upstream of -+ junction attaches to left  broken - cds in next window
    ##              - cds upstream of -+ junction attaches to left  broken + cds in next window
    ##
    ##
    ##

                                        # to achieve this create a graphs on elements of this.tx.span
                                        # connecting elements i and j if a window pair k k+1 in (the corresponding) input grl
                                        # produces a connection with the correct orientation

                                        # match up on the basis of the above factors to determine edges in graph

    ##

    edges = merge(
        data.table(
            i = 1:length(this.tx.span),
            key1 = walks.u$grl.ix[this.tx.span$subject.id],
            key2 = walks.u$grl.iix[this.tx.span$subject.id],
            key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
            key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$right.broken, this.tx.span$left.broken)),
        data.table(
            j = 1:length(this.tx.span),
            key1 = walks.u$grl.ix[this.tx.span$subject.id],
            key2 = walks.u$grl.iix[this.tx.span$subject.id]-1,
            key3 = this.tx.span$cds.sign*this.tx.span$window.sign,
            key4 = ifelse(this.tx.span$window.sign>0, this.tx.span$left.broken, this.tx.span$right.broken)), by = c('key1', 'key2', 'key3', 'key4'), allow.cartesian = TRUE
    )[key4 == TRUE, ]

    ## remove edges that link different transcripts of same gene

    if (filter.splice)
        edges = edges[!(this.tx.span$gene_name[edges$i] == this.tx.span$gene_name[edges$j] & this.tx.span$transcript_id[edges$i] != this.tx.span$transcript_id[edges$j]), ]

    if (nrow(edges)==0)
        return(GRangesList())

    require(Matrix)
    A = sparseMatrix(edges$i, edges$j, x = 1, dims = rep(length(this.tx.span),2))
    sources = which(Matrix::colSums(A!=0)==0)
    sinks = which(Matrix::rowSums(A!=0)==0)

    G = graph.adjacency(A)
    C = igraph::clusters(G, 'weak')
    vL = split(1:nrow(A), C$membership)

    ## collate all paths through this graph
    paths = do.call('c', mclapply(1:length(vL), function(i) {
        if (verbose & (i %% 10)==0)
            cat(i, ' of ', length(vL), '\n')
        x = vL[[i]]
        tmp.source = setdiff(match(sources, x), NA)
        tmp.sink = setdiff(match(sinks, x), NA)
        tmp.mat = A[x, x, drop = FALSE]!=0
        if (length(x)<=1)
            return(NULL)
        if (length(x)==2)
            list(x[c(tmp.source, tmp.sink)])
        else if (all(Matrix::rowSums(tmp.mat)<=1) & all(Matrix::colSums(tmp.mat)<=1))
            get.shortest.paths(G, from = intersect(x, sources), intersect(x, sinks))$vpath
        else
        {
            if (exhaustive)
                lapply(all.paths(A[x,x, drop = FALSE], source.vertices = tmp.source, sink.vertices = tmp.sink, verbose = FALSE)$paths, function(y) x[y])
            else
            {
                out = do.call('c', lapply(intersect(x, sources),
                                          function(x, sinks) suppressWarnings(get.shortest.paths(G, from = x, to = sinks)$vpath), sinks = intersect(x, sinks)))
                out = out[sapply(out, length)!=0]
                if (length(out)>0)
                    out = out[!duplicated(sapply(out, paste, collapse = ','))]
                return(out)
            }
        }
    }, mc.cores = mc.cores))

    fus.sign = this.tx.span$cds.sign * this.tx.span$window.sign
    paths = lapply(paths, function(x) if (fus.sign[x][1]<0) rev(x) else x) ## reverse "backward paths" (i.e. those producing fusions in backward order)
    paths.first = sapply(paths, function(x) x[1])
    paths.last = sapply(paths, function(x) x[length(x)])
    paths.broken.start = ifelse(as.logical(strand(this.tx.span)[paths.first] == '+'), this.tx.span$left.broken[paths.first], this.tx.span$right.broken[paths.first])
    paths.broken.end = ifelse(as.logical(strand(this.tx.span)[paths.last] == '+'), this.tx.span$right.broken[paths.last], this.tx.span$left.broken[paths.last])
    paths.u = unlist(paths)
    paths.i = unlist(lapply(1:length(paths), function(x) rep(x, length(paths[[x]]))))
    tmp.gr = this.tx.span[paths.u]

    ## annotate steps of walk with out of frame vs in frame (if cds is a grl)

    ## left and right exon frame
    paths.u.lec = tmp.gr$left.coord
    paths.u.rec = tmp.gr$right.coord
    paths.u.lef = tmp.gr$left.frame
    paths.u.ref = tmp.gr$right.frame
    paths.u.str = as.character(strand(tmp.gr))
    paths.u.lcds = tmp.gr$left.feat == 'cds'
    paths.u.rcds = tmp.gr$right.feat == 'cds'
    paths.u.lout = tmp.gr$left.feat %in% c('promoter', 'utr')
    paths.u.rout = tmp.gr$right.feat %in% c('promoter', 'utr')

    ## a fragment is in frame either if (1) it begins a walk at frame 0 outside of the cds
    ## or if (2) its frame is concordant with previous cds
    paths.u.inframe = rep(NA, length(paths.u))
    paths.u.cdsend = paths.u.cdsstart = rep(FALSE, length(paths.u)) # this keeps track of cds starts and cds ends in the fusion

    outside = TRUE
    for (i in 1:length(paths.u.inframe))
    {
        if (i == 1)
            outside = TRUE
        else if (paths.i[i] != paths.i[i-1] | (paths.u.str[i] == '+' & paths.u.lout[i]) | (paths.u.str[i] == '-' & paths.u.rout[i]))
            outside = TRUE

        if (outside)
        {
            if (paths.u.str[i] == '+')
                paths.u.inframe[i] = paths.u.lout[i] & paths.u.lef[i] == 0
            else
                paths.u.inframe[i] = paths.u.rout[i] & paths.u.ref[i] == 0

            paths.u.cdsstart[i] = paths.u.inframe[i]
            outside = F
        }
        else
        {
            if (paths.u.str[i] == '+' & paths.u.str[i-1] == '+')
                paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.lef[i] == ((paths.u.ref[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.rcds[i-1]
            else if (paths.u.str[i] == '+' & paths.u.str[i-1] == '-')
                paths.u.inframe[i] = paths.u.lec[i] != 1 & paths.u.ref[i]  == ((paths.u.ref[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.rcds[i-1]
            else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '-')
                paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.ref[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.rcds[i] & paths.u.lcds[i-1]
            else if (paths.u.str[i] == '-' & paths.u.str[i-1] == '+')
                paths.u.inframe[i] = paths.u.rec[i] != 1 & paths.u.lef[i] == ((paths.u.lef[i-1]+1) %% 3) & paths.u.lcds[i] & paths.u.lcds[i-1]
        }

        if ((paths.u.str[i] == '+' & paths.u.rout[i]) | (paths.u.str[i] == '-' & paths.u.lout[i]))
        {
            paths.u.cdsend[i] = paths.u.inframe[i]
            outside = T
        }
    }

    tmp.gr$in.frame = paths.u.inframe;
    tmp.gr$cds.start = paths.u.cdsstart
    tmp.gr$cds.end = paths.u.cdsend
    tmp.gr$del5 = ifelse(paths.u.str == '+', tmp.gr$left.exon_del, tmp.gr$right.exon_del)
    tmp.gr$del3 = ifelse(paths.u.str == '+', tmp.gr$right.exon_del, tmp.gr$left.exon_del)
    tmp.gr$first.coord = ifelse(paths.u.str == '+', tmp.gr$left.coord, tmp.gr$right.coord)
    tmp.gr$last.coord  = ifelse(paths.u.str == '+', tmp.gr$right.coord, tmp.gr$left.coord)
    tmp.gr$first.boundary = ifelse(paths.u.str == '+', tmp.gr$left.boundary, tmp.gr$right.boundary)
    tmp.gr$last.boundary  = ifelse(paths.u.str == '+', tmp.gr$right.boundary, tmp.gr$left.boundary)
    tmp.gr$first.exon = ifelse(paths.u.str == '+', tmp.gr$left.exon_id, tmp.gr$right.exon_id)
    tmp.gr$last.exon = ifelse(paths.u.str == '+', tmp.gr$right.exon_id, tmp.gr$left.exon_id)
    fusions = split(tmp.gr, paths.i)

                                        # now annotate fusions

    values(fusions)$walk.id = data.table(wid = walks.u$grl.ix[tmp.gr$subject.id], fid = paths.i)[, wid[1], keyby = fid][, V1]
                                        #    values(fusions)$walk.id = vaggregate(walks.u$grl.ix[tmp.gr$subject.id], by = list(paths.i), FUN = function(x) x[1])

    tmp.g = tmp.gr$gene_name
    tmp.cds = tmp.gr$transcript_id
    tmp.fe = as.numeric(tmp.gr$first.exon)
    tmp.le = as.numeric(tmp.gr$last.exon)
    tmp.fb = tmp.gr$first.boundary
    tmp.lb = tmp.gr$last.boundary
    tmp.fc = tmp.gr$first.coord
    tmp.lc = tmp.gr$last.coord
    tmp.5d = tmp.gr$del5
    tmp.3d = tmp.gr$del3
    tmp.5d[is.na(tmp.5d)] = FALSE
    tmp.3d[is.na(tmp.3d)] = FALSE
    paths.u.cdsend[is.na(paths.u.cdsend)] = FALSE
    paths.u.cdsstart[is.na(paths.u.cdsstart)] = FALSE
    paths.u.inframe[is.na(paths.u.inframe)] = FALSE

    totpaths = max(paths.i)

    if (verbose)
        cat('Populating coordinates\n')

    values(fusions)[, 'coords'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                           split(gr.string(this.tx.span[paths.u], mb = TRUE, round = 1), paths.i), mc.cores = mc.cores)

    if (verbose)
        cat('Populating transcript names\n')
    values(fusions)[, 'transcript_names'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                     split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                     split(values(tx.span)[, 'transcript_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

    if (verbose)
        cat('Populating transcript ids\n')
    values(fusions)[, 'transcript_ids'] = mcmapply(function(x, y) paste(x, ' (', y, ')', sep = '', collapse = '; '),
                                                   split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i),
                                                   split(values(tx.span)[, 'transcript_id'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

    if (verbose)
        cat('Populating gene names\n')
    values(fusions)[, 'genes'] = mcmapply(function(x) paste(unique(x), collapse = '; '),
                                          split(values(tx.span)[, 'gene_name'][this.tx.span$query.id[paths.u]], paths.i), mc.cores = mc.cores)

    if (verbose)
        cat('Populating alteration\n')
    values(fusions)$alteration =  vaggregate(1:length(paths.i), by = list(paths.i),
                                             FUN = function(x)
                                             {
                                                 if (verbose & (x[1] %% 10)==0)
                                                     cat('Path', unique(paths.i[x]), 'of', totpaths, '\n')
                                                 if (length(unique((tmp.cds[x])))==1) ## single transcript event
                                                 {
                                                     out = NULL
                                                     x = x[!is.na(tmp.fe[x]) & !is.na(tmp.le[x])]
                                                     if (length(x)>0)
                                                     {
                                        #                                                   browser()
                                                         ir = IRanges(pmin(tmp.le[x], tmp.fe[x]), pmax(tmp.fe[x], tmp.le[x]))
                                                         if (length(del <- setdiff(IRanges(min(tmp.fe[x]), max(tmp.le[x])), ir))>0)
                                                         {
                                                             del.fc = pmax(tmp.lc[x[match(start(del)-1, tmp.le[x])]]+1, 1, na.rm = TRUE)
                                                             del.lc = pmin(tmp.fc[x[match(end(del)+1, tmp.fe[x])]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                             out = c(out,## some portion deleted
                                                                     ifelse(start(del)==end(del),
                                                                            paste('deletion of exon ', start(del),
                                                                                  ' [', del.fc, '-', del.lc, 'bp]',
                                                                                  sep = '', collapse = ', '),
                                                                            paste('deletion of exons ', start(del), '-', end(del),
                                                                                  ' [', del.fc, '-', del.lc, 'bp]',
                                                                                  sep = '', collapse = ', ')))
                                                         }

                                                         if (length(amp <- IRanges(coverage(ir)>1))>0)
                                                         {
                                                             amp.fc = tmp.lc[x[match(start(amp), tmp.le[x])]]
                                                             amp.lc = tmp.fc[x[match(end(amp), tmp.fe[x])]]
                                                             out = c(out,   ## some portion duplicated
                                                                     ifelse(start(amp)==end(amp),
                                                                            paste('duplication of exon ', end(amp),
                                                                                  '[', amp.fc, '-', amp.lc, 'bp]',
                                                                                  sep = '', collapse = ', '),
                                                                            paste('duplication of exons ', start(amp), '-', end(amp),
                                                                                  ' [', amp.fc, '-', amp.lc, 'bp]',
                                                                                  sep = '', collapse = ', ')))
                                                         }

                                                         if (any(ix <- tmp.5d[x]))
                                                         {
                                                             out = c(out, paste("partial 5' deletion of exon ", tmp.fe[x[ix]],
                                                                                ' [', tmp.fb[x[ix]], '-', tmp.fc[x[ix]], 'bp]',
                                                                                sep = '', collapse = ', '))  ## some portion duplicated
                                                         }
                                                         if (any(ix <- tmp.3d[x]))
                                                         {
                                                             del.fc = pmax(tmp.lc[x[ix-1]] + 1, 1, na.rm = TRUE)
                                                             del.lc = pmin(tmp.fc[x[ix+1]]-1, max(tmp.lc[x]), na.rm = TRUE)
                                                             out = c(out, paste("partial 3' deletion of exon ", tmp.fe[x[ix]],
                                                                                ' [', tmp.lc[x[ix]], '-', tmp.lb[x[ix]], 'bp]',
                                                                                sep = '', collapse = ', '))  ## some portion duplicated
                                                         }
                                                     }

                                                     if (length(out)>0)
                                                         paste(out, collapse = '; ')
                                                     else
                                                         ''
                                                 }
                                                 else
                                                 {
                                                     return(paste(tmp.g[x], ' ', ifelse(paths.u.cdsstart[x], 'S', ''), ifelse(tmp.5d[x],  'tr', ''),
                                                                  ifelse(is.na(tmp.fe[x]), 'UTR',
                                                                  ifelse(tmp.le[x]==tmp.fe[x],
                                                                         paste('exon ', tmp.le[x], sep = ''),
                                                                         paste('exons ', tmp.fe[x], '-', tmp.le[x], sep = ''))),
                                                                  ' [', tmp.fc[x], '-', tmp.lc[x], 'bp]',
                                                                  ifelse(tmp.3d[x],  'tr', ''), ifelse(paths.u.cdsend[x], 'E', ''), ' ',
                                                                  ifelse(c(paths.u.inframe[x[-1]], FALSE), '-',
                                                                  ifelse((1:length(x))!=length(x), '-X', '')), sep = '', collapse = '-> '))
                                                 }
                                             })

    values(fusions)$max.inframe = vaggregate(paths.u.inframe, by = list(paths.i),
                                             FUN = function(x) return(max(c(0, rle(x)$lengths[which(rle(x)$values == T)]))))
    values(fusions)$num.win = vaggregate(paths.u.inframe, by = list(paths.i), length)
    values(fusions) = cbind(values(walks)[values(fusions)$walk.id, ], values(fusions))

    fusions = fusions[nchar(values(fusions)$alteration)>0, ]
    return(fusions)
}

#########################################
#' @name get.tile.id
#' Create the name vector of a skew-symmetric node set
#'
#' @param x a strand-specific \code{GRanges}, where both strand of the same range must present
#'
#' @return a \code{numeric} vector of the same length, where a pair of opposite values
#' indicate the two strands of the same range
#'
#' @details
#'
#' @export
########################################
get.tile.id = function(segs){
    if (!inherits(segs, "GRanges")){
        stop("Only takes GRanges as input for now.")
    }
    hb = hydrogenBonds(segs = segs)
    if (hb[, any(is.na(from) | is.na(to))]){
        stop("Not fully strand paired.")
    }
    hb[, pix := ifelse(from <= to, paste(from, to), paste(to, from))]
    hb[, tile.id := as.numeric(as.factor(pix))]
    hb[, tile.id := ifelse(as.logical(strand(segs)[from]=="+"), tile.id, -tile.id)]
    setkey(hb, "from")
    return(hb[.(seq_along(segs)), tile.id])
}

#########################################
#' @name ul
#' Upper left corner of a matrix
#'
#' @param x a \code{matrix} or \code{Matrix} object
#' @param n the number of rows and cols to show
#'
#' @return the top left corner matrix
#' @export
########################################
ul = function(x, n=6){
    n = pmin(min(dim(x)), n)
    if (n==0) {
        return(NULL)
    }
    return(x[1:n, 1:n])
}

############################################
#' @name jab2json
#' @title jab2json
#'
#' @description
#'
#' Dumps JaBbA graph into json
#'
#' @param jab input jab object
#' @param filename output json file
#' @author Marcin Imielinski
###########################################
jab2json = function(jab,
                    filename,
                    maxcn = 100,
                    maxweight = 100){
    ## ++ = RL
    ## +- = RR
    ## -+ = LL
    ## -- = LL
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

    writeLines(out, filename)
}

########################################################
#' @name gr2json
#' @title gr2json
#'
#' @description
#'
#' Dumps GRanges into JSON with metadata features as data points in  "intervals"
#'
#'
#' @param GRange input jab object
#' @param filename output json file
#' @author Marcin Imielinski
########################################################
gr2json = function(intervals,
                   filename,
                   y = rep("null", length(intervals)),
                   labels = '',
                   maxcn = 100,
                   maxweight = 100){
    ## ++ = RL
    ## +- = RR
    ## -+ = LL
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
    if (is.null(node.dt$type)){
        node.dt$type = 'interval'
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

    writeLines(out, filename)
    return(out)
}

#############################
#' @name levapply
#' @title levapply
#'
#' @description
#' Applies FUN locally to levels of x and returns vector of length()
#' (eg can do a "local" order within levels)
#'
#' @param x input vector of data
#' @param by length(x) vector of categorical labels
#' @param FUN function that takes a length k vector and outputs a length k vector, used for processing each "level" of by
#' @return length(x) vector of outputs, the results of applying FUN to each "by" defined level of x
#' @author Marcin Imielinski
#############################
levapply = function(x,
                    by,
                    FUN = 'order'){
    if (!is.list(by)){
        by = list(by)
    }

    f = factor(do.call('paste', c(list(sep = '|'), by)))
    ixl = split(1:length(x), f);
    ixv = lapply(ixl, function(y) x[y])
    res = structure(unlist(lapply(ixv, FUN)), names = unlist(ixl))
    out = rep(NA, length(x))
    out[as.numeric(names(res))] = res;
    return(out)
}



#############################
#' get.ploidy
#' We define ploidy as the width-weighted mean of copy number. In other words, how many copies
#' of a unique set of genomic ranges are in the input.
#'
#' @param segs a \code{GRanges} object holding \code{gGraph} node data with copy number annotated
#'
#' @return \code{numeric} scalar of the ploidy value
#'
#' @export
#############################
get.ploidy = function(segs){
    if (!inherits(segs, "GRanges")){
        stop("Error: Not a GRanges!")
        return(NULL)
    }
    ## NOTE: doesn't have to be disjoint
    ## if (!isDisjoint(segs)) {
    ##     warning("Must be disjoint!")
    ##     segs = gr.disjoin(segs)
    ## }

    ## MARCIN COMMENT: WHAT IF THERE IS TWO COLUMNS HERE MATCHING CN???
    if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case=T))==0){
        warning("No copy number (cn) column!")
        return(NULL)
    }

    ## MARCIN COMMENT: WHAT IF THERE IS TWO COLUMNS HERE MATCHING
    cn = mcols(segs)[, cnix[1]]
    wd = width(segs)
    good.ix = which(!is.na(cn))

    pl = weighted.mean(cn[good.ix], wd[good.ix], na.rm=T)
    return(pl)
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
#'
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
#'
#' @author Marcin Imielinski
######################################################
mmatch = function(A, B, dir = 1)
{
    SEP = ' ';
    Atxt = apply(A, dir, function(x) paste(x, collapse = SEP))
    Btxt = apply(B, dir, function(x) paste(x, collapse = SEP))

    return(match(Atxt, Btxt))
}

############################################
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
###########################################
alpha = function(col, alpha)
{
    col.rgb = col2rgb(col)
    out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
}

##########################################
#' @name e2j
#' @title given node and edges set, return GRangesList of the junction locations
#' @description
#' Any adjacency between DNA segments can be defined as a pair of locations on the reference genome,
#' each with the orientation specifying the fused side of the molecule. This data is easily represented
#' by a \code{GRangesList}.
#'
#' @param segs the GRanges of nodes
#' @param es the data.table of edges
#' @param etype \code{character} scalar, specifying the type of edges to convert
#'
#' @details
#' \code{etype} can be any of the following values, "all", "aberrant", "reference", "loose". It can also be
#' regular expressions that match these values.
#'
#' @return \code{junctions} object specifying the breakpoint pair and orientations of selected edges
#' @export
###########################################
e2j = function(segs, es, etype="aberrant"){
    if (verbose <- getOption("gGnome.verbose")){
        message("Return the junctions based on edges in this graph.")
    }

    segs = copy(segs)
    es = copy(es)

    strmap = setNames(c("+", "-"), c("-", "+"))

    if (!is.element("type", colnames(es)) |
        !is.element("loose", colnames(segs))){
        tmp = etype(segs, es, force=T, both=TRUE)
        es = tmp$es
        segs = tmp$segs
    }

    es[, eix := 1:.N]

    if (etype=="all"){
        etype = c("aberrant", "reference", "loose")
    }

    abe = es[grepl(pattern = etype, type)]
    if (nrow(abe)==0){
        empty.out = junctions()
        return(empty.out)
    }

    if (any(! c("fromChr", "fromStr", "fromStart", "fromEnd",
                "toChr", "toStr", "toStart", "toEnd") %in%
            colnames(abe))){
        if (verbose){
            message("Redo the important metadata gathering.")
        }

        abe[, fromStr := ":="(fromChr = as.vector(seqnames(segs[from])),
                              fromStr = as.vector(strand(segs[from])),
                              fromStart = start(segs[from]),
                              fromEnd = end(segs[from]),
                              toChr = as.vector(seqnames(segs[to])),
                              toStr = as.vector(strand(segs[to])),
                              toStart = start(segs[to]),
                              toEnd = end(segs[to]))]
    }

    if (any(!c("eid", "reid") %in% colnames(abe))){
        hb = hydrogenBonds(segs)
        hb.map = hb[, setNames(from, to)]
        abe[, ":="(eid = paste(from, to),
                   reid = paste(hb.map[as.character(to)],
                                hb.map[as.character(from)]))]
    }

    abe[, ":="(ix = 1:.N,
               rix = match(reid, eid))]
    abe[, unique.ix := ifelse(rix>=ix,
                              paste(ix, rix),
                              paste(rix, ix))]
    abe[, eclass := as.numeric(as.factor(unique.ix))]
    abe[, iix := 1:.N, by=eclass]
    setkeyv(abe, c("eclass", "iix"))

    jdt = abe[iix==1, .(eclass, from, to,
                        fromChr, fromStr, fromStart, fromEnd,
                        toChr, toStr, toStart, toEnd)]

    bp1 = dt2gr(jdt[, .(seqnames = fromChr,
                        strand = strmap[fromStr],
                        start = ifelse(fromStr=="+", fromEnd, fromStart-1),
                        end = ifelse(fromStr=="+", fromEnd, fromStart-1),
                        eclass)])

    bp2 = dt2gr(jdt[, .(seqnames = toChr,
                        strand = toStr,
                        start = ifelse(toStr=="+", toStart-1, toEnd),
                        end = ifelse(toStr=="+", toStart-1, toEnd),
                        eclass)])

    bps = gr.fix(GRangesList(bp1, bp2), segs)
    junc = junctions(grl.pivot(bps))
    values(junc)$eclass = bp1$eclass
    values(junc)$type = abe[.(values(junc)$eclass, 1), type]
    if ("cn" %in% colnames(abe)){
        values(junc)$cn = abe[iix==1][.(values(junc)$eclass), cn]
    }
    values(junc)$from1 = abe[.(values(junc)$eclass, 1), from]
    values(junc)$to1 = abe[.(values(junc)$eclass, 1), to]
    values(junc)$from2 = abe[.(values(junc)$eclass, 2), from]
    values(junc)$to2 = abe[.(values(junc)$eclass, 2), to]

    return(junc)
}

##########################################
#' @name etype
#' @title infer edge type based on node coordinates and orientation
#'
#' @param segs the GRanges of nodes
#' @param es the data.table of edges
#' @param force logical, whether to overwrite existing node or edge type
#' @param both logical, if TRUE, return a list of updated segs and es
#'
#' @return es with "type" column.
#' @export
###########################################
etype = function(segs, es, force=FALSE, both=FALSE){
    if (is.null(segs)){
        warning("Empty nodes.")
        return(NULL)
    }
    if (is.null(es)){
        es = data.table(from = numeric(0),
                        to = numeric(0))
    }
    if (!inherits(segs, "GRanges")){
        stop("Error:segs must be GRanges")
    }
    if (!inherits(es, "data.frame")){
        stop("Error:es must be data.frame")
    }
    if (!all(c("from", "to") %in% colnames(es))){
        stop("Error: 'from' & 'to' must be in es!")
    }
    if (!inherits(es, "data.table")){
        es = data.table::as.data.table(es)
    }

    if ("type" %in% colnames(es) & nrow(es)>0){
        if (all(es$type %in% c("reference", "aberrant", "loose")) & force==FALSE){
            return(es)
        }
    }

    es2 = copy(es)
    if (nrow(es)==0){
        segs$loose = FALSE
        es2$type = character(0)
    } else {
        ## as definition of graph, segs and es must be both sets
        ## so no redundancy

        es2[, type := "unknown"]

        ## then dedup edges
        es2[, ":="(eid = paste(from, to))]
        if (es2[, anyDuplicated(eid)]){
            if ("cn" %in% colnames(es2)){
                ## if edges comes with CN field, add them together
                es2[, .(from, to, cn=sum(cn), type), by=eid]
            }
            else {
                ## otherwise just ignore it
                es2 = es2[!duplicated(eid),]
            }
        }

        ## Start to determine edge types
        es2[, ":="(fromChr = as.vector(seqnames(segs[from])),
                   fromStr = as.vector(strand(segs[from])),
                   fromStart = start(segs[from]),
                   fromEnd = end(segs[from]),
                   toChr = as.vector(seqnames(segs[to])),
                   toStr = as.vector(strand(segs[to])),
                   toStart = start(segs[to]),
                   toEnd = end(segs[to]))]

        if (!"loose" %in% colnames(values(segs)) | force){
            ## a loose end is a degree 1, width 1 node
            ## that intersect the incident end of its only neighbor
            ## segs$terminal = !seq_along(segs) %in% es2[, from] | !seq_along(segs) %in% es2[, to]
            segs$out.degree = as.vector(es2[, table(from)][as.character(seq_along(segs))])
            segs$out.degree[which(is.na(segs$out.degree))] = 0
            segs$in.degree = as.vector(es2[, table(to)][as.character(seq_along(segs))])
            segs$in.degree[which(is.na(segs$in.degree))] = 0
            segs$degree = segs$in.degree + segs$out.degree

            ## find the loose ends with 1 out degree
            which.loose.src = which(width(segs)==1 &
                                    segs$degree==1 &
                                    segs$out.degree==1)
            ## find the loose ends with 1 in degree
            which.loose.sink = which(width(segs)==1 &
                                     segs$degree==1 &
                                     segs$in.degree==1)

            ## find their partners
            term.map = c(es2[from %in% which.loose.src, setNames(to, from)], es2[to %in% which.loose.sink, setNames(from, to)])

            ## compare source loose ends to its partner's 5' end
            loose.src.candidates = segs[which.loose.src]
            loose.src.partners = segs[term.map[as.character(which.loose.src)]]
            final.loose.src = which(loose.src.candidates[, c()] == gr.start(loose.src.partners, ignore.strand=FALSE)[, c()])

            ## compare sink loose ends to its partner's 3' end
            loose.sink.candidates = segs[which.loose.sink]
            loose.sink.partners = segs[term.map[as.character(which.loose.sink)]]
            final.loose.sink = which(loose.sink.candidates[, c()] == gr.end(loose.sink.partners, ignore.strand=FALSE)[, c()])

            which.loose = c(which.loose.src[final.loose.src],
                            which.loose.sink[final.loose.sink])
            segs$loose = seq_along(segs) %in% which.loose
        }

        ## if we know who are the loose ends:
        which.loose = setNames(segs$loose,
                               as.character(seq_along(segs)))
        es2[, ":="(fromLoose = which.loose[as.character(from)],
                   toLoose = which.loose[as.character(to)])]
        ## any edge involves a loose node must be loose
        es2[fromLoose==T | toLoose==T, type:="loose"]

        ## interchr or interstrand edges must be aberrant
        es2[(fromChr!=toChr | fromStr!=toStr), type := "aberrant"]

        ## identify reference edges
        eps=1e-9
        es2[fromChr==toChr & fromStr==toStr &
            fromStr=="+" & abs(toStart-fromEnd-1) < eps,
            type := "reference"]
        es2[fromChr==toChr & fromStr==toStr &
            fromStr=="-" & abs(fromStart-toEnd-1) < eps,
            type := "reference"]

        ## the rest is same strand jumping events
        ## "deletion bridges"
        es2[type=="unknown", type := "aberrant"]
    }

    if (both==TRUE) {
        return(list(segs = segs, es = es2))
    } else {
        return(es2)
    }
}

#####################################
#' @name setxor
#' @title XOR operation on sets
#'
#' @param A set A
#' @param B set B
#'
#' @return The set of elements belong to either A or B, but not both.
#' @author Marcin Imielinski
####################################s
setxor = function (A, B)
{
    return(setdiff(union(A, B), intersect(A, B)))
}

############################################
#' @name write.tab
#' @title wrapper around write.table
#' @author Marcin Imielinski
###########################################
write.tab = function (x, ..., sep = "\t", quote = F, row.names = F)
{
    if (!is.data.frame(x)){
        x = as.data.frame(x)
    }
    write.table(x, ..., sep = sep, quote = quote, row.names = row.names)
}

################################
#' @name dedup
#' @title dedup
#'
#' @description
#' relabels duplicates in a character vector with .1, .2, .3
#' (where "." can be replaced by any user specified suffix)
#'
#' @param x input vector to dedup
#' @param suffix suffix separator to use before adding integer for dups in x
#' @return length(x) vector of input + suffix separator + integer for dups and no suffix for "originals"
#' @author Marcin Imielinski
################################
dedup = function(x, suffix = '.')
{
    dup = duplicated(x);
    udup = setdiff(unique(x[dup]), NA)
    udup.ix = lapply(udup, function(y) which(x==y))
    udup.suffices = lapply(udup.ix, function(y) c('', paste(suffix, 2:length(y), sep = '')))
    out = x;
    out[unlist(udup.ix)] = paste(out[unlist(udup.ix)], unlist(udup.suffices), sep = '');
    return(out)
}

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
#' @author Marcin Imielinski
#'
#############################################################
munlist = function(x, force.rbind = F, force.cbind = F, force.list = F)
{
    if (!any(c(force.list, force.cbind, force.rbind)))
    {
        if (any(sapply(x, function(y) is.null(dim(y))))){
            force.list = T
        }
        if (length(unique(sapply(x, function(y) dim(y)[2]))) == 1){
            force.rbind = T
        }
        if ((length(unique(sapply(x, function(y) dim(y)[1]))) == 1)){
            force.cbind = T
        }
    }
    else{
        force.list = T
    }

    if (force.list){
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, length(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (length(x[[y]])>0) 1:length(x[[y]]) else NULL)),
                     unlist(x)))
    }
    else if (force.rbind){
        return(cbind(ix = unlist(lapply(1:length(x), function(y) rep(y, nrow(x[[y]])))),
                     iix = unlist(lapply(1:length(x), function(y) if (nrow(x[[y]])>0) 1:nrow(x[[y]]) else NULL)),
                     do.call('rbind', x)))
    }
    else if (force.cbind){
        return(t(rbind(ix = unlist(lapply(1:length(x), function(y) rep(y, ncol(x[[y]])))),
                       iix = unlist(lapply(1:length(x), function(y) if (ncol(x[[y]])>0) 1:ncol(x[[y]]) else NULL)),
                       do.call('cbind', x))))
    }
}

################################################
#' read_vcf: utility function to read VCF into GRanges object
#'
#' @name read_vcf
#' @importFrom VariantAnnotation readVcf
#'
###############################################
read_vcf = function (fn, gr = NULL, hg = "hg19", geno = NULL, swap.header = NULL,
                     verbose = FALSE, add.path = FALSE, tmp.dir = "~/temp/.tmpvcf",
                     ...)
{
    in.fn = fn
    if (verbose)
        cat("Loading", fn, "\n")
    if (!is.null(gr)) {
        tmp.slice.fn = paste(tmp.dir, "/vcf_tmp", gsub("0\\.",
                                                       "", as.character(runif(1))), ".vcf", sep = "")
        cmd = sprintf("bcftools view %s %s > %s", fn, paste(gr.string(gr.stripstrand(gr)),
                                                            collapse = " "), tmp.slice.fn)
        if (verbose)
            cat("Running", cmd, "\n")
        system(cmd)
        fn = tmp.slice.fn
    }
    if (!is.null(swap.header)) {
        if (!file.exists(swap.header))
            stop(sprintf("Swap header file %s does not exist\n",
                         swap.header))
        system(paste("mkdir -p", tmp.dir))
        tmp.name = paste(tmp.dir, "/vcf_tmp", gsub("0\\.", "",
                                                   as.character(runif(1))), ".vcf", sep = "")
        if (grepl("gz$", fn))
            system(sprintf("zcat %s | grep '^[^#]' > %s.body",
                           fn, tmp.name))
        else system(sprintf("grep '^[^#]' %s > %s.body", fn,
                            tmp.name))
        if (grepl("gz$", swap.header))
            system(sprintf("zcat %s | grep '^[#]' > %s.header",
                           swap.header, tmp.name))
        else system(sprintf("grep '^[#]' %s > %s.header", swap.header,
                            tmp.name))
        system(sprintf("cat %s.header %s.body > %s", tmp.name,
                       tmp.name, tmp.name))
        vcf = readVcf(tmp.name, hg, ...)
        system(sprintf("rm %s %s.body %s.header", tmp.name, tmp.name,
                       tmp.name))
    }
    else vcf = readVcf(fn, hg, ...)
    out = granges(vcf)
    if (!is.null(values(out)))
        values(out) = cbind(values(out), info(vcf))
    else values(out) = info(vcf)
    if (!is.null(geno)) {
        if (geno)
            for (g in names(geno(vcf))) {
                geno = names(geno(vcf))
                warning(sprintf("Loading all geno field:\n\t%s",
                                paste(geno, collapse = ",")))
            }
        gt = NULL
        for (g in geno) {
            m = as.data.frame(geno(vcf)[[g]])
            names(m) = paste(g, names(m), sep = "_")
            if (is.null(gt))
                gt = m
            else gt = cbind(gt, m)
        }
        values(out) = cbind(values(out), as(gt, "DataFrame"))
    }
    if (!is.null(gr))
        system(paste("rm", tmp.slice.fn))
    if (add.path)
        values(out)$path = in.fn
    return(out)
}

############################################################
#' ra_breaks: parse junction data from various common formats
#'
#' @name ra_breaks
#'
#' @description Parsing various formats of structural variation data into junctions.
#'
#' @usage ra_breaks(rafile,
#' keep.features = T,
#' seqlengths = hg_seqlengths(),
#' chr.convert = T,
#' geno=NULL,
#' flipstrand = FALSE,
#' swap.header = NULL,
#' breakpointer = FALSE,
#' seqlevels = NULL,
#' force.bnd = FALSE,
#' skip = NA)
#'
#' @param rafile path to the junctions file. See details for the compatible formats.
#' @param keep.features \code{logical}, if TRUE preserve meta data from the input
#' @param seqlengths a named \code{numeric} vector containing reference contig lengths
#' @param chr.convert \code{logical}, if TRUE strip "chr" prefix from contig names
#' @param geno \code{logical}, whether to parse the 'geno' fields of VCF
#' @param flipstrand \code{logical}, if TRUE will flip breakpoint strand
#' @param swap.header path to the alternative VCF header file
#' @param breakpointer \code{logical}, if TRUE will parse as breakpointer output
#' @param seqlevels vector for renaming the chromosomes
#' @param force.bnd if TRUE overwrite all junction "type" to "BND"
#' @param skip \code{numeric} lines to skip
#'
#' @details
#' A junction is a unordered pair of strand-specific genomic locations (breakpoints). Within a given
#' reference genome coordinate system, we call the direction in which coordinates increase "+". A breakpoint
#' is a width 1 (\code{start==end})genomic range with \code{strand} specified, and "+" means the side with larger
#' coordinate is fused with the other breakpoint in a junction.
#'
#' \code{rafile} must be one of the following formats:
#' 1) Some VCF (variant call format). We currently support the VCF output
#' from a number of structural variation detection methods, namely
#' SvABA (https://github.com/walaj/svaba),
#' DELLY (https://github.com/dellytools/delly),
#' LUMPY (https://github.com/arq5x/lumpy-sv),
#' novoBreak (https://sourceforge.net/projects/novobreak/). In theory,
#' VCF defined with BND style should be compatible but be cautious
#' when using the output from other methods since
#' no universal data definition is adopted by the community yet.
#' 2) BEDPE (http://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
#' 3) Textual output from Breakpointer
#' (http://archive.broadinstitute.org/cancer/cga/breakpointer)
#' 4) R serialized object storing junctions (.rds)
#'
#' @section Warning:
#' We assume the orientation definition in the input is consistent with ours. Check with
#' the documentation of your respective method to make sure. If the contrary, use
#' \code{flipstrand=TRUE} to reconcile.
#'
#' @return a \code{GRangesList} of the junctions
#'
#' @importFrom VariantAnnotation readVcf
#' @import data.table
#'
#' @export
###########################################################
ra_breaks = function(rafile,
                     keep.features = T,
                     seqlengths = hg_seqlengths(),
                     chr.convert = T,
                     geno=NULL,
                     flipstrand = FALSE,
                     swap.header = NULL,
                     breakpointer = FALSE,
                     seqlevels = NULL,
                     force.bnd = FALSE,
                     skip = NA,
                     get.loose = FALSE){
    ## if TRUE will return a list with fields $junctions and $loose.ends
    if (is.character(rafile))
    {
        if (grepl('.rds$', rafile)){
            ra = readRDS(rafile)
            ## validity check written for "junctions" class
            return(junctions(ra))
        }
        else if (grepl('(.bedpe$)', rafile)){
            ra.path = rafile
            cols = c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'name', 'score', 'str1', 'str2')

            ln = readLines(ra.path)
            if (is.na(skip))
            {
                nh = min(c(Inf, which(!grepl('^((#)|(chrom))', ln))))-1
                if (is.infinite(nh)){
                    nh = 1
                }
            }
            else{
                nh = skip
            }

            if ((length(ln)-nh)==0){
                ## if (get.loose){
                ##     return(list(junctions = GRangesList(GRanges(seqlengths = seqlengths))[c()], loose.ends = GRanges(seqlengths = seqlengths)))
                ## }
                ## else{
                return(GRangesList(GRanges(seqlengths = seqlengths))[c()])
                ## }
            }

            if (nh ==0){
                rafile = fread(rafile, header = FALSE)
            }
            else
            {

                rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh), error = function(e) NULL)
                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = '\t'), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    rafile = tryCatch(fread(ra.path, header = FALSE, skip = nh, sep = ','), error = function(e) NULL)
                }

                if (is.null(rafile)){
                    stop('Error reading bedpe')
                }
            }
            ## this is not robust enough! there might be mismatching colnames
            setnames(rafile, 1:length(cols), cols)
            rafile[, str1 := ifelse(str1 %in% c('+', '-'), str1, '*')]
            rafile[, str2 := ifelse(str2 %in% c('+', '-'), str2, '*')]
        }
        else if (grepl('(vcf$)|(vcf.gz$)', rafile)){
            require(VariantAnnotation)
            vcf = readVcf(rafile, Seqinfo(seqnames = names(seqlengths), seqlengths = seqlengths))

            ## vgr = rowData(vcf) ## parse BND format
            vgr = read_vcf(rafile, swap.header = swap.header, geno=geno)
            mc = data.table(as.data.frame(mcols(vgr)))

            if (!('SVTYPE' %in% colnames(mc))) {
                warning('Vcf not in proper format.  Is this a rearrangement vcf?')
                return(GRangesList());
            }

            if (any(w.0 <- (width(vgr)<1))){
                warning("Some breakpoint width==0.")
                ## right bound smaller coor
                ## and there's no negative width GR allowed
                vgr[which(w.0)] = gr.start(vgr[which(w.0)]) %-% 1
            }

            ## BND format doesn't have duplicated rownames
            if (any(duplicated(names(vgr)))) names(vgr) = NULL

            ## no events
            if (length(vgr) == 0)
                return (GRangesList())

            ## local function that turns old VCF to BND
            .vcf2bnd = function(vgr){
                if (!"END" %in% colnames(values(vgr)))
                    stop("Non BND SV should have the second breakpoint coor in END columns!")

                if (!"CHR2" %in% colnames(values(vgr)) | any(is.na(vgr$CHR2)))
                    vgr$CHR2 = as.character(seqnames(vgr))

                bp2 = data.table(as.data.frame(mcols(vgr)))
                bp2[, ":="(seqnames=CHR2, start=as.numeric(END), end=as.numeric(END))]
                bp2.gr = dt2gr(bp2)
                mcols(bp2.gr) = mcols(vgr)

                if (!is.null(names(vgr)) & !anyDuplicated(names(vgr))){
                    jid = names(vgr)
                } else {
                    jid = seq_along(vgr)
                }
                names(vgr) = paste(paste0("exp", jid), "1", sep=":")
                names(bp2.gr) = paste(paste0("exp", jid), "2", sep=":")

                vgr=resize(c(vgr, bp2.gr), 1)

                if (all(grepl("[_:][12]$",names(vgr)))){
                    ## row naming same with Snowman
                    nm <- vgr$MATEID <- names(vgr)
                    ix <- grepl("1$",nm)
                    vgr$MATEID[ix] = gsub("(.*?)(1)$", "\\12", nm[ix])
                    vgr$MATEID[!ix] = gsub("(.*?)(2)$", "\\11", nm[!ix])
                    vgr$SVTYPE="BND"
                }
                return(vgr)
            }

            ## TODO: Delly and Novobreak
            ## fix mateids if not included
            if (!"MATEID" %in% colnames(mcols(vgr))) {
                ## TODO: don't assume every row is a different junction
                ## Novobreak, I'm looking at you.
                ## now delly...
                ## if SVTYPE is BND but no MATEID, don't pretend to be
                if (length(fake.bix <- which(values(vgr)$SVTYPE=="BND"))!=0){
                    values(vgr[fake.bix])$SVTYPE = "TRA"
                }

                ## add row names just like Snowman
                if (all(names(vgr)=="N" | ## Novobreak
                        is.null(names(vgr)) |
                        all(grepl("^DEL|DUP|INV|BND", names(vgr)))) ## Delly
                    ){
                    ## otherwise if all "N", as Novobreak
                    ## or starts with DEL|DUP|INV|BND, as Delly
                    ## expand and match MATEID
                    vgr=.vcf2bnd(vgr)
                }
            } else if (any(is.na(mid <- as.character(vgr$MATEID)))){
                ## like Lumpy, the BND rows are real BND but blended with non-BND rows
                ## treat them separately
                if (is.null(vgr$CHR2)){
                    vgr$CHR2 = as.character(NA)
                }

                names(vgr) = gsub("_", ":", names(vgr))
                vgr$MATEID = sapply(vgr$MATEID, function(x) gsub("_", ":", x))

                values(vgr) = data.table(as.data.frame(values(vgr)))

                ## break up the two junctions in one INV line!
                if ("STRANDS" %in% colnames(mc) & any(ns <- sapply(vgr$STRANDS, length)>1)){
                    ## first fix format errors, two strand given, but not comma separeted
                    ## so you'd have taken them as single
                    if (any(fuix <- sapply(vgr[which(!ns)]$STRANDS, str_count, ":")>1)){
                        which(!ns)[fuix] -> tofix
                        vgr$STRANDS[tofix] = lapply(vgr$STRANDS[tofix],
                                                    function(x){
                                                        strsplit(gsub("(\\d)([\\+\\-])", "\\1,\\2", x), ",")[[1]]
                                                    })
                        ns[tofix] = TRUE
                    }

                    ## for the one line two junction cases
                    ## split into two lines
                    vgr.double = vgr[which(ns)]
                    j1 = j2 = vgr.double
                    st1 = lapply(vgr.double$STRANDS, function(x)x[1])
                    st2 = lapply(vgr.double$STRANDS, function(x)x[2])
                    j1$STRANDS = st1
                    j2$STRANDS = st2
                    vgr.double = c(j1, j2)
                    names(vgr.double) = dedup(names(vgr.double))
                    vgr = c(vgr[which(!ns)], vgr.double)
                }

                mid <- as.logical(sapply(vgr$MATEID, length))
                vgr.bnd = vgr[which(mid)]
                vgr.nonbnd = vgr[which(!mid)]

                vgr.nonbnd = .vcf2bnd(vgr.nonbnd)

                mc.bnd = data.table(as.data.frame(values(vgr.bnd)))
                mc.nonbnd = data.table(as.data.frame(values(vgr.nonbnd)))
                mc.bnd$MATEID = as.character(mc.bnd$MATEID)

                vgr = c(vgr.bnd[,c()], vgr.nonbnd[,c()])
                values(vgr) = rbind(mc.bnd, mc.nonbnd)
            }

            ## sanity check
            if (!any(c("MATEID", "SVTYPE") %in% colnames(mcols(vgr))))
                stop("MATEID or SVTYPE not included. Required")

            vgr$mateid = vgr$MATEID
            ## what's this???
            vgr$svtype = vgr$SVTYPE

            if (!is.null(info(vcf)$SCTG))
                vgr$SCTG = info(vcf)$SCTG

            if (force.bnd)
                vgr$svtype = "BND"

            if (sum(vgr$svtype == 'BND')==0)
                warning('Vcf not in proper format.  Will treat rearrangements as if in BND format')

            if (!all(vgr$svtype == 'BND')){
                warning(sprintf('%s rows of vcf do not have svtype BND, treat them as non-BND!',
                                sum(vgr$svtype != 'BND')))

            }

            bix = which(vgr$svtype == "BND")
            vgr = vgr[bix]
            alt <- sapply(vgr$ALT, function(x) x[1])

            ## Determine each junction's orientation
            if ("CT" %in% colnames(mcols(vgr))){
                message("CT INFO field found.")
                if ("SVLEN" %in% colnames(values(vgr))){
                    ## proceed as Novobreak
                    ## ALERT: overwrite its orientation!!!!
                    del.ix = which(vgr$SVTYPE=="DEL")
                    dup.ix = which(vgr$SVTYPE=="DUP")
                    vgr$CT[del.ix] = "3to5"
                    vgr$CT[dup.ix] = "5to3"
                }

                ## also, Delly is like this
                ori = strsplit(vgr$CT, "to")
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                orimap = setNames(c("+", "-"), c("5", "3"))
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd
                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            }
            else if ("STRANDS" %in% colnames(mcols(vgr))){
                ## TODO!!!!!!!!!!!!!!!
                ## sort by name, record bp1 or bp2
                message("STRANDS INFO field found.")
                iid = sapply(strsplit(names(vgr), ":"), function(x)as.numeric(x[2]))
                vgr$iid = iid
                vgr = vgr[order(names(vgr))]
                iid = vgr$iid

                ## get orientations
                ori = strsplit(substr(unlist(vgr$STRANDS), 1, 2), character(0))
                orimap = setNames(c("+", "-"), c("-", "+"))

                ## map strands
                strd = orimap[sapply(seq_along(ori), function(i) ori[[i]][iid[i]])]
                strand(vgr) = strd

                vgr.pair1 = vgr[which(iid==1)]
                vgr.pair2 = vgr[which(iid==2)]
            }
            else if (any(grepl("\\[", alt))){
                message("ALT field format like BND")
                ## proceed as Snowman
                vgr$first = !grepl('^(\\]|\\[)', alt) ## ? is this row the "first breakend" in the ALT string (i.e. does the ALT string not begin with a bracket)
                vgr$right = grepl('\\[', alt) ## ? are the (sharp ends) of the brackets facing right or left
                vgr$coord = as.character(paste(seqnames(vgr), ':', start(vgr), sep = ''))
                vgr$mcoord = as.character(gsub('.*(\\[|\\])(.*\\:.*)(\\[|\\]).*', '\\2', alt))
                vgr$mcoord = gsub('chr', '', vgr$mcoord)

                ## add extra genotype fields to vgr
                geno(vcf)
                values(vgr)
                if (all(is.na(vgr$mateid)))
                    if (!is.null(names(vgr)) & !any(duplicated(names(vgr))))
                    {
                        warning('MATEID tag missing, guessing BND partner by parsing names of vgr')
                        vgr$mateid = paste(gsub('::\\d$', '', names(vgr)),
                        (sapply(strsplit(names(vgr), '\\:\\:'), function(x) as.numeric(x[length(x)])))%%2 + 1, sep = '::')
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
                    else{
                    stop('Error: MATEID tag missing')
                }

                vgr$mix = as.numeric(match(vgr$mateid, names(vgr)))

                pix = which(!is.na(vgr$mix))

                vgr.pair = vgr[pix]

                if (length(vgr.pair)==0){
                    stop('Error: No mates found despite nonzero number of BND rows in VCF')
                }

                vgr.pair$mix = match(vgr.pair$mix, pix)

                vix = which(1:length(vgr.pair)<vgr.pair$mix)
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
            }

            ra = grl.pivot(GRangesList(vgr.pair1[, c()], vgr.pair2[, c()]))

            ## ALERT: vgr has already been subsetted to only include BND rows
            ## bix is the original indices, so NOT compatible!
            ## this.inf = values(vgr)[bix[pix[vix]], ]
            if (exists("pix") & exists("vix")) this.inf = values(vgr)[pix[vix], ]
            if (exists("iid")) this.inf = values(vgr[which(iid==1)])

            if (is.null(this.inf$POS)){
                this.inf = cbind(data.frame(POS = ''), this.inf)
            }
            if (is.null(this.inf$CHROM)){
                this.inf = cbind(data.frame(CHROM = ''), this.inf)
            }

            if (is.null(this.inf$MATL)){
                this.inf = cbind(data.frame(MALT = ''), this.inf)
            }

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

            if (is.null(values(ra)$TIER)){
                ## baseline tiering of PASS vs non PASS variants
                ## ALERT: mind the naming convention by diff programs
                ## TODO: make sure it is compatible with Delly, Novobreak, Meerkat
                ## Snowman/SvABA uses "PASS"
                ## Lumpy/Speedseq uses "."
                values(ra)$tier = ifelse(values(ra)$FILTER %in% c(".", "PASS"), 2, 3)
            }
            else {
                values(ra)$tier = values(ra)$TIER
            }

            ## ra = ra.dedup(ra)

            if (!get.loose | is.null(vgr$mix)){
                return(ra)
            }
            else
            {
                npix = is.na(vgr$mix)
                vgr.loose = vgr[npix, c()] ## these are possible "loose ends" that we will add to the segmentation

                ## NOT SURE WHY BROKEN
                tmp =  tryCatch( values(vgr)[bix[npix], ],
                                error = function(e) NULL)
                if (!is.null(tmp)){
                    values(vgr.loose) = tmp
                }
                else{
                    values(vgr.loose) = cbind(vcf@fixed[bix[npix], ], info(vcf)[bix[npix], ])
                }

                return(list(junctions = ra, loose.ends = vgr.loose))
            }
        }
        else{
            rafile = read.delim(rafile)
        }
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

    if (flipstrand) ## flip breaks so that they are pointing away from junction
    {
        rafile$str1 = ifelse(rafile$strand1 == '+', '-', '+')
        rafile$str2 = ifelse(rafile$strand2 == '+', '-', '+')
    }

    if (!is.null(seqlevels)) ## convert seqlevels from notation in tab delim file to actual
    {
        rafile$chr1 = seqlevels[rafile$chr1]
        rafile$chr2 = seqlevels[rafile$chr2]
    }


    if (is.null(rafile$str1)){
        rafile$str1 = rafile$strand1
    }

    if (is.null(rafile$str2)){
        rafile$str2 = rafile$strand2
    }

    if (!is.null(rafile$pos1) & !is.null(rafile$pos2))
    {
        if (breakpointer)
        {
            rafile$pos1 = rafile$T_BPpos1
            rafile$pos2 = rafile$T_BPpos2
        }

        if (!is.numeric(rafile$pos1)){
            rafile$pos1 = as.numeric(rafile$pos1)
        }

        if (!is.numeric(rafile$pos2)){
            rafile$pos2 = as.numeric(rafile$pos2)
        }

        ## clean the parenthesis from the string

        rafile$str1 <- gsub('[()]', '', rafile$str1)
        rafile$str2 <- gsub('[()]', '', rafile$str2)

        ## goal is to make the ends point <away> from the junction where - is left and + is right
        if (is.character(rafile$str1) | is.factor(rafile$str1)){
            rafile$str1 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str1))))
        }

        if (is.character(rafile$str2) | is.factor(rafile$str2)){
            rafile$str2 = gsub('0', '-', gsub('1', '+', gsub('\\-', '1', gsub('\\+', '0', rafile$str2))))
        }


        if (is.numeric(rafile$str1)){
            rafile$str1 = ifelse(rafile$str1>0, '+', '-')
        }

        if (is.numeric(rafile$str2)){
            rafile$str2 = ifelse(rafile$str2>0, '+', '-')
        }

        rafile$rowid = 1:nrow(rafile)

        bad.ix = is.na(rafile$chr1) | is.na(rafile$chr2) | is.na(rafile$pos1) | is.na(rafile$pos2) | is.na(rafile$str1) | is.na(rafile$str2) | rafile$str1 == '*'| rafile$str2 == '*' | rafile$pos1<0 | rafile$pos2<0

        rafile = rafile[which(!bad.ix), ]

        if (nrow(rafile)==0){
            return(GRanges())
        }

        seg = rbind(data.frame(chr = rafile$chr1, pos1 = rafile$pos1, pos2 = rafile$pos1, strand = rafile$str1, ra.index = rafile$rowid, ra.which = 1, stringsAsFactors = F),
                    data.frame(chr = rafile$chr2, pos1 = rafile$pos2, pos2 = rafile$pos2, strand = rafile$str2, ra.index = rafile$rowid, ra.which = 2, stringsAsFactors = F))

        if (chr.convert){
            seg$chr = gsub('chr', '', gsub('25', 'M', gsub('24', 'Y', gsub('23', 'X', seg$chr))))
        }

        out = seg2gr(seg, seqlengths = seqlengths)[, c('ra.index', 'ra.which')];
        out = split(out, out$ra.index)
    }
    else if (!is.null(rafile$start1) & !is.null(rafile$start2) & !is.null(rafile$end1) & !is.null(rafile$end2))
                         {
                             ra1 = gr.flipstrand(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
                             ra2 = gr.flipstrand(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
                             out = grl.pivot(GRangesList(ra1, ra2))
                         }



    if (keep.features){
        values(out) = rafile[, ]
    }

    ## if (!is.null(pad)){
    ##     out = ra.dedup(out, pad = pad)
    ## }

    if (!get.loose){
        return(out)
    }
    else{
        return(list(junctions = out, loose.ends = GRanges()))
    }

    return(new("junctions", out))
}



##########################
#' @name chr2num
#' @title chr2num
#' @description
#'
#' Convert from chrXX to numeric format
#'
#' @param x factor, Rle or character vector with chromosome names
#' @param xy Flag to convert M to 25, Y to 24 and X to 23. Default FALSE
#' @return character vector with xy=FALSE, or numeric vector with xy=TRUE
##########################
chr2num = function(x, xy = FALSE)
{
    if (inherits(x, 'factor') | inherits(x, 'Rle')){
        x = as.character(x)
    }

    out = gsub('chr', '', x);

    if (!xy){
        out = as.numeric(gsub('M', '25', gsub('Y', '24', gsub('X', '23', out))))
    }

    return(out)
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
#' @param label.edges boolean Flag for etc.
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
############################################
karyograph = function(junctions, ## this is a grl of breakpoint pairs (eg output of ra_breaks(dranger.df) where dranger is df of dranger output)
                      tile = NULL, ## pre-existing set of intervals on top of which to build a graph (eg endpoints from a copy number based segmentation)
                      label.edges = FALSE)
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

        if (any(as.logical(strand(bp1) == '*') | as.logical(strand(bp2) == '*'))){
            stop('Error: bp1 and bp2 must be signed intervals (i.e. either + or -)')
        }

        if (length(bp1) != length(bp2)){
            stop('Error: bp1 and bp2 inputs must have identical lengths')
        }

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
    else{
        if (is.null(tile)){
            tile = si2gr(junctions)
            if (length(tile)==0){
                warning('Empty input given, producing empty output')
                return(NULL)
            }
            A = sparseMatrix(1,1, x = 0, dims = rep(length(tile), 2))
            return(list(tile = tile, adj = A, G = graph.adjacency(A), ab.adj = A != 0, ab.edges = NULL, junctions = junctions))
        }

        junctions = GRangesList()
        bp1 = bp2 = GRanges()
    }

    if (!is.null(tile)){

        ## find disjoint union of tile and join with gaps
        tile = gr.fix(tile)
        tile = gr.fix(tile, bp1)
        bp1 = gr.fix(bp1, tile) ## argh argh argh .. more pain avoiding hacks
        strand(tile) = '+'
        tile = disjoin(tile)
        tile = sort(c(tile, gaps(tile)))

        ## make sure seqlevels / seqinfo are identical
        if (!identical(sort(seqlevels(tile)), seqlevels(junctions))){
            tile = gr.fix(tile, junctions)
            junctions = gr.fix(junctions, tile)
        }

        if(length(junctions)>0){
            tbp = setdiff(gr.stripstrand(gr.trim(tile, 1)), gr.stripstrand(grbind(bp1, bp2)))
            bp1 = gr.fix(bp1, tbp)
            bp2 = gr.fix(bp2, tbp) ## seqlengths pain
            tbp = gr.fix(tbp, bp1)
        }
        else{
            tbp = gr.stripstrand(gr.trim(tile, 1))
        }

        tbp = tbp[start(tbp)!=1]

        if (length(tbp)>0){
            tbp$seg.bp = TRUE
        }
    }
    else{
        tbp = NULL;
    }

    if (length(junctions)>0){

        if (length(tbp)>0){
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()], tbp[, c()]))))
        }
        else{
            g = gaps(gr.stripstrand(sort(c(bp1[, c()], bp2[, c()]))))
        }
    }
    else{
        g = gaps(gr.stripstrand(sort(tbp)));
    }

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
    if (length(junc.bp)>0){
        junc.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(junc.bp), start(junc.bp)))
    }

    ## make sure all seqlenths are compatible (so effing annoying)
    tile = gr.fix(tile, bp1)
    tile = gr.fix(tile, bp2)
    bp1 = gr.fix(bp1, tile)
    bp2 = gr.fix(bp2, tile)


    ## also keep track of tbp associatd bp.ix
    all.bp = grbind(bp1, bp2, tbp)
    all.bpix = numeric()

    if (length(all.bp)>0){
        all.bpix = which(paste(seqnames(tile), end(tile)) %in% paste(seqnames(all.bp), start(all.bp)))
    }

    ## now to build the graph, we would like to fuse all the bp associated intervals with their previous interval
    ## UNLESS they are preceded by another bp associated interval
    ##
    if (length(all.bpix>0)){
        to.fuse = all.bpix[which(all.bpix>1 & !((all.bpix-1) %in% all.bpix))]
        end(tile)[to.fuse-1] = end(tile)[to.fuse-1]+1
        tile = tile[-to.fuse]
    }

    if (length(junc.bpix)>0){

        ## we have a partition of genomic segments flanked by tile endpoints and/or ra junctions
        ##
        ## Input junction syntax is interpreted as follows:
        ## a- b+ junctions connect seg ending with position a to seg starting with b+1
        ## a- b- junctions connect seg ending with position a to seg ending with position b (on neg strand)
        ## a+ b+ junctions connect seg starting with position a+1 (on negative strand) to seg starting with position b+1
        ## a+ b- junctions connect seg starting with position a+1 (on negative strand) to seg ending with position b (on neg strand)

        ## collect all pairwise adjacencies implied by breakpoints
        ## eg imagine a|bp1|b
        ##            c|bp2|d
        ## "+" bp point to the right (eg b or d), "-" bp point to the left (a or c)

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

        ## clean up adj pairs
        ## remove any that have crossed a chromosome boundary from their breakpoint
        ## this will occur in cases of badly formed breakpoint input (eg breakpoints that point outward
        ## from their telomeres)
        edge.id = rep(1:nrow(ab.pairs), 2)
        ab.pairs = rbind(ab.pairs, cbind(-ab.pairs[,2], -ab.pairs[,1]));
        ab.pairs.bpid = c(ab.pairs.bpid, ab.pairs.bpid)

        ## build "aberrant" adjacency matrix representing directed graph of edges connecting
        ## <signed> nodes.
        ## note: indices of matrix represent edge labels
        adj.ab = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
                        dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
        tmp.ix = cbind(match(as.character(ab.pairs[,1]), rownames(adj.ab)),
                       match(as.character(ab.pairs[,2]), colnames(adj.ab)))
        adj.ab[tmp.ix[!duplicated(tmp.ix), , drop = F]] = ab.pairs.bpid[!duplicated(tmp.ix)]
    }
    else{
        ab.pairs.bpid = edge.id = c()
        ab.pairs = matrix(nrow = 0, ncol = 2);
        adj.ab = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                        dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
    }

    ## build reference adjacency matrix (representing consecutive segments on the reference genome)
    ## note: indices of matrix represent edge labels
    seg.ix = 1:length(tile)
    ref.pairs = cbind(seg.ix[1:(length(seg.ix)-1)], seg.ix[2:(length(seg.ix))])
                                        # ref.pairs = ref.pairs[ref.pairs[,1]>0 & ref.pairs[,2]!=length(tile), ]
    ref.pairs = ref.pairs[which(as.character(seqnames(tile[ref.pairs[,1]])) == as.character(seqnames(tile[ref.pairs[,2]]))), ]

    if (nrow(ref.pairs)>0){
        edge.id = c(edge.id, max(edge.id) + rep(1:nrow(ref.pairs), 2))
        ref.pairs = rbind(ref.pairs, cbind(-ref.pairs[,2], -ref.pairs[,1])) # reverse ref pairs
        adj.ref = Matrix(0, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
        adj.ref[cbind(match(as.character(ref.pairs[,1]), rownames(adj.ref)),
                      match(as.character(ref.pairs[,2]), colnames(adj.ref)))] = nrow(ab.pairs)+1:nrow(ref.pairs)
    }
    else{
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

    if (length(ab.pairs.bpid)>0){
        E(G)$bp.id[ab.ix] = ab.pairs.bpid[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
    }
    E(G)$eid = NA; ## what is edge ID??? how is different from edge.ix?
    E(G)$eid[ab.ix] = edge.id[adj.ab[cbind(E(G)$from[ab.ix], E(G)$to[ab.ix])]]
    E(G)$eid[!ab.ix] = edge.id[adj.ref[cbind(E(G)$from[!ab.ix], E(G)$to[!ab.ix])]]
    values(tile) = values(tile)[, c('tile.id', 'is.tel')]
    tile$ab.source = 1:length(tile) %in% E(G)$from[ab.ix]
    tile$ab.target = 1:length(tile) %in% E(G)$to[ab.ix]

    ## important: map input ra to aberrant graph edges, i.e. ab.edges matrix with $from $to and $edge.ix columns
    ## and one row for each aberrant edge
    ab.edges = array(NA, dim = c(length(junctions), 3, 2), dimnames = list(NULL, c('from', 'to', 'edge.ix'), c('+', '-')))
    dupped = duplicated(ab.pairs.bpid)
    ab.edges[,1:2,1] = cbind(match(ab.pairs[!dupped,1], names(tile)), match(ab.pairs[!dupped,2], names(tile)))
    ab.edges[,1:2,2] = cbind(match(ab.pairs[dupped,1], names(tile)), match(ab.pairs[dupped,2], names(tile)))
    ab.edges[,3, 1] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to)) ## must be easier way to perform this taks
    ab.edges[,3, 2] = match(paste(ab.edges[,1,1], '|', ab.edges[,2,1]), paste(E(G)$from, '|', E(G)$to))

    if (label.edges & nrow(ab.edges)>0){
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
#' and length |E| vector of edge copy numbers (eclass), length |E| vector of edge equivalence
#' classes (both outputs of jbaMIP.process) and computes most likely karyotypes that fit the
#' edge copy number profile subject to some prior likelihood over the k extreme paths
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
###############################################################
karyoMIP = function(K,
                    e,
                    eclass = 1:length(e),
                    kclass = NULL,
                    prior = rep(0, ncol(K)),
                    cpenalty = 1,
                    tilim = 100,
                    epgap = 1,
                    nsolutions = 50,
                    objsense = 'max',
                    gurobi = TRUE,
                    cplex = !gurobi,
                    ...){
    M = 1e7;
    K = as(K, 'sparseMatrix')

    if (length(prior)!=ncol(K)){
        stop('Error: prior must be of the same length as number of columns in K')
    }
    v.ix = 1:ncol(K)
    M.ix = max(v.ix) + (1:ncol(K))
    n = max(M.ix);

    ## add big M constraints
    ## upper bound is infinity if indicator is positive
    Zero = sparseMatrix(1, 1, x = 0, dims = c(n, n))
    Amub = Zero[1:length(M.ix), ]
    Amub[cbind(1:length(M.ix), v.ix)] = 1
    Amub[cbind(1:length(M.ix), M.ix)] = -M

    ## lower bound > 0 if indicator is positive
    Amlb = Zero[1:length(M.ix), ]
    Amlb[cbind(1:length(M.ix), v.ix)] = 1
    Amlb[cbind(1:length(M.ix), M.ix)] = -0.1

    if (is.null(kclass)){
        kclass = .e2class(K, eclass)
    }

    kclass.counts = table(kclass)
    ## any equiv i.e. strand flipped contig pairs? then make sure they appear
    ## in solutions togethrer
    if (any(kclass.counts>1))
    {
        bikclass = which(kclass.counts>1)
        Ac = Zero[1:length(bikclass), , drop=F]
        pairs = matrix(unlist(split(1:length(kclass), kclass)[as.character(bikclass)]), ncol = 2, byrow = T)
        Ac[cbind(1:nrow(pairs), pairs[,1])] = 1
        Ac[cbind(1:nrow(pairs), pairs[,2])] = -1
    }
    else{
        Ac = Zero[1,,drop = FALSE]
        ## combine constraints
    }
    A = rBind(cBind(K, Zero[rep(1, nrow(K)), M.ix]), Amub, Amlb, Ac);
    b = c(e, rep(0, nrow(Amlb)*2), rep(0, nrow(Ac)));
    sense = c(rep('E', nrow(K)), rep('L', nrow(Amlb)), rep('G', nrow(Amlb)), rep('E', nrow(Ac)))
    vtype = c(rep('I', length(v.ix)), rep('B', length(M.ix)))

    cvec = c(rep(0, length(v.ix)), prior-cpenalty*rep(1, length(M.ix)))

    if (cplex){

        sol = Rcplex::Rcplex(cvec = cvec,
                             Amat = A,
                             bvec = b,
                             sense = sense,
                             Qmat = NULL,
                             lb = 0,
                             ub = Inf,
                             n = nsolutions,
                             objsense = objsense,
                             vtype = vtype,
                             control = c(list(...), list(tilim = tilim, epgap = epgap)))

        if (!is.null(sol$xopt)){
            sol = list(sol)
        }

        sol = lapply(sol, function(x)
        {
            x$kcn = round(x$xopt[v.ix])
            x$kclass = kclass
            x$mval= round(x$xopt[M.ix])
            return(x)
        })


    } else {

        sense.map = setNames(c("=", "<=", ">="),
                             c("E", "L", "G"))
        model = list(A = A,
                     obj = cvec,
                     sense = setNames(sense.map[sense],NULL),
                     rhs = b,
                     ## lb is 0 by default
                     ## ub is Inf by default
                     vtype = vtype,
                     modelsense = objsense)
        params = list()
        sol = gurobi::gurobi(model, params)

        if (names(sol)[1]=="status"){
            sol$kcn = round(sol$x[v.ix])
            sol$kclass = kclass
            sol$mval = round(sol$x[M.ix])
        }
    }
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
##############################################################
karyoMIP.to.path = function(sol,
                            K,
                            e,
                            gr = NULL,
                            mc.cores = 1,
                            verbose = TRUE){
    contigs = which(sol$kcn!=0)
    c1 =  contigs[!duplicated(sol$kclass[contigs])]
    c2 = setdiff(contigs, c1)
    c2 = c2[match(sol$kclass[c2], sol$kclass[c1])]
    contigs = c1
    contigs2 = c2

    nm.gr = names(gr)
    names(gr) = NULL

    if (is.null(nm.gr)){
        nm.gr  = 1:length(gr)
    }

    if (any(duplicated(nm.gr))){
        nm.gr = 1:length(gr)
    }

    if (!is.character(e)){
        e = matrix(as.character(e), ncol = 2)
    }

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
    out$paths = mclapply(1:length(contigs), function(i){
        if (verbose)
            cat('contig', i, 'of', length(contigs), '\n')

        k = K[, i]
        v.all = setdiff(as.vector(e[k!=0,]), NA)
        ## v.all = rownames(B)[which(rowSums(abs(B) %*% k)>0)]  ## vertices associated with edges in path / cycle  k

        ## this is a slack to slack path involving 1 node
        if (length(v.all)==1){
            return(v.all)
        }

        ## make subgraph corresponding to edges in this path / cycle
        ##       B.tmp = B[, which(!is.slack)[k[!is.slack]!=0], drop = F] ##
        ##       so = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x<0))]
        ##       si = rownames(B.tmp)[apply(B.tmp, 2, function(x) which(x>0))]
        ##       sG = graph(rbind(so, si))
        ##       sG = graph(rbind(so, si))

        tmp.e = e[k!=0, ,drop = F]
        tmp.e = tmp.e[rowSums(is.na(tmp.e))==0,,drop = F]
        sG = graph(t(tmp.e))

        if (out$is.cyc[i]){
            p.fwd = names(get.shortest.paths(sG, v.all[1], v.all[pmin(length(v.all), 2)])$vpath[[1]])
            p.bwd = names(get.shortest.paths(sG, v.all[pmin(length(v.all), 2)], v.all[1])$vpath[[1]])
            return(unique(unlist(c(p.fwd, p.bwd))))
        }
        else{
            io = as.numeric(B[, !is.slack, drop = F] %*% k[!is.slack])
            v.in = rownames(B)[io<0][1]
            v.out = rownames(B)[io>0][1]
            return(names(get.shortest.paths(sG, v.in, v.out)$vpath[[1]]))
        }
    }, mc.cores = mc.cores)

    if (!is.null(gr)){
        if (is.null(nm.gr)){
            nm.gr = names(B)
        }
        names(gr) = NULL
        out$grl = do.call('GRangesList', lapply(out$paths, function(x) gr[match(x, nm.gr), c()]))  ## match non-slack vertices
        names(out$grl) = paste('Contig ', out$kix, ' (CN = ', out$cn, ')', sep = '')
        values(out$grl)$is.cycle = out$is.cyc
    }

    return(out)
}

####################################################################
#' jbaMIP.process
#'
#' process jbaMIP solution "sol" given original graph "g" (karyograph() list output)
#' into JaBbA object
#'
#' output is
#'
#' @param sol JaBbA object
#' @param allelic logical flag specifying whether object is allelic
#' @return
#' list with items:
#' $B incidence matrix of augmented graph (including slack vertices) (vertices x edges)
#' rownames of $B are vertex names of $G and colnames of B are named with character version of their $G indices
#' (i.e. column order of B  respects the original edge order in the solution)
#'
#' $e edge constraints for downstream karyoMIP, i.e the copy numbers at the edges
#' $e.ij numedges x 2 vertex pair matrix denoting what are the vertex pairs corresponding to the cols of $B and entries of $e, $eclass, $etype etc
#' $eclass id for each unique edge / anti-edge equivalence class
#' $etype specifies whether edge is slack or nonslack
###################################################################
jbaMIP.process = function(sol,
                          allelic = FALSE){
    if (allelic){
        sol = list(segstats = sol$asegstats, adj = sol$aadj)
    }

    if (!all(c('segstats', 'adj') %in% names(sol))){
        stop('sol must be output of jbaMIP()')
    }

    if (is.null(sol$segstats$tile.id)){
        stop('sol$segstats must be populated with tile.id')
    }
    else{
        if (!all(table(sol$segstats$tile.id)==2)){
            stop('sol$segstats$tile.id are malformed, there should be exactly two instances of each tile.id in sol$segstats, one for the positive and one for the negative strand of the same interval')
        }

        tmp = lapply(split(1:length(sol$segstats$tile.id), sol$segstats$tile.id), rev)

        recip.ix = rep(NA, length(sol$segstats))
        recip.ix[order(sol$segstats$tile.id)] = unlist(tmp)
    }

    if (is.null(sol$segstats$eslack.in)){
        sol$segstats$eslack.in = sol$segstats$slack.in
    }

    if (is.null(sol$segstats$eslack.out)){
        sol$segstats$eslack.out = sol$segstats$slack.out
    }

    ed.ij = which(sol$adj!=0, arr.ind = T)

    ## B is vertices x edges (i.e. signed incidence matrix)
    B = sparseMatrix(c(ed.ij[,1], ed.ij[,2]),
                     rep(1:nrow(ed.ij), 2),
                     x = rep(c(-1.00001, 1),
                             each = nrow(ed.ij)),
                     dims = c(nrow(sol$adj), nrow(ed.ij)))

    rownames(B) = 1:nrow(B)

    tmp.ix = which(abs(B)>=1)
    B[tmp.ix] = round(B[tmp.ix]) ## "0.00001" hack to take care of eclass matching below, these are length 1 self loop edge cases

    ix.tel.5 = which(Matrix::colSums(sol$adj!=0)==0)  ## make fake slacks for telomeres
    sol$segstats$eslack.in[ix.tel.5] = sol$segstats$cn[ix.tel.5]

    ix.tel.3 = which(Matrix::rowSums(sol$adj!=0)==0)
    sol$segstats$eslack.out[ix.tel.3] = sol$segstats$cn[ix.tel.3]  ## make fake slacks for telomeres

    ix.eslack.out = which(sol$segstats$eslack.out!=0);
    names(ix.eslack.out) = paste('out slack', ix.eslack.out)
    ix.eslack.in = which(sol$segstats$eslack.in!=0);
    names(ix.eslack.in) = paste('in slack', ix.eslack.in)

    names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3] = paste(names(ix.eslack.in)[ix.eslack.in %in% ix.tel.3], 'tel')
    names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5] = paste(names(ix.eslack.out)[ix.eslack.out %in% ix.tel.5], 'tel')

    ## we add "slack edges" and "slack nodes" to incidence matrix
    Zero = sparseMatrix(1, 1, x = 0, dims = c(length(ix.eslack.in) + length(ix.eslack.out), ncol(B)))

    if (nrow(Zero)>0){
        rownames(Zero) = c(paste('slack in', 1:length(ix.eslack.in)), paste('slack out', 1:length(ix.eslack.out)))
    }

    Bs = rBind(B, Zero)
    ed.ij = rbind(ed.ij, cbind(ix.eslack.out, NA), cbind(NA, ix.eslack.in))

    Is = Diagonal(n = nrow(Bs), rep(1, nrow(Bs)))

    Bs = cBind(Bs, -Is[, ix.eslack.out], Is[, ix.eslack.in])
    colnames(Bs) = c(as.character(1:ncol(B)), names(ix.eslack.out), names(ix.eslack.in))

    ## map new "slack nodes" to their reciprocals
    recip.ix = c(recip.ix,
                 nrow(B) + length(ix.eslack.out) +  match(recip.ix[ix.eslack.out], ix.eslack.in),
                 nrow(B) + match(recip.ix[ix.eslack.in], ix.eslack.out)
                 )

    ## match matrix against its reverse complement (i.e. rotation) to find reciprocal edges
    erecip.ix = mmatch(t(Bs), t(-Bs[recip.ix, ])) ## maps edges to their reciprocals

    tmp.na = which(is.na(erecip.ix))
    ## fix the self loops so that they match
    if (length(tmp.na)>0){
        erecip.ix[tmp.na] = tmp.na[mmatch(t(Bs[1:nrow(Bs), tmp.na]), t(Bs[recip.ix,tmp.na, drop = F]))]
    }

    ## now use this mapping to define edge equivalence classes
    rmat = t(apply(cbind(erecip.ix, erecip.ix[erecip.ix]), 1, sort)) ## length(erecip.ix) x 2 matrix of edge ids and their reciprocal, sorted

    ## eclass will map length(erecip.ix) edges to length(erecip.ix)/2 edge equivalence class ids
    eclass = mmatch(rmat, rmat[!duplicated(rmat), ])

    Bs = round(Bs) ## remove the 0.0001 dummy coefficients (i.e. for self loops)

    ## e will store observed copy states corresponding to edges (i.e. columns of Bs)
    e = c(sol$adj[which(sol$adj!=0)], sol$segstats$eslack.out[ix.eslack.out],  sol$segstats$eslack.in[ix.eslack.in])

    return(list(e = e, e.ij = ed.ij, B = Bs, eclass = eclass, etype = c(ifelse(grepl('slack', colnames(Bs)), 'slack', 'nonslack'))))
}



################################
#' @name seg.fill
#' Supplement the other strand if missing from input.
#'
#' @export
################################
seg.fill = function(segs, verbose=FALSE){
    if (length(segs)==0){
        return(segs)
    }
    segs = segs[!duplicated(segs)]

    ## collapse +/- strand
    ss = unique(gr.stripstrand(segs))
    idss = match(gr.stripstrand(segs), ss)
    if (all(table(idss)==2)){
        ## if the variety matches, check the copy numbers
        ## TODO: back to this later
    }
    else {
        if (verbose){
            warning("Some segments do not have matching strand.")
        }
        fill = tofill = segs[which(idss %in% which(table(idss)==1))]
        strmap = setNames(c("+", "-"), c("-", "+"))
        strand(fill) = strmap[as.character(strand(tofill))]
        mcols(fill) = mcols(tofill)
        segs = c(segs, fill)
    }
    return(segs)
}



## #' @name hbonds
## #' @title get the hydrogen bonds between elements that form a dsDNA
## #' @param gr a \code{GRanges} object
## #' @return a \code{data.table} with columns "from", "to", and "type"
## hbonds = function(gr){
##     if (!inherits(gr, "GRanges")){
##         gr = GRanges(gr)
##     }
##     if (!is.element("prop", slotNames(gr))){
##         return(hydrogenBonds(gr))
##     } else {
##         ## do hb in a general way
##         ## need
##     }
## }

## FOR gGnome v2.0:
## TODO: matching up different alleles of the same segment
## ask for columns to define uniqueness
## if not given, return all possible pairs of hydrogen bonds
##########################################
#' @name hydrogenBonds
#' Return a edge data.table connecting two input segments that are two strands of the same range
#' @param segs GRanges
#' @export
##########################################
hydrogenBonds = function(segs){
    ## MARCIN EDIT: fix to take care of situations where loose ends happen
    ## to exactly overlap a seg
    ## causing error here
    if (is.logical(segs$loose)){
        segss = paste(gUtils::gr.string(segs[, c("loose")]), segs$loose, sep="_")
        rsegss = gsub("ZZ", "-_", gsub("\\-_", "+_", gsub("\\+_", "ZZ", segss)))
        hb = match(segss, rsegss)
    } else {
        hb = match(segs, gr.flipstrand(segs))
    }
    hydrogenBs = data.table(from = seq_along(segs),
                            to = hb,
                            type = "hydrogen")
    return(hydrogenBs)
}



## what is this
sparse_subset = function (A, B, strict = FALSE, chunksize = 100, quiet = FALSE)
{
    nz = Matrix::colSums(A != 0, 1) > 0
    if (is.null(dim(A)) | is.null(dim(B))){
        return(NULL)
    }
    C = sparseMatrix(i = c(), j = c(), dims = c(nrow(A), nrow(B)))
    for (i in seq(1, nrow(A), chunksize)) {
        ixA = i:min(nrow(A), i + chunksize - 1)
        for (j in seq(1, nrow(B), chunksize)) {
            ixB = j:min(nrow(B), j + chunksize - 1)
            if (length(ixA) > 0 & length(ixB) > 0 & !quiet){
                cat(sprintf("\t interval A %s to %s (%d) \t interval B %d to %d (%d)\n",
                            ixA[1], ixA[length(ixA)], nrow(A), ixB[1],
                            ixB[length(ixB)], nrow(B)))
            }
            if (strict){
                C[ixA, ixB] = (sign((A[ixA, , drop = FALSE] !=
                                     0)) %*% sign(t(B[ixB, , drop = FALSE] != 0))) *
                    (sign((A[ixA, , drop = FALSE] == 0)) %*% sign(t(B[ixB,
                                                                    , drop = FALSE] != 0)) > 0)
            }
            else C[ixA, ixB] = (sign(A[ixA, nz, drop = FALSE] !=
                                     0) %*% sign(t(B[ixB, nz, drop = FALSE] == 0))) ==
                     0
        }
    }
    return(C)
}

#################################################
#' @name draw.paths.y
#' Determine the Y axis elevation of segments in a walk
draw.paths.y = function(grl, path.stack.x.gap=0, path.stack.y.gap=1){
    ## if grl is not named
    if (is.null(names(grl))){
        names(grl) = seq_along(grl)
    }

    if (any(is.na(names(grl)))){
        names(grl) = seq_along(grl)
    }

    grl.props = cbind(data.frame(group = names(grl), stringsAsFactors = F),
                      as.data.frame(values(grl)))

    gr = tryCatch(grl.unlist(grl),
                  error = function(e){
                      gr = unlist(grl);
                      if (length(gr)>0)
                      {
                          tmpc = textConnection(names(gr));
                          cat('budget .. \n')
                          gr$grl.ix = read.delim(tmpc, sep = '.', header = F)[,1];
                          gr$grl.iix = data.table::data.table(ix = gr$grl.ix)[
                                                     , iix := 1:length(ix), by = ix][, iix]
                          close(tmpc)
                      }
                      return(gr)
                  })

    gr$group = grl.props$group[gr$grl.ix]
    gr$group.ord = gr$grl.iix
    gr$first = gr$grl.iix == 1

    gr$last = iix = NULL ## NOTE fix, what is this??
    if (length(gr)>0){
        gr$last = data.table::data.table(
                                  iix = as.numeric(gr$grl.iix),
                                  ix = gr$grl.ix)[
                                , last := iix == max(iix), by = ix][, last]
    }
    grl.props$group = as.character(grl.props$group)
    S4Vectors::values(gr) =
        cbind(as.data.frame(values(gr)),
              grl.props[match(values(gr)$group, grl.props$group),
                        setdiff(colnames(grl.props),
                                c(colnames(values(gr)), 'group', 'labels')),
                        drop = FALSE])

    seqlevels(gr) = seqlevels(gr)[seqlevels(gr) %in% as.character(seqnames(gr))]
    windows = as(GenomicRanges::coverage(gr), 'GRanges'); ## Too deeply recursion
    windows = windows[values(windows)$score!=0]
    windows = reduce(windows, min.gapwidth = 1);

    win.gap = mean(width(windows))*0.2

    ## add 1 bp to end for visualization .. ranges avoids weird width < 0 error
    if (length(gr)>0)
    {
        IRanges::ranges(gr) =
            IRanges::IRanges(
                         start(gr),
                         pmax(end(gr),
                              ##                              pmin(end(gr)+1,
                              pmin(end(gr), ## FIXED BY MARCIN, above was causing needless stacking
                                   GenomeInfoDb::seqlengths(gr)[as.character(seqnames(gr))],
                                   na.rm = T),
                              na.rm = T)) ## jeremiah commented
    }

    suppressWarnings(end(windows) <- end(windows) + 1) ## shift one needed bc gr.flatmap has continuous convention, we have categorical (e.g. 2 bases is width 2, not 1)
    mapped = gr.flatmap(gr, windows, win.gap);

    grl.segs = mapped$grl.segs
    window.segs = mapped$window.seg

    dt = data.table(grl.segs)
    mx = dt[, max(c(as.numeric(pos1), as.numeric(pos2)))]
    int.mx = as.double(.Machine$integer.max)

    grl.segs$pos1 = round(as.double(grl.segs$pos1)/as.double(mx)*int.mx)
    grl.segs$pos2 = round(as.double(grl.segs$pos2)/as.double(mx)*int.mx)
    window.segs$start = round(as.double(window.segs$start)/as.double(mx)*int.mx)
    window.segs$end = round(as.double(window.segs$end)/as.double(mx)*int.mx)

    ix.l = lapply(split(1:nrow(grl.segs), grl.segs$group),
                  function(x) x[order(grl.segs$group.ord[x])])
    grl.segs$y.relbin = NA

    ## we want to layout paths so that we prevent collissions between different paths
    grl.segs$y.relbin[unlist(ix.l)] = unlist(lapply(ix.l, function(ix)
    {
        if (length(ix)>1)
        {
            iix = 1:(length(ix)-1)
            concordant = ((grl.segs$pos1[ix[iix+1]] >= grl.segs$pos2[ix[iix]]
                & grl.segs$strand[ix[iix+1]] != '-' & grl.segs$strand[ix[iix]] != '-') |
                (grl.segs$pos1[ix[iix+1]] <= grl.segs$pos2[ix[iix]]
                    & grl.segs$strand[ix[iix+1]] == '-' & grl.segs$strand[ix[iix]] == '-'))
            return(c(0, cumsum(!concordant)))
        }
        else{
            return(0)
        }
    }))

    contig.lim = data.frame(
        group = names(vaggregate(formula = y.relbin ~ group, data = grl.segs, FUN = max)),
        pos1  = vaggregate(formula = pos1 ~ group, data = grl.segs, FUN = min),
        pos2  = vaggregate(formula = pos2~ group, data = grl.segs, FUN = max),
        height = vaggregate(formula = y.relbin ~ group, data = grl.segs, FUN = max)
    );
    contig.lim$width = contig.lim$pos2 - contig.lim$pos1
    contig.lim$y.bin = 0;

    contig.lim = contig.lim[order(-contig.lim$width), ]

    if (nrow(contig.lim)>1){
        for (i in 2:nrow(contig.lim))
        {
            ir1 = IRanges::IRanges(contig.lim[1:(i-1), 'pos1'], contig.lim[1:(i-1), 'pos2'])
            ir2 = IRanges::IRanges(contig.lim[i, 'pos1'], contig.lim[i, 'pos2'])
            clash = which(ir1 %over% (ir2 + path.stack.x.gap))
            pick = clash[which.max(contig.lim$y.bin[clash] + contig.lim$height[clash])]
            contig.lim$y.bin[i] = c(contig.lim$y.bin[pick] + contig.lim$height[pick] + path.stack.y.gap, 0)[1]
        }
    }

    grl.segs$y.bin = contig.lim$y.bin[match(grl.segs$group, contig.lim$group)] + grl.segs$y.relbin + 1

    m.y.bin = max(grl.segs$y.bin)
    ylim = c(1, m.y.bin) + c(-0.5*m.y.bin, 0.5*m.y.bin)

    ## squeeze y coordinates into provided (or inferred) ylim
    tmp.ylim = ylim

    ## provide bottom and top padding of y.bin
    y.pad = 1/(m.y.bin+1)/2
    y.pad = pmin(1/(m.y.bin+1)/2, 0.125)
    tmp.ylim = tmp.ylim + c(1, -1)*y.pad*diff(tmp.ylim);

    ## make final y coordinates by squeezing y.bin into tmp.ylim
    grl.segs$y = affine.map(grl.segs$y.bin, tmp.ylim)

    ## MARCIN EDIT: grl.segs are not in order of paths
    ## but in coordinate order and so the order of the ys will be misintepreted
    ## down the line as being aligned to the order of segs in each path
    ## which will cause a mixup in the graphics

    grl.segs = grl.segs[order(grl.segs$group, grl.segs$group.ord), ]

    return(split(grl.segs$y, grl.segs$group)[names(grl)])
}

#' @name gr.flatmap
gr.flatmap = function(gr,
                      windows,
                      gap = 0,
                      strand.agnostic = TRUE,
                      squeeze = FALSE,
                      xlim = c(0, 1)){
    if (strand.agnostic){
        GenomicRanges::strand(windows) = "*"
    }

    ## now flatten "window" coordinates, so we first map gr to windows
    ## (replicating some gr if necessary)
                                        #    h = findOverlaps(gr, windows)

    h = gr.findoverlaps(gr, windows);

    window.segs = gr.flatten(windows, gap = gap)

    grl.segs = BiocGenerics::as.data.frame(gr);
    grl.segs = grl.segs[values(h)$query.id, ];
    grl.segs$query.id = values(h)$query.id;
    grl.segs$window = values(h)$subject.id
    grl.segs$start = start(h);
    grl.segs$end = end(h);
    grl.segs$pos1 = pmax(window.segs[values(h)$subject.id, ]$start,
                         window.segs[values(h)$subject.id, ]$start + grl.segs$start - start(windows)[values(h)$subject.id])
    grl.segs$pos2 = pmin(window.segs[values(h)$subject.id, ]$end,
                         window.segs[values(h)$subject.id, ]$start + grl.segs$end - start(windows)[values(h)$subject.id])
    grl.segs$chr = grl.segs$seqnames

    if (squeeze)
    {
        min.win = min(window.segs$start)
        max.win = max(window.segs$end)
        grl.segs$pos1 = affine.map(grl.segs$pos1, xlim = c(min.win, max.win), ylim = xlim)
        grl.segs$pos2 = affine.map(grl.segs$pos2, xlim = c(min.win, max.win), ylim = xlim)
        window.segs$start = affine.map(window.segs$start, xlim = c(min.win, max.win), ylim = xlim)
        window.segs$end = affine.map(window.segs$end, xlim = c(min.win, max.win), ylim = xlim)
    }

    return(list(grl.segs = grl.segs, window.segs = window.segs))
}

#' @name affine.map
affine.map = function(x,
                      ylim = c(0,1),
                      xlim = c(min(x), max(x)),
                      cap = F,
                      cap.min = cap,
                      cap.max = cap,
                      clip = T,
                      clip.min = clip,
                      clip.max = clip){
    if (xlim[2]==xlim[1]){
        y = rep(mean(ylim), length(x))
    }
    else{
        y = (ylim[2]-ylim[1]) / (xlim[2]-xlim[1])*(x-xlim[1]) + ylim[1]
    }

    if (cap.min){
        y[x<min(xlim)] = ylim[which.min(xlim)]
    }
    else if (clip.min){
        y[x<min(xlim)] = NA;
    }

    if (cap.max){
        y[x>max(xlim)] = ylim[which.max(xlim)]
    }
    else if (clip.max){
        y[x>max(xlim)] = NA;
    }

    return(y)
}

## setting CPLEX runtime parameters
cplex_customparams = function(out.file, numthreads = 0, nodefileind = NA, treememlim = NA)
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

## LANDFILL
## =================== functions whose fate to be determined
## ####################################################
## #' jabba.walk
## #'
## #' Computes walks around all aberrant edges in JABbA object
## #'
## #' Takes in JaBbA solution and computes local
## #' reconstructions around all aberrant edges (default).  Reconstructions (i.e. Huts) consists
## #' of collections of walks, each walk associated with a copy number, and a given
## #' region (collection of genomic windows).  The interval sum of walks in a given region, weighted
## #' by copy numbers will recapitulate the marginal copy profile (as estimated by JaBbA).
## #' The reconstruction is chosen to maximize parsimony.
## #'
## #' Optional flags allow making huts around specific junctions or specified loci (GRangesList)
## #'
## #' Walks are reconstructed locally within "clustersize" nodes of each aberrant edge, where
## #' clustersize is measured by the number of total edges.  Larger cluster sizes may fail to be
## #' computationally tractable, i.e. with a highly rearranged genome in an area of dense interconnectivity.
## #'
## #' @param sol JaBbA object
## #' @param outdir output directory
## #' @param junction.ix junction indices around which to build walks (default is all junctions)
## #' @param loci  loci around which to build walks (over-rides junction.ix), alternatively can be a list of  "all.paths" objects (i.e. each a list utput of initial all.paths = TRUE run  +/- field $prior for walk to re-eval a given all.paths combo
## #' @param clustersize size of the cluster to output around the locus or junction of interest
## #' @param trim logical flag whether trim in neighborhood of junction (only applicable if loci = NULL, default = TRUE)
## #' @param trim.w integer width to which to trim to
## #' @param prune flag whether to prune trivial walks for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
## #' @param prune.d1 local distance threshold for walk pruning
## #' @param prune.d2 referenc distance threshold for walk pruning
## #' @param mc.cores number of cores to use, default 1
## #' @param genes character vector of gene symbols with which to annotate walk (eg cancer genes)
## #' @param verbose logical flag
## #' @return list of walk set around each locus or junction that is inputted to analysis, each list item is a list with the following fields
## #' $win = input locus of interest, $grl = GRangesList of walks, $grs is a collapsed footprint of all walks in the walk list for this locu
## #' $td gTrack of of the output, additional outputs for debugging: $sol, $K, $Bc, $eix, $vix, $h
## ####################################################
## jabba.walk = function(sol, kag = NULL, digested = TRUE, outdir = 'temp.walk', junction.ix = NULL, loci = NULL, clustersize = 100,
##                       trim = FALSE, ## whether to trim around junction (only applicable when loci = NULL)
##                       trim.w = 1e6, ## how far to trim in neighborhood of junction (only applicable when loci = NULL
##                       prune = FALSE, ## whether to prune trivial walks i.e. those for whom a path can be drawn from first to last interval in a graph linking intervals with pairwise distance < d1 on the walk or distance < d2 on the reference
##                       prune.d1 = 1e5, ## local distance threshold for walk pruning
##                       prune.d2 = 1e5, ## reference distance threshold for walk pruning
##                       maxiterations = Inf,
##                       mc.cores = 1,
##                       genes = read.delim('~/DB/COSMIC/cancer_gene_census.tsv', strings = FALSE)$Symbol,
##                       verbose = TRUE,
##                       max.threads = 4,
##                       customparams = TRUE,
##                       mem = 6,
##                       all.paths = FALSE,
##                       nomip = FALSE,
##                       tilim = 100,
##                       nsolutions = 100,
##                       cb.interval = 1e4,
##                       cb.chunksize = 1e4,
##                       cb.maxchunks = 1e10)
## {
##     system(paste('mkdir -p', outdir))
##     ## awkward workaround to limit the number of processors Cplex will gobble up
##     ##

##     if (customparams)
##     {
##         out.file = paste(outdir, 'tmp.prm', sep = '/')
##         max.threads = Sys.getenv("LSB_DJOB_NUMPROC")
##         if (nchar(max.threads) == 0){
##             max.threads = Inf
##         }
##         else{
##             max.threads = as.numeric(max.threads)
##         }
##         max.threads = min(max.threads, mc.cores)
##         if (is.infinite(max.threads)){
##             max.threads = 0
##         }

##         param.file = paste(out.file, '.prm', sep = '')
##         cplex_customparams(param.file, max.threads, treememlim = mem * 1e3)

##         Sys.setenv(ILOG_CPLEX_PARAMETER_FILE = normalizePath(param.file))
##         print(Sys.getenv('ILOG_CPLEX_PARAMETER_FILE'))
##     }


##     if (is.null(sol)){
##         sol = kag
##     }

##     if (is.null(sol$segstats)){
##         sol$segstats = sol$tile
##         sol$segstats$cn = 2
##         sol$segstats$eslack.out = 0
##         sol$segstats$eslack.in = 0
##     }

##     if (is.null(kag)){
##         kag = sol
##     }


##     out = list()
##     tmp.adj = sol$adj

##     ## if input is already "digested", then don't need to bother with slacks
##     if (digested){
##         sol$segstats$eslack.in = 0
##         sol$segstats$eslack.out = 0
##         G = sol$G
##     }
##     else ## soon to be deprecated
##     {
##         ix = which(sol$segstats$eslack.in!=0 | sol$segstats$eslack.out!=0)
##         tmp.adj[ix, ix] = 0
##         pos.ix = which(as.logical(strand(sol$segstats)=='+'))
##         sol$segstats$tile.id = match(gr.stripstrand(sol$segstats), gr.stripstrand(sol$segstats[pos.ix]))
##         G = graph.adjacency(tmp.adj!=0)
##     }

##     h = jbaMIP.process(sol)

##     if (verbose){
##         cat(paste('Finished processing JaBbA, getting ready to construct walks\n'))
##     }

##     if (is.null(junction.ix) & is.null(loci)){
##         junction.ix = 1:nrow(kag$ab.edges)
##     }

##     if (!is.null(junction.ix)){
##         if (is.null(names(junction.ix))){
##             names(junction.ix) = 1:length(junction.ix)
##         }
##     }


##     if (is.null(loci)) ## junction.ix should be not null here
##     {
##         loci = do.call('GRangesList', mclapply(junction.ix, function(i){
##             if (verbose){
##                 cat(paste('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\nDefining subgraph around junction', i, '\n'))
##             }
##             vix = vix.i = setdiff(kag$ab.edges[i, 1:2, ], NA)
##             if (length(vix)==0){
##                 return(GRanges())
##             }
##             k = 0
##             last.clustersize = 0
##             while (length(vix)<clustersize & k < maxiterations & length(vix)>last.clustersize){
##                 k = k + 1
##                 last.clustersize = length(vix)
##                 vix = unique(unlist(neighborhood(G, vix.i, order = k)))
##             }
##             if (verbose){
##                 cat(paste('Outputting', length(vix), 'vertices around junction', i, '\n'))
##             }

##             return(kag$segstats[vix])
##         }, mc.cores = mc.cores))

##         names(loci) = names(junction.ix)
##         loci = loci[sapply(loci, length)>0]
##     }
##     else ## if loci are provided (i.e. not junction centric) then we will not trim or prune
##     {
##         trim = F
##         prune = F
##     }

##     if (verbose){
##         cat(paste('Finished defining subgraphs\n'))
##     }

##     starts = gr.start(sol$segstats, ignore.strand = F)
##     ends = gr.end(sol$segstats, ignore.strand = F)

##     names(sol$segstats) = 1:length(sol$segstats)

##     if (is.null(names(loci))){
##         lnames =  paste('locus', 1:length(loci), sep = '')
##     }
##     else{
##         lnames = names(loci)
##     }

##     all.junc.pair = c(paste(sol$ab.edges[, 1, 1], sol$ab.edges[, 2, 1], sep = ','), paste(sol$ab.edges[, 1, 2], sol$ab.edges[, 2, 2], sep = ','))
##     names(all.junc.pair) = c(1:nrow(sol$ab.edges), -c(1:nrow(sol$ab.edges)))

##     if (length(loci)>0){

##         out = mclapply(1:length(loci), function(i){

##             label = lnames[i]
##             outfile.rds = sprintf('%s/%s.rds', outdir, label)
##             outfile.pdf = sprintf('%s/%s.pdf', outdir, label)
##             outfile.txt = sprintf('%s/%s.txt', outdir, label)
##             outfile.allpaths.txt = sprintf('%s/%s.allpaths.txt', outdir, label)
##             if (inherits(loci[[i]], 'GRanges')){
##                 vix = which(gr.in(kag$segstats, loci[[i]]))
##                 cat('Number of vertices:', length(vix), '\n')
##                 eix = which((h$e.ij[,1] %in% vix | h$e.ij[,2] %in% vix) & h$e>0)
##                 Bc = as.matrix(h$B)[vix, eix]
##                 K = tryCatch(convex.basis(Bc, interval = cb.interval, chunksize = cb.chunksize, verbose = T, maxchunks = cb.maxchunks), error = function(e) as.character(e))
##                 if (is.character(K)){
##                     return(list(README = K))
##                 }
##                 prior = rep(1, ncol(K))
##             }
##             else ## assume we are re-heating a previous all.paths = TRUE output (and presumably adding a prior)
##             {
##                 K = loci[[i]]$K
##                 h = loci[[i]]$h
##                 eix = loci[[i]]$eix
##                 Bc = loci[[i]]$Bc
##                 vix = loci[[i]]$vix
##                 prior = rep(1, ncol(K))
##                 if (!is.null(loci[[i]]$prior)){
##                     prior[c(values(loci[[i]]$allpaths.og)$kix,values(loci[[i]]$allpaths.og)$kix2)]  = loci[[i]]$prior
##                 }
##                 loci[[i]] = loci[[i]]$win
##             }

##             is.cyc = Matrix::colSums(K[h$etype[eix] == 'slack', ])==0 & Matrix::colSums((h$B[, eix, drop = F] %*% K)!=0)==0
##             ## MEAT
##             karyo.sol = karyoMIP(K, h$e[eix], h$eclass[eix], nsolutions = nsolutions, tilim = tilim, cpenalty = 1/prior)
##             kag.sol = karyo.sol[[1]]
##             p = karyoMIP.to.path(kag.sol, K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores))
##             values(p$grl)$cn = p$cn
##             values(p$grl)$is.cyc = p$is.cyc
##             ##          td.rg$stack.gap = 5e6

##             if (!is.null(kag$junctions))
##             {
##                 values(kag$junctions)$lwd = sol$adj[kag$ab.edges[,1:2, 1]]
##                 values(kag$junctions)$lty = 1
##                 values(kag$junctions)$label = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, sol$adj[kag$ab.edges[,1:2, 1]], '')
##                 values(kag$junctions)$col = ifelse(sol$adj[kag$ab.edges[,1:2, 1]]>0, alpha('red', 0.3), alpha('white', 0))
##             }
##             win = streduce(sol$segstats[vix], 1e4)

##             y1 = max(sol$segstats$cn[gr.in(sol$segstats, win)], na.rm = T)*1.1
##             pdf(outfile.pdf, height = 30, width = 24)
##             grs = gr.simplify(grl.unlist(p$grl), 'grl.ix', split = T)
##             values(grs) = values(p$grl)
##             names(grs) = names(p$grl)

##             if (!is.null(sol$td)){
##                 td.seg = sol$td
##                 td.seg$y1 = y1
##                 td = td.seg
##                 ## td = c(td.seg, td.rg)
##             }
##             else{
##                 td.seg = gTrack(sol$segstats, y.field = 'cn', angle = 0, col ='black', height = 6, labels.suppress = T, y1 = y1)

##                 ## td = c(gTrack(grs, draw.paths = T, path.cex.arrow = 0, border = NA, angle = 0, ywid = 0.5, path.stack.x.gap = 1e6, height = 20, labels.suppress.gr = T),

##                 gt.walk = gTrack(grs, draw.paths = T, border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
##                 gt.walk$path.cex.arrow = 0
##                 gt.walk$path.stack.x.gap = 1e6
##                 td = c(gt.walk, td.seg)

##                 plot(td, windows = win, links = kag$junctions)

##                 dev.off()
##             }

##             df = data.frame(label = label, cn = p$cn, walk = sapply(grs, function(x) paste(gr.string(x, mb = F), collapse = ',')), widths = sapply(grs, function(x) paste(width(x), collapse = ',')), width = sapply(grs, function(x) sum(width(x))), numpieces = sapply(grs, length), type = 'walk')
##             df = rbind(data.frame(label = label, cn = NA, walk = paste(gr.string(win, mb = F), collapse = ','), widths = paste(width(win), collapse = ','), width = sum(width(win)), type = 'window', numpieces = length(win)), df)
##             write.tab(df, outfile.txt)
##             out = list(
##                 win = win, grl = p$grl, grls = grs, td = td, sol = karyo.sol,
##                 K = K, Bc = Bc, eix = eix, vix = vix, h = h,
##                 README = 'win=windows, grl = raw granges list corresponding to paths, grls = simplified granges list corresponding to paths, td = gTrack object plotting walks, sol = solution object from karyoMIP of local walks, K = incidence matrix input to karyomip, Bc = input to convex.basis, eix = eix input to karyomip, vix = vix input corresponding to rows of Bc, h = h input to karyomip')

##             if (all.paths){

##                 outfile.allpaths.pdf = sprintf('%s/%s.allpaths.pdf', outdir, label)

##                 if (verbose){
##                     cat('Generating all walks\n')
##                 }

##                 ## repurpose karyoMIP.to.path to generate all paths using "fake solution" i.e. all 1 weights,  to karyoMIP as input
##                 pallp = karyoMIP.to.path(list(kcn = kag.sol$kcn*0 + 1, kclass = kag.sol$kclass), K, h$e.ij[eix, ], sol$segstats, mc.cores = pmin(4, mc.cores), verbose = verbose)
##                 allp = pallp$grl

##                 allps = gr.simplify(grl.unlist(allp), 'grl.ix', split = T)
##                 allps[values(allp)$is.cycle] = do.call('GRangesList', lapply(which(values(allp)$is.cycle), function(x) c(allps[[x]], allps[[x]])))
##                 allps.og = allps; ## save for later
##                 values(allps.og)$kix = pallp$kix
##                 values(allps.og)$kix2 = pallp$kix2

##                 ## text encoding of junctions
##                 if (!is.null(junction.ix)){
##                     junc.pair = paste(sol$ab.edges[junction.ix[i], 1, ], sol$ab.edges[junction.ix[i], 2, ], sep = ',')
##                 }
##                 ## junction.ix should be not null here (i.e. they were provided as input or loci = NULL)
##                 if (trim | prune){
##                     allps.u = grl.unlist(allps)
##                     allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
##                     allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
##                     allps = split(allps.u, allps.u$grl.ix)
##                     allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
##                     allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
##                     allps.w = split(width(allps.u), allps.u$grl.ix)
##                     allps.endc = split(levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum), allps.u$grl.ix)

##                     ## only include windows around the junction of interest
##                     if (trim){
##                         ## allps.ix.pairs tells us what junction indices are present in a walk collection
##                         allps.ix.pairs = mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F)
##                         ## first, which windows contain the junction

##                         wix = which(sapply(allps.ix.pairs, length)>0)
##                         allps = allps[wix]

##                         if (length(allps)>0){
##                             allps.ixs = allps.ixs[wix] ## start interval id of kth interval in ith walk
##                             allps.ixe = allps.ixe[wix] ## end interval id of kth interval in ith walk
##                             allps.endc = allps.endc[wix] ## end walk coordinate of kth interval in ith walk
##                             allps.w = allps.w[wix]
##                             allps.ix.pairs = allps.ix.pairs[wix]

##                             ## start window for trimming
##                             values(allps)$allps.junc.first =
##                                             pmax(0, mapply(function(x, y) y[x[1]], allps.ix.pairs, allps.endc)) ## walk position of first junction
##                             values(allps)$allps.junc.last =
##                                             pmax(0, mapply(function(x, y) y[x[length(x)]], allps.ix.pairs, allps.endc)) ## walk position of last junction

##                             ## check for any quasi-palindromic walks that contain both orientations of a junction
##                             ## split each of these into two so we can maintain the width limit
##                             pal.wix = which(values(allps)$allps.win.firstix != values(allps)$allps.win.lastix)
##                             if (length(pal.wix)>0){
##                                 allps.dup = allps[pal.wix]
##                                 values(allps.dup)$allps.junc.first = values(allps)$allps.junc.last
##                                 allps = c(allps, allps.dup)
##                                 allps.endc = c(allps.endc, allps.endc[pal.wix])
##                                 allps.w = c(allps.w, allps.w[pal.wix])
##                             }

##                             values(allps)$allps.win.first =
##                                             pmax(0, values(allps)$allps.junc.first - trim.w) ## walk coordinate of new window start
##                             values(allps)$allps.win.last =
##                                             pmin(sapply(allps.endc, function(x) x[length(x)]), values(allps)$allps.junc.first + trim.w) ## walk coordinate of new window end
##                             values(allps)$allps.win.firstix = ## first walk interval to trim to
##                                             mapply(function(x, y) setdiff(c(which(x>y)[1], 1), NA)[1], allps.endc, values(allps)$allps.win.first)
##                             values(allps)$allps.win.lastix = ## last walk interval to trim to
##                                             mapply(function(x, y) setdiff(c(which(x>y)[1], length(x)), NA)[1], allps.endc, values(allps)$allps.win.last)
##                             values(allps)$allps.win.first.keep =
##                                             mapply(function(p,e,i) e[i] - p, values(allps)$allps.win.first, allps.endc, values(allps)$allps.win.firstix)
##                             values(allps)$allps.win.last.keep =
##                                             mapply(function(p,e,i,w) w[i] - (e[i] - p), values(allps)$allps.win.last, allps.endc, values(allps)$allps.win.lastix, allps.w)
##                             ## apply trimming
##                             ## we are trimming walks so that they are within trim.w bases of junction
##                             allps.u = grl.unlist(allps)
##                             iix = mapply(function(x,y) y %in% values(allps)$allps.win.firstix[x]:values(allps)$allps.win.lastix[x], allps.u$grl.ix, allps.u$grl.iix)
##                             allps.u = allps.u[iix]
##                             allps.u$keep.end = mapply(function(x, y)
##                                 ifelse(y == values(allps)$allps.win.firstix[x], values(allps)$allps.win.first.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)
##                             allps.u$keep.start = mapply(function(x, y)
##                                 ifelse(y == values(allps)$allps.win.lastix[x], values(allps)$allps.win.last.keep[x], NA), allps.u$grl.ix, allps.u$grl.iix)

##                             ## we keep the end of the first segment
##                             if (any(tmp.ix <- !is.na(allps.u$keep.start))) {
##                                 allps.u[tmp.ix] = gr.start(allps.u[tmp.ix], allps.u$keep.start[tmp.ix], ignore.strand = F)
##                             }

##                             ## we keep the beginning of the last segment
##                             if (any(tmp.ix <- !is.na(allps.u$keep.end))){
##                                 allps.u[tmp.ix] = gr.end(allps.u[tmp.ix], allps.u$keep.end[tmp.ix], ignore.strand = F)
##                             }

##                             ## if there are multiple walks with the same aberrant junction set, then pick the longest of these

##                             ## first need to find the aberrant walks in each set
##                             ij = paste(allps.u$ix.e[-length(allps.u)], allps.u$ix.s[-1], sep = ',') ## indices of all walk adjacent interval pairs
##                             names(ij) = 1:length(ij)
##                             ij = ij[diff(allps.u$grl.ix)==0] ## only pick intra-walk interval pairs
##                             ij.ix = names(all.junc.pair)[match(ij, all.junc.pair)]
##                             ## then compute the width of each walk

##                             allps = split(allps.u, allps.u$grl.ix)
##                             ij.ix.l = split(ij.ix, allps.u$grl.ix[as.numeric(names(ij))])[names(allps)]
##                             values(allps)$ab.junc = lapply(ij.ix.l, paste, collapse = ',')
##                             values(allps)$wid = vaggregate(width(allps.u), by = list(allps.u$grl.ix), FUN = sum)[names(allps)]
##                             ix.w = order(-values(allps)$wid)
##                             allps = allps[ix.w[which(!duplicated(values(allps)$ab.junc[ix.w]))]] ## only keep the longest non-duplicate walks
##                         }
##                     }

##                     ## now dedup and trim contigs to locus (mainly useful if loci was provided as argument)
##                     if (length(allps)>0){
##                         win = reduce(gr.stripstrand(loci[[i]]))
##                         allps.u = grl.unlist(allps)

##                         ## trim to locus
##                         ix = gr.match(allps.u, win)
##                         allps.u = allps.u[!is.na(ix)]
##                         ix = ix[!is.na(ix)]
##                         start(allps.u) = pmax(start(allps.u), start(win)[ix])
##                         end(allps.u) = pmin(end(allps.u), end(win)[ix])

##                         allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
##                         allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)

##                         ## remove dups
##                         allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of intervals
##                         allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of intervals

##                         allps.u = allps.u[allps.u$grl.ix %in% which(!duplicated(paste(sapply(allps.ixs, paste, collapse = ','), sapply(allps.ixe, paste, collapse = ','))))]
##                         allps = split(allps.u, allps.u$grl.ix)
##                     }

##                     ## this is to prune pseudo-aberrant walks that basically consist of short insertions of non-reference
##                     ## sequences in a big reference chunk
##                     if (prune & length(allps)>0){

##                         ## for each walk create graph of intervals by determining whether pair ij is BOTH near on the walk (<= d1)
##                         ## and near on the refernce (<= d2)
##                         allps.u = grl.unlist(allps)

##                         ## what are the ij pairs we want to test from this collapsed list
##                         ij = merge(cbind(i = 1:length(allps.u), ix = allps.u$grl.ix), cbind(j = 1:length(allps.u), ix = allps.u$grl.ix))[, c('i', 'j')]

##                         tmp = levapply(width(allps.u), by = list(allps.u$grl.ix), FUN = cumsum)
##                         allps.u.ir = IRanges(tmp - width(allps.u) + 1, tmp)

##                         ## distance on the walk
##                         D1 = sparseMatrix(ij[, 'i'],ij[, 'j'],
##                                           x = suppressWarnings(
##                                               distance(IRanges(start = end(allps.u.ir[ij[,'i']]), width = 1),
##                                                        IRanges(start(allps.u.ir[ij[,'j']]), width = 1))) + 1e-5, dims = rep(length(allps.u.ir), 2))

##                         ## distance on the reference
##                         D2 = sparseMatrix(ij[, 'i'],ij[, 'j'],
##                                           x = suppressWarnings(
##                                               distance(gr.end(allps.u[ij[,'i']], ignore.strand = F),
##                                                        gr.start(allps.u[ij[,'j']], ignore.strand = F))) + 1e-5, dims = rep(length(allps.u.ir), 2))

##                         D1 = pmin(as.matrix(D1), as.matrix(t(D1)))
##                         D2 = pmin(as.matrix(D2), as.matrix(t(D2)))

##                         tmp = D1>0 & D1<prune.d1 & D2>0 & D2<prune.d2
##                         tmp[which(is.na(tmp))] = FALSE
##                         G = graph.adjacency(tmp)
##                         cl = clusters(G, 'weak')$membership ## clusters based on this adjacency relationship
##                         cls = split(1:length(cl), cl)
##                         lens = sapply(allps, length)

##                         ## check if there any clusters that contain both the first and last member  of a walk
##                         cls.fl = cls[mapply(function(x) all(c(1,lens[allps.u$grl.ix[x[1]]]) %in% allps.u$grl.iix[x]), cls)]

##                         if (length(cls.fl)>0){
##                             toprune = allps.u$grl.ix[sapply(cls.fl, function(x) x[1])]
##                             if (length(toprune)>0){
##                                 cat('Pruning', length(toprune), 'walks\n')
##                             }
##                             allps = allps[-toprune]
##                         }
##                     }
##                 }

##                 if (length(allps)>0){
##                     win = streduce(unlist(allps), 0)
## ### win = streduce(unlist(allps), sum(width(unlist(allps)))*0)
##                 }

##                 values(allps) = NULL
##                 out$allpaths = allps
##                 out$allpaths.og = allps.og ## untouched all.paths if we want to reheat eg after computing 10X support
##                 gt.walk = gTrack(out$allpaths, draw.paths = T,border = NA, angle = 0, ywid = 0.5, height = 20, labels.suppress.gr = T)
##                 gt.walk$path.cex.arrow = 0
##                 gt.walk$path.stack.x.gap = 1e6
##                 out$td.allpaths = c(gt.walk, td.seg)
##                 pdf(outfile.allpaths.pdf, height = 30, width = 24)
##                 plot(out$td.allpaths, windows = win, links = kag$junctions)
##                 dev.off()
##                 out$README = paste(out$README, 'allpaths= all paths through windows (not just optimal ones), td.allpaths = gTrack object of plot of all paths')
##             }

##             ## if junction.ix was specified then label which positions in the walks represent the rearrangement junction
##             if (!is.null(junction.ix) & length(out$allpaths)>0){
##                 allps = out$allpaths
##                 allps.u = grl.unlist(allps)
##                 allps.u$ix.s = gr.match(gr.start(allps.u, ignore.strand = F), starts, ignore.strand = F)
##                 allps.u$ix.e = gr.match(gr.end(allps.u, ignore.strand = F), ends, ignore.strand = F)
##                 allps.ixs = split(allps.u$ix.s, allps.u$grl.ix) ## start indices of walk intervals in sol$segstats
##                 allps.ixe = split(allps.u$ix.e, allps.u$grl.ix) ## end indices of walks intervals in sol$segstats
##                 allps.ix.pairs = sapply(mapply(function(x,y) if (length(x)<=1) NULL else which(paste(x[-length(x)], y[-1], sep = ',') %in% junc.pair), allps.ixe, allps.ixs, SIMPLIFY = F), paste, collapse = ',')
##                 values(allps)$junction.id = names(junction.ix)[i]
##                 values(allps)$junction.ix = allps.ix.pairs
##                 out$allpaths = allps
##             }

##             if (length(out$allpaths)>0){
##                 values(out$allpaths)$string = grl.string(out$allpaths)
##                 values(out$allpaths)$wid = sapply(out$allpaths, function(x) sum(width(x)))
##                 values(out$allpaths)$wids = sapply(out$allpaths, function(x) paste(width(x), collapse = ','))
##                 write.tab(as.data.frame(values(out$allpaths)), outfile.allpaths.txt)
##             }

##             saveRDS(out, outfile.rds)
##             return(out)
##         }, mc.cores = mc.cores)
##     }

##     ## awkward workaround to limit the number of processors Cplex will gobble up
##     if (customparams){
##         system(paste('rm', param.file))
##         Sys.setenv(ILOG_CPLEX_PARAMETER_FILE='')
##         cat('Finished\n')
##     }

##     return(out)
## }
####################################
#' @name .e2class
#'
#' @title edge to contig class conversion
#'
#' @description given matrix K of k contigs over e edges, each belonging to cardinality 1 or cardinality 2 equivalence classes,
#' assigns id's to equivalent contigs
#'
####################################
.e2class = function(K, eclass){
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
        R = Matrix::diag(!(eclass %in% biclasses));  ## edges belonging to classes of cardinality 1 are on the diagonal

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
## #####################################################
## #' gtf2json
## #' Turning a GTF format gene annotation into JSON
## #' @importFrom rtracklayer import.gff
## ####################################################
## gtf2json = function(gtf=NULL,
##                     gtf.rds=NULL,
##                     gtf.gr.rds=NULL,
##                     filename="./gtf.json",
##                     genes=NULL,
##                     grep=NULL,
##                     grepe=NULL,
##                     chrom.sizes=NULL,
##                     include.chr=NULL,
##                     gene.collapse=TRUE,
##                     verbose = TRUE){

##     if (!is.null(gtf.gr.rds)){
##         message("Using GRanges from rds file.")
##         infile = gtf.gr.rds
##         gr = readRDS(gtf.gr.rds)
##         dt = gr2dt(gr)
##     } else if (!is.null(gtf.rds)){
##         message("Using GTF data.table from rds file.")
##         infile = gtf.rds
##         dt = as.data.table(readRDS(gtf.rds))
##     } else if (!is.null(gtf)){
##         message("Using raw GTF file.")
##         infile = gtf

##         gr = rtracklayer::import.gff(gtf)
##         dt = gr2dt(gr)
##     } else {
##         warning("No input gene annotation. Use the built-in GENCODE v19 in gUtils package")
##         require(skidb)
##         gr = read_gencode()
##         infile = "default"
##         dt = gr2dt(gr)
##     }

##     if (verbose){
##         message("Finished reading raw data, start processing.")
##     }

##     ## get seqlengths
##     if (is.null(chrom.sizes)){
##         message("No ref genome seqlengths given, use default.")
##         ## chrom.sizes = system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
##         ## system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
##         Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "hg19.regularChr.chrom.sizes", package="gUtils"))
##     }

##     sl = hg_seqlengths(include.junk=TRUE)

##     if (!is.null(include.chr)){
##         sl = sl[include.chr]
##     }
##     chrs = data.table(seqnames = names(sl), seqlengths=sl)

##     ## meta data field
##     require(RColorBrewer)
##     qw = function(x) paste0('"', x, '"') ## quote

##     meta.json =paste(paste0('\t',qw("metadata"),': [\n'),
##                      chrs[, paste("\t\t{",
##                                   qw("chromosome"),":", qw(seqnames),
##                                   ",", qw("startPoint"),":", 1,
##                                   ",", qw("endPoint"), ":", seqlengths,
##                                   ",", qw("color"),
##                                   ":", qw(substr(tolower(brewer.master( max(.I), 'BrBG' )), 1, 7)), " }",
##                                   collapse=",\n",
##                                   sep="")],
##                      '\n]')

##     if (verbose){
##         message("Metadata fields done.")
##     }

##     ## reduce columns: seqnames, start, end, strand, type, gene_id, gene_name, gene_type, transcript_id
##     ## reduce rows: gene_status, "KNOWN"; gene_type, not "pseudo", not "processed transcript"
##     dtr = dt[gene_status=="KNOWN" & !grepl("pseudo", gene_type) &
##              gene_type != "processed_transcript",
##              .(chromosome=seqnames, startPoint=start, endPoint=end, strand,
##                title = gene_name, gene_name, type, gene_id, gene_type,
##                transcript_id, transcript_name)]

##     if (!is.null(genes)){
##         dtr = dtr[title %in% genes]
##     }
##     else if (!is.null(grep) | !is.null(grepe)) {
##         if (!is.null(grep)){
##             dtr = dtr[grepl(grep, title)]
##         }
##         if (!is.null(grepe)){
##             dtr = dtr[!grepl(grepe, title)]
##         }
##     }

##     if (nrow(dtr)==0){
##         stop("Error: No more data to present.")
##     }

##     if (gene.collapse){
##         ## collapse by gene
##         dtr[, hasCds := is.element("CDS", type), by=gene_id]
##         dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","gene")],
##                     dtr[hasCds==FALSE][type %in% c("exon", "gene")])
##         ## dedup
##         dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, gene_id))]
##         dtr[, title := gene_name]
##         dtr = dtr[type != "transcript"]

##         ## group id
##         dtr[, gid := as.numeric(as.factor(gene_id))]
##         if (verbose){
##             message("Intervals collapsed to gene level.")
##         }
##     }
##     else {
##         ## collapse by transcript
##         dtr[, hasCds := is.element("CDS", type), by=transcript_id]
##         dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","transcript")],
##                     dtr[hasCds==FALSE][type %in% c("exon","transcript")])
##         ## dedup
##         dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, transcript_id))]
##         dtr[, title := transcript_name]
##         dtr = dtr[type != "gene"]

##         ## group id
##         dtr[, gid := as.numeric(as.factor(transcript_id))]
##         if (verbose){
##             message("Intervals collapsed to transcript level.")
##         }
##     }

##     dtr[, iid := 1:nrow(dtr)]

##     ## processing intervals
##     intervals.json = dtr[, paste0(
##         c(paste0(qw("intervals"),": ["),
##           paste(
##               "\t{",
##               qw("iid"), ":", iid,
##               ",", qw("chromosome"), ":", chromosome,
##               ",", qw("startPoint"), ":", startPoint,
##               ",", qw("endPoint"), ":", endPoint,
##               ",", qw("y"), ":", 0,
##               ",", qw("title"), ":", qw(title),
##               ",", qw("group_id"), ":", qw(gid),
##               ",", qw("type"), ":", qw(type),
##               ",", qw("strand"), ":", qw(strand),
##               "}",
##               sep = "",
##               collapse = ',\n'),
##           "]"),
##         collapse = '\n')
##         ]

##     ## assembling the JSON
##     out = paste(c("var dataInput = {", paste(
##                                            c(meta.json,
##                                              intervals.json),
##                                            collapse = ',\n'
##                                        ),"}"),
##                 sep = "")

##     writeLines(out, filename)
##     message(sprintf('Wrote JSON file of %s to %s', infile, filename))
##     return(filename)
## }
## read.js = function(file){
##     if (!file.exists(file)){
##         stop("File not found.")
##     }

##     require(data.table)
##     require(jsonlite)
##     js = read_json(file)

##     browser()
##     if (all(c("intervals", "connections") %in% names(js))){
##         intervals = rbindlist(js$intervals, fill=TRUE)
##         connections = rbindlist(js$connections, fill=TRUE)
##     } else {
##         stop("This is not a gGraph.json file.")
##     }

##     if ("walks" %in% names(js)){
##         walks.ls = js$walks
##         mc = rbindlist(lapply(walks.ls, function(x) x[1:4]))
##         mc = mc[, ":="(is.cycle = type=="cycle", str = strand)][, .(pid, cn, is.cycle, str)]

##         paths = lapply(walks.ls,
##                        function(x){

##                        })
##     }
## }
## ######################################################
## #' @name jabba.gwalk
## #' @title jabba.gwalk
## #' @description
## #'
## #' Computes greedy collection (i.e. assembly) of genome-wide walks (graphs and cycles) by finding shortest paths in JaBbA graph.
## #'
## #' @param jab JaBbA object
## #'
## #' @return GRangesList of walks with copy number as field $cn, cyclic walks denoted
## #' as field $is.cycle == TRUE, and $wid (width) and $len (segment length) of walks
## #' as additional metadata
## #'
## #' @import igraph
## #'
## #' @author Marcin Imielinski
## #' @author Xiaotong Yao
## #'
## #######################################################
## jabba.gwalk = function(jab, verbose = FALSE, return.grl = TRUE)
## {
##     segs = jab$segstats
##     hb = hydrogenBonds(segs)
##     hb.map = hb[, setNames(from, to)]
##     cn.adj = jab$adj
##     adj = as.matrix(cn.adj)
##     adj.new = adj*0
##     ## ALERT!!! see below
##     adj[which(adj!=0, arr.ind = TRUE)] = width(segs)[which(adj!=0, arr.ind = TRUE)[,2]] ## make all edges a large number by default
##     ## adj[which(adj!=0, arr.ind = TRUE)] = width(segs)[which(adj!=0, arr.ind = TRUE)[,1]] ## make all edges a large number by default
##     if (verbose){
##         ## ALERT!!! I'm gonna switch to source node width for default weight of edges
##         message('Setting edge weights to destination widths for reference edges and 1 for aberrant edges')
##         ## message('Setting default edge weights to SOURCE widths for edges and 1% less for aberrant edges')
##     }

##     ## ab.edges = rbind(jab$ab.edges[,1:2, 1], jab$ab.edges[,1:2, 2])
##     ## ab.edges = ab.edges[Matrix::rowSums(is.na(ab.edges))==0, ]
##     ab.edges = jab$edges[type=="aberrant", cbind(from, to)]
##   ## ALERT!!!
##   if (nrow(ab.edges)>0)
##   {
##     adj[ab.edges] = sign(cn.adj[ab.edges]) ## make ab.edges = 1
##   }
##     ## adj[ab.edges] = adj[ab.edges] * 0.99 ## make ab.edges 1 bp shorter than ref!
##     adj[is.na(adj)] = 0
##     cn.adj[which(is.na(cn.adj))] = 0

##     ## ALERT!!! major change
##     ## adjj = adj/as.matrix(cn.adj)
##     ## adjj[which(is.nan(adjj))] = 0
##     ## adjj[which(adjj<0)] = 0
##     ## G = graph.adjacency(adjj, weighted = 'weight')
##     ## esl = which(adj != 0, arr.ind=T)
##     ## eids = paste(esl[,1], esl[,2])
##     ## weights = adj[esl]
##     ## eclasses = ed[.(eids), eclass]
##     G = graph.adjacency(adj, weighted = 'weight')
##     ## G = make_graph(t(esl), )

##     ## DD = shortest.paths(G, mode="out")
##     ## IJ = which(!is.infinite(DD), arr.ind=T)

##     ## define ends not using degree (old method) but using either telomeres or loose ends
##     ## (otherwise lots of fake ends at homozygous deleted segments)
##     ss = gr2dt(segs)[ , vid:= 1:length(seqnames)]
##     ss[loose == TRUE, is.end := TRUE]
##     ss[loose == FALSE, is.end := 1:length(loose) %in% c(which.min(start), which.max(end)), by = list(seqnames, strand)]
##     ends = which(ss$is.end)
##     src = (Matrix::colSums(adj)[ends]==0) ## indicate which are sources

##     ## sanity check
##     unb = which(!ss$is.end & Matrix::rowSums(jab$adj, na.rm = TRUE) != Matrix::colSums(jab$adj, na.rm = TRUE))

##     if (length(unb)>0)
##     {
##         message(sprintf('JaBbA model not junction balanced at %s non-ends! Adding these to "ends"', length(unb)))
##         ends = c(ends, unb)         ## shameless HACK ... TOFIX
##     }

##     ## ends = which(degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
##     i = 0
##     ## adjust weight just before creating D
##     ## assign lighter weight to higher copy
##     ## D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

##     ## D records distance from ends to every node
##     D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]

##     ## sort shortest paths
##     ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]

##     ## ij only record end to end distance
##     ## ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[col %in% ends, ][, dist := D[cbind(row, col)]][, row := ends[row]][row != col, ][order(dist), ]

##     maxrow = length(ends)*max(cn.adj[ends, ends], na.rm = TRUE)
##     vpaths = rep(list(NA), maxrow)
##     epaths = rep(list(NA), maxrow)
##     cns = rep(NA, maxrow)
##     palindromic.path = rep(FALSE, maxrow)
##     palindromic.cycle = rep(FALSE, maxrow)

##     nb.all = which(Matrix::rowSums(cn.adj) != Matrix::colSums(cn.adj))
##     cn.adj0 = cn.adj
##     G0 = G
##     D0 = D

##     #' first peel off "simple" paths i.e. zero degree
##     #' ends with >0 copy number
##     psimp =  which(degree(G, mode = 'out')==0 & degree(G, mode = 'in')==0 & segs$cn>0)
##     i = 0
##     if (length(psimp)>0)
##     {
##         vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
##         epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA)) ## there is no "edge" associated with a zero total degree node
##         cns[1:length(psimp)] = segs$cn[psimp]
##         i = length(psimp)
##     }

##     ## now iterate from shortest to longest path
##     ## peel that path off and see if it is still there ..
##     ## and see if it is still there
##     ## peel off top path and add to stack, then update cn.adj

##     segs$tile.id = segs$tile.id + as.numeric(segs$loose)*0.5

##     tile.map =
##         gr2dt(segs)[, .(id = 1:length(tile.id),
##                         tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
##     rtile.map =
##         gr2dt(segs)[, .(id = 1:length(tile.id),
##                         tile.id = ifelse(strand == '+', 1, -1)*tile.id)]
##     setkey(tile.map, id)
##     setkey(rtile.map, tile.id)

##     ## unique pair of edge ids: rev comp of a foldback edge will be identical to itself!!!
##     ed = data.table(jab$edges)[cn>0, .(from, to , cn)]

##     if (nrow(ed)==0) ## make trivial gwalk of reference segments
##     {
##         ss = segs[segs$loose == FALSE, ]
##         paths = split(ss[, 'tile.id'], 1:length(ss))
##         values(paths)$ogid = 1:length(paths)
##         values(paths)$cn = ss$cn
##         values(paths)$label = paste('CN=', ss$cn, sep = '')
##         values(paths)$is.cycle = FALSE
##         values(paths)$numsegs = elementNROWS(paths)
##         values(paths)$num.ab = 0
##         values(paths)$wid = width(ss)
##     }
##     else {

##         ed[, ":="(fromss = tile.map[ .(from), tile.id],
##                   toss = tile.map[ .(to), tile.id]),
##            by = 1:nrow(ed)]
##         ed[, weight :=  adj[cbind(from, to)]]
##         print(ed)
##         ed[fromss*toss > 0, eclass := ifelse(fromss>0, paste(fromss, toss), paste(-toss, -fromss))]
##         ed[fromss*toss < 0, eclass := ifelse(abs(fromss)<=abs(toss),
##                                              paste(fromss, toss), paste(-toss, -fromss))]
##         ed[, eclass := as.numeric(as.factor(eclass))]
##         ed[, eid := paste(from, to)]
##         setkey(ed, "eid")
##         eclass.cn = ed[!duplicated(eclass), setNames(cn, eclass)]

##         cleanup_mode = FALSE


##         while (nrow(ij)>0)
##     {
##         if (verbose)
##             message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' end-pairs to resolve' )
##         i = i+1
##         ## swap this
##         ##        vpaths[[i]] = p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])

##         p = get.constrained.shortest.path(cn.adj, G, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

##         if (is.null(p)){
##             message('Came up empty!')
##             i = i -1
##             ij = ij[-1, , drop = FALSE]
##         }
##         else
##         {
##             ## Don't forget to update ed here
##             ed$cn = cn.adj[cbind(ed$from, ed$to)]

##             vpaths[[i]] = p
##             epaths[[i]] = cbind(p[-length(p)], p[-1])
##             eids = paste(epaths[[i]][,1], epaths[[i]][,2])
##             cns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palinrdomic edges by 1/2

##             rvpath = rtile.map[list(tile.map[list(vpaths[[i]]), -rev(tile.id)]), id]
##             repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
##             plen = length(rvpath)
##             hplen = floor(length(rvpath)/2)

##             ## (awkward) check for palindromicity for odd and even length palindromes
##             ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
##             if (ed[eids, any(table(eclass)>1)])
##                 palindromic.path[i] = TRUE
##             ## else
##             ## {
##             vpaths[[i+1]] = rvpath
##             epaths[[i+1]] = repath
##             cns[i+1] = cns[i]
##             palindromic.path[i+1] = TRUE
##             ## }
##             ##        palindromic = TRUE ## set to true while we "figure things out"


##             #' so now we want to subtract that cn units of that path from the graph
##             #' so we want to update the current adjacency matrix to remove that path
##             #' while keeping track of of the paths on the stack
##             cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]

##             ## if (!palindromic) ## update reverse complement unless palindromic
##             cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]
##             if (any(is.na(cn.adj))){
##                 browser()
##             }
##             if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
##             {
##                 message('backtracking ...') ## maybe we got stuck in a quasi-palindrome and need to backtrack
##                                         #            browser()
##                 cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]
##                 ## if (!palindromic) ## update reverse complement unless palindromic
##                 cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
##                 i = i-1
##                 ij = ij[-1, , drop = FALSE]
##             }
##             else ## continue, reduce
##             {
##                 adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
##                 ## if (!palindromic)
##                 adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

##                 ## ## make sure I didn't overuse any edge
##                 ## if (nrow(overdue <- which((as.matrix(jab$adj)-adj.new)<0, arr.ind=T))>0) {
##                 ##     print("Edge copy deficit!")
##                 ##     browser()
##                 ## }

##                 ## intermediate check
##                 ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
##                 ##     browser()

##                 to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,drop = FALSE]
##                 ## if (!palindromic) ## update reverse complement
##                 to.rm = rbind(to.rm, epaths[[i+1]][which(cn.adj[epaths[[i+1]]]==0), ,drop = FALSE])

##                 if (nrow(to.rm)>0)
##                 {
##                     adj[to.rm] = 0
##                     ## ALERT!!! major change
##                     ## adjj = adj/as.matrix(cn.adj)
##                     ## adjj[which(is.nan(adjj))] = 0
##                     ## adjj[which(adjj<0)] = 0
##                     G = graph.adjacency(adj, weighted = 'weight')
##                     ## G = graph.adjacency(adjj, weighted = 'weight')
##                     new.ends = setdiff(which(
##                     (degree(G, mode = 'out')==0 | degree(G, mode = 'in')==0)
##                     & degree(G)>0), ends)

##                     ## ## check if cn.adj out of balance
##                     ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
##                     ##     print("Junction OUT OF BALANCE!")
##                     ##     browser()
##                     ## }

##                     ## ## should be no new ends
##                     ## if (length(new.ends)>0){
##                     ##     print("Please, no new ends!")
##                     ##     browser()
##                     ## }

##                     ## remain = as.matrix(jab$adj) - adj.new
##                     ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
##                     ## if (any(!is.element(nb, nb.all)))
##                     ##     browser()

##                     D = shortest.paths(G, v = ends, mode = 'out', weight = E(G)$weight)[, ends]
##                     ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
##                 }
##                 else
##                     ij = ij[-1, , drop = FALSE]

##                 ## if (!palindromic) ## increase extra counter to account for reverse complement
##                 ## TOFIX: just update counter by 2 above, since we are just doing every path and its rc
##                 i = i+1
##             }
##         }


##         ## DEBUG DEBUG DEBUG
##         seg.ix = which(as.character(strand(segs))=='+'); seg.rix = which(as.character(strand(segs))=='-');


##         if (nrow(ij)==0 & cleanup_mode == FALSE)
##         {
##             message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR PATHS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
##             ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
##             cleanup_mode = TRUE
##         }
##     }
##         if (verbose)
##             message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )

##         ## ## record G, D, remaining edges at the end of path peeling
##         ## G1 = G
##         ## D1 = D
##         ## remain1 = remain

##         vpaths = vpaths[1:i]
##         epaths = epaths[1:i]
##         cns = cns[1:i]
##         palindromic.path = palindromic.path[1:i]

##         vcycles = rep(list(NA), maxrow)
##         ecycles = rep(list(NA), maxrow)
##         ccns = rep(NA, maxrow)

##         csimp = which(diag(cn.adj)!=0)
##         ipath = i
##         i = 0
##         if (length(csimp)>0)
##         {
##             vcycles[1:length(csimp)] = split(csimp, 1:length(csimp))
##             ecycles[1:length(csimp)] = lapply(csimp, function(x) cbind(x, x))
##             ccns[1:length(csimp)] = diag(cn.adj)[csimp]
##             cn.adj[cbind(csimp, csimp)] = 0
##             adj[cbind(csimp, csimp)] = 0
##             i = length(csimp)

##             for (j in 1:length(csimp))
##                 adj.new[ecycles[[j]]] = adj.new[ecycles[[j]]] + ccns[j]
##         }

##         ## sort shortest paths and find which connect a node to its ancestor (i.e. is a cycle)
##         .parents = function(adj)
##         {
##             tmp = apply(adj, 2, function(x) which(x!=0))
##             ix = which(sapply(tmp, length)>0)
##             if (length(ix)>0)
##             {
##                 parents = rbindlist(lapply(ix, function(x) data.table(x, tmp[[x]])))
##                 setnames(parents, c('node', 'parent'))
##                 setkey(parents, node)
##             } else {
##                 parents = data.table(node = 0, parent = NA)
##                 setkey(parents, node)
##             }
##         }

##         parents = .parents(adj)

##         #' then find paths that begin at a node and end at (one of its) immediate upstream neighbors
##         #' this will be a path for whom col index is = parent(row) for one of the rows
##         ## ALERT!!! major change
##         ## adjj = adj/as.matrix(cn.adj)
##         ## adjj[which(is.nan(adjj))] = 0
##         ## adjj[which(adjj<0)] = 0
##         G = graph.adjacency(adj, weighted = 'weight')
##         ## G = graph.adjacency(adjj, weighted = 'weight')
##         D = shortest.paths(G, mode = 'out', weight = E(G)$weight)

##         ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]


##         ## now iterate from shortest to longest path
##         ## peel that path off and see if it is still there ..
##         ## and see if it is still there

##         ## peel off top path and add to stack, then update cn.adj

##         cleanup_mode = FALSE
##         while (nrow(ij)>0)
##     {
##         if (verbose)
##             message('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
##         i = i+1
##                                         #        p = as.numeric(get.shortest.paths(G, ij[1, 1], ij[1, 2], mode = 'out', weight = E(G)$weight)$vpath[[1]])

##         p = get.constrained.shortest.path(cn.adj, G, allD = D, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

##         if (is.null(p)){
##             message('Came up empty!')
##             i = i -1
##             ij = ij[-1, , drop = FALSE]
##         } else
##         {

##             ed$cn = cn.adj[cbind(ed$from, ed$to)]
##             vcycles[[i]] = p
##             ecycles[[i]] = cbind(p, c(p[-1], p[1]))
##             eids = paste(ecycles[[i]][,1], ecycles[[i]][,2])
##             ccns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palindromic edges by 1/2

##             rvcycle = rtile.map[list(tile.map[list(vcycles[[i]]), -rev(tile.id)]), id]
##             recycle = cbind(rvcycle, c(rvcycle[-1], rvcycle[1]))
##             clen = length(rvcycle)
##             hclen = floor(length(rvcycle)/2)
##             ## (awkward) check for palindromicity for odd and even length palindromes

##             ## if (all((vcycles[[i]]==rvcycle)[c(1:hclen,(clen-hclen+1):clen)]))
##             if (ed[eids, any(table(eclass)>1)])
##                 palindromic.cycle[i] = TRUE
##             ## else
##             ## {
##             vcycles[[i+1]] = rvcycle
##             ecycles[[i+1]] = recycle
##             ccns[i+1] = ccns[i]
##             palindromic.cycle[i+1] = TRUE
##             ##     palindromic = FALSE
##             ## }
##             ##        palindromic = TRUE ## set to true while we "figure things out"

##             #' so now we want to subtract that cn units of that path from the graph
##             #' so we want to update the current adjacency matrix to remove that path
##             #' while keeping track of of the cycles on the stack
##             cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
##             ## if (!palindromic) ## update reverse complement unless palindromic
##             cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]
##             if (any(is.na(cn.adj))){
##                 browser()
##             }
##             if (!all(cn.adj[ecycles[[i]]]>=0))
##             {
##                 message('backtracking')
##                 ## browser()
##                 cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]+ccns[i]
##                 ## if (!palindromic) ## update reverse complement unless palindromic
##                 cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]+ccns[i+1]
##                 i = i-1
##                 ij = ij[-1, , drop = FALSE]
##             }
##             else
##             {
##                 adj.new[ecycles[[i]]] = adj.new[ecycles[[i]]] + ccns[i]

##                 ## ## if (!palindromic)
##                 ##     adj.new[ecycles[[i+1]]] = adj.new[ecycles[[i+1]]] + ccns[i]

##                 ## ## ## make sure I didn't overuse any edge
##                 ## ## if (length(overdue <- which((as.matrix(jab$adj)-adj.new)<0))) {
##                 ## ##     print("Edge copy deficit!")
##                 ## ##     browser()
##                 ## ## }

##                 ## ## ## intermediate cross check
##                 ## ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
##                 ## ##     browser()

##                 to.rm = ecycles[[i]][which(cn.adj[ecycles[[i]]]==0), ,drop = FALSE]

##                 ## if (!palindromic) ## update reverse complement
##                 to.rm = rbind(to.rm, ecycles[[i+1]][which(cn.adj[ecycles[[i+1]]]==0), ,drop = FALSE])

##                 if (nrow(to.rm)>0)
##                 {
##                     adj[to.rm] = 0
##                     parents = .parents(adj)
##                     ## G = graph.adjacency(adj, weighted = 'weight')

##                     ## ALERT!!! major change
##                     ## adjj = adj/as.matrix(cn.adj)
##                     ## adjj[which(is.nan(adjj))] = 0
##                     ## adjj[which(adjj<0)] = 0
##                     G = graph.adjacency(adj, weighted = 'weight')
##                     ## G = graph.adjacency(adjj, weighted = 'weight')

##                     ## if (any((Matrix::colSums(cn.adj)*Matrix::rowSums(cn.adj) != 0) & (Matrix::colSums(cn.adj) != Matrix::rowSums(cn.adj)))){
##                     ##     print("Junction OUT OF BALANCE!")
##                     ##     browser()
##                     ## }

##                     ## remain = as.matrix(jab$adj) - adj.new
##                     ## nb <- which(Matrix::colSums(remain) != Matrix::rowSums(remain))
##                     ## if (any(!is.element(nb, nb.all)))
##                     ##     browser()

##                     D = shortest.paths(G, mode = 'out', weight = E(G)$weight)
##                     ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]
##                 }
##                 else
##                     ij = ij[-1, ,drop = FALSE]

##                 ## if (!palindromic) ## increase extra counter to account for reverse complement
##                 i = i+1
##             }
##         }

##         if (nrow(ij)==0 & cleanup_mode == FALSE)
##         {
##             message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR CYCLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
##             ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row %in% parents$parent & row != col, ][order(dist), ][, is.cycle := parents[list(row), col %in% parent], by = row][is.cycle == TRUE, ]

##             cleanup_mode = TRUE
##         }
##     }

##         if (verbose)
##             message('Cycle peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )


##         if (i>0)
##         {
##             vcycles = vcycles[1:i]
##             ecycles = ecycles[1:i]
##             ccns = ccns[1:i]
##         }
##         else
##         {
##             vcycles = NULL
##             ecycles = NULL
##             ccns = NULL
##         }

##         vall = c(vpaths, vcycles)
##         eall = c(epaths, ecycles)
##         ecn = c(cns, ccns)

##         ## ## record G, D, remaining edges at the end of cycle peeling
##         ## G2 = G
##         ## D2 = D
##         ## remain2 = remain
##         remain = as.matrix(jab$adj) - adj.new
##         remain.ends = which(Matrix::colSums(remain)*Matrix::rowSums(remain)==0 & Matrix::colSums(remain)-Matrix::rowSums(remain)!=0)
##         if (length(remain.ends)>0){
##             if (verbose)
##                 message(length(remain.ends), "ends were not properly assigned a path. Do them.")
##         }

##         tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
##         ix = which(Matrix::rowSums(is.na(tmp[, 1:2]))==0)

##         if (length(ix)>0)
##             adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
##         else
##             adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
##         vix = munlist(vall)

##         segs$node.id = 1:length(segs)
##         pathsegs = segs[vix[,3]]
##         pathsegs$grl.ix = vix[,1]
##         browser()
##         ## abjuncs =  as.data.table(ab.edges)[,id := rep(1:(nrow(ab.edges)/2),2)*
##         ##                                         rep(c(1, -1), each = nrow(ab.edges)/2)][
##         ##     !is.na(from), ]
##         abjuncs = as.data.table(ab.edges)
##         abjuncs[, ":="(eid = paste(from, to),
##                        reid = paste(hb.map[as.character(to)],
##                                     hb.map[as.character(from)]))]
##         abjuncs[, ":="(ix = 1:.N,
##                   rix = match(reid, eid))]
##         abjuncs[, unique.ix := ifelse(rix>=ix, paste(ix, rix), paste(rix, ix))]
##         abjuncs[, eclass := as.numeric(as.factor(unique.ix))]
##         abjuncs[, iix := 1:.N, by=eclass]
##         abjuncs[, id := -(iix-1.5)/0.5*eclass]
##         abjuncs = abjuncs[, tag := structure(paste(from, to), names = id)]
##         setkey(abjuncs, tag)

##         ## annotate ab.id (if any) following each segment in each path
##         pathsegs$ab.id = gr2dt(pathsegs)[ , ab.id := c(abjuncs[paste(node.id[-length(node.id)], node.id[-1]), id], NA), by = grl.ix][, ab.id]

##         paths = split(pathsegs, vix[,1] )
##         values(paths)$ogid = 1:length(paths)
##         values(paths)$cn = ecn[as.numeric(names(paths))]
##         values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
##         values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
##         values(paths)$numsegs = elementNROWS(paths)
##         values(paths)$num.ab = sapply(paths, function(x) sum(!is.na(x$ab.id)))
##         values(paths)$wid = sapply(lapply(paths, width), sum)

##         check = which((adj.new - jab$adj) !=0, arr.ind = TRUE)

##         if (length(check)>0)
##             stop('Alleles do not add up to marginal copy number profile!')
##         else if (verbose)
##             message('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')
##     }

##     ## match up paths and their reverse complements
##     psig = lapply(paths, function(x) ifelse(as.logical(strand(x)=='+'), 1, -1)*x$tile.id)
##     psig.flip = sapply(psig, function(x) -rev(x))

##     unmix = data.table(
##         ix = 1:length(paths),
##         mix = match(sapply(psig, paste, collapse = ','), sapply(psig.flip, paste, collapse = ',')))[, pos := 1:length(mix)<mix][order(!pos), ]
##     setkey(unmix, ix)
##     unmix[is.na(mix), pos := TRUE] ## if we have paths with no reverse complement i.e. NA mix, then call "+" for now

##     remix = rbind(
##         unmix[pos == TRUE, ][, id := 1:length(ix)],
##         unmix[list(unmix[pos == TRUE, mix]), ][, id := 1:length(ix)][!is.na(ix), ]
##     )

##     paths = paths[remix$ix]
##     names(paths) = paste(remix$id, ifelse(remix$pos, '+', '-'), sep = '')
##     values(paths)$id = remix$id
##     values(paths)$str = ifelse(remix$pos, '+', '-')

##     if (length(setdiff(values(paths)$ogid, 1:length(paths))))
##         message('Warning!!! Some paths missing!')

##     ## for gGnome compatibiliity
##     if (!return.grl)
##     {
##         tmp.dt = as.data.table(copy(paths))[, pid := group_name][, nix := 1:.N, by =pid]
##         setkeyv(tmp.dt, c('pid', 'nix'))

##         ## mark nodes that precede a reference junction
##         tmp.dt[, d.to.next := c((start-data.table::shift(end))[-1], NA), by = pid]
##         tmp.dt[, d.to.next.neg := c((data.table::shift(start)-end)[-1], NA), by = pid]
##         tmp.dt[, same.strand := c((strand==data.table::shift(strand))[-1], NA), by = pid]
##         tmp.dt[, same.chrom := c((as.character(seqnames)==data.table::shift(as.character(seqnames)))[-1], NA), by = pid]
##         tmp.dt[, last.node := 1:.N == .N, by = pid]
##         tmp.dt[, before.ref :=
##                      (((d.to.next<=1 & d.to.next>=0 & strand == '+') |
##                        (d.to.next.neg<=1 & d.to.next.neg>=0 & strand == '-')
##                      ) & same.strand & same.chrom)]
##         tmp.dt[is.na(before.ref), before.ref := FALSE]

##         ## label reference runs of nodes then collapse
##         .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), as.integer(NA))
##         tmp.dt[, ref.run := .labrun(before.ref), by = pid]
##         tmp.dt[, ref.run.last := data.table::shift(ref.run), by = pid]
##         tmp.dt[is.na(ref.run) & !is.na(ref.run.last), ref.run := ref.run.last]
##         tmp.dt[!is.na(ref.run), ref.run.id := paste(pid, ref.run)]

## ### TODO: store ab.ids in walks
##                                         #tmp.dt[loose == TRUE, ref.run.id := NA] ## make sure loose ends stay ungrouped
##         if (any(!is.na(tmp.dt$ref.run.id)))
##         {
##             collapsed.dt = tmp.dt[!is.na(ref.run.id), .(
##                                                           nix = nix[1],
##                                                           pid = pid[1],
##                                                           seqnames = seqnames[1],
##                                                           start = min(start),
##                                                           end = max(end),
##                                                           loose = FALSE,
##                                                           strand = strand[1]
##                                                       ), by = ref.run.id]
##             tmp.dt = rbind(
##                 tmp.dt[is.na(ref.run.id),
##                        .(pid, nix, seqnames, start, end, strand, loose)],
##                 collapsed.dt[, .(pid, nix, seqnames, start, end, strand, loose)])

##         }

##         ## concatenate back with nodes that precede a non reference junctiono
##         setkeyv(tmp.dt, c('pid', 'nix'))

##         tmp.gr = dt2gr(tmp.dt)
##         tmp.segs = unique(tmp.gr)
##         tmp.gr$seg.id = match(tmp.gr, tmp.segs)
##         tmp.paths = split(tmp.gr$seg.id, tmp.gr$pid)
##         tmp.vals = as.data.frame(values(paths[names(tmp.paths)]))

##         names(tmp.paths) = ifelse(grepl('\\-', names(tmp.paths)), -1, 1)*as.numeric(gsub('\\D', '', names(tmp.paths)))

##         ##      gw = gGnome::gWalks$new(segs=tmp.segs,
##         gw = gWalks$new(segs=tmp.segs,
##                         paths=tmp.paths,
##                         metacols=tmp.vals)
##         return(gw)
##     }
##     return(paths)
## }
