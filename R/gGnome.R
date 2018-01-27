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
#' @import parallel
#' @import IRanges
#' @import GenomicRanges
#' @import data.table
#' @import igraph
#' @import S4Vectors
#' @import Matrix
#' @import gUtils
#' @import gTrack
#' @import Rcplex
#'
NULL

#' TODO:
#'
#' simplify
#' decouple
#' add
#' subtract
#' gg$junctions, gg$ab.edges
#' junctions: ra.merge, ra.dedup, ra.dist, ra.equal
#' gwalks: as.gGraph, write.json
#' jgraph
#' find.fusion
#' proximity
#'
#' Naming:
#' -- use S3 to overload exported function names with pure lowercase
#' -- the exposed fields should have lower case names too
#' -- Okay to keep the camel cases of object or internal method names for now,
#' to replace gradually later
#' -- arithmetics of graphs
#' -- ## TODO: develop the most efficient way to r/w GFA1 format
#'
#' Structure:
#' -- use S4 to extend GRL class to junctions, don't use R6
#'
#' Documentation:
#' -- define exported functions and for each come up with a short use case
#'
#' Final destination: I want something like
#' -- jab = read.jab("jabba.simple.rds"); and jab is a bGraph unless told otherwise
#' -- gw = gwalk(jab); and gw is a gWalk unless told otherwise
#' -- plot(jab) plots the default gTrack and returns in situ
#' -- write.json(jab); write.json(gw); saves the JSON format for viz
#' -- hood(jab, win, d=1e6); returns the +/-1Mb neighborhood of the
#' -- subgraph(jab, expr); returns the subgraph where the nodes evaluate to TRUE in expr
#' --

##############################
## junctions
##############################
#' junctions
#' S4 wrapper around GRangesList to store junction info
#'
#' @import methods
#' @import gUtils
#' @export
junctions = setClass("junctions",
                     contains="GRangesList")
## validity test when intializing
setValidity("junctions",
            function(object){
                if (length(object)==0){
                    message("Empty junction set.")
                    return(TRUE)
                } else if (!all(elementNROWS(object)==2))
                    "Each element must be length 2."
                else if (!all(strand(unlist(object)) %in% c("+", "-")))
                    "All strand info must be present."
                else
                    TRUE
            })
## explicit coercion and that's it!
setAs("GRangesList", "junctions", function(from){new("junctions", from)})
## now extend S4 methods special for "junctions"
## size?
## set operations!!!
## union, setdiff, xor, union




#' gGraph
#'
#' the central class for rearrangement graphs
#'
#' @import R6
#' @import data.table
#' @import Matrix
#' @import igraph
#' @import gUtils
#' @import gTrack
#'
#' @export
gGraph = R6Class("gGraph",
                 public = list(
                     ## public fields
                     ## name = NULL,
                     refG = "GENOME", ## seqinfo of ref genome

                     ## constructor
                     initialize = function(tile=NULL, junctions=NULL, cn = FALSE,
                                           jabba=NULL, weaver=NULL, prego=NULL,
                                           segs=NULL, es=NULL, ploidy=NULL, purity=NULL,
                                           regular=TRUE, rescue.balance=FALSE){
                         ## control how to construct
                         if (!is.null(segs) & !is.null(es)){
                             private$gGraphFromScratch(segs, es,
                                                       junctions, ploidy, purity)
                         } 
                         else if (!is.null(tile) | !is.null(junctions)) {
                             message("Initializing with 'tile' and 'junctions'")
                             self$karyograph(tile, junctions, cn = cn)
                         } 
                         else if (!is.null(jabba)) {
                             ## message("only use 'jabba' or 'weaver' field, not both")
                             message("Reading JaBbA output")
                             self$jabba2gGraph(jabba)
                             if (rescue.balance){
                                 self$rescueBalance()
                             }
                         } 
                         else if (!is.null(weaver)) {
                             ## message("only use 'jabba' or 'weaver' field, not both")
                             message("Reading Weaver output")
                             self$weaver2gGraph(weaver)
                             if (rescue.balance){
                                 self$rescueBalance()
                             }
                         } 
                         else if (!is.null(prego)) {
                             message("Reading Prego output")
                             self$prego2gGraph(prego)
                             if (rescue.balance){
                                 self$rescueBalance()
                             }
                         } 
                         else {
                             ## generate null graph
                             self$nullGGraph(regular)
                         }
                     },

                     ## initialize from global ref genome seqinfo
                     ## This is actually a diploid graph
                     ## null graph should be all zero
                     nullGGraph = function(regular=TRUE){
                         private$segs = GRanges()
                         private$g = make_empty_graph()
                         private$junction = new("junctions", GRangesList())
                         private$es = data.table(from=integer(0),
                                                 to=integer(0),
                                                 type=character(0))
                         return(self)
                     },

                     ## start as diploid chromosomes
                     dipGraph = function(genome = NULL, chr=FALSE, regular=TRUE){
                         ## default seqlengths
                         sl = hg_seqlengths(genome = genome, chr = chr)
                         tmp = si2gr(sl)
                         names(tmp) = NULL

                         if (regular){
                             tmp = tmp %Q% (seqnames %in% regularChr)
                         }

                         private$segs = c(tmp, gr.strandflip(tmp)) ## null segs are ref
                         private$segs$cn = private$.ploidy = 2 ## diploid
                         private$segs$loose = FALSE ## all non-loose end

                         private$es = data.table(from = integer(0),
                                                 to = integer(0),
                                                 type = character(0),
                                                 cn = numeric(0))

                         if (!regular){
                             ## close the circle on mitochondrial DNA
                             sinfo = as.data.frame(seqinfo(get(self$refG)))
                             sinfo = data.table(seqnames=rownames(sinfo), sinfo)
                             circChr = sinfo[isCircular==T, seqnames]

                             ## QUESTION: why Rle doesn't allow match function???
                             circIx = which(as.vector(seqnames(private$segs)) %in% circChr)
                             if ( length(circIx)>0 ){
                                 ## constructing edges: 5 required columns
                                 ## from, to, cn, type, weight (len o' source node)
                                 private$es = rbindlist(
                                     list(private$es, list(
                                                          from=circIx,
                                                          to=circIx,
                                                          type=rep("reference", length(circIx)),
                                                          cn=private$segs$cn[circIx],)
                                          )
                                 )

                                 ## Q: constructing an igraph easier than modifying?
                                 private$g = add_edges(private$g,
                                                       t(as.matrix(private$es[,.(from,to)])),
                                                       ## from-to-from-to-...
                                                       attr = as.list(private$es))
                             }
                         }

                         return(self)
                     },

                     ## initialize from segmenatation AND/OR rearrangement junctions
                     addJuncs = function(junc){
                         ## DONE: populate abEdges while adding new junction!!!!
                         ## NOTE: the bps in junc must be width 2
                         ## TODO: what if junctions come with a CN?
                         ## ALERT: convention of junction orientation!!!
                         "Given a GRL of junctions add them plainly to this gGraph."
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
                             junc = junctions(junc)
                         }

                         if (is.null(private$junction)){
                             tmp = self$e2j()
                             rm(tmp); gc(verbose=FALSE)
                         }

                         junctions = junc
                         jadd = data.table(jix = seq_along(junctions)) ## determine what to add
                         ## save the junctions in the object
                         ## DONE: what if I am adding some existing junctions that are just not
                         ## incorporated???
                         j.ov = na.omit(ra.overlaps(junctions, private$junction))
                         j.exist = data.table(ra1.ix=numeric(0), ra2.ix=numeric(0))
                         if (nrow(na.omit(j.ov))>0 &
                             "ra1.ix" %in% colnames(na.omit(j.ov)) &
                             "ra2.ix" %in% colnames(na.omit(j.ov))){
                             if (length(new.jix <- setdiff(seq_along(junc), j.exist[, ra1.ix]))>0){
                                 private$junction = c(private$junction, junc[new.jix])
                             }
                             jadd[j.exist[, ra1.ix], exist := j.exist[, ra2.ix]]
                         } else {
                             jadd[, exist := as.numeric(NA)]
                         }


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
                         jadd[, j.in := jIn]

                         if ("cn" %in% colnames(values(junctions))){
                             jadd[, cn := values(junctions)$cn]
                         } else {
                             jadd[, cn := 1]
                         }

                         if (jadd[, !any(cn>0 & j.in==TRUE)]){
                             return(self)
                         }

                         ## for existing junctions modify the copy number
                         tomod = jadd[, which(j.in==TRUE & cn>0 & !is.na(exist))]
                         values(private$junction)$cn[jadd[tomod, exist]] =
                                                    values(private$junction)$cn[jadd[tomod, exist]] + values(junctions)$cn[tomod]

                         ## resize to width 1, left
                         jUl = grl.unlist(junctions)
                         if (!all(width(jUl))==1){
                             jUl = gr.start(jUl)
                         }
                         names(jUl) = NULL

                         ## start processing
                         ## DONE: write as JaBbA::karyograph() with modifications
                         ## e.g. (30, 2) --> pivot (2, 30)
                         jadd = jadd[j.in==TRUE & cn>0, jix]
                         bp.p = split(jUl %Q% (grl.ix %in% jadd), rep(1:2, length(jadd)))
                         bp.p = gr.fix(bp.p, get(self$refG))
                         juncTile = c(bp.p[[1]], bp.p[[2]])
                         ## BP 1 and 2, retaining strand-orientation info
                         ## bp1 = gr.start(bp.p[[1]], ignore.strand = T)
                         ## bp2 = gr.start(bp.p[[2]], ignore.strand = T)
                         ## tmpBps = c(bp1, bp2)
                         self$addSegs(juncTile) ## DONE: addSegs is giving dup edges!!!

                         ## sanity check: at this point, all bps should map to only right bound
                         if (!all(unlist(bp.p) %^% gr.end(private$segs))){
                             ## browser()
                             stop("Something went wrong when breaking up the segs!")
                         }

                         hb =self$hydrogenBonds()
                         hb = hb[, c(setNames(from, to), setNames(to, from))]
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

                         abEs = rbind(
                             data.table(from = from1, to = to2, edge.ix = seq_along(anc1)),
                             data.table(from = from2, to = to1, edge.ix = seq_along(anc2)))
                         ## weight is assigned to the width of source node
                         abEs[, cn := c(bp.p[[1]]$cn, bp.p[[2]]$cn)]
                         abEs = abEs[,.(from, to, cn, type="aberrant")]

                         ## make edges
                         private$es = rbind(private$es, abEs)

                         ## connecting the junctions
                         private$g = add_edges(graph = private$g,
                                               edges = as.vector(t(as.matrix(abEs[, .(from, to)]))),
                                               attr = as.list(abEs[,.(cn, type)]))

                         tmp = self$e2j()
                         return(self)
                     },

                     addSegs = function(tiles){
                         ## Given a GRanges obj of a segmentation (complete or not),
                         ## break the gGraph at their ends.
                         ## extract breakpoints
                         ## bps = reduce(c(gr.start(tile), gr.end(tile)))
                         bps = tiles
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
                         if (is.null(private$es)){
                             ## NOTE: don't understand why es is NULL sometimes
                             private$es = data.table(from = integer(0),
                                                     to = integer(0),
                                                     type = character(0),
                                                     cn = numeric(0))
                         }

                         if (!"cn" %in% colnames(private$es)){
                             private$es[, cn := as.numeric(NA)] ## TODO: will this work?
                         }

                         private$es[, .(from = tmpDt[, which(qid %in% from & isTail==T)],
                                        to = tmpDt[, which(qid %in% to & isHead==T)],
                                        cn, type)] -> newEs

                         ## introduce ref edges between new breakpoints
                         refEs = tmpDt[, .(from=.I[isTail==F], to=.I[isHead==F]), by=qid]
                         refEs[, ":="(cn=tmpDt[c(from, to), min(cn)], type="reference")]
                         refEs = refEs[!is.na(from) & !is.na(to), -c("qid"), with=F]

                         newEs = rbindlist(list(newEs, refEs)) ## combine the two parts
                         newEs[!duplicated(newEs)]

                         ## update: es, g
                         private$es = newEs
                         private$g = graph_from_data_frame(
                             newEs, directed = TRUE,
                             vertices=data.frame(name = as.integer(rownames(tmpDt)),
                                                 cn = tmpDt[, cn]))
                         ## reset tmpSegs
                         private$tmpSegs = NULL

                         return(self)
                     },

                     ## karograph: initialize `nullGGraph()`,
                     ## add junctions to it, then add tiles to it
                     karyograph = function(tile=NULL, juncs=NULL, cn=FALSE, regular=FALSE){
                         ## TODO: make this compatible with JaBbA!!
                         self$dipGraph(regular = regular)

                         ## no tile, no cn
                         if (is.null(tile)){
                             cn = FALSE
                         } else if (!"cn" %in% colnames(values(tile))){
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
                             self$addSegs(c(tile[,c()], gr.stripstrand(unlist(juncs[jadd])[,c()])))
                             self$addJuncs(juncs)
                             if (cn == TRUE) {
                                 private$segs = private$segs %$% tile
                                 ## TODO: if anything drops below edge CN sum,
                                 ## tune down the edge CN too
                                 node.cn = data.table(id=seq_along(private$segs),
                                                      cn=private$segs$cn)
                             }
                         } else if (!is.null(tile) & length(tile)>0){
                             self$addSegs(tile)
                         } else if (!is.null(juncs) & length(juncs)>0){
                             ## if empty, ignore these GRanges lists
                             self$addJuncs(juncs)
                         }
                         return(self)
                     },

                     simplify = function(){
                         ## if two or more segment are only connected by ref edges
                         ## and they have the same copy number
                         ## merge them into one node
                         "Merge all pairs of noded only connected by reference edge."
                         verbose = getOption("gGnome.verbose")
                         browser()
                         ## MOMENT
                         ## get the part of the graph where the nodes are
                         ## those at least one side is connecting a single reference edge
                         node.dt = data.table(nix = seq_along(private$segs))
                         private$es[]

                     },

                     decouple = function(){
                         "When there's overlapping nodes, break them down and reconnect."
                         if (isDisjoint(private$segs %Q% (strand=="+" & loose==FALSE)))
                             return(self)
                         ## MOMENT

                     },

                     add = function(gg){
                         "Simply put two gGraphs together."
                         verbose = getOption("gGnome.verbose")

                         if (!is(gg, "gGraph"))
                             stop("Can only deal with addition of two gGraph objects.")

                         node = private$segs
                         if (!"cn" %in% colnames(values(node))){
                             node$cn = 1
                             if (verbose) warning("Arbitrarilly set node copy numbers to 1.")
                         }

                         node2 = gg$segstats
                         if (!"cn" %in% colnames(values(node2))){
                             node2$cn = 1
                             if (verbose) warning("Arbitrarilly set node copy numbers to 1.")
                         }

                         new.segs = c(node, node2) ## ALERT: dedup later
                         ## LATER
                     },

                     ## initialize from JaBbA output
                     jabba2gGraph = function(jabba, regular.only=F){
                         ## ptm = proc.time()
                         if (is.list(jabba)) {
                             if (all(is.element(c("segstats", "adj", "ab.edges", "edges", "G",
                                                  "td", "purity", "ploidy", "junctions"),
                                                names(jabba))))
                                 jabba = jabba
                         } else if (is.character(jabba) & grepl(".rds$", jabba)){
                             if (file.exists(jabba))
                                 jabba = readRDS(jabba)
                         } else {
                             stop("Input must be either JaBbA list output or the RDS file name that contains it!")
                         }

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

                         ## MARCIN DEBUG: THIS MISSES WHOLE CHROMOSOMES WHICH
                         ## HAVE BOTH ZERO IN AND ZERO OUT DEGREE!!
                         ## whichTerminal = private$es[, setxor(from, to)]
                         ## private$segs$terminal = seq_along(private$segs) %in% whichTerminal
                         private$segs$terminal  = !(1:length(private$segs) %in% private$es$from) |
                             !(1:length(private$segs) %in% private$es$to)

                         ## private$abEdges = jabba$ab.edges
                         private$.ploidy = jabba$ploidy
                         private$.purity = jabba$purity

                         if (regular.only==T){
                             self$trim(gr.stripstrand(si2gr(hg_seqlengths()[regularChr])))
                         }
                         ## cat("subgraph done")
                         ## print(proc.time() - ptm)
                         ## cat("\n")

                         return(self)
                     },

                     ## initialize from Weaver result
                     weaver2gGraph = function(weaver){
                         ## DONE: get Weaver done!!!! GEt it done@!!!
                         ## input weaver: directory that contains three and only three files
                         ## named: REGION_CN_PHASE, SV_CN_PHASE, and SNP_CN_PHASE
                         if (!dir.exists(weaver)) stop("Invalid input weaver directory!")

                         if (!all(
                                  is.element(c("SV_CN_PHASE", "REGION_CN_PHASE"),
                                             dir(weaver)))
                             ){
                             stop('Need "SV_CN_PHASE" and "REGION_CN_PHASE".')
                         }

                         sl = fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]

                         require(data.table)
                         region = data.table(read.delim(
                             paste(weaver, "REGION_CN_PHASE", sep="/"),
                             header = FALSE, sep = "\t"))

                         sv.fn = paste(weaver, "SV_CN_PHASE", sep="/")
                         if (file.size(sv.fn)>0){
                             sv = data.table(read.delim(sv.fn, header = FALSE, sep = "\t"))
                             names(sv) = c("chr1", "pos1", "side1", "allele1",
                                           "chr2", "pos2", "side2", "allele2",
                                           "cn", "unknown1", "unknown2", "timing", "class")[1:ncol(sv)]
                         } else {
                             sv = NULL
                         }

                         ## define the columns
                         names(region) = c("seqnames", "start", "end", "acn", "bcn")
                         region[, cn := acn + bcn]
                         ## names(snp) = c("seqnames", "pos", "ref", "alt", "acn", "bcn")

                         ## ALERT, hardcoded!
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
                             if (any(!bps %^% ss.ends))
                                 warning("Eligible SVs not matching segment ends!")

                             ## create junctions
                             junc = grl.pivot(split(bps, bps$ii))
                             toget = intersect(c("allele1", "allele2", "cn", "unknown1", "unknown2", "timing", "class"), colnames(sv))
                             values(junc) = sv[, toget, with=F]
                         } else {
                             junc = NULL
                         }

                         ## edges and graph
                         ## ALERT!! ALERT!!
                         ## TODO: still can't use addxxx() functions in a chain
                         ## doing these two steps apart will result in breakpoint missing from
                         ## self$nullGGraph()$addSegs(ss)$addJuncs(junc)
                         ## private$abEdges = self$makeAbEdges()
                         self$karyograph(tile = ss, juncs = junc, cn = TRUE)
                         return(self)
                     },

                     ## initialize from Prego result
                     prego2gGraph = function(fn){
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
                             warning("No mapping seqnames info, will throw out all non 1:24 values.")
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
                                                    } else {
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
                                            strand = "+", cn = res[[1]]$cn,
                                            left.tag = res[[1]]$node1,
                                            right.tag = res[[1]]$node2,
                                            loose=FALSE)
                         segstats = gr.fix(c(segstats, gr.strandflip(segstats)), sl)
                         neg.ix = which(strand(segstats) == "-")
                         tag1 = segstats$right.tag
                         tag1[neg.ix] = segstats$left.tag[neg.ix]
                         tag2 = segstats$left.tag
                         tag2[neg.ix] = segstats$right.tag[neg.ix]
                         ## private$segs = segstats

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
                         ## private$es = ed

                         ## create g
                         if (nrow(ed)>0){
                             g = make_directed_graph(
                                 t(as.matrix(ed[,.(from,to)])))
                         } else {
                             g = igraph::make_empty_graph(n=length(segstats))
                         }

                         private$g = g

                         ## junctions, many of them are copy 0
                         ve = data.table(res$`variant edges`)
                         if (nrow(ve)>0){
                             bp1 = dt2gr(ve[, .(seqnames = chr1, start = pos1, end = pos1)])
                             bp2 = dt2gr(ve[, .(seqnames = chr2, start = pos2, end = pos2)])
                             ## strand of breakpoint: matching left of interval, +, right, -
                             ss = gr.stripstrand(segstats %Q% (strand=="+"))
                             strand(bp1) = ifelse(is.na(match(bp1, gr.end(ss))), "+", "-")
                             strand(bp2) = ifelse(is.na(match(bp2, gr.end(ss))), "+", "-")
                             ## ALERT: don't forget to move + bp 1 nucleotide left
                             bp1 = do.call(gUtils::`%-%`, list(bp1, as.numeric(strand(bp1)=="+")))
                             bp2 = do.call(gUtils::`%-%`, list(bp2, as.numeric(strand(bp2)=="+")))
                             ## assemble the grl
                             grl = grl.pivot(GRangesList(list(bp1, bp2)))
                             mc = ve[, -c("node1", "chr1", "pos1", "node2", "chr2", "pos2"), with=F]
                             values(grl) = mc
                             ## private$junction = junctions$new(grl)
                         } else {
                             grl = junctions$new()
                         }

                         self$karyograph(tile = segstats, juncs = grl, cn=TRUE)
                         private$segs$loose = FALSE ## ALERT, no loose ends in Prego, not even implicit

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
                         cat('Edge counts:\n')
                         print(private$es[, table(type)/2])
                     },

                     ## TODO: find better default settings
                     plot = function(pad=1e3, colorful=FALSE, ...){
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
                         plot(td, win, links = private$junction)
                     },

                     ## TODO: find better default settings
                     layout = function(){
                         ## TODO: return the plot value
                         ## TODO: decide best visual parameters depend on the size of the graph!!!
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

                     ## TODO: make it informative
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

                     gGraph2gTrack = function(seg.col, ...){
                         "Create gTrack for static genome browser-style viz."
                         ## DONE: allow users to define extra fields to annotate segs or edges!!!
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
                             if (is.list(pid)) pid = unlist(pid)
                             ss$cn[lid] = ss$cn[pid]*1.2
                         }
                         ## col, border, ywid
                         ss$col = ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5))
                         ss$border = ifelse(ss$loose, ss$col, alpha("black", 0.5))
                         ss$ywid = ifelse(ss$loose, 0.001, 0.8)

                         gt = gTrack(ss, y.field="cn", edges=ed, name="CN", angle=0, ...)
                         return(gt)
                     },

                     json = function(filename='.',
                                     maxcn=100,
                                     maxweight=100,
                                     ## trim will only output seqnames
                                     ## that are relevant to the plot
                                     trim = TRUE){
                         self$gGraph2json(filename, maxcn, maxweight, trim)
                     },

                     html = function(filename='.',
                                     maxcn=100,
                                     maxweight=100,
                                     ## trim will only output seqnames that are relevant to the plot
                                     trim = TRUE){
                         if (grepl('\\.json$', filename)){
                             stop("Please refrain from naming directory with .json suffix.")
                         }

                         ## if filename not changed by user, make a ./out to dump things in
                         if (filename=="."){
                             filename = "./out"
                             system("mkdir -p ./out")
                         }
                         self$gGraph2json(filename, maxcn, maxweight, trim, all.js=TRUE)
                     },

                     gGraph2json = function(filename='.',
                                            maxcn=100,
                                            maxweight=100,
                                            ## trim will only output relevant seqnames
                                            trim = TRUE, ## ALERT trim only compatible with hg19 now
                                            all.js = FALSE,
                                            save = TRUE,
                                            ignore.strand=TRUE){
                         ## TODO: if do not ignore.strand
                         ## we will over plot pairs of edges and intervals
                         if (save){
                             if (grepl('\\.js(on)*$', filename))
                                 ## if json path was provided
                                 basedir = dirname(filename)
                             else if (filename==".") {
                                 ## default path was provided
                                 basedir = './'
                                 filename = "data.js"
                             } else {
                                 ## a directory was provided
                                 basedir = filename
                                 filename = paste(filename, 'data.json', sep = '/')
                             }

                             if (!file.exists(basedir)) {
                                 message('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }

                             if (all.js){
                                 if (!file.exists(
                                          system.file("extdata",
                                                      "gTrack.js/complete-genome-interval-graph",
                                                      package = 'gGnome')))
                                     stop("No file to copy!!")
                                 ## copy the structure of the viz system
                                 system(sprintf(
                                     'cp -r %s %s',
                                     paste0(system.file("extdata",
                                                        "gTrack.js/complete-genome-interval-graph",
                                                        package = 'gGnome'), '/*'),
                                     paste0(basedir, '/')
                                 ))
                             }
                         }

                         "Create json file for interactive visualization."
                         qw = function(x) paste0('"', x, '"') ## quote

                         ## range of CN
                         ymin=0
                         ymax=maxcn

                         ## ALERT: for a clean viz for now, only contain regular chromosomes

                         ## MARCIN EDIT: WHY IS regularChr missing in so many places? causing errors
                         ## Also should not hardcode chromosome names - makes it only applicable
                         ## to human / hg19
                         ## and if you must hardcode it, shouldn't have to be redefined so many times
                         ## what ensures that these definitions are harmonized, what if you have to change
                         ## the definition
                         regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs
                         regsegs.ix = which(as.character(seqnames(private$segs)) %in% regularChr)

                         ## processing nodes
                         ## reduce strand
                         ## remove loose nodes
                         oid = gr2dt(private$segs)[, which(strand=="+" &
                                                           loose==F &
                                                           !is.na(cn) &
                                                           seqnames %in% regularChr)]
                         ## ori ind of rev comps
                         rid = gr2dt(private$segs)[, which(strand=="-" &
                                                           loose==F &
                                                           !is.na(cn) &
                                                           seqnames %in% regularChr)]
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
                             ## keep track of gGraph node ids
                             title = paste0(seq_along(nodes), ' (', oid, '|', rid, ')'),
                             type = "interval",
                             y = pmin(maxcn, nodes$cn)
                         )

                         ## processing edges
                         ed = private$es

                         if (nrow(ed)>0){
                             ## TMPFIX: remove NA edges .. not clear where these are coming from
                             ## but likely the result of trimming / hood
                             ed = ed[!is.na(from) & !is.na(to) &
                                     from %in% regsegs.ix & to %in% regsegs.ix, ]

                             ## ALERT: bc strandlessness, I only retained half of the edges
                             ## for gwalks, we will need strandedness, so will retain everything
                             ed[,":="(soStr = as.character(strand(private$segs[from])),
                                      siStr = as.character(strand(private$segs[to])))]
                             edByType = by(ed, ed$type, function(x) x)
                             ## see which of the ab edges are "+"
                             abe = edByType$aberrant
                             ## put 3 edge types back together
                             if (is.null(edByType$loose)){
                                 ed = rbindlist(list(edByType$reference[soStr=="+"],
                                                     abe))
                             } else {
                                 ed = rbindlist(list(edByType$reference[soStr=="+"],
                                                     edByType$loose[soStr=="+"],
                                                     abe))
                             }
                             ## if encountered, switch to 0
                             ## mapping from type field to label in json
                             eType = setNames(c("REF", "ALT", "LOOSE"),
                                              c("reference", "aberrant", "loose"))
                             ## processing edges, cont.
                             fmap = node.dt[, .(oid, iid)]; setkey(fmap, oid);
                             rmap = node.dt[, .(rid, iid)]; setkey(rmap, rid);
                             ## edge data.table
                             ed.dt = ed[,.(from,
                                           to,
                                           ## source
                                           so = ifelse(soStr=="+",
                                                       node.dt[oid == from, iid],
                                                       node.dt[rid == from, iid]),
                                           ## sink
                                           si = ifelse(siStr=="+",
                                                       node.dt[oid == to, iid],
                                                       node.dt[rid == to, iid]),
                                           so.str = ifelse(soStr=="+",1,-1),
                                           si.str = ifelse(siStr=="+",1,-1),
                                           ## diff than defined in es field
                                           weight=pmin(maxweight, cn),
                                           title = "",
                                           type = eType[type]),
                                        by=1:nrow(ed)]
                             ## need to flip the negative segs to positive
                             ed.dt[, sig := ifelse(so<si, ## assuming the sorting of segs
                                                   paste0(so * so.str, '_', -si*si.str),
                                                   paste0(-si * si.str, '_', so*so.str))]
                             ed.dt[!duplicated(sig), ][, cid := seq_along(.I)]
                             ed.dt[, cid := 1:length(from)]
                             ed.dt[,":="(so = so*so.str, si = -si*si.str)]

                             ##TMPFIX: quick hack to remove dup edges
                             ed.dt = ed.dt[
                                 -which(duplicated(paste(
                                      apply(cbind(so*so.str, -si*si.str), 1,
                                            function(x) paste(sort(x), collapse = ' '))))), ]
                             ## finally, convert to JSON 'connections' string
                             connections.json = ed.dt[, paste0(
                                 c(paste0(qw("connections"),": ["),
                                   paste(
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

                         ## converting to JSON 'intervals' string
                         intervals.json = node.dt[, paste0(
                             c(paste0(qw("intervals"),": ["),
                               paste(
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
                             ## only retain important stuff
                             chrs = chrs[seqnames %in%
                                         intersect(regularChr,
                                                   unique(as.character(seqnames(private$segs))))]
                         else
                             chrs = chrs[seqnames %in% levels(seqnames(private$segs))]

                         ## converting JSON 'metadata' string
                         meta.json =
                             paste(paste0('\t',qw("metadata"),': [\n'),
                                   chrs[, paste(
                                       "\t\t{",
                                       qw("chromosome"),":", qw(seqnames),
                                       ",", qw("startPoint"),":", 1,
                                       ",", qw("endPoint"), ":", seqlengths,
                                       ",", qw("color"),
                                       ":", qw(substr(
                                                tolower(brewer.master( max(.I), 'BrBG')),
                                                1, 7)
                                               ),
                                       " }",
                                       collapse=",\n",
                                       sep="")],
                                   '\n]')

                         ## assembling the JSON
                         if (nrow(ed)>0){
                             out = paste(c("var dataInput = {",
                                           paste(
                                               c(meta.json,
                                                 intervals.json,
                                                 connections.json),
                                               collapse = ',\n'
                                           ),"}"),
                                         sep = "")
                         } else {
                             message("No edges in the graph.")
                             out = paste(c("var dataInput = {",
                                           paste(
                                               c(meta.json,
                                                 intervals.json),
                                               collapse = ',\n'
                                           ),"}"),
                                         sep = "")
                         }

                         message("Saving JSON to: ", filename)
                         writeLines(out, filename)
                         ## MARCIN COMMENT: NOT SURE WHY ANYONE WOULD NEED THE JSON BACK,
                         ## AND IT CRASHES EMACS
                         ## solution: only saves to file
                         return(self)
                     },

                     ## self-annotating functions
                     hydrogenBonds = function(){
                       ## collapse +/- strand
                       ss = unique(gr.stripstrand(private$segs))
                       idss = match(gr.stripstrand(private$segs), ss)
                       
                       ## MARCIN EDIT: fix to take care of situations where loose ends happen to exactly overlap a seg
                       ## causing error here
                       ss = paste(gr.string(gr.stripstrand(private$segs)), private$segs$loose)
                       uss = unique(ss)
                       idss = match(ss, uss)
                         if (!all(table(idss)==2)){
                             stop("Malformed object. Suggest creation again.")
                         }
                         tmpDt = data.table(ssid = seq_along(ss))
                         tmpDt[, ":="(n1 = which(idss==ssid)[1],
                                      n2 = which(idss==ssid)[2]), by=ssid]
                         hydrogenBs = tmpDt[, .(from = n1, to = n2,
                                                type="hydrogen")]
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

                     ## DONE:
                     ## if na.rm==F, balanced graph's subgraph should always be balanced!!!!!
                     subgraph = function(v=numeric(0), na.rm=T, mod=T){
                         "Given a numeric vector of vertices, \
                         change this gGraph to its subgraph consists of only these vertices."
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

                             if (length(loose.v <- intersect(v, which(private$segs$loose==T)))>0){
                                 warning("Some v is loose end. Ignore!")
                                 v = setdiff(v, loose.v)
                             }

                             ## DONE: also recover v's missing reverse complements
                             hB = self$hydrogenBonds()
                             vid = sort(unique(c(v, hB[from %in% v, to], hB[to %in% v, from])))

                             ## get the subgraph
                             newSegs = private$segs[vid]

                             newId = setNames(seq_along(vid), vid)
                             newEs = private$es[cn>0][from %in% vid & to %in% vid,
                                                      .(from=newId[as.character(from)],
                                                        to=newId[as.character(to)],
                                                        cn, type)]

                             ## DONE: use "fillin" function on the graph if na.rm=F
                             jIdx = which(grl.in(private$junction, newSegs, only=T))
                             newJuncs = private$junction[unique(jIdx)]

                             if (mod==T){
                                 private$gGraphFromScratch(segs=newSegs,
                                                           es=newEs,
                                                           junc=newJuncs,
                                                           ploidy=private$.ploidy,
                                                           purity=private$.purity)
                                 if (na.rm==F) self$fillin()
                                 return(self)
                             } else {
                                 out = gGraph$new(segs=newSegs,
                                                  es=newEs,
                                                  junctions=newJuncs,
                                                  ploidy=private$.ploidy,
                                                  purity=private$.purity)
                                 if (na.rm==F) out$fillin()
                                 return(out)
                             }
                         } else {
                             stop("Invalid input.")
                         }
                     },

                     ## DONE!!!!!!
                     ## the idea of loose end: accesorries, only exist to BALANCE the graph
                     ## make them transient
                     fillin = function(){
                         "Fill in the missing copies of edges to make the graph balanced."
                         ## GOAL: make loose ends a very free thing, add it, remove it, fuse a
                         ## pair of them or convert to a terminal feature.
                         adj = self$getAdj()
                         inSum = Matrix::colSums(adj)
                         outSum = Matrix::rowSums(adj)
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

                         if (!all(inE[, setNames(cn, toid)] == inSum) |
                             !all(outE[, setNames(cn, fromid)] == outSum))
                             stop("Adjacency not matching edges table!")

                         ## next we determine if it is feasible to fill the slacks
                         ## test if inSum>cns | outSum>cns
                         ## TODO: No!! reference edges are given copy 2!!!

                         if (any(inSum>cns | outSum>cns, na.rm = TRUE)){
                             warning("Infeasible graph!!")
                         } else {
                             colnames(inE)[2] = "cn.in"
                             colnames(outE)[2] = "cn.out"
                             ## Now fill in the loose ends
                             node.cn = merge(inE, outE, by.x="toid", by.y="fromid")
                             colnames(node.cn)[1] = "id"
                             node.cn[, cn := cns]
                             node.cn[, terminal := private$segs$terminal]
                             node.cn[, loose.out := ifelse(terminal==T & cn.out==0, 0, cn-cn.out)]
                             node.cn[, loose.in := ifelse(terminal==T & cn.in==0, 0, cn-cn.in)]

                             ## before any action, if nothing to be filled in, then stop
                             if (node.cn[, !any(loose.in>0)] | node.cn[, !any(loose.out>0)])
                                 return(self)

                             ## construct GR for new loose ends required
                             new.loose.in = node.cn[loose.in>0,
                                                    gr.start(private$segs[id], ignore.strand=FALSE)]
                             values(new.loose.in) = NULL
                             values(new.loose.in)$cn = node.cn[loose.in>0, loose.in]
                             values(new.loose.in)$loose = TRUE
                             values(new.loose.in)$terminal = TRUE

                             new.loose.out = node.cn[loose.out>0,
                                                     gr.end(private$segs[id], ignore.strand=FALSE)]
                             values(new.loose.out) = NULL
                             values(new.loose.out)$cn = node.cn[loose.out>0, loose.out]
                             values(new.loose.out)$loose = TRUE
                             values(new.loose.out)$terminal = TRUE

                             new.loose = c(new.loose.in, new.loose.out)
                             segs = private$tmpSegs = private$segs
                             old.n = length(segs)
                             node.cn[loose.out>0,
                                     new.loose.id := old.n+length(new.loose.in)+
                                         seq_along(new.loose.out)]
                             node.cn[loose.in>0,
                                     new.loose.id := old.n+seq_along(new.loose.in)]
                             ## Warning! We throw out things that was in JaBbA output!
                             ## TODO: reconstruct them on demand
                             private$segs = c(segs[,c("cn", "loose", "terminal")], new.loose)
                             newEs = rbind(node.cn[loose.in>0, .(from=new.loose.id,
                                                                 to=id,
                                                                 cn = loose.in,
                                                                 type="loose")],
                                           node.cn[loose.out>0, .(from=id,
                                                                  to=new.loose.id,
                                                                  cn = loose.out,
                                                                  type="loose")])
                             private$es = rbind(private$es[,.(from, to, cn, type)],
                                                newEs[, .(from, to, cn, type)])
                             private$g = make_directed_graph(
                                 t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))
                         }
                         return(self)
                     },

                     trim = function(gr=NULL){
                         ## DONE
                         ## if input gr is super set of private$segs, do nothing!
                         ## Only returning new obj
                         verbose = getOption("gGnome.verbose")

                         "Given a GRanges, return the trimmed subgraph overlapping it."
                         if (is.null(gr))
                             return(self)

                         gr = gr.fix(gr, get(self$refG)) ## TODO: replace the use of refG
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
                         } else {
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
                                 if (verbose) warning("Only 'cn' field is carried over.")
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
                         newSg = gGraph$new(segs=ord.nss,
                                            es=new.es)
                         return(newSg)
                     },

                     getSeqInfo = function(){
                         as.data.table(attributes(seqinfo(get(self$refG))))
                     },

                     makeAbEdges = function(){
                         "Do I even need this as a field?"
                         ## TODO: reimplement, derive junctions from edges, which is much easier
                         ## This function returns 3-d array of matching junctions to edges
                         ## DONE: derive abEdges from junction
                         if (length(private$junction)==0){
                             return(
                                 array(dim=c(0,3,2),
                                       dimnames=list(NULL,
                                                     c("from", "to", "edge.ix"),
                                                     c("+","-")))
                             )
                         } else {
                             ## based on junctions, get
                             junc = private$junction
                             ## remember, there has to be a cn field in junctions here
                             if (!"cn" %in% colnames(values(junc))){
                                 warning("'cn' not found in junction meta cols, use 1 for all.")
                                 values(junc)$cn = 1
                             }
                             jadd = which(values(junc)$cn > 0)
                             junc = junc[jadd]


                             ## TODO: why are some junctions derikved from gw 1 off??

                             abEdges = array(dim=c(length(private$junction),3,2),
                                             dimnames=list(NULL,
                                                           c("from", "to", "edge.ix"),
                                                           c("+","-")))

                             hb = self$hydrogenBonds()[, c(setNames(from, to), setNames(to, from))]
                             ## find coresponding edge.ix for abe
                             ## ASSUMPTION: junctions are width 1, marking the left nt of a bp
                             bps = grl.unlist(junc)
                             ## get one side of the edges firstn
                             seg = private$segs %Q% (loose==FALSE)

                             seg.ix = which(private$segs$loose==FALSE)

                             bps.seg = c(
                             (bps[which(strand(bps)=="-"),c("grl.ix", "grl.iix")] %**%
                              gr.end(seg[,c()]))[, c("query.id", "subject.id", "grl.ix", "grl.iix")],
                             ((bps[which(strand(bps)=="+"),c("grl.ix", "grl.iix")] %+% 1)
                                 %**%
                                 gr.start(seg[,c()]))[, c("query.id", "subject.id", "grl.ix", "grl.iix")]
                             )

                             ## discard the grl.ix that not both breakpoints match end of segments
                             ## also, the ones with cn <= 0
                             jIn = sort(as.numeric(names(which(
                                 table(bps.seg$grl.ix)==2
                             ))))


                             ## MARCIN COMMENT: why would there be cn < 0
                             if ("cn" %in% colnames(values(junc))){
                                 jIn = setdiff(jIn, which(values(junc)$cn <= 0))
                             }
                             bps.seg = bps.seg %Q% (grl.ix %in% jIn)
                             bps.seg = bps.seg %Q% (order(grl.ix, grl.iix))

                             to.node = bps.seg$subject.id
                             from.node = hb[as.character(to.node)]
                             ix1 = which(bps.seg$grl.iix == 1)
                             ix2 = which(bps.seg$grl.iix == 2)
                             ed1 = data.table(from = from.node[ix1],
                                              to = to.node[ix2])
                             ed2 = data.table(from = from.node[ix2],
                                              to = to.node[ix1])
                             jeid = c(ed1[, paste(from, to)], ed2[, paste(from, to)])


                             eids = private$es[, paste(from, to)]
                             abEdges[jadd[jIn],1:2,"+"] = as.matrix(ed1[, .(from, to)])
                             abEdges[jadd[jIn],1:2,"-"] = as.matrix(ed2[, .(from, to)])

                             ## MARCIN COMMENT: struggling to understand the reason for these expressions inside ed1
                             ## ie there is nothing in the expression inside ed1 that is accessing any elements of ed1
                             ## furthermore
                             abEdges[jadd[jIn],3,"+"] = ed1[, which(eids %in% jeid[1:nrow(ed1)] & private$es$type=="aberrant")]
                             abEdges[jadd[jIn],3,"-"] = ed2[, which(eids %in% jeid[(nrow(ed1)+1):length(jeid)] & private$es$type=="aberrant")]

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
                                 ss = grbind(c(private$segs[private$segs$loose == FALSE, c()],
                                               win[, c()]))

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
                             gr1 = c(gr1, gr.strandflip(gr1[ix]))
                         }

                         if (any(ix <- strand(gr2)=='*'))
                         {
                             strand(gr2)[ix] = '+'
                             gr2 = c(gr2, gr.strandflip(gr2[ix]))
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
                     ## NOW TODO
                     proximity = function(query, subject,
                                          verbose=F, mc.cores=1,
                                          max.dist=1e6){

                         ## TODO:
                         adj = self$getAdj()
                         ix = which(adj[private$abEdges[,1:2,1]]>0)
                         if (length(ix)>0) {
                             ra1 = gr.strandflip(
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
                         ## TODO: make karyograph output compatible with Marcin's!!! THis is important and hard...
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
                                 get.shortest.paths(kg$G,
                                                    vix.query[x, 'end'],
                                                    vix.subject[y, 'start'],
                                                    weights = E(kg$G)$weight,
                                                    mode = 'out')$vpath[[1]]
                             else if ((ra.which[x, y]) == 2)
                                 rev(get.shortest.paths(kg$G,
                                                        vix.query[x, 'start'],
                                                        vix.subject[y, 'end'],
                                                        weights = E(kg$G)$weight,
                                                        mode = 'in')$vpath[[1]])
                             else if ((ra.which[x, y]) == 3)
                                 get.shortest.paths(kg$G,
                                                    vix.query[x, 'end.n'],
                                                    vix.subject[y, 'start'],
                                                    weights = E(kg$G)$weight,
                                                    mode = 'out')$vpath[[1]]
                             else if ((ra.which[x, y]) == 4)
                                 rev(get.shortest.paths(kg$G,
                                                        vix.query[x, 'start.n'],
                                                        vix.subject[y, 'end'],
                                                        weights = E(kg$G)$weight,
                                                        mode = 'in')$vpath[[1]])
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

                     },

                     e2j = function(etype="aberrant"){
                         "Return the junctions based on edges in this graph."
                         verbose = getOption("gGnome.verbose")
                         ## if (!"type" %in% colnames(private$es)){
                         ##     private$es = etype(private$segs, private$es)
                         ## } else if (any(!private$es$type %in% c("reference", "aberrant", "loose"))){
                         ##     private$es = etype(private$segs, private$es, force=T)
                         ## }

                         private$es = etype(private$segs, private$es, force=T)
                         es = private$es
                         es[, eix := 1:.N]

                         if (etype=="all") etype = c("aberrant", "reference", "loose")

                         abe = es[type %in% etype]
                         if (nrow(abe)==0){
                             empty.out = private$junction = junctions()
                             return(empty.out)
                         }

                         ## MOMENT
                         ## browser()
                         if (any(! c("fromChr", "fromStr", "fromStart", "fromEnd",
                                     "toChr", "toStr", "toStart", "toEnd")
                                 %in% colnames(abe))){

                             if (verbose)
                                 message("Redo the important metadata gathering.")

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
                             hb = self$hydrogenBonds()
                             hb.map = hb[, c(setNames(from, to),
                                             setNames(to, from))]
                             abe[, ":="(eid = paste(from, to),
                                        reid = paste(hb.map[as.character(to)],
                                                     hb.map[as.character(from)]))]
                         }

                         abe[, ":="(ix = 1:.N,
                                    rix = match(reid, eid))]
                         abe[, unique.ix := ifelse(rix>=ix, paste(ix, rix), paste(rix, ix))]
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

                         bp2 = dt2gr(jdt[, .(seqnames = fromChr,
                                             strand = toStr,
                                             start = ifelse(toStr=="+", toStart-1, toEnd),
                                             end = ifelse(fromStr=="+", toStart-1, toEnd),
                                             eclass)])

                         junc = junctions(grl.pivot(GRangesList(bp1, bp2)))
                         values(junc)$eclass = bp1$eclass
                         values(junc)$type = abe[.(values(junc)$eclass, 1), type]
                         values(junc)$from1 = abe[.(values(junc)$eclass, 1), from]
                         values(junc)$to1 = abe[.(values(junc)$eclass, 1), to]
                         values(junc)$from2 = abe[.(values(junc)$eclass, 2), from]
                         values(junc)$to2 = abe[.(values(junc)$eclass, 2), to]

                         private$junction = junc
                         return(junc)
                     },

                     jGraph = function(){
                         ##TODO: migrate the jGraph function here

                     },

                     ## property constraints
                     isJunctionBalanced = function(){
                         ## ALERT: this is too loose!!!
                         ## TODO: redo this function!!!
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
                         tCsum = Matrix::colSums(adj)[validTerminal]
                         tRsum = Matrix::rowSums(adj)[validTerminal]
                         terminalConSide = ifelse(tCsum==0, tRsum, tCsum)
                         terminalTrue = terminalConSide == private$segstats[validTerminal]$cn
                         return(all(middleTrue) & all(terminalTrue))
                     },

                     ## isDoubleStrand = function(){
                     ##     ## DONE: test if segs come in +/- pairs
                     ##     identical((ss %Q% (strand=="-"))[, c()],
                     ##               gr.strandflip(ss %Q% (strand=="+"))[, c()])
                     ## },

                     getLooseEnds = function(){
                         ## TODO: return all loose ends as a GRanges
                     },

                     walk = function(v = numeric(0), j = numeric(0), peel=FALSE){
                         "Generate gWalks object given by node or edge sequences."
                         ## TODO: if given j, override v
                         ## MOMENT
                         ## test validity (existence and CN) of path defined by j
                         ## if passed, convert to v and recurse
                         browser()
                     },

                     chromoplexy = function(pad = 1e3){
                         "Identifying parts of the graph that are probably produced from chromoplexy events. In gGraph class the method ignores information from CN."
                         
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
                     ## igraph obj representing the graph structure
                     g = NULL,
                     ## temporary segs for backtrace when modifying segs
                     tmpSegs = NULL,
                     ## putative junctions, junctions
                     junction = NULL,
                     abEdges = NULL,
                     ## ploidy is set to 2, only to init null graph,
                     ## otherwise inferred from segs
                     .ploidy = NULL,
                     ## tumor cell proportion
                     .purity = NULL,
                     ## the partition result of 'g'
                     partition = NULL,

                     ## ----- private methods
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
                         "Nodes as GRanges, edges as data.frame or adj matrix."
                         if (!is.null(names(segs))){
                             if (any(duplicated(names(segs)))) names(segs) = NULL
                         }

                         if (length(segs)==0)
                             return(self)

                         ## ALERT: if no "loose" col in segs, default to FALSE
                         if (!"loose" %in% colnames(values(segs))) segs$loose=FALSE

                         private$segs = segs

                         hB = self$hydrogenBonds()
                         map = hB[, c(setNames(from, to), setNames(to, from))]
                         hB[, tile.id := 1:.N]
                         tile.id = c(hB[, setNames(tile.id, from)],
                                     hB[, setNames(tile.id, to)])
                         private$segs$tile.id = tile.id[as.character(seq_along(private$segs))]

                         ## check es input
                         if (is(es, "data.frame")){
                             if (!all(c("from", "to") %in% colnames(es)))
                                 stop("Given edge data must have 'from' and 'to' fields.")

                             es = es[from %in% seq_along(segs) & to %in% seq_along(segs)]
                         } else if (is(es, "matrix") | is(es, "Matrix")){
                             A = es
                             if (!all(dim(A)==length(segs)))
                                 stop("Given adjacency matrix in wrong dimension.")

                             ## we allow numeric or logical values
                             if (is(A[1,1], "numeric")){
                                 ## when it's numeric & non-negative, the value is copy number
                                 if (any(A<0)){
                                     A = A>0
                                 } else {
                                     es = which(A>0, arr.ind=T)
                                     colnames(es) = c("from", "to")
                                     cn = A[es]
                                     es = as.data.table(es)
                                     es[, cn := cn]
                                 }
                             }

                             if (is(A[1,1], "logical")){
                                 es = which(A==TRUE, arr.ind=T)
                                 colnames(es) = c("from", "to")
                                 es = as.data.table(es)
                                 es[, cn := 0]
                             }
                         }

                         ## finished converting adj matrix to data.table
                         ## now check if it is skew-symmetric
                         es[, eid := paste(from, to)]
                         es[, reid := paste(map[as.character(to)],
                                            map[as.character(from)])]
                         ematch = es[, match(eid, reid)]

                         ## ALERT: sometimes there are NAs in ematch!!!
                         ## TODO: how to deal with NA in ematch???
                         if (all(es[ematch, cn]==es[, cn], na.rm=T)){
                             if (as.logical(getOption("gGnome.verbose"))){
                                 message("Edge copies balanced!")
                             }
                         } else {
                             ## TODO: maybe don't try to do too much????
                             ## or help the user with this????
                             stop("Given edge data is not skew-symmetric!!!")
                         }

                         ## when "type" is missing, infer it
                         if (!is.element("type", colnames(es))){
                             es = etype(private$segs, es)
                         } else if (any(!es$type %in% c("reference", "aberrant", "loose"))){
                             es = etype(private$segs, es, force=TRUE)
                         }

                         private$es = es

                         if (private$es[, any(type=="aberrant")]){
                             abEs = private$es[type=="aberrant"]

                             ## TODO: what if no aberrant edge is here?
                             abEs[, ":="(tmp.id = paste(from, to, sep="-"),
                                         tmp.id.r = paste(map[as.character(to)],
                                                          map[as.character(from)],
                                                          sep="-"))]
                             abEs[, first.id := sort(c(tmp.id, tmp.id.r))[1], by=1:nrow(abEs)]
                             ## each junction shows up twice now
                             abEs = abEs[!duplicated(first.id)]
                         }

                         ## relabel the terminals!
                         ## whichTerminal = private$es[, setxor(from, to)]
                         ## private$segs$terminal = seq_along(private$segs) %in% whichTerminal
                         private$g = make_directed_graph(
                             t(as.matrix(private$es[,.(from,to)])), n=length(private$segs))

                         ## if (is.null(junc) & exists("abEs")){
                         ##     warning("Junctions not provided. Inferring from edges.")
                         ##     bp.from =
                         ##         gr.end(private$segs[abEs$from], ignore.strand=FALSE)[,c()]
                         ##     bp.from[which(strand(bp.from)=="-")] =
                         ##         bp.from[which(strand(bp.from)=="-")] %-% 1
                         ##     bp.from = gr.strandflip(bp.from)
                         ##     ## the other end
                         ##     bp.to =
                         ##         gr.start(private$segs[abEs$to], ignore.strand=FALSE)[,c()]
                         ##     bp.to[which(strand(bp.from)=="+")] =
                         ##         bp.to[which(strand(bp.from)=="+")] %-% 1

                         ##     junc = grl.pivot(GRangesList(bp.from, bp.to))
                         ## }
                         ## if (is.null(junc)){
                         ##     private$junction = new("junctions", GRangesList())
                         ## } else if (is(junc, "junctions")){
                         ##     private$junction = junc
                         ## } else if (is(junc, "GRangesList")){
                         ##     private$junction = new("junctions", junc)
                         ## } else if (is.null){
                         ##     stop("Junctions malformed!")
                         ## }

                         ## MARCIN EDIT: changed to allow failures of makeAbEdges (useful for creating gwalks
                         ## with overlapping segments)

                         private$abEdges = tryCatch(self$makeAbEdges(), error = function(e) NULL)
                         private$.ploidy = ploidy
                         private$.purity = purity
                         return(self)
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
                         if (is.null(private$junction))
                             self$e2j()
                         return(private$junction)
                     },
                     igPlot = function(){
                         ## DONE: make igraph plot
                         return(self$layout())
                     },
                     ## ALERT:
                     ## not much sense to use active binding, deprecated for now
                     ## json = function(file='~/public_html/gGraph'){
                     ##     return(self$gGraph2json(filename=file))
                     ## },
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
                     },
                     purity = function(){
                         return(private$.purity)
                     },
                     ploidy = function(){
                         return(private$.ploidy)
                     }
                 )
                 )

## ============= S3 generics of gGraph ============= ##
#'
#'
#'
components <- function (x, ...) {
    UseMethod("components", x)
}
components.igraph <- function(iGraph){
    return(igraph::components(iGraph))
}

#' components
#' strongly connected components, returned as a list of gGraph objects
#'
#' @return a list of gGraph objects representing each partition of the input
components.gGraph <- function(gGraph){
    ## input must be a gGraph!
    if (!is(gGraph, "gGraph")){
        stop("Invalid input.")
    }
    return(gGraph$components())
}

#' @name length
#'
#' @description return the number of strongly connected components of the graph
#' @export
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

#' @name %+%
#'
`%+%.gGraph` <- function(gg1, gg2){
    return(gg1$add(gg2))
}

## ============= exported functions of gGraph ============= ##
setAs("gGraph", "bGraph",
      function(from){
          return(bGraph$new(from))
      })




#' @name gread
#' Parse the outputs from rearrangement graph callers.
#'
#' @param file filename to JaBbA's rds, PREGO's intervalFile, or Weaver's output directory
#' @export
gread = function(file){
    verbose = getOption("gGnome.verbose")

    if (is.list(file)){
        if (all(is.element(c("segstats", "adj", "ab.edges", "edges", "G",
                             "td", "purity", "ploidy", "junctions"),
                           names(jabba))))
            jabba = jabba
    }
    ## MOMENT
    ## decide what output this is
    if (!file.exists(file)) stop("No such file or directory!")

    if (dir.exists(file)){
        if (verbose)
            message("Given a directory, assume it's Weaver.")
        return(gGraph$new(weaver=file))
    } else if (grepl(".rds$", file)){
        if (verbose)
            message("Try reading the RDS.")

        rds = tryCatch(readRDS(file),
                       error=function(e)
                           stop("Given file can't be read as RDS."))

        if (is(rds, "gGraph")) {
            return(rds)
        } else if (is(rds, "list")){
            jab = bGraph$new(jabba = file)
        }
    } else {
        ## prego = bGraph$new(prego = file)
        prego = gGraph$new(prego = file)
    }
}

#' @name gwrite
#' Write the data in various formats into file.
#'
#' @param file filename to JaBbA's rds, PREGO's intervalFile, or Weaver's output directory
#' @param format can be any one of c("json", "gfa", "vcf")
#' 
#'
#' @export
gwrite = function(obj, filename=paste("foo", format, sep="."), format = "json"){
    verbose = getOption("gGnome.verbose")

    ## MOMENT
    ## decide what output format this is
    if (format==json) {
        
    }
}
##############################
## bGraph
##############################
#' Descendant of gGraph class, where junction balance restraint must be met at all times
#'
#' @import R6
#' @import Matrix
#'
#' @export
bGraph = R6Class("bGraph",
                 inherit = gGraph,
                 public = list(
                     ## overwrite constructor: restrict about junction balance
                     initialize = function(gG=NULL, jabba=NULL, prego=NULL){
                         if (!is.null(gG)){
                             if (is(gG, "gGraph")){
                                 ## TODO
                                 private$gGraphFromScratch(gG$segstats, gG$edges, gG$junctions, gG$ploidy, gG$purity)
                                 return(self)
                             } else {
                                 stop("Invalid input gG.")
                             }
                         } else if (!is.null(jabba)) {
                             ## MARCIN EDIT: this will break if jabba is not a character (ie on file.exists)
                             ##  if (is.character(jabba) & file.exists(jabba)) jabba = readRDS(jabba)
                             if (is.character(jabba))
                             {
                                 if (file.exists(jabba))
                                     jabba = readRDS(jabba)
                                 else
                                     stop(paste('file', jabba, 'not found'))
                             }

                             ## allRegChr = all(
                             ##     as.vector(seqnames(unlist(jabba$junctions))) %in% regularChr
                             ## )
                             self$jabba2gGraph(jabba=jabba)
                             if (self$isJunctionBalanced()){
                                 return(self)
                             } else {
                                 stop("Invalid input gG.")
                             }
                         } else if (!is.null(prego)){
                             if (is.character(prego))
                             {
                                 if (file.exists(prego))
                                     self$prego2gGraph(fn=prego)
                                 else
                                     stop(paste('file', jabba, 'not found'))
                             }

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
                             out = bGraph$new(out)
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
                         "Enumerate all the possible multiset of walks that can be represented by this graph."
                         ## TODO: something's wrong here, need redo
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

                         ## TODO: convert karyoMIP solution to gWalks object, with the new gw definition
                         ## is.cyc = Matrix::colSums(K[h$etype == 'slack', ])==0 & Matrix::colSums((Bc %*% K)!=0)==0
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
                         gw = gWalks$new(segs=segs[whichSeg], paths=p$paths, is.cycle=p$is.cyc, cn = p$cn)

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
                     walk2 = function(verbose = FALSE, grl=TRUE){
                         ## TODO: how come the cn.adj have NA values while none when initialized?
                         cn.adj = self$getAdj()
                         adj = as.matrix(cn.adj)
                         adj.new = adj*0
                         adj[which(adj!=0, arr.ind = TRUE)] = width(private$segs)[which(adj!=0, arr.ind = TRUE)[,2]] ## make all edges a large number by default
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
                         ss = gr2dt(private$segs)[ , vid:= 1:length(seqnames)]
                         ss[loose == TRUE, is.end := TRUE]
                         ss[loose == FALSE, is.end := 1:length(loose) %in% c(which.min(start), which.max(end)), by = list(seqnames, strand)]
                         ends = which(ss$is.end)

                         ## sanity check
                         unb = which(!ss$is.end &
                                     Matrix::rowSums(self$getAdj(), na.rm = TRUE) !=
                                     Matrix::colSums(self$getAdj(), na.rm = TRUE))

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
                         ## MARK
                         palindromic.path = rep(FALSE, maxrow)
                         palindromic.cycle = rep(FALSE, maxrow)

                         #' first peel off "simple" paths i.e. zero degree
                         #' ends with >0 copy number
                         psimp =  which(degree(G, mode = 'out')==0 &
                                        degree(G, mode = 'in')==0 &
                                        private$segs$cn>0)
                         i = 0
                         if (length(psimp)>0)
                         {
                             vpaths[1:length(psimp)] = split(psimp, 1:length(psimp))
                             ## there is no "edge" associated with a zero total degree node
                             epaths[1:length(psimp)] = lapply(psimp, function(x) cbind(NA, NA))
                             cns[1:length(psimp)] = private$segs$cn[psimp]
                             i = length(psimp)
                         }

                         ## now iterate from shortest to longest path
                         ## peel that path off and see if it is still there ..
                         ## and see if it is still there
                         ## peel off top path and add to stack, then update cn.adj

                         if (!"tile.id" %in% colnames(values(private$segs))){
                             warning("Creating tile.id on the fly.")
                             ## ALERT!! TODO!! Should we always have tile.id???
                             ## What to do when there is no tile.id?
                         }

                         tile.map = gr2dt(private$segs)[, .(id = 1:length(tile.id),
                                                            tile.id = ifelse(strand == '+',
                                                                             1, -1)*tile.id)]
                         rtile.map = gr2dt(private$segs)[, .(id = 1:length(tile.id),
                                                             tile.id = ifelse(strand == '+',
                                                                              1, -1)*tile.id)]
                         setkey(tile.map, id)
                         setkey(rtile.map, tile.id)

                         ## unique pair of edge ids:
                         ## rev comp of a foldback edge will be identical to itself!!!
                         ed = data.table(private$es)[cn>0, .(from, to , cn)]
                         ed[, ":="(fromss = tile.map[ .(from), tile.id],
                                   toss = tile.map[ .(to), tile.id]),
                            by = 1:nrow(ed)]
                         ed[, weight :=  adj[cbind(from, to)]]

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
                                 message('Path peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left and ', nrow(ij), ' ends to resolve' )
                             i = i+1
                             p = get.constrained.shortest.path(cn.adj, G, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

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
                                 cns[i] = ed[.(eids), if (length(cn)>1) cn/2 else cn, by = eclass][, floor(min(V1))] ## update cn correctly, adjusting constraints for palindromic edges by 1/2


                                 rvpath = rtile.map[list(tile.map[list(vpaths[[i]]), -rev(tile.id)]), id]
                                 repath = cbind(rvpath[-length(rvpath)], rvpath[-1])
                                 plen = length(rvpath)
                                 hplen = floor(length(rvpath)/2)

                                 ## (awkward) check for palindromicity for odd and even length palindromes
                                 ## if (all((vpaths[[i]]==rvpath)[c(1:hplen,(plen-hplen+1):plen)]))
                                 if (ed[eids, any(table(eclass)>1)])
                                     palindromic.path[i] = TRUE
                                 ## else
                                 ## {
                                 vpaths[[i+1]] = rvpath
                                 epaths[[i+1]] = repath
                                 cns[i+1] = cns[i]
                                 palindromic.path[i+1] = TRUE
                                 ## }
                                 ##        palindromic = TRUE ## set to true while we "figure things out"


                                        # so now we want to subtract that cn units of that path from the graph
                                        # so we want to update the current adjacency matrix to remove that path
                                        # while keeping track of of the paths on the stack
                                 cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]-cns[i]

                                 ## if (!palindromic) ## update reverse complement unless palindromic
                                 cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]-cns[i+1]

                                 if (!all(cn.adj[epaths[[i]]]>=0)) ## something wrong, backtrack
                                 {
                                     message('backtracking ...') ## maybe we got stuck in a quasi-palindrome and need to backtrack
                                        #            browser()
                                     cn.adj[epaths[[i]]] = cn.adj[epaths[[i]]]+cns[i]
                                     ## if (!palindromic) ## update reverse complement unless palindromic
                                     cn.adj[epaths[[i+1]]] = cn.adj[epaths[[i+1]]]+cns[i+1]
                                     i = i-1
                                     ij = ij[-1, , drop = FALSE]
                                 }
                                 else ## continue, reduce
                                 {
                                     adj.new[epaths[[i]]] = adj.new[epaths[[i]]] + cns[i]
                                     ## if (!palindromic)
                                     adj.new[epaths[[i+1]]] = adj.new[epaths[[i+1]]] + cns[i]

                                     ## ## make sure I didn't overuse any edge
                                     ## if (nrow(overdue <- which((as.matrix(jab$adj)-adj.new)<0, arr.ind=T))>0) {
                                     ##     print("Edge copy deficit!")
                                     ##     browser()
                                     ## }

                                     ## intermediate check
                                     ## if (length(which(((adj.new + cn.adj) - jab$adj)!=0, arr.ind = TRUE)))
                                     ##     browser()

                                     to.rm = epaths[[i]][which(cn.adj[epaths[[i]]]==0), ,drop = FALSE]
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
                                         ## if (any((colSums(cn.adj)*rowSums(cn.adj) != 0) & (colSums(cn.adj) != rowSums(cn.adj)))){
                                         ##     print("Junction OUT OF BALANCE!")
                                         ##     browser()
                                         ## }

                                         ## ## should be no new ends
                                         ## if (length(new.ends)>0){
                                         ##     print("Please, no new ends!")
                                         ##     browser()
                                         ## }

                                         ## remain = as.matrix(jab$adj) - adj.new
                                         ## nb <- which(colSums(remain) != rowSums(remain))
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
                             seg.ix = which(strand(private$segs)=='+');
                             seg.rix = which(strand(private$segs)=='-');


                             if (nrow(ij)==0 & cleanup_mode == FALSE)
                             {
                                 message('!!!!!!!!!!!!!!!!!!!!!!!!!!STARTING CLEANUP MODE FOR PATHS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
                                 ij = as.data.table(which(!is.infinite(D), arr.ind = TRUE))[, dist := D[cbind(row, col)]][row != col, ][order(dist), ][, row := ends[row]][, col := ends[col]]
                                 cleanup_mode = TRUE
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

                         csimp = which(Matrix::diag(cn.adj)!=0)
                         i = 0
                         if (length(csimp)>0)
                         {
                             vcycles[1:length(csimp)] = split(csimp, 1:length(csimp))
                             ecycles[1:length(csimp)] = lapply(csimp, function(x) cbind(x, x))
                             ccns[1:length(csimp)] = Matrix::diag(cn.adj)[csimp]
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

                         ## peel off top cycles and add to stack, then update cn.adj
                         while (nrow(ij)>0)
                         {
                             if (verbose)
                                 message('Cycle-peeling iteration ', i, ' with ', sum(adj!=0, na.rm = TRUE), ' edges left ', nrow(ij) )
                             i = i+1
                             p = get.constrained.shortest.path(cn.adj, G, allD = D, v = ij[1, 1], to = ij[1, 2], weight = E(G)$weight, edges = ed, verbose = TRUE, mip = cleanup_mode)

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

                                        # so now we want to subtract that cn units of that path from the graph
                                        # so we want to update the current adjacency matrix to remove that path
                                        # while keeping track of of the cycles on the stack
                                 cn.adj[ecycles[[i]]] = cn.adj[ecycles[[i]]]-ccns[i]
                                 ## if (!palindromic) ## update reverse complement unless palindromic
                                 cn.adj[ecycles[[i+1]]] = cn.adj[ecycles[[i+1]]]-ccns[i+1]

                                 ## ALERT!!!!!!! There are NAs in cn.adj at this point
                                 ## TODO!! fix it, temporarily ignored!!!
                                 ## if (!all(cn.adj[ecycles[[i]]]>=0, na.rm=TRUE))
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

                                         ## if (any((colSums(cn.adj)*rowSums(cn.adj) != 0) & (colSums(cn.adj) != rowSums(cn.adj)))){
                                         ##     print("Junction OUT OF BALANCE!")
                                         ##     browser()
                                         ## }

                                         ## remain = as.matrix(jab$adj) - adj.new
                                         ## nb <- which(colSums(remain) != rowSums(remain))
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

                         tmp = cbind(do.call(rbind, eall), rep(ecn, sapply(eall, nrow)), munlist(eall))
                         ix = which(rowSums(is.na(tmp[, 1:2]))==0)

                         if (length(ix)>0)
                             adj.new = sparseMatrix(tmp[ix,1], tmp[ix,2], x = tmp[ix,3], dims = dim(adj))
                         else
                             adj.new = sparseMatrix(1, 1, x = 0, dims = dim(adj))
                         vix = munlist(vall) ## here is the node indices
                         paths = split(private$segs[vix[,3]], vix[,1])

                         values(paths)$ogid = 1:length(paths)
                         values(paths)$cn = ecn[as.numeric(names(paths))]
                         values(paths)$label = paste('CN=', ecn[as.numeric(names(paths))], sep = '')
                         values(paths)$is.cycle = !(as.numeric(names(paths)) %in% 1:length(vpaths))
                         values(paths)$numsegs = elementNROWS(paths)
                         values(paths)$wid = sapply(lapply(paths, width), sum)

                         check = which((adj.new - self$getAdj()) !=0, arr.ind = TRUE)
                         if (length(check)>0)
                             stop('Alleles do not add up to marginal copy number profile!')
                         else if (verbose)
                             message('Cross check successful: sum of walk copy numbers = marginal JaBbA edge set!')

                         ## match up paths and their reverse complements
                         ## TODO: fix this matching
                         psig = lapply(paths,
                                       function(x) {
                                           ifelse(as.logical(strand(x)=='+'), 1, -1)*x$tile.id
                                       })
                         psig.flip = sapply(psig, function(x) -rev(x))

                         unmix = data.table(
                             ix = 1:length(paths),
                             mix = match(sapply(psig, paste, collapse = ','),
                                         sapply(psig.flip, paste, collapse = ','))
                         )[, pos := 1:length(mix)<mix][order(!pos), ]
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

                         ## DONE: construct gWalks object as output

                         if (grl)
                             return(paths)
                         else
                         {
                             ## EDITS BY MARCIN
                             ## simplify grl before sending to gwalks
                             ## i.e. collapse reference adjacent intervals into single intervals among walks / paths
                             ## tmp.dt = as.data.table(paths)[, pid := group_name][, nix := 1:.N, by =pid]
                             ## setkeyv(tmp.dt, c('pid', 'nix'))

                             ## ## mark nodes that precede a reference junction
                             ## tmp.dt[, d.to.next := c((start-shift(end))[-1], NA), by = pid]
                             ## tmp.dt[, d.to.next.neg := c((shift(start)-end)[-1], NA), by = pid]
                             ## tmp.dt[, same.strand := c((strand==shift(strand))[-1], NA), by = pid]
                             ## tmp.dt[, same.chrom := c((as.character(seqnames)==shift(as.character(seqnames)))[-1], NA), by = pid]
                             ## tmp.dt[, last.node := 1:.N == .N, by = pid]
                             ## tmp.dt[, before.ref :=
                             ##              (((d.to.next<=1 & d.to.next>=0 & strand == '+') |
                             ##                (d.to.next.neg<=1 & d.to.next.neg>=0 & strand == '-')
                             ##              ) & same.strand & same.chrom)]
                             ## tmp.dt[is.na(before.ref), before.ref := FALSE]

                             ## ## label reference runs of nodes then collapse
                             ## .labrun = function(x) ifelse(x, cumsum(diff(as.numeric(c(FALSE, x)))>0), as.integer(NA))
                             ## tmp.dt[, ref.run := .labrun(before.ref), by = pid]
                             ## tmp.dt[, ref.run.last := shift(ref.run), by = pid]
                             ## tmp.dt[is.na(ref.run) & !is.na(ref.run.last), ref.run := ref.run.last]
                             ## tmp.dt[!is.na(ref.run), ref.run.id := paste(pid, ref.run)]
                             ## collapsed.dt = tmp.dt[!is.na(ref.run.id), .(nid = paste(nid, collapse = ' '),
                             ##                                             nix = nix[1],
                             ##                                             pid = pid[1],
                             ##                                             seqnames = seqnames[1],
                             ##                                             start = min(start),
                             ##                                             end = max(end),
                             ##                                             strand = strand[1]
                             ##                                             ), by = ref.run.id]

                             ## ## concatenate back with nodes that precede a non reference junction
                             ## tmp.dt = rbind(tmp.dt[is.na(ref.run.id), .(pid, nid = as.character(nid), nix, seqnames, start, end, strand)],
                             ##                collapsed.dt[, .(pid, nid, nix, seqnames, start, end, strand)])
                             ## setkeyv(tmp.dt, c('pid', 'nix'))

                             ## tmp.gr = dt2gr(tmp.dt)
                             ## tmp.segs = unique(tmp.gr)
                             ## tmp.gr$seg.id = match(tmp.gr, tmp.segs)
                             ## tmp.paths = split(tmp.gr$seg.id, tmp.gr$pid)
                             ## tmp.vals = as.data.frame(values(paths[names(tmp.paths)]))
                             gw = gWalks$new(grl = paths)
                             return(gw)
                         }
                     }
                 ),
                 private = list(

                 ),
                 active = list())

setAs("gGraph", "bGraph",
      function(from){
          return(bGraph$new(gG = from))
      })

setAs("gWalks", "bGraph",
      function(from){
          return(as(from$gw2gg(), "bGraph"))
      })

## Utilities
ul = function(x, n=6){
    n = pmin(pmin(dim(x)), n)
    return(x[1:n, 1:n])
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
                                         mip = TRUE
                                         )
{

    if (is.null(allD)) allD = shortest.paths(G, mode="out", weights = weight)

    v = as.numeric(v)
    to = as.numeric(to)

    if (is.infinite(allD[v, to]) | allD[v, to]==0) return(NULL)

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
        if (verbose)
            message('Shortest path is good enough!')
        return(tmp.p)
    }

    if (!mip)
        return(NULL)


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
    c = edges$weight

    res = Rcplex(c, rbind(A[ix,], R), c(b[ix], Rb), sense = c(rep('E', length(ix)), rep('L', length(Rb))),
                 lb = 0, vtype = "B",
                 objsense = 'min')

    if (verbose)
        message('YES WE ARE DOING PROPER MIP!!!!')

    if (res$status!=101)
    {
        if (verbose)
            message('No solution to MIP!')

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
        if (length(overdrafts.eclass)==0)
            message('No overdrafts after MIP')
        else
        {
            message('Still overdraft!')
            browser()
        }
    }

                                        #    browser()
    return(tmp.p)
}

#' gtf2json
#' Turning a GTF format gene annotation into JSON
#'
#' @export
#'
gtf2json = function(gtf=NULL, gtf.rds=NULL, gtf.gr.rds=NULL, filename="./gtf.json",
                    genes=NULL, grep=NULL, grepe=NULL, chrom.sizes=NULL, include.chr=NULL,
                    gene.collapse=TRUE, verbose = TRUE){
    require(data.table)
    require(gUtils)
    if (!is.null(gtf.gr.rds)){
        message("Using GRanges from rds file.")
        infile = gtf.gr.rds
        gr = readRDS(gtf.gr.rds)
        dt = gr2dt(gr)
    } else if (!is.null(gtf.rds)){
        message("Using GTF data.table from rds file.")
        infile = gtf.rds
        dt = as.data.table(readRDS(gtf.rds))
    } else if (!is.null(gtf)){
        message("Using raw GTF file.")
        infile = gtf
        dt = fread(gtf)
        dt = dt[, .(seqnames = V1, start = V4, end = V5,
                    strand = V7, type = V3, tosp = V9)]

        ## split metadata columns
        tosp = strsplit(dt$tosp, ";")

        gene_id = gsub("\"", "",
                       gsub("gene_id \"", "",
                            sapply(tosp, grep, pattern="gene_id", value=T)))
        gene_name = gsub("\"", "",
                         gsub("gene_name \"", "",
                              sapply(tosp, grep, pattern="gene_name", value=T)))
        gene_type = gsub("\"", "",
                         gsub("gene_type \"", "",
                              sapply(tosp, grep, pattern="gene_type", value=T)))
        transcript_id = gsub("\"", "",
                             gsub("transcript_id \"", "",
                                  sapply(tosp, grep, pattern="transcript_id", value=T)))
        transcript_name = gsub("\"", "",
                               gsub("transcript_name \"", "",
                                    sapply(tosp, grep, pattern="transcript_name", value=T)))

        dt = dt[, .(.SD, gene_id = gene_id, gene_name = gene_name, gene_type = gene_type,
                    transcript_id = transcript_id, transcript_name = transcript_name)]
    } else {
        warning("No input gene annotation. Use the built-in GENCODE v19 in gUtils package")
        require(skidb)
        gr = read_gencode()
        infile = "default"
        dt = gr2dt(gr)
    }

    if (verbose) message("Finished reading raw data, start processing.")

    ## get seqlengths
    if (is.null(chrom.sizes)){
        message("No ref genome seqlengths given, use default.")
        ## chrom.sizes = system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
        ## system.file("extdata", "hg19.regularChr.chrom.sizes", package="gGnome")
        Sys.setenv(DEFAULT_BSGENOME=system.file("extdata", "hg19.regularChr.chrom.sizes", package="gUtils"))
    }
    sl = hg_seqlengths(include.junk=TRUE)
    if (!is.null(include.chr)){
        sl = sl[include.chr]
    }
    chrs = data.table(seqnames = names(sl), seqlengths=sl)

    ## meta data field
    require(RColorBrewer)
    qw = function(x) paste0('"', x, '"') ## quote

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

    if (verbose) message("Metadata fields done.")

    ## reduce columns: seqnames, start, end, strand, type, gene_id, gene_name, gene_type, transcript_id
    ## reduce rows: gene_status, "KNOWN"; gene_type, not "pseudo", not "processed transcript"
    dtr = dt[gene_status=="KNOWN" & !grepl("pseudo", gene_type) &
             gene_type != "processed_transcript",
             .(chromosome=seqnames, startPoint=start, endPoint=end, strand,
               title = gene_name, gene_name, type, gene_id, gene_type,
               transcript_id, transcript_name)]

    if (!is.null(genes)){
        dtr = dtr[title %in% genes]
    } else if (!is.null(grep) | !is.null(grepe)) {
        if (!is.null(grep)) dtr = dtr[grepl(grep, title)]
        if (!is.null(grepe)) dtr = dtr[!grepl(grepe, title)]
    }

    if (nrow(dtr)==0){
        stop("No more data to present.")
    }

    if (gene.collapse){
        ## collapse by gene
        dtr[, hasCds := is.element("CDS", type), by=gene_id]
        dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","gene")],
                    dtr[hasCds==FALSE][type %in% c("exon", "gene")])
        ## dedup
        dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, gene_id))]
        dtr[, title := gene_name]
        dtr = dtr[type != "transcript"]

        ## group id
        dtr[, gid := as.numeric(as.factor(gene_id))]
        if (verbose) message("Intervals collapsed to gene level.")
    } else {
        ## collapse by transcript
        dtr[, hasCds := is.element("CDS", type), by=transcript_id]
        dtr = rbind(dtr[hasCds==TRUE][type %in% c("CDS","UTR","transcript")],
                    dtr[hasCds==FALSE][type %in% c("exon","transcript")])
        ## dedup
        dtr = dtr[!duplicated(paste(chromosome, startPoint, endPoint, transcript_id))]
        dtr[, title := transcript_name]
        dtr = dtr[type != "gene"]

        ## group id
        dtr[, gid := as.numeric(as.factor(transcript_id))]
        if (verbose) message("Intervals collapsed to transcript level.")
    }

    dtr[, iid := 1:nrow(dtr)]

    ## processing intervals
    intervals.json = dtr[, paste0(
        c(paste0(qw("intervals"),": ["),
          paste(
              "\t{",
              qw("iid"), ":", iid,
              ",", qw("chromosome"), ":", chromosome,
              ",", qw("startPoint"), ":", startPoint,
              ",", qw("endPoint"), ":", endPoint,
              ",", qw("y"), ":", 0,
              ",", qw("title"), ":", qw(title),
              ",", qw("group_id"), ":", qw(gid),
              ",", qw("type"), ":", qw(type),
              ",", qw("strand"), ":", qw(strand),
              "}",
              sep = "",
              collapse = ',\n'),
          "]"),
        collapse = '\n')
        ]

    ## assembling the JSON
    out = paste(c("var dataInput = {",
                  paste(
                      c(meta.json,
                        intervals.json),
                      collapse = ',\n'
                  ),"}"),
                sep = "")

    writeLines(out, filename)
    message(sprintf('Wrote JSON file of %s to %s', infile, filename))
    return(filename)
}

#' getPloidy
#'
#' @export
getPloidy = function(segs){
    if (!is(segs, "GRanges")) stop("Not a GRanges!")

    ## NOTE: doesn't have to be disjoint
    ## if (!isDisjoint(segs)) {
    ##     warning("Must be disjoint!")
    ##     segs = gr.disjoin(segs)
    ## }

    ## MARCIN COMMENT: WHAT IF THERE IS TWO COLUMNS HERE MATCHING CN???
    if (length(cnix <- grep("CN", colnames(mcols(segs)), ignore.case=T))==0)
        message("No copy number (cn) column!")

    ## MARCIN COMMENT: WHAT IF THERE IS TWO COLUMNS HERE MATCHING
    cn = mcols(segs)[, cnix[1]]
    wd = width(segs)
    good.ix = which(!is.na(cn))

    pl = weighted.mean(cn[good.ix], wd[good.ix], na.rm=T)
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


setClass

##############################
## gwalks
##############################
#' gwalks
#'
#' S4 wrapper around GRangesList to store walks info.
#'
#' @name gwalks-class
#' @rdname gwalks-class
#'
#' @import methods
#' @import gUtils
#'
#' @exportClass gTrack
gwalks = setClass("gwalks",
                  contains="GRangesList")
## validity test when intializing
setValidity("gwalks",
            function(object){
                if (!is(object, "GRangesList")){
                    object = tryCatch(GRangesList(object),
                                      error=function(e) return(NULL))
                    if (is.null(object))
                        return("Input can't be converted into a GRangesList.")
                }
                if (isEmpty(object)){
                    message("Empty gwalks.")
                    return(TRUE)
                } else if (!all(strand(unlist(object)) %in% c("+", "-"))) {
                    "All strand info must be present."
                } else if (!all(c("cn", "str", "is.cycle") %in% colnames(values(object)))) {
                    "Required metadata fields: cn, str, and is.cycle."
                } else return(TRUE)
            })
## explicit coercion and that's it!
setAs("GRangesList", "gwalks", function(from){new("gwalks", from)})

## now extend S4 methods special for "gwalks"
## 0) explicit constructor
gwalks = function(...){
    verbose = getOption("gGnome.verbose")

    grl = tryCatch(GRangesList(...),
                   error = function(e) return(NULL))

    if (is.null(grl)) stop("Input must be convertible to GRangesList class.")

    coln = colnames(values(grl))

    ## helping you with the required fields
    if (!"cn" %in% coln) {
        if (verbose) warning("Arbitrarily setting all walk 'cn' to 1.")
        values(grl)$cn = 1
    }
    if (!"str" %in% coln) {
        if (verbose) warning("Arbitrarily setting all walk 'str' to '+'.")
        values(grl)$str = "+"
    }
    if (!"is.cycle" %in% coln) {
        if (verbose) warning("Arbitrarily setting all walk 'is.cycle' to FALSE.")
        values(grl)$is.cycle = FALSE
    }

    return(new("gwalks", grl))
}

## 1) test if the walks are paired up
#' grl.match
#' Matching the GRanges elements in a GRangesList
#'
#' @param gr1, gr2
#' @param ordered if TRUE the order of elements are considered in comparison
#' @param ignore.strand if TRUE the strand info is ignored
#'
#' @return vector same length as grl1, containing the indices of the first matching object in grl2
grl.match = function(grl1, grl2,
                     ordered=FALSE, ignore.strand=FALSE,
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
#' @export
rev.comp = function(gr){
    strmap = setNames(c("+", "-"), c("-", "+"))
    if (!is(gr, "GRanges")){
        stop("Input must be GRanges.")
    } else if (!all(strand(gr) %in% strmap)) {
        stop("Input must be all strand specific.")
    }
    return(rev(gr.strandflip(gr)))
}

## 2) converting a gwalks to

## ============= R6 gWalks class definition ============= ##
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

                     initialize = function(grl=NULL, segs=NULL, paths=NULL,
                                           is.cycle=NULL, cn=NULL, str=NULL,
                                           metacols = NULL, kh=FALSE){
                         if (!is.null(segs)){
                             private$gwFromScratch(segs, paths, is.cycle,
                                                   cn, str, metacols=metacols)
                         } else if (!is.null(grl)) {
                             self$grl2gw(grl, kh=kh)
                         } else {
                             self$nullGWalks()
                         }
                     },

                     ## TODO: construct null gWalks
                     nullGWalks = function(){
                         return(self)
                     },

                     ## TODO: gw2bg convert to a list of bGraphs
                     gw2gg = function(){
                         verbose = getOption("gGnome.verbose")
                         if (!self$isStrandPaired()){
                             ## MARCIN EDIT: NOT SURE WHY THIS FAILS SOMETIMES
                             ## first check the segs
                             try = tryCatch(seg.fill(private$segs),
                                            error = function(e) NULL) ## TODO

                             if (!is.null(try))
                                 private$segs = try
                         }

                         ## then check the paths
                         hB = hydrogenBonds(private$segs)
                         hbmap = hB[, c(setNames(from, to),
                                        setNames(to, from))]

                         ## cache the original and rev comp paths
                         rpaths = lapply(private$paths,
                                         function(x) {
                                             out = rev(hbmap[as.character(x)])
                                             names(out) = NULL
                                             return(out)
                                         })
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

                         ## summarize the edges
                         es = self$path2edges()

                         ## amplitude of each walk
                         amp = rep(private$metacols$cn, elementNROWS(private$paths))
                         cns = table(rep(unlist(private$paths), amp))

                         private$segs$cn = ifelse(
                             as.character(seq_along(private$segs)) %in% names(cns),
                             as.numeric(cns[as.character(seq_along(private$segs))]),
                             0)

                         ## in case two strands are not both present: fill it in
                         pl = getPloidy(private$segs)

                         ## NOTE: rest assured no seg info is lost!
                         gg = gGraph$new(segs = private$segs,
                                         es = es,
                                         ploidy = pl,
                                         purity = 1)
                         ## private$gg = gg
                         return(gg)
                     },

                     gw2grl = function(ix=NULL){
                         if (is.null(ix))
                             ix = seq_along(private$paths)
                         segs = private$segs
                         ss = gr.stripstrand(segs)
                         ss = ss[!duplicated(ss)] %Q% (order(seqnames, start))
                         ss$tile.id = seq_along(ss)

                         mix = match(gr.stripstrand(segs[,c()]), ss[,c()])
                         segs$tile.id = ss$tile.id[mix]

                         ## segs$tile.id = rep(LETTERS[1:23][1:(length(private$segs)/2)], 2)
                         grl = lapply(private$paths[ix],
                                      function(pt){
                                          return(segs[pt])
                                      })
                         grl = GRangesList(grl)
                         mcols(grl) = private$metacols[ix]
                         private$grl = grl
                         return(grl)
                     },

                     grl2gw = function(grl, kh = FALSE, mc.cores=1){
                         if (length(grl)==0) return(self)

                         ## first of all, both strand of one range must be present
                         ## if not double them
                         grs = grl.unlist(grl)
                         segs = grs[!duplicated(grs)]
                         names(segs) = NULL

                         private$segs = segs
                         if ("loose" %in% colnames(values(segs))){
                             segs = segs %Q% (order(loose, strand, seqnames, start))
                         } else {
                             segs = segs %Q% (order(strand, seqnames, start))
                         }

                         paths = lapply(grl, function(x) match(x, segs))
                         mc = as.data.frame(values(grl))
                         if (ncol(mc)==0) {
                             ## if no meta data columns, default
                             mc = data.table(is.cycle = rep(FALSE, length(paths)),
                                             cn = rep(1,length(paths)),
                                             str=rep("+", length(paths)))
                         } else {
                             mc = as.data.table(mc)
                             if (!"cn" %in% colnames(mc)) mc[,cn:=rep(1,length(paths))]
                             if (!"is.cycle" %in% colnames(mc))
                                 mc[,is.cycle := rep(FALSE, length(paths))]
                             if (!"str" %in% colnames(mc)) mc[,str:=rep("+",length(paths))]
                         }
                         private$gwFromScratch(segs = segs, paths = paths, metacols = mc)
                         return(self)
                     },
                     gw2gTrack = function(ix=NULL, colorful=FALSE, mc.cores=1){
                         grl = self$gw2grl(ix)
                         gts = gTrack(grl, draw.path=T)
                         if (colorful) {
                             ## TODO
                         }
                         return(gts)
                     },
                     gw2json = function(filename = ".",
                                        save=TRUE,
                                        trim=TRUE,
                                        mc.cores=1){
                         ## TODO: match up the cids with the gGraph
                         if (save){
                             if (grepl('\\.js(on)*$', filename))
                                 ## if json path was provided
                                 basedir = dirname(filename)
                             else if (filename==".") {
                                 ## default path was provided
                                 basedir = './'
                                 filename = "gwalk.json"
                             } else {
                                 ## a directory was provided
                                 basedir = filename
                                 filename = paste(filename, 'gwalk.json', sep = '/')
                             }

                             if (!file.exists(basedir)) {
                                 message('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }
                         }

                         ## get the gGraph part of JSON

                         gg = self$gw2gg()
                         grl = self$gw2grl()
                         ys = draw.paths.y(grl)

                         ## no walk, just graph
                         if (length(private$paths)==0){
                             gg.js = gg$gGraph2json(filename=filename)
                             return(gg.js)
                         }

                         ## gather data
                         "Create json file for interactive visualization."
                         qw = function(x) paste0('"', x, '"') ## quote

                         ## ADDED BY MARCIN: define regularChr
                         regularChr = c(as.character(1:22), "X", "Y") ## 24 regular chrs

                         ## ALERT: for a clean viz for now, only contain regular chromosomes
                         regsegs.ix = which(as.character(seqnames(private$segs)) %in% regularChr)

                         ## EDIT BY MARCIN: trimming json since CX said we only need "walks:" field and we don't need "intervals:" or "connections:" in walks JSON
                         ## processing nodes
                         ## reduce strand
                         ## remove loose nodes
                         oid = gr2dt(private$segs)[, which(strand=="+" &
                                                           loose==F &
                                                           !is.na(cn) &
                                                           seqnames %in% regularChr)]
                         ## ori ind of rev comps
                         rid = gr2dt(private$segs)[, which(strand=="-" &
                                                           loose==F &
                                                           !is.na(cn) &
                                                           seqnames %in% regularChr)]
                         nodes = private$segs[c(oid, rid)]
                         ## ori ix of loose nodes
                         loose.id = which(private$segs$loose==T)
                         ## binding into dt
                         node.dt = data.table(
                             ## each row is a non-loose positive segment
                             node.id = c(oid,rid),
                             iid = seq_along(nodes),
                             chromosome = as.character(seqnames(nodes)),
                             startPoint = as.character(start(nodes)), ## smaller coor side
                             endPoint = as.character(end(nodes)),
                             strand = as.character(strand(nodes)),
                             ## keep track of gGraph node ids
                             title = paste0(seq_along(nodes), ' (', oid, '|', rid, ')'),
                             type = "interval",
                             y = ""
                         )
                         ## ## converting to JSON 'intervals' string
                         ## intervals.json = node.dt[, paste0(
                         ##     c(paste0(qw("intervals"),": ["),
                         ##       paste(
                         ##           "\t{",
                         ##           qw("iid"), ":", iid,
                         ##           ",", qw("chromosome"), ":", qw(chromosome),
                         ##           ",", qw("startPoint"), ":", startPoint,
                         ##           ",", qw("endPoint"), ":", endPoint,
                         ##           ## ",", qw("y"), ":", y,
                         ##           ",", qw("title"), ":", qw(title),
                         ##           ",", qw("type"), ":", qw(type),
                         ##           ",", qw("strand"), ":", qw(strand),
                         ##           "}",
                         ##           sep = "",
                         ##           collapse = ',\n'),
                         ##       "]"),
                         ##     collapse = '\n')
                         ##     ]

                         ## paths not empty, there must be edges
                         ## TMPFIX: remove NA edges .. not clear where these are coming from
                         ## but likely the result of trimming / hood
                         ## also from NA nodes
                         ed = gg$edges
                         ed = ed[!is.na(from) & !is.na(to) &
                                 from %in% regsegs.ix & to %in% regsegs.ix, ]
                         ed[, eid := paste(from, to, sep="_")]


                         ##ALERT: bc strandlessness, I only retained half of the edges
                         ##for gwalks, we will need strandedness, so will retain everything
                         ed[,":="(soStr = as.character(strand(private$segs[from])),
                                  siStr = as.character(strand(private$segs[to])))]

                         ## ## mapping from type field to label in json
                         eType = setNames(c("REF", "ALT", "LOOSE"),
                                          c("reference", "aberrant", "loose"))
                         ## ## processing edges, cont.
                         ## fmap = node.dt[, .(oid, iid)]; setkey(fmap, oid);
                         ## rmap = node.dt[, .(rid, iid)]; setkey(rmap, rid);

                         ## edge data.table
                         ed.dt = ed[,.(from,
                                       to,
                                       ## source
                                       so = node.dt[node.id==from, iid],
                                       ## sink
                                       si = node.dt[node.id==to, iid],
                                       so.str = ifelse(soStr=="+",1,-1),
                                       si.str = ifelse(siStr=="+",1,-1),
                                       ## diff than defined in es field
                                       ## weight=ifelse(type=="aberrant",
                                       ##               1L, weight),
                                       title = "",
                                       cn,
                                       type = eType[type]),
                                    by=1:nrow(ed)]

                         ed.dt[, cid := 1:.N]
                         ed.dt[,":="(so = so*so.str, si = -si*si.str)]


                         ## ## finally, convert to JSON 'connections' string
                         ## connections.json = ed.dt[, paste0(
                         ##     c(paste0(qw("connections"),": ["),
                         ##       paste(
                         ##           "\t{",
                         ##           qw("cid"), ":", cid,
                         ##           ifelse(is.na(so), "", paste0(",",qw("source"),":")),
                         ##           ifelse(is.na(so), "", so),
                         ##           ifelse(is.na(si), "", paste0(",",qw("sink"),":")),
                         ##           ifelse(is.na(si), "", si),
                         ##           ",", qw("title"), ":", qw(title),
                         ##           ",", qw("type"), ":", qw(type),
                         ##           ",", qw("weight"), ": ", cn,
                         ##           "}",
                         ##           sep = "",
                         ##           collapse = ',\n'),
                         ##       "]"),
                         ##     collapse = '\n')]

                         ## finally, turn node path to edge path
                         ed.dt[, eid := paste(from, to, sep="-")]
                         setkey(ed.dt, "eid")
                         setkey(node.dt, "node.id")

                         path.dt = do.call(
                             `rbind`,
                             mclapply(seq_along(private$paths),
                                      function(pti){
                                          this.pname = names(private$paths)[pti]
                                          this.npath = private$paths[[pti]]
                                          this.cyc = private$metacols[pti, is.cycle]


                                          ## MARCIN EDIT: FIX TO DEAL WITH LENGTH 1 CYCLES
                                          this.epath.eid = c()
                                          if (length(this.npath)>1)
                                          {
                                              this.epath.eid =
                                                  paste(this.npath[1:(length(this.npath)-1)],
                                                        this.npath[2:length(this.npath)],
                                                        sep="-")
                                          }
                                          this.ys = ys[[this.pname]]
                                          if (this.cyc){
                                              ##if (this.cyc & length(this.npath)>1){
                                              this.epath.eid = c(this.epath.eid,
                                                                 paste(this.npath[length(this.npath)],
                                                                       this.npath[1], sep="-"))
                                          }

                                          this.ndt = data.table(nid = this.npath, y = this.ys)
                                          this.ndt = cbind(this.ndt,
                                                           node.dt[.(this.ndt$nid),
                                                                   .(node.id, iid,
                                                                     chromosome, startPoint, endPoint,
                                                                     strand, title, type)])

                                          ## MARCIN EDIT: RETURN NULL IF NO NODES AFTER regularChr
                                          ## TRIMMING
                                          this.ndt = this.ndt[!is.na(iid), ]
                                          if (nrow(this.ndt)==0)
                                              return(NULL)

                                          this.nids.json =
                                              this.ndt[!is.na(iid),
                                                       paste0(
                                                           c(paste(
                                                               "\t\t{",
                                                               qw("iid"), ":", iid,
                                                               ",", qw("chromosome"), ":", qw(chromosome),
                                                               ",", qw("startPoint"), ":", startPoint,
                                                               ",", qw("endPoint"), ":", endPoint,
                                                               ",", qw("y"), ":", sprintf("%.2f", y),
                                                               ",", qw("title"), ":", qw(title),
                                                               ",", qw("type"), ":", qw(type),
                                                               ",", qw("strand"), ":", qw(strand),
                                                               "}",sep = "",
                                                               collapse = ',\n'
                                                           )),
                                                           collapse="\n"
                                                       )]

                                          ## ALERT: throwing away good edges
                                          ## just bc they are not in ed.dt

                                          ## EDIT BY MARCIN:
                                          ## SOME VALID PATHS WILL HAVE >=1 nodes and NO EDGES
                                          ## if (any(!this.epath.eid %in% ed.dt[, eid])) return(NULL)

                                          ## just remove any edges that are off the grid, should be good enough
                                          this.epath.eid = this.epath.eid[this.epath.eid %in% ed.dt[, eid]]

                                          this.cids.json = ""

                                          if (length(this.epath.eid)>0)
                                          {
                                              this.pdt = data.table(eid = this.epath.eid)
                                              this.pdt = merge(this.pdt, ed.dt, by="eid")

                                              ## if (nrow(this.pdt)==0) return(NULL)

                                              this.cids.json = this.pdt[this.epath.eid,
                                                                        paste0(
                                                                            c(paste(
                                                                                "\t\t{",
                                                                                qw("cid"), ":", cid,
                                                                                ifelse(is.na(so), "", paste0(",",qw("source"),":")),
                                                                                ifelse(is.na(so), "", so),
                                                                                ifelse(is.na(si), "", paste0(",",qw("sink"),":")),
                                                                                ifelse(is.na(si), "", si),
                                                                                ",", qw("title"), ":", qw(title),
                                                                                ",", qw("type"), ":", qw(type),
                                                                                ",", qw("weight"), ": ", cn,
                                                                                "}",
                                                                                collapse = ',\n'
                                                                            )),
                                                                            collapse="\n"
                                                                        )]
                                          }

                                          this.mc = private$metacols[pti,][, pid := as.numeric(this.pname)]

                                          this.mc[, cids.js := this.cids.json]
                                          this.mc[, nid.js := this.nids.json]

                                          return(this.mc)
                                      },

                                      mc.cores=mc.cores))

                         ##EDIT BY MARCIN
                         ##                         path.dt[, pid := 1:.N]
                         if (any(is.na(path.dt$pid)))
                         {
                             warning('Making numeric pathnames')
                             path.dt[, pid := 1:.N]
                         }

                         path.json = path.dt[str=="+",
                                             paste0(
                                                 c(paste0(qw("walks"),": ["),
                                                   paste0("\t{",
                                                          qw("pid"), ":", pid, ",",
                                                          qw("cn"), ":", cn, ",",
                                                          qw("type"), ":", qw(ifelse(is.cycle, "cycle", "path")), ",",
                                                          qw("strand"), ":", qw(str), ",",
                                                          qw("cids"), ":", "[\n", cids.js,"\n\t]",",\n",
                                                          qw("iids"), ":", "[\n", nid.js,"\n\t]",
                                                          "}",
                                                          collapse=",\n"),
                                                   "]"
                                                   ),
                                                 collapse="\n"
                                             )]


### ADDED BY MARCIN: "gwalks: " header to json file
                         out.json = paste0('{"gwalks":\n{',
                                           paste(
                                               c(
### EDIT BY MARCIN: see above, we don't need these sections anymore according to CX
                                        #                                               intervals.json,
                                        #                                               connections.json,
                                                   path.json
                                               ),
                                               collapse = ',\n'
                                           ),"}\n}")

                         if (save) {
                             message("Writing JSON to ", paste(normalizePath(basedir),filename, sep="/"))
                             writeLines(out.json, filename)
                             return(filename)
                         }
                         return(out.json)
                     },
                     v2e = function(mc.cores=1){
                         ## converting default node path into edges paths
                         if (is.null(private$gg)){
                             gg = self$gw2gg()
                         } else {
                             gg = private$gg
                         }
                     },
                     ## TODO: helper function to turn paths into edges
                     path2edges = function(mc.cores=1){
                         ## whenever this function runs, it will assign result to
                         ## private$es, which will be refreshed to NULL whenever
                         ## a modifying action happens
                         es = do.call(
                             'rbind',
                             mclapply(1:length(private$paths),
                                      function(i){
                                          if (private$metacols[i, cn==0])
                                              return(NULL)
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
                                          ## thisWeight = width(private$segs)[thisFrom]
                                          thisEs = data.table(from = thisFrom,
                                                              to = thisTo,
                                                              cn = private$metacols$cn[i],
                                                              type = "unknown",
                                                              ## weight = thisWeight,
                                                              path.ix = i)
                                          if (thisEs[, any(is.na(to) | is.na(from))])
                                              browser()
                                          return(thisEs)
                                      },
                                      mc.cores=mc.cores)
                         )

                         ## if same edges shows up more than once, dedup and populate cn
                         es[, tmp.id := paste(from, to, sep="-")]
                         setkey(es, "tmp.id")
                         old.es = es
                         es[, cn := sum(cn), by=tmp.id]
                         es = es[!duplicated(tmp.id), .(from, to, cn, type)]

                         ## test junction balance
                         by.from = es[, .(out.cn = sum(cn)), by=from]
                         setkey(by.from, "from")
                         by.to = es[, .(in.cn = sum(cn)), by=to]
                         setkey(by.to, "to")

                         es = etype(private$segs, es)
                         private$es = es
                         return(es)
                     },

                     simplify = function(mod=TRUE, reorder=FALSE){
                         verbose = getOption("gGnome.verbose")
                         ## TODO: merge ref connected segments into one big
                         ## MOMENT
                         if (is.null(private$es))
                             es = self$path2edges()
                         else
                             es = private$es

                         es = etype(private$segs, es, force=T)
                         new.segs = private$segs[, c()]
                         setkey(es, eid)

                         new.paths = list()
                         for (e in seq_along(private$paths)){
                             pth = private$paths[[e]]
                             ## find out runs of at least one reference edge
                             ep = data.table(from = shift(pth),
                                              to = pth)[-1, ]
                             ep[, eid := paste(from, to)]
                             ep[, type := es[.(ep$eid), type]]
                             if (!any(ep$type=="reference")){
                                 break
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

                             to.merge = ep[run.ix>0, .(n.ix = c(head(from, 1), to)), by=run.ix]
                             ## create new node
                             new.node =
                                 do.call(`c`,
                                         sapply(unique(to.merge$run.ix),
                                                function(ix){
                                                    ns = reduce(private$segs[to.merge[run.ix==ix, n.ix]])
                                                    ns$run.ix = ix
                                                    return(ns)
                                                })
                                         )

                             if (verbose)
                                 warning("Updating segs, metadata fields are modified. Proceed with care.")

                             ## ALERT: modifying segs
                             new.node$new.ix = new.node$run.ix + length(new.segs)
                             new.segs = c(new.segs, new.node[,c()])
                             new.nix = new.node$new.ix[new.node$run.ix]

                             ## modifying paths
                             ep[, ":="(last.run.ix = shift(run.ix),
                                       next.run.ix = c(tail(run.ix, -1), NA))]

                             ep[run.ix==0 & last.run.ix!=0, from := new.nix[last.run.ix]]
                             ep[run.ix==0 & next.run.ix!=0, to := new.nix[next.run.ix]]
                             new.ep = ep[run.ix==0]
                             if (nrow(new.ep)==0){
                                 new.paths[[e]] = new.nix
                             } else {
                                 new.paths[[e]] = new.ep[, c(head(from, 1), to)]
                             }
                         }

                         ## final step, dedup the segments, relabel the paths
                         if (any(duplicated(new.segs))){
                             nr.segs = new.segs
                             nr.map = match(new.segs, nr.segs)
                             new.paths = lapply(new.paths, function(x) nr.map[x])
                             new.segs = nr.segs
                         }

                         if (mod==TRUE){
                             ## modify required fields THIS instance
                             private$segs = new.segs
                             private$paths = new.paths
                             ## reset all optional fields
                             private$es = NULL
                             private$grl = NULL
                             private$gg = NULL
                             private$rp.map = NULL
                             if (reorder==TRUE){
                                 self$label()
                             }
                             return(self)
                         } else {
                             out = gWalks$new(segs = new.segs, paths = new.paths)
                             if (reorder==TRUE){
                                 out$label()
                             }
                             return(out)
                         }
                     },

                     subset = function(ix){
                         ## TODO:subset and overwrite `[`
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
                         if (!all(table(gr.match(private$segs, private$segs))==2))
                             return(FALSE)
                         else if (any(duplicated(private$segs)))
                             return(FALSE)

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
                             private$grl = NULL
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
                     grl = NULL,
                     gg = NULL,
                     rp.map = NULL,

                     ## ----- private methods
                     ## you can give me columns separately as vectors
                     ## or you can give me data.frame as a whole
                     gwFromScratch = function(segs, paths=NULL, is.cycle=NULL,
                                              cn=NULL, str=NULL, metacols=NULL){
                         ## segs must be a GRanges
                         if (!is(segs, "GRanges")) stop("segs needs to be a GRanges.")

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
                     segstats = function(){
                         return(private$segs)
                     },
                     td = function(){
                         ## default viz, will not show CN==0 or str=="-"
                         ix = private$metacols[, which(cn>0 & str=="+")]
                         return(self$gw2gTrack(ix))
                     },
                     path = function(){
                         return(private$paths)
                     },
                     values = function(){
                         return(self$metaCols())
                     },
                     json = function(fn = "."){
                         return(self$gw2json(fn))
                     }
                 ))

## ============= R6 gWalks exported functions ============= ##
setAs("gWalks", "gwalks",
      function(from){
          return(from$gw2grl())
      })

setAs("gWalks", "GRangesList",
      function(from){
          return(from$gw2grl())
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

## ============= Utility functions ============= ##
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
etype = function(segs, es, force=FALSE, both=FALSE){
    if (!is(segs, "GRanges")) stop("segs must be GRanges")
    if (!is(es, "data.frame")) stop("es must be data.frame")
    if (!all(c("from", "to") %in% colnames(es))) stop("'from' & 'to' must be in es!")
    if (!is(es, "data.table")) es = as.data.table(es)
    if (nrow(es)==0 | length(segs)==0) return(NULL)

    if ("type" %in% colnames(es) & force==FALSE)
        return(es)

    ## as definition of graph, segs and es must be both sets
    ## so no redundancy
    es2 = copy(es)
    es2[, type := "unknown"]

    ## first dedup nodes
    if (any(dup.ix <- duplicated(segs))){
        new.segs = segs[!dup.ix]
        smap = match(segs, new.segs)
        es2[, ":="(from = smap[from],
                   to = smap[to])]
    }

    ## then dedup edges
    es2[, ":="(eid = paste(from, to))]
    if (es2[, anyDuplicated(eid)]){
        if ("cn" %in% colnames(es2)){
            ## if edges comes with CN field, add them together
            es2[, .(from, to, cn=sum(cn), type), by=eid]
        } else {
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

        which.loose.src = which(width(segs)==1 & segs$degree==1 & segs$out.degree==1)
        loose.src.partner = gr.start(es2[from %in% which.loose.src, segs[to, c()]], ignore.strand=FALSE)
        which.loose.src = which.loose.src[which(loose.src.partner==segs[which.loose.src, c()])]

        which.loose.sink = which(width(segs)==1 & segs$degree==1 & segs$in.degree==1)
        loose.sink.partner = gr.end(es2[to %in% which.loose.sink, segs[to, c()]], ignore.strand=FALSE)
        which.loose.sink = which.loose.sink[which(loose.sink.partner==segs[which.loose.sink, c()])]

        which.loose = c(which.loose.src, which.loose.sink)
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

    if (both) return(list(segs = segs, es = es2))

    return(es2)
}

#' @name setxor
#' @title XOR operation on sets
#'
#' @param A set A
#' @param B set B
#'
#' @return The set of elements belong to either A or B, but not both.
#' @author Marcin Imielinski
#' @export
setxor = function (A, B)
{
    return(setdiff(union(A, B), intersect(A, B)))
}

#' @name write.tab
#' @title wrapper around write.table
#' @author Marcin Imielinski
#' @export
write.tab = function (x, ..., sep = "\t", quote = F, row.names = F)
{
    if (!is.data.frame(x))
        x = as.data.frame(x)
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
#' @export
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

###############################
#' @name isInteger
#' @title Testing if a numeric value is integer
#'
#' @param x numeric vector
#' @param eps infinitely small positive number
#'
#' @return logical vector of same length
#' @author Xiaotong Yao
isInterger = function(x, eps = 1e-300){
    if (!is.numeric(x))
        return(FALSE)

    return(x %% 1 < eps)
}

################################


gencode2json = function(gencode=NULL, file="."){
    ## ASSUMPTION: gencode is a GR, presumably read from skidb function
    if (is.null(gencode)){
        require(skidb)
        ## ALERT: if you don't give me anything, I'm only including known genes
        gencode = skidb::read_gencode()
    }
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

#' read_vcf: utility function to read VCF into GRanges object
#'
#' @name read_vcf
#' @import VariantAnnotation
#' @export
read_vcf = function (fn, gr = NULL, hg = "hg19", geno = NULL, swap.header = NULL,
                     verbose = FALSE, add.path = FALSE, tmp.dir = "~/temp/.tmpvcf",
                     ...)
{
    require(VariantAnnotation)
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

#' ra_breaks: parse junction data from various common formats
#'
#' @name ra_breaks
#' @description
#'
#' @return a junctions object
#'
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import data.table
#'
#' @export
ra_breaks = function(rafile, keep.features = T, seqlengths = hg_seqlengths(), chr.convert = T, geno=NULL,
                     snowman = FALSE, swap.header = NULL,  breakpointer = FALSE, seqlevels = NULL, force.bnd = FALSE, skip = NA,
                     get.loose = FALSE, pad = NULL){
    ## if TRUE will return a list with fields $junctions and $loose.ends
    if (is.character(rafile))
    {
        if (grepl('.rds$', rafile)){
            ra = readRDS(rafile)

            ## a few check points
            if (!is(ra, "GRangesList")) stop("Junctions must be GRangesList!")

            if (any(elementNROWS(ra)!=2)) stop("Each element must be length 2!")

            bps = unlist(ra)
            if (any(!(strand(bps) %in% c("+", "-")))) stop("Breakpoints must have orientation!")

            if (any(width(bps)>1)) stop("Breakpoints must be points!")

            return(ra)
        } else if (grepl('(.bedpe$)', rafile)){
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

            ra = ra.dedup(ra)

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
                             ra1 = gr.strandflip(GRanges(rafile$chr1, IRanges(rafile$start1, rafile$end1), strand = rafile$str1))
                             ra2 = gr.strandflip(GRanges(rafile$chr2, IRanges(rafile$start2, rafile$end2), strand = rafile$str2))
                             out = grl.pivot(GRangesList(ra1, ra2))
                         }



    if (keep.features){
        values(out) = rafile[, ]
    }

    if (!is.null(pad)){
        out = ra.dedup(out, pad = pad)
    }

    if (!get.loose){
        return(out)
    }
    else{
        return(list(junctions = out, loose.ends = GRanges()))
    }

    return(new("junctions", out))
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
#' @param GRange input jab object
#' @param file output json file
#' @author Marcin Imielinski
#' @export
gr2json = function(intervals, file, y = rep("null", length(intervals)), labels = '', maxcn = 100, maxweight = 100)
{

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

    writeLines(out, file)
    return(out)
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
            ra1 = gr.strandflip(gr.end(jab$segstats[jab$ab.edges[ix,1,1]], 1, ignore.strand = F))
            ra2 = gr.start(jab$segstats[jab$ab.edges[ix,2,1]], 1, ignore.strand = F)
            ra1 = GenomicRanges::shift(ra1, ifelse(as.logical(strand(ra1)=='+'), -1, 0))
            ra2 = GenomicRanges::shift(ra2, ifelse(as.logical(strand(ra2)=='+'), -1, 0))
            ra = grl.pivot(GRangesList(ra1,ra2))
        }
    }

    if (length(ra)==0){
        return(list())
    }

    if (length(query)==0 | length(subject)==0){
        return(list())
    }

    if (is.null(names(query))){
        names(query) = 1:length(query)
    }

    if (is.null(names(subject))){
        names(subject) = 1:length(subject)
    }

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

    if (length(query)==0 | length(subject)==0){
        return(list())
    }

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
        if (verbose){
            cat('starting interval', i, 'of', length(ix.query), '\n')
        }

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

        if (verbose){
            cat('finishing interval', i, 'of', length(ix.query), ':', paste(round(D.rel[i, ix],2), collapse = ', '), '\n')
        }

        return(list(D.rel = D.rel, D.ref = D.ref, D.ra = D.ra, D.which = D.which))
    }, mc.cores = mc.cores)

    for (i in 1:length(tmp))
    {
        if (class(tmp[[i]]) != 'list'){
            warning(sprintf('Query %s failed', ix.query[i]))
        }
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
            stop('Error: bp1 and bp2 must be signed intervals (i.e. either + or -)')

        if (length(bp1) != length(bp2)){
            stop('Error: bp1 and bp2 inputs must have identical lengths')

                                        #    if (sum(width(reduce(bp1))) != sum(width(bp1)) | sum(width(reduce(bp2))) != sum(width(bp2)))
                                        #      stop('bp1 or bp2 cannot have duplicates / overlaps (with respect to location AND strand)')
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
    else{
        adj.ref = Matrix(FALSE, nrow = 2*length(tile), ncol = 2*length(tile),
                         dimnames = rep(list(as.character(c(1:length(tile), -(1:length(tile))))), 2))
    }

    ## current tile is partition of genome only in positive orientation + dummy intervals for breakpoints
    ## output tile is forward partition and followed by reverse partition
    ## (this is what is currently referenced by adj.ref and adj.ab)
    ## TODO: clean up this part
    tmp.nm = as.character(c(1:length(tile), -(1:length(tile))))
    tile = c(tile, gr.strandflip(tile))
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

    if (length(prior)!=ncol(K)){
        stop('Error: prior must be of the same length as number of columns in K')

                                        # variable indices
    }
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
    else{
        Ac = Zero[1,,drop = FALSE]
        ## combine constraints
    }
    A = rBind(cBind(K, Zero[rep(1, nrow(K)), M.ix]), Amub, Amlb, Ac);
    b = c(e, rep(0, nrow(Amlb)*2), rep(0, nrow(Ac)));
    sense = c(rep('E', nrow(K)), rep('L', nrow(Amlb)), rep('G', nrow(Amlb)), rep('E', nrow(Ac)))
    vtype = c(rep('I', length(v.ix)), rep('B', length(M.ix)))

    cvec = c(rep(0, length(v.ix)), prior-cpenalty*rep(1, length(M.ix)))


    sol = Rcplex(cvec = cvec, Amat = A, bvec = b, sense = sense, Qmat = NULL, lb = 0, ub = Inf, n = nsolutions, objsense = objsense, vtype = vtype, control = c(list(...), list(tilim = tilim, epgap = epgap)))

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





## TODO:
## 1) make it always going upwards
## 2) only change Y when travel through a aberrant junction OR when strand changes (try this 1st)
## 3) OR not change when it goes far away enough

#' @name draw.paths.y
#' Determine the Y axis elevation of segments in a walk
#'
draw.paths.y = function(grl, path.stack.x.gap=0, path.stack.y.gap=1){
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

    gr$last = iix = NULL ## NOTE fix
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
    windows = as(coverage(gr), 'GRanges');
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
#'
gr.flatmap = function(gr, windows, gap = 0, strand.agnostic = TRUE, squeeze = FALSE, xlim = c(0, 1))
{

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
#'
affine.map = function(x, ylim = c(0,1), xlim = c(min(x), max(x)), cap = F, cap.min = cap, cap.max = cap, clip = T, clip.min = clip, clip.max = clip)
{
                                        #  xlim[2] = max(xlim);
                                        #  ylim[2] = max(ylim);

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



#' @name seg.fill
#' Supplement the other strand if missing from input.
#'
#' @export
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




#' @name hydrogenBonds
#' Return a edge data.table connecting two input segments that are two strands of the same range
#'
#' @param segs GRanges
#'
#' @export
hydrogenBonds = function(segs){
    ## collapse +/- strand
    ss = unique(gr.stripstrand(segs))
    idss = match(gr.stripstrand(segs), ss)
    if (!all(table(idss)==2)){
        stop("Error: Malformed object. Suggest creation again.")
    }
    tmpDt = data.table(ssid = seq_along(ss))
    tmpDt[, ":="(n1 = which(idss==ssid)[1],
                 n2 = which(idss==ssid)[2]), by=ssid]
    hydrogenBs = tmpDt[, .(from = n1, to = n2,
                           type="hydrogen")]
    return(hydrogenBs)
}
