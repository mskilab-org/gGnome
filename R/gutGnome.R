## FIXME: REMOVE THESE LIBRARIES
library(gUtils)
library(data.table)
library(R6)
library(gTrack)


#'' @title gGnome
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
                                        #-'    Github: https://github.com/mskilab/gGnome
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
#'              ploidy=NULL, purity=NULL, segs = NULL)
#'
#'   gread(filename)
#'
#' Public fields:
#'
#'   # A GRanges of the nodes
#'   gg$nodes
#'
#'   # A data.table of the edge connectsion
#'   gg$edges
#'
#'   gg$junctions
#'
#'   # A directed igraph of the edge connections
#'   gg$graph
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
#'   gg$addJuncs(juncs)
#'
#'   gg$addSegs(bps)
#'
#'   gg$addVariant(tile, cn = TRUE)
#'
#'   gg$addSNVs(tile, cn = TRUE)
#'
#'   gg$resetSubgraphs(subs = NULL)
#'
#'   gg$addSubgraph(gr, es = NULL, parent = NULL)
#'
#'   gg$deleteSubgraph(index, parent = NULL)
#'
#'   gg$karyograph(tile = NULL, juncs = NULL)
#'
#'   gg$simplify(mod = TRUE)
#'
#'   gg$decouple(mod = TRUE)
#'
#'   gg$add(gg, mod = FALSE, decouple = TRUE)
#'   ## TODO: gg$add(); gg$subtract()
#'
#'   gg$jab2gg(jabba, regular = NULL)
#'
#'   gg$wv2gg(weaver)
#'
#'   gg$pr2gg(prego)
#'
#'   gg$cougar2gg(cougar)
#'
#'   gg$remixt2gg(remixt)
#'
#'   gg$print()
#'
#'   print(gg)
#'
#'   gg$plot(pad = 1000) ## TODO: add node and edge viz configurations
#'
#'   plot(gg)
#'
#'   gg$window(pad = 0)
#'
#'   gg$layout()
#'   ## TODO: rewrite layout(gg); add layout method option
#'
#'   gg$summary()
#'   ## TODO: rewrite summary(gg)
#'
#'   gg$gg2td(seq.col,...)
#'
#'   gg$json(filename = ".", maxcn = 100, maxweight = 100)
#'
#'   gg$html(filename = ".", gGnome.js = Sys.getenv("DEFAULT_GGENOMEJS"), maxcn = 100, maxweight=100)
#'
#'   gg$gg2js(filename = ".", maxcn = 100, maxweight = 100, save = TRUE, setting = NULL, no.y = FALSE)
#'
#'   gg$julie.copy.gg2js(filename = ".", maxcn = 100, maxweight = 100, save = TRUE, setting = NULL, no.y = FALSE)
#'
#'   gg$length()
#'
#'   length(gg)
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
#'   gg$make.balance(mod = TRUE)
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
#' @import prodlim
#' #' @importFrom stringr str_trim
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


## THIS IS MY EDIT


## ================= Node class definition ============= ##
#' @export
Node = setClass("Node")

Node = R6::R6Class("Node",
                   public = list(

                       ## Set up the constructor
                       ## index - snode.id
                       ## graph
                       initialize = function(snode.id, graph)
                       {
                           if (!snode.id %in% graph$fullnodes$snode.id) {
                               stop("Not valid snode.id")
                           }

                           private$graph = graph
                           private$node.id = abs(snode.id)
                           private$snode.id = snode.id
                           
                           lookup = private$graph$queryLookup(snode.id)
                           private$index = lookup[, index]
                           private$rev.index = lookup[, rev.index]

                           private$orientation = ifelse(snode.id > 0, "+", "-")
                       },

                       flip = function() {
                           tmp = private$index
                           private$index = private$rev.index
                           private$rev.index = tmp
                           private$orientation = ifelse(private$orientation == "+", "-", "+")

                           return(self)
                       },

                       subset = function(i)
                       {
                           private$node.id = private$node.id[i]
                           private$index = private$index.id[i]
                           private$rev.index = private$rev.index[i]
                           private$orientation = private$orientation[i]
                       },

                       print = function()
                       {
                           message('Node object of length ', self$length())
                       },
                       
                       length = function()
                       {
                           return(length(private$node.id))
                       }
                   ),
                   ## NODES PRIVATE FIELDS
                   private = list(
                       ## stores the ID of this node in the gGraph
                       node.id = NULL,
                       snode.id = NULL,

                       index = NULL,
                       rev.index = NULL,

                       orientation = NULL,
                       
                       ## Pointer to gGraph
                       graph = NULL
                   ),


                   ## NODES ACTIVE BINDINGS
                   active = list(
                       id  = function() {
                           return(private$node.id)
                       },

                       ## FIMXE: This should be more of a list for each nothing rather than a vector for all nodes
                       left = function() {
                           leftNodes = private$graph$pedges[to %in% private$index, from]
                           return(Node$new(snode.id = private$graph$private$pnodes$snode.id[leftNodes],
                                           private$graph))
                       },

                       ## Returns the nodes connected to the right of the nodes
                       right = function() {
                           rightNodes = private$graph$pedges[from %in% private$index, to]
                           return(Node$new(snode.id = private$graph$private$pnodes$snode.id[rightNodes],
                                           private$graph))
                       },

                       ## Returns a gGraph containing the nodes in this Node Class
                       subgraph = function() {
                           
                       
                           
                       }
                       

                       
                   )
                   )


## ================= Edge class definition ============= ##
#' @export
Edge = setClass("Edge")

Edge = R6::R6Class("Edge",
                   public = list(
                       ## 
                       initialize = function(sedge.id, graph) {
                           if (!sedge.id %in% graph$fulledges[, sedge.id]) {
                               stop("sedge.id is not in the valid list of sedge.id's in the graph")
                           }

                           private$orientation = ifelse(sedge.id > 0, "+", "-")
                           private$graph = graph
                           private$sedge.id = sedge.id
                           private$edge.id = abs(sedge.id)
                           
                           private$edges = private$graph$fulledges[private$edge.id %in% edge.id]
                       }
                   ),

                   
                   private = list(
                       ## data.table of the edges in this class
                       edges = NULL,

                       edge.id = NULL,
                       sedge.id = NULL,
                       
                       ## Defaults to positive, switch to negative if you flip it
                       orientation = NULL,

                       ## Pointer to the gGraph
                       graph = NULL
                   ),


                   active = list(
                       

                       
                   )
                   )


#' @name [.Node
#' @title Node
#' @description
#'
#' Overloads subset operator for nodes
#'
#' @param obj Node object This is the FishHookObject to be subset
#' @param i subset of
#' @return A new FishHook object that contains only the given hypotheses and/or covariates
#' @export
'[.Node' = function(obj, i = NULL){
    ret = obj$clone()
    ret$subset(i)
    return(ret)
}

gGraph = R6::R6Class("gGraph",
                     public = list(

                         ## public fields
                         ## name = NULL,
                         ## refG = "GENOME", ## seqinfo of ref genome
                         ## constructor INIT
                         initialize = function(genome = NULL,
                                               nodes = NULL,
                                               edges = NULL,
                                               looseterm = TRUE)
                         {
                             ## control how to construct
                             verbose = getOption("gGnome.verbose")
                             if (!is.null(nodes)){
                                 private$gGraphFromNodes(nodes, edges, looseterm = looseterm)
                             } else {
                                 private$emptyGGraph(genome)
                             }

                             ## FIXME: MOVE THIS GUY INTO CONSTRUCTORS, TMP FIX
                             private$buildLookupTable()
                         },

                         ## snode.id is a vector of signed node.id's from our gGraph
                         queryLookup = function(snode.id) {
                             return(lookup[.(snode.id)])
                         },

                         length = function()
                         {
                             return(self$nodes$length())
                         },

                         print = function() {
                             cat('\n\n')
                             cat(sprintf('gGraph with %s nodes and %s edges', self$length(), nrow(private$pedges)))
                         }
                     ),


                     private = list( #### PRIVATE GGRAPH
                         ## ----- private fields
                         ## ===== required
                         ## node/vertex, a GRanges obj of strand-specific ranges
                         pnodes = NULL,
                         ## data.table of all edges in g, from, to, cn, type
                         ## type can be ref, aberrant, loose
                         pedges = NULL,

                         ## Lookup table - must be reset with buildLookupTable
                         lookup = NULL,


                         ## ===== optional slots
                         ## ALERT: whenever segs or es changes, all of these need to be reset!
                         ## igraph obj representing the graph structure
                         pgraph = NULL,

                         ## ----- private methods
                         ## reset optional fields
                         reset = function(){
                             private$pgraph = NULL
                         },


                         ## Builds the lookup table for this gGraph which is stored in private$lookup
                         ## Lookup table contains the columns:
                         ##     snode.id - signed node id of the nodes
                         ##     index - index corresponding to the input snode.id
                         ##     rev.index - index corresponding to the complement snode.id
                         buildLookupTable = function() {
                             lt = match(private$pnodes$snode.id, -private$pnodes$snode.id)
                             private$lookup = data.table(snode.id = private$pnodes$snode.id,
                                                         index = seq_along(private$pnodes),
                                                         rev.index = lt)
                             setkey(private$lookup, snode.id)
                         },


                         ## @name emptyGGraph
                         ## @brief Constructor, initializes an empty gGraph object. If the user does not provide genome
                         ##        and using to update this class, this will try to inherit the current gGraph's seqinfo.
                         ##        See class documentation for usage.
                         ## @param genome Optional seqinfo if the user wants to
                         emptyGGraph = function(genome=NULL)
                         {
                             "Create empty gGraph."

                             ## If there is an old private segs and it has seqinfo, inherit that seqinfo
                             ## If the user provides genome, skip this step
                             if (is.null(genome)) {

                                 if(!is.null(private$pnodes) && length(seqinfo(private$pnodes)@seqlengths) > 0) {

                                     if (getOption("gGnome.verbose")){
                                         message("Adopting reference instance's seqinfo. Ignoring input.")
                                     }
                                     genome = seqinfo(private$pnodes)
                                 }
                             }

                             ## Set up the private fields to be empty
                             private$pnodes = GRanges()
                             private$pgraph = igraph::make_empty_graph()
                             private$pjunctions = new("junctions", GRangesList())
                             private$pedges = data.table(from=integer(0),
                                                         to=integer(0),
                                                         type=character(0))

                             ## Try to set the seqinfo to what the user provides
                             self$set.seqinfo(genome = genome)

                             return(self)
                         },

                         ##NODES XXXXX
                         ## granges of intervals, <sign is IGNORED>
                         ## edges is a data.table with required fields n1, n1.side, n2, n2.side representing node 1 and 2 index and side, where side = 0 is left side, side = 1 is right
                         ## to specify loose ends provide edges with NA in <either> n1 or n2 field
                         gGraphFromNodes = function(nodes,
                                                    edges = NULL,
                                                    looseterm = TRUE,
                                                    keepLoose = FALSE)
                         {

                             loose.left = loose.right = c()

                             if (is.null(edges) || nrow(edges) == 0)
                             {
                                 private$pedges = data.table(from = integer(0),
                                                             to = integer(0),
                                                             type = character(0))
                                 edges = data.table()

                                 if (length(nodes)==0){
                                     private$pnodes = GRanges(seqinfo = seqinfo(nodes))
                                     return(self)
                                 }

                                 if (looseterm)
                                 {
                                     loose.right = loose.left = 1:length(nodes)
                                 }
                             }
                             else
                             {
                                 if (!all(c('n1', 'n1.side', 'n2', 'n2.side') %in% names(edges))) {
                                     stop('edges table not in proper format: requires columns n1, n2, n1.side, n2.side, where n1 and n2 index the provided nodes and the side arguments specify right side if 1 and left side if 0')
                                 }
                                 edges = as.data.table(edges)

                                 if (any(edges$n1>length(nodes) | edges$n2>length(nodes), na.rm = TRUE))
                                 {
                                     stop('n1 or n2 fields in edges table indexing out of range nodes')
                                 }

                                 if (any(ix  <- is.na(edges$n1) & is.na(edges$n2)))
                                 {
                                     warning(paste('Removed', sum(ix), 'edges from graph that have NA in both n1 and n2, edges should have either n1 or n2 non AN'))
                                     edges = edges[!ix, ]
                                 }

                                 if (any(ix <- (!is.na(edges$n1) & is.na(edges$n1.side)) | (!is.na(edges$n2) & is.na(edges$n2.side))))
                                 {
                                     stop(paste('All non NA n1 or n2 must have an n1.side or n2.side, but we found', sum(ix), 'edges that violate this requirement'))
                                 }

                                 if (looseterm)
                                 {
                                     loose.right = setdiff(1:length(nodes), c(edges[n1.side==1, n1],  edges[n2.side==1, n2]))
                                     loose.left = setdiff(1:length(nodes), c(edges[n1.side==0, n1],  edges[n2.side==0, n2]))
                                 }
                             }

                             nodes$loose = FALSE

                             if (length(loose.right)>0)
                             {
                                 edges = rbind(edges, data.table(n1 = loose.right, n1.side = 1, n2 = NA, n2.side = NA), fill = TRUE)
                             }

                             if (length(loose.left)>0)
                             {
                                 edges = rbind(edges, data.table(n2 = loose.left, n2.side = 0, n1 = NA, n1.side = NA), fill = TRUE)
                             }

                             ## add all loose end nodes
                             if (nrow(edges)>0)
                             {
                                 edges[, type := "unknown"]
                                 if (any(nix <- is.na(edges$n1)))
                                 {
                                     ## adding loose ends nodes for n1 based on n2 side which is  either left (n2.side == 0) or right (n2.side == 1)
                                     ## if left, then we use gr.start of n2 node, if right we use gr.end of n2 node to make the loose end
                                     loose.nodes = parse.gr(edges[nix, ifelse(n2.side == 1, gr.string(gr.end(nodes[n2])), gr.string(gr.start(nodes[n2])))])
                                     ## now we fill in the n1 field with the index of the newly appended nodes
                                     edges[nix, n1 := length(nodes) + 1:.N]
                                     edges[nix, n1.side := ifelse(n2.side == 1, 0, 1)]
                                     edges[nix, type := 'loose']
                                     nodes = grbind(nodes, loose.nodes)
                                 }

                                 ## now do the same thing for n2 NA nodes, using features of n1 node
                                 if (any(nix <- is.na(edges$n2)))
                                 {
                                     ## adding loose ends nodes for n1 based on n2 side which is  either left (n2.side == 0) or right (n2.side == 1)
                                     ## if left, then we use gr.start of n2 node, if right we use gr.end of n2 node to make the loose end
                                     loose.nodes = parse.gr(edges[nix, ifelse(n1.side == 1, gr.string(gr.end(nodes[n1])), gr.string(gr.start(nodes[n1])))])
                                     ## now we fill in the n1 field with the index of the newly appended nodes
                                     edges[nix, n2 := length(nodes) + 1:.N]
                                     edges[nix, n2.side := ifelse(n1.side == 1, 0, 1)]
                                     edges[nix, type := 'loose']
                                     nodes = grbind(nodes, loose.nodes)
                                 }

                                 ## Set all the new loose nodes to loose
                                 nodes$loose = ifelse(is.na(nodes$loose), TRUE, FALSE)
                             }

                             strand(nodes) = '+'
                             nodes$node.id = 1:length(nodes) ## for reverse compatibility, to get rid
                             names(nodes) = NULL
                             segs = c(nodes, gr.flipstrand(nodes))
                             segs$snode.id = ifelse(as.logical(strand(segs)=='+'), 1, -1)*segs$node.id
                             private$pnodes = segs

                             if (nrow(edges)>0)
                             {
                                 ## n1.side, n2.side --> strand combo
                                 ## 1, 0 --> n1+,n2+ and n2-,n1-
                                 ## 0, 1 --> n1-,n2- and n2+,n1+
                                 ## 1, 1 --> n1+,n2- and n2+,n1-
                                 ## 0, 0 --> n1-,n2+ and n2-,n1+
                                 
                                 map = data.table(pid = private$pnodes$node.id, str = as.character(strand(private$pnodes)), id = 1:length(private$pnodes))
                                 setkeyv(map, c('pid', 'str'))
                                 edges[, jid := 1:.N]

                                 ## See if there is a cn in the edges, in which case we should carry it over
                                 is.cn = 'cn' %in% names(edges)
                                 is.type = 'type' %in% names(edges)

                                 tmp = rbind(
                                     edges[, .(jid, from = ifelse(n1.side == 1, map[.(n1, '+'), id], map[.(n1, '-'), id]),
                                               to = ifelse(n2.side == 1, map[.(n2, '-'), id], map[.(n2, '+'), id]),
                                               sedge.id = 1:.N)],
                                     edges[, .(jid, from = ifelse(n2.side == 1, map[.(n2, '+'), id], map[.(n2, '-'), id]),
                                               to = ifelse(n1.side == 1, map[.(n1, '-'), id], map[.(n1, '+'), id]),
                                               sedge.id = -1*(1:.N))]
                                 )
                                 tmp = merge(tmp, edges, by = "jid")

                                 tmp[, type := 'aberrant']
                                 tmp[, edge.id := abs(sedge.id)]

                                 
                                 ## Drop unused columns
                                 tmp[, c("jid","n1","n2","n1.side","n2.side") := NULL]
                                 private$pedges = tmp
                             }

                             private$pgraph = igraph::make_directed_graph(t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes))
                             return(self)
                         }
                     ),

                     
                     active = list(
                         nodes = function()
                         {
                             return(Node$new(private$pnodes$snode.id, self))
                         },

                         fullnodes = function() {
                             return(private$pnodes)
                         },

                         fulledges = function() {
                             return(copy(private$pedges))
                         }
                         
                     )
                     )
