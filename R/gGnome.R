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
#' gNodes and edges are necessary and sufficient to define a \code{gGraph} instance, while
#' optional metadata fields like copy numbers, edge attributes can be extended.
#'
#' In the following examples \code{gg} is a gGraph object.
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


## ================= gNode class definition ============= ##
#' @export
gNode = setClass("gNode")

gNode = R6::R6Class("gNode",
                   public = list(

                     ## Set up the constructor GNODE
                     ## index - snode.id
                     ## graph
                     initialize = function(snode.id = NULL, graph)
                     {
                       private$pgraph = graph
                       private$porientation = private$prindex = private$pindex = c()                       

                       if (is.null(snode.id))
                         return(self)

                       if (!all(snode.id %in% graph$gr$snode.id)) {
                         stop("Not valid snode.id(s)")
                       }
                       
                       private$pnode.id = abs(snode.id)
                       
                       lookup = private$pgraph$queryLookup(snode.id)
                       private$pindex = lookup[, index]
                       private$prindex = lookup[, rindex]

                       private$porientation = sign(snode.id)
                     },

                     flip = function() {
                       tmp = private$index
                       private$pindex = private$prindex
                       private$prindex = tmp
                       private$porientation = -private$porientation
                       return(self)
                     },                     

                     subset = function(i)
                     {
                       i = with(self$dt, eval(i)) ## allows subsetting based on metadata

                       if (is.logical(i))
                         i = which(i)

                       if (is.numeric(i) | is.integer(i))
                       {
                         if (any(i<0) | max(i, na.rm = TRUE)>self$length())
                           stop('index out of bounds')
                       }
                       
                       private$pnode.id = private$pnode.id[i]
                       private$pindex = private$pindex[i]
                       private$prindex = private$prindex[i]
                       private$porientation = private$porientation[i]

                       return(self)
                     },

                     print = function()
                     {
                       message('gNode object of length ', self$length())
                       print(self$gr)
                     },
                     
                     length = function()
                     {
                       return(length(private$pnode.id))
                     }
                   ),

                   
                   ## NODES PRIVATE FIELDS
                   private = list(
                     ## stores the ID of this node in the gGraph
                     pnode.id = NULL,

                     ## Stores the index and rindex (reverse index) of the pnode.id in this graph
                     pindex = NULL,

                     prindex = NULL,

                     ## Stores the sign of the node associated with index
                     porientation = NULL,
                     
                     ## Pointer to gGraph that this node belongs to 
                     pgraph = NULL
                   ),


                   ## NODES ACTIVE BINDINGS
                   active = list(
                     dt = function() {
                       return(as.data.table(private$pgraph$gr[private$pindex]))
                     },

                     gr = function() {
                       return(private$pgraph$gr[private$pindex])
                     },

                     id  = function()
                     {
                       return(private$pnode.id)
                     },
                     
                     sign = function()
                     {
                       return(private$porientation)
                     },
                     
                     left = function()
                     {
                       
                       ## Since there is only one index, get its left gNodes
                       tmp = private$pgraph$edgesdt[to %in% private$pindex, from]
                       
                       ## Catch for if we are on a loose end
                       leftNodes = gNode$new(snode.id =
                                              private$pgraph$gr$snode.id[tmp],
                                            graph = private$pgraph)
                       return(leftNodes)
                     },                                                      

                     ## Returns the nodes connected to the right of the nodes
                     right = function()
                     {
                       tmp = private$pgraph$edgesdt[from %in% private$pindex, to]
                       rightNodes = gNode$new(snode.id =
                                               sort(private$pgraph$gr$snode.id[tmp]),
                                             graph = private$pgraph)                       
                       return(rightNodes)
                     },                     
                     
                     ## FIXME: Returns a gGraph containing the nodes in this gNode Class
                     subgraph = function()
                     {
                       edges = private$pgraph$edgesdt[to %in% private$pindex & from %in% private$pid]
                       edgeObj = gEdge$new(snode.id = edges[, sedge.id],
                                          graph = private$pgraph)
                       
                       return(gGraph$new(gNodeObj = self, gEdgeObj = edgeObj))
                     }                                            
                   )
                   )


## ================= gEdge class definition ============= ##
#' @export
gEdge = setClass("gEdge")

gEdge = R6::R6Class("gEdge",
                   public = list(
                     ## Constructor GEDGE
                     ## Builds an gEdge Class from a given edge.id 
                     ## graph - the gGraph you want to build this class from 
                     initialize = function(seid, graph)
                     {
                       private$pgraph = graph ## reference to graph that these edges belong to
                       private$porientation = private$pedge.id = private$psedge.id = c()                       

                       if (is.null(seid))
                         return(self)

                       if (any(is.na(private$pgraph$edgesdt[.(seid), edge.id])))
                         stop("one or more provided signed edge ids are out of bounds")

                       private$porientation = 1
                       private$psedge.id = seid ## signed edge id, referring to the edge in the directed graph for which the $left direction is 5' and $right direction is '
                       private$pedge.id = abs(seid) ## unsigned edge id
                       private$pedges = private$pgraph$edgesdt[list(seid), ]
                       return(self)
                     },

                     ## Returns the number of edge pairs in this class
                     length = function()
                     {
                       return(length(private$pedge.id))
                     },

                     subset = function(i)
                     {
                       i = with(self$dt, eval(i)) ## allows subsetting based on metadata

                       if (is.logical(i))
                         i = which(i)

                       if (is.numeric(i) | is.integer(i))
                       {
                         if (any(i<0) | max(i, na.rm = TRUE)>self$length())
                           stop('index out of bounds')
                       }
                       
                       private$pedges = private$pedges[i]
                       private$pedge.id = private$pedge.id[i]
                       private$psedge.id = private$psedge.id[i]
                       private$porientation = private$porientation[i]
                       return(self)
                     },


                     ## Prints out the number of edges and the count of each type
                     print = function()
                     {
                       message("gEdge object with ", self$length(), " edges")                       
                       print(self$dt)
                     }                       
                   ),

                   
                   private = list(
                     ## data.table of the edges in this class
                     pedges = NULL,

                     pedge.id = NULL,
                     psedge.id = NULL,
                     
                     ## Defaults to positive, switch to negative if you flip it
                     porientation = NULL,

                     ## Pointer to the gGraph
                     pgraph = NULL
                   ),


                   active = list(

                     ## returns junctions representation of these edges
                     ## i.e. pair of signed intervals pointing away from the fusion
                     junctions = function()
                     {
                     },

                     ## every inter-edge in the bidirected graph is connection of the form a <-> b
                     ## in the directed graph each inter edge is represented by
                     ## two edges a->b (+ sign), b_->a_ (- sign) versions of that edge
                     ## for the + signed edge a is "left" and b is "right"
                     ## for the - signed edge b_ is "left" and a_ is "right"

                     left = function() {                      
                         tmp = private$pedges[edge.id == x, from]
                         index = which.min(start(private$pgraph$gr[tmp]))
                         leftNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[index],
                                              graph = private$pgraph)                                                
                         return(leftNodes)
                     },

                     
                     ## Returns the nodes connected to the right of the nodes
                     right = function() {                                            
                       rightNodes = private$pedges[, to]
                       return(gNode$new(snode.id = private$pgraph$private$ppnodes$snode.id[rightNodes],
                                       private$pgraph))
                     },
                     
                     ## Returns a data.table edges in this class format (to, from, type, edge.id)
                     dt = function() {
                       return(copy(private$pedges))
                     },

                     ## Returns a grangeslist representation of each of the edges with metadata populated
                     grl = function() {
                       gr1 = gr.flipstrand(gr.end(private$pgraph$gr[private$pedges$from], ignore.strand = FALSE))
                       gr2 = gr.end(private$pgraph$gr[private$pedges$to], ignore.strand = FALSE)
                       grl = split(c(gr1, gr2), rep(1:length(gr1), 2))[as.character(1:length(gr1))]
                       names(grl) = private$edges$sedge.id
                       values(grl) = private$pedges
                       return(grl)
                     },

                     ## Returns the edge.id's of the edges in the table
                     id = function() {
                       return(private$pedge.id)
                     },
                                          
                     ## Returns the subgraph associated with the edges in this class
                     subgraph = function()
                     {
                       
                     }
                   )
                   )


#' @name [.gNode
#' @title gNode
#' @description
#'
#' Overloads subset operator for nodes
#'
#' @param obj gNode object This is the gNode object to be subset
#' @param i  integer, logical, or expression in gNode metadata used to subset gEdges
#' @return A new node object that contains only the given id's
#' @export
'[.gNode' = function(obj, i = NULL){
  nodes = obj$clone()
  nodes$subset(substitute(i))
  return(nodes)
}

#' @name [.gEdge
#' @title gEdge
#' @description
#'
#' Overloads subset operator for edges
#'
#' @param obj gEdge object This is the gEode object to be subset
#' @param i integer, logical, or expression in  gEdge metadata used to subset gEdges
#' @return A new node object that contains only the given id's
#' @export
'[.gEdge' = function(obj, i = NULL){
  edges = obj$clone()
  edges$subset(substitute(i))
  return(edges)
}


gGraph = R6::R6Class("gGraph",
                     public = list(

                       ## public fields
                       ## name = NULL,
                       ## refG = "GENOME", ## seqinfo of ref genome
                       ## constructor INIT GGRAPH
                       initialize = function(genome = NULL,
                                             nodes = NULL,
                                             edges = NULL,
                                             nodesObj = NULL,
                                             edgesObj = NULL,
                                             looseterm = TRUE)
                       {
                         ## control how to construct
                         verbose = getOption("gGnome.verbose")
                         if (!is.null(nodes)){
                           private$gGraphFromNodes(nodes, edges, looseterm = looseterm)
                         } else if (!is.null(nodesObj))
                         {
                           private$gGraphFromNodesObj(nodesObj, edgesObj, looseterm = looseterm)
                         } else
                         {
                           private$emptyGGraph(genome)
                         }

                         ## FIXME: MOVE THIS GUY INTO CONSTRUCTORS, TMP FIX
                         private$buildLookupTable()
                       },

                       ## snode.id is a vector of signed node.id's from our gGraph
                       queryLookup = function(id) {
                         dt = private$lookup[.(id)]
                         return(dt)
                       },

                       length = function()
                       {
                         return(self$nodes$length())
                       },

                       print = function() {
                         cat("gGraph with ", self$length(), " nodes and ", nrow(private$pedges), " edges")
                         if (self$nodes$length()>00)
                           cat(', containing:\n')
                         else
                           cat('\n')                         
                         self$nodes$print()

                         if (length(self$edges)>0)
                         {
                           message()
                           self$edges$print()
                         }
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
                       ##     rindex - index corresponding to the complement snode.id
                       buildLookupTable = function() {
                         private$lookup = data.table()
                         if (length(private$pnodes)==0)
                           return(self)

                         lt = match(private$pnodes$snode.id, -private$pnodes$snode.id)
                         private$lookup = data.table(snode.id = private$pnodes$snode.id,
                                                     index = seq_along(private$pnodes),
                                                     rindex = lt)
                         setkey(private$lookup, snode.id)
                         return(self)
                       },


                       ## @name emptyGGraph
                       ## @brief Constructor, initializes an empty gGraph object. If the user does not provide genome
                       ##        and using to update this class, this will try to inherit the current gGraph's seqinfo.
                       ##        See class documentation for usage.
                       ## @param genome Optional seqinfo if the user wants to
                       emptyGGraph = function(genome=NULL)
                       {
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
                         private$pnodes = GRanges(seqinfo = genome)
                         private$pgraph = igraph::make_empty_graph()
                         private$pedges = data.table(from=integer(0),
                                                     to=integer(0),
                                                     type=character(0))

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

                         if (any(is.na(seqlengths(nodes))))
                           nodes = gUtils::gr.fix(nodes)

                         if (length(loose.right)>0)
                         {
                           edges = rbind(edges, data.table(n1 = loose.right, n1.side = 1, n2 = NA_integer_, n2.side = NA_integer_), fill = TRUE)
                         }

                         if (length(loose.left)>0)
                         {
                           edges = rbind(edges, data.table(n2 = loose.left, n2.side = 0, n1 = NA_integer_, n1.side = NA_integer_), fill = TRUE)
                         }

                         ## add all loose end nodes
                         if (nrow(edges)>0)
                         {
                           edges[, type := 'aberrant']
                           
                           if (any(nix <- is.na(edges$n1)))
                           {
                             ## adding loose ends nodes for n1 based on n2 side which is  either left (n2.side == 0) or right (n2.side == 1)
                             ## if left, then we use gr.start of n2 node, if right we use gr.end of n2 node to make the loose end
                             loose.nodes = parse.gr(edges[nix, ifelse(n2.side == 1, gr.string(gr.end(nodes[n2])), gr.string(gr.start(nodes[n2])))], seqlengths = seqlengths(nodes))
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
                             loose.nodes = parse.gr(edges[nix, ifelse(n1.side == 1, gr.string(gr.end(nodes[n1])), gr.string(gr.start(nodes[n1])))], seqlengths = seqlengths(nodes))
                             ## now we fill in the n1 field with the index of the newly appended nodes
                             edges[nix, n2 := as.numeric(length(nodes) + 1:.N)]
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
                         segs$index = 1:length(segs)
                         names(segs) = segs$snode.id
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

                           tmp[, edge.id := abs(sedge.id)]

                           
                           ## Drop unused columns
                           tmp[, c("jid","n1","n2","n1.side","n2.side") := NULL]
                           private$pedges = tmp
                           setkey(private$pedges, sedge.id)
                         }
                         
                         private$pgraph = igraph::make_directed_graph(t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes))
                         return(self)
                       },

                       
                       ## Takes nodes and edges and converts the edge table for the nodes if they were strandless
                       ## gEdge table will be the appropriate table for if we did:
                       ##      nodes[which(as.logical(strand(nodes)=="+"))]
                       ##      strand(nodes) = "*"
                       convertEdges = function(nodes, edges, metacols = FALSE)
                       {
                         ## Check to make sure we have some edge table, if not return empty
                         if(is.null(edges) || nrow(edges) == 0) {
                           return(data.table(n1 = integer(), n2 = integer(), n1.side = numeric(), n2.side = numeric()))
                         }

                         es = copy(edges)
                         
                         ## Map between to/from and n1.side, n2.side
                         ##    to (+) ---- n2.side = 0
                         ##    to (-) ---- n2.side = 1
                         ##    from (+) -- n1.side = 1
                         ##    from (-) -- n1.side = 0
                         es[, ":="(n1.side = ifelse(as.logical(strand(nodes[from]) == "+"), 1, 0),
                                   n2.side = ifelse(as.logical(strand(nodes[to]) == "+"), 0, 1))]

                         ## Get positive non-loose nodes
                         new.nodes = nodes %Q% (loose == FALSE & strand == "+")

                         ## Create map between old id positions and new positions (pos - pos, neg - pos, loose - NA)
                         ## Assumes no two values are not NA which they shouldnt be if everything is set up right on backend
                         indexMap = pmin(match(nodes$snode.id, new.nodes$snode.id), match(nodes$snode.id, -new.nodes$snode.id), na.rm=T)

                         ## Map the edges to their new location
                         es[, ":="(n1 = indexMap[from], n2 = indexMap[to])]

                         ## Get only the positive non-NA nodes
                         es = es[sedge.id > 0]
                         es = es[!is.na(n1) & !is.na(n2)]
                         
                         ## If the user chooses to keep metacols, just remove to and from, otherwise, only keep the essentials
                         if (metacols) {
                           es[, c("to","from") := NULL]
                           return(es)
                         } else {
                           return(es[, .(n1, n2, n1.side, n2.side,
                                         cn = if('cn' %in% names(es)) cn,
                                         type = if('type' %in% names(es)) type)])
                         }
                       },

                       
                       gGraphFromNodeClass = function(NodeObj, EdgeObj = NULL)
                       {
                         edges = NULL

                         edgesdt = convertEdges(NodeObj$gr, EdgeObj$dt)
                         
                         ## Get the granges associated with our NodeObj
                         indicies = self$queryLookup(NodeObj$id)
                         nodes = private$pnodes[indicies]
                         
                         if (!is.null(EdgeObj)) {
                           edges = EdgeObj$my.edges ##FIXME: change accessor                                                                                                            
                         }                                                  
                       }
                       
                     ),
                     
                     active = list(
                       nodes = function()
                       {
                         return(gNode$new(private$pnodes$snode.id[private$pnodes$snode.id>0], self))
                       },
                       edges = function()
                       {
                         return(gEdge$new(private$pedges$sedge.id[private$pedges$sedge.id>0], self))
                       },
                       
                       gr = function() {
                         return(private$pnodes)
                       },

                       dt  = function() {
                         return(as.data.table(private$pnodes))
                       },
                       
                       edgesdt = function() {
                         return(copy(private$pedges))
                       }                       
                     )
                     )

