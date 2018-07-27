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
library(GenomicRanges)
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


                            if (is.null(snode.id)) {
                                return(self)
                            }
                            
                            if (!all(snode.id %in% graph$gr$snode.id)) {
                                stop("Not valid snode.id(s)")
                            }
                            
                            private$pnode.id = abs(snode.id)
                            
                            lookup = private$pgraph$queryLookup(snode.id)
                            private$pindex = lookup[, index]
                            private$prindex = lookup[, rindex]

                            private$porientation = sign(snode.id)
                        },

                        annotate = function(col1, vec){
                            private$pgraph$annotate(col1, vec, TRUE)
                        },
                        
                        c.gNode= function(...){                            
                            gNode.list=list(...)
                            isg = sapply(gNode.list, function(x) class(x)[1]=='gNode')

                            if(any(!isg)){
                                stop('Error: All inputs must be of class gNode.')
                            }

                            ##Check to make sure they all come from the same graph
                            graphs = lapply(gNode.list, function(x) x$graph)
                            if(length(unique(graphs))>1){
                                stop('Error: All gNodes must come from the same graph')
                            }
                            ##Get all the pnode.id's to create new gNode
                            ids = lapply(gNode.list, function(x) x$id)
                            return (gNode$new(unlist(ids), graphs[[1]]))
                            
                            
                        },

                        ## Flips the current orientation of the Node object (swap index/rindex)
                        flip = function() {
                            tmp = private$index
                            private$pindex = private$prindex
                            private$prindex = tmp
                            private$porientation = -private$porientation
                            return(self)
                        },


                        ## Allows subseting of the Node object using bracket notation
                        subset = function(i)
                        {
                            i = with(self$dt, eval(i)) ## allows subsetting based on metadata

                            if (is.logical(i))
                                i = which(i)


                            if (is.numeric(i) | is.integer(i)) {
                                if (any(i<0) | max(i, na.rm = TRUE)>self$length()) {
                                    stop('index out of bounds')
                                }
                            }
                            
                            private$pnode.id = private$pnode.id[i]
                            private$pindex = private$pindex[i]
                            private$prindex = private$prindex[i]
                            private$porientation = private$porientation[i]
                            
                            return(self)
                        },



                        ## Prints the Node Object
                        print = function()
                        {
                            message('gNode object of length ', self$length())
                            print(self$gr)
                        },


                        ## Returns the number of nodes in the Node Object
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
                        graph = function()
                        {
                            return(private$pgraph)
                        },

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
                        
                        ## Returns our signed node.id
                        sid = function()
                        {
                            return(ifelse(private$porientation == "+", private$pnode.id, -private$pnode.id))
                        },


                        ## Returns the nodes connected to the left of the nodes
                        left = function()
                        {
                            ## Since there is only one index, get its left gNodes
                            tmp = private$pgraph$edgesdt[to %in% private$pindex, from]
                            
                            ## Catch for if we are on a loose end
                            leftNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[tmp],
                                                  graph = private$pgraph)
                            
                            return(leftNodes)
                        },                                                      

                        
                        ## Returns the nodes connected to the right of the nodes
                        right = function()
                        {
                            tmp = private$pgraph$edgesdt[from %in% private$pindex, to]
                            rightNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[tmp],
                                                   graph = private$pgraph)                       
                            return(rightNodes)
                        },
                        
                        
                        ## Returns the edges connected to the left of the nodes in this Node object as an Edge Object
                        eleft = function()
                        {
                            return(Edge$new(sedge.id = private$graph$fulledges[to %in% private$index, sedge.id],
                                            graph = graph))
                        },
                        

                        ## Returns the edges connected to the right of the nodes in this Node object as an Edge Object
                        eright = function()
                        {
                            return(Edge$new(sedge.id = private$graph$fulledges[from %in% private$index, sedge.id],
                                            graph = graph))
                        },


                        ## Creates a subgraph of all the Nodes in this Class and the edges that correspond to these nodes
                        subgraph = function()
                        {
                            ## Get all the edges that connect our nodes and remove the duplicate edge.id's as they will
                            ## be automatically filled in in the Edge Constructor
                            edge.id = private$pgraph$edgesdt[to %in% c(private$pindex, private$prindex)
                                                             & from %in% c(private$pindex, private$prindex), edge.id]
                            edge.id = edge.id[!duplicated(edge.id)]
                            
                            edgeObj = gEdge$new(seid = edge.id,
                                                graph = private$pgraph)
                            return(gGraph$new(nodeObj = self, edgeObj = edgeObj))
                        }
                    )
                    )


## ================ Non-Member Functions for gNode =============== ##

#' @name length
#' The number of nodes in the gNode Object
#'
#' @param gNode a gNode object
#'
#' @return the number of nodes in the gNode Object
#' 
#' @export
`length.gNode` = function(gNode)
{
    if(!inherits(gNode, "gNode")) {
        stop("Error: Invalid input")
    }
    return(gNode$length())
}



#' @name setdiff
#' Returns a new gNode object which is the difference between x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#'
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNode Object containing the difference between x and y
setMethod("setdiff", c("gNode", "gNode"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = setdiff(x$id, y$id)
              return(gNode$new(new.ids, x$graph))
          })



#' @name union
#' Returns a new gNode object which is the union of x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNode Object containing the union of x and y
setMethod("union", c("gNode", "gNode"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = union(x$id, y$id)
              return(gNode$new(new.ids, x$graph))
          })


#' @name intersect
#' Returns a new gNode object which is the intersection of x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNode Object containing the intersection of x and y
setMethod("intersect", c("gNode", "gNode"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = intersect(x$id, y$id)
              return(gNode$new(new.ids, x$graph))
          })


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
                            
                            if (is.null(seid)) {
                                return(self)
                            }

                            if (any(is.na(private$pgraph$edgesdt[.(seid), edge.id]))) {
                                stop("one or more provided signed edge ids are out of bounds")
                            }

                            private$porientation = 1
                            private$psedge.id = seid ## signed edge id, referring to the edge in the directed graph for which the $left direction is 5' and $right direction is '
                            private$pedge.id = abs(seid) ## unsigned edge id
                            private$pedges = private$pgraph$edgesdt[list(seid), ]
                            return(self)
                        },
                        
                        annotate=function(col1, vec){                            
                            private$pgraph$annotate(col1, vec, FALSE)
                        },
                        
                        ## Returns the number of edge pairs in this class
                        length = function()
                        {                            
                            return(length(private$pedge.id))
                        },
                        
                        ## Combine gEdge objects into a single gEdge
                        c.gEdge= function(...){                            
                             gEdge.list=list(...)
                             isg = sapply(gEdge.list, function(x) class(x)[1]=='gEdge')
                             
                            if(any(!isg)){
                                stop('Error: All inputs must be of class gEdge.')
                            }

                            ##Check to make sure they all come from the same graph
                            graphs = lapply(gEdge.list, function(x) x$graph)
                            if(length(unique(graphs))>1){
                                stop('Error: All gEdges must come from the same graph')
                            }
                            ##Get all the pnode.id's to create new gNode
                            ids = lapply(gEdge.list, function(x) x$id)
                            return (gEdge$new(unlist(ids), graphs[[1]]))                                                        
                        },
                        
                        ## Allows for subsetting of the Edge Object using bracket notation
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

                        graph = function()
                        {
                            return(private$pgraph)
                        },
                        
                        left = function()
                        {                      
                           ##tmp = private$pedges[edge.id == x, from]
                            leftNodes = private$pedges[, from]
                            ##index = which.min(start(private$pgraph$gr[tmp]))
##                            leftNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[index],
  ##                                                graph = private$pgraph)                                                
                            leftNodes = gNode$new(snode.id=private$pgraph$dt[, snode.id[leftNodes]], private$pgraph)
                            return(leftNodes)
                        },

                        
                        ## Returns the nodes connected to the right of the nodes
                        right = function()
                        {
                            
                            rightNodes = private$pedges[, to]
                            return(gNode$new(snode.id = private$pgraph$dt[,snode.id[rightNodes]],
                                             private$pgraph))
                        },

                        
                        ## Returns a data.table edges in this class format (to, from, type, edge.id)
                        dt = function()
                        {
                            return(copy(private$pedges))
                        },

                        
                        ## Returns a grangeslist representation of each of the edges with metadata populated
                        grl = function()
                        {                            
                            gr1 = gr.flipstrand(gr.end(private$pgraph$gr[private$pedges$from], ignore.strand = FALSE))
                            gr2 = gr.end(private$pgraph$gr[private$pedges$to], ignore.strand = FALSE)
                            grl = split(c(gr1, gr2), rep(1:length(gr1), 2))[as.character(1:length(gr1))]
                            names(grl) = private$edges$sedge.id
                            values(grl) = private$pedges
                            return(grl)
                        },

                        
                        ## Returns the edge.id's of the edges in the table
                        id = function()
                        {
                            return(private$pedge.id)
                        },

                        
                        ## Returns the subgraph associated with the edges in this class
                        subgraph = function()
                        {
                            ## Get all of the nodes in either to or from
                            nodeIndex = c(private$pedges[,to], private$pedges[,from])

                            ## Get the unsigned node.id's and remove any duplicate nodes
                            node.id = private$pgraph$dt[index %in% nodeIndex, node.id]
                            node.id = node.id[!duplicated(node.id)]

                            nodeObj = gNode$new(snode.id = node.id,
                                                graph = private$pgraph)
                            return(gGraph$new(nodeObj = nodeObj, edgeObj = self))
                        }
                        

                    )
                    )


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


#' @name length
#' The number of edge pairs in the gEdge Object
#'
#' @param gNode a gEdge object
#'
#' @return the number of edge pairs in the gEdge Object
#' 
#' @export
`length.gEdge` = function(gEdge)
{
    if(!inherits(gEdge, "gEdge")) {
        stop("Error: Invalid input")
    }
    return(gEdge$length())
}



#' @name setdiff
#' Returns a new gNode object which is the difference between x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#'
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNodeObject containing the difference between x and y
setMethod("setdiff", c("gEdge", "gEdge"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = setdiff(x$id, y$id)
              return(gEdge$new(new.ids, x$graph))
          })


#' @name union
#' Returns a new gNode object which is the union of x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNode Object containing the union of x and y
setMethod("union", c("gEdge", "gEdge"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = union(x$id, y$id)
              return(gEdge$new(new.ids, x$graph))
          })


#' @name intersect
#' Returns a new gNode object which is the intersection of x and y (id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#'
#' @return new gNode Object containing the intersection of x and y
setMethod("intersect", c("gEdge", "gEdge"),
          function(x, y) {
              ## Make sure that both come from the same graph
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = intersect(x$id, y$id)
              return(gEdge$new(new.ids, x$graph))
          })


gGraph = R6::R6Class("gGraph",
                     public = list(


                         ## public fields
                         ## name = NULL,
                         ## refG = "GENOME", ## seqinfo of ref genome
                         ## constructor INIT GGRAPH
                         initialize = function(genome = NULL,
                                               prego=NULL,
                                               jabba = NULL,
                                               cougar=NULL,
                                               nodes = NULL,                                               
                                               edges = NULL,
                                               nodeObj = NULL,
                                               edgeObj = NULL,
                                               looseterm = TRUE)
                         {
                             ## control how to construct
                             verbose = getOption("gGnome.verbose")
                             if (!is.null(nodes)){
                                 private$gGraphFromNodes(nodes, edges, looseterm = looseterm)
                             }
                             else if (!is.null(nodeObj))
                             {
                                 private$gGraphFromNodeObj(nodeObj, edgeObj, looseterm = looseterm)
                             }                             
                             else if(!is.null(prego)){
                                 if(verbose){
                                     message("Reading Prego output")

                                 }
                                 private$pr2gg(prego, looseterm)
                                 }                             
                             
                             else if (!is.null(jabba))
                             {
                                 private$jab2gg(jabba)
                             }
                             else if (!is.null(cougar))
                             {
                                 private$cougar2gg(cougar)
                             }
                             else
                             {
                                 private$emptyGGraph(genome)
                             }
                             
                             return(self)
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
                         },
                         annotate = function(col1, vec, node=FALSE){                                                          
                             if (node==TRUE){                                 
                                 gr.dt=gr2dt(private$pnodes)
                                 gr.dt[, paste(col1):=vec]
                                 private$pnodes=makeGRangesFromDataFrame(gr.dt, TRUE)
                             }
                             else{                                
                                 private$pedges[, paste(col1:=vec)]
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

                             if (length(private$pnodes) == 0) {
                                 return(self)
                             }
                             
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
                             
                             ## FIXME: MOVE THIS GUY INTO CONSTRUCTORS, TMP FIX
                             private$buildLookupTable()

                             return(self)
                         },

                         
                         ## Takes nodes and edges and converts the edge table for the nodes if they were strandless
                         ## gEdge table will be the appropriate table for if we did:
                         ##      nodes[which(as.logical(strand(nodes)=="+"))]
                         ##      strand(nodes) = "*"
                         ## pre: Nodes have snode.id and index column (indicates index)
                         ##      edges have sedge.id 
                         ## 
                         convertEdges = function(nodes, edges, metacols = FALSE)
                         {                            
                             ## Check to make sure we have some edge table, if not return empty
                             if(is.null(edges) || nrow(edges) == 0) {
                                 return(data.table(n1 = integer(), n2 = integer(), n1.side = numeric(), n2.side = numeric()))
                             }

                             es = copy(edges)
                             nodedt = gr2dt(nodes)
                             setkey(nodedt, index)
                             
                             ## Map between to/from and n1.side, n2.side
                             ##    to (+) ---- n2.side = 0
                             ##    to (-) ---- n2.side = 1
                             ##    from (+) -- n1.side = 1
                             ##    from (-) -- n1.side = 0
                             
                             es[, ":="(n1.side = ifelse(nodedt[.(es[,from]), strand] == "+", 1, 0),
                                       n2.side = ifelse(nodedt[.(es[,to]), strand] == "+", 0, 1))]


                             ## Get positive non-loose nodes
                             new.nodes = nodes %Q% (loose == FALSE & strand == "+")
                             
                             ## Create map between old id positions and new positions (pos - pos, neg - pos, loose - NA)
                             ## Assumes no two values are not NA which they shouldnt be if everything is set up right on backend
                             indexMap = pmin(match(nodes$snode.id, new.nodes$snode.id), match(nodes$snode.id, -new.nodes$snode.id), na.rm=T)
                             indexMap = data.table(index = nodes$index, new.index = indexMap)
                             setkey(indexMap, index)
                             
                             ## Map the edges to their new location
                             es[, ":="(n1 = indexMap[.(from), new.index],
                                       n2 = indexMap[.(to), new.index])]

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


                         ## Adds the proper snode.id and index to nodes, adds proper sedge.id/edge.id to edges
                         ## Returns list(nodes, edges) with updated fields and miscellaneous metadata removed
                         ## pre: loose column is in nodes, all nodes/edges are paired up
                         ## USE THIS FUNCTION WITH CAUTION - gGraphFromNodes will require that sedge.id/edge.id is removed before running (FIXME: don't make it do that)
                         pairNodesAndEdges = function(nodes, edges)
                         {
                             edges = as.data.table(edges)
                             not.loose = which(!nodes$loose)
                             ix = not.loose[match(gr.stripstrand(nodes[not.loose]), gr.stripstrand(nodes[not.loose]))]
                             nodes$snode.id = NA
                             nodes$snode.id[not.loose] = ifelse(as.logical(strand(nodes)[not.loose] == '+'), ix, -ix)
                             nodes$index = 1:length(nodes)
                             
                             edges[, tag := paste(nodes$snode.id[to], nodes$snode.id[from])]
                             edges[, rtag := paste(-nodes$snode.id[from], -nodes$snode.id[to])]
                             edges[, id := 1:.N]
                             edges[, rid := match(tag, rtag)]                                 
                             edges[, edge.id := clusters(graph.edgelist(cbind(id, rid)), 'weak')$membership]
                             edges[, sedge.id := ifelse(duplicated(edge.id), -edge.id, edge.id)]
                             edges = edges[, .(from, to,
                                               cn = if("cn" %in% names(edges)) cn,
                                               type = if("type" %in% names(edges)) type,
                                               edge.id, sedge.id)]

                             return(list(nodes, edges))
                         },


                         ## Operates only on positive strand, treats all nodes in NodeObj as strandless -- so if a node is marked
                         ## negative it will be treated as a duplicate positive node
                         ## Possibly remove them idk
                         gGraphFromNodeObj = function(NodeObj, EdgeObj = NULL, looseterm = TRUE)
                         {
                             ## Make sure that NodeObj and EdgeObj are gNode and gEdge respectively
                             if(!inherits(NodeObj, "gNode")) {
                                 stop("NodeObj is not a gNode Object")
                             }

                             if(!is.null(EdgeObj) && !inherits(EdgeObj, "gEdge")) {
                                 stop("EdgeObj is not a gEdge Object")

                             }
                             if(!identical(NodeObj$graph, EdgeObj$graph)) {
                                 stop("NodeObj and EdgeObj do not point to the same graph")
                             }

                             ## If we have no edges, just call the constructor on the Nodes Object that aren't loose
                             if(is.null(EdgeObj) || EdgeObj$length() == 0) {
                                 private$gGraphFromNodes(nodes = NodeObj[loose == FALSE]$gr, looseterm = looseterm)
                                 return(self)
                             }

                             ## If we have Node Obj but no EdgeObj, throw an error
                             if(NodeObj$length() == 0) {
                                 stop("Error: Edges cannot non-empty if nodes are empty")
                             }
                             
                             nodes = c(NodeObj$gr, NodeObj$clone()$flip()$gr)


                             ## Validate EdgeObj, no indices not present in index column of 
                             edges = EdgeObj$dt
                             if(!all(nix <- c(edges[,to], edges[,from]) %in% nodes$index)) {
                                 stop(paste0("Edge Object contains indices not in NodeObj. Indicies to remove are: ",
                                            paste(c(edges[,to], edges[,from])[!nix], collapse = " ")))

                             }

                             ## FIXME: this only works because of how our nodes are set up but might fail later, need to use $index within convertEdges
                             edges = private$convertEdges(nodes, edges)
                             nodes = NodeObj[loose == FALSE]$gr
                             
                             private$gGraphFromNodes(nodes = nodes, edges = edges, looseterm = looseterm)
                         },
cougar2gg = function(cougar){
                             "Convert the cougar output directory to gGraph."
                             if (!dir.exists(cougar)){
                                 stop("Error: invalid input CouGaR directory!")
                             }

                             if (!dir.exists(paste(cougar, 'solve',sep = '/'))){
                                 stop("No CouGaR solutions found in the input directory!")
                             }

                             .parsesol = function(this.sol)
                             {
                                 verbose = getOption("gGnome.verbose")
                                 tmp = unlist(.parseparens(this.sol[2]))
                                 tmp2 = as.data.table(
                                     matrix(tmp[nchar(stringr::str_trim(tmp))>0], ncol = 3, byrow = TRUE))
                                 segs = cbind(
                                     as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames = V1, start = V2)],
                                     data.table(end = as.numeric(sapply(strsplit(tmp2$V2, ' '), '[', 2)), strand = '+'),
                                     as.data.table(matrix(unlist(strsplit(stringr::str_trim(tmp2$V3), ' ')),
                                                          ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)])
                                 segs = suppressWarnings(dt2gr(segs))
                                 segs$id = 1:length(segs)
                                 nodes = c(segs, gr.flipstrand(segs))
                                 nodes$nid = ifelse(as.logical(strand(nodes) == '+'), 1, -1)*nodes$id
                                 nodes$ix = 1:length(nodes)
                                 nodes$rix = match(-nodes$nid, nodes$nid)
                                 adj = array(0, dim = rep(length(nodes),2))
                                 adj = sparseMatrix(length(nodes),length(nodes), x = 0)

                                 tmp = unlist(.parseparens(this.sol[3]))
                                 if (length(tmp)>0) ## are there any somatic edges?
                                 {
                                     tmp2 = as.data.table(matrix(tmp[nchar(str_trim(tmp))>0], ncol = 3, byrow = TRUE))
                                     abadj = cbind(
                                         as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames1 = V1, pos1 = V2)],
                                         as.data.table(matrix(unlist(strsplit(tmp2$V2, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames2 = V1, pos2 = V2)],
                                         as.data.table(matrix(unlist(strsplit(str_trim(tmp2$V3), ' ')),
                                                              ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)]
                                     )
                                     abadj$strand1 = ifelse(abadj$type %in% c(0,2), '+', '-')
                                     abadj$strand2 = ifelse(abadj$type %in% c(0,3), '+', '-')
                                     abadj$start.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), start(segs)))
                                     abadj$end.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), end(segs)))
                                     abadj$start.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), start(segs)))
                                     abadj$end.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), end(segs)))

                                     ## if strand1 == '+' then end match
                                     ## if strand1 == '-' then start match
                                     ## if strand2 == '+' then start match
                                     ## if strand2 == '-' then end match
                                     
                                     abadj[, match1 := ifelse(strand1 == '+', end.match1, -start.match1)]
                                     abadj[, match2 := ifelse(strand2 == '+', start.match2, -end.match2)]

                                     
                                     abadj[, nmatch1 := match(match1, nodes$nid)]
                                     abadj[, nmatch2 := match(match2, nodes$nid)]

                                     abadj[, nmatch1r := match(-match1, nodes$nid)]
                                     abadj[, nmatch2r := match(-match2, nodes$nid)]
                                     
                                     adj[cbind(abadj$nmatch1, abadj$nmatch2)] = abadj$cn
                                     adj[cbind(abadj$nmatch2r, abadj$nmatch1r)] = abadj$cn
                                 }

                                 ## how many node copies are unaccounted for by aberrant edges on left and right
                                 node.diff.in = nodes$cn - colSums(adj)
                                 node.diff.out = nodes$cn - rowSums(adj)

                                 norm.adj = as.data.table(cbind(1:length(segs), match(gr.end(segs), gr.start(segs))))[!is.na(V2), ]
                                 norm.adj = rbind(norm.adj, norm.adj[, .(V2 = -V1, V1 = -V2)])[, nid1 := match(V1, nodes$nid)][, nid2 := match(V2, nodes$nid)]

                                 ## now add non-aberrant edge copy numbers that are the minimum of the unaccounted
                                 ## for copy number going <out> of the source node and going <in> to the sink node
                                 adj.old = adj
                                 ## ALERT: extremely hacky solution
                                 adj[as.matrix(norm.adj[, .(nid1, nid2)])] =
                                     pmax(pmin(node.diff.out[norm.adj[, nid1]],
                                          node.diff.in[norm.adj[, nid2]]), 0)

                                 nodes$eslack.in = nodes$cn - colSums(adj)
                                 nodes$eslack.out = nodes$cn - rowSums(adj)


                                 if (sum(adj!=0)>0)
                                 {
                                     if (!identical(adj[which(adj>0)],
                                                    adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = nodes$rix[col], col = nodes$rix[row])])]))
                                     {
                                         stop('reciprocality violated')
                                     }
                                 }
                                 end(nodes) = end(nodes)-1

                                 return(list(nodes, as(adj, 'Matrix')))
                             }

                             .parseparens = function(str)
                             {
                                 cmd = gsub(',$', '',
                                            gsub(',\\)', ')',
                                                 gsub('\\)', '),',
                                                      gsub('\\(', 'list(',
                                                           gsub('([^\\(^\\[^\\]^\\)]+)', '"\\1",', perl = TRUE, gsub('\\]', ')', gsub('\\[', '\\(', str)))))))
                                 eval(parse(text = cmd))
                             }

                             sols = lapply(dir(dir(paste(cougar, 'solve',sep = '/'), full = TRUE)[1], '^g_', full = TRUE), readLines)
                             if (length(sols)==0){
                                 if (verbose){
                                     
                                 }
                                 return(self$nullGGraph())
                             }

                             ## parse cougar graphs
                             graphs = lapply(sols, .parsesol)

                             ## concatenate nodes and block diagonal bind adjacency matrices
                             segs = do.call('c', lapply(graphs, '[[', 1))
                             segs$id = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$id, sep = '.')
                             segs$nid = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$nid, sep = '.')
                             segs$ix = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$ix, sep = '.')
                             segs$rix = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$rix, sep = '.')
                             segs$rix = match(segs$rix, segs$ix)
                             segs$ix = 1:length(segs)
                             adj = do.call('bdiag', lapply(graphs, '[[', 2))


                             ## final double check for identicality
                             if (!(identical(adj[which(adj>0)], adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = segs$rix[col], col = segs$rix[row])])])))
                             {
                                 stop('Reciprocality check failed!')
                             }

                             ## TODO: figure out why there are negative CN edges in CouGaR!!!
                             private$gGraphFromScratch(segs = segs, es = adj)$fillin() ## IT MAY BE NOT BALANCED
                             return(self)                         
                         },                         


 
                         pr2gg = function(fn, looseterm = TRUE)
                         {
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
                              
                              ## turn into our nodes
                              posNodes = GRanges(res[[1]]$chr1,
                                                 IRanges(res[[1]]$pos1,
                                                 res[[1]]$pos2),
                                                 strand = "+",
                                                 cn = res[[1]]$cn,
                                                 left.tag = res[[1]]$node1,
                                                 right.tag = res[[1]]$node2,
                                                 loose=FALSE)
                              
                              ## Set up our new nodes
                              posNodes$snode.id = 1:length(posNodes)
                              nodes = gr.fix(c(posNodes, gr.flipstrand(posNodes)), sl)
                              neg.ix = which(as.logical(strand(nodes) == "-"))
                              nodes$snode.id[neg.ix] = -1 * nodes$snode.id[neg.ix]
                              nodes$index=1:length(nodes)
                              ## tag1 is the 3' end
                              tag1 = nodes$right.tag
                              tag1[neg.ix] = nodes$left.tag[neg.ix]
                              ## tag2 is the 5' end
                              tag2 = nodes$left.tag
                              tag2[neg.ix] = nodes$right.tag[neg.ix]
                              
                              ## adjacency in copy number
                         adj.cn = matrix(0, nrow = length(nodes), ncol = length(nodes),
                                         dimnames = list(tag1, tag2))
                              adj.cn[cbind(res[[2]]$node1, res[[2]]$node2)] = res[[2]]$cn
                              adj.cn[cbind(res[[2]]$node2, res[[2]]$node1)] = res[[2]]$cn
                              adj.cn[cbind(res[[3]]$node1, res[[3]]$node2)] = res[[3]]$cn
                              adj.cn[cbind(res[[3]]$node2, res[[3]]$node1)] = res[[3]]$cn
                              
                              ## create es
                              ed = as.data.table(which(adj.cn>0, arr.ind=T))
                              colnames(ed) = c("from", "to")
                              ed[, ":="(cn = adj.cn[cbind(from, to)])]
                              edges=private$pairNodesAndEdges(nodes, ed)[[2]]
                              edges = private$convertEdges(nodes, edges)                              
                              ## FIXME: update constructor header with looseterm
                              private$gGraphFromNodes(nodes = granges(posNodes), edges = edges, looseterm = looseterm)
                              
                              return(self)
                          },
                         
                         
                         ## @name jab2gg
                         ## @brief Constructor for creating a gGraph object from a jabba output
                         jab2gg = function(jabba, regular=FALSE, looseterm = TRUE)
                         {
                             ## Validate our input
                             if (is.list(jabba)) {
                                 if (!all(is.element(c("segstats", "adj",
                                                       "purity", "ploidy"),
                                                     names(jabba))))
                                     stop("The input is not a JaBbA output.")
                             } else if (is.character(jabba) & grepl(".rds$", jabba)){
                                 if (file.exists(jabba)){
                                     jabba = readRDS(jabba)
                                 }
                             } else {
                                 stop("Error: Input must be either JaBbA list output or the RDS file name that contains it!")
                             }

                             
                             nodes = jabba$segstats
                             edges = NULL
                             
                             ## If we have edges, set them up for constructor
                             if ("edges" %in% names(jabba)){
                                 ## Make sure edges is a data.table
                                 edges = as.data.table(jabba$edges)[type != "loose", .(from, to, cn, type)]
                                 
                                 paired = private$pairNodesAndEdges(nodes, edges)
                                 nodes = paired[[1]]
                                 edges = paired[[2]]
                                 
                                 ## Convert edges to the proper form (n1,n2,n1.side,n2.side
                                 edges = private$convertEdges(nodes, edges)
                             }
                             
                             ## We want to get only the positive strand of nodes
                             nodes = nodes %Q% (loose == FALSE & strand == "+")
                             
                             ## Need to get seqinfo, if all of the seqlenths are missing, fill them in
                             if (any(is.na(seqlengths(nodes)))){
                                 
                                 ## Try to use junction seqlenths
                                 if (!any(is.na(seqlengths(jabba$junctions)))){
                                     
                                     if (getOption("gGnome.verbose")){
                                         warning("No valid seqlengths found in $segstats, force to the same as $junctions.")
                                     }
                                     nodes = gUtils::gr.fix(nodes, jabba$junctions)

                                 } else {
                                     
                                     if (getOption("gGnome.verbose")){
                                         warning("No valid seqlengths found anywhere in input, force to DEFAULT.")
                                     }
                                     default.sl = data.table::fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]

                                     nodes = gUtils::gr.fix(nodes, default.sl)
                                 }
                             }
                             
                             ## Set up gGraph using other constructor, standardizes our input although not totally necessary
                             private$gGraphFromNodes(nodes, edges, looseterm=TRUE)
                             
                             ## If regular, remove non-regular nodes from nodes by trimming
                             if (regular==T) {

                                 if (getOption("gGnome.verbose")){
                                     warning("Forcing regular chromosomes. Will try default. See `Sys.getenv('DEFAULT_REGULAR_CHR')`.")
                                 }
                                 regularChr = si2gr(data.table::fread(Sys.getenv('DEFAULT_REGULAR_CHR'))[, setNames(V2, V1)])
                                 self$trim(regularChr, mod=T)
                             }

                             return(self)
                         }                         

                     ),
                     
                     
                     
                     active = list(
                         ## Returns a Node Object of the nodes in the graph
                         nodes = function()
                         {
                             return(gNode$new(private$pnodes$snode.id[private$pnodes$snode.id>0], self))
                         },



                         ## Returns an Edge Object of the edges in the graph
                         edges = function()
                         {
                             return(gEdge$new(private$pedges$sedge.id[private$pedges$sedge.id>0], self))
                         },

                         
                         ## Returns all of the nodes in the graph as a GRanges
                         gr = function() {
                             return(private$pnodes)
                         },

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
                                 gmessage(chrm.map.fn)
                                 gmessage("Seqnames mapping found.")
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

cougar2gg = function(cougar){
                             "Convert the cougar output directory to gGraph."
                             if (!dir.exists(cougar)){
                                 stop("Error: invalid input CouGaR directory!")
                             }

                             if (!dir.exists(paste(cougar, 'solve',sep = '/'))){
                                 stop("No CouGaR solutions found in the input directory!")
                             }

                             .parsesol = function(this.sol)
                             {
                                 verbose = getOption("gGnome.verbose")
                                 tmp = unlist(.parseparens(this.sol[2]))
                                 tmp2 = as.data.table(
                                     matrix(tmp[nchar(stringr::str_trim(tmp))>0], ncol = 3, byrow = TRUE))
                                 segs = cbind(
                                     as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames = V1, start = V2)],
                                     data.table(end = as.numeric(sapply(strsplit(tmp2$V2, ' '), '[', 2)), strand = '+'),
                                     as.data.table(matrix(unlist(strsplit(stringr::str_trim(tmp2$V3), ' ')),
                                                          ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)])
                                 segs = suppressWarnings(dt2gr(segs))
                                 segs$id = 1:length(segs)
                                 nodes = c(segs, gr.flipstrand(segs))
                                 nodes$nid = ifelse(as.logical(strand(nodes) == '+'), 1, -1)*nodes$id
                                 nodes$ix = 1:length(nodes)
                                 nodes$rix = match(-nodes$nid, nodes$nid)
                                 adj = array(0, dim = rep(length(nodes),2))
                                 adj = sparseMatrix(length(nodes),length(nodes), x = 0)

                                 tmp = unlist(.parseparens(this.sol[3]))
                                 if (length(tmp)>0) ## are there any somatic edges?
                                 {
                                     tmp2 = as.data.table(matrix(tmp[nchar(str_trim(tmp))>0], ncol = 3, byrow = TRUE))
                                     abadj = cbind(
                                         as.data.table(matrix(unlist(strsplit(tmp2$V1, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames1 = V1, pos1 = V2)],
                                         as.data.table(matrix(unlist(strsplit(tmp2$V2, ' ')), ncol = 2, byrow = TRUE))[, .(seqnames2 = V1, pos2 = V2)],
                                         as.data.table(matrix(unlist(strsplit(str_trim(tmp2$V3), ' ')),
                                                              ncol = 4, byrow = TRUE))[, .(type = V1, cn = as.numeric(V2), ncov = V3, tcov  = V4)]
                                     )
                                     abadj$strand1 = ifelse(abadj$type %in% c(0,2), '+', '-')
                                     abadj$strand2 = ifelse(abadj$type %in% c(0,3), '+', '-')
                                     abadj$start.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), start(segs)))
                                     abadj$end.match1 = match(abadj[, paste(seqnames1, pos1)], paste(seqnames(segs), end(segs)))
                                     abadj$start.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), start(segs)))
                                     abadj$end.match2 = match(abadj[, paste(seqnames2, pos2)], paste(seqnames(segs), end(segs)))

                                     ## if strand1 == '+' then end match
                                     ## if strand1 == '-' then start match
                                     ## if strand2 == '+' then start match
                                     ## if strand2 == '-' then end match
                                     
                                     abadj[, match1 := ifelse(strand1 == '+', end.match1, -start.match1)]
                                     abadj[, match2 := ifelse(strand2 == '+', start.match2, -end.match2)]

                                     
                                     abadj[, nmatch1 := match(match1, nodes$nid)]
                                     abadj[, nmatch2 := match(match2, nodes$nid)]

                                     abadj[, nmatch1r := match(-match1, nodes$nid)]
                                     abadj[, nmatch2r := match(-match2, nodes$nid)]
                                     
                                     adj[cbind(abadj$nmatch1, abadj$nmatch2)] = abadj$cn
                                     adj[cbind(abadj$nmatch2r, abadj$nmatch1r)] = abadj$cn
                                 }

                                 ## how many node copies are unaccounted for by aberrant edges on left and right
                                 node.diff.in = nodes$cn - colSums(adj)
                                 node.diff.out = nodes$cn - rowSums(adj)

                                 norm.adj = as.data.table(cbind(1:length(segs), match(gr.end(segs), gr.start(segs))))[!is.na(V2), ]
                                 norm.adj = rbind(norm.adj, norm.adj[, .(V2 = -V1, V1 = -V2)])[, nid1 := match(V1, nodes$nid)][, nid2 := match(V2, nodes$nid)]

                                 ## now add non-aberrant edge copy numbers that are the minimum of the unaccounted
                                 ## for copy number going <out> of the source node and going <in> to the sink node
                                 adj.old = adj
                                 ## ALERT: extremely hacky solution
                                 adj[as.matrix(norm.adj[, .(nid1, nid2)])] =
                                     pmax(pmin(node.diff.out[norm.adj[, nid1]],
                                          node.diff.in[norm.adj[, nid2]]), 0)

                                 nodes$eslack.in = nodes$cn - colSums(adj)
                                 nodes$eslack.out = nodes$cn - rowSums(adj)


                                 if (sum(adj!=0)>0)
                                 {
                                     if (!identical(adj[which(adj>0)],
                                                    adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = nodes$rix[col], col = nodes$rix[row])])]))
                                     {
                                         stop('reciprocality violated')
                                     }
                                 }
                                 end(nodes) = end(nodes)-1

                                 return(list(nodes, as(adj, 'Matrix')))
                             }

                             .parseparens = function(str)
                             {
                                 cmd = gsub(',$', '',
                                            gsub(',\\)', ')',
                                                 gsub('\\)', '),',
                                                      gsub('\\(', 'list(',
                                                           gsub('([^\\(^\\[^\\]^\\)]+)', '"\\1",', perl = TRUE, gsub('\\]', ')', gsub('\\[', '\\(', str)))))))
                                 eval(parse(text = cmd))
                             }

                             sols = lapply(dir(dir(paste(cougar, 'solve',sep = '/'), full = TRUE)[1], '^g_', full = TRUE), readLines)
                             if (length(sols)==0){
                                 if (verbose){
                                     
                                 }
                                 return(self$nullGGraph())
                             }

                             ## parse cougar graphs
                             graphs = lapply(sols, .parsesol)

                             ## concatenate nodes and block diagonal bind adjacency matrices
                             segs = do.call('c', lapply(graphs, '[[', 1))
                             segs$id = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$id, sep = '.')
                             segs$nid = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$nid, sep = '.')
                             segs$ix = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$ix, sep = '.')
                             segs$rix = paste(rep(1:length(graphs), sapply(lapply(graphs, '[[', 1), length)), segs$rix, sep = '.')
                             segs$rix = match(segs$rix, segs$ix)
                             segs$ix = 1:length(segs)
                             adj = do.call('bdiag', lapply(graphs, '[[', 2))


                             ## final double check for identicality
                             if (!(identical(adj[which(adj>0)], adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = segs$rix[col], col = segs$rix[row])])])))
                             {
                                 stop('Reciprocality check failed!')
                             }

                             ## TODO: figure out why there are negative CN edges in CouGaR!!!
                             private$gGraphFromScratch(segs = segs, es = adj)$fillin() ## IT MAY BE NOT BALANCED
                             return(self)                         
                         },

                         remixt2gg = function(remixt){
                             if (!dir.exists(remixt)){
                                 stop("Input ReMixT directory not found.")
                             } else if (length(rmt.out <- dir(remixt, "cn.tsv$|brk.tsv$", full.names=TRUE)) != 2){
                                 stop("Required output files cn.tsv$ and brk.tsv$ cannot be located.")
                             }
                             rmt.seg = fread(grep("cn.tsv", rmt.out, value=TRUE))
                             rmt.seg[, ":="(start = shift(end)), by=chromosome]
                             rmt.seg[, start := start + 1]
                             rmt.seg[is.na(start), start:=1]
                             rmt.seg[, cn := major_1 + minor_1]
                             rmt.tile = dt2gr(rmt.seg)
                             rmt.bks = fread(grep("brk.tsv", rmt.out, value=TRUE))
                             if (nrow(rmt.bks)>0){
                                 strmap = setNames(c("+", "-"), c("-", "+"))
                                 rmt.bks[, cn := cn_1] ## only consider major clone right now
                                 bp1 = dt2gr(rmt.bks[, .(seqnames=chromosome_1, start=position_1, end=position_1, strand=strmap[strand_1])])
                                 bp2 = dt2gr(rmt.bks[, .(seqnames=chromosome_2, start=position_2, end=position_2, strand=strmap[strand_2])])
                                 juncs = grl.pivot(GRangesList(list(bp1, bp2)))
                                 values(juncs) = rmt.bks[, .(prediction_id, cn, cn_0, cn_1, cn_2, n_1, side_1, n_2, side_2)]
                             } else {
                                 juncs = NULL
                             }

                             rmt.tile$loose = FALSE
                             out = self$initialize(tile = rmt.tile,
                                                   juncs = juncs,
                                                   cn=TRUE)$simplify()$fillin()

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

                             gmessage('None')
                         } else if (nrow(private$es)==0){
                             gmessage('None')

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
                                         vertex.frame.color=NA, vertex.label.color = "black"
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
                         ## ## DONE
                         ## if (length(private$segs)==0){
                         ##     return(0L)
                         ## }
                         ## if (is.null(private$partition)){
                         ##     private$partition = self$components()
                         ## }
                         ## return(private$parition$no)
                         ## changing the definition
                         if (!"loose" %in% colnames(private$segs)){
                             return(length(private$segs %Q% (strand=="+")))
                         } else {
                             return(length(private$segs %Q% (loose==FALSE & strand=="+")))
                         }                         
                     },
                     ##

                     gg2td = function(seg.col, ...){
                         if (verbose <- getOption("gGnome.verbose")){
                             gmessage("Create gTrack for static genome browser-style viz.")
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

                         if ("orig" %in% colnames(values(ss))){
                             cn.map = setNames(ss[ss$orig]$cn, as.character(ss[ss$orig]$qid))
                             ss$count = table(gr.string(ss))[gr.string(ss)]
                             ss$cn = cn.map[as.character(ss$qid)]
                             ss$id = gr.match(ss, unique(ss))
                             ss = dt2gr(gr2dt(ss)[, cn := cn + seq(0, by=0.5, length.out = count[1]) - 0.5 * ((count[1]-1)/2),by=id])
                             gt = gTrack(ss, y.field="cn", edges=ed, name="CN", angle=0, gr.colorfield="base", ...)
                         } else if ("cn" %in% colnames(values(ss))){
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
                             gmessage("No gGnome.js repository found on your system.")
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
                             gmessage("Writing your JSON file to:", filename)
                         }

                         ## copy the whole directory
                         system(paste("cp",
                                      system.file("extdata", "gGnome.js", package="gGnome"),
                                      basedir
                                      ))

                         ## generating the JSON
                         if (verbose <- getOption("gGnome.verbose")){
                             gmessage("Writing your JSON file to:", filename)
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
                                      settings = list(y_axis = list(title = "copy number",
                                                                    visible = TRUE)),
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
                                 gmessage('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }
                         }

                         if (verbose <- getOption("gGnome.verbose")){
                             gmessage("Create json file for interactive visualization.")
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

                         hb = hydrogenBonds(private$segs, private$id.column)

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

                         ## Allow the code to work if there is no cn field
                         if(!(is.null(private$segs$cn) | any(is.na(private$segs$cn)))) {
                             node.dt[, y := private$segs$cn[oid]]
                         } else {
                             node.dt[, y := 1]
                         }

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

                             ## Will eclass break SNVs ?
                             
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
                                       weight = 1, #temporary FIXME
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
                                 gmessage("Saving JSON to: ", filename)
                             }
                             jsonlite::write_json(gg.js, filename,
                                                  pretty=TRUE, auto_unbox=TRUE, digits=4)
                             return(normalizePath(filename))
                         } else {
                             return(gg.js)
                         }
                     },

                     julie.copy.gg2js = function(filename='.',
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
                                 gmessage('Creating directory ', basedir)
                                 system(paste('mkdir -p', basedir))
                             }
                         }

                         if (verbose <- getOption("gGnome.verbose")){
                             gmessage("Create json file for interactive visualization.")
                         }


                         
                         ## Returns all of the nodes in the graph as a data.table
                         dt = function() {
                             return(as.data.table(private$pnodes))
                         },

                         

                         ## Returns all the edges in the graph as a data.table
                         edgesdt = function() {
                             return(copy(private$pedges))
                         }                       
                     ),
                     
                     )

