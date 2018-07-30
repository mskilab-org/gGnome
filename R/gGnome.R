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
library(skitools)
library(gUtils)
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
                            browser()
                            i = with(self$dt, eval(i)) ## allows subsetting based on metadata

                            if (is.logical(i)){
                                i = which(i)
                                }

                            if (is.numeric(i) | is.integer(i)) {
                                if (any(i<0) | max(i, na.rm = TRUE)>self$length()) {
                                    stop('index out of bounds')
                                }
                                
                            private$pnode.id = private$pnode.id[i]
                            private$pindex = private$pindex[i]
                            private$prindex = private$prindex[i]
                            private$porientation = private$porientation[i]

                            }
                            else
                            {
                                i=as.numeric(i)
                                dt=gr2dt(private$pgraph$gr)
                                private$pnode.id = private$pnode.id[dt[snode.id==i, index]]
                                private$pindex = private$pindex[dt[snode.id==i, index]]
                                private$prindex = private$prindex[dt[snode.id==i, index]]
                                private$porientation=private$porientation[dt[snode.id==i, index]]
                            }
                            
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
                            new.dt=copy(self$dt)
                            new.dt[, start:=self$graph$dt[from==start, snode.id]]
                            new.dt[, end:=self$graph$dt[to==end, snode.id]]
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
                         
                         get.adj = function(flat=FALSE){
                             if (is.null(private$pgraph)){
                                 self$get.g()
                             }
                             adjMat = igraph::as_adj(private$pgraph)
                             if (flat) {
                                 return(adjMat)
                             } else {
                                 if (is.element("cn", colnames(private$pedges))){
                                     adjMat[private$pedges[,cbind(from, to)]] = private$pedges$cn
                                 }
                                 return(adjMat)
                             }
                         },
                         
                          get.g = function(force=FALSE){
                           if (!is.null(private$pgraph) & !force){
                               return(private$pgraph)
                           } else {
                               private$pgraph = igraph::make_directed_graph(
                                   t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes))
                               return(private$pgraph)
                           }
                       },
                         
                         simplify = function(mod=FALSE)
                       {
                           ## if two or more segment are only connected by ref edges
                           ## and they have the same copy number
                           ## merge them into one node
                           if (verbose <- getOption("gGnome.verbose")){
                               message("Merge all pairs of nodes only connected by reference edge.")
                           }

                           if (length(private$pnodes)==0 || is.null(private$pedges) || nrow(private$pedges) == 0) {
                               if (verbose){
                                   warning("Empty graph/no edges. Cannot be simpler.")
                               }
                               return(self)
                           }


                           ## Get the non-loose nodes and assign them their index
                           ## We can only consider positive strand because ref edges go from pos-pos or neg-neg
                           node.dt = gr2dt(private$pnodes)
                           node.dt[, id := seq_along(private$pnodes)]

                           ## Get the number of edges into and out of each node
                           node.dt[, ":="(in.d = private$pedges[, table(to)[as.character(1:nrow(node.dt))]],
                                          out.d = private$pedges[, table(from)[as.character(1:nrow(node.dt))]])]
                           node.dt[is.na(in.d), in.d:=0]
                           node.dt[is.na(out.d), out.d := 0]

                           ## Get the edges that are reference and get the ones with only 1 edge on either side
                           ref.es = private$pedges[type=="reference"]
                           ref.es[, ":="(to.in.d = node.dt[to, in.d],
                                         from.out.d = node.dt[from, out.d])]
                           ref.es = ref.es[to.in.d==1 | from.out.d==1]

                           ## Get the number of ref edges going into and out of each node, if there are any merge them
                           node.dt[, ":="(ref.in.d = ref.es[, table(to)[as.character(1:nrow(node.dt))]],
                                          ref.out.d = ref.es[, table(from)[as.character(1:nrow(node.dt))]])]
                           node.dt[is.na(ref.in.d), ref.in.d := 0]
                           node.dt[is.na(ref.out.d), ref.out.d := 0]
                           

                           ## THIS IS PROBABLY THE PROBLEM
                           tos = node.dt[(in.d == 1 & ref.in.d == 1), id]
                           froms = node.dt[(out.d == 1 & ref.out.d == 1), id]

                           ref.es = ref.es[to %in% tos & from %in% froms]

                           ## now with candidate ref edges, merge the two nodes one by one
                           nix.to.simplify = ref.es[, sort(unique(c(from, to)))]

                           ref.A = Matrix::sparseMatrix(i = ref.es[, from],
                                                        j = ref.es[, to],
                                                        x = 1,
                                                        dims=dim(self$get.adj()))
                           ref.G = igraph::graph_from_adjacency_matrix(ref.A)

                           ## find all weak components here
                           ref.comps = igraph::components(ref.G, "weak")

                           ## within each cluster, identify the start and end of the Hamilton path
                           node.dt[, map := ref.comps$membership]

                           ## Here we want to mark nodes that weren't changed
                           node.dt[, left.bound := !id %in% nix.to.simplify[ nix.to.simplify %in% ref.es[,to] ]]
                           node.dt[, right.bound := !id %in% nix.to.simplify[ nix.to.simplify %in% ref.es[,from] ]]

                           ## check point: every weakly connected component should have one head and one tail!!!!!
                           if (!node.dt[, .(sum(left.bound)==1, sum(right.bound)==1), by = map][, all(V1) & all(V2)]){
                               stop("You found a bug in gg$simplify!")
                           }

                           ## start merging
                           segs = private$pnodes
                           segs$map = node.dt[,map]
                           new.segs = gr.simplify(segs, field="map")
                           new.tail = node.dt[left.bound==TRUE, setNames(map, id)]
                           new.head = node.dt[right.bound==TRUE, setNames(map, id)]

                           if ("cn" %in% colnames(private$pedges)){
                               new.es = private$pedges[, .(from = new.tail[as.character(from)],
                                                           to = new.head[as.character(to)],
                                                           type, cn)] ## assume there is always "type" in es
                           } else {
                               new.es = private$pedges[, .(from = new.tail[as.character(from)],
                                                           to = new.head[as.character(to)],
                                                           type)] ## assume there is always "type" in es
                           }
                           new.es = new.es[!is.na(from) & !is.na(to)]

                           ## Set up nodes and edges for the constructor                           
                           new.es=private$pairNodesAndEdges(new.segs, new.es)[[2]]
                           new.es = private$convertEdges(new.segs, new.es)
                           new.segs = new.segs %Q% (loose == FALSE & strand == "+")

                           if (mod){
                               private$reset()
                               private$gGraphFromNodes(new.segs, new.es)
                               return(self)
                           } else {
                               out = gGraph$new(nodes = new.segs, edges = new.es)
                               return(out)
                           }
                       },

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
                           if (length(setdiff(streduce(private$pnodes), win))==0){
                               return(self)
                           }

                           ## overlapping window and segs, removing loose ends
                           interGr = gr.findoverlaps(private$pnodes, win, ignore.strand=ignore.strand)
                           lid = which(private$pnodes$loose==T)
                           interGr = interGr %Q% (!query.id %in% lid)
                           qix = interGr$query.id

                           if (is.null(k)){
                               ## no k, use distance
                               if (is.null(d) | d < 0){
                                   stop("Must provide either valid k or d.")
                               }

                               ## blend window with segs
                               win = gr.fix(win, private$pnodes)## fix seqinfo
                               ss = tryCatch(c(private$pnodes[private$pnodes$loose == F, c()],
                                               win[, c()]), error = function(e) NULL)

                               if (is.null(ss)){
                                   ss = grbind(c(private$pnodes[private$pnodes$loose == FALSE, c()],
                                                 win[, c()]))
                               }

                               if (ignore.strand){
                                   ss = gr.stripstrand(ss)
                               }

                               ## TODO: can we take overlapping segs into the picture??
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
                               browser()
                               return(self$trim(hoodRange, mod=mod))
                           }
                           else {
                               ## with k, go no more k steps
                               kNeighbors = unique(unlist(igraph::ego(private$pgraph, order=k, nodes=qix)))
                               return(self$subgraph(v = kNeighbors, mod=mod)) ## not garanteed size to scale
                           }
                       },

                        dist = function(gr1,
                                       gr2 = NULL,
                                       matrix=T,
                                       EPS=1e-9,
                                       include.internal=TRUE, ## consider bp within feature "close"
                                       directed=FALSE, ## if TRUE, only consider gr1-->gr2 paths
                                       verbose=FALSE)
                       {
                           if (verbose <- getOption("gGnome.verbose")){
                               message("Given two GRanges, return pairwise shortest path distance.")
                           }
                           
                           if (is.null(gr2)){
                               gr2 = gr1 ## self distance
                           }
                           ## DONE
                           if (verbose)
                               now = Sys.time()

                           ##Make simpler
                           simpleNodes = self$simplify()

                           intersect.ix = gr.findoverlaps(gr1, gr2, ignore.strand = FALSE)

                           ngr1 = length(gr1)
                           ngr2 = length(gr2)

                           tiles = simpleNodes$nodes$gr

                           G = self$get.g()

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
                               gr1 = gr1[, 'id'] %**% tiles
                               ## BUG here, why does this overlap command give NA seqnames output?
                               gr2 = gr2[, 'id'] %**% tiles
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
                                   igraph::shortest.paths(G, uix1, uix2, weights = E(G)$weight,
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
                                       epathw = sapply(epaths,
                                                       function(x,w) if (length(x)==0) { Inf } else { sum(w[x]) }, E(G)$weight) ## calculate the path weights
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
                        trim = function(gr=NULL, mod=FALSE){
                           "Return a trimmed subgraph that all the nodes are wihtin the range defined by gr."
                           if (verbose <- getOption("gGnome.verbose")){
                               message("Given a GRanges, return the trimmed subgraph overlapping it.")
                           }

                           ## Catch cases of empty self and trim and return self
                           ## Return copy of self if mod is FALSE
                           if ((length(private$pnodes)==0 || is.null(gr) || length(gr)==0) & mod) {
                               return(self)
                           } else if (length(private$pnodes)==0 || is.null(gr) || length(gr)==0) {
                             return(self$clone())
                           }

                           ## FIXME: Should have a catch here for if we don't need to trim

                           es = private$convertEdges(private$pnodes, private$pedges, metacols=T)

                           gr = gr.fix(gr, private$pnodes)
                           gr = streduce(gr)

                           ## Get positive overlaps
                           segs = private$pnodes %Q% (loose == FALSE & strand == "+")
                           new.segs = gr.findoverlaps(segs, gr)

                           segs = segs[new.segs$query.id]

                           ## Find which nodes had their ends changed
                           ## Edges only will be kept if the end of a node was not trimmed
                           new.segs$left = start(new.segs)==start(segs)
                           new.segs$right = end(new.segs)==end(segs)
                           mcols(new.segs) = cbind(mcols(new.segs), mcols(segs)) ## carry over the metadata

                           ## map the edges
                           if (nrow(private$pedges)==0){
                               new.es = private$pedges
                           }
                           else {
                               ## Get edges that have n1 or n2 == query.id
                               validNodes = new.segs$query.id
                               new.es = es[n1 %in% validNodes & n2 %in% validNodes]

                               ## Want to remove edges that went into nodes that were trimmed on one side
                               ## Find query.id's that have all FALSE in left or right - might still be multiple trims but one side is trimmed
                               sg = gr2dt(new.segs)[, .(mean(left)>0, mean(right)>0), by=query.id]
                               leftRemove = sg[V1 == FALSE, query.id]
                               rightRemove = sg[V2 == FALSE, query.id]

                               ## Remove edges that had one of their sides trimmed
                               new.es = new.es[!(n1 %in% leftRemove & n1.side == 0) & !(n2 %in% leftRemove & n2.side == 0)]
                               new.es = new.es[!(n1 %in% rightRemove & n1.side == 1) & !(n2 %in% rightRemove & n2.side == 1)]

                               ## Now we have to remap the edges
                               ## map left==TRUE to n1 or n2 side == 0
                               ## map right==TRUE to n1 or n2 side == 1
                               map = data.table(old = c(new.segs[new.segs$left]$query.id, new.segs[new.segs$right]$query.id),
                                                side = c(rep(0, length(which(new.segs$left))), rep(1, length(which(new.segs$right)))),
                                              new = c(which(new.segs$left), which(new.segs$right)))
                               setkeyv(map, c("old", "side"))
                               new.es[, ":="(n1 = map[.(n1,n1.side), new], n2 = map[.(n2,n2.side), new])]
                           }

                           ## Remove miscellaneous metatcols added in this function
                           new.segs$left = NULL
                           new.segs$right = NULL
                           new.segs$query.id = NULL
                           new.segs$subject.id = NULL
                           new.es[, sedge.id:=NULL]
                           ## finally, recreate the trimmed graph
                           if(mod) {
                               private$gGraphFromNodes(nodes = new.segs,
                                                       edges = new.es,
                                                       looseterm=T)
                           } else {
                               newSg = gGraph$new(nodes = new.segs,
                                                  edges = new.es,
                                                  looseterm=T)
                               ## FIXME: This should be overriden in bGraph then
                               ## if (!newSg$isBalance() & "cn" %in% colnames(newSg$edges)){
                               ##     ## NOTE: why do I have to reassign it here??
                               ##     newSg = newSg$fillin(mod=TRUE)
                               ## }
                               ## if (newSg$isBalance()){
                               ##     if (inherits(self, "bGraph")){
                               ##         newSg = as(newSg, "bGraph")
                               ##     }
                               ## }
                               return(newSg)
                           }
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
                             browser()
                             if (node==TRUE){                                 
                                 gr.dt=gr2dt(private$pnodes)
                                 gr.dt[, paste(col1):=vec]
                                 private$pnodes=makeGRangesFromDataFrame(gr.dt, TRUE)
                             }
                             else{                                
                                 private$pedges[, col1:=vec]
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
                             segs$index = 1:length(segs)
                             snodes=list()
                             segdt=gr2dt(segs)                             
                             x=1                            
                             for (i in 1:length(segs)){                                                                  
                                if (any(segdt[index<i, id]==segdt[index==i, id])){
                                    old.ix=which(segdt[index<i, id]==segdt[index==i, id])
                                    snodes[i]=snodes[[old.ix]]*-1                           
                                }
                                else {
                                    if (segdt[index==i, strand]=="+"){
                                        snodes[i]=x
                                        x=x+1
                                    }
                                    else{
                                        snodes[i]=-1*x
                                        x=x+1
                                    }
                                }
                             }
                             
                             snodes=unlist(snodes)
                             segs$snode.id=snodes
                             adj = do.call('bdiag', lapply(graphs, '[[', 2))
                             
                             
                             ## final double check for identicality
                             if (!(identical(adj[which(adj>0)], adj[as.matrix(as.data.table(which(adj!=0, arr.ind = TRUE))[, .(row = segs$rix[col], col = segs$rix[row])])])))
                             {
                                 stop('Reciprocality check failed!')
                             }
                             
                             ## TODO: figure out why there are negative CN edges in CouGaR!!!                            
                             summ=summary(adj)
                             adj=data.table(from= summ$i, to = summ$j, cn = summ$x)                             
                             segs$loose=FALSE
                             private$pairNodesAndEdges(segs, adj)
                             private$convertEdges(segs, adj)
                             private$gGraphFromNodes(nodes = segs, edges = adj) ## IT MAY BE NOT BALANCED
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

