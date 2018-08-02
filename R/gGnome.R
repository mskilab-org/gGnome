##REMOVE THESE EVENTUALLY
library(data.table)
library(Matrix)
library(jsonlite)
library(gUtils)
library(igraph)
library(gTrack)
library(GenomicRanges)
library(skitools)
library(gUtils)

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



## ================= gNode class definition ================== ##
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

                        
                      mark = function(...)
                      {
                        args = list(...)
                        for (nm in names(args))
                          private$pgraph$annotate(nm, args[[nm]],
                                                  private$pnode.id, "node")
                      },
                        
                        ## Flips the current orientation of the Node object (swap index/rindex)
                        flip = function()
                        {
                            tmp = private$index
                            private$pindex = private$prindex
                            private$prindex = tmp
                            private$porientation = -private$porientation
                            return(self)
                        },


                        ## Allows subseting of the Node object using bracket notation
                        subset = function(i)
                        {                      
                            if (is.logical(i)){
                                i = which(i)
                            }

                            if (is.numeric(i) | is.integer(i)) {
                                if (max(abs(i), na.rm = TRUE)>self$length()) {
                                  stop('index out of bounds')                              
                                }
                            }

                            ## if character is provided as the query, we will
                            ## use signed node id to look up signed node.ids in this object
                            ## *** we will find the first match of the absolute value of the queried
                            ## snode id and flip that node.id if necessary
                            if (is.character(i)) 
                            {
                              i = as.integer(i)
                              i = sign(i)*match(abs(i), abs(private$pnode.id))
                              
                              if (any(is.na(i)))
                                stop('one or more signed node ids are not found in this object')
                            }

                            private$pnode.id = private$pnode.id[abs(i)]
                            new.pindex = ifelse(i>0, private$pindex[abs(i)], private$prindex[abs(i)])
                            new.prindex = ifelse(i>0, private$prindex[abs(i)], private$pindex[abs(i)])
                            private$pindex = new.pindex
                            private$prindex = new.prindex
                            private$porientation = private$porientation[abs(i)]
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

                      ## returns all the edges associated with this node
                      edges = function()
                      {
                        return(union(self$eleft, self$eright))
                      },

                        ## Returns the edges connected to the left of the nodes in this Node object as an Edge Object
                        eleft = function()
                        {
                            sedge.id = private$pgraph$edgesdt[to %in% private$pindex, sedge.id]
                            return(gEdge$new(sedge.id, graph = private$pgraph))
                        },
                        

                        ## Returns the edges connected to the right of the nodes in this Node object as an Edge Object
                        eright = function()
                        {
                            sedge.id = private$pgraph$edgesdt[from %in% private$pindex, sedge.id]
                            return(gEdge$new(sedge.id, graph = private$pgraph))
                        },


                        lleft = function()
                        {
                            return(private$pgraph$loose[private$pgraph$eloose[to %in% private$pindex, from] - length(private$pgraph$gr)])
                        },


                        lright = function()
                        {
                            return(private$pgraph$loose[private$pgraph$eloose[from %in% private$pindex, to] - length(private$pgraph$gr)])
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


## ================== Non-Member Functions for gNode ================== ##

#' @name length
#' @title length
#' @description
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


#' @name c
#' @title c.gNode
#' @description
#' Concatenates gEdge objects by id's
#'
#' @param gEdge objects
#'
#' @return a new concatenated gEdge Object
#' @export
`c.gNode` = function(...)
{                            
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
}


#' @name setdiff
#' @title setdiff.gNode
#' @description
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
#' @title union
#' @description
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
#' @title intersect
#' @description
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
#' @return A new gNode object that contains only the given id's
#' @export
#' 
'[.gNode' = function(obj, i = NULL){
  nodes = obj$clone()
  inew = tryCatch(with(nodes$dt, eval(parse(text = lazyeval::expr_text(i)))), error = function(e) NULL)
  if (!is.null(inew))
    i = inew
  done = nodes$subset(i)
  return(nodes)
}


## ================= gEdge class definition ================== ##
#' @export
gEdge = setClass("gEdge")

gEdge = R6::R6Class("gEdge",
                    public = list(
                        ## Constructor GEDGE
                        ## Builds an gEdge Class from a given edge.id 
                        ## graph - the gGraph you want to build this class from 
                        initialize = function(seid = NULL, graph)
                        {
                            private$pgraph = graph ## reference to graph that these edges belong to
                            private$porientation = private$pedge.id = private$psedge.id = c()                       
                            
                            if (is.null(seid)) {
                                return(self)
                            }
                            
                            if (any(is.na(private$pgraph$edgesdt[.(seid), edge.id]))) {
                                stop("one or more provided signed edge ids are out of bounds")
                            }
                            
                            private$porientation = rep(1, length(seid))
                            private$psedge.id = seid ## signed edge id, referring to the edge in the directed graph for which the $left direction is 5' and $right direction is '
                            private$pedge.id = abs(seid) ## unsigned edge id
                            private$pedges = private$pgraph$edgesdt[list(seid), ]
                            return(self)
                        },
                        
                        mark = function(...)
                        {
                          args = list(...)
                          for (nm in names(args))
                            private$pgraph$annotate(nm, args[[nm]],
                                                    c(private$psedge.id, -private$psedge.id), "edge")

                        },
                        
                        ## Returns the number of edge pairs in this class
                        length = function()
                        {                            
                            return(length(private$pedge.id))
                        },
                        
                        ## Allows for subsetting of the gEdge Object using bracket notation
                        subset = function(i)
                        {

                            if (is.logical(i))
                                i = which(i)

                            if (is.numeric(i) | is.integer(i))
                            {
                                if (max(i, na.rm = TRUE)>self$length())
                                    stop('index out of bounds')
                            }

                            if (is.character(i))
                            {
                              i = as.integer(i)
                              i = sign(i)*match(abs(i), abs(private$psedge.id))
                            }

                            private$psedge.id = sign(i)*private$psedge.id[abs(i)]
                            private$pedges = private$pgraph$edgesdt[.(private$psedge.id), ]
                            private$pedge.id = private$pedges$edge.id
                            private$porientation = sign(i)*private$porientation[abs(i)]

                            return(self)
                        },
                        
                        
                        ## Prints out the number of edges and the count of each type
                        print = function()
                        {
                          message("gEdge object with ", self$length(), " edges")
                          sides = c('left', 'right')
                          print(cbind(convertEdges(private$pgraph$gr, self$dt, metacols = TRUE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
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
                            leftNodes = private$pedges[, from]
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

                      nodes = function()
                      {                      
                        return(union(self$left, self$right))
                      },

                      ## Returns a data.table edges in this class format (to, from, type, edge.id)
                      dt = function()
                      {
                        return(copy(private$pedges))
                      },
                      
                      grl = function()
                      {
                        gr1 = gr.flipstrand(gr.end(private$pgraph$gr[private$pedges$from], ignore.strand = FALSE))
                        gr2 = gr.start(private$pgraph$gr[private$pedges$to], ignore.strand = FALSE)
                        grl = split(c(gr1, gr2), rep(1:length(gr1), 2))[as.character(1:length(gr1))]
                        names(grl) = private$edges$sedge.id
                        meta = cbind(private$pedges, data.table(bp1 = gr.string(gr1), bp2 = gr.string(gr2)))
                        values(grl) = meta[, unique(c("edge.id", "sedge.id", "from", "to", "bp1", "bp2", colnames(meta))), with = FALSE]
                        return(grl)
                      },
                      
                      ## Returns a Junctions Object representation of each of the edges with metadata populated
                      junctions = function()
                      {                        
                        return(Junction$new(self$grl))
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


## ================== Non-Member Functions for gEdge ================== ##

#' @name [.gEdge
#' @title gEdge
#' @description
#'
#' Overloads subset operator for edges
#'
#' @param obj gEdge object This is the gEdge object to be subset
#' @param i integer, logical, or expression in gEdge metadata used to subset gEdges
#' @return A new gEdge object that contains only the given id's
#' @export
'[.gEdge' = function(obj, i = NULL){
  edges = obj$clone()
  inew = tryCatch(with(edges$dt, eval(parse(text = lazyeval::expr_text(i)))), error = function(e) NULL)
  if (!is.null(inew))
    i = inew
  edges$subset(i)
  return(edges)
}


#' @name length
#' @title length
#' @description 
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


#' @name c
#' @title c.gNode
#' @description
#' Concatenates gNode objects by id's
#'
#' @param gNode objects
#'
#' @return a new concatenated gNode Object
#' @export
`c.gEdge` = function(...)
{                            
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
}


#' @name setdiff
#' @title setdiff
#' @description
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
#' @title union
#' @description
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
#' @title intersect
#' @description
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

#' @name union
#' @title union
#' @description
#' Returns a new Junction object which is the union of x and y.
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#'
#' @return new Junction Object containing the union of x and y
setMethod("union", c("Junction", "Junction"),
          function(x, y) {                            
              newJunc=c(x, y)
              newJunc$removeDups()
              return(newJunc)
          })

#' @name setdiff
#' @title setdiff
#' @description
#' Returns a new Junction object which is the difference between x and y (id's).
#'
#' @param x a Junction Object
#' @param y a Junction Object
#'
#' @return new Junction containing the difference between x and y
setMethod("setdiff", c("Junction", "Junction"),
          function(x, y) {
              ## Make sure that both come from the same graph                           
              juncs1=x             
              juncs2=y                          
              difs = lapply(1:length(x), function(x){                  
                  z=TRUE                  
                  for (i in 1:length(juncs2$juncs)){
                      if (identical(juncs1$juncs[[x]], juncs2$juncs[[i]]))
                      {
                          z=FALSE
                      }
                  }
                  if(z==TRUE){
                      return(juncs1$juncs[[x]])
                  }
              })                   
              difs=GRangesList(plyr::compact(difs))              
              return(Junction$new(difs))
          })

#' @name intersect
#' @title intersect
#' @description
#' Returns a new gNode object which is the intersection of x and y (id's).
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#'
#' @return new Junction Object containing the intersection of x and y
setMethod("intersect", c("Junction", "Junction"),
          function(x, y) {              
              diff=c(setdiff(x, y), setdiff(y, x))
              all=c(x, y)
              all$removeDups()
              intersect=setdiff(all, diff)
              return(intersect)
          })


## ================== Junction class definition ================== ##
#' @export
Junction = setClass("Junction")

Junction = R6::R6Class("Junction",
                       public = list( ##JUNCTION
                         ## Builds junction class, grl must be GRangesList with each GRanges of length 2
                         ## Empty junctions are removed
                         ## If grl is empty GRangesList, Junctions class is empty

                           initialize = function(grl)
                           {
                               ## Check to make sure input is GRangesList with each element of length 2
                               ## Will allow elements to be empty and will just remove them
                               if (!inherits(grl, "GRangesList")) {
                                   stop("Input is not a GRangesList")
                               }
                               
                               widths = lengths(grl)
                               empty = widths == 0
                               
                               ix = widths != 2
                               if(any(ix[!empty])) {
                                   stop(paste0("grl contains junctions that are of improper length. Indices are: ", paste(ix[!empty], collapse=" ")))
                               }

                               grl.unlisted=unlist(grl)
                               if (any(width(grl.unlisted)!=1)){
                                   stop(paste0("grl contains GRanges that are of improper length."))
                               }
                               
                               private$pjuncs = grl[!empty]

                               return(self)
                           },

                           intersect=function(junc1, junc2){
                               
                           },
                           ## Allows subseting of the Junction object using bracket notation
                           subset = function(i)
                           {                               
                               if (is.logical(i))
                                   i = which(i)
                               
                               if (is.numeric(i) | is.integer(i)) {
                                   if (max(i, na.rm = TRUE)>self$length()) {
                                       stop('index out of bounds')
                                   }

                                   if (any(i<0))
                                     {
                                       if (!all(i<0))
                                         stop('cannot mix positive with negative subscripts for Junction object')

                                       i = setdiff(1:self$length(), abs(i))
                                     }
                               }
                               
                               private$pjuncs = private$pjuncs[i]
                               
                               return(self)
                           },
                           
                           
                           removeDups = function()
                           {                               
                               list.gr=lapply(1:length(private$pjuncs), function(x){                                   
                                   return(private$pjuncs[[x]])
                               })                               
                               private$pjuncs=GRangesList(unique(list.gr))
                               return(self)
                           },
                           
                           ## Prints the Junctions
                           print = function()
                           {
                               message("Junction Object with ", self$length(), " junctions\n")
                               print(self$dt)
                           },


                           ## Returns the number of junctions
                           length = function()
                           {                               
                               return(length(private$pjuncs))
                           }
                       ),
                       
                       private = list(
                           ## GRangesList with 2 pairs per
                           pjuncs = NULL
                       ),

                       active = list(
                           ## Returns a GRangesList of the junctions
                           grl = function()
                           {
                               return(private$pjuncs)
                           },

                           dt = function()
                           {
                               return(as.data.table(values(private$pjuncs)))
                           },
                           
                           ## Returns a GRanges representing the spots a genome needs to break for these junctions to be connected
                           ## Will remove duplicate breakpoints
                           breakpoints = function()
                           {
                               bps = granges(unlist(private$pjuncs))
                               names(bps) = NULL
                               
                               ## Deal with shift at value one
                               gr = bps %Q% (strand == "+")
                               gr = GenomicRanges::shift(gr, -1)
                               
                               gr = c(bps %Q% (strand == "-"), gr)
                               gr = gr[!duplicated(gr)]
                               
                               return(gr)
                           },
                           
                           ## Returns a gGraph created from this junctions Object
                           graph = function()
                           {
                               return(gGraph$new(juncs = self))
                           }
                       )
                       )


## ================== Non-Member Functions for Junction ================== ##

#' @name length
#' @title length
#' @description
#' The number of junctions in the Junction Object
#'
#' @param Junction a Junction object
#'
#' @return the number of junctions in the Junction Object
#' 
#' @export
`length.Junction` = function(Junction)
{
    if(!inherits(Junction, "Junction")) {
        stop("Error: Invalid input")
    }
    return(Junction$length())
}


#' @name c
#' @title c
#' @description
#' Concatenates Junction objects
#'
#' @param Junction object
#'
#' @return a new concatenated Junction Object
#' @export
`c.Junction` = function(...)
{                            
    juncs.list=list(...)
    isg = sapply(juncs.list, function(x) class(x)[1]=='Junction')
    
    if(any(!isg)){
        stop('Error: All inputs must be of class Junction.')
    }

    ## Get all the pjuncs to create new Junction Object
    grls = lapply(juncs.list, function(x) x$grl)
    return (Junction$new(BiocGenerics::Reduce(c, grls)))
}


#' @name [.Junction
#' @title Junction
#' @description
#'
#' Overloads subset operator for junctions
#'
#' @param obj Junction object This is the Junction object to be subset
#' @param i integer, logical, or expression in Junction metadata used to subset Junction
#' #' @return A new Junction object that contains only the given id's
#' @export
'[.Junction' = function(obj, i = NULL){  
  juncs = obj$clone()
  inew = tryCatch(with(juncs$dt, eval(parse(text = lazyeval::expr_text(i)))), error = function(e) NULL)
  if (!is.null(inew))
    i = inew
  juncs$subset(i)
  return(juncs)
}


## ================== gGraph class definition ================== ##
gGraph = setClass("gGraph")

gGraph = R6::R6Class("gGraph",
                     public = list(


                         ## public fields
                         ## name = NULL,
                         ## refG = "GENOME", ## seqinfo of ref genome
                         ## constructor INIT GGRAPH
                         initialize = function(genome = NULL,
                                               tile = NULL,
                                               juncs = NULL,
                                               prego = NULL,
                                               jabba = NULL,
                                               cougar = NULL,
                                               weaver = NULL,
                                               remixt = NULL,
                                               nodes = NULL,
                                               edges = NULL,
                                               nodeObj = NULL,
                                               edgeObj = NULL)
                         {
                             ## control how to construct
                             verbose = getOption("gGnome.verbose")
                             if (!is.null(nodes)){
                                 private$gGraphFromNodes(nodes, edges)
                             }
                             else if (!is.null(nodeObj))
                             {
                                 private$gGraphFromNodeObj(nodeObj, edgeObj)
                             }
                             else if(!is.null(tile) || !is.null(juncs))
                             {
                                 private$karyograph(tile, juncs)
                             }
                             else if(!is.null(prego))
                             {
                                 private$pr2gg(prego)
                             }
                             else if (!is.null(jabba))
                             {
                                 private$jab2gg(jabba)
                             }
                             else if (!is.null(cougar))
                             {
                                 private$cougar2gg(cougar)
                             }
                             else if (!is.null(weaver))
                             {
                                 private$wv2gg(weaver)
                             }
                             else if (!is.null(remixt))
                             {
                                 private$remixt2gg(remixt)
                             }
                             else
                             {
                                 private$emptyGGraph(genome)
                             }
                             
                             return(self)
                         },
                         
                         ## snode.id is a vector of signed node.id's from our gGraph
                         ## Returns a data.table of the index and reverse index in the graph for this snode.id
                         queryLookup = function(id) {
                             dt = private$lookup[.(id)]
                             return(dt)
                         },

                         
                         ## Returns the number of nodes in the graph
                         length = function()
                         {
                             return(self$nodes$length())
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

                             ## Assign them their index
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
                             new.es = private$pairNodesAndEdges(new.segs, new.es)[[2]]
                             new.es = convertEdges(new.segs, new.es)
                             new.segs = new.segs %Q% (strand == "+")

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
                             lid = 1:length(private$pnodes)
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

                         
                         print = function()
                         {
                             cat("gGraph with ", self$length(), " nodes and ", nrow(private$pedges)/2, " edges")
                             if (self$nodes$length() > 0)
                                 message(', containing:\n')
                             else
                                 message('\n')
                             self$nodes$print()
                             
                             if (nrow(private$pedges))
                             {
                                 message()
                                 self$edges$print()
                             }
                         },


                         ## id - unsigned node.id/edge.id, etc.
                         ## class - character containing the name of the class
                         annotate = function(colName, data, id, class)
                         {
                             if (class == "node") {
                                 
                                 id = self$queryLookup(id)
                                 index = c(id[, index], id[, rindex])
                                 
                                 gr.dt = gr2dt(private$pnodes)
                                 gr.dt[index, paste(colName) := rep(data, 2)]
                                 values(private$pnodes)[[colName]] = gr.dt[[colName]]
                                 
                             } else if (class == "edge") {
                                 private$pedges[.(id), paste(colName) := data]
                                 
                             } else {
                                 stop("Not sure how we got to this error at all, we should never be here")
                             }

                             return (self)
                         },

                         
                         window = function(pad = 0)
                         {
                             return(streduce(gr.stripstrand(private$pnodes + pad)))
                         },
                         
                       
                       gtrack = function(y.field = 'cn', name = 'gGraph', labels.suppress = TRUE, stack.gap = 1e5, ...)
                       {
                         ss = private$pnodes
                         ed = private$pedges

                         if (is.null(ss$loose))
                           ss$loose = FALSE

                         if (!is.null(ed))
                         {
                           
                           ## set edge apperances
                           ## lwd, lty, col, cex.arrow, v, not.flat, h, dangle.w
                           if (!is.element(y.field, colnames(ed))) {
                             ed[, y := 1]
                           } else
                           {                               
                             ed$y = ed[[y.field]]
                           }

                           if (!("col" %in% names(ed)))
                           {
                             ed$col = as.character(NA)
                           }

                           if (!("lty" %in% names(ed)))
                           {
                             ed$lty = as.numeric(NA)
                           }
                           
                           if (!("v" %in% names(ed)))
                           {
                             ed$v = as.numeric(NA)
                           }

                           if (!("h" %in% names(ed)))
                           {
                             ed$h = as.numeric(NA)
                           }

                           if (!("cex.arrow" %in% names(ed)))
                           {
                             ed$cex.arrow = as.numeric(NA)
                           }

                           if (!("lwd" %in% names(ed)))
                           {
                             ed$lwd = as.numeric(NA)
                           }

                           if (!("dangle.w" %in% names(ed)))
                           {
                             ed$dangle.w = as.numeric(NA)
                           }

                           if (!("not.flat" %in% names(ed)))
                           {
                             ed$not.flat = as.logical(NA)
                           }

                           ed$col = as.character(ed$col)
                           ed$border = as.character(ed$border)
                           
                           ed[, ":="(lwd = ifelse(is.na(lwd), ifelse(type=="aberrant", log2(0.2*y+2)+1, 1),lwd),
                                     lty = ifelse(is.na(lty), ifelse(type=='loose', 3, 1),lty),
                                     col = ifelse(is.na(col), ifelse(is.na(col), ifelse(type=="aberrant", 
                                                              ifelse(y>0,
                                                                     alpha("red", 0.4),
                                                                     alpha("purple", 0.3)),
                                                              ifelse(type=="loose",
                                                                     alpha("blue",0.6),
                                                                     alpha("grey",0.5))), col), col),
                                     cex.arrow = ifelse(is.na(cex.arrow), 0, cex.arrow),
                                     not.flat = ifelse(is.na(not.flat), type=="aberrant", not.flat),
                                     v = ifelse(is.na(v), ifelse(type=="aberrant", 2, 1),v),
                                     h = ifelse(is.na(h), ifelse(type=="aberrant", 2, 1),h),
                                     dangle.w = ifelse(is.na(dangle.w), 0.5, dangle.w))]

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
                           ss$y[lid] = ss$y[pid]*1.2
                         }
                         
                         ## col, border, ywid
                         if (!("col" %in% names(values(ss))))
                           ss$col = as.character(NA)
                         
                         ## col, border, ywid
                         if (!("border" %in% names(values(ss))))
                           ss$border = as.character(NA)

                         ss$border = as.character(ss$border)
                         ss$col = as.character(ss$col)

                         if (!("lwd.border" %in% names(values(ss))))
                           ss$lwd.border = as.numeric(NA)
                         
                         if (!("ywid" %in% names(values(ss))))
                           ss$ywid = as.numeric(NA)

                         ss$col = ifelse(is.na(ss$col), ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5)), ss$col)
                         ss$border = ifelse(is.na(ss$border), ifelse(ss$loose, ss$col, alpha("black", 0.5)), ss$border)
                         ss$lwd.border = ifelse(is.na(ss$lwd.border), 1, ss$lwd.border)
                         ss$ywid = ifelse(is.na(ss$ywid), ifelse(ss$loose, 0.001, 0.8), ss$ywid)

                         if (y.field %in% colnames(values(ss))){
                           ss$y = values(ss)[, y.field]
                           gt = gTrack::gTrack(unname(ss), y.field=y.field, edges=ed, labels.suppress = labels.suppress,
                                               name=name, angle=0, ...)
                         } else {
                           ## stack node pairs via stack.gap
                           tmp.ss = ss[ss$snode.id>0]
                           tmp.ss$y = disjointBins(tmp.ss+stack.gap)
                           ss$y = tmp.ss$y[match(ss$node.id, tmp.ss$node.id)]
                           gt = gTrack::gTrack(ss, y.field = 'y', yaxis = FALSE, edges=ed, labels.suppress = labels.suppress,
                                               name=name, angle=0, ...)
                         }
                         return(gt)
                       },
                       
                         
                         ## Trim
                         ## Returns a new graph trimmed around a provided GRanges
                         ## tile - granges to trim the graph to
                         trim = function(tile, mod=F)
                         {
                             ## Some quick catch cases
                             if(length(tile) == 0) {
                                 stop("tile cannot contain no nodes, nothing to trim around")
                             }
                             if(length(self) == 0) {
                                 return(self)
                             }

                             es = convertEdges(private$pnodes, private$pedges, metacols=T)
                             
                             tile = gr.fix(tile, private$pnodes)
                             tile = streduce(tile)
                             
                             ## Get positive overlaps
                             nodes = unname(private$pnodes) %Q% (strand == "+")
                             new.nodes = gr.findoverlaps(nodes, tile)

                             ## If the new trim has no nodes, return empty graph
                             if(length(new.nodes) == 0) {
                                 return(gGraph$new())
                             }
                             
                             nodes = nodes[new.nodes$query.id]
                             
                             ## Find which nodes had their ends changed
                             ## Edges only will be kept if the end of a node was not trimmed
                             new.nodes$left = start(new.nodes)==start(nodes)
                             new.nodes$right = end(new.nodes)==end(nodes)
                             mcols(new.nodes) = cbind(mcols(new.nodes), mcols(nodes)) ## carry over the metadata
                             
                             ## map the edges
                             if (nrow(private$pedges)==0){
                                 edges = private$pedges
                             }
                             else {
                                 ## Get edges that have n1 or n2 == query.id
                                 validNodes = new.nodes$query.id
                                 new.es = es[n1 %in% validNodes & n2 %in% validNodes]
                                 
                                 ## Want to remove edges that went into nodes that were trimmed on one side
                                 ## Find query.id's that have all FALSE in left or right - might still be multiple trims but one side is trimmed
                                 sg = gr2dt(new.nodes)[, .(mean(left)>0, mean(right)>0), by=query.id]
                                 leftRemove = sg[V1 == FALSE, query.id]
                                 rightRemove = sg[V2 == FALSE, query.id]
                                 
                                 ## Remove edges that had one of their sides trimmed
                                 new.es = new.es[!(n1 %in% leftRemove & n1.side == 0) & !(n2 %in% leftRemove & n2.side == 0)]
                                 new.es = new.es[!(n1 %in% rightRemove & n1.side == 1) & !(n2 %in% rightRemove & n2.side == 1)]
                                 
                                 ## Now we have to remap the edges
                                 ## map left==TRUE to n1 or n2 side == 0
                                 ## map right==TRUE to n1 or n2 side == 1
                                 map = data.table(old = c(new.nodes[new.nodes$left]$query.id, new.nodes[new.nodes$right]$query.id),
                                                  side = c(rep(0, length(which(new.nodes$left))), rep(1, length(which(new.nodes$right)))),
                                                  new = c(which(new.nodes$left), which(new.nodes$right)))
                                 setkeyv(map, c("old", "side"))
                                 new.es[, ":="(n1 = map[.(n1,n1.side), new], n2 = map[.(n2,n2.side), new])]
                             }
                             
                             ## Remove miscellaneous metatcols added in this function
                             new.nodes$left = NULL
                             new.nodes$right = NULL
                             new.nodes$query.id = NULL
                             new.nodes$subject.id = NULL

                             if(mod) {
                                 private$gGraphFromNodes(nodes = new.nodes,
                                                         edges = new.es)
                                 return(self)
                             } else {
                                 return(gGraph$new(nodes = new.nodes,
                                                   edges = new.es))
                             }
                         },


                         ## @name mergeOverlaps
                         ## @brief Takes all overlapping nodes and collapses them by breaking the
                         ##        nodes below them and merging.
                         ## @param mod TRUE if we want to modify this graph by reference
                         ## @return Decoupled gGraph object, but only if mod=FALSE
                         mergeOverlaps = function(mod=TRUE)
                         {
                             ## ASSUMPTION: nodes of a gGraph are always skew-symmetric
                             ## Check to make sure the graph has overlaps
                             if (isDisjoint(unname(private$pnodes) %Q% (strand=="+"))){
                                 print("Overlaps already merged")
                                 return(self)
                             }

                             ##FIXME: have extra edge copies because aberrant edges can be all of the edge types
                             
                             ## If we only have aberrant edges, don't use the nodes in tile
                             nodes = disjoin(unname(private$pnodes) %Q% (strand=="+"))
                             
                             ## Turn the edges that aren't reference into junctions, constructor will
                             ## remake the reference edges
                             all.j = self$edges[type == "aberrant"]$junctions ##FIXME: want to do here $removeDuplicates()
                             
                             if (mod==T) {
                                 private$karyograph(tile = nodes, juncs = all.j, genome = seqinfo(nodes))
                                 return(self)
                             } else {
                                 out = gGraph$new(tile = nodes, juncs = all.j)
                                 return(out)
                             }
                         },
                         
                         
                         ## Takes a set of nodes and adds them to our graph. Changes this graph by reference. May provide
                         ## optional table of edges to connect these nodes (in the form n1,n2,n1.side,n2.side).
                         ## Treats gr as strandless
                         addNodes = function(gr, edges = NULL, mod = TRUE)
                         {
                             if(length(gr) == 0) {
                                 return(self)
                             }
                             
                             gg = gGraph$new(nodes = gr, edges = edges)
                             self$mergeGraphs(gg = gg, mod = mod)

                             return(self)
                         },


                         ## @name mergeGraphs
                         ## @brief Merges two gGraphs together into a new gGraph. By default, this function creates
                         ##        a new gGraph and merges overlapping nodes via decouple.
                         ## @param mod If TRUE, sets this gGraph equal to the merged result.
                         ## @param decouple If TRUE, combines overlapping nodes using decouple().
                         ## @return Merged gGraph object, but only returns if mod is FALSE (this is the default)
                         mergeGraphs = function(gg, mod = FALSE)
                         {
                             ## Make sure that gg is a gGraph
                             if (!inherits(gg, "gGraph"))
                                 stop("Error: Can only deal with addition of two gGraph objects.")

                             ## If gg is a null gGraph or has no nodes, don't do anything
                             if (gg$length() == 0) {
                                 return (self)
                             }

                             ## If self is a null gGraph or has no nodes, return gg
                             if (self$length() == 0) {
                                 return (gg)
                             }

                             ## bare GRanges
                             new.segs = c(granges(self$gr), granges(gg$gr))

                             ## Get metadata columns that both graphs have and bind them
                             common.segs.mc = intersect(colnames(values(self$gr)),colnames(values(gg$gr)))

                             if (length(common.segs.mc)>0){
                                 values(new.segs) = rbind(values(self$gr)[, common.segs.mc],
                                                          values(gg$gr)[, common.segs.mc])
                             }

                             new.es = convertEdges(self$gr, self$edgesdt)
                             gg.edges = convertEdges(gg$gr, gg$edgesdt)

                             common.es.mc = setdiff(intersect(colnames(new.es), colnames(gg.edges)), c("from", "to"))
                             
                             ## If both edge tables are non-empty, bind them together
                             ## If this graph has empty edge table, swap new.es with gg.edges
                             if (dim(gg.edges)[1] != 0 && dim(new.es)[1] != 0) {
                                 new.es = rbind(new.es[,.(n1,n2,n1.side,n2.side,
                                                          cn = if('cn' %in% common.es.mc) cn,
                                                          type = if('type' %in% common.es.mc) type)],
                                                gg.edges[,.(n1 = n1 + nrow(self$dt[strand == "+"]),
                                                            n2 = n2 + nrow(self$dt[strand == "+"]),
                                                            n1.side, n2.side,
                                                            cn = if('cn' %in% common.es.mc) cn,
                                                            type = if('type' %in% common.es.mc) type)]
                                                )
                             } else if (dim(gg.edges)[1] != 0) {
                                 new.es = gg.edges
                             }

                             names(new.segs) = NULL
                             new.segs = new.segs %Q% (strand == "+")

                             ##browser()
                             
                             if (mod){
                                 private$gGraphFromNodes(nodes = new.segs,
                                                         edges = new.es)
                                 return(self)
                             } else {
                                 return(gGraph$new(nodes = new.segs, edges = new.es))
                             }
                         },


                         ## adds a GRangesList of junctions or a Junction Object to this object via breakpoints
                         ## If the desired junctions to add do not overlap with our graph, do not add them
                         ## juncs - GRangesList() or Junction Object of junctions to add
                         ## mod - TRUE if we to change this graph
                         addJuncs = function(juncs, mod = TRUE)
                         {
                             ## Potential problems: What if some junctions aren't within the range of our current graph

                             ## If juncs is not Junction Object, try to convert it
                             if(!is.null(juncs) && !inherits(juncs, "Junction")) {
                                 juncs = tryCatch(Junction$new(juncs),
                                                  error = function(e) {
                                                      NULL
                                                  })
                                 if (is.null(juncs)) {
                                     stop("Input is not a valid junctions set.")
                                 }
                             }
                             
                           nodes = private$pnodes
                           edges = convertEdges(private$pnodes, private$pedges)
                           names(nodes) = NULL
                           nodes = unname(nodes) %Q% (strand == "+")


                             ## Behavior: juncs not contained in gGraph are removed, juncs with half in the graph
                             ## will break at those locations but will not add an edge from the junctions

                             ## Want to break out genome at all the junction locations
                             ## Want to check our junctions and remove from the unlisted one all the pairs of junctions where either of them is not present in the overlaps
                             ## Then add the junctions to the edge table
                             
                             ## Break our current genome
                             bps = juncs$breakpoints
                             strand(bps) = "*"
                             nodes = gr.breaks(bps, nodes)
                             
                             juncsGRL = unlist(juncs$grl)

                             ## Remove junctions that aren't needed
                             tmp = juncsGRL
                             strand(tmp) = "+"
                             foundJuncs = findOverlaps(tmp,nodes)@from
                             
                             fullJuncs = which(table(ceiling(foundJuncs/2)) > 1)
                             juncsGRL = juncsGRL[c(fullJuncs, 2*fullJuncs)]
                             
                             ## Find which node each junction overlaps
                             tmp = juncsGRL
                             strand(tmp) = "+"
                             index = findOverlaps(tmp, nodes)@to
                             
                             startEdge = juncsGRL[seq(1, length(index), by=2)]
                             endEdge = juncsGRL[seq(2, length(index), by=2)]
                             
                             ## FIXME: still have this labeling edges problem happening here/cn problem
                             new.edges = data.table(n1 = index[seq(1, length(index), by=2)], n2 = index[seq(2, length(index), by=2)],
                                                    n1.side = ifelse(as.logical(strand(startEdge) == "-"), 1, 0),
                                                    n2.side = ifelse(as.logical(strand(endEdge) == "-"), 1, 0),
                                                    type = "aberrant",
                                                    cn = if ("cn" %in% names(edges)) NA)
                             
                             edges = rbind(edges, new.edges)
                             
                             if (mod) {
                                 private$gGraphFromNodes(nodes = nodes, edges = edges)
                                 return(self)
                             } else {
                                 return(gGraph$new(nodes = nodes,
                                                   edges = edges))
                             }
                         },
                         
                         
                         ## Creates a json file for active visualization using gGnome.js
                         json = function(filename='.',
                                         maxcn=100,
                                         maxweight=100,
                                         save = TRUE,
                                         settings = list(y_axis = list(title = "copy number",
                                                                       visible = TRUE)),
                                         no.y = FALSE)
                         {
                             ## Make sure the our nodes are not empty before visualizing
                             if (is.null(private$pnodes) || length(private$pnodes) == 0) {
                                 stop("Cannot have empty graph and visualize it")
                             }

                             if (save){
                                 if (grepl('\\.js(on)*$', filename)){
                                     ## if json path was provided
                                     basedir = dirname(filename)
                                 }
                                 else if (filename==".") {
                                     ## default path was provided
                                     basedir = './'
                                     filename = "data.json"
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
                             regsegs.ix = which(as.character(seqnames(private$pnodes))
                                                %in% names(regular.sl))

                             ## Add loose nodes/edges
                             
                             loose.ix = which(private$pnodes$loose==TRUE)

                             ed = copy(private$pedges) ## otherwise change by reference!
                             ## construct intervals
                             node.dt = data.table(oid = which(as.logical(strand(private$pnodes)=="+")))
                             ## node.dt[, rid := seq_along(private$pnodes)[-oid][match(private$pnodes[-oid],
                             ##                                                      gUtils::gr.flipstrand(
                             ##                                                                  private$pnodes[oid]
                             ##                                                              ))]]

                             private$buildLookupTable()
                             lookup = copy(private$lookup)
                             setkey(lookup, index)

                             ##browser()
                             
                             node.dt[, rid := lookup[.(oid), rindex]]

                             node.dt = node.dt[oid %in%
                                               which(private$pnodes$loose==FALSE &
                                                     as.character(seqnames(private$pnodes))
                                                     %in% names(regular.sl))]

                             node.dt[, iid := 1:.N]
                             setkey(node.dt, "iid")
                             node.dt[, ":="(chr = as.character(seqnames(private$pnodes[oid])),
                                            start = start(private$pnodes[oid]),
                                            end = end(private$pnodes[oid])
                                            ),
                                     ]
                             node.map = node.dt[, c(setNames(iid, oid),
                                                    setNames(iid, rid))]

                             ## Allow the code to work if there is no cn field
                             if (!is.null(private$pnodes$cn)) {
                                 node.dt[, y := private$pnodes$cn[oid]]
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
                                           reid = paste(lookup[.(to), rindex],
                                                        lookup[.(from), rindex],
                                                        sep="-"))]

                                 ## ALERT: bc strandlessness, I only retained half of the edges
                                 ## to map edges in gwalks, we will need strandedness,
                                 ## so will retain everything
                                 ed[,":="(soStr = as.character(strand(private$pnodes[from])),
                                          siStr = as.character(strand(private$pnodes[to])))]

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
                                           weight = if("cn" %in% names(ed)) { cn } else { 1 },
                                           cid = eclass)]

                                 ed[, ":="(source = so*so.str,
                                           sink = -si*si.str)]
                                 
                                 ed[, weight := ifelse(is.na(weight), 1, weight)]
                                 
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


                         ## vectorized function, connects with an edge index1 -> index2 and then creates the other edge (snode.id? possibly?)
                         ## FIXME: Need to find a way to label edges
                         connectNodes = function(snode.id1, snode.id2)
                         {
                             if(!all(c(snode.id1,snode.id2) %in% private$pnodes$snode.id)) {
                                 stop("Some snode.id's are not present in this graph")
                             }
                             if(length(snode.id1) != length(snode.id2)) {
                                 stop("snode.id1 and snode.id2 must have the same length")
                             }
                             if(length(snode.id1) == 0 || length(snode.id2) == 0) {
                                 stop("snode.id1 or snode.id2 are empty")
                             }

                             ## Build the lookup table again incase it wasn't reset
                             private$buildLookupTable()
                             
                             index1.dt = self$queryLookup(snode.id1)
                             index2.dt = self$queryLookup(snode.id2)

                             edges = copy(private$pedges)
                             new.edges = rbind(data.table(to = index2.dt[, index],
                                                          from = index1.dt[, index],
                                                          edge.id = max(edges[, edge.id]) + 1:nrow(index1.dt),
                                                          sedge.id = max(edges[, edge.id]) + 1:nrow(index1.dt),
                                                          type = "aberrant",
                                                          cn = if("cn" %in% names(edges)) 1),
                                               data.table(to = index1.dt[, rindex],
                                                          from = index2.dt[, rindex],
                                                          edge.id = max(edges[, edge.id]) + 1:nrow(index1.dt),
                                                          sedge.id = -1*(max(edges[, edge.id]) + 1:nrow(index1.dt)),
                                                          type = "aberrant",
                                                          cn = if("cn" %in% names(edges)) 1))
                             
                             ## FIXME: Possible metadata problem
                             edges = rbind(edges, new.edges)
                             
                             private$pedges = edges
                             private$pgraph = igraph::make_directed_graph(
                                 t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes))
                             
                             return(self)
                         }
                     ),
                     

                     private = list( #### PRIVATE GGRAPH
                         ## ----- private fields
                         ## ===== required
                         ## node/vertex, a GRanges obj of strand-specific ranges
                         pnodes = NULL,
                         ## data.table of all edges in g, from, to, cn, type
                         ## type can be ref, aberrant
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
                                                    edges = NULL)
                         {
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
                                 
                                 ## Make sure to remove these columns as we will add them again later
                                 if ("sedge.id" %in% names(edges)) {
                                     edges[, sedge.id := NULL]
                                 }
                                 if("edge.id" %in% names(edges)) {
                                     edges[, edge.id := NULL]
                                 }
                             }

                             if (any(is.na(seqlengths(nodes))))
                                 nodes = gUtils::gr.fix(nodes)
                             
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
                                 ## FIXME: Need to change this to make sure it labels edges
                                 if(!"type" %in% names(edges)) {
                                     edges[, type := 'aberrant']
                                 }
                                 
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
                             
                             private$buildLookupTable()

                             return(self)
                         },

                         ## Adds the proper snode.id and index to nodes, adds proper sedge.id/edge.id to edges
                         ## Returns list(nodes, edges) with updated fields and miscellaneous metadata removed
                         ## USE THIS FUNCTION WITH CAUTION - gGraphFromNodes will require that sedge.id/edge.id is removed before running (FIXME: don't make it do that)
                         pairNodesAndEdges = function(nodes, edges)
                         {                            
                             if('loose' %in% names(values(nodes))) {
                                 not.loose = which(!nodes$loose)
                             } else {
                                 not.loose = seq_along(nodes)
                             }

                             edges = as.data.table(edges)

                             ix = not.loose[match(gr.stripstrand(nodes[not.loose]), gr.stripstrand(nodes[not.loose]))]
                             nodes$snode.id = NA
                             nodes$snode.id[not.loose] = ifelse(as.logical(strand(nodes[not.loose]) == '+'), ix, -ix)
                             nodes$index = 1:length(nodes)
                             
                             edges[, tag := paste(nodes$snode.id[to], nodes$snode.id[from])]
                             edges[, rtag := paste(-nodes$snode.id[from], -nodes$snode.id[to])]
                             edges[, id := 1:.N]
                             edges[, rid := match(tag, rtag)]                                 
                             edges[, edge.id := igraph::clusters(graph.edgelist(cbind(id, rid)), 'weak')$membership]
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

                             ## If we have no edges, just call the constructor on the Nodes Object
                             if(is.null(EdgeObj) || EdgeObj$length() == 0) {
                                 private$gGraphFromNodes(nodes = NodeObj$gr)
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
                           nodes = nodes[!duplicated(nodes$index), ] ## MARCIN FIX .. happens with subgraphs from node / edge neibhbors

                             edges = convertEdges(nodes, edges)
                             nodes = NodeObj$gr
                             
                             private$gGraphFromNodes(nodes = nodes, edges = edges)

                             return(self)
                         },


                         ## Builds a gGraph by breaking the reference genome at the points specified in tile
                         ## Treats tile as nothing but breakpoints, no metadata, nothing special
                         ## Juncs can be either a GRangesList of junctions in proper format or a Junction Object
                         ## If tile has length 0, this function creates a simple graph (1 node per chromosome, each one full length)
                         ## If genome is specified, it will try to use that genome instead - if it is invalid it wil use the default
                         karyograph = function(tile = NULL,
                                               juncs = NULL,
                                               genome = NULL)
                         {
                             ## Make sure user entered some input
                             if (is.null(tile) & is.null(juncs)) {
                                 stop("Cannot have both tile and juncs be NULL, must be some input")
                             }

                             ## If juncs is not Junction Object, try to convert it
                             if(!is.null(juncs) && !inherits(juncs, "Junction")) {
                                 juncs = tryCatch(Junction$new(juncs),
                                                  error = function(e) {
                                                      NULL
                                                  })
                                 if (is.null(juncs)) {
                                     stop("Input is not a valid junctions set.")
                                 }
                             }

                             ## Validate tile as GRanges
                             if (!is.null(tile) && !inherits(tile, "GRanges")){
                                 tile = tryCatch(GRanges(tile),
                                                 error=function(e){
                                                     NULL
                                                 })
                                 if (is.null(tile)){
                                     stop("Input cannot be converted into a GRanges object.")
                                 }
                             }

                             ## If the user provided a genome, check if it is valid. If it isn't, set to NULL
                             if (!is.null(genome)) {
                                 tmp = tryCatch(si2gr(genome), error=function(e) NULL)
                                 if (is.null(tmp)) {
                                     warning("Input 'genome' cannot be converted to a seqinfo. Set to NULL.")
                                     genome = NULL
                                 }
                             }
                             
                             ## Build a GRanges from the default genome
                             ## FIXME: if using misc chr, definitely need to specify genome                           
                             sl = gUtils::hg_seqlengths(genome=genome)

                           nodes = si2gr(sl)
                           edges = data.table(n1 = numeric(0), n2 = numeric(0), n1.side = numeric(0), n2.side = numeric(0))
                             
                             ## Break the genome based on whether there is tile or juncs
                             if (!is.null(tile) && length(tile) > 0) {
                                 nodes = gUtils::gr.breaks(tile, nodes)
                             }
                             if (!is.null(juncs) && length(juncs) > 0) {
                                 bps = juncs$breakpoints
                                 strand(bps) = "*"
                                 nodes = gUtils::gr.breaks(bps, nodes)
                             }

                             ##browser()
                             
                             ## If we broke from junctions, add their edges
                             if (!is.null(juncs) && length(juncs) > 0) {
                               juncsGRL = gr.fix(unlist(juncs$grl), nodes)
                               nodes = gr.fix(nodes, juncsGRL)

                                 ## Remove junctions that aren't in the genome selected
                                 tmp = juncsGRL
                                 strand(tmp) = "+"
                                 foundJuncs = findOverlaps(tmp,nodes)@from
                                 fullJuncs = which(table(ceiling(foundJuncs/2)) > 1)
                                 juncsGRL = juncsGRL[c(fullJuncs, 2*fullJuncs)]

                                 ## Find which node each junction overlaps
                                 tmp = juncsGRL
                                 strand(tmp) = "+"
                                 index = findOverlaps(tmp, nodes)@to
                                 
                                 startEdge = juncsGRL[seq(1, length(index), by=2)]
                                 endEdge = juncsGRL[seq(2, length(index), by=2)]

                                 ## Map the Junctions to an edge table
                                 ## Mapping is as follows in terms of junctions a -> b
                                 ##         a- => leaves/enters left side of base a
                                 ##         a+ => leaves/enters right side of base a
                                 ##  a- <-> a+ => leaves left side of base a, enters right side of base a (NOT THE BASE NEXT TO a)
                                 new.edges = data.table(n1 = index[seq(1, length(index), by=2)], n2 = index[seq(2, length(index), by=2)],
                                                        n1.side = ifelse(as.logical(strand(startEdge) == "-"), 1, 0),
                                                        n2.side = ifelse(as.logical(strand(endEdge) == "-"), 1, 0))

                                 edges = rbind(edges, new.edges)
                             }
                             
                             ## If we broke from tile, add their edges
                             if (!is.null(tile) && length(tile) > 0) {
                                 ## BIG ASSUMPTION HERE THAT ALL TILES WERE INVOLVED IN BREAKING
                                 strand(tile) = "*"

                                 ## Find starting tile points that are past the first base
                                 ## FIXME: might have problem with last base or off the seqnames bases etc.
                                 ## NOTE: gr.breaks treats single base GRanges as endpoints so dont count them in starts
                                 startTile = tile[lengths(tile) > 1]

                                 ## Gets the end of every chromosome that corresponds to the tile, lets us make sure no ref edges are marked when the tile reaches the end of the chromosome
                                 endCut = seqlengths(tile)[as.character(rep(seqnames(tile)@values, seqnames(tile)@lengths))]
                                 
                                 endIndex = findOverlaps(gr.start(startTile[start(startTile) > 1]), nodes)@to
                                 startIndex = findOverlaps(gr.end(tile[end(tile) > endCut]), nodes)@to

                                 new.edges = data.table(n1 = c(endIndex-1, startIndex),
                                                        n2 = c(endIndex, startIndex+1),
                                                        n1.side = rep(1, length(endIndex) + length(startIndex)),
                                                        n2.side = rep(0, length(endIndex) + length(startIndex)))
                                 
                                 edges = rbind(edges, new.edges)
                             }
                             
                             nodes$qid = NULL
                             names(nodes) = NULL
                             
                             private$gGraphFromNodes(nodes, edges)
                             
                             return(self)
                         },
                         
                         
                         ## Generates a gGraph from a prego output file
                         pr2gg = function(fn)
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
                                                right.tag = res[[1]]$node2)
                             
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
                             
                             edges = private$pairNodesAndEdges(nodes, ed)[[2]]
                             edges = convertEdges(nodes, edges)
                             private$gGraphFromNodes(nodes = granges(posNodes), edges = edges)
                             
                             return(self)
                         },
                         
                         
                         ## @name jab2gg
                         ## @brief Constructor for creating a gGraph object from a jabba output
                         jab2gg = function(jabba, regular=FALSE)
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
                                 edges = convertEdges(nodes, edges)
                             }
                            
                             ## We want to get only the positive strand of nodes
                             nodes = unname(nodes) %Q% (strand == "+")
                             
                             ## Need to get seqinfo, if all of the seqlenths are missing, fill them in
                             if (any(is.na(seqlengths(nodes)))) {
                                 
                                 ## Try to use junction seqlenths
                                 if (!any(is.na(seqlengths(jabba$junctions)))){
                                     
                                     nodes = gUtils::gr.fix(nodes, jabba$junctions)

                                 } else {
                                     
                                     default.sl = data.table::fread(Sys.getenv("DEFAULT_BSGENOME"))[, setNames(V2, V1)]

                                     nodes = gUtils::gr.fix(nodes, default.sl)
                                 }
                             }
                             
                             ## Set up gGraph using other constructor, standardizes our input although not totally necessary
                             private$gGraphFromNodes(nodes, edges)
                             
                             ## If regular, remove non-regular nodes from nodes by trimming
                             if (regular==T) {

                                 if (getOption("gGnome.verbose")){
                                     warning("Forcing regular chromosomes. Will try default. See `Sys.getenv('DEFAULT_REGULAR_CHR')`.")
                                 }
                                 regularChr = si2gr(data.table::fread(Sys.getenv('DEFAULT_REGULAR_CHR'))[, setNames(V2, V1)])
                                 self$trim(regularChr, mod=T)
                             }
                             
                             return(self)
                         },


                         wv2gg = function(weaver)
                         {
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

                             ## FIXME: When I update karyograph for cn
                             private$karyograph(tile = ss, juncs = junc)
                             
                             ## ALERT: Weaver out put is not balanced!
                             
                             return(self)
                         },

                         ## Initialize a gGraph from RemixT output
                         ## FIXME: Need to simplify the output here but that has no meaning right now
                         remixt2gg = function(remixt)
                         {
                             if (!dir.exists(remixt)){
                                 stop("Input ReMixT directory not found.")
                             } else if (length(rmt.out <- dir(remixt, "cn.tsv$|brk.tsv$", full.names=TRUE)) != 2){
                                 stop("Required output files cn.tsv$ and brk.tsv$ cannot be located.")
                             }
                             rmt.seg = fread(grep("cn.tsv", rmt.out, value=TRUE))
                             rmt.seg[, ":="(start = data.table::shift(end)+1)]
                             rmt.seg[, start := ifelse(start > end, 1, start)]
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
                             private$karyograph(tile = rmt.tile, juncs = juncs) ## FIXME: Need to incorporate cn here
                             return(self)
                         }
                         
                     ),

                     active = list(
                         ## Returns a Node Object of the nodes in the graph
                         nodes = function()
                         {
                           if (length(self$gr)>0)
                             tmp = private$pnodes$snode.id[private$lookup[.(1:(length(private$pnodes)/2)), index]]
                           else
                             tmp = NULL
                           return(gNode$new(tmp, self))
                         },

                         ## Returns an Edge Object of the edges in the graph
                         edges = function()
                         {
                            if (nrow(self$edgesdt)>0)
                             tmp = private$pedges[.(1:(nrow(private$pedges)/2)), sedge.id]
                           else
                             tmp = NULL
                            return(gEdge$new(tmp, self))
                         },

                                                                       ## Returns an Edge Object of the edges in the graph
                         junctions = function()
                         {
                             return(self$edges$junctions)
                         },

                         ## Returns all the loose nodes in the graph - dynamically calculated
                         ## Only essential metadata - snode.id, node.id, index will be kept on these nodes
                       
                         loose = function()
                         {
                             nodes = unname(private$pnodes) %Q% (strand == "+")
                             
                             loose.left = setdiff(nodes$index, private$pedges[, to])
                             loose.right = setdiff(nodes$index, private$pedges[, from])
                             loose.nodes = granges(c(gr.start(nodes[loose.left]), gr.end(nodes[loose.right])))
                             
                             loose.nodes$node.id = (max(nodes$node.id) + 1):(max(nodes$node.id) + length(loose.nodes))
                             loose.nodes = c(loose.nodes, gr.flipstrand(loose.nodes))
                             loose.nodes$snode.id = ifelse(as.logical(strand(loose.nodes) == "+"), loose.nodes$node.id, -loose.nodes$node.id)
                             loose.nodes$index = (length(private$pnodes) + 1):(length(private$pnodes) + length(loose.nodes))
                                                                             
                             return (loose.nodes)
                         },

                       #' returns adjacency matrix of graph as sparseMatrix 
                       adj = function(){

                         adjMat = igraph::as_adj(self$G)

                         if (is.element("cn", colnames(private$pedges))){
                           adjMat[private$pedges[,cbind(from, to)]] = private$pedges$cn
                         }
                                 return(adjMat)
                       },
                                              
                       G = function(force=FALSE){
                         if (!is.null(private$pgraph) & !force){
                           return(private$pgraph)
                         } else {
                           return(igraph::make_directed_graph(
                                                      t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes)))
                         }
                       },
                       
                         ## Returns a data.table of the edges needed to append if $loose was run and appended to private$pnodes
                         eloose = function()
                         {
                             nodes = unname(private$pnodes) %Q% (strand == "+")

                             loose.left = setdiff(nodes$index, private$pedges[, to])
                             loose.right = setdiff(nodes$index, private$pedges[, from])
                             
                             eleft = data.table(from = (length(private$pnodes)+1):(length(private$pnodes) + length(loose.left)),
                                                to = loose.left,
                                                type = "loose",
                                                edge.id = (max(private$pedges[, edge.id]) + 1):(max(private$pedges[, edge.id]) + length(loose.left)))
                             eright = data.table(from = loose.right,
                                                 to = (max(eleft[, from]) + 1):(max(eleft[, from]) + length(loose.right)),
                                                 type = "loose",
                                                 edge.id = (max(eleft[, edge.id]) + 1):(max(eleft[, edge.id]) + length(loose.right)))
                             
                             releft = data.table(from = self$queryLookup(nodes[loose.left]$snode.id)[, rindex],
                                                 to = eleft[, from] + length(loose.left) + length(loose.right),
                                                 type = "loose",
                                                 edge.id = (max(eright[, edge.id]) + 1):(max(eright[, edge.id]) + length(loose.left)))
                             
                             reright = data.table(from = eright[, to] + length(loose.left) + length(loose.right),
                                                  to = self$queryLookup(nodes[loose.right]$snode.id)[, rindex],
                                                  type = "loose",
                                                  edge.id = (max(releft[, edge.id]) + 1):(max(releft[, edge.id]) + length(loose.right)))
                             
                             
                             loose.edges = rbindlist(list(eleft[, snode.id := edge.id], eright[, snode.id := edge.id],
                                                          releft[, snode.id := -edge.id], reright[, snode.id := -edge.id]))
                             
                             return (loose.edges)
                         },
                         
                         
                         ## Returns the igraph associated with this graph
                         graph = function()
                         {
                             if(!is.null(private$pgraph)) {
                                 return(private$pgraph)
                             } else {
                                 private$pgraph = igraph::make_directed_graph(
                                     t(as.matrix(private$pedges[,.(from,to)])), n=length(private$pnodes))
                                 return(private$pgraph)
                             }
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
                         },

                         
                         ## Returns a gTrack
                         gt = function()
                         {
                             self$gtrack(y.field = 'cn', stack.gap = 1e5)
                         },


                         ## Returns the footprint of this graph
                         footprint = function()
                         {
                             return(self$window())
                         }
                     )
)



#' @name [.gGraph
#' @title [.gGraph
#' @description
#'
#' Overloads subset operator for gGraph
#'
#' @param obj gWalk object this is the gGraph
#' @param i integer, logical, or expression in gNode metadata used to subset gNodes
#' @param j integer, logical, or expression in gEdge metadata used to subset gEdges
#' @return A new gGraph object that only contains the given nodes and edges
#' @export
'[.gGraph' = function(obj, i = NULL, j = NULL){

  if (deparse(substitute(j)) != "NULL")
  {
    edges = obj$edges[j]
    if (deparse(substitute(i)) == "NULL")
      nodes = edges$nodes
    else
      nodes = union(edges$nodes, obj$nodes[i])
  } else
  {
    if (deparse(substitute(i)) == "NULL")
      nodes = obj$nodes
    else
      nodes = obj$nodes[i]

    edges = nodes$edges
  }

  return(gGraph$new(nodeObj = nodes, edgeObj = edges))
}


## ================== Non-Member Functions for gGraph ================== ##

#' @name length
#' @title length
#' @description
#' The number of nodes in the gGraph
#'
#' @param gGraph a \code{gGraph} object
#'
#' @return the number of strongly connected components in this graph
#' @export
`length.gGraph` = function(gGraph){
    ## input must be a gGraph!
    if (!inherits(gGraph, "gGraph")){
        stop("Error: Invalid input.")
    }
    return(gGraph$length())
}


#' @name dim
#' @title dim
#' @description
#' The number of nodes and edges in the gGraph
#'
#' @param gGraph a \code{gGraph} object
#'
#' @return the number of nodes and edges in the gGraph
#' @export
`dim.gGraph` = function(gGraph){
    ## input must be a gGraph!
    if (!inherits(gGraph, "gGraph")){
        stop("Error: Invalid input.")
    }
    return(c(gGraph$length(), nrow(gGraph$edgesdt)/2))
}



#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param gGraph a gGraph object
#'
#' @return the seqinfo of this graph
setMethod("seqinfo", c("gGraph"),
          function(x) {
              return(x$seqinfo)
          })


#' @name %+%
#' @title gGraphPlus
#' @description
#' Adding two \code{gGraph} instances
#'
#' @param gg1, gg2 instances of the \code{gGraph} class
#'
#' @return a copy of a new \code{gGraph} object that is the simple sum of the two inputs
setGeneric("%+%", function(gg, x) standardGeneric("%+%"))

#' Adds together a gGraph and a gGraph
setMethod("%+%", c("gGraph","gGraph"),
          function(gg, x) {
              return(gg$mergeGraphs(x))
          })

#' Adds together a gGraph and gNode Object (helpful when gNode does not point to the same graph)
setMethod("%+%", c("gGraph", "gNode"),
          function(gg, x) {
              gg1 = x$subgraph
              return(gg$mergeGraphs(gg1))
          })

#' Adds together a gGraph and gEdge Object (helpful when gEdge does not point to the same graph)
setMethod("%+%", c("gGraph", "gEdge"),
          function(gg, x) {
              gg1 = x$subgraph
              return(gg$mergeGraphs(gg1))
          })

#' Adds together a gGraph and Junction Object - breaks gGraph at junction points and adds edges
setMethod("%+%", c("gGraph", "Junction"),
          function(gg, x) {
              return(gg$addJuncs(x, mod=F))
          })



## #' @name %Q%
## #'
## #' @param gg a gGraph Object to query
## #' @param query a String to query in the following format (node_queries %
## setGeneric("%Q%", function(gg, ...) standardGeneric("%Q%"))
## setMethod("%Q%", signature(gg = "gGraph"),
##           function(gg, y) {
##               condition_call  = substitute(y)
##               browser()
##               ## serious R voodoo gymnastics .. but I think finally hacked it to remove ghosts
##               ## create environment that combines the calling env with the granges env
##               env = as(c(as.list(parent.frame(2)), as.list(as.data.frame(gg$gr))), 'environment')
##               parent.env(env) = parent.frame()
##               ix = tryCatch(eval(condition_call, env), error = function(e) NULL)
##               if (is.null(ix))
##               {
##                   condition_call  = substitute(y)
##                   ix = eval(condition_call, GenomicRanges::as.data.frame(x))
##               }
##               return(x[ix])
##           })



#' @name refresh
#' @title refresh
#' @description
#' Updates gGraph object to reflect changes in source code
#' 
#' @param gGraph object
#'
#' @return gGraph object
refresh = function(gg)
{
    if(!inherits(gg, "gGraph")) {
        stop("gg is not a gGraph object")
    }
    return(gGraph$new(nodeObj = gg$nodes,
                      edgeObj = gg$edges))
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
alpha = function(col, alpha)
{
    col.rgb = col2rgb(col)
    out = rgb(red = col.rgb['red', ]/255, green = col.rgb['green', ]/255, blue = col.rgb['blue', ]/255, alpha = alpha)
    names(out) = names(col)
    return(out)
}



#' @name convertEdges
#' @title convertEdges
#' @description
#' Takes nodes and edges and converts the edge table for the nodes if they were strandless
#' gEdge table will be the appropriate table for if we did:
#'      nodes[which(as.logical(strand(nodes)=="+"))]
#'      strand(nodes) = "*"
#' pre: Nodes have snode.id and index column (indicates index)
#'      edges have sedge.id 
#'
#' @author Joe DeRose
#' @param nodes GRanges of signed nodes (ie $pnodes in gGraph object)
#' @param edges data.table of signed edges (ie $pedges in gGraph object)
#' @keywords internal
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


  ## Get positive nodes
  new.nodes = unname(nodes) %Q% (strand == "+")
  
  ## Create map between old id positions and new positions (pos - pos, neg - pos)
  ## Assumes no two values are not NA which they shouldnt be if everything is set up right on backend
  indexMap = pmin(match(nodes$snode.id, new.nodes$snode.id), match(nodes$snode.id, -new.nodes$snode.id), na.rm=T)
  indexMap = data.table(index = nodes$index, new.index = indexMap)
  setkey(indexMap, index)
  
  ## Map the edges to their new location
  es[, ":="(n1 = indexMap[.(from), new.index],
            n2 = indexMap[.(to), new.index])]

  ## Get only the positive non-NA nodes
  ## remove duplicate edges ie those that are present in both positive and negative orientations hence have
  ## abs(sedge.id) present twice
  es = es[!duplicated(abs(sedge.id)), ] 
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
}



## ================= gWalk class definition ================== ##
#' @export
gWalk = setClass("gWalk")

gWalk = R6::R6Class("gWalk", ## GWALKS
                     public = list(
                         ## snode.id is a list of snode.id vectors
                         ## sedge.id is a list of sedge.id vectors
                         ## grl is an arbitrary GRangesList representing a walk
                         ## graph is a pointer to the graph that snode.id/sedge.id comes from if being used
                         ## Only one of these can be non-NULL when initializing this object (graph must be used with snode.id and sedge.id)
                         initialize = function(snode.id = NULL,
                                               sedge.id = NULL,
                                               grl = NULL,
                                               graph = NULL,
                                               meta = NULL,
                                               collapse = FALSE ## flag whether to collapse ie disjoin the intervals
                                               )
                         {
                           if (all(is.null(snode.id), is.null(sedge.id), is.null(grl))) {
                             ## initializing empty gWalk

                             if (is.null(graph))
                               graph = gGraph$new()
                             
                             ## data.table of node.ids in the walk
                             private$pnode = gGraph$new()$nodes
                             
                             ## data.table of edge.ids in the walk
                             private$pedge = gGraph$new()$edges
                             
                             ## data.table of walk metadata 
                             private$pmeta = data.table()
                             
                             ## Pointer to the graph the snode.id/sedge.ids come from
                             ## If initialized with grl, this is the graph from the gWalk
                             private$pgraph = NULL

                             return(self)
                           } else if (sum(c(!is.null(snode.id), !is.null(sedge.id), !is.null(grl))) > 1) {
                             stop("More than one of snode.id, sedge.id and grl cannot be non-NULL")
                           }


                             if (!is.null(grl)) {

                               ## will create the graph first using the unsigned nodes / edges constructor
                               ## doing it differently depending on whether disjoin is flagged

                               grlu = grl.unlist(grl)
                               if (collapse) ## here we disjoin all nodes first and then build the graph
                               {
                                 nodes = disjoin(gr.stripstrand(grlu))
                               }
                               else ## otherwise we keep all nodes separate 
                               {
                                 nodes = gr.stripstrand(grlu[, c()])
                               }

                               tmpdt = as.data.table(grlu[, c('grl.ix', 'grl.iix')] %*% nodes)[, strand := as.character(strand(grlu)[query.id])]

                               ## pos will take of duplicated created in the "disjoin" case where a single walk interval comprises multiple nodes
                               ## i.e. has teh same grl.ix and grl.iix but different pos
                               ## we want to make sure that these get ordered properly in the key call below
                               ## so that on negative strands the leftward segments come after the rightward in the same grl.iix
                               tmpdt[, pos := (sign(strand=='+')-0.5)*start]
                               setkeyv(tmpdt, c("grl.ix", "grl.iix", 'pos'))

                               edges = tmpdt[, .(n1 = subject.id[-.N], n2 = subject.id[-1],
                                                 n1.side = sign(strand[-.N] == '+'),
                                                 n2.side = sign(strand[-1] == '-')), by = grl.ix]

                               edges$grl.ix = NULL
                               graph = gGraph$new(nodes = nodes, edges = edges)

                               snode.id = split(ifelse(tmpdt$strand=='+', 1, -1)*tmpdt$subject.id, tmpdt$grl.ix)
                             }


                           if ((!is.null(snode.id) | !is.null(sedge.id)) & is.null(graph))
                           {
                             stop('graph must be non NULL if snode.id or sedge.id are provided')
                           }

                           if (!is.null(snode.id)) {
                             ## List of snode.id's
                             ## Get list of sedge.ids with only snode.ids from each element
                             private$pmeta = data.table(walk.id = 1:length(snode.id))
                             if (!is.null(meta))
                             {
                               if (nrow(meta) != length(snode.id))
                               {
                                 stop('data.table of meta.data must be same length and order as node, edge, or GRangesList list input to gWalk constructor')
                               }
                               private$pmeta = cbind(data.table(walk.id = 1:length(snode.id)), meta)
                             }

                             ## always faster to unlist than lapply
                             pnode = dunlist(snode.id)[, .(walk.id = as.integer(listid), snode.id = V1)]

                             ## first need to check if the edges corresponding to the consecutive
                             ## node pairs in the input lists exist
                             pedge = pnode[, .(from = snode.id[-.N], to = snode.id[-1]), by = walk.id]
                             edgesdt = copy(graph$edgesdt)[, from := graph$gr$snode.id[from]][, to := graph$gr$snode.id[to]]
                             setkeyv(edgesdt, c('from', 'to'))
                             pedge$sedge.id = edgesdt[.(pedge$from, pedge$to), sedge.id]
                             if (any(is.na(pedge$sedge.id)))
                               stop('One or more provided walks refers to a non-existent edge')

                             private$pnode = pnode
                             private$pedge = pedge
                             private$pgraph = graph
                             
                           } else if (!is.null(sedge.id)) {

                             private$pmeta = data.table(walk.id = 1:length(sedge.id))
                             if (!is.null(meta))
                             {
                               if (nrow(meta) != length(sedge.id))
                               {
                                 stop('data.table of meta.data must be same length and order as node, edge, or GRangesList list input')
                               }
                               private$pmeta = cbind(private$pmeta, meta)
                             }

                             pedge = dunlist(sedge.id)[, .(walk.id = as.integer(listid), sedge.id = V1)]
                             tmpe = graph$edgesdt[.(pedge$sedge.id), .(from, to)]

                             if (any(is.na(tmpe$to)))
                               stop('one or more provided sedge.ids are out of bounds')

                             pedge = cbind(pedge, tmpe)

                             ## first need to check if consecutive edges in the provided
                             ## walk share sink and source nodes
                             pedge[, matches.next := c(to[-.N] == from[-1], NA), by = walk.id]

                             if (any(!pedge$matches.next, na.rm = TRUE))
                               stop('at least one pair of consecutive edges i and i+1 in the walk do not share a sink and source node, respectively')

                             pedge[, to := graph$dt$snode.id[to]]
                             pedge[, from := graph$dt$snode.id[from]]

                             pnode = pedge[, .(snode.id = c(from, to[.N])), by = walk.id]
                             private$pedge = pedge[, .(walk.id, sedge.id, to, from)]
                             private$pnode = pnode
                             private$pgraph = graph
                           } 
             
                             setkey(private$pmeta, walk.id)
                             setkey(private$pnode, walk.id)                             
                             setkey(private$pedge, walk.id)
                             
                             return (self)
                         },

                       ## Allows for subsetting of the gWalk Object using bracket notation
                       subset = function(i)
                       {
                         if (is.logical(i))
                           i = which(i)
                         
                         if (is.numeric(i) | is.integer(i))
                         {
                           if (max(i, na.rm = TRUE)>self$length())
                             stop('index out of bounds')
                         }


                         tmp.node = rbind(private$pnode, copy(private$pnode)[,
                                                 ":="(walk.id = -walk.id, snode.id = -snode.id)][rev(1:.N), ])
                         tmp.edge = rbind(private$pedge, copy(private$pedge)[,
                                                 ":="(walk.id = -walk.id, sedge.id = -sedge.id)][rev(1:.N), ])
                         tmp.meta = private$pmeta[.(abs(i)), ]

                         setkey(tmp.node, walk.id)
                         setkey(tmp.edge, walk.id)

                         new.node = merge(data.table(query.id = 1:length(i), walk.id = i), tmp.node, by = 'walk.id', allow.cartesian = TRUE)
                         new.edge = merge(data.table(query.id = 1:length(i), walk.id = i), tmp.edge, by = 'walk.id',, allow.cartesian = TRUE)
                         new.meta = merge(data.table(query.id = 1:length(i), walk.id = abs(i)), private$pmeta, by = 'walk.id', , allow.cartesian = TRUE)

                         new.node[, walk.id := query.id]
                         new.edge[, walk.id := query.id]
                         new.meta[, walk.id := query.id]

                         private$pnode = new.node[,-2]
                         private$pedge = new.edge[,-2]
                         private$pmeta = new.meta[, -2]

                         setkey(private$pnode, walk.id)
                         setkey(private$pedge, walk.id)
                         setkey(private$pmeta, walk.id)

                         return(self)
                       },

                       length = function() {
                         return(nrow(private$pmeta))
                       },

                       print = function()
                       {
                         TOOBIG = 5
                         TEASER = 2
                         .tease = function(stuff, toobig = TOOBIG, teaser = TEASER)
                         {
                           if (length(stuff)>toobig)
                             paste(paste(stuff[1:teaser], collapse = ' -> '), '-> ... ->', paste(stuff[length(stuff)-TEASER:1+1], collapse = ' -> '))
                           else
                             paste(stuff, collapse = ' -> ')
                         }

                         message('gWalk object with ', self$length(), ' walks')
                         if (self$length()>0)
                           {
                             ix = seq_along(self$length())
                             tmp.node = data.table(walk.id = private$pnode$walk.id,  nix = private$pgraph$queryLookup(private$pnode$snode.id)$index)
                             out = cbind(private$pmeta[.(ix), ], private$pnode[, .(nodes = .tease(snode.id)), keyby = walk.id][.(ix),][, -1],
                                         tmp.node[, .(gr = .tease(gr.string(private$pgraph$gr[nix], mb = FALSE))), keyby = walk.id][.(ix),][, -1],
                                         private$pedge[, .(edges = .tease(sedge.id)), keyby = walk.id][.(ix),][, -1])

                             print(out[1:min(nrow(out), TOOBIG)])
                             if (self$length()>TOOBIG)
                             {
                               message('\n ... \n(', self$length()-TOOBIG, ' more walks )')
                             }
                           }
                       },

                       mark = function(...)
                       {
                         self$nodes$mark(...)
                         self$edges$mark(...)
                       },
                       
                       gtrack = function(name = 'gWalk', stack.gap = 1e5, ...)
                       {
                         gTrack(self$grl, draw.paths = TRUE, stack.gap = stack.gap, name = name,  ...)
                       }
                     ),

                     private = list(
                       ## data.table of node.ids in the walk
                       pnode = NULL,
                       
                       ## data.table of edge.ids in the walk
                       pedge = NULL,

                       ## data.table of walk metadata 
                       pmeta = NULL,
                       
                       ## Pointer to the graph the snode.id/sedge.ids come from
                       ## If initialized with grl, this is the graph from the gWalk
                       pgraph = NULL
                     ),

                    active = list(
                      ## Returns a GRangesList of walks in the graph
                      grl = function()
                      {
                        ## Does not get both strands only 1 strand
                        nix = private$pgraph$queryLookup(private$pnode$snode.id)$index
                        return(split(private$pgraph$gr[nix], private$pnode$walk.id))
                      },

                      nodes = function()
                      {
                        return(private$pgraph$nodes[private$pnode$snode.id])
                      },

                      edges = function()
                      {
                        return(private$pgraph$edges[private$pedge$sedge.id])
                      },

                      footprint = function()
                      {
                        return(sort(reduce(gr.stripstrand(self$nodes$gr))))
                      },

                      gt = function()
                      {
                        return(self$gtrack())
                      },
                      
                      ## Allows the user to "port" the graph and then access it (if generated from grl, allows us to get the dummy graph we made)
                      graph = function()
                      {
                        return(private$pgraph)
                      },
                      
                      dt = function() {
                        ix = 1:self$length()
                        out = cbind(private$pmeta[.(ix), ],
                                    private$pnode[, .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),][, -1],
                                    private$pedge[, .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),][, -1])
                        return(out)
                      }
                    )
                    )


#' @name c
#' @title c
#' @description
#' 
#' Concatenates gWalk objects by id's
#'
#' @param gWalk objects
#'
#' @return a new concatenated gWalk Object
#' @export
`c.gWalk` = function(...)
{                            
    gWalk.list=list(...)
    isg = sapply(gWalk.list, function(x) class(x)[1]=='gWalk')

    if(any(!isg)){
        stop('Error: All inputs must be of class gNode.')
    }
    
    ##Check to make sure they all come from the same graph
    graphs = lapply(gWalk.list, function(x) x$graph)
    if(length(unique(graphs))>1){
        stop('Error: All gNodes must come from the same graph')
    }
    
    ##Get all the pnode.id's to create new gNode
    sedge.id = do.call('c', lapply(gWalk.list, function(x) x$dt$sedge.id))

    return (gWalk$new(sedge.id = sedge.id, graph = gWalk.list[[1]]$graph))
}



#' @name [.gWalk
#' @title gWalk
#' @description
#'
#' Overloads subset operator for walks
#'
#' @param obj gWalk object This is the gWalk object to be subset
#' @param i integer, logical, or expression on gWalk metadata used to subset gWalk
#' @return A new gWalk object that contains only the given subset of walks
#' @export
'[.gWalk' = function(obj, i = NULL){
  walks = obj$clone()
  inew = tryCatch(with(walks$dt, eval(parse(text = lazyeval::expr_text(i)))), error = function(e) NULL)
  if (!is.null(inew))
    i = inew
  walks$subset(i)
  return(walks)
}


#' @name dunlist
#' @title dunlist
#'
#' @description
#' unlists a list of vectors, matrices, data.tables into a data.table indexed by the list id
#' $listid
#'
#' does fill = TRUE in case the data.tables inside the list do not have compatible column names 
#' 
#' @param x list of vectors, matrices, or data frames
#' @return data.frame of concatenated input data with additional fields $ix and $iix specifying the list item and within-list index from which the given row originated from
#' @author Marcin Imielinski
#' @export
#' @keywords internal
#' @noRd 
#############################################################
dunlist = function(x)
{
  if (is.null(names(x)))
    names(x) = 1:length(x)
  tmp = lapply(x, as.data.table)
  
  out = cbind(data.table(listid = rep(names(x), elementNROWS(x)), rbindlist(tmp, fill = TRUE)))
  setkey(out, listid)
  return(out)
}
