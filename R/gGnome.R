#' @name  gGnome
#' @description
#' Reference-based graph representation of structurally altered genome
#' employing GenomicRanges framework.
#'
#'
#' Copyright (C) 2018 Xiaotong Yao, Joseph DeRose, Marcin Imielinski
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
#' @importFrom reshape2 melt
#' @import gUtils
#' @import gTrack
"_PACKAGE"


## ================= gNode class definition ================== ##
#' 
#' @export
gNode = setClass("gNode")
gNode = R6::R6Class("gNode",
                    public = list(

                        #' @name initialize
                        #' @description
                        #' Set up the constructor GNODE
                        #' index - snode.id
                        #' graph
                        #' @param snode.id signed node id
                        #' @param graph gGraph associated with this sNode
                        #' @author Joe DeRose
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

                        #' @name mark
                        #' @description
                        #'
                        #' Marks the nodes pointed at by this gNode by adding metadata in the gGraph they came from.
                        #' Metadata names can be arbitrary, such as "lwd = 7" or "highlight = TRUE".
                        #' gGraph nodes have metadata specified appended to them
                        #'
                        #' @param ... metadata names = data to store in metadata columns
                        #' 
                        #' @usage
                        #'
                        #' gg$nodes[1:5]$mark(col = "purple")
                        #' gg$nodes$mark(changed = FALSE)
                        #' @param  ... name = value pairs of scalar or vector (length edges in graph) arguments
                        #' @author Marcin Imielinski
                        mark = function(...)
                        {
                            NONO.FIELDS = c('node.id', 'snode.id', 'index', 'loose.left', 'loose.right')
                            args = list(...)

                            if (any(names(args) %in% NONO.FIELDS)){
                                stop(paste('Cannot modify these reserved gNode metadata fields:', paste(NONO.FIELDS, collapse = ', ')))
                            }

                            for (nm in names(args)){
                                private$pgraph$annotate(nm, args[[nm]],
                                                        private$pnode.id, "node")
                            }
                        },

                        #' @name subset
                        #' @description
                        #' Allows subseting of the Node object using bracket notation
                        #'
                        #' @param i integer or logical index
                        #' @return subset of nodes indexed by i 
                        #' @author Joe DeRose
                        subset = function(i)
                        {

                            if (is.null(i)){
                                i = integer()
                            }
                            if (is.logical(i)){
                                i = which(i)
                            }
                            
                            if (length(i)>0 && (is.numeric(i) | is.integer(i))) {
                                if (max(abs(i), na.rm = TRUE) > self$length || max(abs(i), na.rm = TRUE) == 0) {
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
                                {
                                    stop('one or more signed node ids are not found in this object')
                                }
                            }
                            
                            private$pnode.id = private$pnode.id[abs(i)]
                            new.pindex = ifelse(i>0, private$pindex[abs(i)], private$prindex[abs(i)])
                            new.prindex = ifelse(i>0, private$prindex[abs(i)], private$pindex[abs(i)])
                            private$pindex = new.pindex
                            private$prindex = new.prindex
                            private$porientation = private$porientation[abs(i)]*sign(i)
                            return(self)                            
                        },



                        #' @name ego
                        #' @description
                        #'
                        #' returns all nodes with k degrees of sepraation (i.e. "order")
                        #' from these
                        #'
                        #' @param order integer edge distance from these (seed) nodes to return
                        #' @param mindist minimum distance from these (seed) nodes
                        #' @return gNode object comprising nodes within the ego of this node
                        #' @author Marcin Imielinski
                        ego = function(order = 0, mindist = 0)
                        {
                            G = self$graph$igraph
                            nodes = c(private$pindex, private$prindex)
                            egoi = unique(unlist(ego(G, order = order, nodes = nodes, mindist = mindist)))
                            egoini = unique(self$graph$gr[egoi]$node.id)

                            return(self$graph$nodes[egoini])
                        },
                        
                        #' @name print
                        #' @description
                        #' 
                        #' Prints out the gNode Object. Prints the length and the GRanges of the nodes.
                        print = function()
                        {
                            message('gNode object of length ', self$length)
                            print(self$gr)
                        }
                    ),
                    

                    private = list(
                        ## stores the unsigned node.id's of the nodes in private$pgraph in this gNode
                        pnode.id = NULL,
                        
                        ## Stores the index and reverse index of the pnode.id's in this gNode
                        pindex = NULL,
                        prindex = NULL,
                        
                        ## Stores the signs of the nodes associated with indices in private$pindex
                        porientation = NULL,
                        
                        ## Pointer to gGraph that this gNode belongs to
                        pgraph = NULL
                    ),


                    active = list(
                        #' @name length
                        #' @description
                        #' 
                        #' Returns the number of nodes in this gNode Object.
                        #' 
                        #' @return Number of nodes in this gNode Object
                        length = function()
                        {
                            return(length(private$pnode.id))
                        },

                        copy = function() self$clone(),

                        #' @name graph
                        #' @description
                        #'
                        #' Returns the gGraph this gNode points at.
                        #'
                        #' @return gGraph this gNode points at
                        graph = function()
                        {
                            return(private$pgraph)
                        },


                        #' @name dt
                        #' @description
                        #'
                        #' Returns the GRanges of the nodes in the gNode as a data.table.
                        #'
                        #' @return data.table GRanges of the nodes coverted to a data.table
                        dt = function() {
                            return(as.data.table(private$pgraph$gr[private$pindex]))
                        },
                        
                        
                        #' @name gr
                        #' @description
                        #'
                        #' Returns a GRanges of the nodes in the gNode.
                        #'
                        #' @return GRanges of the nodes
                        gr = function() {
                            return(private$pgraph$gr[private$pindex])
                        },
                        

                        #' @name id
                        #' @description
                        #'
                        #' Returns the node.id's of the nodes in this gNode
                        #'
                        #' @return Vector of the node.id's in this gNode
                        id  = function()
                        {
                            return(private$pnode.id)
                        },


                        #' @name sid
                        #' @description
                        #'
                        #' Returns the snode.id's of the nodes in this gNode
                        #'
                        #' @return Vector of the snode.id's in this gNode
                        sid = function()
                        {
                            return(ifelse(private$porientation == 1, private$pnode.id, -private$pnode.id))
                        },                       

                        ## returns flipped version of this node
                        flip = function()
                        {
                            sid = self$dt$snode.id
                            return(self$graph$nodes[-sid])
                        },

                        ## Returns the nodes connected to the left of the nodes

                        #' @name left
                        #' @description
                        #'
                        #' Returns a gNode containing the nodes connected to the left of the nodes in this
                        #' gNode. If there are no nodes to the left, returns an empty gNode.
                        #'
                        #' @return gNode Nodes connected to the left of the nodes in this gNode Object
                        left = function()
                        {
                            dt = copy(private$pgraph$sedgesdt)
                            setkeyv(dt, "to")
                            
                            tmp = dt[.(private$pindex), from]
                            tmp = tmp[!is.na(tmp)]

                            leftNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[tmp],
                                                  graph = private$pgraph)                       
                            return(leftNodes)
                        },


                        #' @name right
                        #' @description
                        #'
                        #' Returns a gNode containing the nodes connected to the right of the nodes in this
                        #' gNode. If there are no nodes to the right, returns an empty gNode.
                        #'
                        #' @return gNode Nodes connected to the right of the nodes in this gNode
                        right = function()
                        {
                            dt = copy(private$pgraph$sedgesdt)
                            setkeyv(dt, "from")
                            
                            tmp = dt[.(private$pindex), to]
                            tmp = tmp[!is.na(tmp)]
                            
                            rightNodes = gNode$new(snode.id = private$pgraph$gr$snode.id[tmp],
                                                   graph = private$pgraph)                       
                            return(rightNodes)
                        },

                        
                        #' @name edges
                        #' @description
                        #'
                        #' Returns a gEdge containing the all edges connected to the nodes in this gNode.
                        #' If no edges are associated with this gNode, returns an empty gEdge.
                        #'
                        #' @return gEdge Edges connected to the nodes in this gNode
                        edges = function()
                        {
                            ed = c(self$eleft, self$eright)
                            if (length(ed)>0)
                                ed =ed[!duplicated(edge.id)]

                            return(ed)
                        },

                        #' @name eleft
                        #' @description
                        #'
                        #' Returns a gEdge containing the all edges connected to the left of the nodes in this gNode.
                        #' If no edges are left of this gNode, returns an empty gEdge.
                        #'
                        #' @return gEdge Edges connected to the left of the nodes in this gNode
                        eleft = function()
                        {
                            sedge.id = private$pgraph$sedgesdt[to %in% private$pindex, sedge.id]
                            return(gEdge$new(sedge.id, graph = private$pgraph))
                        },
                        

                        #' @name eright
                        #' @description
                        #'
                        #' Returns a gEdge containing the all edges connected to the right of the nodes in this gNode.
                        #' If no edges are right of this gNode, returns an empty gEdge.
                        #'
                        #' @return gEdge Edges connected to the left of the nodes in this gNode
                        eright = function()
                        {
                            sedge.id = private$pgraph$sedgesdt[from %in% private$pindex, sedge.id]
                            return(gEdge$new(sedge.id, graph = private$pgraph))
                        },

                        
                        #' @name loose
                        #' @description
                        #' 
                        #' Returns a gNode containing the loose ends connected to the nodes in this
                        #' gNode. If there are no loose ends, returns an empty gNode.
                        #'
                        #' @return GRanges Loose ends connected to the nodes in this gNode Object
                        loose = function()
                        {
                            return(unique(c(self$lleft, self$lright)))
                        },

                        #' @name loose.left
                        #' @description
                        #' sets or accesses a logical vector specifying whether the given node has a loose.end on its left 
                        #'
                        #' @param value logical value to replace loose end value by
                        #' @author Marcin Imielinski
                        loose.left = function(value)
                        {
                            if (!missing(value))
                            {
                                if (!is.logical(value)){
                                    stop('replacement must be logical vector')}

                                if (any(is.na(value))){
                                    stop('replacement cannot contain NAs')}

                                ldeg = self$ldegree

                                value = value | ldeg==0

                                ix = private$porientation>0
                                if (any(ix)){                        
                                    private$pgraph$annotate('loose.left', value[ix],
                                                            private$pnode.id[ix], "node")}
                                
                                if (any(!ix)){
                                    private$pgraph$annotate('loose.right', value[!ix],
                                                            private$pnode.id[!ix], "node")}
                                ##FIXME: R6 doesn't like something wrong w setting active binding here

                            } else
                            {
                                return(ifelse(private$porientation>0,
                                              self$gr$loose.left, self$gr$loose.right))
                            }
                        },


                        #' @name loose.right
                        #' @description
                        #' sets or accesses a logical vector specifying whether the given node has a loose.end on its right
                        #'
                        #' @param value logical value to replace loose end value by
                        #' @author Marcin Imielinski
                        loose.right = function(value)
                        {
                            if (!missing(value))
                            {
                                if (!is.logical(value)){
                                    stop('replacement must be logical vector')
                                }
                                if (any(is.na(value))){
                                    stop('replacement cannot contain NAs')
                                }
                                rdeg = self$rdegree

                                value = value | ldeg==0

                                ix = private$porientation>0
                                if (any(ix)){                
                                    private$pgraph$annotate('loose.right', value[ix],
                                                            private$pnode.id[ix], "node")
                                }
                                if (any(!ix)){
                                    private$pgraph$annotate('loose.left', value[!ix],
                                                            private$pnode.id[!ix], "node")
                                }
                                ##FIXME: R6 doesn't like something wrong w setting active binding here

                            } else
                            {
                                return(ifelse(private$porientation>0,
                                              self$gr$loose.right, self$gr$loose.left))
                            }
                        },

                        degree = function()
                        {
                            ldeg = private$pgraph$edgesdt[n2.side == 'left', .N, keyby = n1][private$pnode.id, N] + private$pgraph$edgesdt[n2.side == 'left', .N, keyby = n2][private$pnode.id, N]
                            
                            rdeg = private$pgraph$edgesdt[n2.side == 'right', .N, keyby = n1][private$pnode.id, N] + private$pgraph$edgesdt[n2.side == 'right', .N, keyby = n2][private$pnode.id, N]

                            return(rowSums(cbind(ldeg, rdeg, 0), na.rm = TRUE))
                        },

                        terminal = function()
                        {
                            if (self$length==0)
                                return(GRanges(seqlengths = seqlengths(self)))
                            ix = which(self$ldegree==0 | self$rdegree==0)
                            return(self[ix]$loose)
                        },

                        ldegree = function()
                        {
                            ldeg = rowSums(cbind(private$pgraph$edgesdt[, sum(n1.side == 'left'), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == 'left'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)

                            rdeg = ldeg
                            if (any(private$porientation<0)){
                                rdeg = rowSums(cbind(private$pgraph$edgesdt[ , sum(n1.side == 'right'), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == 'right'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)}

                            return(pmax(0, ifelse(private$porientation>0, ldeg, rdeg), na.rm = TRUE))
                        },

                        rdegree = function()
                        {
                            rdeg = rowSums(cbind(private$pgraph$edgesdt[, sum(n1.side == 'right'), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == 'right'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)
                            ldeg = rdeg
                            if (any(private$porientation<0)){
                                ldeg = rowSums(cbind(private$pgraph$edgesdt[, sum(n1.side == 'left'), keyby = n1][.(private$pnode.id), V1],  private$pgraph$edgesdt[, sum(n2.side == 'left'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)}

                            return(pmax(0, ifelse(private$porientation>0, rdeg, ldeg), na.rm = TRUE))
                        },                        

                        #' @name lleft
                        #' @description
                        #' 
                        #' Returns a gNode containing the loose ends connected to the left of the nodes in this
                        #' gNode. If there are no loose ends to the left, returns an empty gNode.
                        #'
                        #' @return GRanges Loose ends connected to the left of the nodes in this gNode Object
                        lleft = function()
                        {
                            return(gr.start(self$gr[self$gr$loose.left>0]))
                        },
                        
                        
                        #' @name lright
                        #' @description
                        #' 
                        #' Returns a gNode containing the loose ends connected to the right of the nodes in this
                        #' gNode. If there are no loose ends to the right, returns an empty gNode.
                        #'
                        #' @return GRanges Loose ends connected to the right of the nodes in this gNode Object
                        lright = function()
                        {
                            return(gr.end(self$gr[self$gr$loose.right>0]))
                        },

                        #' @name dist
                        #' @description
                        #' 
                        #' returns a distance matrix of all node pairs in this gGnode
                        #' representing base pairs separating each node pair
                        #'
                        #' @author Marcin Imielinski
                        dist = function()
                        {
                            return(private$pgraph$dist(private$pnode.id))
                        },

                        #' @name subgraph
                        #' @description
                        #'
                        #' Generates a new gGraph from the gGraph pointed at by this gNode containing only the nodes
                        #' in this gNode and the edges that connect them. The only edges included in this subgraph are
                        #' those that only contain the nodes in this gNode.
                        #'
                        #' @return gGraph Subgraph of all of the nodes in this gNode and the edges containing them
                        subgraph = function()
                        {
                            ## Get all the edges that connect our nodes and remove the duplicate edge.id's as they will
                            ## be automatically filled in in the Edge Constructor
                            edge.id = private$pgraph$sedgesdt[to %in% c(private$pindex, private$prindex)
                                                              & from %in% c(private$pindex, private$prindex), edge.id]
                            edge.id = edge.id[!duplicated(edge.id)]
                            
                            edgeObj = gEdge$new(seid = edge.id,
                                                graph = private$pgraph)
                            return(gGraph$new(nodeObj = self, edgeObj = edgeObj, meta = private$pgraph$meta))
                        }
                    )
                    )


## ================== Non-Member Functions for gNode ================== ##

#' @name length
#' @title length.gNode
#' @description
#'
#' The number of nodes in the gNode Object
#'
#' @param gNode a gNode object
#' @return the number of nodes in the gNode Object
#' @author Joe DeRose
#' @export
`length.gNode` = function(gNode)
{
    if(!inherits(gNode, "gNode")) {
        stop("Error: Invalid input")
    }
    return(gNode$length)
}


#' @name c
#' @title c.gNode
#' @description
#' 
#' Concatenates gNode objects by node.id's stored in them. All arguments
#' must point at the same graph or an error will be thrown.
#'
#' @param gNode objects
#' @return a new concatenated gNode
#' @author Joe DeRose
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
#' 
#' Returns a new gNode which is the difference between x and y (node.id's).
#' All arguments must point at the same graph or an error will be thrown.
#'
#' @param x a gNode Object
#' @param y a gNode Object
#' @return new gNode Object containing the difference between x and y
#' @author Joe DeRose
#' @export
setMethod("setdiff", c("gNode", "gNode"),
          function(x, y)
          {  
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = setdiff(x$id, y$id)
              return(gNode$new(new.ids, x$graph))
          })

#' @name union
#' @title union.gNode
#' @description
#' 
#' Returns a new gNode which is the union of x and y (node.id's). All
#' arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#' @return new gNode Object containing the union of x and y
#' @author Joe DeRose
#' @export
#'

setMethod("union", c("gNode", "gNode"), function(x, y, ...)
{
    if(!identical(x$graph, y$graph)) {
        stop("Arguments do not point to the same graph")
    }

    new.ids = union(x$id, y$id)
    return(gNode$new(new.ids, x$graph))
})

#' @name intersect
#' @title intersect.gNode
#' @description
#' 
#' Returns a new gNode which is the intersection of x and y (node.id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gNode Object
#' @param y a gNode Object
#' @return new gNode Object containing the intersection of x and y

#' @author Joe DeRose
#' @export
setMethod("intersect", c("gNode", "gNode"),

          function(x, y)
          {             
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
#' Overloads subset operator for gNode. Allows subsetting of gNode via index or
#' snode.id (as a key). Also, allows for data.table style queries on gNode GRanges
#' and corresponding metadata.
#'
#' @param obj This is the gNode object to be subset
#' @param i integer, logical, or expression in gNode metadata used to subset gNodes
#' @return A new gNode that contains only the given id's
#' @author Joe DeRose
#' @export
'[.gNode' = function(obj, i = NULL, ...)
{
  nodes = obj$clone()
  ## yes I finally figured it out!!!! grrrrrr
  inew = tryCatch(
      eval(eval(parse(text = substitute(deparse(substitute(i)))),
                parent.frame()),
           nodes$dt, parent.frame(2)),
      error = function(e) NULL)
  if (is.null(inew)){
      inew = i ## just give up
  }
  nodes$subset(inew)
  return(nodes)
}


## ================= gEdge class definition ================== ##
#'
#' @export
gEdge = setClass("gEdge")
gEdge = R6::R6Class("gEdge",
                    public = list(

                        #' @name gEdge constructor
                        #' @description
                        #' Constructor GEDGE
                        #' Builds an gEdge Class from a given edge.id 
                        #' graph - the gGraph you want to build this class from
                        #' @param seid vector of signed integer edge ids
                        #' @param graph gGraph that this edge belongs to
                        #' @author Joe DeRose
                        initialize = function(seid = NULL, graph)
                        {
                            private$pgraph = graph ## reference to graph that these edges belong to
                            private$porientation = private$pedge.id = private$psedge.id = c()                       
                            
                            if (is.null(seid)) {
                                return(self)
                            }
                            
                            if (any(is.na(private$pgraph$sedgesdt[.(seid), edge.id]))) {
                                stop("one or more provided signed edge ids are out of bounds")
                            }
                            
                            private$porientation = rep(1, length(seid))
                            private$psedge.id = seid ## signed edge id, referring to the edge in the directed graph for which the $left direction is 5' and $right direction is '
                            private$pedge.id = abs(seid) ## unsigned edge id
                            private$pedges = private$pgraph$sedgesdt[list(seid), ]
                            return(self)
                        },


                        #' @name mark
                        #' @description
                        #'
                        #' Marks the nodes pointed at by this gEdge by adding metadata in the gGraph they came from.
                        #' Metadata names can be arbitrary, such as "lwd = 7" or "highlight = TRUE".
                        #' gGraph nodes have metadata specified appended to them
                        #' 
                        #' @param ... metadata names = data to store in metadata columns
                        #' 
                        #' @usage
                        #'
                        #' gg$edges[1:5]$mark(col = "purple")
                        #' gg$edges$mark(changed = FALSE)
                        #' @param  ... name = value pairs of scalar or vector (length edges in graph) arguments
                        #' @author Marcin Imielinski
                        mark = function(...)
                        {
                            NONO.FIELDS = c('from', 'to', 'sedge.id', 'edge.id', 'type')
                            args = list(...)
                            
                            if (any(names(args) %in% NONO.FIELDS)){
                                stop(paste('Cannot modify the following reserved gEdge metadata fields:', paste(NONO.FIELDS, collapse = ', ')))
                            }
                            
                            for (nm in names(args)){
                                private$pgraph$annotate(nm, args[[nm]],
                                                        c(private$psedge.id, -private$psedge.id), "edge")                          }
                        },


                        #' @name subset
                        #' @description
                        #' Allows subseting of the gNode object using bracket notation
                        #'
                        #' @param i integer or logical index
                        #' @return subset of nodes indexed by i 
                        #' @author Joe DeRose
                        #'  Allows for subsetting of the gEdge Object using bracket notation
                        subset = function(i)
                        {
                            if (is.null(i)){
                                i = integer()
                            }
                            if (is.logical(i)){
                                i = which(i)
                            }
                            if ((length(i)>0) && (is.numeric(i) | is.integer(i)))
                            {
                                if (max(i, na.rm = TRUE)>self$length){
                                    stop('index out of bounds')
                                }
                            }
                            
                            if (is.character(i))
                            {
                                i = as.integer(i)
                                i = sign(i)*match(abs(i), abs(private$psedge.id))
                            }
                            
                            private$psedge.id = sign(i)*private$psedge.id[abs(i)]
                            private$pedges = private$pgraph$sedgesdt[.(private$psedge.id), ]
                            private$pedge.id = private$pedges$edge.id
                            private$porientation = sign(i)*private$porientation[abs(i)]
                            
                            return(self)
                        },
                        

                        #' @name print
                        #' @description
                        #' 
                        #' Prints out the gEdge. Prints the length and the data.table of the edges.
                        #' Prints the edges in the form n1, n2, n1.side, n2.side instead of to, from
                        print = function()
                        {
                            message("gEdge object with ", self$length, " edges")
                            print(self$dt)
                        }
                    ),


                    private = list(
                        ## data.table of the edges in this gEdge
                        pedges = NULL,

                        ## edge.id's and sedge.id's of the edges in this gEdge
                        pedge.id = NULL,
                        psedge.id = NULL,
                        
                        ## Sign of our sedge.id's
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

                        #' @name graph
                        #' @description
                        #'
                        #' Returns the gGraph this gEdge points at.
                        #'
                        #' @return gGraph this gEdge points at
                        graph = function()
                        {
                            return(private$pgraph)
                        },


                        copy = function() self$clone(),

                        #' @name length
                        #' @description
                        #' 
                        #' Returns the number of edge pairs in this gEdge
                        #' 
                        #' @return Number of edge pairs in this gEdge
                        length = function()
                        {                            
                            return(length(private$pedge.id))
                        },

                        #' @name span
                        #' @description
                        #' 
                        #' Returns a vector of spans of the junctions referred to by these
                        #' edges (i.e. reference distance between the breakpoints)
                        #' 
                        #' @return Number of edge pairs in this gEdge
                        span = function()
                        {                            
                            return(self$junctions$span)
                        },
                        
                        #' @name sign
                        #' @description
                        #' 
                        #' Returns a vector of signs of the junctions referred to by these
                        #' edges (i.e. the product of the numeric interpretation of the strands)
                        #' 
                        #' @return Number of edge pairs in this gEdge
                        sign = function()
                        {                            
                            return(self$junctions$sign)
                        },
                        
                        #' @name left
                        #' @description
                        #'
                        #' Returns a gNode containing the nodes connected to the left (from's) of the edges in this
                        #' gEdges. If there are no nodes to the left, returns an empty gNode.
                        #'
                        #' @return gNode Nodes connected to the left of the edges in this gNode Object
                        left = function()
                        {                      
                            leftNodes = private$pedges[, from]
                            leftNodes = gNode$new(snode.id=private$pgraph$dt[, snode.id[leftNodes]], private$pgraph)
                            return(leftNodes)
                        },


                        #' @name right
                        #' @description
                        #'
                        #' Returns a gNode containing the nodes connected to the right (to's) of the edges in this
                        #' gNode. If there are no nodes to the right, returns an empty gNode.
                        #'
                        #' @return gNode Nodes connected to the right of the nodes in this gNode
                        ## Returns the nodes connected to the right of the nodes
                        right = function()
                        {
                            rightNodes = private$pedges[, to]
                            return(gNode$new(snode.id = private$pgraph$dt[,snode.id[rightNodes]],
                                             private$pgraph))
                        },


                        #' @name nodes
                        #' @description
                        #'
                        #' Returns a gNode containing the all nodes connected to the edges in this gEdges.
                        #' If no nodes are associated with this gEdge, returns an empty gNode.
                        #'
                        #' @return gNode Nodes connected to the nodes in this gEdge
                        nodes = function()
                        {
                            no = c(self$left, self$right)
                            if (length(no)>0)
                                no = no[!duplicated(node.id)]
                            return(no)
                        },

                        #' @name sdt
                        #' @description
                        #'
                        #' Returns a data.table of the signed edges in this gEdge in backend format (to, from, type, edge.id)
                        #'
                        #' @return data.table of the signed edges in this gEdge in to,from,type,edge.id form
                        sdt = function()
                        {
                            return(copy(private$pedges))
                        },


                        #' @name dt
                        #' @description
                        #'
                        #' Returns a data.table of the unsigned edges in this gEdge in front end format 
                        #'
                        #' @return data.table of the unsigned edges in this gEdge 
                        dt = function()
                        {
                            sides = c('left', 'right')
                            return(copy(convertEdges(private$pgraph$gr, private$pedges, metacols = TRUE, cleanup = FALSE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
                        },

                        #' @name junctions
                        #' @description
                        #'
                        #' Coverts the edges in this gEdge to junctions and returns these as a Junctions Object.
                        #' See Junctions Object documentation to see how junctions are defined in this package.
                        junctions = function()
                        {
                            grl = self$grl
                            return(Junction$new(grl))
                        },
                        

                        grl = function()
                        {
                            if (length(self)==0)
                            {
                                ## awkward but best constructor for empty grl with proper seqinfo
                                return(GRangesList(self$graph$gr[c()])[c()])
                            }

                            gr1 = gr.flipstrand(gr.end(private$pgraph$gr[private$pedges$from], ignore.strand = FALSE))
                            gr2 = gr.start(private$pgraph$gr[private$pedges$to], ignore.strand = FALSE)
                            grl = split(c(gr1, gr2), rep(1:length(gr1), 2))[as.character(1:length(gr1))]
                            names(grl) = private$edges$sedge.id
                            meta = cbind(private$pedges, data.table(bp1 = gr.string(gr1), bp2 = gr.string(gr2)))
                            values(grl) = meta[, unique(c("edge.id", "sedge.id", "from", "to", "bp1", "bp2", colnames(meta))), with = FALSE]

                            return(grl)
                        },
                        
                        #' @name id
                        #'
                        #' Returns the edge.id's of the edges in this gEdge
                        #'
                        #' @return edge.id's of the edges in this gEdge
                        id = function()
                        {
                            return(private$pedge.id)
                        },


                        #' @name subgraph
                        #' @description
                        #'
                        #' Generates a new gGraph from the gGraph pointed at by this gEdge containing only the edges
                        #' in this gEdgges and the nodes that connect them. The only nodes included in this subgraph are
                        #' those that only contain the edges in this gEdge.
                        #'
                        #' @return gGraph Subgraph of all of the edges in this gEdges and the nodes containing them
                        subgraph = function()
                        {
                            ## Get all of the nodes in either to or from
                            nodeIndex = c(private$pedges[,to], private$pedges[,from])
                            
                            ## Get the unsigned node.id's and remove any duplicate nodes
                            node.id = private$pgraph$dt[index %in% nodeIndex, node.id]
                            node.id = node.id[!duplicated(node.id)]

                            nodeObj = gNode$new(snode.id = node.id,
                                                graph = private$pgraph)
                            return(gGraph$new(nodeObj = nodeObj, edgeObj = self, meta = private$pgraph$meta))
                        }
                        

                    )
                    )


## ================== Non-Member Functions for gEdge ================== ##
#' @name [.gEdge
#' @title gEdge
#' @description
#'
#' Overloads subset operator for gEdge. Allows subsetting of gEdge via index or
#' sedge.id (as a key). Also, allows for data.table style queries on gEdge and
#' corresponding metadata.
#'
#' @param obj This is the gEdge object to be subset
#' @param i integer, logical, or expression in gEdge metadata used to subset gEdge
#' @return A new gEdges that contains only the given id's
#' @author Joe DeRose
#' @export
'[.gEdge' = function(obj, i = NULL, ...)
{
    edges = obj$clone()  
    ## yes I finally figured it out grrrrrr
    inew = tryCatch(
        eval(eval(parse(text = substitute(deparse(substitute(i)))),
                  parent.frame()),edges$dt, parent.frame(2)),
        error = function(e) NULL)
    if (is.null(inew)){
        inew = i ## just give up
    }
    edges$subset(inew)
    return(edges)
}


#' @name length
#' @title length.gEdge
#' @description
#' 
#' The number of edge pairs in the gEdge Object
#'
#' @param gNode a gEdge object
#' @return the number of edge pairs in the gEdge Object
#' @author Joe DeRose
#' @export
`length.gEdge` = function(gEdge)
{
    if(!inherits(gEdge, "gEdge")) {
        stop("Error: Invalid input")
    }
    return(gEdge$length)
}


#' @name c
#' @title c.gEdge
#' @description
#'
#' Concatenates gEdge by edge.id's stored in them. All arguments
#' must point at the same graph or an error will be thrown.
#'
#' @param gEdge objects
#' @return a new concatenated gEdge Object
#' @author Joe DeRose
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
#' @title setdiff.gEdge
#' @description
#' 
#' Returns a new gEdge which is the difference between x and y (edge.id's).
#' All arguments must point at the same graph or an error will be thrown.
#'
#' @param x a gEdge Object
#' @param y a gEdge Object
#' @export
#' @author Joe DeRose
#' @return new gEdge containing the difference between x and y
setMethod("setdiff", c("gEdge", "gEdge"),
          function(x, y)
          {
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }
              
              new.ids = setdiff(x$id, y$id)
              return(gEdge$new(new.ids, x$graph))
          })

#' @name intersect
#' @title intersect.gEdge
#' @description
#'
#' Returns a new gEdge which is the intersection of x and y (edge.id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gEdge Object
#' @param y a gEdge Object
#' @export
#' @author Joe DeRose
#' @return new gEdge containing the intersection of x and y
setMethod("intersect", c("gEdge", "gEdge"),
          function(x, y)
          {
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }

              new.ids = intersect(x$id, y$id)
              return(gEdge$new(new.ids, x$graph))
          })


#' @name union
#' @title union.gEdge
#' @description
#' 
#' Returns a new gEdge which is the union of x and y (edge.id's).
#' All arguments must point at the same graph or an error will be thrown.
#' 
#' @param x a gEdge Object
#' @param y a gEdge Object
#' @export
#' @return new gEdge containing the union of x and y
setMethod("union", c("gEdge", "gEdge"), function(x, y, ...)
{
    if(!identical(x$graph, y$graph)) {
        stop("Arguments do not point to the same graph")
    }  
    new.ids = union(x$id, y$id)
    return(gEdge$new(new.ids, x$graph))
})


## ================== Junction class definition ================== ##
#' 
#' @exportClass Junction
#' @export
Junction = setClass("Junction")
Junction = R6::R6Class("Junction",
                       public = list(
                           #' @name Junction class constructor
                           #' @description 
                           #' Builds junction class, grl must be GRangesList with each GRanges of length 2 JUNCTION
                           #' Empty junctions are removed
                           #' If grl is empty GRangesList, Junctions class is empty
                           #' @param grl GRangesList of signed breakpoint pairs
                           #' @author Joe DeRose and Rick Mortensen                         
                           initialize = function(grl = NULL, ...)
                           {
                               ## Check to make sure input is GRangesList with each element of length 2
                               if (is.character(grl))
                                   grl = read.juncs(grl, ...)

                               if (is.null(grl))
                                   grl = GRangesList()

                               ## Will allow elements to be empty and will just remove them
                               if (!inherits(grl, "GRangesList")) {
                                   stop("Input is not a GRangesList")
                               }
                               
                               widths = elementNROWS(grl)
                               empty = widths == 0
                               
                               ix = widths != 2
                               if(any(ix[!empty])) {
                                   stop(paste0("grl contains junctions that are of improper length. Indices are: ", paste(ix[!empty], collapse=" ")))
                               }
                               
                               private$pjuncs = grl[!empty]
                               
                               return(self)
                           },
                           
                           #' @name removeDups
                           #' @description
                           #'
                           #' Removes duplicate junctions in this Junction Object
                           #'
                           #' this Junction Object has no duplicate junctions
                           #' @author Rick Mortensen
                           removeDups = function()

                           {
                               tmp = gr2dt(sort(grl.unlist(private$pjuncs)))
                               tmp = tmp[, paste(seqnames, start, end, strand, collapse = ','), by = grl.ix]
                               uix = tmp[!duplicated(V1), grl.ix]

                               private$pjuncs = private$pjuncs[uix]
                               return(self)
                           },

                           #' @name subset
                           #' @description 
                           #' Allows subseting of the Junction object using bracket notation
                           #' @param i integer or self$length logical vector specifying subset                        
                           subset = function(i)
                           {
                               if (is.null(i)){
                                   i = integer()
                               }

                               if (is.logical(i)){
                                   i = which(i)
                               }
                               if (length(i)>0 && (is.numeric(i) | is.integer(i))) {
                                   if (max(i, na.rm = TRUE)>self$length) {
                                       stop('index out of bounds')
                                   }
                                   
                                   if (any(i<0))
                                   {
                                       if (!all(i<0)){
                                           stop(paste('cannot mix positive ',
                                                      'with negative subscripts',
                                                      ' for Junction object'))
                                       }
                                       
                                       i = setdiff(1:self$length, abs(i))
                                   }
                               }
                               
                               private$pjuncs = private$pjuncs[i]
                               
                               return(self)
                           },
                           
                           
                           #' @name print
                           #' @description
                           #' 
                           #' Prints out the Junction Object.
                           #' Prints the length and the GRangesList of the junctions.
                           print = function()
                           {
                               message("Junction Object with ", self$length, " junctions\n")
                               if (self$length>0)
                               {
                                   HEAD = 1:pmin(4, self$length)
                                   bp = grl.unlist(self$grl[HEAD])
                                   bpstr = gr.string(bp)
                                   bpdt = data.table(str = bpstr,
                                                     ix = bp$grl.ix)[
                                     , .(junc = paste(str, collapse = ' <-> ')),
                                       keyby = ix][.(HEAD), ]

                                   if (ncol(self$dt)>0)
                                   {
                                       bpdt = cbind(bpdt[, "junc", with = FALSE],
                                                    self$dt[HEAD, ])
                                   }
                                   print(bpdt)
                                   more = self$length-HEAD[length(HEAD)]
                                   if (more>0)
                                       message('... (', more,' additional junctions)')
                               }
                           }
                       ),

                       
                       private = list(
                           ## GRangesList of junctions, each element has length 2
                           pjuncs = NULL
                       ),

                       
                       active = list(

                           #' @name grl
                           #' @description
                           #'
                           #' Returns the GRangesList of the junctions in this Junction Object
                           #'
                           #' @return GRangesList of the junctions in this Junction Object
                           grl = function()
                           {
                               return(private$pjuncs)
                           },

                           #' @name length
                           #' @description
                           #' 
                           #' Returns the number of junctions in this Junction Object.
                           #' 
                           #' @return Number of junctions in this Junction Object
                           length = function()
                           {                               
                               return(length(private$pjuncs))
                           },

                           copy = function() self$clone(),

                           #' @name dt
                           #' @description
                           #'
                           #' Returns the GRangesList of the junctions in the Junction Object as a data.table.
                           #'
                           #' @return data.table GRangesList of the junctions coverted to a data.table
                           dt = function()
                           {
                               return(as.data.table(values(private$pjuncs)))
                           },


                           #' @name breakpoints
                           #' @description
                           #'
                           #' Returns a GRanges represnting the spots a genome needs to be broken for these
                           #' junctions to be added to the genome. If there is a pad on a junctions or
                           #' there is an approximate junction, the left end of a "-" strand and the right
                           #' most edge of a "+" strand will be used. Will remove duplicate breakpoints.
                           #'
                           #' @return GRanges of breakpoints generated from the junctions in this Junction Object
                           breakpoints = function()
                           {
                               ## Get the breakpoints, left most end of - strand, right most end of + strand
                               bps = gr.end(granges(unlist(private$pjuncs)), ignore.strand = FALSE)
                               names(bps) = NULL
                               
                               ## Deal with shift at value one
                               gr = bps %Q% (strand == "+")
                               gr = GenomicRanges::shift(gr, -1)
                               
                               gr = c(bps %Q% (strand == "-"), gr)
                               gr = gr[!duplicated(gr)]
                               return(gr)
                           },
                           #' @name span
                           #' @description
                           #'
                           #' Returns the distance between breakpoint pairs on the genome
                           #'
                           #' @return vector of spans of all junctions in junction object
                           span = function()
                           {
                               grl.eval(private$pjuncs, ifelse(seqnames[1]==seqnames[2], as.numeric(abs(diff(start))), Inf))
                           },

                           #' @name flip
                           #' @description
                           #'
                           #' Flips strand of junctions in junction object, i.e.
                           #' returning the reciprocal juncion
                           #'
                           #' @return vector of spans of all junctions in junction object
                           flip = function()
                           {

                           },

                           #' @name sign
                           #' @description
                           #'
                           #' Sign is the product of strands, i.e. ++ and -- is 1
                           #' and "+-" and "-+" is -1
                           #'
                           #' @return vector of orientations of each junction in junction object
                           sign = function()
                           {
                               return(grl.eval(private$pjuncs, sign((sign(strand[1]=='+')-0.5)*(sign(strand[2]=='+')-0.5))))
                           },
                           

                           #' @name graph
                           #' @description
                           #'
                           #' Returns the gGraph this Junction Object points at.
                           #'
                           #' @return gGraph this Junction Object points at
                           graph = function()
                           {
                               return(gGraph$new(juncs = self))
                           }
                       )
                       )


## ================== Non-Member Functions for Junction ================== ##

#' @name length
#' @title length.Junction
#' @description
#' 
#' The number of junctions in this Junction Object
#'
#' @param Junction a Junction Object
#' @return the number of junctions in the Junction Object
#' @export
`length.Junction` = function(Junction)
{
    return(Junction$length)
}


#' @name c
#' @title c.Junction
#' @description
#'
#' Concatenates Junction Objects
#'
#' @param Junction object
#'
#' @return a new concatenated Junction Object
#' @author Rick Mortensen
#' @export
`c.Junction` = function(...)
{                            
    juncs.list=list(...)

    isg = sapply(juncs.list, inherits, 'Junction')
    
    if(any(!isg)){
        stop('Error: All inputs must be of class Junction.')
    }

    ## Get all the pjuncs to create new Junction Object

    grll = lapply(juncs.list, function(x) x$grl)
    grlm = lapply(grll, function(x) as.data.table(values(x)))
    newgrl = dodo.call(grl.bind, lapply(grll, function(x) {values(x) = NULL; x}))
    values(newgrl) = rbindlist(grlm, fill = TRUE)
    return (Junction$new(newgrl))
}

#' @name union
#' @title union.Junction
#' @description
#' #o'
#' Returns a new Junction Object which is the union of x and y.
#'
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the union of x and y
#'
#' @export
#' setMethod("union", c('Junction', "Junction"), function(x, y, pad = 0, ignore.strand = FALSE, ...)
#' {
#'   newJunc=c(x, y)
#'   return(unique(newJunc, pad, ignore.strand))
#' })
#'
#' @name unique
#' @title unique.Junction
#' @description
#'
#' Returns the subset of Junction object that it is unique
#'
#' @param x a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the union of x and y
#' @export
"unique.Junction" = function(x, pad = 0, ignore.strand = FALSE)
{
    if (pad==0)
        return(x$copy$removeDups())
    else
        return(x[!ra.duplicated(x$grl)])
}
setMethod("unique", c('Junction'), unique.Junction)

#' @name setdiff
#' @title setdiff.Junction
#' @description
#'
#' Returns a new Junction Object which is the difference between x and y.
#'
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the difference between x and y
#' @export
setMethod("setdiff", c('Junction', "Junction"), function(x, y, pad = 0, ...)
{
    ov = ra.overlaps(x$grl, y$grl, pad = pad)
    ix = setdiff(1:length(x$grl), ov[,1])
    return(x[ix])
})

#' @name intersect
#' @title intersect.Junction
#' @description
#'
#' Returns a new Junction Object which is the intersection of x and y.
#'
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the intersection of x and y
#' @export
setMethod("intersect", c('Junction', 'Junction'), function(x, y, pad = 0, ...) {
    ov = ra.overlaps(x$grl, y$grl, pad = pad)
    return(unique(x[ov[, 'ra1.ix']], pad = pad))
})

###############
#' @name union
#' @title union.Junction
#' @description
#' Returns a new Junction Object which is the union of x and y.
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the union of x and y
#' @export
## setMethod("union", c('Junction', "Junction"), function(x, y, pad = 0, ignore.strand = FALSE, ...)
`union.Junction` = function(x, y, pad = 0, ignore.strand = FALSE, ...)
{
    newJunc=c(x, y)
    return(unique(newJunc, pad, ignore.strand))
}
## )

#' @name unique
#' @title unique.Junction
#' @description
#' 
#' Returns the subset of Junction object that it is unique
#' 
#' @param x a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the union of x and y
#' @export
"unique.Junction" = function(x, pad = 0, ignore.strand = FALSE)
{
    if (pad==0)
        return(x$copy$removeDups())
    else
        return(x[!ra.duplicated(x$grl)])
}
## setMethod("unique", c('Junction'), unique.Junction)


#' @name setdiff
#' @title setdiff.Junction
#' @description
#' 
#' Returns a new Junction Object which is the difference between x and y.
#'
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the difference between x and y
#' @export
setMethod("setdiff", c('Junction', "Junction"), function(x, y, pad = 0, ...)
{  
    ov = ra.overlaps(x$grl, y$grl, pad = pad)
    ix = setdiff(1:length(x$grl), ov[,1])
    return(x[ix])
})
## `setdiff.Junction` = function(x, y, pad = 0, ...)
## 
## )


#' @name intersect
#' @title intersect.Junction
#' @description
#' 
#' Returns a new Junction Object which is the intersection of x and y.
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the intersection of x and y
#' @export
`intersect.Junction` = function(x, y, pad = 0, ...) {
    ov = ra.overlaps(x$grl, y$grl, pad = pad)
    return(unique(x[ov[, 'ra1.ix']], pad = pad))
}
## setMethod("intersect", c('Junction', 'Junction'), function(x, y, pad = 0, ...) {
## )


#' @name [.Junction
#' @title Junction
#' @description
#'
#' Overloads subset operator for Junction. Allows subsetting of Junction via index.
#' Also, allows for data.table style queries on Junction GRangesList and
#' corresponding metadata.
#'
#' @param obj Junction object This is the Junction object to be subset
#' @param i integer, logical, or expression in Junction metadata used to subset Junction
#' @return A new Junction object that contains only the given id's
#' @export
'[.Junction' = function(obj, i, j){

    juncs = obj$clone()
    if (!missing(i))
    {
        inew = tryCatch(
            eval(eval(parse(text = substitute(deparse(substitute(i)))),
                      parent.frame()),juncs$dt, parent.frame(2)),
            error = function(e) NULL)
        if (is.null(inew))
            inew = i ## just give up      
        juncs$subset(inew)
    }

    if (!missing("j"))
    {
        grl = juncs$grl
        values(grl) = values(grl)[, j, drop = FALSE]
        juncs = jJ(grl)
    }
    return(juncs)
}


#' @name +.Junction
#' @description
#'
#' Allows padding of junctions. The rpad will be added to the left of a "+" junction
#' and to the right of "-" junction.
#'
#' @param jj Junction Object
#' @param pad Positive number representing amount to pad the Junction Object.
#' @return a new Junction class with the padding applied
#' @export
'+.Junction' = function(jj, pad)
{
    new.grl = resize(jj$grl, pad+1, fix="end")
    values(new.grl) = values(jj$grl)
    return(Junction$new(new.grl))
}


#' @name refresh
#' @description
#' 
#' Updates Junction object to reflect changes in source code
#' 
#' @param Junction object
#' @return Junction object
#' @export
setGeneric("refresh", function(x) standardGeneric("refresh"))
setMethod("refresh", "Junction",
          function(x) {
              return(Junction$new(x$juncs))
          })


## ================== gGraph class definition ================== ##
#' @export
gGraph = setClass("gGraph")
gGraph = R6::R6Class("gGraph",
                     public = list(
                         ## public fields

                         ## constructor INIT GGRAPH
                         #' @name gGraph constructor 
                         #' @description
                         #' All purpose constructor of gGraphs from
                         #' nodes, edges, junctions or various input formats (JaBbA, Weaver, etc)
                         #' @param genome Seqinfo or object coercible to seqinfo
                         #' @param breaks GRanges whose endpoints specify breakpoints in the genome
                         #' @param juncs Junction object or GRangesList coercible to Junction
                         #' @param jabba JaBbA graph rds file
                         #' @param cougar CouGar output directory path
                         #' @param weaver Weaver output directory path
                         #' @param remixt RemiXT output directory path 
                         #' @param prego  PREGO output directory path
                         #' @author Joe DeRose
                         initialize = function(genome = NULL,
                                               breaks = NULL,
                                               juncs = NULL,
                                               prego = NULL,
                                               jabba = NULL,
                                               cougar = NULL,
                                               weaver = NULL,
                                               remixt = NULL,
                                               nodes = NULL,
                                               edges = NULL,
                                               nodeObj = NULL,
                                               edgeObj = NULL,
                                               meta = NULL,
                                               verbose = FALSE
                                               )
                         {
                             ## control how to construct
                             meta = c(meta, list())

                             if (is.null(nodes) & is.null(edges))
                             {
                                 if (!is.null(nodeObj))
                                 {
                                     private$gGraphFromNodeObj(nodeObj, edgeObj)
                                     meta = meta[!duplicated(names(meta))]
                                     private$pmeta = meta
                                     return(self)
                                 }
                                 else if(!is.null(breaks) || !is.null(juncs))
                                 {                                                           
                                     ne = breakgraph(breaks, juncs)
                                 }
                                 else if(!is.null(prego))
                                 {
                                     ne = pr2gg(prego)
                                     meta[['y.field']] = 'cn'
                                     meta[['name']] = 'PREGO'
                                     meta[['y0']] = 0
                                     meta[['by']] = 'cn'
                                 }
                                 else if (!is.null(jabba))
                                 {
                                     ne = jab2gg(jabba)
                                     meta[['y.field']] = 'cn'
                                     meta[['name']] = 'JaBbA'
                                     meta[['y0']] = 0
                                     meta[['by']] = 'cn'
                                 }
                                 else if (!is.null(cougar))
                                 {
                                     ne = cougar2gg(cougar)
                                     ne$nodes = inferLoose(ne$nodes, ne$edges)
                                     meta[['name']] = 'CouGAR'
                                     meta[['y.field']] = 'cn'
                                     meta[['y0']] = 0
                                     meta[['by']] = 'cn'
                                 }
                                 else if (!is.null(weaver))
                                 {
                                     ne = wv2gg(weaver)
                                     ne$nodes = inferLoose(ne$nodes, ne$edges)
                                     meta[['name']] = 'Weaver'
                                     meta[['y.field']] = 'cn'
                                     meta[['y0']] = 0
                                     meta[['by']] = 'cn'
                                 }
                                 else if (!is.null(remixt))
                                 {
                                     ne = remixt2gg(remixt)
                                     ne$nodes = inferLoose(ne$nodes, ne$edges)
                                     meta[['name']] = 'ReMiXt'
                                     meta[['y.field']] = 'cn'
                                     meta[['y0']] = 0
                                     meta[['by']] = 'cn'
                                 }
                                 else
                                 {
                                     private$emptyGGraph(genome)
                                     return(self)
                                 }

                                 nodes = ne$nodes
                                 edges = ne$edges
                             }

                             private$gGraphFromNodes(nodes, edges)

                             ## set gGraph metadata
                             ## (mostly gTrack settings at this point)
                             meta = meta[!duplicated(names(meta))]
                             private$pmeta = meta

                             return(self)
                         },
                         
                         #' @name walks
                         #' @description
                         #' Exhaustively generates walks (if greedy = FALSE)
                         #' or otherwise applies a greedy heuristic (greedy = TRUE)
                         #' returns a gWalk object tied to this graph
                         #' warning: greedy = FALSE will not scale to large graphs
                         #' (i.e. may even hang on a couple of hundred nodes, depending
                         #' on the topology)
                         #' @param greedy logical scalar specifying whether to generate greedy walks
                         #' @param verbose logical scalar
                         #' @author Marcin Imielinski
                         walks = function(greedy = FALSE, verbose = FALSE)
                         {
                             if (greedy==TRUE){
                                 stop('Greedy walks not yet implemented - please stay tuned')
                             }
                             A = self$adj
                             colnames(A) = rownames(A) = self$gr$snode.id


                             lleft = which(self$nodes$loose.left)
                             lright = which(self$nodes$loose.right)

                             llefti = self$queryLookup(lleft)$index
                             lrighti = self$queryLookup(lright)$index

                             lleftir = self$queryLookup(-lleft)$index
                             lrightir = self$queryLookup(-lright)$index

                             sources = c(llefti, lrightir)
                             sinks = c(lrighti, lleftir)

                             ap = all.paths(A, sources = sources, sinks = sinks, verbose = verbose)

                             ## first make paths from any nodes that are both
                             ## sources and sinks
                             tmp = self$gr$snode.id[intersect(sources, sinks)]
                             paths = split(tmp, seq_along(tmp))

                             paths = c(paths, lapply(ap$paths, function(x) self$gr$snode.id[x]))
                             cycles = lapply(ap$cycles, function(x) self$gr$snode.id[x])

                             ## dedup reciprocals
                             ## gets tricky for cycles since may be
                             ## out of phase, however the (unordered) set of vertices / antivertices
                             ## will be unique to that path / cycle and it's reciprocal
                             ## so we don't need to worry about the sequence to match these up
                             ## can just sort and match


                             ## label each path / cycle
                             pstr = sapply(paths, function(x) paste(sort(x), collapse = ' '))
                             cstr = sapply(cycles, function(x) paste(sort(x), collapse = ' '))

                             ## label reverse complement version
                             rpstr = sapply(paths, function(x) paste(sort(-x), collapse = ' '))
                             rcstr = sapply(cycles, function(x) paste(sort(-x), collapse = ' '))


                             ## match up paths and their reverse complement
                             ped = cbind(seq_along(rpstr), match(pstr, rpstr))
                             pcl = igraph::clusters(igraph::graph.edgelist(ped), 'weak')$membership

                             ## match up cycles and their reverse complement
                             ## use vgrep here since cycles generated by all paths
                             ## may be "out of phase" with their reverse complements
                             ced = cbind(seq_along(rcstr), match(cstr, rcstr))
                             ccl = igraph::clusters(igraph::graph.edgelist(ced), 'weak')$membership

                             paths = paths[!duplicated(pcl)]
                             cycles = cycles[!duplicated(ccl)]

                             circular = c(rep(FALSE, length(paths)),
                                          rep(TRUE, length(cycles)))

                             return(gWalk$new(snode.id = c(paths, cycles), graph = self, circular = circular))
                         },                      

                         #' @name set
                         #' @description
                         #'
                         #' set metadata of gGraph
                         #' right now mainly useful for gTrack defaults
                         #' such as "name" or "colormaps", and also used
                         #' for setting default "by" field for simplify
                         #' but can be used for configuring other settings in the future
                         #' (note that this is graph level and not gNode or gEdge level
                         #' metadata)
                         #' @param ... name value pairs
                         #' @author Marcin Imielinski
                         set = function(...)
                         {
                             args = list(...)

                             for (arg in names(args)){
                                 private$pmeta[[arg]] = args[[arg]]
                             }
                         },

                         #' @name queryLookup
                         #' @description
                         #'
                         #' Returned a data.table of the provided snode.ids, their indicies and the indicies of
                         #' their reverse complements in the graph. data.table is keyed on snode.id.
                         #'
                         #' @param id snode.ids to look up
                         #' @return data.table of snode.ids, indicies and reverse complement indicies
                         #' @usage
                         #'
                         #' id = c(1,3,-4,2)
                         #' gg$queryLookup(id)
                         #' @param id signed node ids in graph
                         #' @author Joe DeRose
                         queryLookup = function(id) {
                             dt = private$lookup[.(id)]
                             return(dt)
                         },

                         #' @name disjoin
                         #' @description
                         #' disjoins (i.e. collapses) all overlapping nodes in graph (subject to "by" argument), and aggregates node and edge
                         #' metadata among them using FUN
                         #' modifies the current graph
                         #' optional input gr will first concatenate a reference graph with GRanges gr prior to disjoining
                         #' collapse argument (if TRUE) will output a graph where there is a single
                         #' node per reference interval
                         #' and if collapse = FALSE will only disjoin all the nodes in the graph but keep all overlapping nodes separate
                         #' i.e. so that overlapping graphs are composed of a common set of disjoint intervals, but we allow there 
                         #' to be several instances of a given interval among the different graphs
                         #' @param gr GRanges around which to disjoin the graph
                         #' @param by metadata field of current graph around which to limit disjoining
                         #' @param collapse logical scalar specifying whether to collapse graph nodes after disjoining
                         #' @param na.rm logical scalar specifying whether to remove NA's when aggregating metadata after collapsing
                         #' @param avg logical scalar specifying whether to average (if TRUE) or sum (if FALSE) numeric metadata during aggregation (default = FALSE)
                         #' @param FUN function which should take (numeric or character) x and na.rm = TRUE and return a scalar value
                         #' @author Marcin Imielinski
                         disjoin = function(gr = NULL,
                                            by = NULL,
                                            collapse = TRUE,
                                            na.rm = TRUE,
                                            avg = FALSE,
                                            sep = ',',
                                            FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                         {
                             this = self

                             if (length(self)==0) 
                                 return(self)

                             if (!is.null(gr)){
                                 ## adding the breaks but no edges in there?
                                 this = c(self, gG(breaks = gr))
                             }

                             ## now we disjoin all the nodes and merge them back with their parents
                             ## XT debug: this dnodes is NOT disjoint!
                             dnodes = disjoin(this$nodes$gr) %*% this$nodes$gr

                             ## nmap maps our original node ids to the new node ids
                             ## (but we throw out all "internal nodes" created by the disjoiin ## here since we will use nmap to map edges onto the new disjoined nodes
                             nmap = gr2dt(dnodes)[, .(subject.id)][, dnode.id := 1:.N]
                             nmap[, left := start(this$nodes$gr)[subject.id] == start(dnodes)]
                             nmap[, right := end(this$nodes$gr)[subject.id] == end(dnodes)]

                             ## clean up loose ends dnodes
                             ## i.e. internal nodes can't be loose ends
                             dnodes$loose.left[nmap[left == FALSE, dnode.id]] = FALSE
                             dnodes$loose.right[nmap[right == FALSE, dnode.id]] = FALSE

                             ## remove all non left  / non right (i.e. "internal") nodes
                             ## since we don't need them to remap edges
                             nmap = nmap[left | right, ]

                             ## nmap will now be a one to 1/2 mapping
                             setkey(nmap, subject.id)
                             
                             ## update edge table to new dnodes n1 and n2
                             edges = copy(this$edges$dt)

                             ## first merge n1 and n2 with nmap
                             ## converting n1 / n2 edge pairs
                             ## into candidate  dnode.id.x / dnode.id.y edge pairs

                             ## since we've already labeled nodes as left or right (or both)
                             ## then we just need to subset on the merged dnode pairs
                             ## that are abutting the correct side of the parent node

                             dedges = data.table(
                                 n1 = integer(),
                                 n1.side = character(),
                                 n2 = integer(),
                                 n2.side = character(),
                                 type = character(),
                                 sedge.id = integer())

                             if (nrow(edges)>0)
                             {
                                 dedges =  merge(edges[, .(sedge.id, n1, n1.side, n2, n2.side, type)],
                                                 nmap, by.x = 'n1', by.y = 'subject.id', allow.cartesian = TRUE)
                                 setnames(dedges, c('left', 'right'), c('n1.left', 'n1.right'))
                                 dedges = merge(dedges, nmap, by.x = 'n2', by.y = 'subject.id', allow.cartesian = TRUE)
                                 setnames(dedges, c('left', 'right'), c('n2.left', 'n2.right'))
                                 
                                 ## remove edges to internal nodes
                                 dedges = dedges[((n1.side == 'left' & n1.left) |
                                                  (n1.side == 'right' & n1.right)) &
                                                 ((n2.side == 'left' & n2.left) |
                                                  (n2.side == 'right' & n2.right)), ]
                                 
                                 dedges = dedges[, .(n1 = dnode.id.x, n2 = dnode.id.y, n1.side, n2.side, sedge.id, type)]
                             }                         

                             ## we now add new reference
                             ## adjacent edges connecting disjoined nodes
                             ## in sequence to each other 
                             tmp = gr2dt(dnodes)[, dnode.id := 1:.N][, .(dnode.id, seqnames, start, end, subject.id)]
                             setkeyv(tmp, c("subject.id", "start"))
                             missing = tmp[, .(
                                 n1 = dnode.id[-.N], 
                                 n1.side = 'right', 
                                 n2 = dnode.id[-1],
                                 n2.side = 'left'), by = subject.id][!is.na(n1) & !is.na(n2), ][, type := 'REF'][, sedge.id := NA]
                             
                             ## add these internal reference edges to new edges 
                             newedges = rbind(dedges, missing[, -1, with = FALSE])
                             
                             ## identify unique dnodes in the new graph
                             ## (+/- applying by = argument)
                             if (is.null(by))
                             {
                                 udnodes = dnodes[!duplicated(dnodes$query.id), ]
                             }
                             else
                             {
                                 tmp = gr2dt(dnodes)[, dup := duplicated(query.id), by = eval(by)]
                                 udnodes = dt2gr(tmp[dup == FALSE, ][, -ncol(tmp), with = FALSE], seqlengths = seqlengths(dnodes))
                             }

                             ## now we want to map, then project the newedges onto the the udnodes
                             ## applying the aggregation function 
                             dnodes$udnode.id = match(dnodes$query.id, udnodes$query.id)


                             ## now project
                             colsn = c('node.id', 'snode.id', 'index', 'loose.left', 'loose.right')
                             colse = c('sedge.id', 'n1', 'n1.side', 'n2', 'n2.side', 'type')
                             metacolsn = setdiff(names(values(this$gr)), colsn)
                             metacolse = setdiff(names(this$edgesdt), colse)

                             final.nodes = gr2dt(dnodes[, c('udnode.id', 'loose.left', 'loose.right')])
                             if (length(metacolsn)>0){
                                 final.nodes = cbind(final.nodes, as.data.table(values(this$gr))[dnodes$subject.id, metacolsn, with = FALSE])
                             }

                             final.edges = gr2dt(newedges[, colse, with = FALSE])
                             if (nrow(final.edges)>0)
                             {
                                 if (length(metacolse)>0)
                                 {
                                     final.edges = cbind(final.edges, edges[.(newedges$sedge.id), metacolse, with = FALSE])
                                 }
                                 ## XT debug: if original edge is empty, final.edges here won't have a key!!
                                 if (!is.null(key(edges))){
                                     final.edges[, type := ifelse(is.na(edges[.(newedges$sedge.id), type]),
                                                                  type, edges[.(newedges$sedge.id), type])]
                                 }
                             }
                             
                             if (!collapse)
                             {
                                 private$gGraphFromNodes(dt2gr(final.nodes, seqlengths = seqlengths(dnodes)), final.edges)
                                 return(self)
                             }

                             final.edges[, ":="(n1 = final.nodes$udnode.id[n1], n2 = final.nodes$udnode.id[n2])]

                             ## fix edge type

                             if (is.null(FUN)){
                                 FUN = function(x) x[1]
                             }
                             id.col = c("seqnames", "start", "end", "strand", "width", "udnode.id")

                             other.col = c(id.col, metacolsn)
                             loose.col  =  c(id.col, 'loose.left', 'loose.right')

                             ## merge loose ends separately
                             ## i.e. for gGraph we take their logical OR
                             ## in bGraph, loose.agg.fun will be replaced
                             ## by addition
                             tmp.nodes = 
                                 final.nodes[, loose.col, with = FALSE][, lapply(.SD, private$loose.agg.fun), by = .(seqnames, start, end, strand, width, udnode.id)]

                             ## merge / aggregate metadata columns
                             ## if any exist
                             if (length(metacolsn)>0)
                             {
                                 final.nodes = skrub(final.nodes)
                                 
                                 meta = final.nodes[, other.col, with = FALSE][, lapply(.SD, FUN), by = .(seqnames, start, end, strand, width, udnode.id)][, metacolsn, with = FALSE]

                                 tmp.nodes = cbind(tmp.nodes, 
                                                   meta)
                             }
                             final.nodes = tmp.nodes

                             final.edges = skrub(final.edges)

                             ## here we will aggregate metadata for "identical"
                             ## edges ... i.e. those that connect the same node / sides
                             ## here we need to dedup based on those edge pairs that are identical
                             ## but have n1 and n2 flipvped
                             ## to do this we standardized the edge notation
                             ## so that n1 < n2 

                             tmp.n1 = final.edges$n1
                             tmp.n2 = final.edges$n2
                             tmp.n1.side = final.edges$n1.side
                             tmp.n2.side = final.edges$n2.side

                             final.edges[, n1 := pmin(tmp.n1, tmp.n2)]
                             final.edges[, n2 := pmax(tmp.n1, tmp.n2)]

                             final.edges[, n1.side := ifelse(tmp.n1<tmp.n2, tmp.n1.side, tmp.n2.side)]
                             final.edges[, n2.side := ifelse(tmp.n1>=tmp.n2, tmp.n1.side, tmp.n2.side)]

                             final.edges = final.edges[, lapply(.SD, FUN), by = .(n1, n1.side, n2, n2.side)]

                             private$gGraphFromNodes(dt2gr(final.nodes, seqlengths = seqlengths(dnodes)), final.edges)

                             return(self)
                         },


                         #' @name simplify
                         #' @title simplify
                         #' @description
                         #' Simplifies gGraph by collapsing reference adjacent nodes
                         #' that lack a junction or loose end (ignore.loose = FALSE)
                         #' between them.
                         #'
                         #' Takes an optional "by" column. If by is not NULL simplify
                         #' will only collapse adjacent nodes if they share metadata
                         #' in the columns specified by "by"
                         #' @param by metadata field of current graph around which to limit simplification
                         #' @param na.rm logical scalar specifying whether to remove NA's when aggregating metadata after collapsing
                         #' @param avg logical scalar specifying whether to average (if TRUE) or sum (if FALSE) numeric metadata during aggregation (default = FALSE)
                         #' @param FUN function which should take (numeric or character) x and na.rm = TRUE and return a scalar value
                         #' @author Marcin Imielinski
                         simplify = function(by = private$pmeta$by, na.rm = TRUE, avg = TRUE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg), ignore.loose = FALSE)
                         {
                             if (length(self)==0) 
                                 return(self)

                             edges = copy(self$edgesdt)
                             nodes = self$nodes$gr
                             ## let's figure out reference adjacent dnode pairs
                             nodes.dt = gr2dt(nodes)[, c("node.id", "seqnames", "start", "end", by), with = FALSE]
                             nodes.dt[, endnext := end +1]
                             
                             refadj = merge(nodes.dt, nodes.dt, by.x = 'endnext', by.y = 'start', 
                                            allow.cartesian = TRUE)[, .(n1 = node.id.x, n1.side = 'right', n2 = node.id.y, n2.side = 'left', qsedge.id = 1:.N)]
                             refadjr = copy(refadj)
                             refadjr$qsedge.id = -refadjr$qsedge.id
                             setnames(refadjr, c('n2', 'n2.side', 'n1', 'n1.side', 'qsedge.id'))
                             refadj = rbind(refadj, refadjr)
                             
                             ## query dedges for nodes that have reference non-adjacent edges to the left or to the right
                             setkeyv(edges, c("n1", "n1.side", "n2", "n2.side"))
                             edges[refadj, ref := TRUE]
                             abed = edges[is.na(ref), ]

                             ## now mark nodes that abut a non-reference adjacency
                             nodes.dt[, ab.right := node.id %in% c(abed[n1.side == 'right', n1], abed[n2.side == 'right', n2])]
                             nodes.dt[, ab.left := node.id %in% c(abed[n1.side == 'left', n1], abed[n2.side == 'left', n2])]

                             if (!ignore.loose)
                             {
                                 nodes.dt[, ab.left := ab.left |  self$nodes$loose.left]
                                 nodes.dt[, ab.right := ab.right | self$nodes$loose.right]
                             }

                             ## create fake column to sort and group intervals on
                             by = unique(c("seqnames", by))

                             ## label groups of intervals to collapse
                             setkeyv(nodes.dt, c(by, "seqnames", "start"))

                             ## two consecutive nodes belong to different groups if (1) don't have 
                             ## adjacent coordinates and (2) there is either an aberrant junction
                             ## or loose end on the right of the first node or the left of the second
                             nodes.dt[, diff := c(end[-.N],NA)+1 != c(start[-1], NA) |
                                            (c(ab.right[-.N], NA) | c(ab.left[-1], NA)), by = eval(paste(by))]

                             ## label the interval after each diff == TRUE with a new group
                             nodes.dt[, group := cumsum(1:.N %in% (which(diff == TRUE)+1)), by = eval(paste(by))]

                             ## start a fresh dt just in case we overwrote
                             ## some user metadata columns
                             nodes.dt2 = gr2dt(nodes)
                             setkey(nodes.dt, node.id)
                             nodes.dt2$simplify.group = nodes.dt[.(nodes.dt2$node.id), group]
                             metadata.cols = setdiff(names(nodes.dt2),
                                                     c('seqnames', 'start', 'end', 'strand', 'width',
                                                       'snode.id', 'node.id', 'index', 'loose.left',
                                                       'loose.right',  'simplify.group', by))

                             ## now start to aggregate / merge using FUN
                             by = unique(c(by, 'seqnames', 'simplify.group'))

                             if (is.null(FUN)) ## if no fun specified just take the first 
                                 FUN = function(x) x[1]

                             setkeyv(nodes.dt2, c(by, 'start'))
                             nodes.out = 
                                 nodes.dt2[, .(start = start[1], end = end[.N], node.id.left = node.id[1], node.id.right = node.id[.N], loose.left = loose.left[1], loose.right = loose.right[.N]), by = eval(paste(by))]

                             if (length(metadata.cols)>0)
                             {
                                 nodes.out = cbind(nodes.out, skrub(nodes.dt2)[, lapply(.SD, FUN), .SDcols = metadata.cols, by = eval(paste(by))][, metadata.cols, with = FALSE])
                             }

                             ## keep track of which edges connect to the node.id.left
                             ## and node.id.right nodes,
                             ## every edge being removed from the graph
                             ## should be a reference adjacency
                             edges$n2.present = edges$n1.present = FALSE
                             setkeyv(edges, c('n1', 'n1.side'))
                             edges[.(nodes.out$node.id.left, rep('left', nrow(nodes.out))), n1.present := TRUE]
                             edges[.(nodes.out$node.id.right, rep('right', nrow(nodes.out))), n1.present := TRUE]

                             setkeyv(edges, c('n2', 'n2.side'))
                             edges[.(nodes.out$node.id.left, rep('left', nrow(nodes.out))), n2.present := TRUE]
                             edges[.(nodes.out$node.id.right, rep('right', nrow(nodes.out))), n2.present := TRUE]

                             
                             edges.out = edges[(n1.present == TRUE & n2.present == TRUE), ]
                             nodes.out[, node.id := 1:.N]
                             edges.out$n1.new = edges.out$n2.new = as.integer(NA)
                             edges.out[n1.side == 'left', n1.new := match(n1, nodes.out$node.id.left)]
                             edges.out[n1.side == 'right', n1.new := match(n1, nodes.out$node.id.right)]
                             edges.out[n2.side == 'left', n2.new := match(n2, nodes.out$node.id.left)]
                             edges.out[n2.side == 'right', n2.new := match(n2, nodes.out$node.id.right)]
                             edges.out[, n1 := n1.new]
                             edges.out[, n2 := n2.new]

                             ## cleanup
                             ## FIXME: what if user specified these columns as edge metadata
                             ## they will have to adjust for now
                             edges.out$n1.new = NULL
                             edges.out$n2.new = NULL
                             edges.out$n1.present = NULL
                             edges.out$n2.present = NULL

                             private$gGraphFromNodes(dt2gr(nodes.out, seqlengths = seqlengths(nodes)), edges.out)
                             return(self)

                         },

                         #' @name reduce
                         #' @title reduce
                         #' @description
                         #'
                         #' Reduces graph which is $disjoin() followed by a simplify()$
                         #' i.e. collapsing overlapping nodes, then merging
                         #' adjacent ones subject to (optional) matching on some
                         #' metadata field
                         #' @param by metadata field of current graph around which to limit reduction (i.e. disjoining and simplification)
                         #' @param na.rm logical scalar specifying whether to remove NA's when aggregating metadata after collapsing
                         #' @param avg logical scalar specifying whether to average (if TRUE) or sum (if FALSE) numeric metadata during aggregation (default = FALSE)
                         #' @param FUN function which should take (numeric or character) x and na.rm = TRUE and return a scalar value
                         #' @author Marcin Imielinski
                         reduce = function(by = private$pmeta$by, na.rm = TRUE, avg = FALSE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                         {
                             self$disjoin(by = NULL, na.rm = na.rm, avg = avg, sep = sep, FUN = FUN)
                             self$simplify(by = by, na.rm = na.rm, avg = avg, sep = sep, FUN = FUN)
                         },

                         #' @name subgraph
                         #' @description
                         #' compute subgraph within a certain distance or degree of separation 
                         #' of (all nodes) intersection given GRanges "seed" window win
                         #' @param seed GRanges around which to subgraph
                         #' @param d distance in bp around which to subgraph
                         #' @param k order in number of edges which to subgraph
                         #' @param pad positive integer scalar padding to add to seed
                         #' @param ignore.strand logical scalar specifying whether to ignore.strand
                         #' @param verbose logical scalar
                         #' @author Xiaotong Yao
                         subgraph = function(seed = si2gr(self),
                                             d=NULL,
                                             k=0,
                                             bagel=FALSE,
                                             mod = FALSE,
                                             ignore.strand=T,
                                             verbose=FALSE)
                         {
                             if (verbose){
                                 message("Get the trimmed subgraph around a given GRanges within a distance on the graph.")
                             }
                             

                             win = seed;

                             if (is.character(win))
                                 win = parse.gr(win)

                             if (ignore.strand){
                                 win = gr.stripstrand(win)
                             }

                             win = gr.fix(win, private$pnodes, drop = TRUE)

                             ## DONE: what to do when win is larger than segs?????
                             ## ans: return self
                             if (length(setdiff(streduce(private$pnodes), win))==0){
                                 return(self)
                             }

                             ## overlapping window and segs, removing loose ends
                             interGr = gr.findoverlaps(private$pnodes, win, ignore.strand=ignore.strand)                         
                             qix = interGr$query.id

                             if (!is.null(d)){
                                 ## no k, use distance
                                 if (is.null(d) | d < 0){
                                     stop("Must provide either valid k or d.")
                                 }

                                 ## blend window with segs
                                 win = gr.fix(win, private$pnodes)## fix seqinfo
                                 ss = tryCatch(c(private$pnodes[, c()],
                                                 win[, c()]), error = function(e) NULL)

                                 if (is.null(ss)){
                                     ss = grbind(c(private$pnodes[, c()],
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

                                 ## FIXME: we need a purely graph and graph seed distance-based subsetting here
                                 ## in the case of graphs with overlapping nodes  
                                 ## probably just requires a couple of lines

                                 hoodRange = streduce(out)

                                 return(self$trim(hoodRange))
                             }
                             else {
                                 ## with k, go no more k steps
                                 kNeighbors = unique(unlist(igraph::ego(self$igraph, order=k, nodes=qix)))
                                 ix = unique(private$pnodes[kNeighbors]$node.id)
                                 return(self$nodes[ix]$subgraph)
                             }
                         },

                         #' @name clusters
                         #' @description
                         #'
                         #' Marks nodes in graph with metadata field $cluster
                         #' based on one of several algorithms, selected by mode
                         #'
                         #' @param weak character scalar that can take one of the following possible values - "weak" or "strong" specifying weakly or strongly connected components, walktrap specifying cluster_walktrap community detection
                         #' @author Marcin Imielinski
                         clusters = function(mode = 'weak')
                         {
                             algos = c('weak','strong','walktrap')
                             mode = as.character(factor(mode,algos))
                             if (is.na(mode))
                             {
                                 stop(sprintf('mode argument must be one of the following possibilities', paste(allowed.levels, collapse = ',')))
                             }

                             G = self$igraph

                             if (mode %in% c("strong", "weak")){
                                 membership = igraph::clusters(G, mode)$membership
                             }
                             else if (mode== 'walktrap'){
                                 membership = igraph::cluster_walktrap(G)$membership
                             }

                             ## note that membership may be different
                             ## (for certain algorithms) for a node and its
                             ## reverse complement, so we keep track of both
                             names(membership) = self$gr$snode.id

                             ## "positive membership"
                             pmembership = membership[as.character(self$nodes$dt$node.id)]

                             ## "reverse membership" (will be different unless there is
                             ## a palindromic community / component
                             rmembership = membership[as.character(-self$nodes$dt$node.id)]

                             self$annotate('cluster', data = pmembership, id = self$nodes$dt$node.id,
                                           class = 'node')
                             
                             self$annotate('rcluster', data = rmembership, id = self$nodes$dt$node.id,
                                           class = 'node')
                         },

                         jclusters = function(d1 = 1e4,
                                              d2 = 1e6,
                                              mc.cores = 1){
                             self$edges$mark(og.eid = self$edges$dt$edge.id)
                             nodes = self$nodes$dt
                             altedges = self$edges[type == "ALT", ]
                             altg = self[, type=="ALT"]
                             bp = grl.unlist(altedges$grl)

                             bp.dt = gr2dt(bp)[
                               , .(grl.ix,
                                   grl.iix,
                                   strand,
                                   og.eid,
                                   node.id,
                                   snode.id,
                                   bp1 = as.character(bp1),
                                   bp2 = as.character(bp2))]
                             bp.dt[, bp.str := ifelse(grl.iix==1,
                                                      bp1,
                                                      bp2)]
                             bp.dt[, ":="(bp1 = NULL,
                                          bp2 = NULL)]
                             gr = self$gr
                             bp.dt[, ig.ix := match(snode.id, gr$snode.id)]

                             ## if two aberrant junctions share a node/bp
                             neighbors =
                                 unlist(adjacent_vertices(self$igraph,
                                                          bp.dt$ig.ix,
                                                          "total"))
                             next.n =
                                 data.table(
                                     ig.ix.i = as.numeric(
                                         sapply(strsplit(names(neighbors), "\\."),
                                                function(x) x[1])),
                                     ig.ix.j = neighbors)

                             bp.dt[,, by="og.eid"]

                             j.neighbors =
                                 do.call(
                                     `rbind`,
                                     mclapply(altedges$dt$og.eid,
                                              function(eid){
                                                  message(eid)
                                                  bps = bp.dt[og.eid==eid, c(ig.ix)]
                                                  locals =
                                                      unlist(
                                                          neighbors[as.character(bps)]
                                                      )
                                                  next.e = bp.dt[
                                                      node.id %in% gr[locals]$node.id,
                                                      setdiff(og.eid, eid)]
                                                  return(
                                                      data.table(ji = eid,
                                                                 jj = next.e)
                                                  )
                                              },
                                              mc.cores = mc.cores,
                                              mc.preschedule = FALSE))

                             ## 2) pair-wise distance between bp
                             ## three types of distances
                             bp.refd = gr.dist(bp, bp, ignore.strand=FALSE)
                             bp.altd = altg$dist(bp, bp, ignore.strand=TRUE)
                             bp.ggd = self$dist(bp, bp, ignore.strand=TRUE)

                             ## first, the close enough pairs of bps
                             close.ij = data.table(
                                 which(bp.refd<d1, arr.ind=T)
                             )[row!=col, .(i = row, j = col)]
                             close.ij[, refd := bp.refd[cbind(i, j)]]

                             near.ij = data.table(
                                 which(bp.altd<d2, arr.ind=T)
                             )[row!=col, .(i = row, j = col)]
                             near.ij[, altd := bp.altd]

                             prox.ij = data.table(
                                 which(bp.ggd<d2, arr.ind=T)
                             )[row!=col, .(i = row, j = col)]

                             ## build the distance matrix between junctions
                             return(self)
                         },
                         
                         #' @name eclusters
                         #' @description
                         #' Marks ALT edges belonging (quasi) reciprocal cycles 
                         #' @param juncs GRangesList of junctions
                         #' @param mc.cores parallel
                         #' @param ignore.strand usually TRUE
                         #' @return numerical vector of the same length, Inf means they r not facing each other
                         #' @author Marcin Imielinski
                         eclusters = function(thresh = 1e3,
                                              paths = TRUE,
                                              mc.cores = 1,
                                              verbose = FALSE,
                                              chunksize = 1e30)
                         {
                             altedges = self$edges[type == "ALT", ]

                             bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]

                             ix = split(1:length(bp),
                                        ceiling( runif(length(bp)) * ceiling(length(bp)/chunksize) )
                                        )
                             ixu = unlist(ix)
                             eps = 1e-9
                             ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
                             adj = Matrix::sparseMatrix(1, 1, x = FALSE, dims = rep(length(bp), 2))

                             if (verbose){
                                 message(sprintf(
                                     paste0('Computing junction graph across ',
                                            '%s ALT edges with distance threshold %s'),
                                     length(altedges), thresh))
                             }

                             browser()
                             ## matrix of (strand aware) reference distances
                             ## between breakpoint pairs
                             adj[ixu, ] =
                                 do.call(rbind,
                                         mclapply(ix,
                                                  function(iix)
                                                  {
                                                      if (verbose>1)
                                                          cat('.')
                                                      tmpm =
                                                          gr.dist(bp[iix],
                                                                  gr.flipstrand(bp),
                                                                  ignore.strand = FALSE)+eps
                                                      tmpm[is.na(tmpm)] = 0
                                                      tmpm[tmpm>thresh] = 0
                                                      tmpm = as(tmpm>0, 'Matrix')
                                                  },
                                                  mc.cores = mc.cores))

                             ## check bp pairs to see if they are actually
                             ## reference connected (ignore.strand = TRUE)
                             ## on the given graphs ...
                             ## which if we have many graphs overlapping
                             ## athe same intervals
                             ## may not actually be the case
                             ## we only check connectivity using ref edges

                             ## compute reference graph distance and
                             ## remove any bp pairs that are farther away
                             ## on the reference graph than on the
                             ## linear reference
                             ##FIX ME: can't handle when there are no reference edges
                             refg = self[, type == 'REF']
                             bpp = Matrix::which(adj!=0, arr.ind = TRUE)

                             dref = pdist(bp[bpp[,1]], bp[bpp[,2]])
                             drefg = diag(refg$dist(bp[bpp[,1]], bp[bpp[,2]]))
                             
                             ix = drefg>dref
                             if (any(ix)) 
                                 adj[bpp[ix,, drop = FALSE]] = FALSE
                             if (verbose>1)
                                 cat('\n')

                             adj = adj | t(adj) ## symmetrize
                             

                             ## bidirected graph --> skew symmetric directed graph conversion
                             ## split each junction (bp pair) into two nodes, one + and -
                             ## arbitrarily call each bp1-->bp2 junction is "+" orientation
                             ## then all odd proximities adjacent to bp1 will enter the "+"
                             ## version of that junction and exit the "-" version

                             ## new matrix will be same dimension as adj
                             ## however the nodes will represents + and -
                             ## orientation of junctions
                             ## using the foollowing conversion

                             ## i.e.
                             ## bp2 --> bp1 + +
                             ## bp2 --> bp2 + -
                             ## bp1 --> bp1 - +
                             ## bp1 --> bp2 - -

                             ## we'll use the same indices just to keep things confusing
                             junpos = bp1 = bp$grl.iix == 1
                             junneg = bp2 = bp$grl.iix == 2
                             ## clear out adj for new skew symmetric version
                             adj2 = adj & FALSE 
                             adj2[junpos, junpos] = adj[bp2, bp1]
                             adj2[junpos, junneg] = adj[bp2, bp2]
                             adj2[junneg, junpos] = adj[bp1, bp1]
                             adj2[junneg, junneg] = adj[bp1, bp2]

                             if (verbose)
                                 message(sprintf('Created basic junction graph using distance threshold of %s', thresh))

                             ## strongly connected components consists of
                             ## (possibly nested) cycles
                             cl = split(1:length(bp),
                                        igraph::clusters(graph.adjacency(adj2),
                                                         'strong')$membership)

                             ## choose only clusters with length > 1
                             cl = cl[S4Vectors::elementNROWS(cl)>1]
                             cl = cl[order(S4Vectors::elementNROWS(cl))]


                             jcl = lapply(cl, function(x) unique(sort(bp$grl.ix[x])))
                             jcls = sapply(jcl, paste, collapse = ' ')
                             ## bc same pair of junctions might show twice from either side?
                             jcl = jcl[!duplicated(jcls)]
                             adj3 = adj2
                             if (length(jcl)>0)
                             {
                                 dcl = dunlist(unname(jcl))[, listid := paste0('c', listid)]
                                 altedges[dcl$V1]$mark(ecycle = dcl$listid)
                                 altedges[dcl$V1]$mark(ecluster = dcl$listid)
                             } else {
                                 self$edges$mark(ecycle = NA, ecluster = NA)
                             }

                             if (verbose){
                                 message(sprintf(
                                     'Annotated %s junction cycles in edge field $ecycle',
                                     length(jcl)))
                             }
                             
                             if (paths)
                             {
                                 if (verbose)
                                     message('Analyzing paths')

                                 ## remove all cycles and enumerate remaining paths > 1
                                 if (length(jcl)>0)
                                 {
                                     adj3[unlist(jcl), unlist(jcl)] = FALSE
                                 }
                                 sinks = Matrix::which(Matrix::rowSums(adj3)==0)
                                 sources = Matrix::which(Matrix::colSums(adj3)==0)
                                 
                                 cl2 = split(1:length(bp), igraph::clusters(graph.adjacency(adj3), 'weak')$membership)
                                 cl2 = cl2[S4Vectors::elementNROWS(cl2)>1]
                                 
                                 if (any(ix <- S4Vectors::elementNROWS(cl2)>2))
                                 { ## only need to do this for connected components that have 3 or more junctions
                                     cl3 = do.call(c, mclapply(cl2[ix], function(x)
                                     {
                                         tmp.adj = adj3[x, x]
                                         lapply(all.paths(tmp.adj, sources = sources, sinks = sinks)$paths, function(i) x[i])
                                     }, mc.cores = mc.cores))
                                     
                                     cl2 = c(cl2[!ix], cl3)
                                 }
                                 jcl2 = lapply(cl2, function(x) unique(sort(bp$grl.ix[x])))
                                 jcls2 = sapply(jcl2, paste, collapse = ' ')
                                 jcl2 = jcl2[!duplicated(jcls2)]

                                 if (length(jcl2)>0)
                                 {
                                     dcl2 = dunlist(unname(jcl2))[, listid := paste0('p', listid)]
                                     altedges[dcl2$V1]$mark(epath = dcl2$listid)

                                     ## also mark ecluster, though they may have
                                     ## overlapping edges
                                     self$edges$mark(ecluster =
                                                         ifelse(
                                                             is.na(self$edges$dt$ecycle) &
                                                             is.na(self$edges$dt$epath), NA,
                                                             paste0(
                                                                 ifelse(is.na(self$edges$dt$ecycle),
                                                                        '',
                                                                        self$edges$dt$ecycle),
                                                                 ifelse(is.na(self$edges$dt$epath),
                                                                        '',
                                                                        self$edges$dt$epath))))
                                 }

                                 if (verbose)
                                     message(sprintf('Annotated %s paths in edge field $epath',
                                                     length(jcl2)))
                             }

                             ## XT added 9/6/18
                             ## the `ecluster` column points to "edge" between cycle and paths
                             ## we want the connected component of all cycles and paths
                             if (!all(is.element(c("ecycle", "epath"), colnames(self$edges$dt)))){
                                 return(invisible(self))
                             }
                             
                             if (self$edgesdt[, !any(is.na(ecycle) & is.na(epath))]){
                                 return(invisible(self))
                             }
                             ecyc = self$edgesdt[!is.na(ecycle),
                                                 setNames(unique(ecycle), unique(ecycle))]
                             epath = self$edgesdt[!is.na(epath),
                                                  setNames(unique(epath), unique(epath))]
                             ec = c(ecyc, epath)

                             ec.adj = Matrix(FALSE, nrow=length(ec), ncol=length(ec))
                             rownames(ec.adj) = colnames(ec.adj) = ec
                             ec.adj[
                                 self$edgesdt[!is.na(ecycle) & !is.na(epath)][
                                     !duplicated(ecluster),
                                     .(which(ec==ecycle), which(ec==epath)),
                                     by=ecluster][
                                   , rbind(cbind(V1, V2), cbind(V2, V1))]] = TRUE
                             ec.ig = igraph::graph.adjacency(ec.adj, "undirected")
                             ec.cl = components(ec.ig)
                             ec.dt = data.table(ec = names(ec.cl$membership),
                                                cl = ec.cl$membership)[order(cl)]
                             ec.dt[, cl.size := .N, by=cl]
                             setkeyv(ec.dt, "ec")

                             qix = ifelse(is.na(self$edgesdt$ecycle),
                                          self$edgesdt$epath,
                                          self$edgesdt$ecycle)
                             ec.cls = ec.dt[, setNames(cl, ec)]
                             self$annotate("ec.cl",
                                           ec.cls[qix],
                                           self$edgesdt$edge.id,
                                           "edge")
                             return(invisible(self))
                         },

                         #' @name paths
                         #' @description
                         #' Returns shortest paths from query
                         #' to subject in the form of gWalks
                         #'
                         #' Each output path is a  gWalk that connects query-subject on the genome
                         #' described by gGraph gg.  Each gWalk is  annotated by the metadata of the
                         #' corresponding query-subject GRanges pair as well as fields "altdist" and "refdist"
                         #' specifying the "alternate and "reference" gGraph distance of the
                         #' query-subject pair.  The gWalk metadata field "reldist" specifies 
                         #' the relative distance (i.e. ratio of altdist to refdist) for that walk.
                         #'
                         #' NOTE: this operation can be quite expensive for large combinations of
                         #' of query and subject, so max.dist parameter will by default only compute
                         #' paths for query-subject pairs that are less then max.dist apart (default 1MB).
                         #' That default is chosen for large queries (eg >10K on each side), however
                         #' for smaller queries (eg length <100) the user may want to set max.dist = Inf
                         #' 
                         #' @param query GRanges of "intervals of interest" eg regulatory elements
                         #' @param subject GRanges of "intervals of interest" eg genes
                         #' @param ref gGraph of the "reference genome", by default is the reference genome but can be any gGraph
                         #' @param ignore.strand whether to ignore strand of input GRanges
                         #' @param verbose logical flag
                         #' @param mc.cores how many cores (default 1)
                         #' @param max.dist maximum genomic distance to store and compute (1MB by default) should the maximum distance at which biological interactions may occur
                         #' @return gWalk object each representing a query-subject shortest path (if any exist)
                         #' @author Marcin Imielinski
                         path = function(query, subject = NULL, ref = NULL, reduce = TRUE, ignore.strand = TRUE, verbose = F, mc.cores = 1)
                         {
                             if (length(query)!=1 | length(subject)!=1)
                                 stop('$path only support a single query subject pair.  For larger queries consider using proximity function')
                             
                             proximity(self, query = query, subject = subject,
                                       ref = ref, reduce = reduce, ignore.strand = ignore.strand,
                                       verbose = verbose, mc.cores = mc.cores, max.dist = Inf)
                         },
                         
                         #' @name dist
                         #' @description
                         #' Computes a distance matrix of query and subject
                         #' intervals (in base pairs) on the gGraph
                         #' between any arbitrary pairs of granges gr1 and gr2.
                         #' @param gr1 GRanges query 
                         #' @param gr2 GRanges query (if NULL, will set to gr1)
                         #' @param include.internal logical flag whether to allow
                         #' paths that begin or end inside teh query or subject
                         #' @author Marcin Imielinski
                         dist = function(query,
                                         subject,
                                         ignore.strand = TRUE,
                                         include.internal = TRUE,
                                         verbose=FALSE)
                         {
                             ## FIXME: (may) need a strand aware graph distance 
                             if (!ignore.strand)
                                 stop('strand aware distance TBD')

                             ## first try evaluating query and subject as node indices or filters
                             queryf = tryCatch(self$nodes[query], error = function(e) NULL)

                             if (!is.null(queryf))
                             {
                                 if (missing('subject'))
                                     subjectf = queryf
                                 else
                                     subjectf = tryCatch(self$nodes[subject], error = function(e) NULL)                           

                                 if (is.null(subjectf))
                                     stop('Error in parsing query or subject filters')

                                 query = queryf$gr
                                 subject = subjectf$gr
                                 
                                 if (length(query)==0 | length(subject)==0)
                                     return(matrix()[c(), c()])

                                 query$id = 1:length(query)
                                 subject$id = 1:length(subject)

                                 query.s = query.e = query
                                 subject.s = subject.e = subject

                                 simpleg = self ## don't need to create a new simpleg
                                 grds = simpleg$nodes$gr[, 'snode.id']
                                 values(grds) = cbind(values(grds),
                                                      simpleg$queryLookup(grds$snode.id)[, .(index, rindex)])

                                 grds$is.query = abs(grds$snode.id) %in% abs(query$snode.id)
                                 grds$is.subject = abs(grds$snode.id) %in% abs(subject$node.id)

                                 query.ix = grds$index[grds$is.query]
                                 subject.ix = grds$index[grds$is.subject]
                                 
                                 ## reverse complement
                                 query.rix = grds$rindex[grds$is.query]
                                 subject.rix = grds$rindex[grds$is.subject]

                                 subject.sweep = width(subject)-1
                             }
                             else
                             {                         
                                 ## means that query and subject are GRanges
                                 ## so we have to chop up graph

                                 if (missing("subject"))
                                     subject = query
                                 
                                 if (length(query)==0 || length(subject)==0)
                                     return(matrix()[c(), c()])

                                 if (!any(query %^% self$nodes$gr) | !any(subject %^% self$nodes$gr))
                                     return(matrix(Inf, nrow = length(query), ncol = length(subject)))

                                 ## strip metadata from query and subject
                                 values(query) = NULL
                                 values(subject) = NULL
                                 
                                 ## copy current graph without metadata
                                 ## and simplify                         
                                 query$id = 1:length(query)
                                 subject$id = 1:length(subject)
                                 
                                 query.og = query
                                 subject.og = subject

                                 ## subject.adjustment                  
                                 subject.sweep = rep(0, length(subject)) 
                                 
                                 simpleg = gG(nodes = self$nodes$gr[, c()],
                                              edges = self$edgesdt[,.(n1,n2,n1.side,n2.side,type)])$simplify()                            
                                 simpleg$nodes$mark(simpleg = TRUE)

                                 ## XT debug: `simpleg` may be edge empty
                                 
                                 
                                 ## split query and subject if we want to allow
                                 ## paths that originate or end inside a node
                                 ## to constitute a "distance"
                                 if (include.internal)
                                 {
                                     query = query %*% simpleg$nodes$gr[, c()]
                                     subject = subject %*% simpleg$nodes$gr[, c()]
                                 }

                                 query.s = gr.start(query)
                                 query.e = gr.end(query)
                                 subject.s = gr.start(subject)
                                 subject.e = gr.end(subject)
                                 
                                 grd = unique(gr.stripstrand(grbind(query.s, query.e, subject.s, subject.e)))

                                 ## this will disjoin by grd without collapsing
                                 ## and only keep the nodes and edges that were previously in simpleg
                                 simpleg$disjoin(gr = grd, collapse = FALSE)
                                 simpleg = simpleg[which(simpleg==TRUE), ]

                                 grd$is.query = grd %^% grbind(query.s, query.e)
                                 grd$is.subject = grd %^% grbind(subject.s, subject.e)

                                 grds = grd %*% simpleg$nodes$gr[, 'snode.id']
                                 values(grds) = cbind(values(grds),
                                                      simpleg$queryLookup(grds$snode.id)[, .(index, rindex)])
                                 
                                 ## forward
                                 query.ix = grds$index[grds$is.query]
                                 subject.ix = grds$index[grds$is.subject]
                                 
                                 ## reverse complement
                                 query.rix = grds$rindex[grds$is.query]
                                 subject.rix = grds$rindex[grds$is.subject]
                             }

                             ## get igraph and populate edge weights
                             G = simpleg$igraph
                             edG = simpleg$sedgesdt[.(igraph::E(G)$sedge.id), ]
                             E(G)$weight = width(simpleg$gr)[edG$to]

                             ## since the edge weights are the target widths
                             ## and all our query and subject interval are width 1
                             ## then we just need to subtract one base from all the
                             ## distances (ie the number of bases separating our
                             ## query and target
                             Dff = igraph::shortest.paths(G, query.ix, subject.ix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             Dfr = igraph::shortest.paths(G, query.ix, subject.rix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             Drr = igraph::shortest.paths(G, query.rix, subject.rix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             Drf = igraph::shortest.paths(G, query.rix, subject.ix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             ## we pmax by zero since for edge case of self-self
                             ## node pairs we overshoot by 1
                             D = pmin(Dff, Dfr, Drr, Drf)

                             ## now we melt D and merge back with input
                             ## using overlaps of query / subject with grds
                             Dt = as.data.table(melt(D))[!is.infinite(value), ]
                             setkeyv(Dt, c("Var1", "Var2"))

                             query.ends = grbind(query.s, query.e)[, 'id'] %*% grds[, 'index']
                             subject.ends = grbind(subject.s, subject.e)[, 'id'] %*% grds[, 'index']

                             out.dt = merge(gr2dt(query.ends)[, .(id, qindex = index)],
                                            merge(gr2dt(subject.ends)[, .(id, sindex = index)],
                                                  Dt, by.x = 'sindex', by.y = 'Var2',
                                                  allow.cartesian = TRUE),
                                            by.x = 'qindex', by.y = 'Var1',
                                            allow.cartesian = TRUE)

                             setkey(out.dt, value) ## sorts by value i.e. distance
                             out.dt = unique(out.dt, by = c("id.x", "id.y"))

                             D = matrix(Inf, nrow = length(query.og), ncol = length(subject.og))
                             D[cbind(out.dt$id.x, out.dt$id.y)] = out.dt$value

                             ## sweep out subject (will be 0 if we used granges and made simpleg)
                             ## pmax only matters in case of self-self connections where we
                             ## "oversweep"
                             D = pmax(sweep(D, 2, subject.sweep, '-'), 0)

                             return(D)
                         },
                         
                         #' @name gGraph$toString
                         #' @description
                         #' Convert this gGraph to string
                         toString = function(){
                             nl = length(self$loose)
                             nt = length(self$terminal)
                             nalt = sum(private$pedges$type == "ALT")/2
                             nref = sum(private$pedges$type == "REF")/2
                             ne = nalt+nref
                             str1 = sprintf(
                                 paste('gGraph with %s nodes, ',
                                       '%s loose ends (%s terminal and %s internal),',
                                       ' and %s edges (%s REF and %s ALT)'),
                                 self$length,
                                 nl, nt, nl-nt,
                                 ne, nref, nalt)
                             ## node str
                             str2 = ifelse(length(self$nodes)>0,
                                           '\n  comprising\n',
                                           '\n')
                             str2 = paste(str2,
                                          paste(
                                              utils::capture.output(
                                                  print(self$nodes)
                                              ), collapse="\n")
                                          )
                             ## edge str
                             str3 = ifelse(length(self$edges)>0,
                                           paste(
                                               utils::capture.output(
                                                   print(self$edges)),
                                               collapse="\n"),
                                           '')
                             str = paste(str1, str2, "\n", str3, collapse="\n")
                             return(str)
                         },
                         
                         #' @name gGraph$print
                         #' @description
                         #'
                         #' Prints out this gGraph. Prints number of nodes and edges, the gNode associated
                         #' with this gGraph and the gEdge associated with this gGraph
                         print = function()
                         {
                             message(self$toString(),
                                     appendLF = TRUE)

                             ## if (self$nodes$length > 0){
                             ##     message(' comprising:\n')
                             ## }
                             ## else
                             ##     message('\n')
                             ## self$nodes$print()
                             
                             ## if (nrow(private$pedges))
                             ## {
                             ##     message()
                             ##     self$edges$print()
                             ## }
                         },


                         #' @name annotate
                         #' @description
                         #'
                         #' Used by the mark() functions in gNode, gEdge and gWalks to alter the metadata
                         #' associated with the nodes and edges in this gGraph. Not recommended to use this
                         #' function. It is much safer to use mark.
                         #'
                         #' FYI
                         #' id for nodes is the node.id (not snode.id)
                         #' id for edges is the edge.id (not sedge.id)
                         annotate = function(colName, data, id, class)
                         {
                             ## TODO: if id not given, why not label all the elements
                             if (class == "node") {
                                 NONO.FIELDS = c('node.id', 'snode.id', 'index', 'loose.left', 'loose.right', 'loose.left', 'loose.right')
                                 if (colName %in% NONO.FIELDS)
                                     stop(paste('Cannot alter these protected gNode fields: ', paste(NONO.FIELDS, collapse = ', ')))

                                 if (is.null(data))
                                 {
                                     values(private$pnodes)[[colName]] = NULL
                                     return(self)
                                 }

                                 id = self$queryLookup(id)
                                 id$data = data

                                 index = c(id[, index], id[, rindex])
                                 data = c(id[, data], id[, data])
                                 
                                 gr.dt = gr2dt(private$pnodes)
                                 gr.dt[index, paste(colName) := data]
                                 values(private$pnodes)[[colName]] = gr.dt[[colName]]
                                 
                             } else if (class == "edge") {
                                 NONO.FIELDS = c('from', 'to', 'sedge.id', 'edge.id', 'type', 'n1', 'n2', 'n1.side', 'n2.side')
                                 if (colName %in% NONO.FIELDS)
                                     stop(paste('Cannot alter these protected gEdge fields: ', paste(NONO.FIELDS, collapse = ', ')))
                                 
                                 if (is.null(data))
                                 {
                                     private$pedges[[colName]] = NULL
                                     return(self)
                                 }

                                 id2 = id ## for some reason data.table needs this to assign correctly (SCARY!!) ... must be a substitution / promise / NSE issue
                                 private$pedges[.(id2), paste(colName) := data]                           
                             } else {
                                 stop("Not sure how we got to this error at all, we should never be here")
                             }
                             return (self)
                         },


                         #' @name window
                         #' @description
                         #'
                         #' Returns the region this gGraph spans as a GRanges
                         #'
                         #' @param pad A positive amount to pad the window by
                         #' @return GRanges of the region this gGraph covers
                         window = function(pad = 0)
                         {
                             return(streduce(gr.stripstrand(private$pnodes + pad)))
                         },


                         gtrack = function(y.field = NULL, ...)
                         {
                             ss = private$pnodes
                             ed = private$pedges

                             if (is.null(ss$loose))
                             {
                                 lleft = self$nodes$lleft[, c('snode.id', 'node.id')]
                                 lright = self$nodes$lright[, c('snode.id', 'node.id')]

                                 if ((length(lleft) + length(lright))>0)
                                 {
                                     last = length(ss)
                                     lleft$index = seq_along(lleft) + last

                                     last = last + length(lleft)
                                     lright$index = seq_along(lright) + last

                                     loose.ed = rbind(
                                         data.table(from = lleft$index,
                                                    to = self$queryLookup(lleft$snode.id)$index),
                                         data.table(from = self$queryLookup(lright$snode.id)$index,
                                                    to = lright$index)
                                     )[, type := 'loose']

                                     loose =  gr.flipstrand(grbind(lleft, lright))
                                     loose$loose = TRUE
                                     ss$loose = FALSE
                                     ss = grbind(ss, loose)
                                     ed = rbind(ed, loose.ed, fill = TRUE)
                                 }
                             }


                             if (!is.null(ed))
                             {

                                 ## set edge apperances
                                 ## lwd, lty, col, cex.arrow, v, not.flat, h, dangle.w
                                 if (is.null(y.field) || !is.element(y.field, colnames(ed))) {
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

                                 ed[, ":="(lwd = ifelse(is.na(lwd), ifelse(type=="ALT", log2(0.2*pmax(0, y, na.rm = TRUE)+2)+1, ifelse(type == 'loose', 1.5, 1)),lwd),
                                           lty = ifelse(is.na(lty), ifelse(type=='loose', 1, 1),lty),
                                           col = ifelse(is.na(col), ifelse(is.na(col), ifelse(type=="ALT",
                                                                                       ifelse(y>0,
                                                                                              alpha("red", 0.4),
                                                                                              alpha("purple", 0.3)),
                                                                                       ifelse(type=="loose",
                                                                                              alpha("blue",0.6),
                                                                                              alpha("grey",0.5))), col), col),
                                           cex.arrow = ifelse(is.na(cex.arrow), 0, cex.arrow),
                                           not.flat = ifelse(is.na(not.flat), type=="ALT", not.flat),
                                           v = ifelse(is.na(v), ifelse(type=="ALT", 2, 1),v),
                                           h = ifelse(is.na(h), ifelse(type=="ALT", 2, 1),h),
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
                             if (!("col" %in% names(values(ss)))){
                                 ss$col = as.character(NA)
                             }

                             ## col, border, ywid
                             if (!("border" %in% names(values(ss)))){
                                 ss$border = as.character(NA)
                             }
                             ss$border = as.character(ss$border)
                             ss$col = as.character(ss$col)

                             if (!("lwd.border" %in% names(values(ss)))){
                                 ss$lwd.border = as.numeric(NA)
                             }
                             if (!("ywid" %in% names(values(ss)))){
                                 ss$ywid = as.numeric(NA)
                             }
                             ss$col = ifelse(is.na(ss$col), ifelse(ss$loose, alpha("white", 0), alpha("grey", 0.5)), ss$col)
                             ss$border = ifelse(is.na(ss$border), ifelse(ss$loose, ss$col, alpha("black", 0.5)), ss$border)
                             ss$lwd.border = ifelse(is.na(ss$lwd.border), 1, ss$lwd.border)
                             ss$ywid = ifelse(is.na(ss$ywid), ifelse(ss$loose, 0.001, 0.8), ss$ywid)
                             gt.args = private$pmeta[intersect(c('colormap',
                                                                 'gr.colorfield',
                                                                 'y.field',
                                                                 'col',
                                                                 'border',
                                                                 'name',
                                                                 "height", 
                                                                 'y1',
                                                                 'y0', 
                                                                 'lwd.border',
                                                                 'labels.suppress',
                                                                 'labels.suppress.gr',
                                                                 'label.suppress.grl',
                                                                 'yaxis'), names(private$pmeta))]
                             
                             args = list(...)
                             for (arg in names(args))
                             {
                                 gt.args[[arg]] = args[[arg]]
                             }

                             if (!is.null(y.field)){
                                 gt.args[['y.field']] = y.field
                             }
                             stack.gap = gt.args$stack.gap
                             y.field = gt.args$y.field

                             if (is.null(stack.gap)){
                                 stack.gap = 1e5
                             }
                             if (is.null(gt.args$angle)){
                                 gt.args[['angle']] = 0
                             }

                             if (is.null(gt.args$ylab)){
                                 gt.args[['ylab']] = y.field
                             }

                             if (is.null(gt.args$col)){
                                 gt.args$col = NA                             
                             }
                             gt.args$edges = ed

                             if (is.null(gt.args$name)){
                                 gt.args$name = 'gGraph'
                             }
                             if (!is.null(y.field) && y.field %in% colnames(values(ss))){
                                 ss$y = values(ss)[, y.field]
                                 gt.args[['data']] = unname(ss)                               
                                 gt = do.call(gTrack::gTrack, gt.args)
                             } else {
                                 gt.args$y.field = "y"
                                 gt.args$yaxis = FALSE
                                 ## stack node pairs via stack.gap
                                 tmp.ss = ss[ss$snode.id>0]
                                 tmp.ss$y = suppressWarnings(disjointBins(tmp.ss+stack.gap))
                                 ss$y = tmp.ss$y[match(ss$node.id, tmp.ss$node.id)] 

                                 ## let loose ends inherit their y value from their parent nodes
                                 if (any(lix <- ss$loose))
                                     ss$y[lix] = ss$y[ss$node.id[lix]] +
                                         sign(ss$loose[lix])*0.3

                                 ## stack parent graphs
                                 ## if the data here is the result of a
                                 ## concatenation
                                 if (!is.null(ss$parent.graph))
                                 {
                                     ## fill in any blanks using igraph
                                     if(any(is.na(ss$parent.graph)))
                                     {
                                         cl = igraph::clusters(self$igraph, 'weak')$membership
                                         tmp = data.table(pg = ss$parent.graph,
                                                          cl = cl[ss$node.id])

                                         ## use weak components to impute parent.graph
                                         pg.map = tmp[!is.na(pg), pg[1], keyby = cl]
                                         ss$parent.graph = pg.map[.(tmp$cl), V1]
                                     }

                                     mx = gr2dt(ss)[, max(y), keyby = .(seqnames, parent.graph)]
                                     mx[, offset := c(0, cumsum(V1)[-.N]), by = .(seqnames)]

                                     setkeyv(mx, c('seqnames', 'parent.graph'))
                                     ss$y = ss$y + mx[.(as.character(seqnames(ss)), ss$parent.graph), offset]
                                 }

                                 gt.args[['data']] = unname(ss)
                                 gt = do.call(gTrack::gTrack, gt.args)
                             }
                             return(gt)
                         },
                         
                         
                         #' @name trim
                         #' @description
                         #'
                         #' Trims the current gGraph to the provided GRanges and returns this as a new
                         #' gGraph.
                         #'
                         #' @param tile GRanges to trim on
                         #' @param mod Defaults to FALSE, set to TRUE to modify this gGraph
                         #' @return new gGraph trimmed to tile, unless mod is set to TRUE
                         #' @usage
                         #'
                         #' gr = c(GRanges("1", IRanges(10000,100000), "+"), GRanges("2", IRanges(10000,100000), "+"))
                         #' new.gg = gg$trim(gr)
                         #' @param tile interval around which to trim the gGraph
                         #' @author Joe DeRose
                         trim = function(tile)
                         {
                             
                             ## Some quick catch cases
                             if(length(tile) == 0) {
                                 stop("tile cannot contain no nodes, nothing to trim around")
                             }
                             if(length(self) == 0) {
                                 return(self)
                             }

                             if (is.character(tile))
                                 tile = parse.gr(tile)
                             
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
                                 es = private$pedges
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

                                 es = new.es
                             }


                             ## Remove miscellaneous metatcols added in this function
                             new.nodes$left = NULL
                             new.nodes$right = NULL
                             new.nodes$query.id = NULL
                             new.nodes$subject.id = NULL

                             pmeta = copy(private$pmeta);
                             private$gGraphFromNodes(nodes = new.nodes,
                                                     edges = es)
                             private$pmeta = pmeta
                             return(self)
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

                       #' adds a GRangesList of junctions or a Junction Object to this object via breakpoints
                       #' If the desired junctions to add do not overlap with our graph, do not add them
                       #' juncs - GRangesList() or Junction Object of junctions to add
                       #' mod - TRUE if we to change this graph
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
                                                  type = "ALT"
                                                  ## cn = if ("cn" %in% names(edges)) NA
                                                  )
                           
                           edges = rbind(edges, new.edges)
                           
                           if (mod) {
                               private$gGraphFromNodes(nodes = nodes, edges = edges)
                               return(self)
                           } else {
                               return(gGraph$new(nodes = nodes,
                                                 edges = edges))
                           }
                       },


                       #' @name json
                       #' @description 
                       #' Creates a json file for active visualization using gGnome.js
                       #' @author Marcin Imielinski
                       json = function(filename='.',
                                       maxcn=100,
                                       maxweight=100,
                                       save = TRUE,
                                       verbose = FALSE,
                                       seqlevels = c(1:22, 'X', 'Y'),
                                       settings = list(y_axis = list(title = "copy number",
                                                                     visible = TRUE)),
                                       no.y = FALSE)
                       {
                           ## Make sure the our nodes are not empty before visualizing
                           if (is.null(private$pnodes) || length(private$pnodes) == 0) {
                               stop("Cannot have empty graph and visualize it")
                           }

                           if (save){
                               if (file.exists(filename) && (file.info(filename)$isdir))
                               {
                                   basedir = filename
                                   filename = paste(basedir, "data.json", sep = '/')
                               } else
                               {
                                   if (grepl('\\.js(on)*$', filename)){                              
                                       basedir = dirname(filename)
                                   } else
                                   {
                                       stop('filename must be directory in which data.json will be created or .json filepath')
                                   }
                               }
                               
                               if (!file.exists(basedir)) {
                                   message('Creating directory ', basedir)
                                   system(paste('mkdir -p', basedir))
                               }
                           }

                           if (is.null(seqlevels))
                               seqlevels = seqlevels(seqinfo(self))

                           ## range of CN
                           ymin=0
                           ymax=maxcn
                           
                           node.json = gr2dt(self$nodes$gr[, "snode.id"])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, y = 1)]
                           
                           ed = copy(private$pedges)[sedge.id>0, .(sedge.id, from, to, type)] ## otherwise change by reference!
                           ed$from = private$pnodes$snode.id[ed$from]
                           ed$to = -private$pnodes$snode.id[ed$to]

                           yf = NULL
                           if (!no.y && !is.null(yf <- self$meta$y.field) && yf %in% names(values(self$nodes$gr)))
                           {
                               settings$y.axis = yf
                               node.json[['y']] = pmin(pmax(ymin, values(self$nodes$gr)[, yf]), ymax)
                               node.json[['y']] = ifelse(is.na(node.json[['y']]), 0, node.json[['y']])

                               if (yf %in% names(private$pedges))
                                   ed$weight = private$pedges[.(ed$sedge.id), yf, with = FALSE]
                           }
                           else
                           {
                               no.y = TRUE
                           }

                           ## make loose end nodes
                           if (length(self$loose)>0)
                           {
                               lleft = self$nodes$lleft[, c('snode.id')]
                               lright = self$nodes$lright[, c('snode.id')]

                               last = length(self)
                               lleft$index = seq_along(lleft) + last
                               
                               last = last + length(self)
                               lright$index = seq_along(lright) + last
                               
                               lleft$weight = lright$weight = as.numeric(NA)
                               if (!is.null(yf))
                               {
                                   values(lleft)$weight = values(self$nodes[lleft$snode.id]$gr)[[yf]]
                                   values(lright)$weight = values(self$nodes[lright$snode.id]$gr)[[yf]]
                               } 

                               loose.ed = rbind(
                                   data.table(from = -lleft$snode.id,
                                              to = NA, weight = lleft$weight),
                                   data.table(from = lright$snode.id,
                                              to = NA, weight = lright$weight)
                               )[, type := 'LOOSE']

                               ## remove zero or NA weight loose ends
                               loose.ed = loose.ed[!is.na(weight), ][weight>0, ]
                               loose.ed[, sedge.id := 1:.N + nrow(ed)]

                               ed = rbind(ed, loose.ed, fill = TRUE)
                           }

                           ## remove nodes and edges that are off the map
                           good.nodes = node.json[chromosome %in% seqlevels, iid]
                           node.json = node.json[iid %in% good.nodes, ]
                           ed[, good.count := rowSums(cbind(abs(from) %in% good.nodes, abs(to) %in% good.nodes), na.rm = TRUE)]
                           ed[, node.count := rowSums(cbind(!is.na(from), !is.na(to)))]
                           ed = ed[good.count == node.count, ]


                           ## TODO: do not assume things are paired up
                           ## do not assume the cn field in the segs is correct
                           node.json$title = as.character(node.json$iid)
                           node.json[, type := "interval"]
                           node.json[, strand := "*"]

                           ## if any edge left, process
                           if (nrow(ed)>0){
                               ## EDGE.JSON
                               ed[is.na(weight), weight := 0]
                               ed.json = ed[, 
                                            .(cid = sedge.id,
                                              source = from,
                                              sink = to,
                                              title = as.character(sedge.id),
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

                           gg.js = list(intervals = node.json, connections = ed.json)

                           if (no.y){
                               settings$y_axis = list(visible=FALSE)
                           }

                           if (!is.null(settings)){
                               gg.js = c(list(settings = settings), gg.js)
                           }

                           if (save){
                               if (verbose)
                               {
                                   message("Saving JSON to: ", filename)
                               }
                               jsonlite::write_json(gg.js, filename,
                                                    pretty=TRUE, auto_unbox=TRUE, digits=4)
                               return(normalizePath(filename))
                           } else {
                               return(gg.js)
                           }
                       },

                       plot = function(){
                           plot(self$gtrack(), self$footprint)
                       }

                       

                       ),


                     private = list( #### PRIVATE GGRAPH
                         ## ----- private fields
                         ## ===== required
                         ## node/vertex, a GRanges obj of strand-specific ranges
                         ## pnodes GRanges has these reserved metadata fields:
                         ## $index $loose.left $loose.right $node.id, $snode.id
                         pnodes = NULL,

                         ## data.table of all edges in g, from, to, cn, type
                         ## type can be REF, ALT
                         ## pedges data.table has these reserved columns
                         ## $from $to sedge.id $type $edge.id
                         pedges = NULL,
                         
                         ## Lookup table - must be reset with buildLookupTable
                         lookup = NULL,
                         
                         ## ===== optional slots
                         ## ALERT: whenever segs or es changes, all of these need to be reset!
                         ## igraph obj representing the graph structure DEPRECATED
                         ## pgraph = NULL,
                         
                         ## list of metadata including gTrack default settings
                         ## (e.g. colormaps) which can be saved with the object
                         pmeta = list(),

                         loose.agg.fun = function(x) as.logical(sum(x)),

                         ## ----- private methods                      

                         #' @name buildLookupTable
                         #' Builds the lookup table for this gGraph which is stored in private$lookup
                         #' Lookup table contains the columns:
                         #'     snode.id - signed node id of the nodes
                         #'     index - index corresponding to the input snode.id
                         #'     rindex - index corresponding to the complement snode.id
                         #' @param genome Optional seqinfo if the user wants to
                         #' @author Joe DeRose 
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

                         #' @name emptyGGraph
                         #' @brief Constructor, initializes an empty gGraph object. If the user does not provide genome
                         #'        and using to update this class, this will try to inherit the current gGraph's seqinfo.
                         #'        See class documentation for usage.
                         #' @param genome Optional seqinfo if the user wants to
                         #' @author Joe DeRose 
                         emptyGGraph = function(genome=NULL)
                         {
                             ## If there is an old private segs and it has seqinfo, inherit that seqinfo
                             ## If the user provides genome, skip this step
                             if (is.null(genome)) {

                                 if(!is.null(private$pnodes) && length(seqinfo(private$pnodes)@seqlengths) > 0) {
                                     genome = seqinfo(private$pnodes)
                                 }
                             }

                             if (inherits(genome, 'GRanges') | inherits(genome, 'GRangesList'))
                                 genome = seqinfo(genome)

                             if (!is.null(genome) && !inherits(genome, 'Seqinfo'))
                                 genome = Seqinfo(seqnames = names(genome), seqlengths = genome)

                             ## Set up the private fields to be empty
                             private$pnodes = GRanges(seqinfo = genome)
                             private$pedges = data.table(from=integer(0),
                                                         to=integer(0),
                                                         type=character(0))

                             return(self)
                         },
                         

                         #' @name gGraphFromNodes
                         #' @description
                         #'
                         #' This is the "master" gGraph constructor
                         #' all other constructors funnel into this one. 
                         #'
                         #' @param nodes granges of intervals, <sign is IGNORED>, optional metadata column $loose.left and $loose.right is a logical argument specifying whether the given side of the interval is a "loose end".  (Note: all terminal node "sides" will be automatically labeled as loose ends. 
                         #' @param edges optional
                         #' data.table with required fields n1, n1.side, n2, n2.side
                         #' representing node 1 and 2 index and side, where side = 0
                         #' is left side, side = 1 is right
                         #'
                         #' sets private fields $pnodes and $pedges field of gGraph
                         #' object.
                         #' @param nodes GRanges, strand is ignored
                         #' @param edges data.table with fields n1, n2, n1.side, and n2.side (see description for details)
                         #' @author Joe DeRose, Marcin Imielinski
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

                             if (any(is.na(seqlengths(nodes)))){
                                 nodes = gUtils::gr.fix(nodes)
                             }
                             if (is.null(nodes$loose.left)){
                                 nodes$loose.left = FALSE
                             }
                             if (is.null(nodes$loose.right)){
                                 nodes$loose.right = FALSE                         
                             }
                             if (nrow(edges)>0)
                             {
                                 if (is.character(edges$n1.side)) ## convert to numeric
                                 {
                                     edges[, n1.side := sign(n1.side == 'right')]
                                     edges[, n2.side := sign(n2.side == 'right')]
                                 }
                                 nodes$loose.left = nodes$loose.left | !(1:length(nodes) %in% union(edges$n1[edges$n1.side==0], edges$n2[edges$n2.side==0]))
                                 nodes$loose.right = nodes$loose.right | !(1:length(nodes) %in% union(edges$n1[edges$n1.side==1], edges$n2[edges$n2.side==1]))
                             } else
                             {
                                 nodes$loose.left = TRUE
                                 nodes$loose.right = TRUE
                             }

                             strand(nodes) = '+'
                             nodes$node.id = 1:length(nodes) ## for reverse compatibility, to get rid
                             names(nodes) = NULL
                             segs = c(nodes, gr.flipstrand(nodes))
                             segs$snode.id = ifelse(as.logical(strand(segs)=='+'), 1, -1)*segs$node.id
                             segs$index = 1:length(segs)
                             names(segs) = segs$snode.id

                             segs$loose.left = ifelse(is.na(segs$loose.left), FALSE, segs$loose.left)
                             segs$loose.right = ifelse(is.na(segs$loose.right), FALSE, segs$loose.right)

                             private$pnodes = segs
                             
                             if (nrow(edges)>0)
                             {
                                 ## FIXME: Need to change this to make sure it labels edges
                                 if(!"type" %in% names(edges)) {
                                     edges[, type := 'ALT']
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

                                 ## if not REF then ALT
                                 private$pedges[, type := ifelse(type == 'REF', 'REF', 'ALT')]

                                 setkey(private$pedges, sedge.id)
                             }

                             private$buildLookupTable()

                             return(self)
                         },

                         #' @name gGraphFromNodeObj
                         #' @description
                         #' Operates only on positive strand, treats all nodes in NodeObj as strandless -- so if a node is marked
                         #' negative it will be treated as a duplicate positive node
                         #' Possibly remove them idk
                         #' @author Joe DeRose
                         #' @param NodeObj gNode object e.g. one obtained from an existing gGraph
                         #' @param EdgeObj gEdge object e.g. one obtained from an existing gGraph
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
                             if(is.null(EdgeObj) || EdgeObj$length == 0) {
                                 private$gGraphFromNodes(nodes = NodeObj$gr)
                                 return(self)
                             }

                             ## If we have Node Obj but no EdgeObj, throw an error
                             if(NodeObj$length == 0) {
                                 stop("Error: Edges cannot non-empty if nodes are empty")
                             }
                             
                             nodes = c(NodeObj$gr, NodeObj$clone()$flip$gr)

                             ## Validate EdgeObj, no indices not present in index column of 
                             edges = EdgeObj$sdt
                             if(!all(nix <- c(edges[,to], edges[,from]) %in% nodes$index)) {
                                 stop(paste0("Edge Object contains indices not in NodeObj. Indicies to remove are: ",
                                             paste(c(edges[,to], edges[,from])[!nix], collapse = " ")))

                             }

                             ## FIXME: this only works because of how our nodes are set up but might fail later, need to use $index within convertEdges
                             nodes = nodes[!duplicated(nodes$index), ] ## MARCIN FIX .. happens with subgraphs from node / edge neibhbors


                             ## identify any nodes that previously had edges
                             ## on left or right and make these loose on that side
                             eid = EdgeObj$dt$edge.id
                             nid = NodeObj$dt$node.id
                             lleft = NodeObj$eleft$dt[!(edge.id %in% eid),
                                                      intersect(c(n1, n2), nid)]
                             lright = NodeObj$eright$dt[!(edge.id %in% eid),
                                                        intersect(c(n1, n2), nid)]
                             
                             edges = convertEdges(nodes, edges, metacols = TRUE)
                             nodes = NodeObj$gr
                             ## assign any new loose ends caused by breaking
                             ## these nodes from their previous partners
                             nodes$loose.left = nodes$loose.left | (nodes$node.id %in% lleft)
                             nodes$loose.right = nodes$loose.right | (nodes$node.id %in% lright)
                             private$gGraphFromNodes(nodes = nodes, edges = edges)

                             return(self)
                         }                         
                     ),

                     active = list(
                         ## Returns a Node Object of the nodes in the graph
                         nodes = function()
                         {
                             if (length(self$gr)>0){
                                 tmp = private$pnodes$snode.id[private$lookup[.(1:(length(private$pnodes)/2)), index]]
                             }
                             else{
                                 tmp = NULL
                             }
                             return(gNode$new(tmp, self))
                         },


                         #' @name length
                         #' @description
                         #'
                         #' Returns the number of nodes in this gGraph (pair of nodes counts as 1 node).
                         #'
                         #' @return Number of nodes in this gGraph
                         length = function()
                         {
                             return(self$nodes$length)
                         },

                         copy = function() self$clone(),
                         
                         ## Returns an Edge Object of the edges in the graph
                         edges = function()
                         {
                             if (nrow(self$sedgesdt)>0){
                                 tmp = private$pedges[.(1:(nrow(private$pedges)/2)), sedge.id]
                             }
                             else{
                                 tmp = NULL
                             }
                             return(gEdge$new(tmp, self))
                         },
                         
                         ## Returns all the loose nodes in the graph
                         ## calls gNode $loose which uses $loose.left and $loose.right
                         ## $pnodes metadata to determine the nodes that are loose on
                         ## the left and loose on the right
                         loose = function()
                         {
                             return(self$nodes$loose)
                         },

                         ## returns all terminal nodes in the graph
                         ## (this is a subset of loose ends that have either
                         ## 0 left or right degree)
                         terminal = function()
                         {
                             return(self$nodes$terminal)
                         },

                         #' returns adjacency matrix of graph as sparseMatrix 
                         adj = function(){

                             adjMat = igraph::as_adj(self$igraph)

                             ## if (is.element("cn", colnames(private$pedges))){
                             ##   adjMat[private$pedges[,cbind(from, to)]] = private$pedges$cn
                             ## }
                             return(adjMat)
                         },                                              

                         meta = function()
                         {
                             return(private$pmeta)
                         },


                         ## Returns the igraph associated with this graph
                         igraph = function()
                         {
                             if (length(self$edges)>0){
                                 ed = as.data.frame(self$sedgesdt)[
                                   , c("from", "to", "sedge.id", "edge.id", "type")]
                             } else {
                                 ed = data.frame(from = numeric(0),
                                                 to = numeric(0),
                                                 sedge.id = numeric(0),
                                                 edge.id = numeric(0),
                                                 type = character(0))
                             }
                             v = as.data.frame(
                                 values(self$gr)[, c('index', 'node.id', 'snode.id')])
                             return(igraph::graph_from_data_frame(ed,
                                                                  vertices = v,
                                                                  directed = TRUE))
                         },
                         
                         ## Returns all of the nodes in the graph as a GRanges
                         gr = function() {
                             return(private$pnodes)
                         },                                                   

                         
                         ## Returns all of the nodes in the graph as a data.table
                         dt = function() {
                             return(as.data.table(private$pnodes))
                         },
                         


                         ## Returns an Junction Object of the edges in the graph
                         junctions = function()
                         {
                             return(self$edges$junctions)
                         },

                         ## Returns all the edges in the graph as a data.table
                         sedgesdt = function() {
                             return(copy(private$pedges))
                         },

                         ## Returns all the edges in the graph as a data.table
                         edgesdt = function() {
                             sides = c('left', 'right')
                             return(copy(
                                 convertEdges(self$gr,
                                              private$pedges,
                                              metacols = TRUE)[
                                   , n1.side := sides[n1.side+1]][
                                   , n2.side := sides[n2.side+1]]
                             ))
                         },

                         
                         ## Returns a gTrack
                         gt = function()
                         {
                             stack.gap = private$pmeta$stack.gap
                             if (is.null(stack.gap)){
                                 stack.gap = 1e5
                             }
                             self$gtrack(y.field = private$pmeta$y.field, stack.gap = stack.gap)
                         },


                         ## Returns the footprint of this graph
                         footprint = function()
                         {
                             return(self$window())
                         }
                     )
                     )


#' @name c.gGraph
#' @title c.gGraph
#' @description
#'
#' Concatenates gGraphs without doing any merging or aggregating
#'
#' @param ... set of gGraph arguments (names will be incorporated as graph metadata $parent.graph)
#' @return A new gGraph object that is the union of the nodes and edges in the input gGraphs
#' @author Marcin Imielinski
#' @export
'c.gGraph' = function(...){
    args = list(...)
    args = args[sapply(args, inherits, 'gGraph')]
    if (is.null(names(args))){
        names(args) = paste0('gGraph', 1:length(args))
    }
    all.nodes = lapply(names(args), function(nm)
    {
        gr = args[[nm]]$nodes$gr
        gr$parent.graph = nm
        values(gr) = values(gr)[, setdiff(names(values(gr)), c('snode.id', 'index', 'node.id'))]
        return(gr)
    })

    ## need to offset node and edge ids in edge vector so they don't collide in the union graph 
    ns = elementNROWS(all.nodes)
    noffset = c(0, cumsum(ns))[1:length(args)]

    all.edges = mapply(function(nm, noffset)
    {
        gr = args[[nm]]$gr
        edt = args[[nm]]$sedgesdt
        edt$parent.graph = nm
        edges = convertEdges(gr, edt, metacols = TRUE)
        edges$n1 = edges$n1 + noffset
        edges$n2 = edges$n2 + noffset
        return(edges)
    }, names(args), noffset, SIMPLIFY = FALSE)

    meta = do.call('c', lapply(names(args), function(nm) args[[nm]]$meta))
    meta = meta[!duplicated(names(meta))]

    nodes = do.call(grbind, all.nodes)
    edges = rbindlist(all.edges, fill = TRUE)

    return(gGraph$new(nodes = nodes, edges = edges, meta = meta))
}


#' @name gG
#' @title create gGraph
#' @description
#'
#' Wrapper that instantiates a gGraph object from a variety of different inputs.
#'
#' Only select parameters combos need to be / should be specified simultaneously.
#'
#' @examples
#' mygenome = c('chr1' = 1e6, 'chr2' = 1e7)
#' myjuncs = system.file('extdata', "delly.final.vcf.gz", package = "gGnome") ## junctions
#' mybreaks = gr.tile(mygenome, 10000) ## 1kb tiles across my genome
#' myweaver = system.file('extdata', 'weaver', package='gGnome')
#' myremixt = system.file('extdata', 'remixt', package='gGnome')
#' myprego = system.file('extdata', 'intervalFile.results', package='gGnome')
#' myjab = system.file('extdata', 'jabba.simple.rds', package="gGnome")
#' 
#' gG(genome = mygenome)
#'
#' ## create a graph from junctions
#' gG(breaks = mybreaks, juncs = myjuncs)
#'
#' ## import a genome graph from popular callers
#' gG(jabba = myjab)
#' gG(weaver = myweaver)
#' gG(prego = myprego)
#' gG(remixt = myremixt)
#'
#'
#' ## hard code a graph from GRanges and edges
#' nodes = c(GRanges("1",IRanges(1,100),"*"), GRanges("1",IRanges(101,200),"*"),
#'                GRanges("1",IRanges(201,300),"*"), GRanges("1",IRanges(301,400),"*"),
#'                GRanges("1",IRanges(401,500),"*"))
#' edges = data.table(n1 = c(3,2,4,1,3),
#'                    n2 = c(3,4,2,5,4),
#'                    n1.side = c(1,1,0,0,1),
#'                    n2.side = c(0,0,0,1,0))
#' gg = gG(nodes = nodes, edges = edges, meta = list(gr.colorfield = 'type'))
#'
#' ## mostly for developer use 
#' gG(nodeObj = gg$nodes, edgeObj = gg$edges)
#' 
#' @param genome seqlengths or seqinfo object around which to build an empty gGraph
#' @param breaks  GRanges around which to build a reference gGraph
#' @param juncs BND vcf or bedpe file path, gGraph::Junctions object, or GRangesList specifying pairs of locations to reconnect (can be used in conjunction with breaks)
#' @param jabba path to JaBbA .rds output file, if specified instantiates a gGraph object from JaBbA output with y.field "cn" specifying copy number 
#' @param prego path to PREGO output file. if specified, instantiates a gGraph object from PREGO output with y.field "cn" specifying copy number
#' @param weaver path to Weaver output. if specified, instantiates a gGraph object from Weaver output with y.field "cn" specifying copy number
#' @param remixt path to RemiXT output, if specified, instantiates a gGraph object from ReMiXT output with y.field "cn" specifying copy number@param prego Instantiates a gGraph object from PREGO output with y.field "cn" specifying copy number
#' @param nodes GRanges of unsigned intervals to be rejoined in gGRaph, used in conjunction with edges argument below
#' @param edges data.table with field n1, n2, n1.side, n2.side with each row specifying an edge, i.e. which sides ("left" vs "right") of which integer node indices to connect in the gGRaph
#' @param nodeObj gNode object to create a gGraph around (similar to nodes input above), except generated via the $nodes accessor of an existing gGraph object, used in cojunction with gEdge input to create a new gGnome object from an existing one
#' @param edgeeObj gEdge object to create a gGraph around (similar to edges input above), except generated via the $edges accessor of an existing gGraph object
#' @param meta list of metadata to associate with this gGraph object, for now mostly used to populate gTrack visualization parameters
#' @return A new gGraph object that is the sum of the component gGraphs
#' @export
gG = function(genome = NULL,
              breaks = NULL,
              juncs = NULL,
              prego = NULL,
              jabba = NULL,
              cougar = NULL,
              weaver = NULL,
              remixt = NULL,
              nodes = NULL,
              edges = NULL,
              nodeObj = NULL,
              edgeObj = NULL,
              meta = NULL
              )
{
    return(gGraph$new(genome = genome,
                      breaks = breaks,
                      juncs = juncs,
                      prego = prego,
                      jabba = jabba,
                      cougar = cougar,
                      weaver = weaver,
                      remixt = remixt,
                      nodes = nodes,
                      edges = edges,
                      nodeObj = nodeObj,
                      edgeObj = edgeObj,
                      meta))
}



#' @name sum.gGraph
#' @title sum.gGraph
#' @description
#'
#' Sums gGraphs using default aggregation operation for metadata (sum for numeric and paste(..., collapse = ',') for character, and c() for list)
#'
#' @param ... set of gGraph arguments (names will be incorporated as graph metadata $parent.graph)
#' @param y gGraph
#' @return A new gGraph object that is the sum of the component gGraphs
#' @export
'sum.gGraph' = function(..., by = NULL, na.rm = TRUE){
    cgg = c(...)
    return(cgg$disjoin(by = NULL, na.rm = na.rm))
}

#' @name +.gGraph
#' @title +.gGraph
#' @description
#'
#' Adds two gGraphs using default aggregation operation for metadata (sum for numeric and paste(..., collapse = ',') for character, and c() for list)
#'
#' @param x gGraph
#' @param y gGraph
#' @return A new gGraph object that is the sum of the component gGRaphs
#' @export
'+.gGraph' = function(x,y){
    return(sum(x,y))
}

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
#' @author Marcin Imielinski
#' @export
'[.gGraph' = function(obj, i = NULL, j = NULL, ...){
    if (deparse(substitute(j)) != "NULL")
    {
        edges = obj$edges[j]
        if (deparse(substitute(i)) == "NULL"){
            nodes = edges$nodes
        }
        else
        {
            nodes = obj$nodes[i]
            nids = c(nodes$dt$index, nodes$flip$dt$index);
            eid = edges$sdt[to %in% nids & from %in% nids, sedge.id]    
            edges = edges[eid]
        }

    } else
    {
        if (deparse(substitute(i)) == "NULL"){
            nodes = obj$nodes
        }
        else{
            nodes = obj$nodes[i]
        }
        nids = c(nodes$dt$index, nodes$flip$dt$index);
        eid = nodes$edges$sdt[to %in% nids & from %in% nids, sedge.id]    
        edges = obj$edges[eid]
    }

    return(gGraph$new(nodeObj = nodes, edgeObj = edges, meta = obj$meta))
}


## ================== Non-Member Functions for gGraph ================== ##

#' @name length
#' @title length
#' @description
#' The number of nodes in the gGraph
#'
#' @param gGraph a \code{gGraph} object
#'
#' @return the number of nodes in the gGraph
#' @export
`length.gGraph` = function(gGraph){
    return(gGraph$length)
}


#' @name length
#' @title length
#' @description
#' The number of walks in the gWalk
#'
#' @param gGraph a \code{gWalk} object
#'
#' @return the number of nodes in the gWalk
#' @export
`length.gWalk` = function(gWalk){
    return(gWalk$length)
}


#' @name length
#' @title length
#' @description
#' The number of nodes in the gNode
#'
#' @param gGraph a \code{gNode} object
#'
#' @return the number of nodes in the gNode
#' @export
`length.gNode` = function(gNode){
    return(gNode$length)
}

#' @name length
#' @title length
#' @description
#' The number of edges in the gEdge
#'
#' @param gGraph a \code{gEdge} object
#'
#' @return the number of edges in the gEdge
#' @export
`length.gEdge` = function(gEdge){
    return(gEdge$length)
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
    return(c(gGraph$length, nrow(gGraph$sedgesdt)/2))
}

#' @name refresh
#' @title refresh
#' @description
#' Updates gGraph object to reflect changes in source code
#' 
#' @param gGraph object
#'
#' @return gGraph object
#' @export
setMethod("refresh", "gGraph",
          function(x) {
              return(gGraph$new(nodeObj = x$nodes,
                                edgeObj = x$edges))
          })

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
convertEdges = function(nodes, edges, metacols = FALSE, cleanup = TRUE)
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

    ## remove duplicate edges ie those that are present in both positive and negative orientations hence have
    ## abs(sedge.id) present twice

    if (cleanup)
    {
        es = es[!duplicated(abs(sedge.id)), ] 
        es = es[!is.na(n1) & !is.na(n2)]
    }

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

#' @name is.jb
#' @description
#' test if a gGraph is junction balanced
#' @param gg gGraph object to be tested
#' @return logical scalar
#' @export
is.jb = function(gg){
    ## reserve name "cn"
    if (!"cn" %in% colnames(gg$nodes$dt)){
        msg = "'cn' must be a metadata col of nodes"
        message(msg)
        return(FALSE)
    }
    if (length(gg$edges)>0 & !"cn" %in% colnames(gg$edges$dt)){
        msg = "'cn' must be a metadata col of edges"
        message(msg)
        return(FALSE)
    }
    ## test edge skew symmetric-ness
    ## judging from adj
    if (any(rc.adj(gg$adj, gg$gr$snode.id) !=
            gg$adj)){
        msg = paste("Some edge doesn't have the",
                    "same copy number as its",
                    "reverse complement edge")
        message(msg)
        return(FALSE)
    }
    in.fl = Matrix::colSums(gg$adj)
    if (any(in.fl!=gg$gr$cn)){
        msg = paste("Some node cn differs from",
                    "the sum of incoming edge cn")
        message(msg)
        return(FALSE)
    }
    out.fl = Matrix::rowSums(gg$adj)
    if (any(out.fl!=gg$gr$cn)){
        msg = paste("Some node cn differs from",
                    "the sum of outgoing edge cn")
        message(msg)
        return(FALSE)
    }
    return(TRUE)
}

## ================== bGraph class definition ================== ##
#' @export
bGraph = setClass("bGraph")
bGraph = R6::R6Class("bGraph",
                     inherit = gGraph,
                     public = list(
                         initialize = function(gg){
                             if (!inherits(gg, "gGraph")){
                                 stop("Input is not a gGraph object")
                             }
                             ## reserve name "cn"
                             if (!"cn" %in% colnames(gg$nodes$dt)){
                                 msg = "'cn' must be a metadata col of nodes"
                                 stop(msg)
                             }
                             if (length(gg$edges)>0 & !"cn" %in% colnames(gg$edges$dt)){
                                 msg = "'cn' must be a metadata col of edges"
                                 stop(msg)
                             }
                             ## test edge skew symmetric-ness
                             ## judging from adj
                             if (any(rc.adj(gg$adj, gg$gr$snode.id) !=
                                     gg$adj)){
                                 msg = paste("Some edge doesn't have the",
                                             "same copy number as its",
                                             "reverse complement edge")
                                 stop(msg)
                             }
                             in.fl = Matrix::colSums(gg$adj)
                             if (any(in.fl!=gg$gr$cn)){
                                 msg = paste("Some node cn differs from",
                                             "the sum of incoming edge cn")
                             }
                             out.fl = Matrix::rowSums(gg$adj)
                             if (any(out.fl!=gg$gr$cn)){
                                 msg = paste("Some node cn differs from",
                                             "the sum of outgoing edge cn")
                             }
                             ## good, copy over the data
                             ## super$initialize(nodes = gg$nodes$dt,
                             ##                  edges = gg$edges$d)t
                             super$initialize(nodeObj = gg$nodes,
                                              edgeObj = gg$edges)
                             return(self)
                         },
                         print = function(){
                             str = gsub("gGraph", "bGrpah", super$toString())
                             message(str)
                         }
                     ),
                     private = list(),
                     active = list())


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
                      ## circular is a logical vector of length walks, specifying which walk should be interpreted as a circular contig, in which case an edge is implied from the last interval / node / edge in the walk to the first (i.e. do not include this in the walk definition)
                      initialize = function(snode.id = NULL,
                                            sedge.id = NULL,
                                            grl = NULL,
                                            graph = NULL,
                                            meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
                                            circular = NULL,  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first
                                            disjoin = FALSE ## flag whether to collapse ie disjoin the intervals
                                            )
                      {
                        dgraph = NULL ## only used if grl and graph are both non NULL

                        if (all(is.null(snode.id), is.null(sedge.id), is.null(grl))) {
                          ## initializing empty gWalk

                          if (is.null(graph)){
                            graph = gGraph$new()
                          }
                            ## data.table of node.ids in the walk
                            private$pnode = gGraph$new()$nodes
                            
                            ## data.table of edge.ids in the walk
                            private$pedge = gGraph$new()$edges
                            
                            ## data.table of walk metadata 
                            private$pmeta = data.table()
                            
                            ## Pointer to the graph the snode.id/sedge.ids come from
                            ## If initialized with grl, this is the graph from the gWalk
                            private$pgraph = graph

                            return(self)
                        } else if (sum(c(!is.null(snode.id), !is.null(sedge.id), !is.null(grl))) > 1) {
                            stop("More than one of snode.id, sedge.id and grl cannot be non-NULL")
                        }


                        if (!is.null(grl)) {

                            ## will create the graph first using the unsigned nodes / edges constructor
                            ## doing it differently depending on whether disjoin is flagged

                            grlu = grl.unlist(grl)

                            ## we disjoin by default if a graph is already provided
                            disjoin = disjoin | !is.null(graph)

                            if (length(grlu)>0)
                            {
                                if (disjoin) 
                                {
                                    ## here we disjoin all nodes first and then build the graph                        
                                    nodes = disjoin(gr.stripstrand(grlu))
                                    tmpdt = as.data.table(grlu[, c('grl.ix', 'grl.iix')] %*% nodes)[, strand := as.character(strand(grlu)[query.id])]
                                }
                                else ## otherwise we keep all nodes separate 
                                {
                                    nodes = gr.stripstrand(grlu[, c()])
                                    tmpdt = cbind(as.data.table(grlu[, c('grl.ix', 'grl.iix')]),
                                                  data.table(query.id = 1:length(nodes), subject.id = 1:length(nodes)))
                                    
                                }

                                ## pos will take of duplicated created in the "disjoin" case where a single walk interval comprises multiple nodes
                                ## i.e. has teh same grl.ix and grl.iix but different pos
                                ## we want to make sure that these get ordered properly in the key call below
                                ## so that on negative strands the leftward segments come after the rightward in the same grl.iix
                                tmpdt[, pos := (sign(strand=='+')-0.5)*start]
                                setkeyv(tmpdt, c("grl.ix", "grl.iix", 'pos'))

                                edges = tmpdt[, .(n1 = subject.id[-.N], n2 = subject.id[-1],
                                                  n1.side = sign(strand[-.N] == '+'),
                                                  n2.side = sign(strand[-1] == '-')), by = grl.ix]

                                if (!is.null(circular))
                                {
                                    if (!is.logical(circular)){
                                        stop('circular must be a logical vector')
                                    }
                                    circular = data.table(walk.id = 1:length(grl),
                                                          circular = circular)$circular                            
                                }
                                else {
                                    circular = rep(FALSE, length(grl))
                                }

                                if (any(circular))
                                {
                                    ix = which(circular)
                                    setkeyv(tmpdt, c("grl.ix", "grl.iix"))
                                    edges = rbind(edges, tmpdt[grl.ix %in% ix, .(n1 = subject.id[.N],
                                                                                 n1.side = sign(strand[.N] == '+'),
                                                                                 n2 = subject.id[1],
                                                                                 n2.side = sign(strand[1] == '-')
                                                                                 ), by = grl.ix])
                                }

                                edges$grl.ix = NULL
                                snode.id = split(ifelse(tmpdt$strand=='+', 1, -1)*tmpdt$subject.id, tmpdt$grl.ix)[as.character(1:length(grl))]

                                if (!is.null(names(grl)))
                                    names(snode.id) = names(grl)


                                if (!is.null(graph)) ## disjoin the input graph over the provided grls
                                {
                                    dgraph = graph$clone()
                                    ## first we disjoin this graph over the intervals in ours
                                    gr = disjoin(unlist(grl))
                                    dgraph$disjoin(gr = gr)
                                }

                                edges[, type := 'REF'] ## we will let the graph determine which edges are alt
                                if (ncol(values(grl))>0)
                                {
                                    if (is.null(meta))
                                    {
                                        meta = as.data.table(values(grl))
                                    }
                                    else
                                    {
                                        meta = cbind(meta, as.data.table(values(grl)))
                                    }
                                }

                                graph = gGraph$new(nodes = nodes, edges = edges)
                            }
                            else
                            {
                                snode.id = list()
                                graph = gGraph$new(genome = grl)
                            }
                        }

                        if ((!is.null(snode.id) | !is.null(sedge.id)) & is.null(graph))
                        {
                            stop('graph must be non NULL if snode.id or sedge.id are provided')
                        }

                        if (!is.null(snode.id)) {
                            walk.names = names(snode.id)
                            names(snode.id) = seq_along(snode.id)
                            private$gWalkFromNodes(snode.id = snode.id,
                                                   graph = graph,
                                                   circular = circular,
                                                   meta
                                                   )
                        } else if (!is.null(sedge.id)) {
                            walk.names = names(sedge.id)
                            names(sedge.id) = seq_along(sedge.id)
                            private$gWalkFromEdges(sedge.id = sedge.id,
                                                   graph = graph,
                                                   circular = circular,
                                                   meta = meta)
                        }

                        if (nrow(private$pmeta)>0)
                        {
                            setkey(private$pmeta, walk.id)
                            setkey(private$pnode, walk.id)                             
                            setkey(private$pedge, walk.id)
                        }

                        if (!is.null(dgraph)) ## in this case both grl and graph were provided, so now we will try to reconcile
                        {
                            ## the trick is to label all edges of dgraph and then see
                            ## if the walks use any edges that are <not> inherited from dgraph
                            ## in which case we will break
                            dgraph$edges$mark(dgraph = TRUE)
                            self$disjoin(graph = dgraph)
                            if (any(is.na(self$graph$edgesdt$dgraph)))
                            {
                                stop('Some grl edges are not present in provided graph, please check inputs')
                            }
                            self$graph$edges$mark(dgraph = NULL)
                        }

                        ## fill in remaining metadata
                        if (self$length>0)
                        {
                            tmp = private$pnode[, .(walk.id, snode.id)]                   
                            tmp$wid = width(self$graph$nodes$gr)[abs(tmp$snode.id)]
                            if (is.null(walk.names))
                                walk.names = 1:self$length
                            private$pmeta$name = walk.names
                            private$pmeta$wid = tmp[, .(wid = sum(as.numeric(wid))), keyby = walk.id][.(private$pmeta$walk.id), wid]
                        }


                        return (self)
                      },

                      #' sets metadata of gWalk object
                      #' (accessible through $dt accessor)
                      set = function(...)
                      {
                          args = list(...)

                          NONO.FIELDS = c('walk.id', 'circular', 'sedge.id', 'snode.id', 'length', 'wid')
                          args = list(...)
                          
                          if (any(names(args) %in% NONO.FIELDS)){
                              stop(paste('Cannot modify the following reserved gWalk metadata fields:', paste(NONO.FIELDS, collapse = ', ')))                  
                          }
                          for (arg in names(args)){
                              private$pmeta[[arg]] = args[[arg]]
                          }
                      },

                      #' @name gWalk.subset
                      #' @description 
                      #' Allows for subsetting of the gWalk Object using bracket notation
                      #' @author Marcin Imielinski                       
                      subset = function(i)
                      {
                        if (is.null(i)){
                          i = integer()
                        }
                        if (is.logical(i)){
                            i = which(i)
                        }
                        if (length(i)>0 && (is.numeric(i) | is.integer(i)))
                        {
                          if (max(i, na.rm = TRUE)>self$length){
                              stop('index out of bounds')
                              }
                        }

                        if (length(i)==0)
                        {
                          private$pnode = data.table()
                          private$pedge = data.table()
                          private$pmeta = data.table()
                          return(self)
                        }

                        tmp.node =  rbind(                        
                            private$pnode[.(abs(i)), ],
                            copy(private$pnode[.(abs(i)), ])[,
                           ":="(walk.id = -walk.id, snode.id = -snode.id)][rev(1:.N), ])

                        new.edge = tmp.edge = data.table()
                        if (nrow(private$pedge)>0)
                          {
                            tmp.edge = rbind(
                              private$pedge[.(abs(i)), ],
                              copy(private$pedge[.(abs(i)), ])[, 
                                                               ":="(walk.id = -walk.id, sedge.id = -sedge.id)][rev(1:.N), ])
                            setkey(tmp.edge, walk.id)
                          }
                        tmp.meta = private$pmeta[.(abs(i)), ]

                        setkey(tmp.node, walk.id)                        

                        new.node = merge(data.table(query.id = 1:length(i), walk.id = i), tmp.node, by = 'walk.id', allow.cartesian = TRUE)
                        if (nrow(tmp.edge)>0)
                          {
                            new.edge = merge(data.table(query.id = 1:length(i), walk.id = i), tmp.edge, by = 'walk.id',, allow.cartesian = TRUE)
                          }
                        new.meta = merge(data.table(query.id = 1:length(i), walk.id = abs(i)), private$pmeta, by = 'walk.id', , allow.cartesian = TRUE)

                        new.node[, walk.id := query.id]
                        new.meta[, walk.id := query.id]

                        if (nrow(new.edge)>0)
                        {
                          new.edge[, walk.id := query.id]
                          private$pedge = new.edge[,-2]
                        }
                        else
                          private$pedge = data.table()

                        private$pnode = new.node[,-2]                                             
                        private$pmeta = new.meta[, -2]

                        setkey(private$pnode, walk.id)
                        setkey(private$pmeta, walk.id)
                        if (nrow(new.edge)>0)
                          {
                            setkey(private$pedge, walk.id)
                          }

                        return(self)
                      },

                      dts = function(ix = 1:self$length, makelists = TRUE)
                      {
                        out = private$pmeta[.(ix), ]
                        if (makelists == FALSE) ## just return raw pmeta subset 
                          return(out)
                        node.sum = private$pnode[.(ix), .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),][, -1]
                        out = cbind(out, node.sum)
                        if (nrow(private$pedge)>0)
                          out = cbind(out, private$pedge[.(ix), .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),][, -1])
                        return(out)
                      },

                      print = function()
                      {
                        if (self$length==0)
                        {
                          message('empty gWalk object')
                          return()
                        }

                        TOOBIG = 5
                        TEASER = 3
                        .tease = function(stuff, toobig = TOOBIG, teaser = TEASER)
                        {
                          if (length(stuff)>toobig){
                              paste(paste(stuff[1:teaser], collapse = ' -> '), '-> ... ->', paste(stuff[length(stuff)-TEASER:1+1], collapse = ' -> '))
                              }
                          else{
                              paste(stuff, collapse = ' -> ')
                              }
                        }

                        message('gWalk object with ', self$length, ' walks (',
                                sum(!self$circular), ' linear and ',
                                sum(self$circular), ' circular)')

                        if (self$length>0)
                        {
                          ix = 1:min(self$length, TOOBIG)
                          tmp.node = private$pnode[.(ix), .(walk.id, snode.id)]
                          tmp.node[, nix := private$pgraph$queryLookup(snode.id)$index]

                          setkey(tmp.node, walk.id)
                          cols = setdiff(colnames(private$pmeta), c('snode.id', 'sedge.id'))
                          out = cbind(self$dts(ix)[ , cols, with = FALSE],
                                      tmp.node[.(ix), .(gr = .tease(gr.string(private$pgraph$gr[nix], mb = FALSE))), keyby = walk.id][.(ix),][, -1])
 
                          print(out)
                          if (self$length>TOOBIG)
                          {
                            message('\n ... \n(', self$length-TOOBIG, ' more walks )')
                          }
                        }
                      },
                      #' @name gWalk eval
                      #' @description
                      #'
                      #' Evaluates an expression on gWalk node or edge metadata
                      #' and returns a vector of length length(self) one per walk.
                      #' 
                      #' Ideally the expression should generate a scalar value for each walk.id, though will dedup and re-order
                      #' even if it does not. 
                      #' @param x data.table style expression on node or edge metadata (will try each)
                      #' @param node data.table style expresion on node metadata
                      #' @param edge data.table style expresion on edge metadata
                      #' @author Marcin Imielinski
                      eval = function(x, node, edge)
                      {
                        if (self$length==0)
                          return(NULL)

                        out = NULL

                        if (missing('node') & !missing('x'))
                          delayedAssign("node", x)

                        if (missing('edge') & !missing('x'))
                          delayedAssign("edge", x)

                        if (!missing("node"))
                        {
                          out = tryCatch({
                            tmpdt = merge(private$pnode, private$pgraph$dt, by = 'snode.id')
                            out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(x), ", keyby = walk.id]")))
                            out = unique(out, by = "walk.id")
                            out[.(1:self$length), V1]
                          }, error = function(e) NULL)
                        }

                        if (is.null(out) & !missing('edge'))
                        {
                          out = tryCatch({
                            tmpdt = merge(private$pedge, private$pgraph$sedgesdt, by = 'sedge.id')
                            out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(x), ", keyby = walk.id]")))
                            out = unique(out, by = "walk.id")
                            out[.(1:self$length), V1]
                          }, error = function(e) NULL)
                        }

                        if (is.null(out))
                          stop('error with gWalk eval expression, check against node and edge metadata to see that all referenced fields exist')

                        return(out)
                      },

                      #' @name gWalk disjoin
                      #' @description
                      #' similar to gGraph disjoin but also updates according gWalk
                      #' and allows disjoining around a set of gRanges, essentially
                      #' adding "breaks" to the walks and allowing the ranges in those
                      #' walks to inherit metadata from overlapping GRanges
                      disjoin = function(gr = NULL, graph = NULL, by = NULL, na.rm = TRUE, avg = FALSE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                      {
                        if (self$length==0){
                          return(self)
                          }
                        onode.id = private$pgraph$nodes$dt$node.id

                        ## mark old nodes in curreng raph
                        old.gr = private$pgraph$nodes$gr

                        ## disjoin current graph
                        tmpg = private$pgraph$clone()

                        if (!is.null(graph))
                          tmpg = c(tmpg, graph)

                        if (!is.null(gr))
                        {
                          tmpg = c(tmpg, gG(breaks = gr))
                        }

                        tmpg$disjoin(by = by, na.rm = na.rm,
                                               avg = avg, FUN = FUN)

                        ## we will redo overlaps here
                        ## since we know tmpg$nodes$gr are disjoint
                        ## and we need coordinates to properly reorder
                        ## any disjoint walk nodes that came from a parent
                        ## node

                        old.gr$snode.id.old = old.gr$snode.id
                        map = gr2dt(sort(old.gr[,'snode.id.old'] %*% tmpg$nodes$gr[, 'snode.id']))
                        map = map[, .(query.id, subject.id, snode.id.old, snode.id.new = snode.id)]
                        map[, query.iid := 1:.N, by = query.id]
                        setkey(map, snode.id.old)

                        pn = copy(private$pnode)[, walk.iid := 1:.N, by = walk.id]
                        pn[, node.id := abs(snode.id)]

                        pn.new = merge(pn, map, by.x = 'node.id',
                                     by.y = 'snode.id.old', allow.cartesian = TRUE)
                        pn.new[, snode.id.new := sign(snode.id)*snode.id.new]
                        pn.new[, query.iid := sign(snode.id)*query.iid]

                        ## sorting on (now signed) query.iid
                        ## should cause the more rightwards nodes
                        ## in a query.id to be later in the path
                        setkeyv(pn.new, c("walk.id", "walk.iid", "query.iid"))
                        pn.new = pn.new[, .(walk.id, snode.id = snode.id.new)]

                        private$gWalkFromNodes(
                          snode.id = split(pn.new$snode.id, pn.new$walk.id),
                          graph = tmpg,
                          meta = copy(private$pmeta),
                          circular = private$pmeta$circular)

                        setkey(private$pnode, walk.id)
                        setkey(private$pedge, walk.id)
                        setkey(private$pmeta, walk.id)
                      },

                      #' @name simplify
                      #' @description
                      #'
                      #' gWalk simplify will merge reference adjacent nodes in walk +/- "by field"
                      #' (i.e. which requires the nodes to additional match on by field)
                      #' and simplify the underlying graph.
                      #'
                      #' Note that this only simplifies walks;
                      #' the associated graph may not be completely simplified 
                      #' since we retain all walk endpoints
                      #' 

                      simplify = function(by = NULL, na.rm = TRUE, avg = FALSE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                      {
                        tmpgr = private$pgraph$nodes$gr
                        byval = tmpgr %N% private$pgraph$nodes[private$pnode$snode.id]$gr
                        if (!is.null(by)) ## add user provided by to byval if exists
                          byval = paste(byval, values(tmpgr)[, by])

                        ## we will first simplify tmpg
                        ## then rematch the walk nodes to tmpg
                        tmpg = private$pgraph$copy
                        tmpg$nodes$mark(byval = byval)
                        tmpg$nodes$mark(node.id.og = 1:length(tmpg))
                        tmpg$simplify(by = "byval", FUN = NULL)

                        
                        ## we do simplify here just to keep track of nodes
                        ## will do a "real" simplify downstream
                        ## if the user wants to aggregate metadata from the walks

                        ## now we just need to match up walk nodes 
                        ## to simplified graph
                        ## some of our current walk nodes will become
                        ## "internal" so this just requires doing a simple
                        ## match up and then providing an updated set
                        ## snode.ids after removing NAs

                        newnode = copy(private$pnode)
                        newnode$snode.id = sign(newnode$snode.id)*
                          match(abs(newnode$snode.id), tmpg$nodes$dt$node.id.og)

                        newnode = newnode[!is.na(snode.id), ]

                        if (!is.null(FUN))
                        {
                          newgraph = private$pgraph$copy
                          newgraph$nodes$mark(byval = byval)
                          newgraph$simplify(by = "byval", FUN = FUN, na.rm = na.rm)
                        }
                        else
                        {
                          newgraph = tmpg
                          newgraph$mark(node.id.og = NULL)
                        }

                        newgraph$nodes$mark(byval = NULL)

                        private$gWalkFromNodes(
                                  snode.id = split(newnode[, snode.id], newnode[, walk.id]),
                                  graph = newgraph,
                                  meta = private$pmeta)

                        setkey(private$pnode, walk.id)
                        setkey(private$pedge, walk.id)
                        setkey(private$pmeta, walk.id)

                        return(self)
                      },

                      mark = function(...)
                      {
                        self$nodes$mark(...)
                        self$edges$mark(...)
                      },
                      
                      gtrack = function(name = NULL, stack.gap = 1e5, ...)
                      {
                        ## inherit gTrack arguments from gGraph
                        gt.args = private$pgraph$meta[intersect(c('colormap',
                                                                  'gr.colorfield',
                                                                  'col',
                                                                  'border',
                                                                  'name',
                                                                  'y1',
                                                                  'y0', 
                                                                  'lwd.border',
                                                                  'labels.suppress',
                                                                  'labels.suppress.gr',
                                                                  'label.suppress.grl',
                                                                  'yaxis'), names(private$pgraph$meta))]
                        
                        args = list(...)
                        for (arg in names(args))
                        {
                          gt.args[[arg]] = args[[arg]]
                        }

                        if (!is.null(name)){
                          gt.args[["name"]] = name
                          }
                        if (is.null(gt.args[['name']])){
                          gt.args[['name']] = 'gWalk'
                        }
                        gt.args[['draw.paths']] = TRUE


                        tmp.grl = self$grl
                        if (length(tmp.grl)>0){
                          values(tmp.grl)$is.cycle = self$dt$circular
                        }

                        gt.args[['data']] = tmp.grl

                        gt.args[['grl.labelfield']] = "name"

                        do.call(gTrack, gt.args)
                      },

                      ## TODO
                      ## experimental: this should be a special case of the gw2js function
                      bam2js = function(fn = "./"){
                          ## assume each walk in this object is a read pair
                          browser()
                          node.dt = self$nodes$dt
                          setkey(node.dt, "snode.id")
                          bam.dt = self$dt
                          edge.dt = self$edges$dt
                          ## hard flip the n2 side!!!
                          sides = setNames(c("left", "right"),
                                           c("right", "left"))
                          edge.dt[, n2.side := sides[n2.side]]
                          setkey(edge.dt, "sedge.id")
                          path.json =
                              lapply(bam.dt[, unique(walk.id)],
                                     function(ix){
                                         message(ix)
                                         dt = bam.dt[
                                             .(ix),
                                             .(pid = walk.id,
                                               cn = 1,
                                               type = ifelse(circular, "cycle", "path"),
                                               strand = "+")]
                                         snid = bam.dt[.(ix), unlist(snode.id)]
                                         gr = node.gr[as.character(snid)]
                                         ndt = node.dt[
                                             .(snid),
                                             .(iid = node.id,
                                               chromosome = as.character(seqnames),
                                               startPoint = start,
                                               endPoint = end,
                                               y = 0,
                                               title = paste0("read_", ix),
                                               type = "interval",
                                               strand)]
                                         seid = bam.dt[.(ix), unlist(sedge.id)]
                                         edt = edge.dt[
                                             .(seid),
                                             .(cid = edge.id,
                                               `source` = ifelse(n1.side=="left",
                                                                 -n1, n1),
                                               sink = ifelse(n2.side=="left",
                                                             -n2, n2),
                                               title = "",
                                               type = type,
                                               weight = 1)]
                                         out = as.list(dt)
                                         out$cids = edt
                                         out$iids = ndt
                                         return(out)
                                     })
                          jsonlite::write_json(list(walks = path.json),
                                               path = fn,
                                               pretty = TRUE,
                                               simplifyVector = TRUE,
                                               flatten = TRUE)
                          return(fn)
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
                      pgraph = NULL,
                     
                      gWalkFromNodes = function(snode.id,
                                              graph,
                                              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
                                              circular = NULL)  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first                                            
                      {

                        if (length(snode.id)==0)
                        {
                          private$pnode = data.table()
                          private$pedge = data.table()
                          private$pmeta = data.table()
                          private$pgraph = graph
                          return()
                        }

                          ## List of snode.id's
                          ## Get list of sedge.ids with only snode.ids from each element

                          if (!is.null(circular))
                          {
                            if (!is.logical(circular)){
                              stop('circular must be a logical vector')
                              }
                            }
                          else {
                            circular = rep(FALSE, length(snode.id))
                          }

                        if (!is.null(meta)) ## need this for some reason to get rid of pass by reference stickiness (doesn't work higher in call stack)
                        {
                          if (!is.data.table(meta))
                            meta = as.data.table(meta)
                          meta = copy(meta)
                        }
                        

                        if (is.null(names(snode.id)))
                          names(snode.id) = seq_along(snode.id)

                        private$pmeta = data.table(walk.id = seq_along(snode.id), name = names(snode.id), length = elementNROWS(snode.id), wid = NA, 
                                                   circular = circular)
                        
                        if (!is.null(meta))
                        {                       
                          if (nrow(meta) != length(snode.id))
                          {
                            stop('data.table of meta.data must be same length and order as node, edge, or GRangesList list input to gWalk constructor')
                          }
                          
                          GW.NONO = c('walk.id', 'circular', 'snode.id', 'sedge.id', 'wid', 'length', 'name', 'subject.id', 'query.id', 'i')
                          good.cols = setdiff(names(meta), GW.NONO)
                          
                          if (length(good.cols)>0)
                          {                            
                            private$pmeta = cbind(private$pmeta, meta[, good.cols, with = FALSE])
                          }
                        }

                        ## always mfaster to unlist than lapply
                        pnode = dunlist(unname(snode.id))[, .(walk.id = listid, snode.id = V1)]
                        
                        ## first need to check if the edges corresponding to the consecutive
                        ## node pairs in the input lists exist
                        pedge = pnode[, .(from = snode.id[-.N], to = snode.id[-1]), by = walk.id]
                        sedgesdt = copy(graph$sedgesdt)[, from := graph$gr$snode.id[from]][, to := graph$gr$snode.id[to]]
                        
                        if (any(private$pmeta$circular))
                        {
                          ix = private$pmeta$circular
                          wid = private$pmeta$walk.id[ix]
                            last = sapply(snode.id[ix], function(x) x[length(x)])
                          first = sapply(snode.id[ix], "[", 1)
                          pedge = rbind(pedge,
                                        data.table(
                                          walk.id = wid,
                                          from = last, to = first))
                        }
                                                
                        setkeyv(sedgesdt, c('from', 'to'))

                        ## just in case there is more than one edge
                        ## connecting a node pair (possible)
                        sedgesdt = unique(sedgesdt, by = c("from", "to"))

                        if (nrow(pedge)>0){
                          pedge$sedge.id = sedgesdt[.(pedge$from, pedge$to), sedge.id]
                          if (any(is.na(pedge$sedge.id))){
                            stop('One or more provided walks refers to a non-existent edge')
                          }
                        }                        
                        pnode[, walk.iid := 1:.N, by = walk.id]
                        pedge[, walk.iid := 1:.N, by = walk.id]

                        private$pnode = pnode
                        private$pedge = pedge
                        private$pgraph = graph
                      },


                      gWalkFromEdges = function(sedge.id,
                                              graph = NULL,
                                              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
                                              circular = NULL)  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first                                            
                      {
                        if (length(sedge.id)==0)
                        {
                          private$pnode = data.table()
                          private$pedge = data.table()
                          private$pmeta = data.table()
                          private$pgraph = graph
                          return()
                        }

                        if (!is.null(circular))
                        {
                          if (!is.logical(circular))
                            stop('circular must be a logical vector')
                        }
                        else {
                          circular = rep(FALSE, length(sedge.id))
                        }
                        
                        
                        private$pmeta = data.table(walk.id = seq_along(sedge.id), name = names(sedge.id), length = elementNROWS(sedge.id), wid = NA,
                                                   circular = circular)

                        if (!is.null(meta)) ## need this for some reason to get rid of pass by reference stickiness (doesn't work higher in call stack)
                        {
                          if (!is.data.table(meta))
                            meta = as.data.table(meta)
                          meta = copy(meta)
                          
                          if (nrow(meta) != length(sedge.id))
                          {
                            stop('data.table of meta.data must be same length and order as node, edge, or GRangesList list input')
                          }

                          GW.NONO = c('walk.id', 'circular', 'snode.id', 'sedge.id', 'wid', 'length', 'name', 'query.id', 'subject.id', 'i')
                          good.cols = setdiff(names(meta), GW.NONO)
                          
                            if (length(good.cols)>0)
                              {
                                private$pmeta = cbind(private$pmeta,
                                                      meta[, good.cols, with = FALSE])
                              }
                          }
                        
                        
                        pedge = dunlist(unname(sedge.id))[, .(walk.id = listid, sedge.id = V1)]
                        tmpe = graph$sedgesdt[.(pedge$sedge.id), .(from, to)]

                        if (any(is.na(tmpe$to))){
                          stop('one or more provided sedge.ids are out of bounds')
                        }
                        pedge = cbind(pedge, tmpe)
                        
                        ## first need to check if consecutive edges in the provided
                        ## walk share sink and source nodes
                        pedge[, matches.next := c(to[-.N] == from[-1], NA), by = walk.id]

                        if (any(!pedge$matches.next, na.rm = TRUE)){
                          stop('at least one pair of consecutive edges i and i+1 in the walk do not share a sink and source node, respectively')
                        }
                        pedge[, to := graph$dt$snode.id[to]]
                        pedge[, from := graph$dt$snode.id[from]]


                        pnode = pedge[, .(snode.id = c(from, to[.N])), by = walk.id]
                        private$pedge = pedge[, .(walk.id, sedge.id, to, from)]

                        pnode[, walk.iid := 1:.N, by = walk.id]
                        pedge[, walk.iid := 1:.N, by = walk.id]
                        private$pnode = pnode
                        private$pgraph = graph
                      }
                    ),
                    active = list(
                      ## Returns a GRangesList of walks in the graph
                      grl = function()
                      {
                        if (self$length == 0){
                            return(GRangesList(self$graph$gr[c()])[c()])
                        }
                        ## Does not get both strands only 1 strand
                        nix = private$pgraph$queryLookup(private$pnode$snode.id)$index
                        out = split(private$pgraph$gr[nix], private$pnode$walk.id)[as.character(1:self$length)]
                        values(out) = private$pmeta
                        return(out)
                      },

                      meta = function()
                      {
                        return(private$pmeta)
                      },

                      nodes = function()
                      {
                        return(private$pgraph$nodes[private$pnode$snode.id])
                      },

                      copy = function() self$clone(),

                      ## returns a length(self) logical vector specifying whether
                      ## each walk is circular or not
                      circular = function()
                      {
                        return(private$pmeta$circular)
                      },

                      length = function() {
                        return(nrow(private$pmeta))
                      },

                      edges = function()
                      {
                        if (self$length == 0){
                          return(gEdge$new(graph = private$pgraph))
                        }
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
                      
                      graph = function()
                      {
                        return(private$pgraph)
                      },
                      
                      dt = function() {
                        if (self$length==0){
                            return(data.table())
                            }
                        ix = 1:self$length
                        return(self$dts(ix))              
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
  sedge.id = do.call('c', lapply(gWalk.list, function(x) x$sdt$sedge.id))

  return (gWalk$new(sedge.id = sedge.id, graph = gWalk.list[[1]]$graph))
}


#' @name gW
#' @title create gWalk
#' @description
#'
#' Wrapper that instantiates a gWalk object from a variety of different inputs.  If a gGraph is provided as input (in conjunction
#' with sedge.id, snode.id, or grl  input) then will check for existence of the provided walks in the provided graph, and
#' will error out if those walks do not exist).
#'
#' 
#' @examples
#'
#' ## read in GRangesList of walks
#' walks.grl = readRDS(system.file('extdata', 'walks.grl.rds', package = "gGnome"))
#'
#' ## read in gGraph
#' gg = readRDS(system.file('extdata', 'walks.gg.rds', package="gGnome"))
#'
#' ## create gWalks from GRangesList - creates graph in which walks are disconnected, i.e. phased,
#' ## so there may be several intervals per reference region
#' gW(grl = walks.grl)
#'
#' ## disjoin graph from GRangesList, these refer to a single disjoint graph
#' ## where there is one interval per reference region
#' gW(grl = walks.grl, disjoin = TRUE)
#' 
#' ## thread GRangesList onto existing graph
#' gW(grl = walks.grl, graph = gg)
#'
#' ## from lists of signed node ids
#' nid.list = list(c(-871, 89, 90), c(1105, 911, 912, 912))
#' gW(snode.id = nid.list, graph = gg)
#'
#' ## from lists of signed edge ids
#' eid.list = list(c(-5266, -9765, -9764, -5263), c(-5272, -9317, -5270, -6638, -9370, -5267))
#' gW(sedge.id = eid.list, graph = gg)
#' 
#' @export
#'
#' @param snode.id list of (possibly negative) integer signed node ids in a gGraph
#' @param sedge.id list of (possibly negative) integer signed edge ids in a gGraph, 
#' @param grl  GRangesList of walks in the genome, each walk is a (stranded) GRangesList, if graph is NULL will build a gGraph from the provided GRangesLists
#' @param graph an existing graph from which to build a graph (used in conjunction with snode.id, sedge.id, or grl), if used with grl will attempt to match up the provided grl's with edges in the graph, if the GRL walks non-existent edges, then will error out
#' @param disjoin relevant for grl input when graph is NULL. If TRUE (default FALSE) will create a graph taking the GRanges disjoin of the GRanges comprising the provided GRL
#' @return gWalk object of provided walks, with pointers back to the graph to which they refer.
gW = function(snode.id = NULL,
              sedge.id = NULL,
              grl = NULL,
              graph = NULL,
              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
              circular = NULL,  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first
              disjoin = FALSE ## flag whether to collapse ie disjoin the intervals
              )
{
  return(gWalk$new(snode.id = snode.id,
                    sedge.id = sedge.id,
                    grl = grl,
                    graph = graph,
                    meta = meta,
                    circular = circular,
                    disjoin = disjoin))
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
'[.gWalk' = function(obj, i = NULL, ...){
  walks = obj$clone()
  inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()),obj$dts(makelists = FALSE), parent.frame(2)), error = function(e) NULL)
  if (is.null(inew))
    inew = i ## just give up      
  walks$subset(inew)
  return(walks)
}



#' @name %&%'
#' @title subset x on y edges
#' @description
#'
#'
#' gEdge1 %&% gEdge2/Junction returns the subsets of gEdge1 that overlaps gEdge2/Junction
#'
#' @return subset of gEdge that overlaps 
#' @exportMethod %&%
#' @aliases %&%, gEdge Method
#' @author Rick Mortenson
#setGeneric('%&%', function(x, ...) standardGeneric('%&%'))
setMethod("%&%", signature(x = 'gNode'), function(x, y) {
    if (is.character(y)){
        y = parse.gr(y)
    }
    return(x$gr[gr.in(x$gr, y)])
})

edge.queries = function(x, y) {   
    if (is.character(y)){        
        y = parse.gr(y)
    }
    gr = unlist(x$gr, use.names = FALSE)
    if (inherits(y, "gEdge")) {       
        overlaps=ra.overlaps(x$junctions$grl, y$junctions$grl)       
        good.ix=overlaps[, "ra1.ix"]                      
        good.ix=unique(good.ix)
        if (is.na(good.ix)){
            return("No overlaps")
        }
        else{
            return(x[good.ix])
        }
    }
    
 else if (inherits(y, "Junction")) {             
    juncs=x$junctions$grl
    overlaps=ra.overlaps(juncs, y$grl)
    good.ix=overlaps[, "ra1.ix"]                      
    good.ix=unique(good.ix)
    if (is.na(good.ix)){
        return("No overlaps")
    }
    else{
            return(x[good.ix])
    }
 }
}


#' @name merge.Junction
#' @title merge junctions by overlaps with padding
#'
#' Merges a set of junctions and keeps "seen.by" metadata of junction origin
#' using the argument names to this function
#'
#' @examples
#'
#' ## wil output a Junction object with metadata seen.by.svaba etc.
#' ## will pad with 500 bases prior to merging
#'
#' svaba = jJ(system.file('extdata', "HCC1143.svaba.somatic.sv.vcf", package = "gGnome"))
#' delly = jJ(system.file('extdata', "delly.final.vcf.gz", package = "gGnome"))
#' novobreak = jJ(system.file('extdata', "novoBreak.pass.flt.vcf", package = "gGnome"))
#' 
#' ## merge(svaba = svaba, delly = delly, caller3 = novobreak, pad = 500)
#' 
#' @param ... GRangesList representing rearrangements to be merged
#' @param pad non-negative integer specifying padding
#' @param ind  logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing
"merge.Junction" = function(..., pad = 0, ind = FALSE)
{
  list.args = list(...)
  Junction$new(do.call(ra.merge, c(lapply(list.args, function(x) x$grl), list(pad = pad, ind = ind))))
}

setMethod("%&%", signature(x = 'gEdge'), edge.queries)
    
#' @name default.agg.fun.generator
#' @title default.agg.fun.generator)
#' @description
#'
#' Generates aggregation function for gGraph reduce and disjoin.
#' 
#' Returns a function that takes x and sums (+/- na.rm) if non-character, pastes with given sep if numeric, and concatenates if list
#' (note this returns a function, it is not the actual functioon)
#' 
#' @author Marcin Imielinski
#' @return a function with the given characteristics
#' @keywords internal
#' @noRd
default.agg.fun.generator = function(na.rm = TRUE, avg = FALSE, sep = ',')
{
  function(x)
  {
    if (length(x) == 1){
        out = x
        }
    else if (all(is.na(x))){
        out = x[1]
        }
    else if (is.numeric(x))
    {
      if (avg){
        out = mean(x, na.rm = na.rm)      
        }
      else{
          out = sum(x, na.rm = na.rm)
          }
      }
    else if (is.integer(x))
    {
      if (avg){
        out = as.integer(median(x, na.rm = na.rm))
        }
      else{
        out = as.integer(sum(x, na.rm = na.rm))
        }
      }
    else if (is.logical(x)){
        out = as.logical(sum(x, na.rm = na.rm))
        }
    else if (is.character(x) | is.factor(x)){
          out = paste(unique(x[!is.na(x)]), collapse = sep)
        }
    else if (is.list(x))
      {
        out = do.call(c, x)
      }
    else
      {
        stop('gGraph default aggregation failed for unknown meta data type (numeric, integer, logical, character, list)')
      }


    return(out)

  }
}

#' @name lengths
#' @title lengths
#' @description
#'
#' @param x a gWalk object
#'
#' @return the number of nodes in each gWalk
setMethod("lengths", c("gWalk"),
          function(x) {
            return(lengths(x$grl))
          })

#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gGraph object
#'
#' @return the seqinfo of this graph
setMethod("seqinfo", c("gGraph"),
          function(x) {
            return(seqinfo(x$gr))
          })


#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gNode object
#'
#' @return the seqinfo of this gNode
setMethod("seqinfo", c("gNode"),
          function(x) {
            return(seqinfo(x$gr))
          })


#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gEdge object
#'
#' @return the seqinfo of this gEdge
setMethod("seqinfo", c("gEdge"),
          function(x) {
            return(seqinfo(x$grl))
          })

#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a Junction
#'
#' @return the seqinfo of this Junction
setMethod("seqinfo", c("Junction"),
          function(x) {
            return(seqinfo(x$grl))
          })

#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gWalk
#'
#' @return the seqinfo of this gWalk
setMethod("seqinfo", c("gWalk"),
          function(x) {
            return(seqinfo(x$grl))
          })


#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gGraph object
#'
#' @return the seqlengths of this graph
setMethod("seqlengths", c("gGraph"),
          function(x) {
            return(seqlengths(x$gr))
          })



#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gNode object
#'
#' @return the seqlengths of this gNode
setMethod("seqlengths", c("gNode"),
          function(x) {
            return(seqlengths(x$gr))
          })

#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gEdge object
#'
#' @return the seqlengths of this gEdge
setMethod("seqlengths", c("gEdge"),
          function(x) {
            return(seqlengths(x$grl))
          })


#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a Junction
#'
#' @return the seqlengths of this Junction
setMethod("seqlengths", c("Junction"),
          function(x) {
            return(seqlengths(x$grl))
          })

#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gWalk
#'
#' @return the seqlengths of this gWalk
setMethod("seqlengths", c("gWalk"),
          function(x) {
            return(seqlengths(x$grl))
          })


#' @name width
#' @title width
#' @description
#'
#' @param x a gWalk
#'
#' @return the width of each walk in this gWalk
setMethod("width", c("gWalk"),
          function(x) {
            return(x$dt$wid)
          })


#' @name width
#' @title width
#' @description
#'
#' @param x a gNode
#'
#' @return the width of each gNode
setMethod("width", c("gNode"),
          function(x) {
            return(width(x$gr))
          })


#' @name refresh
#' @title refresh
#' @description
#' Updates gWalk object 
#' 
#' @param gWalk object
#'
#' @return gWalk object
#' @export
setMethod("refresh", "gWalk",
          function(x) {
              return(gWalk$new(snode.id = x$dts()$snode.id,
                               meta = x$meta,
                               graph = x$graph))
          })


#' @name gNode subset
#' @title subset gNode on overlaps
#' @description
#'
#' @param x a gNode
#'
#' @return the subset of gNodes overlapping with the query
setMethod("%&%", signature(x = 'gNode'), function(x, y) {

  if (inherits(y, 'gNode'))
    y = y$gr

  if (inherits(y, 'gEdge'))
    y = unlist(y$grl)

  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGrangesList'))
    y = unlist(y$grl)

  if (is.character(y)){
    y = parse.gr(y)
  }
  return(x[gr.in(x$gr, y)])
})

#' @name gNode over
#' @title subset gNode on overlaps
#' @description
#'
#' @param x a gNode
#'
#' @return logical vector of length(x) specifying overlap 
setMethod("%^%", signature(x = 'gNode'), function(x, y) {

  if (inherits(y, 'gNode'))
    y = y$gr

  if (inherits(y, 'gEdge'))
    y = unlist(y$grl)
    
  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGrangesList'))
    y = unlist(y)

  if (is.character(y)){
    y = parse.gr(y)
  }
  return(gr.in(x$gr, y))
})



#' @name gEdge.in
#' @title subset gEdge on overlaps
#' @description
#'
#' @param x a gEdge
#'
#' @return the subset of gEdges overlapping with the query
setMethod("%&%", signature(x = 'gEdge'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = y$grl
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }


  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGRangesList'))
  {
    ix = setdiff(ra.overlaps(x$grl, y)[,1], NA)
    return(x[ix])
  }

   if (is.character(y)){
     y = parse.gr(y)
   }

  return(x[grl.in(x$grl, y)])
})

#' @name gEdge.over
#' @title gEdge over
#' @description
#'
#' @param x a gEdge
#'
#' @return logical vector of length(x) which is TRUE if an overlap exists
setMethod("%^%", signature(x = 'gEdge'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = y$grl
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }

  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGRangesList'))
  {
    ix = setdiff(ra.overlaps(x$grl, y)[,1], NA)
    return(1:length(x) %in% ix)
  }

   if (is.character(y)){
     y = parse.gr(y)
   }

  return(grl.in(x$grl, y))
})



#' @name gWalk.subset
#' @title subset gWalk on overlaps
#' @description
#'
#' @param x a gWalk
#'
#' @return the subset of gWalks overlapping with the query
setMethod("%&%", signature(x = 'gWalk'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = unlist(y$grl)
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }

  if (is.character(y)){
    y = parse.gr(y)
  }
  
  return(x[grl.in(x$grl, y)])
})



#' @name gWalk.over
#' @title subset gWalk on overlaps
#' @description
#'
#' @param x a gWalk
#'
#' @return logical vector of length(x) which is TRUE if an overlap exists
setMethod("%^%", signature(x = 'gWalk'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = unlist(y$grl)
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }

  if (is.character(y)){
    y = parse.gr(y)
  }
  
  return(grl.in(x$grl, y))
})




#' @name Junction.subset
#' @title subset gEdge on overlaps
#' @description
#'
#' @param x a Junction
#' @param y a Junction, gEdge, GenomicRanges, string
#' @return the subset of gEdges overlapping with the query
setMethod("%&%", signature(x = 'Junction'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = y$grl
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }

  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGRangesList'))
  {
    ix = setdiff(ra.overlaps(x$grl, y)[,1], NA)
    return(x[ix])
  }

   if (is.character(y)){
     y = parse.gr(y)
   }

  return(x[grl.in(x$grl, y)])
})

#' @name Junction.over
#' @title Junction over
#' @description
#'
#' @param x a Junction
#' @param y a Junction, gEdge, GenomicRanges, string
#' @return logical vector of length(x) which is TRUE if an overlap exists
setMethod("%^%", signature(x = 'Junction'), function(x, y) {
  if (inherits(y, 'Junction') | inherits(y, 'gEdge'))
  {
    y = y$grl
  }

  if (inherits(y, 'gNode'))
  {
    y = y$gr
  }

  if (inherits(y, 'GRangesList') | inherits(y, 'CompressedGRangesList'))
  {
    ix = setdiff(ra.overlaps(x$grl, y)[,1], NA)
    return(1:length(x) %in% ix)
  }

   if (is.character(y)){
     y = parse.gr(y)
   }

  return(grl.in(x$grl, y))
})



#' @name jJ
#' @title parse junctions from a variety of formats
#'
#' @description Parsing various formats of structural variation data into junctions.
#'
#' @usage read.junctions(rafile,
#' keep.features = T,
#' seqlengths = NULL,
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
#' @return a \code{Junction} object
#'
#' @export
#' @author Xiaotong Yao
#' @importFrom VariantAnnotation readVcf
#' @import data.table
jJ = function(rafile,
              keep.features = T,
              seqlengths = NULL,
              chr.convert = T,
              geno=NULL,
              flipstrand = FALSE,
              swap.header = NULL,
              breakpointer = FALSE,
              seqlevels = NULL,
              force.bnd = FALSE,
              skip = NA,
              get.loose = FALSE){
  return(Junction$new(
                    rafile,
                    keep.features = keep.features,
                    seqlengths = seqlengths,
                    chr.convert = chr.convert,
                    geno=geno,
                    flipstrand = flipstrand,
                    swap.header = swap.header,
                    breakpointer = breakpointer,
                    seqlevels = seqlevels,
                    force.bnd = force.bnd,
                    skip = skip,
                    get.loose = get.loose)
         )
}



