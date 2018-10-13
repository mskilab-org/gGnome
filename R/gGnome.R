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
                        private$ptimestamp = graph$timestamp ## inherit timestamp from graph object
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
                        self$check
                        if (length(self)>0)
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
                          }
                      },

                      #' @name rep
                      #' @description
                      #'
                      #' Creates "bubbles" in the graph by replicating the nodes in the argument.  Node replication
                      #' replicates edges going in and out of all replicated nodes.  If an edge connects a pair of replicated nodes
                      #' that edge will be replicated across all pairs of those replciated nodes.   Walk replication will create "longer bubbles"
                      #' with fewer edges getting replicated i.e. it will only replicate intra-walk edges within each walk replicate (but not between
                      #' separate walk replicates).
                      #'
                      #' New graph keeps track of the parent node and edge ids in the original graph using node metadata parent.node.id
                      #' and edge metadata parent.edge.id i.e. the replicated nodes will be connected to the sources of the original nodes
                      #' and if replicated nodes connect to each other, then there will exist an edge connecting
                      #' all of their instances to each other.
                      #'
                      #' @param nodes = gNode object must point to a node in the graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                      #' @param times  scalar or vector of length self$length specifying how many times to replicate each of the nodes.
                      #' @author Marcin Imielinski
                      #' @return Returns a pointer to the new nodes
                      rep = function(times)
                      {
                        if (length(self)>0)
                        {
                          private$pgraph$rep(self, times)
                        }
                        return(invisible(self))
                      },

                      #' @name swap
                      #' @description
                      #'
                      #' Swap nodes with granges, grl, or Gwalks.
                      #' Provided replacement vector must be the same length as the inputted nodes, resulting in each node being "swapped" by the provided 
                      #' interval, node, grl (representing a walk), or gWalk.  The replacement will inherit left and right edges for the removed node.
                      #' If the replacement is a walk, then the left side of the first node in the walk will inherit the edges that were
                      #' previously to the left of the node being replaced, and right side of the last node of the walk will inherit the edges
                      #' that were previously to the right of the node being replaced.  
                      #'
                      #' Note: these replacement obey the orientation of the arguments.  So if the node to be replaced is flipped (- orientation with
                      #' respect to the reference, then it's "left" is to the right on the reference.  Similarly for walks whose first interval
                      #' is flipped with respect to the reference, the left edges will be attached to the right of the node on the reference.
                      #' 
                      #' @param replacement  GRanges, GRangesList, or gWalk object whose length is the length(nodes)
                      #' @author Marcin Imielinski
                      #' @return gGraph (also modified in place) with nodes annotated with parent.node.id, parent.rep
                      swap = function(replacement)
                      {
                        if (length(self)>0)
                          {
                            private$pgraph$swap(self, replacement)
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
                        self$check
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

                      #' @name clusters
                      #' @description
                      #'
                      #' Marks nodes in graph with metadata field $cluster
                      #' based on one of several algorithms, selected by mode.
                      #'
                      #' Unlike the gGraph version of this function, this usage enables computation of clusters
                      #' based on a subset of nodes in the graph (i.e. ignoring certain nodes)
                      #' while still marking the original graph (and leaving the excluded graph with an NA
                      #' value for "cluster"
                      #' 
                      #'
                      #' @param weak character scalar that can take one of the following possible values - "weak" or "strong" specifying weakly or strongly connected components, walktrap specifying cluster_walktrap community detection
                      #' @author Marcin Imielinski
                      clusters = function(mode = 'weak')
                      {
                        graph = self$graph$copy
                        graph$nodes$mark(parent.node.id = graph$nodes$dt$node.id)

                        sgraph = graph[self$dt$node.id, ]

                        sgraph$clusters(mode)

                        self$graph$nodes$mark(cluster = as.numeric(NA))
                        self$graph$nodes[sgraph$nodes$dt$parent.node.id]$mark(cluster = sgraph$nodes$dt$cluster)
                        return(invisible(self$graph))
                      },
                      

                      #' @name ego
                      #' @description
                      #'
                      #' returns all nodes with k degrees of separation (i.e. "order")
                      #' from these
                      #'
                      #' @param order integer edge distance from these (seed) nodes to return
                      #' @param mindist minimum distance from these (seed) nodes
                      #' @return gNode object comprising nodes within the ego of this node
                      #' @author Marcin Imielinski
                      ego = function(order = 0, mindist = 0)
                      {
                        self$check
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
                        self$check
                        message('gNode object of length ', self$length)

                        if (self$length>0)
                        {
                          HEAD = 1:pmin(4, self$length)                         
                          print(self$gr[HEAD, ])
                          more = self$length-HEAD[length(HEAD)]
                          if (more>0)
                          {
                            message('... (', more,' additional nodes)')
                          }
                        }
                      }
                    ),
                    

                    private = list(
                      ## stores the unsigned node.id's of the nodes in private$pgraph in this gNode
                      pnode.id = NULL,

                      ptimestamp = NULL,
                      
                      ## Stores the index and reverse index of the pnode.id's in this gNode
                      pindex = NULL,
                      prindex = NULL,
                      
                      ## Stores the signs of the nodes associated with indices in private$pindex
                      porientation = NULL,
                      
                      ## Pointer to gGraph that this gNode belongs to
                      pgraph = NULL
                    ),


                    active = list(

                      ## object is stale if the recorded timestamp of the gGraph
                      ## != timestamp the actual gGraph pointed to by pgraph 
                      ## suggesting that the indices are no longer valid
                      stale = function() private$ptimestamp != private$pgraph$timestamp,


                      ## checks if object is stale i.e.
                      check = function() if (self$stale) stop('object is stale, underlying gGraph has changed. You will need to re-instantiate.'),

                      #' @name length
                      #' @description
                      #' 
                      #' Returns the number of nodes in this gNode Object.
                      #' 
                      #' @return Number of nodes in this gNode Object
                      length = function()
                      {
                        self$check
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
                        self$check
                        return(private$pgraph)
                      },


                      #' @name dt
                      #' @description
                      #'
                      #' Returns the GRanges of the nodes in the gNode as a data.table.
                      #'
                      #' @return data.table GRanges of the nodes coverted to a data.table
                      dt = function() {
                        self$check
                        return(as.data.table(private$pgraph$gr[private$pindex]))
                      },
                      
                      
                      #' @name gr
                      #' @description
                      #'
                      #' Returns a GRanges of the nodes in the gNode.
                      #'
                      #' @return GRanges of the nodes
                      gr = function() {
                        self$check
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
                        self$check
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
                        self$check
                        return(ifelse(private$porientation == 1, private$pnode.id, -private$pnode.id))
                      },                       

                      ## returns flipped version of this node
                      flip = function()
                      {
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
                        ed = c(self$eleft, self$eright)
                        if (length(ed)>0)
                          {
                            ed =ed[!duplicated(edge.id)]
                          }

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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
                        deg = rowSums(cbind(
                          private$pgraph$edgesdt[n1.side == 'left', .N, keyby = n1][.(private$pnode.id), N],
                          private$pgraph$edgesdt[n2.side == 'left', .N, keyby = n2][.(private$pnode.id), N],
                          private$pgraph$edgesdt[n1.side == 'right', .N, keyby = n1][.(private$pnode.id), N],
                          private$pgraph$edgesdt[n2.side == 'right', .N, keyby = n2][.(private$pnode.id), N]),
                          na.rm = TRUE)

                        return(deg)
                      },

                      terminal = function()
                      {
                        self$check
                        if (self$length==0)
                          {
                            return(GRanges(seqlengths = seqlengths(self)))
                          }
                        ix = which(self$ldegree==0 | self$rdegree==0)
                        return(self[ix]$loose)
                      },

                      ldegree = function()
                      {
                        self$check
                        ldeg = rowSums(cbind(private$pgraph$edgesdt[, sum(n1.side == 'left'), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == 'left'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)

                        rdeg = ldeg
                        if (any(private$porientation<0)){
                          rdeg = rowSums(cbind(private$pgraph$edgesdt[ , sum(n1.side == 'right'), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == 'right'), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)}

                        return(pmax(0, ifelse(private$porientation>0, ldeg, rdeg), na.rm = TRUE))
                      },

                      rdegree = function()
                      {
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
'[.gNode' = function(obj, i = NULL, with = TRUE, ...)
{
  nodes = obj$clone()
  if (with)
    {
      ## yes I finally figured it out!!!! grrrrrr
      inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()), nodes$dt, parent.frame(2)), error = function(e) NULL)
      if (is.null(inew))
      {
        inew = i ## just give up
      }
    }
  else
  {
    inew = i
  }
  nodes$subset(inew)
  return(nodes)
}


## ================= gEdge class definition ================== ##
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
                        private$ptimestamp = private$pgraph$timestamp
                        
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
                        self$check
                        if (length(self)>0)
                          {
                            NONO.FIELDS = c('from', 'to', 'sedge.id', 'edge.id')
                            args = list(...)
                            
                            if (any(names(args) %in% NONO.FIELDS)){
                              stop(paste('Cannot modify the following reserved gEdge metadata fields:', paste(NONO.FIELDS, collapse = ', ')))
                            }
                            
                            for (nm in names(args)){
                              private$pgraph$annotate(nm, args[[nm]],
                                                      c(private$psedge.id, -private$psedge.id), "edge")
                            }
                          }
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
                        self$check
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
                        self$check
                        message("gEdge object with ", self$length, " edges")
                        if (self$length>0)
                        {
                          HEAD = 1:pmin(4, self$length)
                          tmpdt = self$dt[HEAD, ]
                          more = self$length-HEAD[length(HEAD)]
                          print(tmpdt)
                          if (more>0)
                            message('... (', more,' additional edges)')
                        }


                      }
                    ),

                    private = list(
                      ## data.table of the edges in this gEdge
                      pedges = NULL,

                      ## edge.id's and sedge.id's of the edges in this gEdge
                      pedge.id = NULL,
                      psedge.id = NULL,

                      ptimestamp = NULL,
                      
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
                        self$check
                        return(private$pgraph)
                      },


                      copy = function() self$clone(),


                      ## object is stale if the recorded timestamp of the gGraph
                      ## != timestamp the actual gGraph pointed to by pgraph 
                      ## suggesting that the indices are no longer valid
                      stale = function() private$ptimestamp != private$pgraph$timestamp,


                      ## checks if object is stale i.e.
                      check = function() if (self$stale) stop('object is stale, underlying gGraph has changed. You will need to re-instantiate.'),

                      #' @name length
                      #' @description
                      #' 
                      #' Returns the number of edge pairs in this gEdge
                      #' 
                      #' @return Number of edge pairs in this gEdge
                      length = function()
                      {
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
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
                        self$check
                        sides = c('left', 'right')
                        return(copy(convertEdges(private$pgraph$gr, private$pedges, metacols = TRUE, cleanup = FALSE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
                      },

                      #' @name class
                      #' @description
                      #' 
                      #' recomputes the "class" of each gEdge (note: automatically done during
                      #' gGRaph instantiation and added to edge metadata table as $class column
                      #' using the standard taxonomy of "DUP", "DEL", "TRA", and "INV" but
                      #' based on the orientation and reference location of junctions
                      #' but appending suffix "-like" to enunciate the fact that these are only
                      #' preliminary classifications of these junctions
                      class = function()
                      {
                        self$check
                        nodes = self$graph$nodes$gr
                        ed = self$dt[, .(n1, n1.side, n2, n2.side)]
                        ed[, intrachr :=  as.logical(seqnames(nodes)[n1] ==  seqnames(nodes)[n2])]
                        ed[intrachr == TRUE,
                           ":="(
                             n1.coord = ifelse(n1.side == 'right', end(nodes)[n1], start(nodes)[n1]),
                             n2.coord = ifelse(n2.side == 'right', end(nodes)[n2], start(nodes)[n2])
                           )]
                        ed[intrachr == TRUE,  span12 :=  n1.coord-n2.coord]
                        ed[, class := ifelse(intrachr,
                                      ifelse(n1.side == n2.side, 'INV-like',
                                      ifelse(span12<0 &  n1.side == 'right' | span12>0 & n2.side == 'right',
                                      ifelse(abs(span12)==1, 'REF', 'DEL-like'),
                                      'DUP-like')
                                      ), 'TRA-like')]
                        return(ed$class)
                      },

                      #' @name junctions
                      #' @description
                      #'
                      #' Coverts the edges in this gEdge to junctions and returns these as a Junctions Object.
                      #' See Junctions Object documentation to see how junctions are defined in this package.
                      junctions = function()
                      {
                        self$check
                        grl = self$grl
                        return(Junction$new(grl))
                      },
                      

                      grl = function()
                      {
                        self$check
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
                        self$check
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
                        self$check
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
'[.gEdge' = function(obj, i = NULL, from = NULL, to = NULL, with = TRUE, ...)
{
  edges = obj$clone()

  if (!is.null(from) | !is.null(to))
  {

    tmpdt = edges$graph$sedgesdt
    sid = tmpdt$sedge.id

    if (!is.null(from))
    {
      n1 = edges$graph$nodes[from]
      setkeyv(tmpdt, c("from"))
      sid = intersect(sid, tmpdt[.(n1$dt$index), sedge.id])
    }
    
    if (!is.null(to))
    {
      n2 = edges$graph$nodes[to]
      setkeyv(tmpdt, c("to"))
      sid = intersect(sid, tmpdt[.(n2$dt$index), sedge.id])      
    }

    sid = setdiff(sid, NA)

    edges$subset(sid)
  }

  if (length(edges)==0)
    return(edges)


  if (deparse(substitute(i)) == "NULL")
    return(edges)

  if (with)
  {
    ## yes I finally figured it out grrrrrr
    inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()),edges$dt, parent.frame(2)), error = function(e) NULL)
    if (is.null(inew))
      inew = i ## just give up          
  } else
  {
    inew = i
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



## ================== Junction class definition ================== ##
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
                                   stop('cannot mix positive with negative subscripts for Junction object')
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
                         #' Prints out the Junction Object. Prints the length and the GRangesList of the junctions.
                         print = function()
                         {
                           message("Junction Object with ", self$length, " junctions\n")
                           if (self$length>0)
                             {
                               HEAD = 1:pmin(4, self$length)
                               bp = grl.unlist(self$grl[HEAD])
                               bpstr = gr.string(bp)
                               bpdt = data.table(str = bpstr, ix = bp$grl.ix)[, .(junc = paste(str, collapse = ' <-> ')), keyby = ix][.(HEAD), ]

                               if (ncol(self$dt)>0)
                               {
                                 bpdt = cbind(bpdt[, "junc", with = FALSE], self$dt[HEAD, ])
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
                         #' Returns a GRanges representing the spots a genome needs to be broken for these
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


                         #' @name left
                         #' @description
                         #'
                         #' Returns a GRanges of "left" ends of this junction object.  Left / right is just an arbitrary
                         #' designation (left refers to the "first" endpoint in each junction
                         #' 
                         #' @return GRanges of left end points of this junction object
                         left = function()
                         {
                           ## Get the breakpoints, left most end of - strand, right most end of + strand
                           return(grl.pivot(private$pjuncs)[[1]])
                         },

                         #' @name right
                         #' @description
                         #'
                         #' Returns a GRanges of "right" ends of this junction object.  Left / right is just an arbitrary
                         #' designation (right refers to the "second" endpoint in each junction)
                         #' 
                         #' @return GRanges of right end points of this junction object
                         right = function()
                         {
                           ## Get the breakpoints, left most end of - strand, right most end of + strand
                           return(grl.pivot(private$pjuncs)[[2]])
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
                         #' @return flipped Junction object
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
#o' 
#' Returns a new Junction Object which is the union of x and y.
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortensen
#' @return new Junction Object containing the union of x and y

#' @export
setMethod("union", c('Junction', "Junction"), function(x, y, pad = 0, ignore.strand = FALSE, ...)
{
  newJunc=c(x, y)
  return(unique(newJunc, pad, ignore.strand))
})

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
'[.Junction' = function(obj, i, j, with = TRUE){
  juncs = obj$clone()
  if (!missing(i))
  {
    if (with)
      {
        inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()),juncs$dt, parent.frame(2)), error = function(e) NULL)
        if (is.null(inew))
          inew = i ## just give up
      }
    else
    {
      inew = i
    }
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



#' @name union
#' Returns a new Junction object which is the union of x and y.
#' 
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortenson
#' @return new Junction Object containing the union of x and y
setMethod("union", c("Junction", "Junction"),
          function(x, y) {
              newJunc=c(x, y)
              newJunc$removeDups()
              return(newJunc)
          })

#' @name setdiff
#' Returns a new Junction object which is the difference between x and y (id's).
#'
#' @param x a Junction Object
#' @param y a Junction Object
#' @author Rick Mortenson
#' @return new Junction containing the difference between x and y
setMethod("setdiff", c("Junction", "Junction"),
          function(x, y) {
              ## Make sure that both come from the same graph                                         
              overlaps=ra.overlaps(x$grl, y$grl)
              overlaps=overlaps[, "ra1.ix"]
              all.ix=c(1:length(x$grl))
              dif.ix=setdiff(all.ix, overlaps)             
              return(Junction$new(x$grl[dif.ix]))
              
          })

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
                             private$stamp()
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
                       disjoin = function(gr = NULL, by = NULL, collapse = TRUE, na.rm = TRUE, avg = FALSE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                       {
                         this = self

                         if (length(self)==0) 
                           return(invisible(self))

#                         if (!is.null(gr))
 #                          this = c(self, gG(breaks = gr))

                         dgr = this$nodes$gr

                         if (!is.null(gr))
                         {
                           dgr = gr.stripstrand(grbind(dgr, gr %&% reduce(gr.stripstrand(dgr))))
                         }

                         ## now we disjoin all the nodes and merge them back with their parents
                         dnodes = disjoin(dgr) %*% this$nodes$gr

                         ## nmap maps our original node ids to the new node ids
                         ## (but we throw out all "internal nodes" created by the disjoin
                         ## here since we will use nmap to map edges onto the new disjoined nodes
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
                         if (nrow(final.edges)>0 & nrow(edges)>0)
                           {
                             if (length(metacolse)>0)
                             {
                               final.edges = cbind(final.edges, edges[.(newedges$sedge.id), metacolse, with = FALSE])
                             }
                             
                             final.edges[, type := ifelse(is.na(edges[.(newedges$sedge.id), type]), type, edges[.(newedges$sedge.id), type])]
                           }
                             
                         if (!collapse)
                         {
                           private$gGraphFromNodes(dt2gr(final.nodes, seqlengths = seqlengths(dnodes)), final.edges)
                           return(invisible(self))
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
                         ## but have n1 and n2 flipped (ie they are the identical edge
                         ## but our aggregation will not recognize it)
                         ## to do this we standardized the edge notation
                         ## so that n1 < n2, where the "order" is arbitrary but standard / unambiguous
                         ## (ordering also includes the "side")

                         ## tmp.n1 = final.edges$n1
                         ## tmp.n2 = final.edges$n2
                         ## tmp.n1.side = final.edges$n1.side
                         ## tmp.n2.side = final.edges$n2.side

                         old.final.edges = copy(final.edges)
                         ## factor then convert to integer so we can apply an (arbitrary) order
                         n1.tag = final.edges[, paste(n1, n1.side)]
                         n2.tag = final.edges[, paste(n2, n2.side)]

                         ulev = union(n1.tag, n2.tag)
                         n1.tag = as.integer(factor(n1.tag, ulev))
                         n2.tag = as.integer(factor(n2.tag, ulev))

                         ## final.edges[, n1 := pmin(tmp.n1, tmp.n2)]
                         ## final.edges[, n2 := pmax(tmp.n1, tmp.n2)]

                         ## final.edges[, n1.side := ifelse(tmp.n1<tmp.n2, tmp.n1.side, tmp.n2.side)]
                         ## final.edges[, n2.side := ifelse(tmp.n1>=tmp.n2, tmp.n1.side, tmp.n2.side)]

                         n1 = ifelse(n1.tag<n2.tag, final.edges$n1, final.edges$n2)
                         n1.side = ifelse(n1.tag<n2.tag, final.edges$n1.side, final.edges$n2.side)

                         n2 = ifelse(n1.tag<n2.tag, final.edges$n2, final.edges$n1)
                         n2.side = ifelse(n1.tag<n2.tag, final.edges$n2.side, final.edges$n1.side)

                         final.edges$n1 = n1
                         final.edges$n2 = n2
                         final.edges$n1.side = n1.side
                         final.edges$n2.side = n2.side

                         ## final merging of edges
                         final.edges = final.edges[, lapply(.SD, FUN), by = .(n1, n1.side, n2, n2.side)]

                         private$gGraphFromNodes(dt2gr(final.nodes, seqlengths = seqlengths(dnodes)), final.edges)

                         return(invisible(self))
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
                           return(invisible(self))

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
                         
                         ## label ref edges in current edges using refadj
                         setkeyv(edges, c("n1", "n1.side", "n2", "n2.side"))
                         edges[, ref := FALSE]
                         edges[refadj, ref := TRUE]

                         ## collect non ref edges
                         abed = edges[is.na(ref), ]

                         edcount = rbind(edges[, .(node = n1, side = n1.side)],
                                         edges[, .(node = n2, side = n2.side)]
                                         )[, .(count = .N), keyby = .(node, side)]

                         abcount = rbind(edges[ref==FALSE, .(node = n1, side = n1.side)],
                                         edges[ref==FALSE, .(node = n2, side = n2.side)]
                                         )[, .(count = .N), keyby = .(node, side)]



                         ## tabulate all edges (including "ref") 

                         ## now mark node sides that abut a non-reference adjacency
                         ## or have more than one adjacency emanating from them
                         nodes.dt[, ab.right := abcount[.(node.id, 'right'), ifelse(is.na(count), 0, count)]>0 |
                                      edcount[.(node.id, 'right'), ifelse(is.na(count), 0, count)]>1]

                         nodes.dt[, ab.left := abcount[.(node.id, 'left'), ifelse(is.na(count), 0, count)]>0 |
                                      edcount[.(node.id, 'left'), ifelse(is.na(count), 0, count)]>1]


                         ## also mark nodes that abut a loose end 
                         if (!ignore.loose)
                         {
                           nodes.dt[, ab.left := ab.left |  self$nodes$loose.left]
                           nodes.dt[, ab.right := ab.right | self$nodes$loose.right]
                         }
                         setkey(nodes.dt, node.id)

                         ## keep track of interval sides
                         sides = rbind(nodes.dt[, .(node.id, side = 'right', ab = ab.right)],
                                       nodes.dt[, .(node.id, side = 'left', ab = ab.left)])
                                                  
                         setkeyv(sides, c('node.id', 'side'))

                         ## an internal edge is a ref edge that connects two non-ab "sides"
                         edges[, n1.ab := sides[.(n1, n1.side), ab]]
                         edges[, n2.ab := sides[.(n2, n2.side), ab]]
                         edges[, internal := ref & !n1.ab & !n2.ab]

                         ## if "by" provided, then we take this also into account
                         ## to place additional constrain on internal edges
                         if (!is.null(by) && (by %in% names(nodes.dt)))
                           {
                             edges$n1.by = nodes.dt[.(edges$n1), ][[by]]
                             edges$n2.by = nodes.dt[.(edges$n2), ][[by]]
                             edges[, internal := internal & n1.by == n2.by]
                           }

                         ## now make a quick graph from only internal edges and find clusters
                         ## and their "ends", of which there is guaranteed to be one left
                         ## and one right due to our definitions above
                         tmp.gr = self$nodes$gr
                         tmp.gr$node.id.og = tmp.gr$node.id
                         igg = gG(nodes = tmp.gr, edges = edges[internal == TRUE, ])
                         igg$clusters('weak')
                         nodemap = data.table(node.id = igg$nodes$gr$node.id.og,
                                              new.node.id = igg$nodes$gr$cluster, key = 'node.id',
                                              left.end = igg$nodes$gr$loose.left,
                                              right.end = igg$nodes$gr$loose.right
                                              )


                         ## start a fresh dt just in case we overwrote
                         ## some user metadata columns
                         nodes.dt2 = gr2dt(nodes)
                         metadata.cols = setdiff(names(nodes.dt2),
                                                 c('seqnames', 'start', 'end', 'strand', 'width',
                                                   'snode.id', 'node.id', 'index', 'loose.left',
                                                   'new.node.id', 'left.end', 'right.end',
                                                   'loose.right'))

                         nodes.dt2 = cbind(nodemap[.(nodes.dt2$node.id), .(new.node.id, left.end, right.end)],
                                           nodes.dt2)


                         if (is.null(FUN)) ## if no fun specified just take the first 
                           FUN = function(x) x[1]

                         new.nodes = 
                           nodes.dt2[, .(seqnames = seqnames[left.end],
                                         start = start[left.end], end = end[right.end],
                                         node.id.left = node.id[left.end],
                                         node.id.right = node.id[right.end],
                                         loose.left = loose.left[left.end],
                                         loose.right = loose.right[right.end]), keyby = new.node.id]

                         if (length(metadata.cols)>0)
                         {
                           new.nodes = cbind(new.nodes, skrub(nodes.dt2)[, lapply(.SD, FUN), .SDcols = metadata.cols, keyby = new.node.id][.(new.node.id), metadata.cols, with = FALSE])
                         }

                         ## now find the old nodes that comprise the "sides" of the new nodes
                         ## and reconnect edges to these new ids
                         new.sides = rbind(
                           nodes.dt2[left.end == TRUE, .(node = node.id, side = 'left', new.node.id)],
                           nodes.dt2[right.end == TRUE, .(node = node.id, side = 'right', new.node.id)])
                         setkeyv(new.sides, c("node", "side"))


                         ## only edges that remain are external edges
                         new.edges = NULL

                         if (nrow(edges)>0)
                           {
                             new.edges = self$edges$dt[abs(sedge.id) %in% abs(edges$sedge.id[!edges$internal])]

                             if (nrow(new.edges)>0)
                             {
                               new.edges$n1 = new.sides[.(new.edges$n1, new.edges$n1.side), new.node.id]
                               new.edges$n2 = new.sides[.(new.edges$n2, new.edges$n2.side), new.node.id]
                             }
                           }

                         ## clean up new.nodes metadata
                         new.nodes$new.node.id = NULL
                         new.nodes$node.id.left = NULL
                         new.nodes$node.id.right = NULL

                         private$gGraphFromNodes(dt2gr(new.nodes, seqlengths = seqlengths(nodes)), new.edges)
                         return(invisible(self))

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
                                           verbose=FALSE
                                           )
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
                           return(invisible(self))
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

                         return(invisible(self))
                       },

                       #' @name eclusters
                       #' @description
                       #' Marks ALT edges belonging (quasi) reciprocal cycles 
                       #' @param juncs GRangesList of junctions
                       #' @param mc.cores parallel
                       #' @param ignore.strand usually TRUE
                       #' @return numerical vector of the same length, Inf means they r not facing each other
                       #' @author Marcin Imielinski
                       eclusters = function(thresh = 1e3, paths = TRUE,  mc.cores = 1, verbose = FALSE, chunksize = 1e30)
                       {                           
                         altedges = self$edges[type == "ALT", ]

                         bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]

                         ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
                         ixu = unlist(ix)
                         eps = 1e-9
                         ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
                         adj = Matrix::sparseMatrix(1, 1, x = FALSE, dims = rep(length(bp), 2))

                         if (verbose)

                           message(sprintf('Computing junction graph across %s ALT edges with distance threshold %s', length(altedges), thresh))                        
                         ## matrix of (strand aware) reference distances between breakpoint pairs
                         adj[ixu, ] = do.call(rbind, mclapply(ix,
                                                              function(iix)
                                                              {
                                                                if (verbose>1)
                                                                  cat('.')
                                                                tmpm = gr.dist(bp[iix], gr.flipstrand(bp), ignore.strand = FALSE)+eps
                                                                tmpm[is.na(tmpm)] = 0
                                                                tmpm[tmpm>thresh] = 0
                                                                tmpm = as(tmpm>0, 'Matrix')
                                                              },
                                                              mc.cores = mc.cores))

                         ## check bp pairs to see if they are actually reference connected (ignore.strand = TRUE)
                         ## on the given graphs ...
                         ## which if we have many graphs overlapping the same intervals
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

                         adj2 = adj & FALSE ## clear out adj for new skew symmetric version
                         adj2[junpos, junpos] = adj[bp2, bp1]
                         adj2[junpos, junneg] = adj[bp2, bp2]
                         adj2[junneg, junpos] = adj[bp1, bp1]
                         adj2[junneg, junneg] = adj[bp1, bp2]

                         if (verbose)
                           message(sprintf('Created basic junction graph using distance threshold of %s', thresh))

                         ## strongly connected components consists of (possibly nested) cycles
                         cl = split(1:length(bp), igraph::clusters(graph.adjacency(adj2), 'strong')$membership)

                         ## choose only clusters with length > 1
                         cl = cl[S4Vectors::elementNROWS(cl)>1]
                         cl = cl[order(S4Vectors::elementNROWS(cl))]


                         jcl = lapply(cl, function(x) unique(sort(bp$grl.ix[x])))
                         jcls = sapply(jcl, paste, collapse = ' ')
                         jcl = jcl[!duplicated(jcls)]
                         adj3 = adj2
                         altedges$mark(ecycle = as.character(NA))
                         if (length(jcl)>0)
                           {
                             dcl = dunlist(unname(jcl))[, listid := paste0('c', listid)]
                             altedges[dcl$V1]$mark(ecycle = dcl$listid)
                             altedges[dcl$V1]$mark(ecluster = dcl$listid)
                           }

                         if (verbose)
                           message(sprintf('Annotated %s junction cycles in edge field $ecycle', length(jcl)))                         
                         
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

                           altedges$mark(epath = as.character(NA))
                           if (length(jcl2)>0)
                             {
                               dcl2 = dunlist(unname(jcl2))[, listid := paste0('p', listid)]
                               altedges[dcl2$V1]$mark(epath = dcl2$listid)

                               ## also mark ecluster, though they may have
                               ## overlapping edges
                               self$edges$mark(ecluster =
                                                 ifelse(
                                                   is.na(self$edges$dt$ecycle) &
                                                   is.na(self$edges$dt$epath), as.character(NA),
                                                   paste0(
                                                     ifelse(is.na(self$edges$dt$ecycle),
                                                            '',
                                                            self$edges$dt$ecycle),
                                                     ifelse(is.na(self$edges$dt$epath),
                                                            '',
                                                            self$edges$dt$epath))))                              
                                 }

                           if (verbose)
                             message(sprintf('Annotated %s paths in edge field $epath', length(jcl2)))
                         }
                       },

                       #' @name paths
                       #' @description
                       #' Returns shortest paths from query gNode to subject gNode in graph in the form of gWalks
                       #' (note: the gNodes must exist in the graph, unlike in the related but more general proximity function)
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
                       #' By default performs a "cartesian" search, i.e. all pairs of query and subject
                       #' but if cartesian is set to FALSE will only search for specified pairs of query
                       #' and subject (then query and subject must be of the same length)
                       #' 
                       #' @param query gNode object or snode.id of gNode in this graph
                       #' @param subject gNode object or snode.id of gNode in this graph
                       #' @param weight edge column to use as the weight (instead of standard weight column)
                       #' @param mc.cores how many cores (default 1)
                       #' @return gWalk object each representing a query-subject shortest path (if any exist)
                       #' @author Marcin Imielinski
                       paths = function(query, subject = query, mc.cores = 1, weight = NULL, meta = NULL,
                                        ignore.strand = TRUE, cartesian = TRUE)
                       {
                         if (!inherits(query, 'gNode'))
                           query = self$nodes[query]

                         if (!inherits(subject, 'gNode'))
                           subject = self$nodes[subject]

                         if (!cartesian && length(query) != length(subject))
                         {
                           stop('If cartesian = FALSE then query and subject must have the same length')
                         }

                         ## get igraph whose edge weights correspond to the width of the target interval
                         G = self$igraph

                         if (length(query)==0 | length(subject)==0)
                           return(gW(graph = self))

                         if (!is.null(weight)) ## replace standard weight with user provided column
                         {
                           weight = weight[1]
                           if (!(weight %in% names(self$sedgesdt)))
                             stop(paste('column', weight, 'not found in edges metadata'))

                           E(G)$weight = self$sedgesdt[.(E(G)$sedge.id),weight, with = FALSE][[1]]
                         }

                         ## only takes a length 1 query, since igraph::get.shortest.paths only supports scalar queries
                         ## for the source, so then we will have to iterate this query on all the input query nodes
                         ## do.trim is there to do.trim the last interval from edge weight is we are using the
                         ## plain vanilla edge weights i.e. where each edge is weighted by the width of its target
                         .path = function(qid, G, query, subject, ignore.strand = TRUE, do.trim = TRUE, cartesian = TRUE)
                         {
                           query = query[qid, with = FALSE]

                           if (!cartesian) ## only choose one of the subject intervals to run
                           {
                             subject = subject[qid, with = FALSE]
                           }
                           
                           .get.shortest.paths = function(G, query, subject, do.trim = TRUE)
                           {
                             snode.id = self$gr$snode.id

                             p = suppressWarnings(igraph::get.shortest.paths(G, query, subject,
                                                                             weight = igraph::E(G)$weight, output = 'both'))
                             
                             ## convert to integer .. and back to snode.id
                             ## if path does not exist get.shortest.paths (for some reason) will return a single vector with the target
                             ## so we clean out all paths that don't have the query                             
                             v  = lapply(p$vpath, function(x) {y = as.integer(x); if (y[1] != query) return(c()) else return(snode.id[y])})
                             ## chop off the last interval to get the "true" distance from query to subject
                             if (do.trim)
                               {
                                 w.all = sapply(p$epath, function(x) sum(x$weight[-length(x)])+1)
                               }
                             else
                               {
                                 w.all = sapply(p$epath, function(x) sum(x$weight))
                               }

                             w.all[elementNROWS(v)==0] = Inf
                             
                             return(data.table(source = snode.id[query], sink = snode.id[subject], path = v, dist = w.all)[!is.infinite(dist), ])
                           }

                           p.dt = .get.shortest.paths(G, query$dt$index, subject$dt$index, do.trim)
                           
                           ## compute the three other flavors of paths - i.e. ++, +-, -+, ++
                           if (ignore.strand)
                           {
                             p.fr = .get.shortest.paths(G, query$dt$index, subject$flip$dt$index, do.trim)
                             p.rr = .get.shortest.paths(G, query$flip$dt$index, subject$flip$dt$index, do.trim)
                             p.rf = .get.shortest.paths(G, query$flip$dt$index, subject$dt$index, do.trim)
                             p.dt = rbind(p.dt, p.fr, p.rr, p.rf)[order(dist), ][!duplicated(abs(sink)), ]
                           }
                           return(p.dt)
                         }

                         do.trim = is.null(weight)

                         
                         walks.dt = rbindlist(mclapply(1:length(query), .path, query = query,
                                                       subject = subject, G = G, ignore.strand = ignore.strand,
                                                       cartesian = cartesian,
                                                       mc.cores = mc.cores, do.trim = do.trim))


                         walks = gW(snode.id = walks.dt$path, graph = self, meta = walks.dt[, setdiff(colnames(walks.dt), 'path'), with = FALSE])
                         return(walks)
                       },
  
                       #' @name dist
                       #' @description
                       #' Computes a distance matrix of query and subject
                       #' intervals (in base pairs) on the gGraph
                       #' between any arbitrary pairs of granges gr1 and gr2.
                       #' @param gr1 GRanges query 
                       #' @param gr2 GRanges query (if NULL, will set to gr1)
                       #' @param weight metadata field of gEdges to use as weight (instead of distance of target node)
                       #' @param dt returns data.table if TRUE, excluding all Inf distances
                       #' @param include.internal logical flag whether to allow
                       #' paths that begin or end inside teh query or subject
                       #' @author Marcin Imielinski
                       dist = function(query,
                                       subject,
                                       weight = NULL,
                                       ignore.strand = TRUE,
                                       include.internal = TRUE,
                                       verbose=FALSE)
                       {
                         if (!missing('query'))
                         {
                           if (inherits(query, 'gNode'))
                             {
                               queryf = query
                             }
                           else
                             {
                               ## first try evaluating query and subject as node indices or filters
                               queryf = tryCatch(self$nodes[query], error = function(e) NULL)
                             }
                         }

                         ## primitive function to compute distance on graph, called by
                         ## both gNode and GRanges queries below 
                         .dist = function(gg, query, subject, weight = NULL, ignore.strand = TRUE)
                         {
                           ## get igraph and (populated with edge weights)
                           G = gg$igraph
                           
                           if (!is.null(weight)) ## replace standard weight with user provided column if provided
                           {
                             weight = weight[1]
                             if (!(weight %in% names(gg$sedgesdt)))
                               stop(paste('column', weight, 'not found in edges metadata'))
                             
                             igraph::E(G)$weight = gg$sedgesdt[.(E(G)$sedge.id),weight, with = FALSE][[1]]
                           }

                           qmap = data.table(
                             query.ix = query$dt$index,
                             query.rix = query$flip$dt$index)
                           
                           smap = data.table(
                             subject.ix = subject$dt$index,
                             subject.rix = subject$flip$dt$index)
                           
                           ## shortest.paths does not like dups so we have to dedup here
                           ## will redup later
                           qmap[, dup := duplicated(query.ix)]
                           smap[, dup := duplicated(subject.ix)]

                           query.ix = qmap[dup == FALSE, query.ix]
                           query.rix = qmap[dup == FALSE, query.rix]
                           
                           subject.ix = smap[dup == FALSE, subject.ix]
                           subject.rix = smap[dup == FALSE, subject.rix]

                           Dff = igraph::shortest.paths(G, query.ix, subject.ix,
                                                        weights = igraph::E(G)$weight,
                                                        mode = 'out')
                           if (ignore.strand)
                           {
                             Dfr = igraph::shortest.paths(G, query.ix, subject.rix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             Drr = igraph::shortest.paths(G, query.rix, subject.rix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')
                             
                             Drf = igraph::shortest.paths(G, query.rix, subject.ix,
                                                          weights = igraph::E(G)$weight,
                                                          mode = 'out')                             
                             D = pmin(Dff, Dfr, Drr, Drf)
                           }
                           else
                             D = Dff

                           ## remap to original values
                           rownames(D) = query.ix
                           colnames(D) = subject.ix

                           D = D[as.character(qmap$query.ix), as.character(smap$subject.ix), drop = FALSE]

                           if (is.null(weight)) ## we only need to sweep for the standard usage, but not for a custom weight
                             {
                               D = pmax(sweep(D, 2, width(subject)-1, '-'), 0)
                             }

                           return(D)
                         }

                         ## we have gNodes or gNode indices as arguments
                         if (!is.null(queryf))
                         {
                           if (missing('subject'))
                             subjectf = queryf
                           else
                           {
                             if (inherits(subject, 'gNode'))
                             {
                                 subjectf = subject
                               }
                             else
                             {
                                 subjectf = tryCatch(self$nodes[subject], error = function(e) NULL)
                               }
                           }

                           if (is.null(subjectf))
                             stop('Error in parsing query or subject filters')

                           D = .dist(self, queryf, subjectf, ignore.strand = ignore.strand, weight = weight)
                           return(D)
                         }
                         ## we have GRanges as arguments
                                                 
                         if (!ignore.strand) ##FIXME
                         {
                           stop('Stranded distance not yet implemented for $dist method with GRanges arguments')
                         }

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

                         ed = self$edgesdt[,.(n1,n2,n1.side,n2.side,type)]

                         if (!is.null(weight))
                           {
                             ed = cbind(ed, self$edgesdt[, weight, with = FALSE])
                           }

                         simpleg = gG(nodes = self$nodes$gr[, c()],
                                      edges = ed)$simplify()
                         
                         simpleg$nodes$mark(simpleg = TRUE)
                         
                         ## split query and subject if we want to allow
                         ## paths that originate or end inside a node
                         ## to constitute a "distance"
                         if (include.internal)
                         {
                           query = query %*% simpleg$nodes$gr[, c()]
                           subject = subject %*% simpleg$nodes$gr[, c()]
                         }
                         
                         grd = unique(gr.stripstrand(grbind(query, subject)))

                         ## this will disjoin by grd without collapsing
                         ## and only keep the nodes and edges that were previously in simpleg
                         simpleg$disjoin(gr = grd, collapse = FALSE)
                         simpleg = simpleg[which(simpleg==TRUE), ]

                         ## if custom weight field is specified then we will have
                         ## created many possibly internal edges that we have to assume
                         ## be zero weight
                         if (!is.null(weight))
                         {
                           wval = simpleg$sedgesdt[[weight]]
                           sed = simpleg$sedgesdt
                           simpleg$annotate(weight, 0,
                                            sed$sedge.id[is.na(wval)], "edge")
                         }

                         grd$is.query = grd %^% query
                         grd$is.subject = grd %^% subject

                         grds = grd %*% simpleg$nodes$gr[, 'snode.id']
                         values(grds) = cbind(values(grds),
                                              simpleg$queryLookup(grds$snode.id)[, .(index, rindex)])

                         qnode = simpleg$nodes[grds$snode.id[grds$is.query]]
                         snode = simpleg$nodes[grds$snode.id[grds$is.subject]]

                         D = .dist(simpleg, qnode, snode,
                                   ignore.strand = TRUE, weight = weight)

                         ## now we melt D and merge back with input
                         ## using overlaps of query / subject with grds
                         Dt = as.data.table(reshape2::melt(D))[!is.infinite(value), ]
                         setkeyv(Dt, c("Var1", "Var2"))

                         queryl = query[, 'id'] %*% grds[, 'index']
                         subjectl = subject[, 'id'] %*% grds[, 'index']

                         out.dt = merge(gr2dt(queryl)[, .(id, qindex = index)],
                                        merge(gr2dt(subjectl)[, .(id, sindex = index)],
                                              Dt, by.x = 'sindex', by.y = 'Var2',
                                              allow.cartesian = TRUE),
                                        by.x = 'qindex', by.y = 'Var1',
                                        allow.cartesian = TRUE)

                         setkey(out.dt, value) ## sorts by value i.e. distance
                         out.dt = unique(out.dt, by = c("id.x", "id.y"))

                         D = matrix(Inf, nrow = length(query.og), ncol = length(subject.og))
                         D[cbind(out.dt$id.x, out.dt$id.y)] = out.dt$value


                         return(D)
                       },

                       #' @name rep
                       #' @description
                       #'
                       #' Creates "bubbles" in the graph by replicating the nodes or gwalks in the argument.  Node replication
                       #' replicates edges going in and out of all replicated nodes.  If an edge connects a pair of replicated nodes
                       #' that edge will be replicated across all pairs of those replciated nodes.   Walk replication will create "longer bubbles"
                       #' with fewer edges getting replicated i.e. it will only replicate intra-walk edges within each walk replicate (but not between
                       #' separate walk replicates).
                       #'
                       #' (note that this changes the current gGraph in place, and thus the input gNode or gWalk will no longer
                       #' apply to the new altered gWalk)
                       #'
                       #' New graph keeps track of the parent node and edge ids in the original graph using node metadata parent.node.id
                       #' and edge metadata parent.edge.id i.e. the replicated nodes will be connected to the sources of the original nodes
                       #' and if replicated nodes connect to each other, then there will exist an edge connecting
                       #' all of their instances to each other.
                       #'
                       #' @param nodes = gNode object must point to a node in the graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                       #' @param times  scalar or vector of length self$length specifying how many times to replicate each of the nodes.
                       #' @author Marcin Imielinski
                       #' @return gGraph (also modified in place) with nodes annotated with parent.node.id, parent.rep
                       rep = function(nodes = NULL, times)
                       {
                         ## by default i.e. if nodes not specified assume times applies to the whole graph
                         if (is.null(nodes))
                           {
                             nodes = self$nodes
                           }
                         ## convert everything to integer list
                         else if (is(nodes, 'gNode'))
                         {
                           if (!(nodes$graph == self))
                             {
                               stop('gWalk input must point to the nodes of this gGraph')
                             }
                           nlist = split(nodes$dt$snode.id, 1:length(nodes))
                         }
                         else if (is(nodes, 'gWalk'))
                         {
                           if (!(nodes$graph == self))
                             {
                               stop('gWalk input must point to this gGraph')
                             }
                           nlist = nodes$snode.id
                         }
                         else if (is.integer(nodes) | is.numeric(nodes))
                         {
                           nlist = split(nodes, 1:length(nodes))
                         }
                         else ## otherwise user provided list of integer node indices as input 
                         {
                           nlist = nodes
                         }

                         if (length(times) != 1 && length(nlist) != length(times))
                         {
                           stop('Length of times must be either length 1 or length of nodes')
                         }

                         times = cbind(1:length(nlist), times)[,2]

                         dtu = dunlist(unname(nlist))
                         
                         setnames(dtu, c('walk.id', 'node.id'))
                         dtu$times = times[dtu$walk.id]
                         dtu[, node.id := abs(node.id)]

                         if (!is.integer(dtu$node.id) && any(is.na(as.integer(dtu$node.id))))
                           {
                             stop('Raw node indices (or lists of node indices) must be integers')
                           }

                         if (any(dtu$V1>length(self)))
                         {
                           stop('One or more provided node indices are outside the bounds of the graph')
                         }

                         gr = self$nodes$gr
                         ed = self$edges$dt
                         ed[, parent.edge.id := edge.id]
                         gr$parent.node.id = gr$node.id

                         
                         other = setdiff(1:length(self), dtu$node.id)

                         ## parent.node.id = is the node.id in the original graph
                         ## walk.id = is the original walk.id
                         ## walk.rep = is a local id unique instance of each walk replicate inside each walk.id
                         nmap = data.table(parent.node.id = c(other, rep(dtu$node.id, dtu$times)),
                                           walk.id = c(rep(NA, length(other)), paste(rep(dtu$walk.id, dtu$times))),
                                           walk.rep = c(rep(NA, length(other)), data.table(ind = 1:nrow(dtu), time = dtu$times)[, .(1:time), by = ind][, V1])
                                           )[, node.id := 1:.N]
                                                  
                         new.gr = gr[nmap$parent.node.id]
                         new.gr$parent.rep = nmap$walk.rep

                         ## now we only replicate edges if parent.node.id match and either
                         ## (1) walk.id of either node is NA
                         ## (2) walk.ids and walk.rep match
                         ## (3) walk.ids are non NA and they don't match
                         ## first we merge everything (may create a ton of edges
                         ed$parent.edge.id = abs(ed$edge.id)
                         new.ed = merge(merge(ed, nmap, by.x = 'n1', by.y = 'parent.node.id', allow.cartesian = TRUE),
                                        nmap, by.x = 'n2', by.y = 'parent.node.id', allow.cartesian = TRUE)

                         new.ed = new.ed[which(is.na(walk.id.x) | is.na(walk.id.y) | ## condition 1
                                               (walk.id.x != walk.id.y) | ## condition 2
                                               (walk.id.x == walk.id.y & walk.rep.x == walk.rep.y) ## condition 3
                                               ), ]
                         new.ed[, n1 := node.id.x]
                         new.ed[, n2 := node.id.y]

                         new.ed$walk.rep.y = new.ed$walk.rep.x = new.ed$walk.id.y = new.ed$walk.id.x = new.ed$node.id.x = new.ed$node.id.y = NULL

                         private$gGraphFromNodes(new.gr, new.ed)
                         return(invisible(self))
                       },

                       #' @name swap
                       #' @description
                       #'
                       #' Swap nodes with granges, grl, or Gwalks.
                       #' Provided replacement vector must be the same length as the inputted nodes, resulting in each node being "swapped" by the provided 
                       #' interval, node, grl (representing a walk), or gWalk.  The replacement will inherit left and right edges for the removed node.
                       #' If the replacement is a walk, then the left side of the first node in the walk will inherit the edges that were
                       #' previously to the left of the node being replaced, and right side of the last node of the walk will inherit the edges
                       #' that were previously to the right of the node being replaced.  
                       #'
                       #' Note: these replacement obey the orientation of the arguments.  So if the node to be replaced is flipped (- orientation with
                       #' respect to the reference, then it's "left" is to the right on the reference.  Similarly for walks whose first interval
                       #' is flipped with respect to the reference, the left edges will be attached to the right of the node on the reference.
                       #' 
                       #' @param nodes = gNode object must point to a node in this graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                       #' @param replacement  GRanges, GRangesList, or gWalk object whose length is the length(nodes)
                       #' @author Marcin Imielinski
                       #' @return gGraph (also modified in place) with nodes annotated with parent.node.id, parent.rep
                       swap = function(nodes, replacement)
                       {
                         if (length(nodes) != length(replacement))
                           stop('nodes to be swapped must be the same length as their replacement')

                         if (is.integer(nodes) | is.numeric(nodes))
                           nodes = self$nodes[nodes]

                         if (nodes$graph != self)
                           stop('nodes must be matched to the current graph')

                         if (inherits(replacement, 'GRangesList'))
                           {
                             replacement = gW(grl = replacement, disjoin = FALSE)

                             ## automatically label all inserted edges as ALT by default
                             replacement$graph$edges$mark(type = 'ALT')
                           }

                         ## at this point should be either gRanges or gWalk
                         if (inherits(replacement, 'GRanges'))
                         {
                           nodes.rep = replacement
                           nodes.rep$parent.node.id = nodes$dt$node.id
                           nodes.rep$parent.rep = 1:length(replacement)
                           ed.rep = data.table() ## there are no edges in this case
                         }
                         else if (inherits(replacement, 'gWalk')) 
                         {
                           nodes.rep = replacement$graph$nodes$gr
                           tmp = dunlist(unname(replacement$snode.id))
                           nodes.rep$parent.node.id[abs(tmp$V1)] = nodes$dt$node.id[tmp$listid]
                           nodes.rep$parent.rep = tmp$listid
                           ed.rep = replacement$graph$edges$dt[, n1 := paste('replacement', n1)][, n2 := paste('replacement', n2)]
                         }

                         ## relabel the nodes that we are swapping into the graph
                         nodes.rep$node.id = paste('replacement', nodes.rep$node.id)


                         ## extract the node GRanges from the current graph, we will swap a subset of these 
                         nodes.old = self$nodes$gr
                         nodes.old$parent.rep = NA

                         ## extract the edges
                         ed = self$edges$dt[, n1 := as.character(n1)][, n2 := as.character(n2)]
                                                  
                         ## these vectors will be all of the same length
                         swap = data.table(
                           toswap = nodes$dt$node.id,
                           firsts = sapply(replacement$snode.id, '[', 1),
                           lasts = sapply(replacement$snode.id, function(x) x[length(x)]),
                           key = 'toswap'
                         )

                         eswap.left = nodes$eleft$dt ## these edges will have the node to be swapped as n2
                         eswap.right = nodes$eright$dt ## these edges will have the node to be swapped as n1

                         ## firsts get left edges  attached to their left, unless they are flipped (i.e. have negative snode.id)
                         eswap.left[, n2.swap := swap[.(n2), paste('replacement', abs(firsts))]]
                         eswap.left[, n2.swap.side := swap[.(n2), ifelse(firsts>0, 'left', 'right')]]

                         ## lasts get right edges attached to their right, unless they are flipped (i.e. have negative snode.id)
                         eswap.right[, n1.swap := swap[.(n1), paste('replacement', abs(lasts))]]
                         eswap.right[, n1.swap.side := swap[.(n1), ifelse(lasts>0, 'right', 'left')]]

                         ## ed1 only has edges in "forward" orientation so need to make sure to swap the correct "end" of each edge
                         ## we will use the sign to guide us - i.e. we swap n2 for eswap.left, unless they are flipped, in which
                         ## case we swap n1
                         ed[.(abs(eswap.left$sedge.id)), ":="(
                                                           n1 = ifelse(eswap.left$sedge.id<0, eswap.left$n2.swap, n1),
                                                           n1.side = ifelse(eswap.left$sedge.id<0, eswap.left$n2.swap.side, n1.side),
                                                           n2 = ifelse(eswap.left$sedge.id>0, eswap.left$n2.swap, n2),
                                                           n2.side = ifelse(eswap.left$sedge.id>0, eswap.left$n2.swap.side, n2.side))]

                         ## for eswap.right we do the reverse, i.e. forward orientation edges get their n1 switched, and n2 if flipped
                         ed[.(abs(eswap.right$sedge.id)), ":="(
                                                           n1 = ifelse(eswap.right$sedge.id>0, eswap.right$n1.swap, n1),
                                                           n1.side = ifelse(eswap.right$sedge.id>0, eswap.right$n1.swap.side, n1.side),
                                                           n2 = ifelse(eswap.right$sedge.id<0, eswap.right$n1.swap, n2),
                                                           n2.side = ifelse(eswap.right$sedge.id<0, eswap.right$n1.swap.side, n2.side))]

                         ## if there are additional internal edges (i.e. inside the provided gWalk) then we add these to the new graph
                         ## along with any metadata 
                         if (nrow(ed.rep)>0) 
                           {
                             ed = rbind(ed, ed.rep, fill = TRUE)
                             ed$sedge.id = ed$edge.id = NULL
                           }

                         ## by now all instances of swap$toswap node.ids should be missing from ed
                         ## and replaced with nodes from n2
                         ## we now concatenate the remaining n1 with n2, find a common
                         nodes.new = grbind(nodes.old[-swap$toswap, ], nodes.rep)
                         ed[, n1 := as.integer(factor(n1, nodes.new$node.id))]
                         ed[, n2 := as.integer(factor(n2, nodes.new$node.id))]

                         ## remove previous node.id and snode.id so as not to confuse instantiator 
                         nodes.new$node.id = NULL
                         nodes.new$snode.id = NULL

                         ## reinstantiate 
                         private$gGraphFromNodes(nodes.new, ed)

                         return(invisible(self))
                       },

                       #' @name connect
                       #' @description
                       #'
                       #' Connect node pairs in the gGraph by adding (optional) edge metadata and (optionally) inserting nodes or grl / walks
                       #' in between the given edge.  Note: the connections are made with respect to the provided node orientation
                       #' so if the node is provided in a "flipped" orientation then it's right direction will point left on the reference.
                       #' 
                       #' @param n1 = gNode object must point to a node in this graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                       #' @param n2 = gNode object must point to a node in this graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                       #' @param n1.side character vector of length length(n1) whose value is either "left" or "right" (default 'right')
                       #' @param n2.side character vector of length length(n1) whose value is either "left" or "right" (default 'left')
                       #' @author Marcin Imielinski
                       #' @return gGraph (also modified in place) with nodes annotated with parent.node.id, parent.rep
                       connect = function(n1, n2, n1.side = 'right', n2.side = 'left', type = 'ALT', meta = NULL, insert = NULL)
                       {
                         if (length(n1) != length(n2))
                           stop('length of n1 and n2 must be the same')

                         if (inherits(n1, 'gNode'))
                         {
                           if (n1$graph != self)
                             stop('n1 must match gGraph')

                           n1 = n1$dt$snode.id
                         }

                         if (inherits(n2, 'gNode'))
                         {
                           if (n2$graph != self)
                             stop('n2 must match gGraph')

                           n2 = n2$dt$snode.id
                         }

                         newedges = data.table(n1 = abs(n1), n2 = abs(n2),
                                               n1.side = ifelse(sign(n1)*sign((n1.side == 'right')-0.5)>0, 'right', 'left'),
                                               n2.side = ifelse(sign(n2)*sign((n2.side == 'right')-0.5)>0, 'right', 'left'))


                         if (!is.null(meta))
                           newedges = cbind(newedges, meta)

                         if (is.null(newedges$type))
                           newedges$type = type

                         newnodes = self$nodes$gr
                         newnodes$parent.node.id = newnodes$node.id
                         newnodes$parent.rep = NA

                         ## if insert is provided then we add more nodes and edges and reconnect them to the provided 
                         if (!is.null(insert))
                         {
                           if (length(insert) != nrow(newedges))
                             stop('length of insert must be equal to the number of inserted edges')

                           if (inherits(insert, 'GRanges'))
                             insert = split(insert, 1:length(insert))
                           
                           if (inherits(insert, 'GRangesList'))
                           {
                             insert = gW(grl = insert, disjoin = FALSE)
                             
                             ## automatically label all inserted edges as ALT by default
                             insert$graph$edges$mark(type = 'ALT')
                           }

                           ## at this point insert is a gWalk

                           insert.gr = insert$graph$nodes$gr
                           tmp = dunlist(unname(insert$snode.id))
                           insert.gr$parent.node.id[abs(tmp$V1)] = insert$dt$node.id[tmp$listid]
                           insert.gr$parent.rep = tmp$listid
                           ed.insert = insert$graph$edges$dt[, n1 := paste('insert', n1)][, n2 := paste('insert', n2)]

                           ## relabel the insert that we are swapping into the graph
                           insert.gr$node.id = paste('insert', insert.gr$node.id)

                           ## these vectors will be all of the same length
                           inodes = data.table(                             
                             firsts = sapply(insert$snode.id, '[', 1),
                             lasts = sapply(insert$snode.id, function(x) x[length(x)])
                           )

                           ## we split newedges into 2
                           newedges1 = copy(newedges)
                           newedges2 = copy(newedges)

                           ## connect n2 side of newedges1 to (the correct side) of firsts
                           ## i.e. if firsts is + then left side otherwise right side
                           ## note: this does not depend on the provided n1.side or n2.side 
                           newedges1[, n2 := insert.gr$node.id[abs(inodes$firsts)]]
                           newedges1[, n2.side := ifelse(inodes$firsts>0, 'left', 'right')]

                           ## connect n1 of newedges2 to (the correct side) of lasts
                           ## i.e. if first is + then right side otherwise left side 
                           newedges2[, n1 := insert.gr$node.id[abs(inodes$lasts)]]
                           newedges2[, n1.side := ifelse(inodes$lasts>0, 'right', 'left')]

                           newedges = rbind(newedges1, newedges2, ed.insert, fill = TRUE)
             
                           ## by now all instances of insert$toinsert node.ids should be missing from ed
                           ## and replaced with nodes from n2
                           ## we now concatenate the remaining n1 with n2, find a common
                           newnodes$node.id = as.character(newnodes$node.id)

                           ## reset the loose end status of the inserted walk node GRanges
                           insert.gr$loose.right = insert.gr$loose.left = FALSE                           
                           newnodes = grbind(newnodes, insert.gr)

                           ## factor across new levels and convert to integers to include the newnodes
                           newedges[, n1 := as.integer(factor(n1, newnodes$node.id))]
                           newedges[, n2 := as.integer(factor(n2, newnodes$node.id))]
                         }
                         
                         
                         if (nrow(self$edges$dt)>0)
                           newedges = rbind(self$edges$dt, newedges, fill = TRUE)

                         newnodes$node.id = newnodes$snode.id = NULL
                         newedges$edge.id = newedges$sedge.id = NULL

                         private$gGraphFromNodes(newnodes, newedges)

                         return(invisible(self))
                       },

                       #' @name toposort
                       #'
                       #' marks node metadata with topo.order after topological sort
                       #' if graph has cycles, will do partial sort, leave some nodes as NA
                       #'
                       #' @return current graph, modified with nodes marked according to topological sort
                       toposort = function()
                       {
                         v = unique(abs(self$gr$snode.id[as.vector(
                                         suppressWarnings(
                                           igraph::topo_sort(self$igraph)))]))

                         self$nodes$mark(topo.sort = NULL)
                         self$nodes[v]$mark(topo.order = 1:length(v))
                         return(invisible(self))
                       },

                       #' @name gGraph$print
                       #' @description
                       #'
                       #' Prints out this gGraph. Prints number of nodes and edges, the gNode associated
                       #' with this gGraph and the gEdge associated with this gGraph
                       print = function()
                       {
                         nl = length(self$loose)
                         nt = length(self$terminal)
                         nalt = sum(private$pedges$type == "ALT")/2
                         nref = sum(private$pedges$type == "REF")/2
                         ne = nalt+nref

                         message(sprintf('gGraph with %s nodes, %s loose ends (%s terminal and %s internal), and %s edges (%s REF and %s ALT)',
                                         self$length, nl, nt, nl-nt, ne, nref, nalt),
                                 appendLF = FALSE)

                         if (self$nodes$length > 0){
                             message(' comprising:\n')
                             }
                         else
                           message('\n')
                         self$nodes$print()
                         
                         if (nrow(private$pedges))
                         {
                           message()
                           self$edges$print()
                         }
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
                         if (class == "node") {
                           NONO.FIELDS = c('node.id', 'snode.id', 'index', 'loose.left', 'loose.right', 'loose.left', 'loose.right')
                           if (colName %in% NONO.FIELDS)
                             stop(paste('Cannot alter these protected gNode fields: ', paste(NONO.FIELDS, collapse = ', ')))

                           if (is.null(data))
                             {
                               values(private$pnodes)[[colName]] = NULL
                               return(invisible(self))
                             }

                           id = self$queryLookup(id)
                           id$data = data

                           index = c(id[, index], id[, rindex])
                           data = c(id[, data], id[, data])
                          
                           gr.dt = gr2dt(private$pnodes)
                           gr.dt[index, paste(colName) := data]
                           values(private$pnodes)[[colName]] = gr.dt[[colName]]
                           
                         } else if (class == "edge") {
                           NONO.FIELDS = c('from', 'to', 'sedge.id', 'edge.id','n1', 'n2', 'n1.side', 'n2.side')
                           if (colName %in% NONO.FIELDS)
                             stop(paste('Cannot alter these protected gEdge fields: ', paste(NONO.FIELDS, collapse = ', ')))
                           
                           if (is.null(data))
                           {
                             private$pedges[[colName]] = NULL
                             return(invisible(self))
                           }

                           if (colName == "type" && (!is.character(data) || !all(data %in% c('REF', 'ALT'))))
                             stop('type is a reserved gEdge metadata field and can only be replaced with values REF and ALT')

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

                           ed[, ":="(lwd = ifelse(is.na(lwd), ifelse(type=="ALT", log2(0.2*pmax(0, y, na.rm = TRUE)+2)+1, ifelse(type == 'loose', 3, 1)),lwd),
                                     lty = ifelse(is.na(lty), ifelse(type=='loose', 1, 1),lty),
                                     col = ifelse(is.na(col), ifelse(is.na(col), ifelse(type=="ALT",
                                                                                 ifelse(is.na(y),
                                                                                        alpha("red", 0.4),
                                                                                        ifelse(y>0,
                                                                                               alpha("red", 0.4),
                                                                                               alpha("purple", 0.3))),
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
                                                             'gr.labelfield', 
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

                         if (is.null(gt.args$stack.gap)){
                           gt.args[['stack.gap']] = 0
                         }

                         stack.gap = gt.args$stack.gap

                         y.field = gt.args$y.field

                         if (is.null(stack.gap)){
                           stack.gap = 1e5
                           }

                         if (is.null(gt.args$angle)){
                           gt.args[['angle']] = 0
                         }

                         if (is.null(gt.args$y0)){
                           gt.args[['y0']] = 0
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
                           gt.args$y.field = 'y'

                           ## let loose ends inherit their y value from their parent nodes
                           if (any(lix <- ss$loose))
                           {
                             ss$y[lix] = ss$y[ss$node.id[lix]] + 0.25
                           }                          
                           gt.args[['data']] = unname(ss)                               
                           gt = do.call(gTrack::gTrack, gt.args)
                         } else {

                           gt.args$y.field = "y"                           
                           gt.args$yaxis = FALSE

                           ## stack node pairs via stack.gap
                           tmp.ss = ss[ss$snode.id>0]
                           tmp.ss$y = suppressWarnings(disjointBins(tmp.ss+stack.gap))-0.5
                           ss$y = tmp.ss$y[match(ss$node.id, tmp.ss$node.id)] 

                           ## stack parent graphs
                           ## if the data here is the result of a
                           ## concatenation
                           if (!is.null(ss$parent.graph) | !is.null(ss$parent.rep))
                           {
                             if (is.null(ss$parent.graph))
                               {
                                 ss$parent.graph = ss$parent.rep
                               }
                             else if (!is.null(ss$parent.rep))
                               {
                                 ss$parent.graph = paste(ss$parent.graph, ss$parent.rep)
                               }

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

                            ## let loose ends inherit their y value from their parent nodes
                           if (any(lix <- ss$loose))
                             ss$y[lix] = ss$y[ss$node.id[lix]] + 0.2
                               
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
                         if (length(tile) == 0) {
                           stop("tile cannot contain no nodes, nothing to trim around")
                         }
                         if(length(self) == 0) {
                           return(invisible(self))
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
                         return(invisible(self))
                       },
                       
                       
                       #' @name add
                       #' @description
                       #'
                       #' Adds GRanges nodes, edges (data.table), or junctions to graph
                       #' Only one of the below parameters can be specified at a time (since
                       #' the graph is modified in place, order matters)
                       #'
                       #' @param nodes GRanges, strand is ignored
                       #' @param edges data.table specifying edges in existing table with field n1, n2, n1.side, n2.side
                       #' @param junctions Junction object or GRangesList coercible to junction object
                       #' @author Marcin Imielinski
                       #' @return current graph modified in place with additional nodes and edges, as specified by user 
                       add = function(nodes = NULL, edges = NULL, junctions = NULL)
                       {
                         if (sum(c(!is.null(nodes), !is.null(edges), !is.null(junctions)))!=1)
                           stop('Exactly one of the following can be added in a given call to add: nodes, edges, junctions')
                         if (!is.null(nodes))
                         {
                           newnodes = grbind(self$nodes$gr, nodes)
                           private$gGraphFromNodes(nodes = newnodes,
                                                   edges = self$edgesdt)
                           return(invisible(self))
                         }
                         else if (!is.null(edges))
                         {
                           if (!is.data.table(edges))
                             edges = as.data.table(edges)

                           if (!all(c('n1', 'n2', 'n1.side', 'n2.side') %in% names(edges)))
                           {
                             stop('edges must be data.table with fields n1, n2, n1.side, n2.side')
                           }

                           if (any(duplicated(names(edges))))
                           {
                             stop('edges has duplicated column names')
                           }

                           if (is.null(edges$type))
                             edges$type = 'ALT'

                           if (any(iix <- is.na(edges$type)))
                             edges[iix, type := 'ALT']

                           newedges = rbind(self$edgesdt, edges, fill = TRUE)
                           private$gGraphFromNodes(nodes = self$nodes$gr,
                                                   edges = newedges)
                           return(invisible(self))
                         }
                         else if (!is.null(junctions))
                         {
                           if (!inherits(junctions, 'Junction'))
                             junctions = jJ(grl = junctions)

                           bp = grl.unlist(junctions$grl)
                           self$disjoin(bp[, c()], collapse = FALSE)
                           ov = bp[, c('grl.ix', 'grl.iix')] %*% self$nodes$gr[, 'node.id']
                           strand(ov) = strand(bp)[ov$query.id]
                           ov1 = ov[ov$grl.iix == 1]
                           ov2 = ov[ov$grl.iix == 2]

                           newedges =
                             merge(gr2dt(ov1),
                                   gr2dt(ov2),
                                   by = 'grl.ix', allow.cartesian = TRUE)[, ":="(
                                                  n1 = node.id.x,
                                                  n1.side = ifelse(strand.x == '+', 'left', 'right'),
                                                  n2 = node.id.y,
                                    n2.side = ifelse(strand.y == '+', 'left', 'right'))][, .(n1, n1.side, n2, n2.side, grl.ix, type = 'ALT')]

                           if (ncol(junctions$dt)>0)
                             {
                               newedges = cbind(newedges, junctions$dt)
                               newedges = newedges[, unique(colnames(newedges)), with = FALSE]
                             }

                           private$gGraphFromNodes(nodes = self$nodes$gr,
                                                   edges = rbind(self$edgesdt, newedges, fill = TRUE))
                           return(invisible(self))                           
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
                       
                       ## timestamp of object helps us identify when the node edge / structure was last updated
                       ptimestamp = NULL,
                       
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

                       stamp = function() private$ptimestamp = tstamp(),

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
                           return(invisible(self))
                         }
                         
                         lt = match(private$pnodes$snode.id, -private$pnodes$snode.id)
                         private$lookup = data.table(snode.id = private$pnodes$snode.id,
                                                     index = seq_along(private$pnodes),
                                                     rindex = lt)
                         setkey(private$lookup, snode.id)
                         
                         return(invisible(self))
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

                         private$stamp()
                         return(invisible(self))
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
                             warning(paste('Removed', sum(ix), 'edges from graph that have NA in both n1 and n2, edges should have either n1 or n2 NA'))
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

                           ## remove reserved fields or else we mess up merge below
                           if (!is.null(edges$from))
                             edges$from = NULL

                           if (!is.null(edges$to))
                             edges$to = NULL

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
                         private$stamp()

                         ## label edges with class
                         self$edges$mark(class = self$edges$class)
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
                         lleft = lright = c()

                         if (length(NodeObj$eleft)>0)
                           {
                             lleft = NodeObj$eleft$dt[!(edge.id %in% eid),
                                                      intersect(c(n1, n2), nid)]
                           }

                         if (length(NodeObj$eright)>0)
                         {
                             lright = NodeObj$eright$dt[!(edge.id %in% eid),
                                                        intersect(c(n1, n2), nid)]
                           }
                         
                         edges = convertEdges(nodes, edges, metacols = TRUE)
                         nodes = NodeObj$gr
                         ## assign any new loose ends caused by breaking
                         ## these nodes from their previous partners
                         nodes$loose.left = nodes$loose.left | (nodes$node.id %in% lleft)
                         nodes$loose.right = nodes$loose.right | (nodes$node.id %in% lright)
                         private$gGraphFromNodes(nodes = nodes, edges = edges)

                         private$stamp()
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

                       timestamp = function() private$ptimestamp,
                       
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
                       ## the edge weights of this igraph will be populated by the width of the target interval
                       igraph = function()
                       {
                         v = as.data.frame(
                           values(self$gr)[, c('index', 'node.id', 'snode.id')])

                         if (nrow(self$sedgesdt)==0)
                         {
                           G = igraph::graph_from_adjacency_matrix(sparseMatrix(1,1, x = 0, dims = rep(nrow(v), 2)))
                           V(G)$index = v$index
                           V(G)$node.id = v$node.id
                           V(G)$snode.id = v$snode.id
                           return(G)
                         }
                         
                         ed = as.data.frame(self$sedgesdt)[, c("from", "to", "sedge.id", "edge.id", "type")]
                         G = igraph::graph_from_data_frame(ed,
                                                           vertices = v,
                                                           directed = TRUE)
                         
                         
                         edG = self$sedgesdt[.(igraph::E(G)$sedge.id), ]
                         E(G)$weight = width(self$gr)[edG$to]
                         
                         return(G)
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
                         return(copy(convertEdges(self$gr, private$pedges, metacols = TRUE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
                       },

                       
                       ## Returns a gTrack
                       gt = function()
                       {
                         stack.gap = private$pmeta$stack.gap
                         if (is.null(stack.gap)){
                           stack.gap = 5e4
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


#' @name ==.gGraph
#' @title ==.gGraph
#' @description
#'
#' Returns TRUE if two graphs are equivalent by value based on having equivalent nodes and edges
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables are identical by reference i.e. point to the exact same location in memory)
#'
#' @param x gGraph
#' @param y gGraph
#' @return TRUE if objects are equivalent
#' @export
'==.gGraph' = function(x,y){
  return(x$nodes == y$nodes & x$edges == y$edges)
}

#' @name ==.gNode
#' @title ==.gNode
#' @description
#'
#' Returns TRUE if two gNode objects are equivalent, based on having identical node grs
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables point to the exact same location in memory)
#'
#' @param x gNode
#' @param y gNode
#' @return TRUE if objects are equivalent
#' @export
'==.gNode' = function(x,y){
  return(identical(x$gr, y$gr))
}


#' @name ==.gEdge
#' @title ==.gEdge
#' @description
#'
#' Returns TRUE if two gEdge objects are  equivalent, based on having identical edge tables
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables point to the exact same location in memory)
#'
#' @param x gNode
#' @param y gNode
#' @return TRUE if objects are equivalent
#' @export
'==.gEdge' = function(x,y){
  return(identical(x$dt, y$dt))
}



#' @name !=.gGraph
#' @title !=.gGraph
#' @description
#'
#' Returns TRUE if two graphs are not equivalent based on having non equivalent nodes and edges
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables point to the exact same location in memory)
#'
#' @param x gGraph
#' @param y gGraph
#' @return TRUE if objects are not equivalent
#' @export
'!=.gGraph' = function(x,y){
  return(x$nodes != y$nodes & x$edges != y$edges)
}

#' @name !=.gNode
#' @title !=.gNode
#' @description
#'
#' Returns TRUE if two gNode objects are not equivalent, based on having non identical node grs
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables point to the exact same location in memory)
#'
#' @param x gNode
#' @param y gNode
#' @return TRUE if objects are not equivalent
#' @export
'!=.gNode' = function(x,y){
  return(!identical(x$gr, y$gr))
}


#' @name !=.gEdge
#' @title !=.gEdge
#' @description
#'
#' Returns TRUE if two gEdge objects are non equivalent, based on having non identical edge tables
#' (note this is different - intentionally - from identical, which will determine whether
#' the two variables point to the exact same location in memory)
#'
#' @param x gNode
#' @param y gNode
#' @return TRUE if objects are not equivalent
#' @export
'!=.gEdge' = function(x,y){
  return(!identical(x$dt, y$dt))
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
'[.gGraph' = function(obj, i = NULL, j = NULL, with = TRUE, ...){
  if (deparse(substitute(j)) != "NULL")
  {
    edges = obj$edges[j, with = with]
    if (deparse(substitute(i)) == "NULL"){
        nodes = edges$nodes
    }
    else
    {
      nodes = obj$nodes[i, with = with]
      nids = c(nodes$dt$index, nodes$flip$dt$index);
      eid = edges$sdt[to %in% nids & from %in% nids, sedge.id]    
      edges = edges[eid, with = with]
    }

  } else
  {
    if (deparse(substitute(i)) == "NULL"){
        nodes = obj$nodes
        }
    else{
      nodes = obj$nodes[i, with = with]
        }
    nids = c(nodes$dt$index, nodes$flip$dt$index);
    eid = nodes$edges$sdt[to %in% nids & from %in% nids, sedge.id]    
    edges = obj$edges[eid, with = with]
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
                              edgeObj = x$edges,
                              meta = x$meta
                              ))
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
                                            disjoin = FALSE, ## flag whether to collapse ie disjoin the intervals
                                            drop = FALSE ## flag whether in the case of snode.id should we just drop any walks with non-existent edges (as opposed to just erroring out)                                            
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
                              NONO.FIELDS = c('node.id', 'snode.id', 'loose.left', 'loose.right', 'index')
                              nodes = gr.stripstrand(grlu[, setdiff(names(values(grlu)), NONO.FIELDS)])
                              nodes$parent.graph = grlu$grl.ix
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
                            
                            ## mark all class != 'REF' edges in the walks as ALT
                            graph$edges$mark(type = ifelse(graph$edges$dt$class == 'REF', 'REF', 'ALT'))
                                             
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
                          if (is.null(walk.names))
                            walk.names = seq_along(snode.id)
                          names(snode.id) = seq_along(snode.id)
                          private$gWalkFromNodes(snode.id = snode.id,
                                         graph = graph,
                                         circular = circular,
                                         drop = drop, 
                                         meta
                                         )
                        } else if (!is.null(sedge.id)) {
                          walk.names = names(sedge.id)
                          if (is.null(walk.names))
                            walk.names = seq_along(sedge.id)
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
                          dgraph$edges$mark(in.dgraph = TRUE)
                          self$graph$edges$mark(pid = as.character(self$graph$edges$dt$sedge.id))
                          cloned = self$copy
                          self$disjoin(graph = dgraph)

                          if (any(is.na(self$graph$edgesdt$in.dgraph)))
                          {
                            stop('Some grl edges are not present in provided graph, please check inputs')
                          }
                          self$graph$edges$mark(in.dgraph = NULL)
                        }

                        ## fill in remaining metadata
                        if (self$length>0)
                        {
                          tmp = private$pnode[, .(walk.id, snode.id)]                   
                          tmp$wid = width(self$graph$nodes$gr)[abs(tmp$snode.id)]
                          if (!is.null(walk.names))
                            private$pmeta$name = walk.names[private$pmeta$walk.id]
                          private$pmeta$wid = tmp[, .(wid = sum(as.numeric(wid))), keyby = walk.id][.(private$pmeta$walk.id), wid]

                          ## we need to fix factors if drop == TRUE and there are missing walk.ids
                          if (drop == TRUE & any(private$pmeta$walk.id>self$length)){
                            lev = private$pmeta$walk.id
                            private$pmeta[, walk.id := as.integer(factor(as.character(walk.id), lev))]
                            private$pnode[, walk.id := as.integer(factor(as.character(walk.id), lev))]
                            setkey(private$pnode, walk.id)
                            setkey(private$pmeta, walk.id)
                            if (nrow(private$pedge)>0)
                            {
                              private$pedge[, walk.id := as.integer(factor(as.character(walk.id), lev))]
                              setkey(private$pedge, walk.id)
                            }                              
                          }
                        }

                        private$ptimestamp = private$pgraph$timestamp
                        return(self)
                      },

                      #' sets metadata of gWalk object
                      #' (accessible through $dt accessor)
                      set = function(...)
                      {
                        self$check
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
                        self$check
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

                        return(invisible(self))
                      },

                      dts = function(ix = 1:self$length, makelists = TRUE)
                      {
                        self$check
                        out = private$pmeta[.(ix), ]
                        if (makelists == FALSE) ## just return raw pmeta subset 
                          return(out)
                        node.sum = private$pnode[.(ix), .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),][, -1]
                        out = cbind(out, node.sum)
                        if (nrow(private$pedge)>0)
                          out = cbind(out, private$pedge[.(ix), .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),][, -1])
                        return(out)
                      },

                      #' @name rep
                      #' @description
                      #'
                      #' Creates "bubbles" in the underlying graph by replicating the gwalks.  Node replication
                      #' replicates edges going in and out of all replicated nodes.  If an edge connects a pair of replicated nodes
                      #' that edge will be replicated across all pairs of those replciated nodes.   Walk replication will create "longer bubbles"
                      #' with fewer edges getting replicated i.e. it will only replicate intra-walk edges within each walk replicate (but not between
                      #' separate walk replicates).
                      #'
                      #' (note that this returns a copy of the gGraph that this node is linked to, i.e. it does not make 
                      #' this gWalk stale)
                      #'
                      #' New graph keeps track of the parent node and edge ids in the original graph using node metadata parent.node.id
                      #' and edge metadata parent.edge.id i.e. the replicated nodes will be connected to the sources of the original nodes
                      #' and if replicated nodes connect to each other, then there will exist an edge connecting
                      #' all of their instances to each other.
                      #'
                      #' @param nodes = gNode object must point to a node in the graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object
                      #' @param times  scalar or vector of length self$length specifying how many times to replicate each of the nodes.
                      #' @author Marcin Imielinski
                      #' @return Returns a pointer to the new nodes
                      rep = function(times)
                      {
                        self$check
                        if (length(self)>0)
                          {
                            self$graph$clone()$rep(self, times)
                          }
                      },
                      
                      print = function()
                      {
                        self$check
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
                        self$check
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
                            out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(node), ", keyby = walk.id]")))
                            out = unique(out, by = "walk.id")
                            out[.(1:self$length), ][[2]]
                          }, error = function(e) NULL)
                        }

                        if (is.null(out) & !missing('edge'))
                        {
                          out = tryCatch({
                            tmpdt = merge(private$pedge, private$pgraph$sedgesdt, by = 'sedge.id')
                            out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(edge), ", keyby = walk.id]")))
                            out = unique(out, by = "walk.id")
                            out[.(1:self$length)][[2]]
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
                        self$check
                        if (self$length==0){
                          return(invisible(self))
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

                        private$ptimestamp = private$pgraph$timestamp
                        setkey(private$pnode, walk.id)
                        setkey(private$pedge, walk.id)
                        setkey(private$pmeta, walk.id)

                        return(invisible(self))
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
                        self$check
                        tmpgr = private$pgraph$nodes$gr

                        ## byval = tmpgr %N% private$pgraph$nodes[private$pnode$snode.id]$gr
                        ## if (!is.null(by)) ## add user provided by to byval if exists
                        ##   byval = paste(byval, values(tmpgr)[, by])

                        ## we will first simplify tmpg
                        ## then rematch the walk nodes to tmpg
                        newgraph = private$pgraph$copy
                        newgraph$nodes$mark(parent.node.id = 1:length(newgraph))
                        newgraph$simplify(FUN = NULL)
                        
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
                          match(abs(newnode$snode.id), newgraph$nodes$dt$parent.node.id)

                        newnode = newnode[!is.na(snode.id), ]
                        
                        private$gWalkFromNodes(
                                  snode.id = split(newnode[, snode.id], newnode[, walk.id]),
                                  graph = newgraph,
                                  meta = private$pmeta)

                        setkey(private$pnode, walk.id)
                        setkey(private$pedge, walk.id)
                        setkey(private$pmeta, walk.id)

                        private$ptimestamp = private$pgraph$timestamp
                        return(invisible(self))
                      },

                      mark = function(...)
                      {
                        self$check
                        self$nodes$mark(...)
                        self$edges$mark(...)
                      },
                      
                      gtrack = function(name = NULL, stack.gap = 1e5, ...)
                      {
                        self$check

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
                      }
                    ),

                    private = list(
                      ## data.table of node.ids in the walk
                      pnode = NULL,
                      
                      ## data.table of edge.ids in the walk
                      pedge = NULL,

                      ## data.table of walk metadata 
                      pmeta = NULL,

                      ## timestamp of underlying graph, used to check staleness
                      ptimestamp = NULL,
                      
                      ## Pointer to the graph the snode.id/sedge.ids come from
                      ## If initialized with grl, this is the graph from the gWalk
                      pgraph = NULL,
                     
                      gWalkFromNodes = function(snode.id,
                                              graph,
                                              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
                                              circular = NULL, ## specifies whether the given walk is circular
                                              drop = FALSE ## logical flag specifying whether to drop provided walks traversing non-existent edges
                                              )  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first
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

                        ## always faster to unlist than lapply
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
                        setkey(pnode, walk.id)

                        ## just in case there is more than one edge
                        ## connecting a node pair (possible)
                        sedgesdt = unique(sedgesdt, by = c("from", "to"))

                        if (nrow(pedge)>0){
                          pedge$sedge.id = sedgesdt[.(pedge$from, pedge$to), sedge.id]
                          setkey(pedge, walk.id)
                          if (!drop) ## default behavior is we stop if the provided node lists refer to non-existent edges
                            {
                              if (any(is.na(pedge$sedge.id))){
                                stop('One or more provided walks refers to a non-existent edge, please check your input or re-instantiate with drop = TRUE to ignore this error and remove those walks')
                              }
                            }
                          else
                          {
                            good.walks = setdiff(pnode$walk.id, pedge[is.na(sedge.id), walk.id])
                            pnode = pnode[.(good.walks), ]
                            pedge = pedge[.(good.walks), ]
                            private$pmeta = private$pmeta[good.walks, ]
                          }

                          if (nrow(pedge)>0)
                          {
                            pedge[, walk.iid := 1:.N, by = walk.id]
                          }
                        }

                        if (nrow(pnode)>0)
                          {
                            pnode[, walk.iid := 1:.N, by = walk.id]
                          }

                        private$pnode = pnode
                        private$pedge = pedge
                        private$pgraph = graph
                        private$ptimestamp = private$pgraph$timestamp
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
                        private$ptimestamp = private$pgraph$timestamp
                      }
                    ),
                    active = list(
                      ## Returns a GRangesList of walks in the graph
                      grl = function()
                      {
                        self$check
                        if (self$length == 0){
                            return(GRangesList(private$pgraph$gr[c()])[c()])
                        }
                        ## Does not get both strands only 1 strand
                        nix = private$pgraph$queryLookup(private$pnode$snode.id)$index
                        out = split(private$pgraph$gr[nix], private$pnode$walk.id)[as.character(1:self$length)]
                        values(out) = private$pmeta
                        return(out)
                      },

                      meta = function()
                      {
                        self$check
                        return(private$pmeta)
                      },

                      ## object is stale if the recorded timestamp of the gGraph
                      ## != timestamp the actual gGraph pointed to by pgraph 
                      ## suggesting that the indices are no longer valid
                      stale = function() if (!is.null(private$ptimestamp)) private$ptimestamp != private$pgraph$timestamp else FALSE,
                                            
                      ## checks if object is stale i.e.
                      check = function() if (self$stale) stop('object is stale, underlying gGraph has changed. You will need to re-instantiate.'),

                      nodes = function()
                      {
                        self$check
                        return(private$pgraph$nodes[private$pnode$snode.id])
                      },

                      nodesdt = function()
                      {
                        private$pnode
                      },

                      edgesdt = function()
                      {
                        private$pedge
                      },

                      snode.id = function()
                      {
                        self$check
                        if (self$length==0)
                          return(list())
                        ix = 1:self$length
                        tmp = private$pnode[.(ix), .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),snode.id]
                        return(tmp)
                      },

                      sedge.id = function()
                      {
                        self$check
                        if (self$length==0)
                          return(list())
                        ix = 1:self$length
                        tmp = private$pedge[.(ix), .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),sedge.id]
                        return(tmp)
                      },

                      copy = function() self$clone(),

                      ## returns a length(self) logical vector specifying whether
                      ## each walk is circular or not
                      circular = function()
                      {
                        self$check
                        return(private$pmeta$circular)
                      },

                      length = function() {
                        self$check
                        return(nrow(private$pmeta))
                      },

                      edges = function()
                      {
                        self$check
                        if (self$length == 0){
                          return(gEdge$new(graph = private$pgraph))
                        }
                        seid = private$pedge$sedge.id
                        seid = seid[!is.na(seid)]
                        return(private$pgraph$edges[seid])
                      },

                      footprint = function()
                      {
                        self$check
                        return(sort(reduce(gr.stripstrand(self$nodes$gr))))
                      },

                      gt = function()
                      {
                        self$check
                        return(self$gtrack())
                      },
                      
                      graph = function()
                      {
                        self$check
                        return(private$pgraph)
                      },
                      
                      dt = function() {
                        self$check
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
`c.gWalk` = function(..., force = FALSE)
{                            
  gWalk.list=list(...)
  isg = sapply(gWalk.list, function(x) class(x)[1]=='gWalk')

  if(any(!isg)){
    stop('Error: All inputs must be of class gWalk.')
  }
  
  ##Check to make sure they all come from the same graph
  if (!force)
  {
    graphs = lapply(gWalk.list, function(x) x$graph)
    if (length(graphs)>1)
    {
      gg = graphs[[1]]
      for (i in 2:length(graphs))
        if (graphs[[i]] != gg)
          stop('All inputs must point to the same graph')
    }
  }     

  ##Get all the pnode.id's to create new gNode
  snode.id = do.call('c', lapply(gWalk.list, function(x) x$snode.id))

  metas = rbindlist(lapply(gWalk.list, function(x) x$meta), fill = TRUE)

  return(gWalk$new(snode.id = snode.id, graph = gWalk.list[[1]]$graph, meta = metas))
}


#' @name gW
#' @title create gWalk
#' @description
#'
#' Wrapper that instantiates a gWalk object from a variety of different inputs.  If a gGraph is provided
#' as input (in conjunction with sedge.id, snode.id, or grl  input) then will check for existence of
#' the provided walks in the provided graph, and will error out if those walks do not exist).
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
#' @param drop logical flag only relevant when snode.id is specified indicating whether to drop any walks that refer to non-existent edges 
#' @return gWalk object of provided walks, with pointers back to the graph to which they refer.
gW = function(snode.id = NULL,
              sedge.id = NULL,
              grl = NULL,
              graph = NULL,
              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
              circular = NULL,  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first
              drop = FALSE, 
              disjoin = FALSE ## flag whether to collapse ie disjoin the intervals
              )
{
  return(gWalk$new(snode.id = snode.id,
                    sedge.id = sedge.id,
                    grl = grl,
                    graph = graph,
                    meta = meta,
                   circular = circular,
                   drop = drop,
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
'[.gWalk' = function(obj, i = NULL, with = TRUE, ...){
  walks = obj$clone()
  if (with)
    {
      inew = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()),obj$dts(makelists = FALSE), parent.frame(2)), error = function(e) NULL)
      if (is.null(inew))
        inew = i ## just give up      
    }
  else
  {
    inew = i
  }

  walks$subset(inew)
  return(walks)
}



#' @name %&%
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
#' @exportMethod lengths
setMethod("lengths", c("gWalk"),
          function(x) {
            return(lengths(x$snode.id))
          })

#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' @param x a gGraph object
#'
#' @return the seqinfo of this graph
#' @exportMethod seqinfo
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
#' @exportMethod seqinfo
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
#' @exportMethod seqinfo
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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
#' @export
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





