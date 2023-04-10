#' @name gGnome
#' @title gGnome
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
#' @importFrom parallel mclapply
#' @importFrom reshape2 melt
#' @importFrom VariantAnnotation readVcf info
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
#' @import fishHook
#' @useDynLib gGnome
"_PACKAGE"

#' @name cbind
#' @title cbind wrapper
#'
#' @description
#' Forcing correct call of cbind
#'
#' @param ... arguments to cbind
#'
#' @return vector of combined arguments
cbind = function(..., deparse.level = 1) {
    lst_ = list(...)
    ## anyS4 = any(vapply(lst_, inherits, FALSE, c("DFrame", "DataFrame", "List")))
    anyS4 = any(vapply(lst_, isS4, FALSE))
    if (anyS4) cbind.DataFrame(...) else BiocGenerics::cbind(..., deparse.level = deparse.level)
}


## ================= gNode class definition ================== ##
#' @name gNode
#' @title gNode
#' @description
#' gNode object.
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
                      #' @details
                      #' ```
                      #' gg$nodes[1:5]$mark(col = "purple")
                      #' gg$nodes$mark(changed = FALSE)
                      #' @param  ... name = value pairs of scalar or vector (length edges in graph) arguments
                      #' ```
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
                      #' is flipped with respect to the reference, the left edges will be attached to the right of the node on the reference.)
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

                        if (!is.null(dim(i)))
                          stop('subscript dimensions malformed')
                        
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

                      #' @name eval
                      #' @description
                      #'
                      #' Evaluates expression over edge metadata for all nodes in self
                      #' with a keyby node.id returning a self$length vector
                      #' @param x expression over edge metadata
                      eval = function(x, right = TRUE, all = FALSE)
                      {
                        self$check
                        if (self$length==0)
                          return(NULL)

                        tmpdt = data.table()
                        if (right)
                          tmpdt = rbind(tmpdt, self$eright$dt[, node.id := n1])

                        if (!right | all)
                          tmpdt = rbind(tmpdt, self$eleft$dt[, node.id := n2])

                        if (nrow(tmpdt)==0)
                          return(rep(NA, self$length))

                        ## remove any dups which will inflate the function
                        ## eg if the same edges happen to be on the left and
                        ## right of the given node
                        tmpdt = unique(tmpdt, by = c("edge.id", "node.id"))

                        out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(x), ", keyby = node.id]")))[.(self$id), V1]

                        return(out)
                      },

                      #' value for "cluster"
                      #' 
                      #'
                      #' @param weak character scalar that can take one of the following possible values - "weak" or "strong" specifying weakly or strongly connected components, walktrap specifying cluster_walktrap community detection
                      #' @author Marcin Imielinski

                      #' @name clusters
                      #' @description
                      #'
                      #' Marks nodes in graph with metadata field $cluster
                      #' based on one of several algorithms, selected by mode.
                      #'
                      #' Unlike the gGraph version of this function, this usage enables computation of clusters
                      #' based on a subset of nodes in the graph (i.e. ignoring certain nodes)
                      #' while still marking the original graph (and leaving the excluded graph with an NA)
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

                        sgraph$clusters(mode = mode)

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
                      },

                      #' @name loose.degree
                      #'
                      #' @param orientation (character) one of 'left' or 'right'
                      #'
                      #' @return number of loose ends with given orientation
                      loose.degree = function(orientation)
                      {
                        if (!(orientation %in% c('right', 'left'))){
                          stop('Bad orientation: "', orientation, '". orientation must be either "right", or "left"')
                        }
                        self$check
                        if (!nrow(private$pgraph$edgesdt))
                          return(0)                        
                        op.orientation = 'right'
                        if (orientation == 'right'){
                            op.orientation = 'left'
                        }

                        deg = rowSums(cbind(private$pgraph$edgesdt[, sum(n1.side == orientation), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == orientation), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)

                        op.deg = deg
                        if (any(private$porientation<0)){
                          op.deg = rowSums(cbind(private$pgraph$edgesdt[ , sum(n1.side == op.orientation), keyby = n1][.(private$pnode.id), V1], private$pgraph$edgesdt[, sum(n2.side == op.orientation), keyby = n2][.(private$pnode.id), V1]), na.rm = TRUE)}

                        return(pmax(0, ifelse(private$porientation>0, deg, op.deg), na.rm = TRUE))
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

                      #' @name copy
                      #' @title copy
                      #' @description
                      #'
                      #' Return a deep copy of the graph
                      #'
                      #' @return copy of the object
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

                      #' @name footprint
                      #' @description
                      #' 
                      #' Returns the reduced genomic footprint of this object
                      #' 
                      #' @return GRanges of footprint of this object
                      footprint = function()
                      {
                        self$check
                        return(sort(reduce(self$gr)))
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

                      #' @name flip
                      #' @title flip
                      #' @description
                      #'
                      #' returns flipped version of this node
                      #'
                      #' @return reverse complemented gNode
                      flip = function()
                      {
                        self$check
                        sid = self$dt$snode.id
                        return(self$graph$nodes[-sid])
                      },

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
                        if (nrow(private$pgraph$sedgesdt)==0)
                          return(gEdge$new(graph = private$pgraph))
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
                        if (nrow(private$pgraph$sedgesdt)==0)
                          return(gEdge$new(graph = private$pgraph))
                        sedge.id = private$pgraph$sedgesdt[from %in% private$pindex, sedge.id]
                        return(gEdge$new(sedge.id, graph = private$pgraph))
                      },

                      
                      #' @name loose
                      #' @description
                      #' Returns a GRanges containing loose ends information
                      #' @details
                      #' Returns info on loose ends connected to the nodes in this gNode. If there are no loose ends, returns an empty GRanges.
                      #'
                      #' @return GRanges Loose ends connected to the nodes in this gNode Object
                      loose = function()
                      {
                        return(c(self$lleft, self$lright))
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
                          ##FIXED: R6 won't let you do this with nested active bindings, i.e.
                          ## gg$nodes$loose.left = TRUE won't work; but nodes = gg$nodes; nodse$loose.left = TRUE will

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

                          value = value | rdeg==0

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
                          ##FIXED: R6 won't let you do this with nested active bindings, i.e.
                          ## gg$nodes$loose.left = TRUE won't work; but nodes = gg$nodes; nodse$loose.left = TRUE will

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

                      #' @name terminal
                      #' @title terminal
                      #' @description
                      #' Get a GRanges containing all terminal loose ends
                      #' @return GRanges
                      terminal = function()
                      {
                        self$check
                        if (self$length==0)
                          {
                            return(GRanges(seqlengths = seqlengths(self)))
                          }
                        ix = which(self$ldegree==0 | self$rdegree==0)
                        terminal.left = self[self$ldegree==0]$lleft
                        terminal.right = self[self$rdegree==0]$lright
                        return(c(terminal.left, terminal.right))
                      },

                      ldegree = function()
                      {
                        return(self$loose.degree(orientation = 'left'))
                      },

                      rdegree = function()
                      {
                        return(self$loose.degree(orientation = 'right'))
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
                        return(gNode.loose(self, orientation = 'left'))
                      },
                      
                      
                      #' @name lright
                      #' @description
                      #' 
                      #' Returns a GRanges containing the loose ends connected to the right of the nodes in this
                      #' gNode. If there are no loose ends to the right, returns an empty GRanges
                      #'
                      #' @return GRanges Loose ends connected to the right of the nodes in this gNode Object
                      lright = function()
                      {
                        return(gNode.loose(self, orientation = 'right'))
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
                        return(gGraph$new(nodeObj = self,
                                          edgeObj = edgeObj,
                                          meta = private$pgraph$meta))
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

#' @name intersect.gNode
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
"intersect.gNode" = function(x, y)
          {             
              if(!identical(x$graph, y$graph)) {
                  stop("Arguments do not point to the same graph")
              }

              new.ids = intersect(x$id, y$id)
              return(gNode$new(new.ids, x$graph))
          }
registerS3method("intersect", "gNode", intersect.gNode, envir = globalenv())

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
#' @name gEdge
#' @title gEdge
#'
#' @description
#' gEdge obejct
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
                      #' @param  ... name = value pairs of scalar or vector (length edges in graph) arguments
                      #'
                      #' @details
                      #' ```
                      #' gg$edges[1:5]$mark(col = "purple")
                      #' gg$edges$mark(changed = FALSE)
                      #' ```
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

                        if (!is.null(dim(i)))
                          stop('subscript dimensions malformed')

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

                      shadow = function() self$junctions$shadow,

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

                      #' @name footprint
                      #' @description
                      #' 
                      #' Returns the reduced genomic footprint of this object
                      #' 
                      #' @return GRanges of footprint of this object
                      footprint = function()
                      {
                        self$check
                        return(sort(reduce(gr.stripstrand(unlist(self$grl)))))
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
                        out = copy(private$pedges)
                        if (nrow(out))
                          {
                            out$n1 = self$graph$dt$snode.id[out$from]
                            out$n2 = self$graph$dt$snode.id[out$to]
                          }
                        return(out)
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

                        if (length(private$psedge.id)==0)
                        {
                          pedges = copy(private$pedges)
                        }
                        else
                        {
                          pedges = private$pgraph$sedgesdt[.(private$psedge.id), ]
                        }
                        return(copy(convertEdges(private$pgraph$gr, pedges, metacols = TRUE, cleanup = FALSE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
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
                        grl = GenomicRanges::split(c(gr1, gr2), rep(1:length(gr1), 2))[as.character(1:length(gr1))]
                        names(grl) = private$edges$sedge.id
                        pedges = private$pgraph$sedgesdt[.(private$psedge.id), ]
                        meta = cbind(pedges, data.table(bp1 = gr.string(gr1), bp2 = gr.string(gr2)))
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


  if (any(deparse(substitute(i)) == "NULL"))
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


#' @name intersect.gEdge
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
"intersect.gEdge" = function(x, y)
{
  if(!identical(x$graph, y$graph)) {
    stop("Arguments do not point to the same graph")
  }
  
  new.ids = intersect(x$id, y$id)
  return(gEdge$new(new.ids, x$graph))
}
registerS3method("intersect", "gEdge", intersect.gEdge, envir = globalenv())



## ================== Junction class definition ================== ##
#' @name Junction
#' @title Junction
#'
#' @description
#' Junction object
#'
#' @details
#' signed adjacency between two genomic loci
#' 
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

                         #' @name gw
                         #' @description
                         #'
                         #' Returns a gWalk representing the left -> right traversal of this junction +- padding
                         #'
                         #' @return GRangesList of the junctions in this Junction Object
                         gw = function(pad = 0)
                         {
                           left = gr.flipstrand(flank(self$left, -pad+1))
                           right = flank(self$right, -pad+1)
                           grl = grl.pivot(GRangesList(left, right))
                           values(grl) = self$meta
                           return(gW(grl = grl))
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

                           if (!is.null(dim(i)))
                             stop('subscript dimensions malformed')

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

                         #' @name set
                         #' @title set
                         #' @description
                         #' set metadata for junction object
                         #'
                         #' @details
                         #' sets metadata to either a scalar or vector
                         #' where the vector is the same length as the junction object
                         #' ```
                         #' jj$set(cn = 7)
                         #' jj$set(cn = c(1,2,3))
                         #' jj$set(col = c("red", "blue", "green"))
                         #' ```
                         #' (after setting, metadata is accessible through $dt accessor)
                         set = function(...)
                         {
                           self$check
                           args = list(...)
                           
                           NONO.FIELDS = c('junc','ix')
                           args = list(...)

                           meta = as.data.table(values(private$pjuncs))

                           if (any(names(args) %in% NONO.FIELDS)){
                             stop(paste('Cannot modify the following reserved gWalk metadata fields:', paste(NONO.FIELDS, collapse = ', ')))                  
                           }
                           for (arg in names(args)){
                             meta[[arg]] = args[[arg]]
                           }
                           values(private$pjuncs) = as.data.frame(meta)
                           return(invisible(self))
                         },

                         #' @name print
                         #' @description
                         #' 
                         #' Prints out the Junction Object. Prints the length and the GRangesList of the junctions.
                         print = function() ## junction::print
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


                         #' @name footprint
                         #' @description
                         #' 
                         #' Returns the reduced genomic footprint of this object
                         #' 
                         #' @return GRanges of footprint of this object
                         footprint = function()
                         {
                           self$check
                           return(sort(reduce(gr.stripstrand(unlist(self$grl)))))
                         },

                         #' @name shadow
                         #' @description
                         #' 
                         #' Returns the "shadow" of this junction object
                         #' which for intrachromosomal junctions is the interval between the breakpoints
                         #' and for interchromosomal junctions is just the locations of the end points
                         #' Note that this is different from footprint: which is only the locations
                         #' of the breakpoints, and the output is not reduced but each junction will produce
                         #' either one or two ranges.  The id can be traced by $id metadata.
                         #' @return GRanges of footprint of this object
                         shadow = function()
                         {
                           self$check
                           if (!length(self))
                             return(unlist(self$grl))
                           gru = gr2dt(grl.unlist(self$grl))
                           shadow = dt2gr(gru[, .(start = min(start), end = max(start)), by = .(seqnames, grl.ix)], seqlengths = seqlengths(self))[, 'grl.ix']
                           names(values(shadow)) = 'id'
                           return(shadow)
                         },


                         #' @name junc
                         #' @description
                         #' short character blurb of this junction
                         #'
                         junc = function()
                         {
                           bp = grl.unlist(self$grl)
                           bpstr = gr.string(bp)
                           data.table(str = bpstr, ix = bp$grl.ix)[, .(junc = paste(str, collapse = ' <-> ')), keyby = ix]$junc
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
                           if (!length(self))
                             return(c())
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
                           tmp = grl.pivot(private$pjuncs)
                           tmp[[1]] = gr.flipstrand(tmp[[1]])
                           tmp[[2]] = gr.flipstrand(tmp[[2]])
                           newjunc = grl.pivot(tmp)
                           values(newjunc) = values(private$pjuncs)
                           private$pjuncs = newjunc
                           return(self)
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
                           if (!self$length)
                             return(c())
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
#' @exportMethod setdiff
#' @return new Junction Object containing the difference between x and y
#' @export
setMethod("setdiff", c('Junction', "Junction"), function(x, y, pad = 0, ...)
{  
  ov = ra.overlaps(x$grl, y$grl, pad = pad)
  ix = setdiff(1:length(x$grl), ov[,1])
  return(x[ix])
})



#' @name intersect.Junction
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
#setMethod("intersect", c('Junction', 'Junction'), function(x, y, pad = 0, ...) {
"intersect.Junction" = function(x, y, pad = 0, ...) {

  ov = ra.overlaps(x$grl, y$grl, pad = pad)
  return(unique(x[ov[, 'ra1.ix']], pad = pad))
}
registerS3method("intersect", "Junction", intersect.Junction, envir = globalenv())


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
#' @exportMethod union
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
#' @exportMethod setdiff
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
#' @exportMethod refresh
#' @export
setGeneric("refresh", function(x) standardGeneric("refresh"))
setMethod("refresh", "Junction",
          function(x) {
            return(Junction$new(x$juncs))
          })


## ================== gGraph class definition ================== ##
#' @name gGraph
#' @title gGraph
#'
#' @description
#' a genome graph object
#' 
#' @export
gGraph = setClass("gGraph")
gGraph = R6::R6Class("gGraph",
                     public = list(
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
                       #' @param rck  RCK output directory path  
                       #' @author Joe DeRose                                              
                       initialize = function(genome = NULL,
                                             breaks = NULL,
                                             juncs = NULL,
                                             alignments = NULL,
                                             prego = NULL,
                                             jabba = NULL,
                                             cougar = NULL,
                                             weaver = NULL,
                                             remixt = NULL,
                                             rck = NULL,
                                             walks = NULL,
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
                           else if (!is.null(alignments))
                           {
                             ne = alignments2gg(alignments, verbose)
                           }
                           else if (!is.null(prego))
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
                             meta$purity = ne$purity
                             meta$ploidy = ne$ploidy
                           }
                           else if (!is.null(cougar))
                           {
                               ne = cougar2gg(cougar)
                               cgg = gG(breaks = ne$breaks, juncs = ne$juncs)
                               ne$nodes = cgg$gr
                               ne$edges = cgg$edgesdt
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
                           else if (!is.null(rck))
                           {
                               ne = rck2gg(rck, haploid = TRUE) ## for now only unphased graphs
                               meta[['name']] = 'RCK'
                               meta[['y.field']] = 'cn'
                               meta[['y0']] = 0
                               meta[['by']] = 'cn'
                           }
                           else if (!is.null(walks))
                           {
                               gn = haplograph(walks, breaks)
                               ne = list(nodes = gn$nodes$gr, edges = gn$edges$dt)
                           }
                           else if(!is.null(breaks) || !is.null(juncs))
                           {                                                           
                               ne = breakgraph(breaks, juncs, genome)
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
                       walks = function(field = NULL, greedy = FALSE, verbose = FALSE)
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
                           return(invisible(self))
                       },

                       #' @name queryLookup
                       #' @description
                       #'
                       #' Returned a data.table of the provided snode.ids, their indicies and the indicies of
                       #' their reverse complements in the graph. data.table is keyed on snode.id.
                       #'
                       #' @param id snode.ids to look up
                       #' @param id signed node ids in graph
                       #' @return data.table of snode.ids, indicies and reverse complement indicies
                       #' @author Joe DeRose
                       queryLookup = function(id) {
                           dt = private$lookup[list(id)]
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
                               ## now we want to map, then project the newedges onto the the udnodes
                               ## applying the aggregation function 
                               dnodes$udnode.id = match(dnodes$query.id, udnodes$query.id)
                           }
                           else
                           {

                               if (!is.character(by))
                                   stop('by must be a character vector specifying one or more node metadata columns')

                               tmp = gr2dt(dnodes)[, dup := duplicated(query.id), by = eval(by)]
                               udnodes = dt2gr(tmp[dup == FALSE, ][, -ncol(tmp), with = FALSE], seqlengths = seqlengths(dnodes))
                               dnodes$udnode.id =
                                   match(do.call('paste', c(as.list(values(dnodes)[, c('query.id', by)]),
                                                            list(sep = ','))),
                                         do.call('paste', c(as.list(values(udnodes)[, c('query.id', by)]),
                                                            list(sep = ','))))

                           }


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

                           ## add annotations from gr if those exist
                           if (length(gr) && ncol(values(gr)))
                           {
                               tmp.nodes = cbind(tmp.nodes, as.data.table(values(gr)[gr.match(dt2gr(tmp.nodes), gr, by = by), ]))
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

                           if (length(self$edges)==0)
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
                           igg$clusters(mode = 'weak')
                           nodemap = data.table(node.id = igg$nodes$gr$node.id.og,
                                                new.node.id = as.integer(factor(igg$nodes$gr$cluster, unique(igg$nodes$gr$cluster))), key = 'node.id',
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
                       reduce = function(...)
                       {
                           self$disjoin(...)
                           self$simplify(...)
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
                           out = GRanges(seqlengths = seqlengths(self))
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
                                 ## out = streduce(c(win[, c()], out[, c()]))
                                 out = streduce(gUtils::grbind(win[, c()], out[, c()]))
                           }

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
                       #' If i and j are specified, graph is first subsetted then 
                       #' clusters computed, then cluster ids are lifted back to mark
                       #' the original graph.
                       #' 
                       #' @param weak character scalar that can take one of the following possible values - "weak" or "strong" specifying weakly or strongly connected components, walktrap specifying cluster_walktrap community detection
                       #' @param i node filter to apply to graph prior to clustering
                       #' @param j edge filter to apply to graph prior to clustering
                       #' @author Marcin Imielinski
                       clusters = function(i = NULL,
                                           j = NULL,
                                           mode = 'weak'
                                           )
                       {
                         filtered = FALSE

                         self$nodes$mark(cluster = NULL)
                         self$nodes$mark(pcluster = NULL)
                         self$nodes$mark(ncluster = NULL)
                         self$nodes$mark(cluster = as.integer(NA))
                         self$nodes$mark(pcluster = as.integer(NA))
                         self$nodes$mark(ncluster = as.integer(NA))
                         graph = self$copy
                         if (length(graph)==0)
                           return(invisible(self))

                         graph$nodes$mark(og.node.id = 1:length(graph))                       
                         if (all(deparse(substitute(i)) != "NULL")
                             | all(deparse(substitute(j)) != "NULL"))
                         {
                           ## NSE voodoo to fill in the parent.frame values before passing on
                           i = deparse(eval(parse(text = substitute(deparse(substitute(i)))), parent.frame()))
                           j = deparse(eval(parse(text = substitute(deparse(substitute(j)))), parent.frame()))
                           if (length(i)>1)
                             i = paste(i, collapse = ' ')

                           if (length(j)>1)
                             j = paste(j, collapse = ' ')

                           graph = eval(parse(text = sprintf('graph[%s, %s]', i, j)))
                                                             
                           graph$clusters(mode = mode)
                           self$nodes[graph$nodes$dt$og.node.id]$mark(
                                                                   cluster = graph$nodes$dt$cluster,
                                                                   pcluster = graph$nodes$dt$pcluster,
                                                                   ncluster = graph$nodes$dt$ncluster
                                                                 )
                           return(invisible(self))
                         }
                            
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
                         names(membership) = as.character(as.integer(self$gr$snode.id))

                         ## "positive membership" membership of the positive strand of the node
                         pmembership = membership[as.character(self$nodes$dt$node.id)]

                         ## "reverse" i.e. negative strand membership (will be different unless there is
                         ## a palindromic community / component
                         nmembership = membership[as.character(-self$nodes$dt$node.id)]

                         ## rename clusters so most popular is first
                         rename = data.table(names(rev(sort(table(c(pmembership, nmembership))))))[, structure(1:.N, names = V1)]

                         ## rename again so that positive clusters are first
                         names(rename) = names(rename)[order(!(names(rename) %in% as.character(pmembership)))]

                         pmembership = rename[as.character(pmembership)]
                         nmembership = rename[as.character(nmembership)]

                         ## final cluster membership of the node is the min of the two clusters
                         ## this will ensure that the node inherits the most popular cluster of its
                         ## two "sides"

                         ## compose strand agnostic cluster name choosing most popular of the signed clusters
                         ## for each node
                         tmp.cluster = pmin(pmembership, nmembership)

                         ## then rename these according to their popularity so they are in decreasing
                         ## order of popularity 
                         tab.cluster = rev(sort(table(tmp.cluster)))
                         rename = structure(1:length(tab.cluster), names = names(tab.cluster))
                         cluster = rename[as.character(tmp.cluster)]

                         self$annotate('cluster', data = cluster, id = self$nodes$dt$node.id,
                                       class = 'node')

                         self$annotate('pcluster', data = pmembership, id = self$nodes$dt$node.id,
                                       class = 'node')
                         
                         self$annotate('ncluster', data = nmembership, id = self$nodes$dt$node.id,
                                       class = 'node')

                         return(invisible(self))
                       },


                       #' @name eclusters
                       #'                        #' @description
                       #' Marks ALT edges belonging (quasi) reciprocal cycles
                       #' @param thresh the distance threshold with which to group nearby quasi-reciprocal junctions - i.e. if thresh=0 then we only consider clusters of exactly reciprocal junctions.
                       #' @param mc.cores parallel
                       #' @param weak logical flag if TRUE will not differentiate between cycles and paths and will return all weakly connected clusters in the graph [FALSE]
                       #' @return numerical vector of the same length, Inf means they r not facing each other
                       #' @author Marcin Imielinski
                       eclusters = function(thresh = 1e3,
                                            range = 1e6,
                                            weak = TRUE,
                                            paths = !weak,
                                            mc.cores = 1,
                                            verbose = FALSE,
                                            chunksize = 1e30,
                                            method = "single")
                       {
                         self$edges$mark(ecluster = as.integer(NA))
                         altedges = self$edges[type == "ALT", ]
                         if (verbose & weak)
                           message('Computing weak eclusters')

                         if (length(altedges)==0){
                           if (verbose){
                             message("No junction in this graph")
                           }
                           return(NULL)
                         }

                         bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix")]
                         bp.dt = gr2dt(bp)

                         ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
                         ixu = unlist(ix)
                         eps = 1e-9
                         ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
                         xt.adj = old.adj = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp), 2))

                         if (verbose){
                           message(sprintf('Computing junction graph across %s ALT edges with distance threshold %s', length(altedges), thresh))
                         }

                         if (!exists(".INF")){
                           .INF = pmax(sum(seqlengths(self)), 1e9)
                         }
                         xt.adj[ixu, ] = do.call(rbind,
                                                 mclapply(ix,
                                                          function(iix)
                                                          {
                                                            if (verbose>1)
                                                              cat('.')
                                                            tmpm =
                                                              gr.dist(bp[iix],
                                                                      gr.flipstrand(bp),
                                                                      ignore.strand = FALSE)+eps
                                                            return(as(tmpm, "Matrix"))
                                                          },
                                                          mc.cores = mc.cores))
                         ## get back to marcin's version
                         adj = xt.adj
                         adj[which(is.na(as.matrix(adj)))] = 0
                         adj[which(as.matrix(adj)>thresh)] = 0

                         ## message(identical(as.logical(adj), as.logical(as.matrix(old.adj))))


                         xt.adj[which(is.na(as.matrix(xt.adj)))] = .INF + 1
                         ## two breakpoints of the same junction should be distance 1
                         bp.pair = t(
                           sapply(unique(bp$grl.ix),
                                  function(ix){
                                    matrix(which(bp$grl.ix==ix), ncol=2, nrow=1)
                                  }))
                         xt.adj[bp.pair] = 1

                         ## do single linkage hierarchical clustering within `range`
                         hcl = stats::hclust(as.dist(xt.adj), method = "single")
                         hcl.lbl = cutree(hcl, h = thresh)
                         bp.dt$hcl = hcl.lbl
                         bp.hcl =
                           bp.dt[,.(hcl.1 = .SD[grl.iix==1, hcl],
                                    hcl.2 = .SD[grl.iix==2, hcl]),
                                 keyby=grl.ix]

                         ## sometimes two breakpoints belong to diff hcl
                         ## merge them!
                         altedges$mark(hcl.1 = bp.hcl[.(seq_along(altedges)), hcl.1])
                         altedges$mark(hcl.2 = bp.hcl[.(seq_along(altedges)), hcl.2])
                         hcl.ig = igraph::graph_from_edgelist(
                                            bp.hcl[, unique(cbind(hcl.1, hcl.2))], directed = FALSE)
                         hcl.comp = components(hcl.ig)
                         altedges$mark(ehcl = as.integer(hcl.comp$membership)[bp.hcl[, hcl.1]])

                         ## connect to MI's code
                         adj[adj>thresh] = 0

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
                         cl = split(1:length(bp), igraph::clusters(graph.adjacency(adj2), ifelse(weak, 'weak', 'strong'))$membership)

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
                           dcl = dunlist(unname(jcl))[, listid := paste0(ifelse(weak, '', 'c'), listid)]
                           if (!weak)
                             altedges[dcl$V1]$mark(ecycle = dcl$listid)
                           altedges[dcl$V1]$mark(ecluster = dcl$listid)
                         }

                         if (verbose)
                           message(sprintf('Annotated %s junction cycles in edge field $ecycle', length(jcl)))

                         if (paths & !weak)
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
                               lapply(all.paths(tmp.adj, sources = sources, sinks = sinks, verbose = verbose)$paths, function(i) x[i])
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
                         return(invisible(self))
                       },

                       #' @name eclusters
                       #' @description
                       #' Marks ALT edges belonging (quasi) reciprocal cycles 
                       #' @param juncs GRangesList of junctions
                       #' @param mc.cores parallel
                       #' @param only_chains TRUE will only pair breakend to its nearest nearest neighbor IFF the nearest neighbor is reciprocal
                       #' @param max.small size below which simple dups and dels are excluded
                       #' @param ignore.strand usually TRUE
                       #' @param weak logical flag if TRUE will not differentiate between cycles and paths and will return all weakly connected clusters in the junction graph [FALSE] 
                       #' @return numerical vector of the same length, Inf means they r not facing each other
                       #' @author Marcin Imielinski
                       eclusters2 = function(thresh = 1e3,
                                            range = 1e6,
                                            weak = TRUE,
                                            paths = !weak,
                                            mc.cores = 1,
                                            verbose = FALSE,
                                            chunksize = 1e30,
                                            method = "single",
                                            return_pairs = FALSE,
                                            ignore.small = TRUE,
                                            max.small = 1e4,
                                            ignore.isolated = TRUE,
                                            strict = c("strict", "one_to_one", "loose"),
                                            min.isolated = max.small,
                                            only_chains = FALSE
                                            )
                       {
                         if (!is.character(strict) ||
                               any(!strict %in% c("strict", "one_to_one", "loose"))) {
                           stop("strict must be one of 'strict', 'one_to_one', or 'loose'")
                         } else if (length(strict) > 1) {
                           strict = "one_to_one"
                         }
                         self$edges$mark(ecluster = as.integer(NA))
                         altedges = self$edges[type == "ALT", ]
                         if (ignore.small &&
                               !is.null(max.small) &&
                               !is.na(max.small) &&
                               max.small > 0) {
                           altedges = altedges[!((class == "DUP-like" |
                                                    class == "DEL-like") &
                                                   altedges$span <= max.small)]
                         }
                         if (verbose & weak)
                           message('Computing weak eclusters')
                         
                         if (length(altedges)==0){
                           if (verbose){
                             message("No junction in this graph")                                     
                           }
                           return(NULL)
                         }

                         deldup = copy3(altedges)[class %in% c("DUP-like", "DEL-like", "INV-like")]

                         if (length(deldup) > 0 && ignore.isolated) {
                           
                           altes = deldup$shadow
                           ## altes$sedge.id = altedges[class %in% c("DUP-like", "DEL-like")]$dt[altes$id]$sedge.id
                           bp = gr.noval(grl.unlist(altedges$grl), keep.col = c("grl.ix", "grl.iix", "class", "sedge.id"))
                           bp$sedge.id.y = bp$sedge.id; bp$sedge.id = NULL
                           addon = deldup$dt[altes$id][, .(sedge.id, class)]
                           altes$sedge.id = addon$sedge.id
                           altes$class = addon$class
                           altes$nbp = altes %N% bp # number of breakpoints of any SV that fall within segment
                           numsum = altedges$shadow %>% gr.sum # using the shadows of all of the SVs not just dels and dups
                           altes = altes %$% numsum
                           iso = ((altes) %Q% (score == 1.0))$id
                           ## rm.edges = unique(altes[iso] %Q% (width < thresh))$sedge.id ## old
                           rm.edges = unique(altes[iso] %Q% (width < min.isolated))$sedge.id
                           rm.dups = S4Vectors::with(altes, sedge.id[class == "DUP-like" & nbp <= 2])
                           rm.dups = c(rm.dups, dedup.cols(gr2dt(altes %*% bp))[sedge.id != sedge.id.y][class == "DUP-like"][, .(all(1:2 %in% grl.iix), class.1 = class.1[1]), by = .(sedge.id, sedge.id.y)][, all(V1 == TRUE) & all(class.1 == "DUP-like"), by = sedge.id][V1 == TRUE]$sedge.id) # removing dups that have only other nested dups
                           rm.inv = S4Vectors::with(altes, sedge.id[class == "INV-like" & nbp <= 2 & width < min.isolated]) # removing inv that are isolated
                           rm.edges = unique(c(rm.edges, rm.dups, rm.inv))
                           keepeid = setdiff(altedges$dt$sedge.id, rm.edges)
                           altedges = altedges[as.character(keepeid)]
                           
                         } # ignoring isolated dup and del edges that are smaller than threshold
                         
                         if (verbose & weak)
                           message("Computing weak eclusters")
                         if (length(altedges) == 0) {
                           if (verbose) {
                             message("No junction in this graph")
                           }
                           return(NULL)
                         }
                         
                         bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix", "class", "edge.id")]
                         bp$m.ix = seq_along(bp)
                         bp.dt = gr2dt(bp)

                         ix = split(1:length(bp), ceiling(runif(length(bp))*ceiling(length(bp)/chunksize)))
                         ixu = unlist(ix)
                         eps = 1e-9
                         ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
                         xt.adj = xt.adj0 = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp), 2))

                         if (verbose){
                           message(sprintf('Computing junction graph across %s ALT edges with distance threshold %s', length(altedges), thresh))
                         }

                         if (!exists(".INF")){
                           .INF = pmax(sum(seqlengths(self)), 1e9)
                         }
                         bp.pair = as.matrix(dcast(data.table(ix = seq_along(bp),
                           grl.ix = bp$grl.ix,
                           grl.iix = bp$grl.iix),
                           grl.ix ~ grl.iix, value.var = "ix")[,2:3])
                         bp.pair = rbind(bp.pair, cbind(bp.pair[,2], bp.pair[,1]))
                         ifun = function(iix, ignore.strand = FALSE,
                                         verbose = FALSE, eps = 1e-9) {
                           if (verbose > 1)
                             cat(".")
                           tmpm = gr.dist(bp[iix], gr.flipstrand(bp), ignore.strand = ignore.strand) +
                             eps
                           return(as(tmpm, "Matrix"))
                         }
                         xt.adj[ixu, ] = do.call(rbind,
                           mclapply(ix, ifun, ignore.strand = FALSE,
                             mc.cores = mc.cores))
                         ## only_chains = TRUE, enforcing that nearest breakpoints are only considered
                         diag(xt.adj) = NA_real_
                         if (only_chains) {
                           
                           ## enforcing that no distances between breakends from the same junction  are considered
                           xt.adj[bp.pair] = NA_real_
                           xt.adj0[ixu, ] = do.call(rbind,
                             mclapply(ix, ifun, ignore.strand = TRUE,
                               mc.cores = mc.cores))
                           ## enforcing that self-to-self breakend distances are not considered
                           diag(xt.adj0) = NA_real_
                           ## enforcing that no distances between breakends from the same junction  are considered
                           xt.adj0[bp.pair] = NA_real_
                           suppressWarnings({
                             nearest_ix = dunlist(lapply(
                               seq_len(nrow(xt.adj)),
                               function(x) which(xt.adj[x,] == min(xt.adj[x,], na.rm = T)))) %>%
                               as.matrix
                             nearest0_ix = dunlist(lapply(
                               seq_len(nrow(xt.adj0)),
                               function(x) which(xt.adj0[x,] == min(xt.adj0[x,], na.rm = T)))) %>%
                               as.matrix
                           })
                           if (nrow(nearest_ix) & nrow(nearest0_ix)) {
                             if (strict == "strict") {
                               nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
                               nearest0_ix = cbind(rowMins(nearest0_ix), rowMaxs(nearest0_ix))
                               nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
                               nearest0_ix = nearest0_ix[duplicated(nearest0_ix),,drop = FALSE]
                               nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
                               nearest0_ix = rbind(nearest0_ix, cbind(nearest0_ix[,2], nearest0_ix[,1]))
                               nearest_ix = as.matrix(merge(nearest_ix, nearest0_ix)) # if there is a breakend closer but in the wrong orientation, that distance will be thrown out downstream, i.e. there is a one to one match of breakend and any nearest breakend that is in the wrong orientation disqualifies the clustser
                             } else if (strict == "one_to_one") {
                               nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
                               nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
                               nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
                               nearest_ix = as.matrix(nearest_ix) # only one-to-one breakends are considered
                             } else if (strict == "loose") {
                               nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
                               nearest_ix = nearest_ix[!duplicated(nearest_ix),,drop = FALSE]
                             }
                           }
                         } else if (!only_chains) {
                           ## two breakpoints of the same junction should be distance 1
                           xt.adj[bp.pair] = 1
                         }
                         rm(xt.adj0)
                         ## get back to marcin's version
                         adj = xt.adj
                         adj[which(is.na(as.matrix(adj)))] = 0
                         adj[which(as.matrix(adj)>thresh)] = 0

                         ## message(identical(as.logical(adj), as.logical(as.matrix(old.adj))))


                         xt.adj[which(is.na(as.matrix(xt.adj)))] = .INF + 1

                         if (only_chains && nrow(nearest_ix)) {
                           tmp = xt.adj[nearest_ix]
                           xt.adj[] = .INF + 1
                           adj[] = 0
                           xt.adj[nearest_ix] = tmp
                           adj[nearest_ix] = tmp
                         }

                         ## below update - grab all the properly oriented links regardless
                         ## of distance threshold - use metadata
                         
                         ## here we are making the reciprocal breakend metadata...
                         
                         ## creating a breakend table of reciprocal distances
                         dt = Matrix::which(xt.adj < thresh, arr.ind = T)
                         dt = unique(data.table(V1 = rowMins(dt), V2 = rowMaxs(dt)))
                         dt[, bp.dist := xt.adj[dt[, cbind(V1, V2)]]]
                         dt$sign = ifelse(strand(bp[dt$V1]) == "+",
                           ifelse(gr.flipstrand(bp[dt$V1]) > bp[dt$V2], -1, 1),
                           ifelse(gr.flipstrand(bp[dt$V1]) < bp[dt$V2], -1, 1))
                         dt2 = dt[,idx := seq_len(.N)] %>%
                           melt(measure.vars = c("V1", "V2"))
                         meta = cbind(gr2dt(bp)[dt2$value], dt2)[order(idx)]
                         meta = rbind(meta,
                           gr2dt(bp)[!paste(grl.ix, grl.iix) %in%
                                       meta[, paste(grl.ix, grl.iix)]],
                           fill = T)
                         meta$thresh = thresh

                         dt = Matrix::which(adj > 0, arr.ind = T)
                         dt = unique(data.table(cbind(rowMins(dt), rowMaxs(dt))))
                         dt[, bp.dist := xt.adj[dt[, cbind(V1, V2)]]]
                         dt$sign = ifelse(strand(bp[dt$V1]) == "+",
                           ifelse(gr.flipstrand(bp[dt$V1]) > bp[dt$V2], -1, 1),
                           ifelse(gr.flipstrand(bp[dt$V1]) < bp[dt$V2], -1, 1))
                         dt2 = dt[,idx := seq_len(.N)] %>% melt(measure.vars = c("V1", "V2"))
                         meta.links = cbind(gr2dt(bp)[dt2$value], dt2)[order(idx)]
                         meta.links = rbind(meta.links,
                           gr2dt(bp)[paste(grl.ix, grl.iix) %nin% meta.links[, paste(grl.ix, grl.iix)]],
                           fill = T)
                         meta.links$thresh = thresh

                         ## do single linkage hierarchical clustering within `range`
                         hcl = stats::hclust(as.dist(xt.adj), method = "single")
                         hcl.lbl = cutree(hcl, h = thresh)
                         bp.dt$hcl = hcl.lbl
                         bp.hcl =
                           bp.dt[,.(hcl.1 = .SD[grl.iix==1, hcl],
                             hcl.2 = .SD[grl.iix==2, hcl]),
                             keyby=grl.ix]

                         ## sometimes two breakpoints belong to diff hcl
                         ## merge them!
                         altedges$mark(hcl.1 = bp.hcl[.(seq_along(altedges)), hcl.1])
                         altedges$mark(hcl.2 = bp.hcl[.(seq_along(altedges)), hcl.2])
                         hcl.ig = igraph::graph_from_edgelist(
                           bp.hcl[, unique(cbind(hcl.1, hcl.2))], directed = FALSE)
                         hcl.comp = components(hcl.ig)
                         altedges$mark(ehcl = as.integer(hcl.comp$membership)[bp.hcl[, hcl.1]])

                         ## connect to MI's code
                         adj[adj>thresh] = 0

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
                         if (length(ix)) 
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
                         cl = split(1:length(bp), igraph::clusters(graph.adjacency(adj2), ifelse(weak, 'weak', 'strong'))$membership)

                         ## choose only clusters with length > 1
                         cl = cl[S4Vectors::elementNROWS(cl)>1]
                         cl = cl[order(S4Vectors::elementNROWS(cl))]

                         jcl.0 = lapply(cl, function(x) unique(sort(bp$edge.id[x])))
                         jcls = sapply(jcl.0, paste, collapse = " ")
                         jcl = jcl.0[!duplicated(jcls)]
                         adj3 = adj2
                         altedges$mark(ecycle = as.character(NA))

                         self$nodes$mark(ecluster = NA_integer_)
                         self$set(recip_event = data.table())
                         
                         if (NROW(jcl)>0)
                         {
                           ## dcl = dunlist(unname(jcl))[, listid := paste0(ifelse(weak, '', 'c'), listid)]
                           dcl.0 = dunlist(unname(jcl.0))
                           setnames(dcl.0, "listid", "jcl.0")
                           dcl = dunlist(unname(jcl))[, `:=`(listid, paste0(ifelse(weak,
                             "", "c"), listid))]
                           dcl.0$listid = dcl[V1 %K% dcl.0$V1]$listid
                           cl.1 = dunlist(unname(cl))
                           setnames(cl.1, "listid", "jcl.0")
                           cl.eclust = merge(cl.1, 
                             unique(dcl.0[, .(jcl.0, listid)]), by = "jcl.0")
                           if (!weak)
                             altedges[as.character(dcl$V1)]$mark(ecycle = dcl$listid)
                           altedges[as.character(dcl$V1)]$mark(ecluster = dcl$listid)
                           ## mark edge id's with ecluster id's
                           meta = merge(meta,
                             altedges$dt[
                              ,.(edge.id,
                                ecluster,
                                hcl.1,
                                hcl.2,
                                ehcl,
                                ecycle)], by = "edge.id")
                           meta.thresh = merge(meta.links[bp.dist < thresh], altedges$dt[, .(edge.id, ecluster, hcl.1, hcl.2, ehcl, ecycle)], by = "edge.id")
                           meta.links = rbind(meta.thresh, meta.links[!na2false(bp.dist < thresh)], fill = T)


                           meta[, `:=`(
                             num_positive = sum(sign[!duplicated(idx)] > 0, na.rm = T),
                             num_negative = sum(sign[!duplicated(idx)] < 0, na.rm = T)
                           ),
                           by = ecluster]
                           ##     self$set(recip_bp = meta)
                           
                           meta[!is.na(ecluster),
                             `:=`(
                               nclust = length(unique(edge.id)),
                               all_positive = all(replace(sign, is.na(sign),  3e9) > 0),
                               all_negative = all(replace(sign, is.na(sign), -3e9) < 0),
                               mixed = {naom = na.omit(sign); any(naom > 0) & any(naom < 0)},
                               bridge = anyNA(bp.dist)
                             ),
                             by = ecluster]
                           self$set(recip_links = meta.links)
                           recip_bp = copy(meta)
                           eclust.edge = self$edges[!is.na(ecluster)]

                           ndt = data.table::melt(eclust.edge$dt, id.vars = c("ecluster"),
                             measure.vars = c("n1", "n2"))[!duplicated(cbind(ecluster, value))]
                           ndt$nodefp = gr.string(gr.stripstrand(self$nodes[as.character(ndt$value)]$gr)) # footprint of the involved node

                           self$nodes[ndt$value]$mark(ecluster = ndt$ecluster)

                           nnodes = ndt[, .(nnodes = length(unique(value))), by = ecluster]
                           recip_bp = merge(recip_bp[!is.na(ecluster)], nnodes, by = "ecluster")


                           recip_bp[intersect((GenomicRanges::order(gr.stripstrand(dt2gr(recip_bp)))),
                             which(!is.na(sign))), `:=`(fstart = start[1], fend = start[2]), by = .(idx)]
                           recip_bp[!is.na(fstart), `:=`(footprint = (unique(paste0(unique(paste0(seqnames, ":", fstart, "-", fend, ifelse(sign > 0, "+", "-"))), collapse = ",")))), by = .(ecluster)] # footprint of the reciprocal bridge
                           recip_bp[, footprint := na.omit(footprint)[1], by = .(ecluster)]

                           recip_event = recip_bp[
                            ,.(njuncs = nclust[1], nnodes = nclust[1],
                              num_positive = num_positive[1], num_negative = num_negative[1],
                              all_positive = all(replace(sign, is.na(sign),  3e9) > 0),
                              all_negative = all(replace(sign, is.na(sign), -3e9) < 0),
                              mixed = unique(na.omit(mixed)),
                              bridge = unique(na.omit(bridge)),
                              footprint = footprint[1]), by = ecluster]

                           edt = copy3(self$edges$dt)[
                           , `:=`(unfused.n1 = n1 + ifelse(n1.side == "right", 1L, -1L),
                             unfused.n2 = n2 + ifelse(n2.side == "right", 1L, -1L))]
                           edt$fp1 = self$nodes[as.character(edt$n1)]$gr %>% gr.stripstrand %>% gr.string
                           edt$fp2 = self$nodes[as.character(edt$n2)]$gr %>% gr.stripstrand %>% gr.string
                           edt$unfused.fp1 = self$nodes[as.character(edt$unfused.n1)]$gr %>% gr.stripstrand %>% gr.string
                           edt$unfused.fp2 = self$nodes[as.character(edt$unfused.n2)]$gr %>% gr.stripstrand %>% gr.string


                           recip_bp = merge(recip_bp,
                             edt,
                             by = "edge.id",
                             all.x = T,
                             suffixes = c("", ".y_remove_"))##  %>%
                             ## dplyr::select(-dplyr::matches(".y_remove_"))
                           rmcols = grepl(".y_remove_", colnames(recip_bp))
                           recip_bp = recip_bp[,!rmcols, with = F]


                           ufp = recip_bp %>% melt(id.vars = c("ecluster"), measure.vars = c("unfused.fp1", "unfused.fp2"))

                           ## ufp2 = with(ufp, parse.gr(value, meta = g2())) %>% gr.spreduce(ecluster) %>% within({unrfp = gr.string(g2())}) %>% gr2df
                           ufp2 = parse.gr(ufp$value, meta = ufp) %>% gr.spreduce(ecluster)
                           ufp2$unrfp = gr.string(ufp2)
                           ufp2 = gr2dt(ufp2)
                           ufp3 = ufp2[!duped(ecluster, unrfp)][, .(unrfp = paste(sort(unique(unrfp)), collapse = ",")), by = .(ecluster)]
                           un = recip_bp %>% melt(id.vars = c("ecluster"), measure.vars = c("unfused.n1", "unfused.n2"))
                           un2 = un[!duped(ecluster, value)][, paste(sort(unique(value)), collapse = ","), by = .(ecluster)]
                           un2 = setcols(un2, "V1", "un.snode.ids")
                           recip_bp = merge.repl(merge.repl(recip_bp, ufp3, by = "ecluster"), un2, by  = "ecluster")
                           self$set(recip_bp = recip_bp)
                           
                           this = ndt[order(ecluster)][GenomicRanges::order(parse.gr(nodefp))]
                           ndredf = parse.gr(this$nodefp, meta = this) %>% gr.spreduce(ecluster)
                           ## ndredf = with(this, parse.gr(nodefp, meta = g2())) %>% gr.spreduce(ecluster) # reduced involved node footprint

                           dt = data.table(ecluster = ndredf$ecluster, wid = width(ndredf), tmp = gr.string(ndredf))[, .(nodefp = paste(tmp, collapse = ","), nrwid = sum(wid), nrfp = .N), keyby = ecluster] # nrfp is the number of reduced footprints
                           
                           recip_event = merge(recip_event, dt, by = "ecluster") %>% merge(this[, .(snode.ids = paste(value, collapse = ",")), keyby = ecluster], by = "ecluster")
                           recip_event = merge.repl(merge.repl(recip_event, ufp3, by = "ecluster"), un2, by  = "ecluster")
                           self$set(recip_event = recip_event)
                         }

                         if (verbose)
                           message(sprintf('Annotated %s junction cycles in edge field $ecycle', length(jcl)))                         
                         
                         if (paths & !weak)
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
                               lapply(all.paths(tmp.adj, sources = sources, sinks = sinks, verbose = verbose)$paths, function(i) x[i])
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
                         return(invisible(self))
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
                       paths = function(query,
                                        subject = query,
                                        mc.cores = 1,
                                        weight = NULL,
                                        meta = NULL,
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
                                   v  = lapply(p$vpath, function(x) {y = as.integer(x); if (length(y) == 0 | is.na(y[1]) | y[1] != query) return(c()) else return(snode.id[y])})
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
                         else
                         {
                           queryf = self$nodes
                         }

                           ## primitive function to compute distance on graph, called by
                           ## both gNode and GRanges queries below 
                           .dist = function(gg, query, subject, weight = NULL, ignore.strand = TRUE)
                           {
                               ## get igraph and (populated with edge weights)
                               G = gg$igraph
                             
                             if (!nrow(gg$sedgesdt))
                             {
                               return(matrix(Inf, nrow = length(gg$nodes), ncol = length(gg$nodes)))
                             }

                               if (!is.null(weight)) ## replace standard weight with user provided column if provided
                               {
                                   weight = weight[1]
                                   if (!(weight %in% names(gg$sedgesdt)))
                                       stop(paste('column', weight, 'not found in edges metadata'))
                                   
                                   igraph::E(G)$weight = gg$sedgesdt[
                                                   .(E(G)$sedge.id),weight, with = FALSE][[1]]
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

                         ## if (nrow(self$edgesdt)>0){
                             ed = self$edgesdt[,.(n1,n2,n1.side,n2.side,type)]
                         ## } else {
                         ##     ed = data.table(n1 = numeric(0), n2 = numeric(0), n1.side = character(0), n2.side = character(0))
                         ## }


                         ## if no edges, then infinite distance 
                         if (!nrow(ed))
                           return(matrix(Inf, nrow = length(query), ncol = length(subject)))

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
                         ## check if graph is a DAG (non-DAGs break topo-sort in newer igraph versions)
                         if (igraph::is_dag(self$igraph)) {
                             v = unique(abs(self$gr$snode.id[as.vector(
                                 suppressWarnings(
                                     igraph::topo_sort(self$igraph)))]))
                         } else {
                             v = unique(abs(self$gr$snode.id[as.vector(
                                 suppressWarnings(
                                     igraph::topo_sort(igraph::mst(self$igraph))))]))
                         }

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
#                           NONO.FIELDS = c('node.id', 'snode.id', 'index', 'loose.left', 'loose.right', 'loose.left', 'loose.right')
                           NONO.FIELDS = c('node.id', 'snode.id', 'index')

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
                          
                           gr.dt = as.data.table(private$pnodes)
                           suppressWarnings(gr.dt[index, paste(colName) := data])
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


                           suppressWarnings(private$pedges[.(id2), paste(colName) := rep(data, length.out = .N)])                           
                         } else {
                           stop("Not sure how we got to this error at all, we should never be here")
                         }
                         return (self)
                       },
                       
                       
                       #' @name maxflow
                       #' @description
                       #'
                       #' Computes the "max flow" between every node pair in self
                       #' for some metadata field.  
                       #'
                       #' The "max flow" for a node pair i, j is the
                       #' maximum value m of node and/or edge metadata for which there
                       #' exists a path p between i and j whose nodes n and/or edges e
                       #' obey field(n)>=m and/or field(e)>=m for all n,e \eqn{\in} p.
                       #' (i.e. m is the maximum lower bound of the value of
                       #' nodes / edges across all paths connecting ij)
                       #'
                       #' The "min version" of this problem (max = FALSE) will 
                       #' determine the min value m for which there exists  p
                       #' whose nodes n and edges e
                       #' obey field(n)>=m and/or field(e)>=m for all n,e \eqn{\in} p.
                       #'
                       #' The user can also do the problem with lower.bound = FALSE
                       #' i.e. where m is the (maximum or minimum) <upper> bound 
                       #' value in each path. 
                       #'
                       #' By default will try to solve problem across both node and
                       #' edge metadata if the field is present in either.  If the field
                       #' is only present in one then will solve for that.  This property
                       #' can be toggled using edges.only and nodes.only parameters.
                       #'
                       #' @param field metadata field to run maxflow on
                       #' @param max logical flag whether to find maximum path or (if max = FALSE) minimum path
                       #' @param walk if TRUE will return the single walk that maximizes the sum of metadata fields 
                       #' @param nfield field to specify a node field to maximize across paths
                       #' @param efield field to specify an edge field to maximize across paths
                       #' @param cfield field to specify a node / edge field that limits / caps  the dosage at nodes / edges
                       #' @param path.only logical flag relevant only if walk = TRUE,  if path is TRUE will only allow path based maxflows (TRUE) ie will not return a solution when the graph contains only cycles
                       #' @param multi logical flag (FALSE) if TRUE will allow the optimization to compute a solution that outputs multiple disjoint paths
                       #' @param  ncopies positive integer representing the number of copies of the flow that we want the graph to support
                       #' @param reverse complement will compute maximum flow between each node i and the reverse complement of node j in a strand specific way 
                       #' @author Marcin Imielinski
                       maxflow = function(field = NA, walk = FALSE, max = TRUE,
                                          lower.bound = TRUE,
                                          nfield = NA,
                                          efield = NA,
                                          cfield = NA, 
                                          path.only = TRUE,
                                          require.nodes = NULL,
                                          multi = FALSE,
                                          ncopies = 1, 
                                          reverse.complement = FALSE,
                                          verbose = FALSE
                                          )
                       {
                         do.edges = FALSE
                         do.nodes = FALSE

                         if (is.na(field))
                           field = 'cn'

                         if (is.na(cfield))
                           cfield = field

                         if (is.na(efield) & is.na(nfield))
                           efield = nfield = field
                           
                         if (!is.na(efield) && efield %in% names(self$edges$dt))
                           do.edges = TRUE
                         
                         if (!is.na(nfield) && nfield %in% names(values(self$nodes$gr)))
                           do.nodes = TRUE

                         if (!do.edges & !do.nodes)
                         {
                           warning(sprintf('field %s not found in node and field %s not found in edge metadata, marking each with dummy values (1)', efield, nfield))
                           nval = eval = list(1)
                           names(eval) = efield
                           names(nval) = nfield
                           do.call(self$nodes$mark, nval)
                           do.call(self$edges$mark, eval)
                         }

                         ## if (walk)
                         ## {
                         ##   if (is.null(self$edges$dt$cn) | is.null(self$nodes$dt$cn))
                         ##   {
                         ##     do.edges = FALSE
                         ##     do.nodes = FALSE

                         ##   if (is.null(self$nodes$dt$loose.cn.left) | is.null(self$nodes$dt$loose.cn.right))
                         ##   {
                         ##     warning('loose cn left and right is currently required for maxflow(walk = TRUE), putting in dummy values using $loose.left and $loose.right node features')

                         ##     if (is.na(field))
                         ##       field = 'cn'

                         ##     if (is.na(cfield))
                         ##       cfield = field

                         ##     if (is.na(efield) & is.na(nfield))
                         ##       efield = nfield = field
                               
                         ##     if (!is.na(efield) && efield %in% names(self$edges$dt))
                         ##       do.edges = TRUE
                             
                         ##     if (!is.na(nfield) && nfield %in% names(values(self$nodes$gr)))
                         ##       do.nodes = TRUE

                         ##     if (!do.edges & !do.nodes)
                         ##     {
                         ##       warning(sprintf('field %s not found in node and field %s not found in edge metadata, marking each with dummy values (1)', efield, nfield))
                         ##       nval = eval = list(1)
                         ##       names(eval) = efield
                         ##       names(nval) = nfield
                         ##       do.call(self$nodes$mark, nval)
                         ##       do.call(self$edges$mark, eval)
                         ##     }

                         if (walk)
                             {
                               if (is.null(self$edges$dt$cn) | is.null(self$nodes$dt$cn))
                               {
                                 warning('cn is required for maxflow(walk = TRUE), putting in dummy values')
                                 
                                 if (is.null(self$edges$dt$cn))
                                   self$edges$mark(cn = 1)
                                 
                                 if (is.null(self$nodes$dt$cn))
                                   self$nodes$mark(cn = 1)
                               }
                               
                               if (is.null(self$nodes$dt$loose.cn.left) | is.null(self$nodes$dt$loose.cn.right))
                               {
                                 warning('loose cn left and right is currently required for maxflow(walk = TRUE), putting in dummy values using $loose.left and $loose.right node features')
                                 
                                 self$nodes$mark(loose.cn.left = self$nodes$dt$cn*sign(self$nodes$dt$loose.left))
                                 self$nodes$mark(loose.cn.right = self$nodes$dt$cn*sign(self$nodes$dt$loose.right))
                               }


                               if (any(is.na(self$edges$dt$cn)))
                                 stop('cn has NA values, cn is a required edge data field for maxflow(walk = TRUE) so either make blank or fill in with non NA values')

                               if (any(is.na(self$nodes$dt$cn)))
                                 stop('cn has NA values, cn is a required node data field for maxflow(walk = TRUE) so either make blank or fill in with non NA values')
                               
                               ed = copy(private$pedges)

                               ## graph needs loose ends to output a walk
                               ## since all walks must begin and end at a loose end
                               if (!nrow(ed) & !length(self$loose))
                                 return(gW(c(), graph = self))
                               
                               ## make incidence matrix nodes x edges + loose ends
                               Inc = sparseMatrix(1, 1, x = 0,
                                                  dims = c(length(self$nodes),
                                                           length(self$edges) + length(self$loose))*2)
                               
                               rownames(Inc) = c(1:length(self$nodes), -(1:length(self$nodes)))
                               colnames(Inc) = 1:ncol(Inc)
                               
                               colnames(Inc)[1:nrow(ed)] = ed$sedge.id
                               Inc[cbind(ed$from, 1:nrow(ed))] = -1
                               Inc[cbind(ed$to,1:nrow(ed))] = 1 + Inc[cbind(ed$to,1:nrow(ed))]
                                                                                       
                               lleft = which(self$nodes$loose.left)
                               lright = which(self$nodes$loose.right)

                               Inc[cbind(match(lleft, rownames(Inc)),
                                         nrow(ed) + 1:length(lleft))] = 1
                               Inc[cbind(match(-lleft, rownames(Inc)),
                                         nrow(ed) + length(lleft) + 1:length(lleft))] = -1
                               Inc[cbind(match(lright, rownames(Inc)),
                                         nrow(ed) + 2*length(lleft)+ 1:length(lright))] = -1
                               Inc[cbind(match(-lright, rownames(Inc)),
                                         nrow(ed) + 2*length(lleft) + length(lright) + 1:length(lright))] = 1

                               ## metadata for edges including loose ends
                               meta = data.table(index = 1:ncol(Inc),
                                                 sedge.id = NA_integer_,
                                                 lledge.id = NA_integer_, ## left loose end
                                                 lredge.id = NA_integer_ ## right loose end
                                                 )
                               meta[1:nrow(ed), sedge.id := ed$sedge.id %>% as.integer]

                               if (length(lleft))
                                 meta[1:(2*length(lleft)) + nrow(ed), lledge.id := c(lleft, -lleft)]

                               if (length(lright))
                                 meta[1:(2*length(lright)) + nrow(ed) + 2*length(lleft), lredge.id := c(lright, -lright)]

                               meta$cn = ed[.(meta$sedge.id), cn]
                               meta[!is.na(lledge.id), cn := self$nodes[lledge.id]$dt$loose.cn.left]
                               meta[!is.na(lredge.id), cn := self$nodes[lredge.id]$dt$loose.cn.right]

                               ## these ids will allow us to group rc edges
                               meta[, id := ifelse(!is.na(sedge.id), abs(sedge.id),
                                            ifelse(!is.na(lledge.id), paste0(abs(lledge.id), 'l'),
                                                   paste0(abs(lredge.id), 'r')))]
                               meta[, id := factor(id) %>% as.integer]
                               cvec = rep(0, ncol(Inc))
                               
                               ## do.nodes determines whether to create the objective from a node or edge metadata
                               if (do.nodes)
                               {
                                 cvec[1:nrow(ed)] = values(private$pnodes)[ed$from,][[nfield]]
                               }
                               else
                               {
                                 cvec[1:nrow(ed)] = ed[[efield]]
                               }

                               ## reward the use of loose ends
                               lix = meta[is.na(sedge.id), index]
                               cvec[lix] = 1

                               Amat = Inc
                               b = data.table(
                                 type = 'flow',
                                 bvec = rep(0, nrow(Amat)),
                                 sense = 'E'
                               )
                               
                               ## if multi we don't place any constraints on the # of loose ends ie paths
                               if (!multi)
                               {
                                 ## add loose end constraint ie require total weight 2 on loose ends
                                 lec = ifelse(1:ncol(Inc) %in% 1:nrow(ed), 0, 1) ## lec = 1 if loose end, 0 otherwise
                                 Amat = rbind(Inc, lec)
                                 if (path.only) ## require single path
                                 {                               
                                   b = rbind(b,
                                             data.table(
                                               type = 'pathonly',
                                               bvec = 2,
                                               sense = 'E'))
                                 }
                                 else ## allow max one path
                                 {
                                   b = rbind(b,
                                             data.table(
                                               type = 'pathonly',
                                               bvec = 2,
                                               sense = 'L')) 
                                 }
                               }

                           ## add basic utilization constraints
                           ## limiting the flow through each signed node to 1
                           ## each edge utilizes half of each node
                           Umat = sparseMatrix(1,1, x = 0,
                                               dims = c(length(self$nodes),
                                                        length(self$edges) + length(self$loose))*2)
                           
                           rownames(Umat) = c(1:length(self$nodes), -(1:length(self$nodes)))
                           colnames(Umat) = 1:ncol(Umat)                           
                           colnames(Umat)[1:nrow(ed)] = ed$sedge.id

                           Umat.so = Umat.si = Umat ## make separate source and sink versions

                           ## Umat.so catalogues the edges for which that node is a source
                           ## Umat.si ... for which that node is a sink
                           Umat.so[cbind(ed$from, 1:nrow(ed))] = 1
                           Umat.si[cbind(ed$to,1:nrow(ed))] = 1

                           Amat = rbind(
                             Amat, 
                             Umat.so,
                             Umat.si
                           )


                           ## this constrains each node to be used at most once
                           ## as a sink or a source
                           b = rbind(b,
                                     data.table(
                                       type = rep(
                                         c('utilization.so', 'utilization.si'),
                                           each = nrow(Umat)),
                                       bvec = rep(1, 2*nrow(Umat)),
                                       sense = 'L')
                                     )

                           ub = rep(1, ncol(Amat))

                           ## add node and edge capacity constraints if cfield is provided
                           if (!is.null(cfield) & cfield %in% names(self$nodes$dt))
                           {
                             ## we want to figure out the
                             ## strand free usage of a node in any flow
                             ## meaning add the number of times that node
                             ## "side" is used in any flow
                             ## meaning that we want to add the usage of
                             ## node sinks and the sources of their reverse complements
                             ## ie have a constraint that
                             ## limits the sum of copy number across the union of
                             ## sink_edges(node_i) and source_edges(rc(node_i))
                             ## to the node copy number
                             Umatc.left = Umat.si[self$nodes$dt$node.id %>% as.character, ] +
                               Umat.so[-self$nodes$dt$node.id %>% as.character, ]
                             Umatc.right = Umat.so[self$nodes$dt$node.id %>% as.character, ] +
                               Umat.si[-self$nodes$dt$node.id %>% as.character, ]
                                                          
                             Amat = rbind(
                               Amat, 
                               Umatc.left,
                               Umatc.right
                             )
                             
                             ## cap the signed node usage to node copy number
                             ## divided by ncopies
                             b = rbind(b,
                                       data.table(
                                         type = 'capacity.left', 
                                         bvec = self$nodes$dt[[cfield]]/ncopies,
                                         sense = 'L'),
                                       data.table(
                                         type = 'capacity.right', 
                                         bvec = self$nodes$dt[[cfield]]/ncopies,
                                         sense = 'L')
                                       )

                             if (!is.null(require.nodes))
                             {
                               Amat = rbind(
                                 Amat, 
                                 Umatc.left[require.nodes,,drop = FALSE] +
                                 Umatc.right[require.nodes,,drop = FALSE]                           
                               )
                               
                               b = rbind(b,
                                         ## data.table(
                                         ##   type = 'require.left', 
                                         ##   bvec = rep(1, length(require.nodes)),
                                         ##   sense = 'G'),
                                         data.table(
                                           type = 'require.right', 
                                           bvec = rep(1, length(require.nodes)),
                                           sense = 'G')
                                         )
                             }
                             
                             ## now strand collapse edges
                             Emat = sparseMatrix(
                               meta$id,
                               meta$index, x = as.numeric(1), dims = c(nrow(meta)/2,ncol(Amat))) %>% as.matrix
       
                             ## this matrix now summarizes the copy number of this flow
                             ## across pairs of stranded edges
                             ## which can't be bigger than capacity constraints
                             ## on the unsigned edge
                             
                             Amat = rbind(
                               Amat, 
                               Emat
                             )
                             
                             ## cap the signed node usage to node copy number
                             ## these capacity constraints are stored in meta
                             b = rbind(b,
                                       data.table(
                                         type = 'ecapacity', 
                                         bvec = meta[match(1:nrow(Emat), id), cn]/ncopies,
                                         sense = 'L')
                                       )

                             ## cap individual signed edge.ids to their capacity
                             ub[1:nrow(ed)] = pmin(ub[1:nrow(ed)], ed[[cfield]])
                           }
                           ## sol = Rcplex::Rcplex(Amat = Amat)
                           sol = Rcplex2(Amat = Amat,
                                          lb = rep(0, ncol(Amat)),
                                          ub = ub,
                                          bvec = b$bvec,
                                          sense = b$sense,
                                          vtype = 'I', 
                                          control = list(trace = verbose > 1),
                                          objsense = ifelse(max, 'max', 'min'),
                                          cvec = cvec)

                           ## find edges used in the solution
                           opted = ed[round(sol$xopt[1:.N])!=0, ]

                           ## now need to convert pile of edges to walk(s)
                           ## (note: there can be more than one walk because of cycles)
                           ## (however the edge pile is constrained because of the above
                           ## MIP to have at most one edge entering and leaving each
                           ## node)
                           ## will use dfs
                           
                           ## pick a source (or random node if none exists) for the DFS
                           ## find source nodes by mining the loose ends that were utilized
                           ## add also additional sources from the used edges
                           lix = meta[is.na(sedge.id), index]
                           so = which(((Inc[, lix]>0) %*% sol$xopt[lix])!=0) %>%                           
                             c(opted$from)%>% unique %>% as.character

                           ## combine with all other utilized vertices to build graph
                           allv = union(so, opted$to)

                           if (!length(so))
                             return(gW(graph = self))

                           ## we have to make adjacency matrix here to take care of edge cases
                           ## where there are graphs with no edges but loose ends
                           adj = sparseMatrix(1, 1, x = 0, dim = rep(length(allv),2), dimnames = rep(list(allv), 2))
                           adj[cbind(match(opted$from, allv), match(opted$to, allv))] = 1
                           G = graph.adjacency(adj)

                           ## we are using dfs to convert this flow into paths / cycles
                           ## we do a separate dfs for each source, prioritizing the true "source" nodes first
                           ## cyclic paths will begin from an arbitrary internal source 
                           ## (this will create some redundant walks, which we prune out below)
                           walksi = lapply(so, function(x) V(G)$name[dfs(G, x, unreachable = FALSE)$order %>% as.integer %>% setdiff(NA)] %>% as.integer)

                           ## for each node, we only keep the first walk that has that node
                           ## (this will remove reverse complements as well as redundant paths from the dfs)
                           snode.id = lapply(walksi, function(x) private$pnodes$snode.id[x])
                           keep = dunlist(snode.id)[!duplicated(abs(V1)), listid] %>% unique
                           walksi = walksi[keep]
                           snode.id = snode.id[keep]

                           ## mark all cycles
                           setkeyv(opted, c('from', 'to'))
                           is.cycle = opted[lapply(walksi, function(x) data.table(x[length(x)], x[1])) %>% rbindlist, !is.na(edge.id)]

                           ## now first let's instantiate a graph around opted so we can get the edge lists
                           tmp.gg = gG(nodes = self$nodes$gr, edges = self$edges[abs(opted$sedge.id)]$dt)
                           tmp.gw = gW(snode.id = snode.id, circular = is.cycle, graph = tmp.gg)
                           
                           sedge.id = lapply(tmp.gw$sedge.id, function(x) if (length(x)) sign(x)*abs(opted$sedge.id)[abs(x)])

                           gw = gW(sedge.id = sedge.id, snode.id = snode.id, circular = is.cycle, graph = self)
                           
                           if (path.only) ## keep only the paths
                           {
                             gw = gw[!gw$circular]
                           }

                           return(gw)
                         } 
                         
                         uval = c()
                         if (do.edges)
                           uval = c(uval, self$edges$dt[[nfield]])
                         if (do.nodes)
                           uval = c(uval, values(self$nodes$gr)[[nfield]])

                         uval = sort(unique(setdiff(uval, NA)))

                         ## if max we will begin with high values
                         ## and determine the
                         ## connectivity of graph nodes at all nodes / edges with
                         ## field>=this.val

                         ## if !max we will begin low values and determine
                         ## the connectivity of all nodes / edges with field<=this.val
                         if (lower.bound)
                           uval = rev(uval)

                         out = matrix(ifelse(max, -Inf, Inf), nrow = length(self$nodes),
                                      ncol = length(self$nodes))

                         if (length(uval)==0)
                           return(out)
                         
                         comparator = ifelse(lower.bound, '>=', '<=')
                         estring = ifelse(do.edges,
                                          sprintf('%s %s this.val', efield, comparator),
                                          '')
                         nstring = ifelse(do.nodes,
                                          sprintf('%s %s this.val', nfield, comparator),
                                          '')
                         cmd = sprintf("graph[%s,%s]", nstring, estring)

                         ## copy self and label og nodes 
                         graph = self$copy
                         graph$nodes$mark(og.node.id = graph$nodes$dt$node.id)
                         for (i in 1:length(uval))
                         {
                           this.val = uval[i]
                           if (verbose)
                             message(sprintf('trying values %s %s (%s of %s)',
                                             comparator, this.val,
                                             i, length(uval)))
                           this.graph = eval(parse(text = cmd))
                           node.ids = this.graph$nodes$dt$og.node.id
                           if (reverse.complement)
                           {
                             nix = this.graph$nodes %>% seq_along
                             mat = pmin(
                               this.graph$dist(nix, -nix, ignore.strand = FALSE),
                               this.graph$dist(-nix, nix, ignore.strand = FALSE))
                             reachable = as.matrix(
                               as.data.table(
                                 which(mat<Inf,
                                       arr.ind = TRUE))[,
                                                        .(nid1 = node.ids[row], nid2 = node.ids[col])])
                           }
                           else
                           {
                             reachable = as.matrix(
                               as.data.table(
                                 which(this.graph$dist()<Inf,
                                       arr.ind = TRUE))[,
                                                        .(nid1 = node.ids[row], nid2 = node.ids[col])])
                           }
                           if (max)
                             out[reachable] = pmax(out[reachable], this.val)
                           else
                             out[reachable] = pmin(out[reachable], this.val)
                         }
                         return(out)
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


                       gtrack = function(y.field = NULL, lwd.loose = 3, col.loose = alpha("blue",0.6), col.alt = alpha("red",0.4), ...)
                       {
                         if (!length(self$nodes))
                           return(gTrack(GRanges(seqlengths = seqlengths(self))))


                         ss = private$pnodes
                         ed = private$pedges
                         
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


                         if (!is.null(ed))
                         {

                           ## set edge apperaances
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

                           ed[, ":="(lwd = ifelse(is.na(lwd), ifelse(type=="ALT", log2(0.2*pmax(0, y, na.rm = TRUE)+2)+1, ifelse(type == 'loose', lwd.loose, 1)),lwd),
                                     lty = ifelse(is.na(lty), ifelse(type=='loose', 1, 1),lty),
                                     col = ifelse(is.na(col), ifelse(is.na(col), ifelse(type=="ALT",
                                                                                 ifelse(is.na(y),
                                                                                        col.alt,
                                                                                        ifelse(y>0,
                                                                                               col.alt,
                                                                                               alpha("purple", 0.3))),
                                                                                 ifelse(type=="loose",
                                                                                        col.loose,
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

                             mx = as.data.table(ss)[, max(y), keyby = .(seqnames, parent.graph)]
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
                       #' @param tile interval around which to trim the gGraph
                       #' @return new gGraph trimmed to tile, unless mod is set to TRUE
                       #'
                       #' @details
                       #' ```
                       #' gr = c(GRanges("1", IRanges(10000,100000), "+"), GRanges("2", IRanges(10000,100000), "+"))
                       #' new.gg = gg$trim(gr)
                       #' ```
                       #' 
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
                           es[, keep := n1 %in% validNodes & n2 %in% validNodes]
                           
                           ## Want to remove edges that went into nodes that were trimmed on one side
                           ## Find query.id's that have all FALSE in left or right - might still be multiple trims but one side is trimmed
                           sg = gr2dt(new.nodes)[, .(mean(left)>0, mean(right)>0), by=query.id]
                           leftRemove = sg[V1 == FALSE, query.id]
                           rightRemove = sg[V2 == FALSE, query.id]
                           
                           ## Remove edges that had one of their sides trimmed                           
                           es = es[, keep := keep & !(n1 %in% leftRemove & n1.side == 0) & !(n2 %in% leftRemove & n2.side == 0)]
                           es = es[, keep := keep & !(n1 %in% rightRemove & n1.side == 1) & !(n2 %in% rightRemove & n2.side == 1)]

                           new.es = es[keep == TRUE, ]
                           ## keep track of node sides that will now acquire loose ends
                           new.loose = unique(es[keep == FALSE, .(n = c(n1, n2), side = c(n1.side, n2.side))][n %in% (new.nodes$node.id), ])

                           if (nrow(new.loose)>0)
                           {
                             new.loose$present = TRUE
                             setkeyv(new.loose, c("n", "side"))
                             new.nodes$loose.right = new.nodes$loose.right |  new.loose[.(new.nodes$node.id, 1), !is.na(present)]
                             new.nodes$loose.left = new.nodes$loose.left |  new.loose[.(new.nodes$node.id, 0), !is.na(present)]
                           }




                           ## map left==TRUE to n1 or n2 side == 0
                           ## map right==TRUE to n1 or n2 side == 1
                           map = data.table(old = c(new.nodes[new.nodes$left]$query.id, new.nodes[new.nodes$right]$query.id),
                                            side = c(rep(0, length(which(new.nodes$left))), rep(1, length(which(new.nodes$right)))),
                                            new = c(which(new.nodes$left), which(new.nodes$right)))
                           setkeyv(map, c("old", "side"))
                           new.es[, ":="(n1 = map[.(n1,n1.side), new], n2 = map[.(n2,n2.side), new])]

                           es = new.es
                         }


                         ## Remove miscellaneous metacols added in this function
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

                       #' @name fix
                       #' @description
                       #'
                       #' Modifies (in place) the current seqlevels of the gGraph, including
                       #' keeping only certain seqlevels, dropping certain seqlevels, and replacing seqlevels.
                       #' 
                       #' Warning: this may modify the graph including getting rid of nodes and edges (i.e. those outside
                       #' the retained seqlevels) and also change coordinates (ie move ranges that were previously on different
                       #' chromosomes to the same chromosome etc.).  Use with caution!
                       #'
                       #' Default behavior is to replace 'chr', with ''.
                       #' 
                       #' @param pattern  character pattern to replace in seqlevels (used in a gsub, can have backreferences)
                       #' @param replacement  character to replace pattern with (used in a gsub, can have backreferences) 
                       #' @param drop  character vector of seqlevels to drop or logical TRUE to drop all seqlevels that are unused  (TRUE)
                       #' @param seqlengths new seqlengths i.e. named integer vector of seqlevels to drop or embed graph into 
                       #' @author Marcin Imielinski
                       #' @return current graph modified in place with additional nodes and edges, as specified by user 
                       fix = function(pattern = NULL, replacement = NULL, drop = TRUE, seqlengths = NULL)
                       {
                         if (is.logical(drop))
                         {
                           if (drop)                             
                             drop = setdiff(seqlevels(self), as.character(seqnames(private$pnodes)))
                           else
                             drop = NULL
                         }

                         if (!is.null(drop) | !is.null(seqlengths))
                         {
                           tmp = as.data.table(self$nodes$dt)

                           if (!is.null(seqlengths))
                             sl.new = seqlengths
                           else
                             sl.new = seqlengths(self)

                           ix = which(tmp$seqnames %in% names(sl.new) & !(tmp$seqnames %in% drop))
                           sl.new = sl.new[setdiff(names(sl.new), drop)]

                           gg = self                           
                           gg = self[ix, ]
                           sl.new = sl.new[setdiff(names(sl.new), drop)]
                           if (length(sl.new))
                             private$pnodes = dt2gr(as.data.table(gg$gr), seqlengths = sl.new)
                           else
                             private$pnodes = GRanges(seqlengths = sl.new)
                           private$pedges = gg$sedgesdt
                         }

                         if (!is.null(pattern) & !is.null(replacement))
                         {
                           tmp = as.data.table(private$pnodes)
                           tmp$seqnames = gsub(pattern, replacement, tmp$seqnames)
                           sl.new = seqlengths(self)
                           names(sl.new) = gsub(pattern, replacement, names(sl.new))
                           sl.new = data.table(sn = names(sl.new), sl = sl.new)[, max(sl), by = sn][, structure(V1, names = sn)]
                           private$pnodes = dt2gr(tmp, seqlengths = sl.new)
                         }

                         private$buildLookupTable() ## essential for graph integrity!!
                         private$stamp()
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
                             junctions = jJ(junctions)

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
                       #' annotations are node / edge features that will be dumped to json
                       #' @param filename character path to save to
                       #' @param save whether to save or return list object representing json contents
                       #' @param nfields which node fields to dump to json (NULL)
                       #' @param efields which edge fields to dump to json (NULL)
                       #' @param annotations which graph annotations to dump to json
                       #' @param settings gGnome.js settings values to add to the output JSON files (list)
                       #' @param cid.field field in the graph edges that should be used for setting the cid values in the JSON (default: 'sedge.id'). 
                       #' This is useful for cases in which there is some unique identifier used across samples to identify 
                       #' identical junctions (for example "merged.ix" field, which is generated by merge.Junction())
                       #' @author Marcin Imielinski
                       json = function(filename='.',
                                       maxcn=100,
                                       maxweight=100,
                                       save = TRUE,
                                       verbose = FALSE,
                                       annotations = NULL,
                                       nfields = NULL, ## node metadata fields to dump
                                       efields = NULL, ## edge metadata fields to dump
                                       settings = list(y_axis = list(title = "copy number",
                                                                     visible = TRUE)),
                                       cid.field = NULL,
                                       no.y = FALSE)
                       {
                         ## Make sure that our nodes are not empty before visualizing
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
                      
                         ## range of CN
                         ymin=0
                         ymax=maxcn
                         
                         node.json = 
                           gr2dt(self$nodes$gr[, "snode.id"])[, .(chromosome = seqnames, startPoint = start, endPoint = end, iid = snode.id, y = 1)]

                         if (length(nfields))
                         {
                           nfields = setdiff(intersect(nfields, names(self$nodes$dt)), names(node.json))
                           node.json = cbind(node.json, gr2dt(self$nodes$gr)[, nfields, with = FALSE])
                         }
                         .dtstring = function(dt)
                           dt[, gsub('\\|+', '|', gsub('\\|+$', '', gsub('^\\|+', '', do.call(paste, c(lapply(names(.SD), function(x) ifelse(!is.na(.SD[[x]]), paste0(x, '=', .SD[[x]]), '')), sep = '|')))))]

                         if (!is.null(annotations))
                         {
                           nodes.overlap.annotations = intersect(annotations, names(values(self$nodes$gr)))
                           if (length(nodes.overlap.annotations) == 0){
                             warning('There is no overlap between the provided annotations and the annotaions available in the nodes in your gGraph.')
                             node.json$annotation = ''
                           } else {
                               node.json = cbind(node.json, data.table(annotation = .dtstring(as.data.table(values(self$nodes$gr))[, intersect(annotations, names(values(self$nodes$gr))), with = FALSE])))                           
                           }
                         }

                         if (is.null(cid.field)){
                             cid.field = 'sedge.id'
                         } else {
                             if (!(cid.field %in% names(private$pedges))){
                                 stop('Invalid cid.field: "', cid.field, '"')
                             }
                             if (any(is.na(get_cids(self, cid.field)))){
                                 stop('Invliad cid.field: "', cid.field)
                             }
                         }

                         ed = data.table()
                         efields = setdiff(intersect(efields, names(private$pedges)),  c("sedge.id", "class", "from", "to", "type", annotations))
                         if (nrow(private$pedges))
                         {
                           ed = copy(private$pedges)[sedge.id>0, intersect(names(private$pedges), c("sedge.id", "class", "from", "to", "type", cid.field, efields, annotations)), with = FALSE] ## otherwise change by reference!

                           # cid.field is sometimes only defines for ALT edges so here we add values for REF (We also need to add this later to LOOSE edges. See below)
                           cid_na = ed[is.na(get(cid.field)), .N]
                           if (cid_na > 0){
                               max.cid = max(ed[, get(cid.field)], na.rm = T)
                               ed[is.na(get(cid.field)), (cid.field) := (max.cid + .I)]
                           }

                           if (!is.null(annotations)){
                             edges.overlap.annotations = ed[, intersect(names(ed), annotations), with = FALSE]
                             if (length(edges.overlap.annotations) == 0){
                               warning('There is no overlap between the provided annotations and the annotaions available in the edges in your gGraph.')
                               ed$annotation = ''
                             } else {
                               ed$annotation = .dtstring(ed[, intersect(names(ed), annotations), with = FALSE])
                             }
                           ed$from = private$pnodes$snode.id[ed$from]
                           ed$to = -private$pnodes$snode.id[ed$to]
                           }
                         }

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
                           ## add cid values for loose ends
                           max.cid = max(ed[, get(cid.field)], na.rm = T)
                           loose.ed[, (cid.field) := 1:.N + max.cid]
                           loose.ed[, class := '']

                           if (!is.null(annotations))
                             loose.ed$annotation = ''
                           ed = rbind(ed, loose.ed, fill = TRUE)
                         }

                         ## remove nodes and edges that are off the map
                         ##  good.nodes = node.json[chromosome %in% seqlevels, iid]
                         ##  node.json = node.json[iid %in% good.nodes, ]
                         
                         ## ed[, good.count := rowSums(cbind(abs(from) %in% good.nodes, abs(to) %in% good.nodes), na.rm = TRUE)]
                       
                         ##  ed[, node.count := rowSums(cbind(!is.na(from), !is.na(to)))]
                         ##  ed = ed[good.count == node.count, ]

                         # KSK: temporary fix for old-style notation (large values instead of negatives) 
                         # of whether edge is going out/coming in from left or right side of node
                         # so that json is rendered correctly in gGnome.js:
                         # max.node.id <- max(self$dt$node.id)
                         # ed[from > max.node.id, from := (from - max.node.id)*-1]
                         # ed[to <= max.node.id, to := -1*to]
                         # ed[to > max.node.id, to := (to - max.node.id)]
                         # KSK end

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
                                        .(cid = get(cid.field),
                                          source = from,
                                          sink = to,
                                          title = class,
                                          type = type,
                                          weight)]             
                        
                           if (!is.null(annotations))
                             ed.json = cbind(ed.json, ed[, "annotation", with = FALSE])


                           ## append list of edge metadata features if efields is specified
                           if (length(efields))
                           {
                             .fix = function(x) as.list(x)[!is.na(unlist(x))]

                             ## list of main features
                             mainl = lapply(split(ed.json, 1:nrow(ed.json)), .fix)

                             metal = lapply(split(ed[, efields, with = FALSE], 1:nrow(ed.json)), .fix)
                             ## combine the two lists one inside the other
                             ed.json = unname(mapply(function(x,y) c(x, metadata = list(y)), mainl, metal, SIMPLIFY = FALSE))
                           }

                         } else {
                           ed.json = data.table(cid = numeric(0),
                                                source = numeric(0),
                                                sink = numeric(0),
                                                title = character(0),
                                                type = character(0),
                                                weight = numeric(0))
                         }

                         ## append list of node metadata features if nfields is specified
                         if (length(nfields))
                         {
                           .fix = function(x) as.list(x)[!is.na(unlist(x))]
                           mainl = lapply(split(node.json[, setdiff(names(node.json), nfields), with = FALSE], 1:nrow(node.json)), .fix)
                           
                           metal = lapply(split(node.json[, nfields, with = FALSE], 1:nrow(node.json)), .fix)
                           ## combine the two lists one inside the other
                           node.json = unname(mapply(function(x,y) c(x, metadata = list(y)), mainl, metal, SIMPLIFY = FALSE))
                         }

                         gg.js = list(intervals = node.json, connections = ed.json)
                         
                          if (no.y){
                           settings$y_axis = list(visible=FALSE)
                         }

                         if (!is.null(self$meta$name) | !is.null(self$meta$description))
                         {
                           name = paste0('<h3>', self$meta$name, '</h3>')
                           description = paste0('', self$meta$description)
                           settings$description = paste(name, description)
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

                       get.diameter = function(weights = NULL){
                          self$check
                          if (length(self$nodes)==0){
                              return(gW())
                          }
                          if (length(self$edges)>0){
                              tmp.ig = igraph::remove.edge.attribute(self$igraph, "weight")
                          } else {
                              tmp.ig = self$igraph
                          }

                          if (!is.null(weights)){
                              if (length(weights)==1 &&
                                  is.character(weights) &&
                                  is.element(weights, colnames(self$sedgesdt))){
                                  tmp.ig = self$igraph
                                  igraph::set.edge.attribute(
                                      tmp.ig, name="weight",
                                      value = self$sedgesdt[get.edge.attribute(tmp.ig)$sedge.id,
                                                            get("weights")])
                              } else if (length(weights)==nrow(self$sedgesdt)){
                                  names(weights) = self$sedgesdt$sedge.id
                                  igraph::set.edge.attribute(
                                      tmp.ig, name="weight",
                                      value = weights[
                                          as.character(get.edge.attribute(tmp.ig)$sedge.id)])
                              }
                          }
                          diam = igraph::get_diameter(tmp.ig, weights = weights)
                          snd = self$gr[as.numeric(diam)]$snode.id
                          return(gW(snode.id = list(snd), graph = self))
                       },

                       circos = function(...){
                           have.circlize = !inherits(try(find.package("circlize")), "try-error")
                           if (!have.circlize){
                               stop("You need the packages 'circlize' to use this function.")
                           }
                           args  = list(...)
                           ## some important pars
                           labels.cex = ifelse(is.null(args$labels.cex), 1, args$labels.cex)
                           bands.height = ifelse(is.null(args$bands.height), 0.1, args$bands.height)
                           cn.height = ifelse(is.null(args$cn.height), 0.1, args$cn.height)
                           link.h.ratio = ifelse(is.null(args$link.h.ratio), 0.75, args$link.h.ratio)
                           junc.col = setNames(c('#e41a1c','#377eb8','#4daf4a','#984ea3'),
                                               c('DUP-like', "DEL-like", "INV-like", "TRA-like"))
                           div.col3 = c('#67a9cf','#f7f7f7','#ef8a62')
                           juncs = gr2dt(gr.chr(grl.unlist(self$edges[type=="ALT"]$grl)))
                           bp1 = juncs[grl.iix==1][order(grl.ix)]
                           bp2 = juncs[grl.iix==2][order(grl.ix)]
                           cn.bed = self$nodes$dt[!is.na(cn)]
                           if (cn.bed[, all(!grepl("chr", seqnames))]){
                               cn.bed[, seqnames := paste0("chr", seqnames)]
                           }
                           circlize::circos.initializeWithIdeogram(
                               track.height = bands.height,
                               labels.cex = labels.cex)
                           ## add CN
                           f = circlize::colorRamp2(
                               breaks = c(0, 2, cn.bed[, max(cn)]),
                               colors = div.col3)
                           circlize::circos.genomicTrackPlotRegion(
                               ylim = cn.bed[!is.na(cn), range(cn)],
                               track.height = cn.height,
                               cn.bed[, .(seqnames, start, end, cn)],
                               bg.border = NA,
                               panel.fun = function(region, value, ...) {
                                   i = circlize::get.cell.meta.data(
                                       "sector.numeric.index")
                                   circlize::circos.genomicRect(
                                       region,
                                       value,
                                       area = TRUE,
                                       border = f(value),
                                       baseline = 0,
                                       col = f(value),
                                       ...)
                               })
                           ## add junctions
                           circlize::circos.genomicLink(
                               bp1[, 1:3, with=FALSE],
                               bp2[, 1:3, with=FALSE],
                               col = setNames(junc.col[as.character(bp1$class)],
                                              rownames(bp1)),
                               border=NA,
                               h.ratio = link.h.ratio)
                           circlize::circos.clear()
                       },
                       
                       split = function(by = "parent.graph"){
                           if (!is.element(by, colnames(self$nodes$dt))){
                               stop(by, " not found among node metadata!")
                           }
                           val = self$nodes$dt[[by]]
                           uval = unique(val)
                           names(uval) = as.character(uval)
                           if (length(na.ix <- which(is.na(uval))) > 0){
                               names(uval)[na.ix] = "NA"
                           }
                           ggs = lapply(uval, function(x){
                               if (!is.na(x)){
                                   return(self$copy[which(val==x),])
                               } else {
                                   return(self$copy[which(is.na(val)),])
                               }
                           })
                           return(ggs)
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
                       #' @description
                       #' Constructor, initializes an empty gGraph object. If the user does not provide genome
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

                           ## HOW TO ADD empty fields to self$edgesdt
                           ## make sure empty edge data table also has the following default fields
                           ## edges = data.table(n1 = numeric(0), n2 = numeric(0), n1.side = character(0), n2.side = character(0))

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
                         if (!is.null(edges) && nrow(edges)>0)
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
                         
                         if (!is.null(edges) && nrow(edges)>0)
                         {
                           ## FIXME: Need to change this to make sure it labels edges
                           type.provided = TRUE
                           if(!"type" %in% names(edges)) {
                             type.provided = FALSE
                             edges[, type := NA_character_]
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
                           if (!is.null(edges[['from']]))
                             edges$from = NULL

                           if (!is.null(edges[['to']]))
                             edges$to = NULL

                           ## See if there is a cn in the edges, in which case we should carry it over
                           ## is.cn = 'cn' %in% names(edges)
                           ## is.type = 'type' %in% names(edges)
                           
                           tmp = rbind(
                             edges[, .(jid, from = ifelse(n1.side == 1, map[.(n1, '+'), id], map[.(n1, '-'), id]),
                                       to = ifelse(n2.side == 1, map[.(n2, '-'), id], map[.(n2, '+'), id]),
                                       sedge.id = 1:.N)],
                             edges[, .(jid, from = ifelse(n2.side == 1, map[.(n2, '+'), id], map[.(n2, '-'), id]),
                                       to = ifelse(n1.side == 1, map[.(n1, '-'), id], map[.(n1, '+'), id]),
                                       sedge.id = -1*(1:.N))]
                           )

                           cols = c('jid', setdiff(names(edges), names(tmp)))
                           tmp = merge(tmp, edges[, cols, with = FALSE], by = "jid")                           
                           tmp[, edge.id := abs(sedge.id)]
                           
                           
                           ## Drop unused columns
                           tmp[, c("jid","n1","n2","n1.side","n2.side") := NULL]
                           private$pedges = tmp

                           ## if not REF then ALT
                           private$pedges[, type := ifelse(is.na(type), type, ifelse(type == 'REF', 'REF', 'ALT'))]

                           
                           setkey(private$pedges, sedge.id)
                         }

                           private$buildLookupTable()
                           private$stamp()

                           if (is.null(edges) || nrow(edges)==0){
                               self$edges$mark(n1 = numeric(0), n2 = numeric(0), n1.side = character(0), n2.side = character(0))
                           }
                         ## label edges with class and type 
                         self$edges$mark(class = self$edges$class)

                         ## label remaining edges with NA type as REF or ALT based on class
                         self$edges$mark(type = ifelse(is.na(self$edges$dt$type),
                                                ifelse(self$edges$class == 'REF',
                                                       'REF', 'ALT'),
                                                self$edges$dt$type)
                                         )
                        
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

                         if (length(self$gr)==0)
                           return(igraph::graph_from_adjacency_matrix(sparseMatrix(1,1, x = 0)))

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
                         if (!nrow(private$pedges))
                           return(data.table(n1 = numeric(0), n2 = numeric(0), n1.side = character(0), n2.side = character(0)))
                         return(copy(convertEdges(self$gr, private$pedges[.(1:(nrow(private$pedges)/2)), ], metacols = TRUE)[, n1.side := sides[n1.side+1]][, n2.side := sides[n2.side+1]]))
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
                       },
                      
                      diameter = function(){
                          self$check
                          return(self$get.diameter())
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
    if (!length(gr))
      return(gr)
    gr$parent.graph = nm
    gr$og.node.id  = gr$node.id
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
    edt$og.edge.id = edt$edge.id
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

#' @name split.gGraph
#' @title split.gGraph
#' @description
#'
#' Split gGraph by a node metadata field
#'
#' @param by field name to split by, default "parent.graph"
#' @return A new gGraph object that is the union of the nodes and edges in the input gGraphs
#' @author Marcin Imielinski
#' @export
'split.gGraph' = function(gg, by = "parent.graph"){
    return(gg$split(by = by))
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
#' @param junctions BND vcf or bedpe file path, gGraph::Junctions object, or GRangesList specifying pairs of locations to reconnect (can be used in conjunction with breaks)
#' @param alignments GRanges with fields $cigar $flag and $qname (e.g. output of bamUtils::read.bam or RSeqLib::BWA specifying alignments to a reference, the graph will then represent the (possibly multimapped) implicit nodes and edges representing an "end to end" walk of all queries lifted to the reference.  Insertions are represented as 0-width nodes and multi-maps yield branches.  The resulting graph can be walked to exhaustively or greedily enumerate all possible linear embeddings of the queries on the reference
#' @param jabba path to JaBbA .rds output file, if specified instantiates a gGraph object from JaBbA output with y.field "cn" specifying copy number 
#' @param prego path to PREGO output file. if specified, instantiates a gGraph object from PREGO output with y.field "cn" specifying copy number
#' @param weaver path to Weaver output. if specified, instantiates a gGraph object from Weaver output with y.field "cn" specifying copy number
#' @param remixt path to RemiXT output, if specified, instantiates a gGraph object from ReMiXT output with y.field "cn" specifying copy number
#' @param rck path to RCK output, if specified, instantiates a gGraph object from RCK output with y.field "cn" specifying total copy number
#' @param nodes GRanges of unsigned intervals to be rejoined in gGRaph, used in conjunction with edges argument below
#' @param walks GRangesList or gWalk to build a "haplograph" around, see ?haplograph
#' @param edges data.table with field n1, n2, n1.side, n2.side with each row specifying an edge, i.e. which sides ("left" vs "right") of which integer node indices to connect in the gGRaph
#' @param nodeObj gNode object to create a gGraph around (similar to nodes input above), except generated via the $nodes accessor of an existing gGraph object, used in cojunction with gEdge input to create a new gGnome object from an existing one
#' @param edgeeObj gEdge object to create a gGraph around (similar to edges input above), except generated via the $edges accessor of an existing gGraph object
#' @param meta list of metadata to associate with this gGraph object, for now mostly used to populate gTrack visualization parameters
#' @return A new gGraph object that is the sum of the component gGraphs
#' @export
gG = function(genome = NULL,
              breaks = NULL,
              junctions = NULL,
              alignments = NULL,
              juncs = NULL,
              prego = NULL,
              jabba = NULL,
              cougar = NULL,
              weaver = NULL,
              remixt = NULL,
              rck = NULL,
              nodes = NULL,
              edges = NULL,
              walks = NULL,
              nodeObj = NULL,
              edgeObj = NULL,
              meta = NULL
              )
{
  if (!is.null(junctions))
    juncs = junctions;

  return(gGraph$new(genome = genome,
                    breaks = breaks,
                    juncs = juncs,
                    alignments = alignments,
                    prego = prego,
                    jabba = jabba,
                    cougar = cougar,
                    weaver = weaver,
                    remixt = remixt,
                    rck = rck,
                    nodes = nodes,
                    edges = edges,
                    walks = walks, 
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
#' @rdname equals.gGraph
#' @title equals.gGraph
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
#' @rdname equals.gNode
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
#' @rdname equals.gEdge
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
#' @rdname notequals.gGraph
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
#' @rdname notequals.gNode
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
#' @rdname notequals.gEdge
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
#' @rdname plus.gGraph
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
#' @rdname subset.gGraph
#' @description
#'
#' Overloads subset operator for gGraph.
#'
#' Subsetting nodes will remove any edges associated with non-existent nodes
#' Subsetting edges will not remove the associated nodes
#'
#' @param obj  gGraph
#' @param i integer, logical, or expression in gNode metadata used to subset gNodes
#' @param j integer, logical, or expression in gEdge metadata used to subset gEdges
#' @return A new gGraph object that only contains the given nodes and edges
#' @author Marcin Imielinski
#' @export
'[.gGraph' = function(obj, i = NULL, j = NULL, with = TRUE, ...){
  if (any(deparse(substitute(j)) != "NULL"))
  {
    edges = obj$edges[j, with = with]
    if (all(deparse(substitute(i)) == "NULL")){
        nodes = obj$nodes
    }
    else
    {
      nodes = obj$nodes[i, with = with]
      nids = c(nodes$dt$index, nodes$flip$dt$index);
      eid = edges$sdt[to %in% nids & from %in% nids, sedge.id]
      if (length(edges)>0)
      {
        eix = edges$dt[, ix := 1:.N][.(eid), ix]
        edges = edges[eix, with = with]
      }
    }
  } else
  {
    if (all(deparse(substitute(i)) == "NULL")){
      nodes = obj$nodes
    }
    else{
      nodes = obj$nodes[i, with = with]
    }
    nids = c(nodes$dt$index, nodes$flip$dt$index);
    if (length(nids) && length(nodes$edges))
    {
      eid = nodes$edges$sdt[to %in% nids & from %in% nids, sedge.id]
      edges = obj$edges[eid, with = with]
      }
    else
      edges = obj$edges[c()]
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
`length.gGraph` = function(x){
  return(x$length)
}


#' @name length
#' @title length
#' @description
#' The number of walks in the gWalk
#'
#' @param gWalk a \code{gWalk} object
#'
#' @return the number of nodes in the gWalk
#' @export
`length.gWalk` = function(x){
  return(x$length)
}


#' @name lengths
#' @title lengths
#' 
#' @description
#' 
#' establish s3 method for lengths
#'
#' @param gWalk a \code{gWalk} object
#'
#' @return the number of nodes per walk in the gWalk
#' @export
lengths = function(x, use.names = T) {
  UseMethod("lengths")
}

#' @name lengths
#' @title lengths
#' @description
#' A vector of walk lengths associated with this walk
#'
#' @param gWalk a \code{gWalk} object
#'
#' @return the number of nodes per walk in the gWalk
#' @export
`lengths.gWalk` = function(x, use.names = FALSE){
  return(x$lengths)
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
`length.gNode` = function(x){
  return(x$length)
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
`length.gEdge` = function(x){
  return(x$length)
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
`dim.gGraph` = function(x){
  ## input must be a gGraph!
  if (!inherits(x, "gGraph")){
    stop("Error: Invalid input.")
  }
  return(c(x$length, nrow(x$sedgesdt)/2))
}

#' @name refresh
#' @title refresh
#' @description
#' Updates gGraph object to reflect changes in source code
#' 
#' @param gGraph object
#'
#' @return gGraph object
#' @exportMethod refresh
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


  if (is.null(es$sedge.id) || any(is.na(es$sedge.id))) ## infer sedge.id if not provided
  {
    ## arbitrary order for edges
    esign = es[, ifelse(n1+0.1*n1.side < n2+0.1*n2.side, 1, -1)]
    nutag = factor(es[, ifelse(esign>0, paste(n1, n1.side, n2, n2.side), paste(n2, n2.side, n1, n1.side))])
    es[, sedge.id := esign*as.integer(nutag)]
  }

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
#' @name gWalk
#' @title gWalk
#' @description
#' gWalk object
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
                      ## however, for circular walks defined through the sedge.id, we specify the edge between the final and the first  node
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
                        } else if (sum(c(!is.null(snode.id) | !is.null(sedge.id), !is.null(grl))) > 1) {
                          stop("More than one of snode.id, sedge.id and grl cannot be non-NULL")
                        }


                        if (!is.null(grl)) {

                          ## will create the graph first using the unsigned nodes / edges constructor
                          ## doing it differently depending on whether disjoin is flagged

                          grlu = grl.unlist(grl)


                          if (any(ix <- as.character(strand(grlu))=='*'))
                            strand(grlu)[ix] = '+'

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
                        }



                        ## overwrite above if edges also provided
                        if (!is.null(sedge.id)) {
                          walk.names = names(sedge.id)
                          if (is.null(walk.names))
                            walk.names = seq_along(sedge.id)
                          names(sedge.id) = seq_along(sedge.id)
                          private$gWalkFromEdges(sedge.id = sedge.id,
                                                 snode.id = snode.id, ## provide node input to over-ride edges and also check if compatible
                                         graph = graph,
                                         circular = circular,
                                         drop = drop,
                                         meta = meta)
                        }

                        if (nrow(private$pmeta)>0)
                          {
                            setkey(private$pmeta, walk.id)
                            if (nrow(private$pnode))
                              setkey(private$pnode, walk.id)

                            if (nrow(private$pedge))
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
                        if (self$length>0 && nrow(private$pnode))
                        {
                          tmp = private$pnode[, .(walk.id, snode.id)]                   
                          tmp$wid = width(self$graph$nodes$gr)[abs(tmp$snode.id)]
                          if (!is.null(walk.names))
                            private$pmeta$name = walk.names[private$pmeta$walk.id]
                          private$pmeta$wid = tmp[, .(wid = sum(as.numeric(wid), na.rm = TRUE)), keyby = walk.id][.(private$pmeta$walk.id), wid]

                          ## we need to fix factors if drop == TRUE and there are missing walk.ids
                          if (drop == TRUE & any(private$pmeta$walk.id>self$length)){
                            lev = private$pmeta$walk.id
                            private$pmeta[, walk.id := as.integer(factor(as.character(walk.id), lev))]
                            private$pnode[, walk.id := as.integer(factor(as.character(walk.id), lev))]
                            if (nrow(private$pnode))
                              setkey(private$pnode, walk.id)

                            if (nrow(private$pmeta))
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
                        return(invisible(self))
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

                        if (!is.null(dim(i)))
                          stop('subscript dimensions malformed')

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

                        ## create nodes and their "mirrors" allows "walk flipping" using negative indices
                        ui = unique(i) ## need unique indices otherwise dups get created below!!
                        tmp.node =  rbind(                        
                            private$pnode[.(abs(ui)), , allow.cartesian = TRUE],
                            copy(private$pnode[.(abs(ui)), ])[,
                           ":="(walk.id = -walk.id, snode.id = -snode.id)][rev(1:.N), ])

                        new.edge = tmp.edge = data.table()
                        if (nrow(private$pedge)>0)
                          {
                            tmp.edge = rbind(
                              private$pedge[.(abs(i)), , allow.cartesian = TRUE],
                              copy(private$pedge[.(abs(i)), , allow.cartesian = TRUE])[, 
                                                               ":="(walk.id = -walk.id, sedge.id = -sedge.id)][rev(1:.N), ])
                            setkey(tmp.edge, walk.id)
                          }
                        tmp.meta = private$pmeta[.(abs(i)), , allow.cartesian = TRUE]

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

                        if (nrow(private$pnode))
                          setkey(private$pnode, walk.id)

                        if (nrow(private$pmeta))
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
                        if (self$length==0)
                          return(data.table())

                        ix = ix[ix %in% 1:self$length]
                        out = private$pmeta[.(ix), ]
                        if (makelists == FALSE) ## just return raw pmeta subset 
                          return(out)

                        node.sum = data.table(
                          snode.id = rep(list(),self$length),
                          walk.id = private$pmeta$walk.id,
                          key = 'walk.id')

                        if (nrow(private$pnode))
                        {
                          ix2 = ix[ix %in% private$pnode$walk.id]
                          node.sum = private$pnode[.(ix2), .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),][, -1]
                        }                        
                        out = cbind(out, node.sum)
                          
                        if (nrow(private$pedge)>0)
                          out = cbind(out, private$pedge[.(ix), .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),][, -1])
                        return(out)
                      },
                      
                       #' @name fix
                       #' @description
                       #'
                       #' Returns a modified version of the gWalk with seqlevels "fixed", including
                       #' keeping only certain seqlevels, dropping certain seqlevels, and replacing seqlevels.
                       #' 
                       #' Warning: this may modify the walks including getting rid of nodes and edges (i.e. those outside
                       #' the retained seqlevels) and also change coordinates (ie move ranges that were previously on different
                       #' chromosomes to the same chromosome etc.).  Use with caution!
                       #'
                       #' Default behavior is to replace 'chr', with ''.
                       #' 
                       #' @param pattern  character pattern to replace in seqlevels (used in a gsub, can have backreferences)
                       #' @param replacement  character to replace pattern with (used in a gsub, can have backreferences) 
                       #' @param drop  character vector of seqlevels to drop or logical TRUE to drop all seqlevels that are unused  (TRUE)
                       #' @param seqlengths new seqlengths i.e. named integer vector of seqlevels to drop or embed graph into 
                       #' @author Marcin Imielinski
                       #' @return current graph modified in place with additional nodes and edges, as specified by user 
                       fix = function(pattern = NULL, replacement = NULL, drop = TRUE, seqlengths = NULL)
                       {
                         if (is.logical(drop))
                         {
                           if (drop)                             
                             drop = setdiff(seqlevels(self), as.character(seqnames(self$nodes$gr)))
                           else
                             drop = NULL
                         }

                         if (!is.null(drop) | !is.null(seqlengths))
                         {
                           if (!is.null(seqlengths))
                             sl.new = seqlengths
                           else
                             sl.new = seqlengths(self)

                           tmp = as.data.table(self$nodes$dt)
                           ix = tmp$seqnames %in% names(sl.new) & !(tmp$seqnames %in% drop)

                           sl.new = sl.new[intersect(names(sl.new), tmp$seqnames[ix])]

                           gg = refresh(self$graph$copy)
                           gg$nodes$mark(og.node.id = gg$nodes$dt$node.id)
                           gg$fix(seqlengths = sl.new, pattern = pattern, replacement = replacement, drop = drop)

                           if (length(sl.new))
                             {
                               map = unique(gg$dt[, .(snode.id = c(og.node.id, -og.node.id), new.snode.id = c(node.id, -node.id))], by = 'snode.id')
                               pnode.new = merge(private$pnode, map, by = 'snode.id', all.x = TRUE)[, .(snode.id, walk.id, walk.iid, new.snode.id)]
                               setkeyv(pnode.new, c('walk.id', 'walk.iid'))

                               ## walks may need to be split or made empty because there are no longer seqlengths associated with them
                               ## basically each NA to non NA should give a new walk id
                               pnode.new[, new.walk := cumsum(c(0, diff(!is.na(new.snode.id))>0)), by = walk.id]

                               ## broken walks ie those with NA lose the right to be circular if they were before
                               pnode.new[, broken := any(is.na(new.snode.id))]
                               pnode.new[, new.walk.id := paste(walk.id, new.walk)]
                               pnode.new[, new.walk.iid := cumsum(!is.na(new.snode.id)), by = new.walk.id]

                               wmap = unique(pnode.new[, .(walk.id, new.walk.id, broken)])

                               meta.new = merge(private$pmeta, wmap, by = 'walk.id')
                               meta.new[, circular := circular & !broken]
                               setkey(meta.new, new.walk.id)                           
                               tmp = pnode.new[!is.na(new.snode.id), ]
                               snode.id = split(tmp$new.snode.id, factor(tmp$new.walk.id, meta.new$new.walk.id))
                             }
                           else
                           {
                             meta.new = private$pmeta
                             snode.id = lapply(meta.new$walk.id, function(x) c())
                             names(snode.id) = meta.new$walk.id
                           }
                           
                           private$gWalkFromNodes(snode.id = snode.id,
                                                  graph = gg,
                                                  circular = meta.new$circular,
                                                  meta = meta.new
                                                  )
                           
                           if (nrow(private$pnode))
                             setkey(private$pnode, walk.id)
                           
                           if (nrow(private$pedge))
                             setkey(private$pedge, walk.id)
                           
                           if (nrow(private$pmeta))
                             setkey(private$pmeta, walk.id)

                           private$ptimestamp = private$pgraph$timestamp
                         }

                         return(invisible(self))
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
                          cols = setdiff(colnames(private$pmeta), c('snode.id', 'sedge.id'))
                          out = self$dts(ix)[ , cols, with = FALSE]
                          if (nrow(private$pnode))
                          {
                            tmp.node = private$pnode[.(ix), .(walk.id, snode.id)]
                            tmp.node[, nix := private$pgraph$queryLookup(snode.id)$index]
                            setkey(tmp.node, walk.id)
                            
                            teaser = tmp.node[.(ix), ][,
                                              .(gr =
                                                  {
                                                   if (any(is.na(nix)))
                                                     ''
                                                   else
                                                     .tease(gr.string(private$pgraph$gr[nix],
                                                                      mb = FALSE))
                                                  }
                                                ), keyby = walk.id][.(ix),][, -1]                        
                            out = cbind(out, teaser)
                          }

                          print(out)
                          if (self$length>TOOBIG)
                          {
                            message('\n ... \n(', self$length-TOOBIG, ' more walks )')
                          }
                        }
                      },
                      

                      #' @name gWalk json
                      #' @description
                      #' Dumps walk (and underlying graph)
                      #' to json file
                      #' @param filename character path to save to
                      #' @param save whether to save or return list object representing json contents
                      #' @param nfields which node fields to dump to json (NULL)
                      #' @param efields which edge fields to dump to json (NULL)
                      #' @param annotations which graph annotations to dump to json
                      #' @param include.graph if set to TRUE, then the output will include the underlying graph for the gWalk (TRUE)
                      #' @param settings for details see gGraph json() module
                      #' @param no.y for details see gGraph json() module
                      #' @param stack.gap this currently does not do anything, but is intended to control the spread of intervals on the y-axis
                      #' @param verbose
                      #' @author Marcin Imielinski
                      json = function(filename = '.',
                                      save = TRUE,
                                      verbose = FALSE,
                                      annotations = NULL,
                                      nfields = NULL,
                                      efields = NULL,
                                      stack.gap = 1e5, 
                                      include.graph = TRUE,
                                      settings = list(y_axis = list(title = "copy number",
                                                                    visible = TRUE)),
                                      cid.field = NULL,
                                      no.y = FALSE)
                      {
                        if (length(self) == 0){
                            warning('This is an empty gWalk so no JSON will be produced.')
                            return(NA)
                        }
                        if (length(self$edges) == 0){
                            warning('There are no edges in this gWalk so no JSON will be produced.')
                            return(NA)
                        }
                        # check if the gWalk includes walks with no ALT edges and if so then remove these and call json again
                        # non.alt.exist = any(self$dt[,sapply(sedge.id, length) == 0])
                        non.alt.exist = any(self$dt[,sapply(sedge.id, length) == 0] | self$dt[,sapply(sedge.id, anyNA) == TRUE])  ## kat: I had to add this bandaid because I was using subsetted gWalk objects, where walks of length one have NAs in place of NULLs for sedge.id... ultimately this should be patched under the hood in the gWalk subset function so that subsetted walk sedge.id output is EXACTLY the same style as original...

                        if (non.alt.exist){
                          # we call json function again but only including the walks that include at least one ALT edge (anyNA condition added to fix subsetting issue as described in the comments just above)
                          return(refresh(self[self$dt[,sapply(sedge.id, length) > 0] & self$dt[,sapply(sedge.id, anyNA) == FALSE]])$json(filename = filename,
                                                                                         save = save,
                                                                                         verbose = verbose,
                                                                                         annotations = annotations,
                                                                                         nfields = nfields,
                                                                                         efields = efields,
                                                                                         stack.gap = stack.gap,
                                                                                         include.graph = include.graph,
                                                                                         settings = settings,
                                                                                         no.y = no.y))
                        }

                        if (include.graph){
                            ## build graph level js 
                            graph.js = refresh(self$graph)$json(filename = NA, save = FALSE, verbose = verbose,
                                                 annotations = annotations, nfields = nfields, efields = efields,
                                                 settings = settings, no.y = no.y)
                        }


                        pids = split(self$dt[, .(pid = walk.id,
                                           strand = '+',
                                           type = ifelse(self$circular, 'cycle', 'path'))], 1:self$length)

                        efields = unique(c('type', efields))
                        protected_efields = c('cid', 'source', 'sink', 'title', 'weight')
                        rejected_efields = intersect(efields, protected_efields)
                        if (length(rejected_efields) > 0){
                            warning(sprintf('The following fields were included in efields: "%s", but since these are conserved fields in the json walks output then they will be not be included in efields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.', paste(rejected_efields, collapse = '" ,"')))
                            efields = setdiff(efields, rejected_efields)
                        }
                        missing_efields = setdiff(efields, names(self$edges$dt))
                        if (length(missing_efields) > 0){
                            warning(sprintf('Invalid efields value/s provided: "%s". These fields were not found in the gWalk and since will be ignored.', paste(missing_efields, collapse = '" ,"')))
                            efields = intersect(efields, names(self$edges$dt))
                        }

                        sedu = dunlist(self$sedge.id)
                        cids = lapply(unname(split(cbind(data.table(cid = sedu$V1,
                                                       source = self$graph$edges[sedu$V1]$left$dt$snode.id,
                                                       sink = -self$graph$edges[sedu$V1]$right$dt$snode.id, # notice that we need to add negative sign here to meet the gGnome.js expectations
                                                       title = "",
                                                       weight = 1),
                                                   self$graph$edges[sedu$V1]$dt[, ..efields]
                                                 ), sedu$listid)),
                                      function(x) unname(split(x, 1:nrow(x))))

                        snu = dunlist(self$snode.id)
                        snu$ys = gGnome:::draw.paths.y(self$grl) %>% unlist

                        protected_nfields = c('chromosome', 'startPoint', 'endPoint',
                                              'y', 'type', 'strand', 'title')
                        rejected_nfields = intersect(nfields, protected_nfields)
                        if (length(rejected_nfields) > 0){
                            warning(sprintf('The following fields were included in nfields: "%s", but since these are conserved fields in the json walks output then they will be not be included in nfields. If these fields contain important metadata that you want included in the json output, then consider renaming these field names in your gWalk object.', paste(rejected_nfields, collapse = '" ,"')))
                            nfields = setdiff(nfields, rejected_nfields)
                        }
                        missing_nfields = setdiff(nfields, names(self$nodes$dt))
                        if (length(missing_nfields) > 0){
                            warning(sprintf('Invalid nfields value/s provided: "%s". These fields were not found in the gWalk and since will be ignored.', paste(missing_nfields, collapse = '" ,"')))
                            nfields = intersect(nfields, names(self$edges$dt))
                        }
                        iids = lapply(unname(split(cbind(
                          data.table(iid = abs(snu$V1)),
                          self$graph$nodes[snu$V1]$dt[,
                                                      .(chromosome = seqnames,
                                                        startPoint = start,
                                                        endPoint = end,
                                                        y = snu$ys,
                                                        type = "interval",
                                                        strand = ifelse(snu$V1 > 0, "+", "-"),
                                                        title = abs(snu$V1))],
                          self$graph$nodes[snu$V1]$dt[,..nfields]), snu$listid)),
                          function(x) unname(split(x, 1:nrow(x))))

                        walks.js = lapply(1:length(self), function(x)
                          c(as.list(pids[[x]]),
                            list(cids = rbindlist(cids[[x]])),
                            list(iids = rbindlist(iids[[x]]))))

                        if (include.graph){
                            out = c(graph.js, list(walks = walks.js))
                        } else {
                            out = list(walks = walks.js)
                        }

                        if (save){
                          if (verbose)
                          {
                            message("Saving JSON to: ", filename)
                          }
                          jsonlite::write_json(out, filename,
                                               pretty=TRUE, auto_unbox=TRUE, digits=4)
                          return(normalizePath(filename))
                        } else {
                          return(out)
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
                            pdt = private$pgraph$dt
                            pdt = pdt[, setdiff(names(pdt), c('walk.id', 'walk.iid')), with = FALSE]
                       
                            tmpdt = merge(private$pnode, pdt, by = 'snode.id')
                            setkeyv(tmpdt, c("walk.id", "walk.iid")) #fix order to match private$pnode
                            tmpdt = tmpdt[.(private$pnode$walk.id, private$pnode$walk.iid), ]
                            out = eval(parse(text = paste("tmpdt[, ", lazyeval::expr_text(node), ", keyby = walk.id]")))
                            out = unique(out, by = "walk.id")
                            out[.(1:self$length), ][[2]]
                          }, error = function(e) NULL)
                        }

                        if (is.null(out) & !missing('edge'))
                        {

                          out = tryCatch({
                            ped = private$pgraph$sedgesdt
                            ped = ped[, setdiff(names(ped), c('walk.id', 'walk.iid')),
                                      with = FALSE]
                            tmpdt = merge(private$pedge, ped, by = 'sedge.id')
                            setkeyv(tmpdt, c("walk.id", "walk.iid")) #fix order to match private$pnode
                            tmpdt = tmpdt[.(private$pedge$walk.id, private$pedge$walk.iid), ]
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
                      disjoin = function(gr = NULL, graph = NULL, by = NULL, collapse = TRUE, na.rm = TRUE, avg = FALSE, sep = ',', FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg))
                      {
                        self$check
                        if (self$length==0){
                          return(invisible(self))
                        }

                        node.id = private$pgraph$nodes$dt$node.id

                                                
                        ## disjoin current graph
                        tmpg = private$pgraph$copy

                        ## mark old nodes in current graph
                        if (!collapse) ## keep track of node ids
                          {
                            tmpg$nodes$mark(node.id.old = 1:length(tmpg))
                          }

                        old.gr = tmpg$nodes$gr

                        if (!is.null(graph))
                        {
                          if (!collapse) ## unmark node.id.old to prevent collisions
                            graph$nodes$mark(node.id.old = NA_integer_)

                          tmpg = c(tmpg, graph)
                        }

                        if (!is.null(gr))
                        {
                          gr$node.id.old = NULL
                          ## grb = gG(breaks = gr)
                          ## if (!collapse)
                          ##   grb$nodes$mark(node.id.old = NA_integer_)
                          ## tmpg = c(tmpg, grb)
                        }

                        tmpg$disjoin(gr = gr, by = by, na.rm = na.rm, collapse = collapse, 
                                               avg = avg, FUN = FUN)

                        ## we will redo overlaps here
                        ## since we know tmpg$nodes$gr are disjoint
                        ## and we need coordinates to properly reorder
                        ## any disjoint walk nodes that came from a parent
                        ## node

                        old.gr$snode.id.old = old.gr$snode.id
                        old.gr$node.id.old = old.gr$node.id
                        if (collapse) ## we match on intervals since many input nodes collapsed into a single node
                          map = gr2dt(sort(gr.findoverlaps(old.gr, tmpg$nodes$gr, qcol = 'snode.id.old', scol = 'snode.id', by = by)))
                        else
                        {
                          old.gr$query.id = 1:length(old.gr)
                          new.gr = tmpg$nodes$gr %>% gr2dt
                          new.gr$subject.id = 1:nrow(new.gr)

                          map = merge(old.gr[, c("query.id", "snode.id.old", "node.id.old")]  %>% gr2dt,
                                      new.gr[, .(subject.id, node.id.old, snode.id)],
                                      by = 'node.id.old')[, .(query.id, subject.id, snode.id.old, snode.id)]
                        }

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

                        if (nrow(private$pnode))
                          setkey(private$pnode, walk.id)

                        if (nrow(private$pedge))
                          setkey(private$pedge, walk.id)

                        if (nrow(private$pmeta))
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

                        ## we will first simplify tmpg
                        ## then rematch the walk nodes to tmpg
                        newgraph = private$pgraph$copy
                        newgraph$nodes$mark(parent.node.id = 1:length(newgraph))
                        newgraph$simplify(FUN = NULL, by = by)
                        
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

                        if (nrow(private$pnode))
                          setkey(private$pnode, walk.id)

                        if (nrow(private$pedge))
                          setkey(private$pedge, walk.id)

                        if (nrow(private$pmeta))
                          setkey(private$pmeta, walk.id)

                        private$ptimestamp = private$pgraph$timestamp
                        return(invisible(self))
                      },


                      #' @name mark
                      #' @description
                      #'
                      #' marks the nodes and edges in the graph corresponding to this walk
                      #' @author Marcin Imielinski
                      mark = function(...)
                      {
                        self$check
                        self$nodes$mark(...)
                        self$edges$mark(...)
                      },

                      #' @name fitcn
                      #' @description
                      #'
                      #' Given gWalks on a graph with node and edge
                      #' metadata "cn" attempts to find a non-negative integer
                      #' combination of those walks that generates those walks
                      #' minimizes the number of unique walks used (L0 norm).
                      #' 
                      #' Modifications to this objecti include (default behavior)
                      #' if min.alt = TRUE then solution will minimizes the number of
                      #' walks that contain an ALT edge.  If weight is provided
                      #' then the objective function be \eqn{\sum_i w_i*I[x_i>0]} where
                      #' w_i is the provided weight for walk i and I[x_i>0] is indicator
                      #' that is 1 if walk i has a nonzero copy number in the solution
                      #' and 0 otherwise. 
                      #'
                      #' Complains if "cn" is not a metadata of nodes and edges
                      #' and if the graph is not "balanced" ie either
                      #' every node cn >= sum(edge cn) for both incoming
                      #' and outgoing edges.
                      #'
                      #' @param logical flag trim whether to "trim" the graph relative to the footprint of the walks
                      #' @param weight either a numeric metadata field of self or a  self$length integer weight
                      #' @param min.alt logical flag whether to make objective function
                      #' @param verbose logical flag whether to output CPLEX output
                      #' @param edgeonly logical
                      #' @param evolve logical
                      #' @param n.sol numeric
                      #' @author Julie Behr
                      fitcn = function(trim=TRUE,                                       
                                       weight=NULL,
                                       obs.mat=NULL,
                                       verbose = FALSE,
                                       min.alt=TRUE,
                                       edgeonly = FALSE, 
                                       evolve=FALSE,
                                       n.sol=2)
                      {
                          gw = self
                          gg = self$graph

                          if (is.null(gw$graph$nodes$dt$cn) |
                              is.null(gw$graph$edges$dt$cn))
                          {
                              stop('cn field is missing from node and edge metadata')
                          }

                          ## quick check to see if graph is junction balanced
                          ## ie that node cn >= sum(junction cn) for each node side
                          ## and node cn == sum(junction cn) on all node sides
                          ## where there is no loose end
                          rcn = gg$nodes$eval(sum(cn));
                          lcn = gg$nodes$eval(sum(cn), FALSE);

                          jbal.left = all((lcn == gg$nodes$dt$cn &
                                           !gg$nodes$dt$loose.left) |
                                          lcn <= gg$nodes$dt$cn &
                                          gg$nodes$dt$loose.left, na.rm = TRUE)

                          jbal.right = all((rcn == gg$nodes$dt$cn &
                                            !gg$nodes$dt$loose.right) |
                                           rcn <= gg$nodes$dt$cn &
                                           gg$nodes$dt$loose.right, na.rm = TRUE)

                          if (!jbal.left | !jbal.right)
                              warning('graph does not appear to be junction balanced, please check inputs')


                          constrain.evolution = function(K, gw, A, b, sense){
                              h = K[gw$edges[type=="ALT"]$dt[!duplicated(edge.id), edge.id],]
                              A = rbind(A, cbind(sparseMatrix(1, 1, x=0, dims=dim(h)), h))
                              b = c(b, rep(1, nrow(h)))
                              sense = c(sense, rep("L", nrow(h)))
                              return(list(A=A, b=b, sense=sense))
                          }

                          ## only affects which x will optimize min(x * c),
                          ## currently doesn't matter because all solutions
                          ## (regardless of how optimal) are returned
                          constrain.observations = function(obs.mat, A, b, cvec, sense, vtype){
                              if(!(ncol(obs.mat)*2)==ncol(A))
                                  stop('input obs.mat contains the wrong number of columns; should match length of gw')
                              p = nrow(obs.mat)
                              w = ncol(obs.mat)
                              Zero = sparseMatrix(1, 1, x=0, dims=c(2*w*p, 2*w*p))
                              A0 = Zero[rep(1, nrow(A)), 1:(2*p)]
                              Ap = cbind(Zero[rep(1, p), 1:w], sign(obs.mat), diag(rep(-1, p)), Zero[rep(1,p), 1:p])
                              Mpub = cbind(Zero[rep(1,p), 1:(2*w)], diag(rep(1, p)), diag(rep(-1e7, p)))
                              Mplb = cbind(Zero[rep(1,p), 1:(2*w)], diag(rep(1, p)), diag(rep(-0.1, p)))
                              Amp = rbind(cbind(A, A0), Ap, Mpub, Mplb)
                              b = c(b, rep(0, 3*p))
                              cvec = c(cvec, rep(0, p), -1*rowMax(obs.mat))
                              sense = c(sense, rep("E", p), rep("L", p), rep("G", p))
                              vtype = c(vtype, rep("I", p), rep("B", p))
                              
                              return(list(A=Amp, b=b, c=cvec, sense=sense, vtype=vtype))
                          }
                          

                          ## walks x edges matrix
                          generate.Ke = function(gw){
                              dt = gw$edgesdt[, c("walk.id", "sedge.id")][, edge.id := abs(sedge.id)]
                              dt$listid = factor(dt$walk.id, 1:length(gw))
                              cdt = dcast(dt[!is.na(sedge.id),], listid ~ edge.id, fun.aggregate=length, value.var="edge.id", drop=FALSE)
                              mat = cdt[, -1]
                              rownames(mat) = cdt$listid
                              return(t(mat))
                          }

                          ## walks x nodes matrix
                          generate.Kn = function(gw){
                              dt = gw$nodesdt[, c("walk.id", "snode.id")][, node.id := abs(snode.id)]
                              dt$listid = factor(dt$walk.id, 1:length(gw))
                              cdt = dcast(dt[!is.na(snode.id),], listid ~ node.id, fun.aggregate=length, value.var="node.id", drop=FALSE)
                              mat = cdt[, -1]
                              rownames(mat) = cdt$listid
                              return(t(mat))
                          }
                          
                          generate.Amat = function(K){
                              M = 1e7
                              K = as(K, 'sparseMatrix')
                              w = ncol(K)

                              ## upper bound is infinity if indicator is positive
                              Zero = sparseMatrix(1, 1, x=0, dims=c(2*w, 2*w))
                              Amub = cbind(diag(rep(1, w)), diag(rep(-M, w)))

                              ## lower bound > 0 if indicator is positive
                              Amlb = cbind(diag(rep(1, w)), diag(rep(-0.1, w)))

                              A = rbind(cbind2(K, Zero[rep(1, nrow(K)), (w+1:w)]), Amub, Amlb)
                              return(A)
                          }

                          generate.bvec = function(e, K){
                              w = ncol(K)
                              bvec = c(e, rep(0, 2*w))
                              return(bvec)
                          }

                          generate.cvec = function(K, weight, min.alt, gw){
                              if(!is.null(weight) | min.alt) weight = prep.weight(weight, min.alt, gw)
                              w = ncol(K)
                              if(is.null(weight)) weight=rep(1, w)
                              cvec = c(rep(0, w), weight)
                              return(cvec)
                          }

                          prep.weight = function(weight, min.alt, gw){
                              if(!is.null(weight)){
                                  if(length(weight)==1 && is.character(weight) && weight %in% colnames(gw$dt)) weight = gw$dt[, weight, with=F]
                                  if(!(is.numeric(weight))){
                                      stop("weight must either be numeric vector of same length as gw or the name of a single numeric annotation in gw")
                                  }
                              }
                              if(min.alt){
                                  if(!is.null(weight)){
                                      warning("modifying input weight to satisfy min.alt=TRUE")
                                  } else weight = rep(1, ncol(K))
                                  numalt = gw$eval(edge = sum(type=='ALT'))
                                  ## walks with no edges or only REF edges will not add to objective function
                                  weight[is.na(numalt) | numalt==0] = 0
                              }
                              return(weight)
                          }
                          
                          generate.vtype = function(K){
                              w = ncol(K)
                              vtype = c(rep("I", w), rep("B", w))
                              return(vtype)
                          }

                          generate.sense = function(K){
                              w = ncol(K)
                              e = nrow(K)
                              sense = c(rep("E", e), rep("L", w), rep("G", w))
                              return(sense)
                          }

                          if(!"cn" %in% colnames(gw$graph$edges$dt)){
                              stop('cn field must be populated in the input graph node and edges metadata')
                                        #                          warning("gw$graph does not have cn field on edges. Defaulting to all cn=1...")
                              if(trim){
                                e = rep(1, length(unique(abs(unlist(gw$sedge.id)))))
                                e2 = rep(1, length(unique(abs(unlist(gw$snode.id)))))
                              } else{
                                e = rep(1, length(gw$graph$edges))
                                e2 = rep(1, length(gw$graph$nodes))
                              }
                          } else{
                              if(trim){
                                e = gw$graph$edges[sedge.id %in% abs(unlist(gw$sedge.id))]$dt$cn
                                e2 = gw$graph$nodes[snode.id %in% abs(unlist(gw$snode.id))]$dt$cn
                              } else{
                                e = gw$graph$edges$dt$cn
                                e2 = gw$graph$nodes$dt$cn
                              }
                          }

                          if (edgeonly)
                            {
                              K = generate.Ke(gw)
                            }
                          else
                            {
                              K = rbind(generate.Ke(gw), generate.Kn(gw))
                              e = c(e, e2)
                            }

                          if(nrow(K) != length(e)){
                              stop("Mismatch between size of A matrix and length of b vector. Some edges in gw$graph are not covered by gw. If this was intended, try trim=TRUE")
                          }
                          A = generate.Amat(unname(K))
                          b = generate.bvec(e, K)
                          c = generate.cvec(K, weight, min.alt, gw)
                          sense = generate.sense(K)
                          vtype = generate.vtype(K)
                          if(evolve){
                              ll = constrain.evolution(K, gw, A, b, sense)
                              A = ll$A
                              b = ll$b
                              sense = ll$sense
                          }
                          if(!is.null(obs.mat)){
                              ll = constrain.observations(obs.mat, A, b, c, sense, vtype)
                              A = ll$A
                              b = ll$b
                              c = ll$c
                              sense = ll$sense
                              vtype = ll$vtype
                          }
                          ## sol = Rcplex::Rcplex(cvec = c,
                          sol = Rcplex2(cvec = c,
                                        Amat = A,
                                        bvec = b,
                                        sense = sense,
                                        Qmat = NULL,
                                        lb = 0,
                                        ub = Inf,
                                        n = n.sol,
                                        objsense = "min",
                                        vtype = vtype,
                                        control = list(
                                            trace = ifelse(verbose>=1, 1, 0),
                                            tilim = 100,
                                            epgap = 1))

                          if (!is.null(sol$xopt)){
                              sol = list(sol)
                          }
                          if(length(sol)==0){
                              stop("No solutions found satisfying given constraints")
                          }

                          rerun=T
                          while(rerun){
                              z = sign(vtype == 'B')
                              P = do.call(rbind, lapply(sol, function(x) x$xopt*z))
                              p = rowSums(P)-1

                              Ahat = rbind(A, P)
                              bhat = c(b, p)
                              sensehat = c(sense, rep('L', length(p)))

                              ## sol.new = Rcplex::Rcplex(cvec = c,
                              sol.new = Rcplex2(cvec = c,
                                                Amat = Ahat,
                                                bvec = bhat,
                                                sense = sensehat,
                                                Qmat = NULL,
                                                lb = 0,
                                                ub = Inf,
                                                n = n.sol,
                                                objsense = "min",
                                                vtype = vtype,
                                                control = list(
                                                    trace = ifelse(verbose>=1, 1, 0),
                                                    tilim = 100,
                                                    epgap = 1))
                              ## this could be optional?
                              ##    sol.new = sol.new[round(sapply(sol.new, function(x) x$obj))==round(sol[[1]]$obj)]
                              if(length(sol.new)==0){
                                  rerun=F
                              } else{
                                  sol = c(sol, sol.new)
                                  if(length(sol)>=n.sol){
                                      sol = sol[1:n.sol]
                                      rerun=F
                                  }
                              }
                          }

                          gw$set(cn = round(sol[[1]]$xopt[1:ncol(K)]))
                          if(length(sol)>1){
                              all = setNames(lapply(sol,
                                                    function(s)
                                                        round(s$xopt[1:ncol(K)])),
                                             paste("cn", 1:length(sol), sep="."))
                              do.call(gw$set, all)
                          }
                          return(invisible(gw))
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


                        tmp.grl = unname(self$grl)
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
                        pnode.empty = as.data.table(list(walk.id = c(), snode.id = c()))
                        pedge.empty = as.data.table(list(walk.id = c(), sedge.id = c()))

                        if (length(snode.id)==0)
                        {
                          private$pnode = pnode.empty
                          private$pedge = pedge.empty
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

                        private$pmeta = data.table(walk.id = seq_along(snode.id), name = names(snode.id), length = elementNROWS(snode.id), wid = 0, 
                                                   circular = circular)
                        
                        if (!is.null(meta))
                        {                       
                          if (nrow(meta) != length(snode.id))
                          {
                            stop('data.table of meta.data must be same length and order as node, edge, or GRangesList list input to gWalk constructor')
                          }
                          
                          GW.NONO = c('walk.id', 'walk.iid', 'circular', 'snode.id', 'sedge.id', 'wid', 'length', 'name', 'subject.id', 'query.id', 'i')
                          good.cols = setdiff(names(meta), GW.NONO)
                          
                          if (length(good.cols)>0)
                          {                            
                            private$pmeta = cbind(private$pmeta, meta[, good.cols, with = FALSE])
                          }
                        }

                        ## always faster to unlist than lapply
                        tmp = dunlist(unname(snode.id))

                        if (nrow(tmp)==0)
                        {
                          private$pnode = pnode.empty
                          private$pedge = pedge.empty
                          private$pgraph = graph
                          return()
                        }

                        pnode = tmp[, .(walk.id = listid, snode.id = V1)]
                        
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
                          setkey(pedge, walk.id)
                          pedge$sedge.id = NA
                          if (nrow(sedgesdt))
                            pedge$sedge.id = sedgesdt[.(pedge$from, pedge$to), sedge.id]

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

                        tmp = copy(pnode)
                        tmp$wid = width(self$graph$nodes$gr)[abs(tmp$snode.id)]
                        private$pmeta$wid = tmp[, .(wid = sum(as.numeric(wid), na.rm = TRUE)), keyby = walk.id][.(private$pmeta$walk.id), wid]
                        if (nrow(private$pmeta) && any(ix <- is.na(private$pmeta$wid)))
                          private$pmeta$wid[ix] = 0


                      },


                      gWalkFromEdges = function(sedge.id,
                                                snode.id = NULL,
                                              graph = NULL,
                                              meta = NULL, ## metadata with one row per walk, can include every metadata field except for $circular, $walk.id
                                              drop = FALSE, 
                                              circular = NULL)  ## logical vector specifying which contigs are circular, i.e. have an implied edge from the final node to the first                                            
                      {
                        pnode.empty = as.data.table(list(walk.id = c(), walk.iid = c(), snode.id = c()))
                        pedge.empty = as.data.table(list(walk.id = c(), walk.iid = c(), sedge.id = c()))

                        if (length(sedge.id)==0)
                        {
                          private$pnode = pnode.empty
                          private$pedge = pedge.empty
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
                        
                        
                        private$pmeta = data.table(walk.id = seq_along(sedge.id), name = names(sedge.id), length = elementNROWS(sedge.id), wid = 0,
                                                   circular = circular)

                        if (!is.null(meta)) ## need this for some reason to get rid of pass by reference stickiness (doesn't work higher in call stack)
                        {
                          if (!is.data.table(meta))
                            meta = as.data.table(meta)
                          meta = copy(meta)
                          
                          if (nrow(meta) != length(sedge.id))
                          {
                            stop('data.table of meta data must be same length and order as node, edge, or GRangesList list input')
                          }

                          GW.NONO = c('walk.id', 'walk.iid', 'circular', 'snode.id', 'sedge.id', 'wid', 'length', 'name', 'query.id', 'subject.id', 'i')
                          good.cols = setdiff(names(meta), GW.NONO)
                          
                          if (length(good.cols)>0)
                          {
                            private$pmeta = cbind(private$pmeta,
                                                  meta[, good.cols, with = FALSE])
                          }
                        }
                        
                        private$pgraph = graph
                        tmp = dunlist(unname(sedge.id))
                        if (nrow(tmp)==0)
                        {
                          private$pnode = pnode.empty
                          private$pedge = pedge.empty
                          return()
                        }
                        setkey(private$pmeta, walk.id)

                        pedge = tmp[, .(walk.id = listid, sedge.id = V1)]

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
                        
                        pnode = merge(pedge, private$pmeta, by = 'walk.id')[, .(snode.id = {if (circular) from else c(from, to[.N])}), by = .(walk.id, circular)][, .(walk.id, snode.id)]

                        if (is.null(snode.id))
                        {
                          snode.id = split(pnode$snode.id, pnode$walk.id)
                        }
                        else ## check to see if provided node ids are compatible for walks with 1 or more edges
                        {
                          .f = function(y) sapply(y, function(x) paste(x, collapse = ', '))
                          ix = names(sedge.id)[base::lengths(sedge.id)>0]
                          if (length(ix))
                          {
                            if (!identical(.f(split(pnode$snode.id, pnode$walk.id)[ix]), .f(snode.id[ix])))
                              stop('plesae check inputs - provided snode.id and sedge.ids are not compatible')
                          }
                          tmp = dunlist(unname(snode.id))
                          ## if so let's populate pnode with the snode.id version
                          pnode = tmp[, .(walk.id = listid, snode.id = V1)]
                        }                          
                        
                        if (nrow(pedge)>0)
                          pedge[, walk.iid := 1:.N, by = walk.id]
                       
                        if (nrow(pnode)>0)
                          pnode[, walk.iid := 1:.N, by = walk.id]

                        setkey(pnode, walk.id)
                        setkey(pedge, walk.id)

                        ## check to see if circular are truly circular
                        ## i.e. does the final edge touch the first node
                        if (any(private$pmeta$circular))
                        {
                          tmp = merge(pedge, private$pmeta[, .(walk.id, circular)], by = 'walk.id') %>%
                            merge(pnode[walk.iid == 1, .(walk.id, snode.id)], by = 'walk.id')                           
                          tmp = tmp[circular == TRUE, .(check = to[.N] == snode.id), by = .(walk.id, snode.id)]

                          if (length(notcircle <- tmp[check == FALSE, walk.id]))
                            stop(paste('Walk', paste(notcircle, collapse = ', '), "not circular (last edge doesn't touch the first node), please check circular flag"))
                        }
                       
                        private$pedge = pedge                        
                        private$pnode = pnode
                        
                        ## run everything through gWalkFromNodes to make sure compatible
                        ## private$gWalkFromNodes(snode.id = snode.id,
                        ##                        graph = graph,
                        ##                        circular = circular,
                        ##                        drop = drop, 
                        ##                        private$pmeta
                        ##                        )
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
                        out = split(private$pgraph$gr[nix], factor(private$pnode$walk.id, private$pmeta$walk.id))[as.character(1:self$length)]
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

                        node.sum = data.table(
                          snode.id = rep(list(),self$length),
                          walk.id = private$pmeta$walk.id,
                          key = 'walk.id')
                        
                        if (nrow(private$pmeta))
                        {
                          ix = private$pmeta$walk.id
                          node.sum = private$pnode[.(ix), .(snode.id = list(c(snode.id))), keyby = walk.id][.(ix),][, -1]
                        }                        

                        return(node.sum$snode.id)
                      },

                      sedge.id = function()
                      {
                        self$check
                        if (self$length==0)
                          return(list())

                        edge.sum = data.table(
                          sedge.id = rep(list(),self$length),
                          walk.id = private$pmeta$walk.id,
                          key = 'walk.id')

                        if (nrow(private$pmeta))
                        {
                          ix = unique(private$pmeta$walk.id)
                          edge.sum = private$pedge[.(ix), .(sedge.id = list(c(sedge.id))), keyby = walk.id][.(ix),][, -1]
                        }                        

                        return(lapply(edge.sum$sedge.id, function(x) x[!is.na(x)]))
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

                      lengths = function() {
                        self$check
                        return(lengths(self$grl))
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

                      #' @name nalt
                      #' @description
                      #' computes number of alt edges in each walk
                      nalt = function()
                      {
                        nalt = pmax(0, self$eval(edge = sum(class != 'REF')), na.rm = TRUE)
                                        #                         private$pmeta[, nalt := nalt]
                        return(nalt)
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
  ##  snode.id = do.call('c', lapply(gWalk.list, function(x) {y = x$snode.id; names(y) = x$dt$name; y}))
  sedge.id = do.call('c', lapply(gWalk.list, function(x) {y = x$sedge.id; names(y) = x$dt$name; y}))
  snode.id = do.call('c', lapply(gWalk.list, function(x) {y = x$snode.id; names(y) = x$dt$name; y}))
  circular = do.call('c', lapply(gWalk.list, function(x) x$circular))
  metas = rbindlist(lapply(gWalk.list, function(x) x$meta), fill = TRUE)

  return(gWalk$new(snode.id = snode.id, sedge.id = sedge.id, circular = circular, graph = gWalk.list[[1]]$graph, meta = metas))
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
#' walks.grl = readRDS(system.file('extdata', 'gw.grl.rds', package = "gGnome"))
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
#' nid.list = list(c(-888, 322, 324, 325, 325), c(-708, 699, -706, 701, 702, 702))
#' gW(snode.id = nid.list, graph = gg)
#'
#' ## from lists of signed edge ids
#' eid.list = list(c(-1462, 1461), c(-133, -132, -134))
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

#' @name gdist
#' @title gdist
#' @description
#'
#' Compute "distance" between two gGraphs using node and edge features.
#'
#' Distance computation involves
#' (1) assigning a value to each node and edge via value / nvalue / evalue expression (default cn metadata field)
#' (2) assigning a weight to each disjoint interval and edge in the union of the two graphs (log10(width) to intervals and 1 to edges)
#' (3) for each disjoint interval computing a ndist = sum(weight*|value1-value2|)/sum(weight), where if the interval does not exist for one of the graphs we use the value "fill" in place
#' (4) for each (ALT) edge present in one or both graphs compute edist = sum(weight*|value1-value2|)/sum(weight), filling in missing edges as above
#' (5) combining node.coef*ndist + edge.coef*edist
#'
#' Variations on these distances include defining alternate weight and value expressions
#' for nodes and edges (based on their metadata), changing the padding to align ALT junctions that are not exactly identical.
#'
#' If graphs are non-haploid (ie have two more intervals for a given reference region)
#' then agg.fun (default min) is used to aggregate the several ndist that map to a
#' that given reference region into a single ndist value for that coordinate. 
#' 
#'
#' @param gg1 first graph (or list of graphs)
#' @param gg2 second graph (or list of graphs)
#' @param value expression on node and edge metadata compute node and edge value
#' @param weight expression on node and edge metadata to compute node and edge weights
#' @param nvalue expression on node metadata to compute node value
#' @param evalue expression on edge metadata to compute edge value
#' @param nweight expression on node metadata to compute node weight
#' @param eweight expression on edge metadata to compute edge weight
#' @param normalize logical flag specifying whether to divide distance by sum(weight) (if TRUE) or not
#' @param fill 0 value to fill value for missing nodes and edges
#' @param na.rm  logical flag specifying whether to remove NA when summing distance (should not be any NA's if value and weight are well specified
gdist = function(gg1, gg2,
                 weight = quote(log10(width)),
                 value = quote(cn),
                 nvalue = NULL, 
                 evalue = NULL,
                 nweight = NULL,
                 eweight = 1,
                 fill = 0,
                 normalize = TRUE,
                 na.rm = TRUE,
                 pad = 1000,
                 agg.fun = min,
                 node.coef = 1,
                 edge.coef = 1,
                 verbose = FALSE,
                 L = 1
                 )
{
  if (L==1)
    lfun = function(x, y) abs(x-y)
  else if (L==0)
    lfun = function(x, y) abs(sign(x-y))
  else if (L==2)
    lfun = function(x, y) (x-y)^2
  else
    stop('L specifies L norm and can be either 0, 1, ,2')
  
  if (deparse(substitute(nvalue)) == "NULL")
    nvalue = substitute(value)

  if (deparse(substitute(evalue)) == "NULL")
    evalue = substitute(value)

  if (deparse(substitute(nweight)) == "NULL")
    nweight = substitute(weight)

  if (deparse(substitute(eweight)) == "NULL")
    eweight = substitute(weight)

  .broken = function(x) is.null(x) || (!is.numeric(x) & !is.integer(x) & !is.logical(x))

  ## compare intervals

  ## make disjoint nodes
  grb = grbind(gg1$nodes$gr, gg2$nodes$gr)
  if (length(grb)==0)
    return(NA)

  dnodes = gr.stripstrand(disjoin(grb))
  values(dnodes) = values(grb)[gr.match(dnodes, grb), ]

  ## weight is a feature of the disjoint granges, if ambiguity will choose the first
  nwt = tryCatch(
    eval(eval(parse(text = substitute(deparse(nweight))),
              parent.frame()), gr2dt(dnodes), parent.frame(2)), error = function(e) NULL)

  if (.broken(nwt)) ## try another way
     nwt = tryCatch(
       eval(eval(parse(text = substitute(deparse(substitute(nweight)))),
                 parent.frame()), gr2dt(dnodes), parent.frame(2)), error = function(e) NULL)
  
   if (.broken(nwt)) ## give up 
  {
    warning('weight expression did not complete for nodes, defaulting to 1')
    nwt = rep(1, length(dnodes))    
  }

  ## disjoin by dnodes
  suppressWarnings(gg1$disjoin(dnodes))
  suppressWarnings(gg2$disjoin(dnodes))

  gr1 = gg1$nodes$gr
  gr2 = gg2$nodes$gr

  ## painful acrobatics to allow expressions as inputs to value
  nval1 = tryCatch(
    eval(eval(parse(text = substitute(deparse(nvalue))),
              parent.frame()), gr2dt(gr1), parent.frame(2)), error = function(e) NULL)

  if (.broken(nval1)) ## try another way
  {
    nval1 = tryCatch(
      eval(eval(parse(text = substitute(deparse(substitute(nvalue)))),
                parent.frame()), gr2dt(gr1), parent.frame(2)), error = function(e) NULL)
  }
  
  nval2 = tryCatch(
    eval(eval(parse(text = substitute(deparse(nvalue))),
              parent.frame()), gr2dt(gr2), parent.frame(2)), error = function(e) NULL)

  if (.broken(nval2)) ## try another way
    nval2 = tryCatch(
      eval(eval(parse(text = substitute(deparse(substitute(nvalue)))),
                parent.frame()), gr2dt(gr2), parent.frame(2)), error = function(e) NULL)


  if (.broken(nval1)) ## give up 
  {
    warning(sprintf('value expression failed for nodes in graph 1, defaulting to 1'))
    nval1 = rep(1, length(gr2))
  }

  if (.broken(nval2)) ## give up 
  {
    warning(sprintf('value expression failed for nodes in graph 2, defaulting to 1'))
    nval2 = rep(1, length(gr2))
  }

  gr1$value = nval1
  gr2$value = nval2

  ## intersect and setdiff node pairs in gg1 and gg2
  nodes.int = gr1[, "value"] %*% gr2[, "value"]
  names(values(nodes.int))[ncol(values(nodes.int))] = 'value.1'
  nodes.sd1 = gr1[!(gr1 %^% gr2)][, "value"]
  nodes.sd2 = gr2[!(gr2 %^% gr1)][, "value"]
  names(values(nodes.sd2)) = 'value.1'

  ## outer + inner join of nodes based on coordinate
  mnodes = grbind(nodes.int, nodes.sd1, nodes.sd2)[, c("value", "value.1")]
  names(values(mnodes)) = c('value', 'value.1')
  mnodes$ix = match(mnodes, dnodes)
  mnodes = gr2dt(mnodes)

  ## fill any NAs for field
  mnodes[is.na(value), value := fill]
  mnodes[is.na(value.1), value.1 := fill]

  ## node distance
  if (normalize)
    ndist = mnodes[, .(diff = agg.fun(nwt[ix]*lfun(value, value.1))), by = ix][, sum(diff, na.rm = na.rm)/sum(nwt[ix], na.rm = na.rm)]
  else
    ndist = mnodes[, .(diff = agg.fun(nwt[ix]*lfun(value, value.1))), by = ix][, sum(diff, na.rm = na.rm)]

  jj1 = gg1$junctions[type == 'ALT']
  jj2 = gg2$junctions[type == 'ALT']

  eval1 = tryCatch(
    eval(eval(parse(text = substitute(deparse(nvalue))),
              parent.frame()), gr2dt(gr1), parent.frame(2)), error = function(e) NULL)

  if (.broken(eval1)) ## try another way 
    eval1 = tryCatch(
      eval(eval(parse(text = substitute(deparse(substitute(nvalue)))),
                parent.frame()), gr2dt(gr1), parent.frame(2)), error = function(e) NULL)
  
  
  eval2 = tryCatch(
      eval(eval(parse(text = substitute(deparse(nvalue))),
              parent.frame()), gr2dt(gr2), parent.frame(2)), error = function(e) NULL)

  if (.broken(eval2))  ## try another way 
    eval2 = tryCatch(
      eval(eval(parse(text = substitute(deparse(substitute(nvalue)))),
                parent.frame()), gr2dt(gr2), parent.frame(2)), error = function(e) NULL)
  
  if (.broken(eval1))  
  {
    warning(sprintf('value expression failed for edges for graph 1, defaulting to 1'))
    eval1 = rep(1, length(jj1))
  }
  
  if (.broken(eval2))  ## give up
  {
    warning(sprintf('value expression failed for edges for graph 2, defaulting to 1'))
    eval2 = rep(1, length(jj2))
  }

  jjm = merge.Junction(jj1, jj2, pad = pad, cartesian = TRUE, all = TRUE)$dt
  jdist = 0
  if (nrow(jjm)>0)
    {
      jjm$value = eval1[jjm$query.id]
      jjm$value.1 = eval1[jjm$subject.id]

      ## some more acrobatics to parse eweight as an expression

      ewt = tryCatch(eval(eval(parse(text = substitute(deparse(eweight))),
                               parent.frame()), jjm, parent.frame(2)), error = function(e) NULL)

      if (.broken(ewt))  ## try another way
        ewt = tryCatch(eval(eval(parse(text = substitute(deparse(substitute(eweight)))),
                                 parent.frame()), jjm, parent.frame(2)), error = function(e) NULL)
      
      if (.broken(ewt))  ## give up 
      {
        warning(sprintf('weight expression failed for edges, defaulting to 1'))
        ewt = rep(1, nrow(jjm))
      }
      
      jjm$weight = ewt
      jjdat = jjm[, .(value, value.1, weight)]
      jjdat[is.na(value), value := fill]
      jjdat[is.na(value.1), value.1 := fill]
      
      if (normalize) ## normalize by the weights
        jdist = jjdat[, diff := weight*lfun(value, value.1)][, sum(diff, na.rm = na.rm)/sum(weight, na.rm = na.rm)]
      else
        jdist = jjdat[, diff := weight*lfun(value, value.1)][, sum(diff, na.rm = na.rm)]
    }

  dist = node.coef*ndist + edge.coef*jdist

  if (L==2)
    dist = sqrt(dist)
  return(dist)
}



#' @name merge
#' @title merge for undefined number of Junction objects
#'
#' @export
merge = function(...) {
    UseMethod("merge")
}


#' @name merge.Junction
#' @title merge junctions by overlaps with padding
#'
#' Merges a set of junctions and keeps "seen.by" metadata of junction origin
#' using the argument names to this function
#'
#' If cartesian = TRUE, can only merge a pair of junction objects but then the output contains the overlapping 
#' junctions and metadata annotated with a $query.id (index into first argument) and $subject.id 
#' (index to into second argument) with deduped outputs,
#' 
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
#' @param cartesian whether to do a pairwise merge of all junction pairs in two junction objects x and y, which can potentially result in more rows than the number of inputs, Note: only works when there are exactly two inputs x and y
#' @param all only applicable if cartesian = TRUE, logical flag specifying whether to keep the junctions and metadata for non-overlapping junction pairs from both x and y inputs, aka "outer join" + "inner join"
#' @param all.x only applicable if cartesian = TRUE, logical flag specifying whether to keep the junctions and metadata for non-overlapping junction pairs from j1 in the output aka "left join" + "inner join""
#' @param all.y only applicable if cartesian = TRUE, logical flag specifying whether to keep the junctions and metadata for non-overlapping junction pairs from j2 in the output aka "right join" + "inner join"
#' @param ind  logical flag (default FALSE) specifying whether the "seen.by" fields should contain indices of inputs (rather than logical flags) and NA if the given junction is missing
#' @export merge.Junction
#' @export
"merge.Junction" = function(..., pad = 0, ind = FALSE, cartesian = FALSE, all = FALSE, all.x = all, all.y = all)
{
  list.args = list(...)
  if (cartesian)
  {
    if (length(list.args)!=2)
      stop('cartesian mode requires exactly two arguments')

    ra1 = list.args[[1]]
    ra2 = list.args[[2]]
    ov = ra.overlaps(ra1$grl, ra2$grl, pad = pad)
    ov = ov[!is.na(ov[,1]), , drop = FALSE]
    if (nrow(ov)==0)
      out.grl = ra1[c()]$grl
    else
      {
        out.grl = ra1$grl[ov[,1]]
        if (ncol(ra2$dt))
          values(out.grl) = cbind(values(out.grl),ra2$dt[ov[,2]])
        values(out.grl)$query.id = ov[,1]
        values(out.grl)$subject.id = ov[,2]
      }

    .cgrl = function(x, y)
    {
      tmp1 = x
      values(tmp1) = NULL
      tmp2 = y
      values(tmp2) = NULL
      out = c(tmp1, tmp2)
      values(out) = rbind(as.data.table(values(x)), as.data.table(values(y)), fill = TRUE)
      out
    }

    if (all.x && length(ra1) && length(ix <- setdiff(1:length(ra1), ov[,1])))
    {
      left.grl = ra1[ix]$grl
      values(left.grl)$query.id = ix
      out.grl = .cgrl(out.grl, left.grl) 
    }
    
    if (all.y && length(ra2) && length(ix <- setdiff(1:length(ra2), ov[,2])))
    {
      right.grl = ra2[ix]$grl
      values(right.grl)$subject.id = ix
      out.grl = .cgrl(out.grl, right.grl) 
    }

    names(values(out.grl)) = dedup(names(values(out.grl)), itemize.all = TRUE)
    Junction$new(out.grl)
  }
  else
  {
    if (all(sapply(list.args, length)==0))
      return(jJ())
    Junction$new(do.call(gGnome::ra.merge, c(lapply(list.args, function(x) x$grl), list(pad = pad))))
  }
}
registerS3method("merge", "Junction", merge.Junction, envir = globalenv())

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


setGeneric("lengths")
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
            return(base::lengths(x$snode.id))
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
#' @exportMethod seqinfo
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
#' @exportMethod seqinfo
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
#' @exportMethod seqinfo
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
#' @exportMethod seqlengths
#' @export
setMethod("seqlengths", c("gGraph"),
          function(x) {
            return(seqlengths(x$gr))
          })


#' @name seqlevels
#' @title seqlevels
#' @description
#'
#' @param x a gGraph object
#'
#' @return the seqlevels of this graph
#' @exportMethod seqlevels
#' @export
setMethod("seqlevels", c("gGraph"),
          function(x) {
            return(seqlevels(x$gr))
          })




#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' @param x a gNode object
#'
#' @return the seqlengths of this gNode
#' @exportMethod seqlengths
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
#' @exportMethod seqlengths
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
#' @exportMethod seqlengths
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
#' @exportMethod seqlengths
#' @export
setMethod("seqlengths", c("gWalk"),
          function(x) {
            return(seqlengths(x$graph))
          })


#' @name width
#' @title width
#' @description
#'
#' @param x a gWalk
#'
#' @return the width of each walk in this gWalk
#' @exportMethod width
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
#' @exportMethod width
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
#' @exportMethod refresh
#' @export
setMethod("refresh", "gWalk",
          function(x) {
            gW(snode.id = x$snode.id, sedge.id = x$sedge.id, circular = x$circular, graph = x$graph, meta = x$meta)
          })

#' @name gGraph.subset
#' @title subset gGraph on overlaps
#' @description
#'
#' @param x a gGraph
#' @param y a GRanges or string coercible to one
#' @return the portion of the gGraph that overlaps y (with loose ends wherever nodes / edges are broken)
#' @export
setMethod("%&%", signature(x = 'gGraph'), function(x, y) {
  ix = x$nodes$gr %^% y
  return(x[ix])
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
    y = unlist(y)

  if (is.character(y)){
    y = parse.gr(y)
  }
  ix = gr.in(x$gr, y)
  return(x[ix])
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

  ix = grl.in(x$grl, y)
  return(x[ix])
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

  ix = grl.in(x$grl, y)
  return(x[ix])
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

  ix = grl.in(x$grl, y)
  return(x[ix])
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
#' @usage jJ(rafile,
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
#' @param hg \code{character}, human genome version
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
jJ = function(rafile = NULL,
              keep.features = T,
              seqlengths = NULL,
              chr.convert = T,
              geno=NULL,
              hg=NULL,
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
                    hg=hg,
                    flipstrand = flipstrand,
                    swap.header = swap.header,
                    breakpointer = breakpointer,
                    seqlevels = seqlevels,
                    force.bnd = force.bnd,
                    skip = skip,
                    get.loose = get.loose)
         )
}


##########
########## KH additional functions
##########

## #' @name eclusters2
## #' @title gGraph R6 public method eclusters2
## #' 
## #' @description
## #' Marks ALT edges belonging (quasi) reciprocal cycles
## #'
## #' 
## #' @param juncs GRangesList of junctions
## #' @param mc.cores parallel
## #' @param only_chains TRUE will only pair breakend to its nearest nearest neighbor IFF the nearest neighbor is reciprocal, see arguments to "strict" for 3 different matching heuristics
## #' @param max.small size below which simple dups and dels are excluded
## #' @param weak logical flag if TRUE will not differentiate between cycles and paths and will return all weakly connected clusters in the junction graph [FALSE]
## #' @param strict Only active if only_chains = TRUE. Can be one of "strict", "one_to_one", or "loose". \cr
## #' "strict": each breakend can only be "monogamously" matched to one other breakend and if the nearest breakend is of the wrong orientation, it is thrown out. \cr
## #' "one_to_one" breakends can only be coupled to a single "monogamous" match but without considering breakends of the wrong orientation. If breakend B is nearest to C, B will only be matched to C. If A's nearest breakend is B but is further away than C, A will not be matched to B. \cr
## #' "loose": the nearest breakend in the correct orientation under the threshold is considered. The same as one_to_one except A will be matched to B, while B will be matched to C. \cr
## #' @param ignore.isolated If TRUE, all simple duplications, duplications nested with only other duplications, and simple deletions without any breakends will be thrown out.
## #' @return gGraph object with edges marked with ecluster id and metadata in gEdge, gNode, and $meta
## #' @author Marcin Imielinski
## eclusters2 = function (thresh = 1000, weak = TRUE, paths = !weak,
##                            mc.cores = 1, verbose = FALSE, chunksize = 1e+30, method = "single",
##                            return_pairs = FALSE, ignore.small = TRUE,
##                            max.small = 1e4, ignore.isolated = TRUE,
##                            strict = c("strict", "one_to_one", "loose"),
##                            min.isolated = max.small,
##                            only_chains = FALSE) {


##   if (!is.character(strict) || any(!strict %in% c("strict", "one_to_one", "loose"))) {
##     stop("strict must be one of 'strict', 'one_to_one', or 'loose'")
##   } else if (length(strict) > 1) {
##     strict = "one_to_one"
##   }
##   self$edges$mark(ecluster = as.integer(NA))
##   altedges = self$edges[type == "ALT", ]
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   if (ignore.small) {
##       altedges = altedges[!((class == "DUP-like" | class == "DEL-like") & altedges$span <= max.small)]
##   }
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   deldup = altedges[class %in% c("DUP-like", "DEL-like")]
##   ## below removes non-nested events below size threshold
##   ## and simple or nested duplications with no subsumed breakends
##   if (length(deldup) > 0 && ignore.isolated) {
##     altes = deldup$shadow
##     ## altes$sedge.id = altedges[class %in% c("DUP-like", "DEL-like")]$dt[altes$id]$sedge.id
##     bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix", "class", "sedge.id")]
##     bp$sedge.id.y = bp$sedge.id; bp$sedge.id = NULL
##     addon = deldup$dt[altes$id][, .(sedge.id, class)]
##     altes$sedge.id = addon$sedge.id
##     altes$class = addon$class
##     altes$nbp = altes %N% bp # number of breakpoints of any SV that fall within segment
##     numsum = altedges$shadow %>% gr.sum # using the shadows of all of the SVs not just dels and dups
##     altes = altes %$% numsum
##     iso = ((altes) %Q% (score == 1.0))$id
##     ## rm.edges = unique(altes[iso] %Q% (width < thresh))$sedge.id ## old
##     rm.edges = unique(altes[iso] %Q% (width < min.isolated))$sedge.id
##     rm.dups = S4Vectors::with(altes, sedge.id[class == "DUP-like" & nbp <= 2])
##     rm.dups = c(rm.dups, dedup.cols(gr2dt(altes %*% bp))[sedge.id != sedge.id.y][class == "DUP-like"][, .(all(1:2 %in% grl.iix), class.1 = class.1[1]), by = .(sedge.id, sedge.id.y)][, all(V1 == TRUE) & all(class.1 == "DUP-like"), by = sedge.id][V1 == TRUE]$sedge.id) # removing dups that have only other nested dups 
##     rm.edges = union(rm.edges, rm.dups)
##     keepeid = setdiff(altedges$dt$sedge.id, rm.edges)
##     altedges = altedges[as.character(keepeid)]
##   } # ignoring isolated dup and del edges that are smaller than threshold
##   if (verbose & weak)
##     message("Computing weak eclusters")
##   if (length(altedges) == 0) {
##     if (verbose) {
##       gmessage("No junction in this graph")
##     }
##     return(NULL)
##   }
##   bp = grl.unlist(altedges$grl)[, c("grl.ix", "grl.iix", "edge.id")]
##   bp$m.ix = seq_along(bp)
##   bp.dt = gr2dt(bp)
##   ix = split(1:length(bp), ceiling(runif(length(bp)) * ceiling(length(bp)/chunksize)))
##   ixu = unlist(ix)
##   eps = 1e-09
##   ## ij = do.call(rbind, split(1:length(bp), bp$grl.ix))
##   xt.adj = xt.adj0 = Matrix::sparseMatrix(1, 1, x = 0, dims = rep(length(bp),
##                                                         2))
##   if (verbose) {
##     message(sprintf("Computing junction graph across %s ALT edges with distance threshold %s",
##                     length(altedges), thresh))
##   }
##   if (!exists(".INF")) {
##     .INF = pmax(sum(seqlengths(self)), 1e+09)
##   }
##   bp.pair = as.matrix(dcast(data.table(ix = bp$m.ix,
##                              grl.ix = bp$grl.ix,
##                              grl.iix = bp$grl.iix),
##                   grl.ix ~ grl.iix, value.var = "ix")[,2:3])
##   bp.pair = rbind(bp.pair, cbind(bp.pair[,2], bp.pair[,1]))
##   ifun = function(iix, ignore.strand = FALSE,
##                   verbose = FALSE, eps = 1e-9) {
##     if (verbose > 1)
##       cat(".")
##     tmpm = gr.dist(bp[iix], gr.flipstrand(bp), ignore.strand = ignore.strand) +
##       eps
##     return(as(tmpm, "Matrix"))
##   }
##   xt.adj[ixu, ] = do.call(rbind,
##                           mclapply(ix, ifun, ignore.strand = FALSE,
##                                    mc.cores = mc.cores))
##   diag(xt.adj) = NA_real_
##   ## only_chains = TRUE, enforcing that nearest breakpoints are only considered
##   if (only_chains) {
    
##     ## enforcing that no distances between breakends from the same junction  are considered
##     xt.adj[bp.pair] = NA_real_
##     xt.adj0[ixu, ] = do.call(rbind,
##                              mclapply(ix, ifun, ignore.strand = TRUE,
##                                       mc.cores = mc.cores)) # this is to find nearest breakends, regardless of orientation
##     ## enforcing that self-to-self breakend distances are not considered
##     diag(xt.adj0) = NA_real_
##     ## enforcing that no distances between breakends from the same junction  are considered
##     xt.adj0[bp.pair] = NA_real_
##     suppressWarnings({
##       nearest_ix = dunlist(lapply(
##         seq_len(nrow(xt.adj)),
##         function(x) which(xt.adj[x,] == min(xt.adj[x,], na.rm = T)))) %>%
##         as.matrix
##       nearest0_ix = dunlist(lapply(
##         seq_len(nrow(xt.adj0)),
##         function(x) which(xt.adj0[x,] == min(xt.adj0[x,], na.rm = T)))) %>%
##         as.matrix
##     })
##     if (nrow(nearest_ix) & nrow(nearest0_ix)) {
##       if (strict == "strict") {
##         nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
##         nearest0_ix = cbind(rowMins(nearest0_ix), rowMaxs(nearest0_ix))
##         nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
##         nearest0_ix = nearest0_ix[duplicated(nearest0_ix),,drop = FALSE]
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest0_ix = rbind(nearest0_ix, cbind(nearest0_ix[,2], nearest0_ix[,1]))
##         nearest_ix = as.matrix(merge(nearest_ix, nearest0_ix)) # if there is a breakend closer but in the wrong orientation, that distance will be thrown out downstream, i.e. there is a one to one match of breakend and any nearest breakend that is in the wrong orientation disqualifies the clustser
##       } else if (strict == "one_to_one") {
##         nearest_ix = cbind(rowMins(nearest_ix), rowMaxs(nearest_ix))
##         nearest_ix = nearest_ix[duplicated(nearest_ix),,drop = FALSE]
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest_ix = as.matrix(nearest_ix) # only one-to-one breakends are considered
##       } else if (strict == "loose") {
##         nearest_ix = rbind(nearest_ix, cbind(nearest_ix[,2], nearest_ix[,1]))
##         nearest_ix = nearest_ix[!duplicated(nearest_ix),,drop = FALSE]
##       }
##     }
##   } else if (!only_chains) {
##     xt.adj[bp.pair] = 1
##   }
##   ## rm(xt.adj0)
##   adj = xt.adj
##   xt.adj[which(is.na(as.matrix(xt.adj)))] = .INF + 1
##   adj[which(is.na(as.matrix(adj)))] = 0
##   adj[which(as.matrix(adj) > thresh)] = 0
##   if (only_chains && nrow(nearest_ix)) {
##     tmp = xt.adj[nearest_ix]
##     xt.adj[] = .INF + 1
##     adj[] = 0
##     xt.adj[nearest_ix] = tmp
##     adj[nearest_ix] = tmp
##   }
##   dt = Matrix::which(xt.adj < thresh, arr.ind = T)
##   dt = unique(data.table(cbind(rowMins(dt), rowMaxs(dt))))
##   dt[, bp.dist := xt.adj[dt[, cbind(V1, V2)]]]
##   dt$sign = ifelse(strand(bp[dt$V1]) == "+",
##             ifelse(gr.flipstrand(bp[dt$V1]) > bp[dt$V2], -1, 1),
##             ifelse(gr.flipstrand(bp[dt$V1]) < bp[dt$V2], -1, 1))
##   dt2 = dt[,idx := seq_len(.N)] %>% melt(measure.vars = c("V1", "V2"))
##   meta = cbind(gr2dt(bp)[dt2$value], dt2)[order(idx)]
##   meta = rbind(meta,
##               gr2dt(bp)[paste(grl.ix, grl.iix) %nin% meta[, paste(grl.ix, grl.iix)]],
##               fill = T)
##   if (return_pairs) {
##     return(meta)
##   }
##   hcl = stats::hclust(as.dist(xt.adj), method = "single")
##   hcl.lbl = cutree(hcl, h = thresh)
##   bp.dt$hcl = hcl.lbl
##   bp.hcl = bp.dt[, .(hcl.1 = .SD[grl.iix == 1, hcl], hcl.2 = .SD[grl.iix ==
##                                                                  2, hcl]), keyby = grl.ix]
##   altedges$mark(hcl.1 = bp.hcl[.(seq_along(altedges)), hcl.1])
##   altedges$mark(hcl.2 = bp.hcl[.(seq_along(altedges)), hcl.2])
##   hcl.ig = igraph::graph_from_edgelist(bp.hcl[, unique(cbind(hcl.1,
##                                                              hcl.2))], directed = FALSE)
##   hcl.comp = components(hcl.ig)
##   altedges$mark(ehcl = as.integer(hcl.comp$membership)[bp.hcl[,
##                                                               hcl.1]])
##   adj[adj > thresh] = 0
##   refg = self[, type == "REF"]
##   bpp = Matrix::which(adj != 0, arr.ind = TRUE)
##   dref = pdist(bp[bpp[, 1]], bp[bpp[, 2]])
##   drefg = diag(refg$dist(bp[bpp[, 1]], bp[bpp[, 2]]))
##   ix = which(drefg > dref)
##   if (length(ix))
##     adj[bpp[ix, , drop = FALSE]] = FALSE
##   if (verbose > 1)
##     cat("\n")
##   adj = adj | t(adj)
##   junpos = bp1 = bp$grl.iix == 1
##   junneg = bp2 = bp$grl.iix == 2
##   adj2 = adj & FALSE
##   adj2[junpos, junpos] = adj[bp2, bp1]
##   adj2[junpos, junneg] = adj[bp2, bp2]
##   adj2[junneg, junpos] = adj[bp1, bp1]
##   adj2[junneg, junneg] = adj[bp1, bp2]
##   if (verbose)
##     message(sprintf("Created basic junction graph using distance threshold of %s",
##                     thresh))
##   cl = split(1:length(bp), igraph::clusters(graph.adjacency(adj2),
##                                             ifelse(weak, "weak", "strong"))$membership)
##   cl = cl[S4Vectors::elementNROWS(cl) > 1]
##   cl = cl[order(S4Vectors::elementNROWS(cl))]
##   ## browser()
##   ## jcl = lapply(cl, function(x) unique(sort(bp$grl.ix[x])))
##   jcl = lapply(cl, function(x) unique(sort(bp$edge.id[x])))
##   jcls = sapply(jcl, paste, collapse = " ")
##   jcl = jcl[!duplicated(jcls)]
##   adj3 = adj2
##   altedges$mark(ecycle = as.character(NA))
##   if (length(jcl) > 0) {
##     dcl = dunlist(unname(jcl))[, `:=`(listid, paste0(ifelse(weak,
##                                                             "", "c"), listid))]
##     if (!weak)
##         ## altedges[dcl$V1]$mark(ecycle = dcl$listid)
##         altedges[as.character(dcl$V1)]$mark(ecycle = dcl$listid)
##     ## altedges[dcl$V1]$mark(ecluster = dcl$listid)
##     altedges[as.character(dcl$V1)]$mark(ecluster = dcl$listid)
##     meta = merge(meta, altedges$dt[, .(edge.id, ecluster, hcl.1, hcl.2, ehcl, ecycle)], by = "edge.id")
##     meta[!is.na(ecluster),
##          `:=`(
##            nclust = length(unique(edge.id)),
##            all_positive = all(replace(sign, is.na(sign),  3e9) > 0),
##            all_negative = all(replace(sign, is.na(sign), -3e9) < 0),
##            mixed = {naom = na.omit(sign); any(naom > 0) & any(naom < 0)},
##            bridge = anyNA(bp.dist)
##            ),
##          by = ecluster]
##     meta[, `:=`(
##         num_positive = sum(sign[!duplicated(idx)] > 0, na.rm = T),
##         num_negative = sum(sign[!duplicated(idx)] < 0, na.rm = T)
##     ),
##     by = ecluster]
##     self$set(recip_bp = meta)
##     if (verbose) {
##         message("Annotated weakly connected junction clusters and added ecluster pair metadata")
##     }
##   }
##   if (verbose)
##     message(sprintf("Annotated %s junction cycles in edge field $ecycle",
##                     length(jcl)))
##   if (paths & !weak) {
##     if (verbose)
##       message("Analyzing paths")
##     if (length(jcl) > 0) {
##       adj3[unlist(jcl), unlist(jcl)] = FALSE
##     }
##     sinks = Matrix::which(Matrix::rowSums(adj3) == 0)
##     sources = Matrix::which(Matrix::colSums(adj3) == 0)
##     cl2 = split(1:length(bp), igraph::clusters(graph.adjacency(adj3),
##                                                "weak")$membership)
##     cl2 = cl2[S4Vectors::elementNROWS(cl2) > 1]
##     if (any(ix <- S4Vectors::elementNROWS(cl2) > 2)) {
##       cl3 = do.call(c, mclapply(cl2[ix], function(x) {
##         tmp.adj = adj3[x, x]
##         lapply(all.paths(tmp.adj, sources = sources,
##                          sinks = sinks, verbose = verbose)$paths, function(i) x[i])
##       }, mc.cores = mc.cores))
##       cl2 = c(cl2[!ix], cl3)
##     }
##     jcl2 = lapply(cl2, function(x) unique(sort(bp$grl.ix[x])))
##     jcls2 = sapply(jcl2, paste, collapse = " ")
##     jcl2 = jcl2[!duplicated(jcls2)]
##     altedges$mark(epath = as.character(NA))
##     if (length(jcl2) > 0) {
##       dcl2 = dunlist(unname(jcl2))[, `:=`(listid, paste0("p",
##                                                          listid))]
##       altedges[dcl2$V1]$mark(epath = dcl2$listid)
##       self$edges$mark(ecluster =
##                         ifelse(is.na(self$edges$dt$ecycle) &
##                                is.na(self$edges$dt$epath),
##                                as.character(NA),
##                                paste0(ifelse(is.na(self$edges$dt$ecycle), "",
##                                              self$edges$dt$ecycle),
##                                       ifelse(is.na(self$edges$dt$epath),
##                                              "", self$edges$dt$epath))))
##     }
##     if (verbose)
##       message(sprintf("Annotated %s paths in edge field $epath",
##                       length(jcl2)))
##   }
##   return(invisible(self))
## }
## ## gGraph$public_methods$eclusters2 = tmpeclustpairs # if replacing a binding

## #' @description
## #' make eclusters
## gGraph$public_methods$eclusters2 = NULL; gGraph$set("public", "eclusters2", eclusters2)



#' @name copy3
#' @title make deep copy, recursively
#'
#' useful for dev
#' makes deep copy of R6 object, S4 object, or anything else really
#'
copy3 = function (x, recurse_list = TRUE) {
    if (inherits(x, "R6")) {
        x2 = rlang::duplicate(x$clone(deep = T))
        for (name in intersect(names(x2$.__enclos_env__), c("private", 
            "public"))) for (nname in names(x2$.__enclos_env__[[name]])) tryCatch({
            x2$.__enclos_env__[[name]][[nname]] = copy3(x2$.__enclos_env__[[name]][[nname]])
        }, error = function(e) NULL)
        return(x2)
    } else if (isS4(x)) {
        x2 = rlang::duplicate(x)
        slns = slotNames(x2)
        for (sln in slns) {
            tryCatch({
                slot(x2, sln) = copy3(slot(x2, sln))
            }, error = function(e) NULL)
        }
        return(x2)
    } else if (inherits(x, c("list"))) {
        x2 = rlang::duplicate(x)
        x2 = rapply(x2, copy3, how = "replace")
        return(x2)
    } else {
        x2 = rlang::duplicate(x)
        return(x2)
    }
}


#' @name rep_len2
#' @title recycle vector along length OR nrow of object
#'
#'
#' @author Kevin Hadi
#' @param x data
#' @param objalong any object to recycle x along if uselen = TRUE, or an actual integer value if uselen = FALSE
#' @return vector
rep_len2 = function(x, objalong, uselen = TRUE) {
    if (uselen)
        rep(x, length.out = NROW(objalong))
    else
        rep(x, length.out = objalong)
}

#' @name seq_along2
#' @title seq along either row of table or length of vector
#'
#'
#' @author Kevin Hadi
#' @param x data
#' @return vector
seq_along2 = function(x)  {
  seq_len(NROW(x))
}


#' @name rleseq
#' @title numbers up within repeating elements of a vector
#'
#' @description
#' returns unique id within each unique element of a vector or set of provided vectors
#' and also a running id within each unique element
#'
#' @param ... Vector(s) to identify with unique id and a running id within each unique id
#' @param clump a logical specifying if duplicates are to be counted together
#' @param recurs a logical that is meant to only be set by the function when using clump = TRUE
#' @return a list of idx and seq
#' @author Kevin Hadi
rleseq = function (..., clump = TRUE, recurs = FALSE, na.clump = TRUE, 
                   na.ignore = FALSE, sep = paste0(" ", rand.string(length = 6), 
                     " "), use.data.table = TRUE) 
{
  force(sep)
  out = if (use.data.table) {
    tryCatch(
    {
      dt = data.table(...)
      setnames(dt, make.names(rep("", ncol(dt)), unique = T))
      ## make.unique
      cmd = sprintf("dt[, I := .I][, .(idx = .GRP, seq = seq_len(.N), lns = .N, I), by = %s]", mkst(colnames(dt), "list"))
      dt = eval(parse(text = cmd))
      setkey(dt, I)[, .(idx, seq, lns)]
    }, error = function(e) structure("data table didn't work...", class = "err"))
  }
  if (!(is.null(out) || class(out)[1] == "err"))
    return(as.list(out))
  rand.string <- function(n = 1, length = 12) {
    randomString <- c(1:n)
    for (i in 1:n) {
      randomString[i] <- paste(sample(c(0:9, letters, LETTERS), 
        length, replace = TRUE), collapse = "")
    }
    return(randomString)
  }
  if (isTRUE(na.clump)) 
    paste = function(..., sep) base::paste(..., sep = sep)
  else paste = function(..., sep) base::paste(stringr::str_c(..., 
    sep = sep))
  lns = base::lengths(list(...))
  if (!all(lns == lns[1])) 
    warning("not all vectors provided have same length")
  fulllens = max(lns, na.rm = T)
  vec = setNames(paste(..., sep = sep), seq_len(fulllens))
  if (length(vec) == 0) {
    out = list(idx = integer(0), seq = integer(0), lns = integer(0))
    return(out)
  }
  if (na.ignore) {
    isnotna = which(rowSums(as.data.frame(lapply(list(...), 
      is.na))) == 0)
    out = list(idx = rep(NA, fulllens), seq = rep(NA, fulllens), 
      lns = rep(NA, fulllens))
    if (length(isnotna)) 
      vec = vec[isnotna]
    tmpout = do.call(rleseq, c(alist(... = vec), alist(clump = clump, 
      recurs = recurs, na.clump = na.clump, na.ignore = FALSE, use.data.table = FALSE)))
    for (i in seq_along(out)) out[[i]][isnotna] = tmpout[[i]]
    return(out)
  }
  if (!isTRUE(clump)) {
    rlev = rle(vec)
    if (isTRUE(recurs)) {
      return(unlist(unname(lapply(rlev$lengths, seq_len))))
    }
    else {
      out = list(idx = rep(seq_along(rlev$lengths), times = rlev$lengths), 
        seq = unlist(unname(lapply(rlev$lengths, seq_len))))
      out$lns = ave(out[[1]], out[[1]], FUN = length)
      return(out)
    }
  }
  else {
    if (!isTRUE(na.clump)) {
      vec = replace2(vec, which(x == "NA"), dedup(dg(x)[dg(x) == 
                                                          "NA"]))
    }
    vec = setNames(vec, seq_along(vec))
    lst = split(vec, factor(vec, levels = unique(vec)))
    ord = as.integer(names(unlist(unname(lst))))
    idx = rep(seq_along(lst), times = base::lengths(lst))
    out = list(idx = idx[order(ord)], seq = rleseq(idx, clump = FALSE, 
      recurs = TRUE, use.data.table = FALSE)[order(ord)])
    out$lns = ave(out[[1]], out[[1]], FUN = length)
    return(out)
  }
}


#' @name match3
#' @title similar to setkey except a general use utility
#'
#' very slow version of keying a la data.table
#' but for general/interactive use
#' @author Kevin Hadi
match3 = function(x, table, nomatch = NA_integer_, old = TRUE, use.data.table = TRUE) {
  out = if (use.data.table) {
    tryCatch({
      dx = data.table(x = x)[, id.x := seq_len(.N)]
      dtb = data.table(table = table)[, id.tb := seq_len(.N)]
      ## setkey(dx, x)[list(dtb$table)]$id.x
      setkey(dtb, table)[list(dx$x)]$id.tb
    }, error = function(e) structure("err", class = "err"))
  }
  if (!is.null(out) && !inherits(out, "err")) return(out)
  if (old) {
    dx = within(data.frame(x = x), {id.x = seq_along(x)})
    dtb = within(data.frame(table = table), {id.tb = seq_along(table)})
    res = merge(dx, dtb, by.x = "x", by.y = "table", all.x = TRUE,
      allow.cartesian = TRUE)
    return(res$id.tb[order(res$id.x)])
  } else  {
    m = match(table,x)
    mat = cbind(m, seq_along(m))
    mat = mat[!is.na(mat[,1]),,drop=FALSE]
    mat = mat[order(mat[,1], na.last = FALSE),,drop = FALSE]
    mat = cbind(mat, seq_len(dim(mat)[1]))
    m2 = match(x,table)
    ix = which(!duplicated(m2) & !is.na(m2))
    mat_rix = unlist(rep(split(mat[,3], mat[,1]), base::tabulate(m2)[m2][ix]))
    ## mat_rix = unlist(rep(split(mat[,3], mat[,1]), base::tabulate(m2)[m2][ix]))
    ix = rep(1, length.out = length(m2))
    ## original line
    ## ix[!is.na(m2)] = base::tabulate(m)[!is.na(m2)]
    ix[!is.na(m2)] = base::tabulate(m)[m][m2][!is.na(m2)]
    out = rep(m2, ix)
    out[!is.na(out)] = mat[mat_rix,,drop=F][,2]
    return(out)
    ## m = match(table, x)
    ## mat = cbind(m, seq_along(m))
    ## mat = mat[!is.na(mat[, 1]), , drop = FALSE]
    ## mat = mat[order(mat[, 1]), , drop = FALSE]
    ## mat = cbind(mat, seq_len(dim(mat)[1]))
    ## m2 = match(x, table)
    ## ix = which(!duplicated(m2))
    ## mat_rix = unlist(rep(split(mat[, 3], mat[, 1]), base::tabulate(m2)[m2][ix]))
    ## mat[mat_rix, , drop = F][, 2]
  }
}


#' @name dedup.cols
#' @title applies dedup to colnames
#'
#' dedup the column names of a data.frame/data.table
#'
#' @return A data.table or data.frame
dedup.cols = function(tbl, remove = FALSE) {
    if (remove) {
        if (!inherits(tbl, "data.table"))
            return(tbl[, match(unique(colnames(tbl)), colnames(tbl))])
        else
            return(tbl[, match(unique(colnames(tbl)), colnames(tbl)), with = FALSE])
    } else {
            colnames(tbl) = base::make.unique(colnames(tbl))
            return(tbl)
    }
}


#' @name rowMins
#' @title rowMins hack
#'
#'
#' @return A vector
rowMins = function(x) {
  do.call(pmin, as.data.frame(x))
}

#' @name rowMaxs
#' @title rowMaxs hack
#'
#'
#' @return A vector
rowMaxs = function(x) {
  do.call(pmax, as.data.frame(x))
}

#' @name gr.noval
#' @title get rid of mcols on GRanges/GRangesLists
#' @description
#'
#' remove all metadata from GRanges or GRangesList
#'
#' @return GRanges or GRangesList
#' @author Kevin Hadi
gr.noval = function(gr, keep.col = NULL, drop.col = NULL) {
    if (is.null(keep.col) & is.null(drop.col)) {
        select_col = NULL
    } else {
        all_col = colnames(gr@elementMetadata)
        if (inherits(gr, "GRangesList")) {
            all_col = c(all_col, colnames(gr@unlistData@elementMetadata))
        }

        if (!is.null(keep.col) & is.null(drop.col)) {
            select_col = intersect(all_col, keep.col)
        } else if (is.null(keep.col) & !is.null(drop.col)) {
            select_col = setdiff(all_col, drop.col)
        } else if (!is.null(keep.col) && !is.null(drop.col)) {
            if (intersect(keep.col, drop.col) > 0) {
                warning("drop.col and keep.col args have overlapping elements\nkeeping the columns that overlap")
                select_col = intersect(setdiff(all_col, setdiff(drop.col, keep.col)), keep.col)
            }
        }
    }
    if (inherits(gr, "GRangesList")) {
        tmp_query = intersect(select_col, colnames(gr@unlistData@elementMetadata))
        gr@unlistData@elementMetadata = gr@unlistData@elementMetadata[,c(tmp_query), drop = FALSE]
    }
    tmp_query = intersect(select_col, colnames(gr@elementMetadata))
    gr@elementMetadata = gr@elementMetadata[,c(tmp_query),drop = FALSE]
    return(gr)
}


#' @name setcols
#' @title convenience function to set columns
#'
#' sets columns of an object
#'
#' @param dt data frame/table or matrix
#' @param old integer or character or logical vector corresponding to current colnames in dt
#' @param new character vector for new column names
#' @return colnamed object
#' @author Kevin Hadi
setcols = function(dt, old, new) {
  if (inherits(dt, c("GRanges", "GRangesList"))) {
    mcols(dt) = setcols(mcols(dt), old, new)
    return(dt)
  }
  cnames = colnames2(dt)
  if (missing(new) || missing(old)) {
    if (missing(old)) {
      old = new
    }
    if (is.character(old) && length(old) == length(cnames)) {
      colnames(dt) = old
      return(dt)
    } else {
      stop("names provided must be same length as ncol(dt)")
    }
  }
  if (is.character(old)) {
    out = merge(data.frame(cnames, seq_along(cnames)), data.frame(cnames = old, new = new),
      allow.cartesian = T)
    cnames[out[[2]]] = out[[3]]
    colnames(dt) = cnames
    return(dt)
  }
  if (is.logical(old)) {
    if (! length(old) == length(cnames)) stop("logical vector must be same length as ncol(dt)")
    old = which(old)
  }
  cnames[old] = new
  colnames(dt) = cnames
  return(dt)
}

#' @name colnames2
#' @title robust colnames
#'
#' gives back character vector same number of columns of input regardless whether named or not
#'
#' @param str a path string
#' @return a string with multiple parentheses replaced with a single parenthesis
#' @author Kevin Hadi
colnames2 = function(x) {
    nm = colnames(x)
    if (is.null(nm))
        return(rep("", length.out = NCOL(x)))
    else
        return(nm)
}


#' @name isNA
#' @title is.na but also tests for "NA" character
#'
#' @description
#'
#' @author Kevin Hadi
isNA = function(x, na.char = c("NA", "NULL", "na", "null")) {
    if (is.character(x)) {
        return(is.na(x) | x %in% na.char)
    } else {
        return(is.na(x))
    }
}

#' @name na2false
#' @title replace logical vector with NA to FALSE
#'
#' @description
#' A convenience function to set a logical vector with NAs to false
#'
#' @return A logical vector with NAs set to FALSE
#' @author Kevin Hadi
na2false = function(v)
{
    v[isNA(v)] = FALSE
    ## mode(v) = "logical"
    v
}

#' @name rand.string
#' @title make a random string
#'
#' @return random string
#' @author Someone from Stackoverflow
rand.string <- function(n=1, length=12)
{
    randomString <- c(1:n)                  # initialize vector
    for (i in 1:n)
    {
        randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                        length, replace=TRUE),
                                 collapse="")
    }
    return(randomString)
}

#' @name gr.spreduce
#' @title reduce based on a field(s) to split by in elementMetadata of GRanges, or given vector
#' @description
#'
#' split and reduce GRanges by field(s)
#' if providing a variable not already within the GRanges,
#' may need to use dynget(variable_name)
#'
#' @return GRanges
#' @author Kevin Hadi
gr.spreduce = function(gr,  ..., ignore.strand = FALSE, pad = 0, return.grl = FALSE, sep = paste0(" ", rand.string(length = 8), " ")) {
  lst = as.list(match.call())[-1]
  ix = which(!names(lst) %in% c("gr", "sep", "pad", "ignore.strand", "return.grl"))
  vars = unlist(sapply(lst[ix], function(x) unlist(sapply(x, toString))))
  if (length(vars) == 1) {
    if (!vars %in% colnames(mcols(gr)))
      vars = tryCatch(unlist(list(...)), error = function(e) vars)
  }
  if (!all(vars %in% colnames(mcols(gr))))
    stop("Must specify valid metadata columns in gr")
  tmpix = do.call(
    function(...) paste(..., sep = sep),
    as.list(mcols(gr)[,vars, drop = F]))
  unix = which(!duplicated(tmpix))
  tmpix = factor(tmpix, levels = tmpix[unix])
  grl = unname(gr.noval(gr) %>% GenomicRanges::split(tmpix))
  grl = GenomicRanges::reduce(grl + pad, ignore.strand = ignore.strand)
  if (return.grl) {
    mcols(grl) = mcols(gr)[unix,vars,drop = F]
    return(grl)
  } else {
    out = unlist(grl)
    mcols(out) = mcols(gr)[rep(unix, times = IRanges::width(grl@partitioning)),
      vars,drop = F]
    return(out)
  }
}

#' @name duped
#' @title duped
#'
#'
#'
#' @param ... vectors to paste by
#' @author Kevin Hadi
duped = function(..., binder = "data.table") {
    duplicated(tryCatch(et(sprintf("%s(...)", binder)),
                        ## error = function(e) do.call(cbind, list(...))))
                        error = function(e) paste(...)))
}

#' @name merge.repl
#' @title merging data tables with collapsing columns with the same name
#'
#' Merge two data tables with various replacing strategies
#' for columns common between x and y that are not used to merge
#' (i.e. not specified in the "by" argument)
#'
#' @param replace_NA logical, only use values in dt.y, any dt.x not in dt.y is clobbered (NA)
#' @param force_y logical, should x and y common columns be merged?
#' @param overwrite_x logical, if force_y = TRUE, should NA values in y replace x?
#' @return A data.table
#' @author Kevin Hadi
merge.repl = function(dt.x,
                      dt.y,
                      sep = "_",
                      replace_NA = TRUE,
                      force_y = TRUE,
                      overwrite_x = FALSE,
                      keep_order = FALSE,
                      keep_colorder = TRUE,
                      keep_factor = TRUE,
                      ...) {
    arg_lst = as.list(match.call())
    by.y = eval(arg_lst[['by.y']], parent.frame())
    by.x = eval(arg_lst[['by.x']], parent.frame())
    by = eval(arg_lst[['by']], parent.frame())
    all.x = eval(arg_lst[['all.x']], parent.frame())
    all.y = eval(arg_lst[['all.y']], parent.frame())
    all = eval(arg_lst[['all']], parent.frame())
    allow.cartesian = eval(arg_lst[['allow.cartesian']])
    key_x = key(dt.x)
    if (is.null(all.x)) {
        all.x = TRUE
    }
    if (is.null(all.y)) {
        all.y = FALSE
    }
    if (!is.null(all) && all) {
        all.y = TRUE
        all.x = TRUE
    }
    if (is.null(allow.cartesian)) {
        allow.cartesian = FALSE
    }
    if (!inherits(dt.x, "data.table")) {
        dt.x = as.data.table(dt.x)
    }
    if (!inherits(dt.y, "data.table")) {
        dt.y = as.data.table(dt.y)
    }
    if (keep_order == TRUE) {
        dt.x[['tmp.2345098712340987']] = seq_len(nrow(dt.x))
    }

    dt.x[['in.x.2345098712340987']] = rep(TRUE, length.out = nrow(dt.x))
    dt.y[['in.y.2345098712340987']] = rep(TRUE, length.out = nrow(dt.y))

    new_ddd_args = list(by = by, by.x = by.x, by.y = by.y, all.x = all.x, all.y = all.y, allow.cartesian = allow.cartesian)

    if (is.null(by.y) & is.null(by.x) & is.null(by)) {

        if (length(attributes(dt.x)[['sorted']]) > 0 &&
            length(attributes(dt.y)[['sorted']]) > 0) {
            k.x = key(dt.x)
            k.y = key(dt.y)
        } else {
            k.y = k.x = intersect(names2(dt.x), names2(dt.y))
            if (length(k.x) == 0)
                stop("no common columns to merge by!")
            message("intersecting by: ", paste(k.x, collapse = ", "))
            new_ddd_args[['by']] = k.x
        }
        if (is.null(k.x) | is.null(k.y) || (k.x != k.y)) {
            stop("neither by.x/by.y  nor by are supplied, keys of dt.x and dt.y must be identical and non NULL")
        }
        x.cols = setdiff(names(dt.x), k.x)
        y.cols = setdiff(names(dt.y), k.y)

    } else if (!is.null(by.x) & !is.null(by.y)) {

        x.cols = setdiff(names(dt.x), by.x)
        y.cols = setdiff(names(dt.y), by.y)
        new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by"))]

    } else if (!is.null(by)) {

        x.cols = setdiff(names(dt.x), by)
        y.cols = setdiff(names(dt.y), by)
        if (! all(by %in% colnames(dt.x)) | ! all(by %in% colnames(dt.y))) {
            stop("column ", by, " does not exist in one of the tables supplied \nCheck the column names")
        }
        new_ddd_args = new_ddd_args[setdiff(names(new_ddd_args), c("by.y", "by.x"))]

    }
    these_cols = intersect(x.cols, y.cols)
    ## if (replace_in_x) {
    if (!replace_NA) {
        dt.x.tmp = copy(dt.x)
        for (this_col in these_cols) {
            data.table::set(dt.x.tmp, i = NULL, j = this_col, value = NULL)
        }
        dt.repl = suppressWarnings(do.call("merge", args = c(list(x = dt.x.tmp, y = dt.y), new_ddd_args)))
        ## dt_na2false(dt.repl, c("in.x.2345098712340987", "in.y.2345098712340987"))
    } else {
        dt.repl = suppressWarnings(do.call("merge", args = c(list(x = dt.x, y = dt.y), new_ddd_args)))
        dt_na2false(dt.repl, c("in.x.2345098712340987", "in.y.2345098712340987"))
        in.x = which(dt.repl[["in.x.2345098712340987"]])
        in.y = which(dt.repl[["in.y.2345098712340987"]])
        this_env = environment()
        for (this_col in these_cols) {
            x_cname = paste0(this_col, ".x")
            y_cname = paste0(this_col, ".y")
            x_col = dt.repl[[x_cname]]
            y_col = dt.repl[[y_cname]]
            xf = inherits(x_col, "factor")
            yf = inherits(y_col, "factor")
            if ( {xf || yf} && keep_factor) {
                if (!xf) { x_col = factor(x_col); xf = TRUE }
                if (!yf) { y_col = factor(y_col); yf = TRUE }
            }
            if (xf && !keep_factor) { x_col = as.character(x_col); xf = FALSE } 
            if (yf && !keep_factor) { y_col = as.character(y_col); yf = FALSE }
            if (force_y) {
                if (!overwrite_x) {
                    ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                    ##     new_col = factor(y_col, forcats::lvls_union(list(y_col, x_col)))
                    ##     new_col[is.na(new_col)] = x_col[is.na(new_col)]
                    ## } else {
                    ##     new_col = ifelse(!is.na(y_col), y_col, x_col)
                    ## }                    
                    if (xf || yf) {
                        new_col = factor(y_col, forcats::lvls_union(list(y_col, x_col)))
                        new_col[is.na(new_col)] = x_col[is.na(new_col)]
                    } else {
                        new_col = ifelse(!is.na(y_col), y_col, x_col)
                    }
                } else {
                    ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                    ##     new_col = factor(x_col, forcats::lvls_union(list(y_col, x_col)))
                    ## } else {
                    ##     new_col = x_col
                    ## }
                    ## new_col[dt.repl[['in.y.2345098712340987']]] = y_col[dt.repl[['in.y.2345098712340987']]]
                    if (xf || yf) {
                        new_col = factor(x_col, forcats::lvls_union(list(y_col, x_col)))
                    } else {
                        new_col = x_col
                    }
                    new_col[in.y] = y_col[in.y]
                }
            } else {
                ## if (inherits(x_col, "factor") & inherits(y_col, "factor")) {
                ##     new_col = factor(x_col, forcats::lvls_union(list(x_col, y_col)))
                ##     new_col[is.na(new_col) & !is.na(y_col)] = y_col[is.na(new_col) & !is.na(y_col)]
                ## } else {
                ##     new_col = ifelse(is.na(x_col) & !is.na(y_col), y_col, x_col)
                ## }
                if (xf | yf) {
                    new_col = factor(x_col, forcats::lvls_union(list(x_col, y_col)))
                    new_col[is.na(new_col) & !is.na(y_col)] = y_col[is.na(new_col) & !is.na(y_col)]
                } else {
                    new_col = ifelse(is.na(x_col) & !is.na(y_col), y_col, x_col)
                }
            }
            data.table::set(dt.repl, j = c(x_cname, y_cname, this_col), value = list(NULL, NULL, this_env[["new_col"]]))
        }
    }
    ## } else if (!replace_in_x & !is.null(suffix)) {
    ##     y.suff.cols = paste0(y.cols, sep, suffix)
    ##     ## dt.y.tmp = copy(dt.y)[, eval(dc(y.suff.cols)) := eval(dl(y.cols))][, eval(dc(y.cols)) := NULL]
    ##     dt.y.tmp = copy(dt.y)
    ##     data.table::set(dt.y, j = y.suff.cols, value = dt.y[, y.cols, with = FALSE])
    ##     data.table::set(dt.y, j = y.cols, value = NULL)
    ##     ## dt.repl = merge(dt.x, dt.y.tmp, all.x = TRUE, ...)
    ##     dt.repl = do.call("merge", args = c(list(x = dt.x, y = dt.y.tmp), new_ddd_args))
    ## }
    if (keep_order == TRUE) {
        data.table::setorderv(dt.repl, "tmp.2345098712340987")
        dt.repl[['tmp.2345098712340987']] = NULL
    }
    data.table::set(dt.repl, j = c("in.y.2345098712340987", "in.x.2345098712340987"),
                    value = list(NULL, NULL))
    if (keep_colorder) {
        x_cols = colnames(dt.x)
        ## get the order of columns in dt.repl in order of X with
        ## additional columns tacked on end
        setcolorder(dt.repl,
                    intersect(union(colnames(dt.x), colnames(dt.repl)),
                              colnames(dt.repl)))
    }
    return(dt.repl)
}


#' @name dt_na2false
#' @title convert columns with NA to false
#'
#' coerce NA in columns of class "logical" to FALSE
#'
#' @param dt data.table
#' @param these_cols NULL by default, will select columns of class logical, otherwise will be specified
#' @return A data.table
#' @author Kevin Hadi
dt_na2false = function(dt, these_cols = NULL) {
    na2false = function(v)
    {
        ## v = ifelse(is.na(v), v, FALSE)
        v[is.na(v)] = FALSE
        as.logical(v)
    }
    if (is.null(these_cols)) {
        these_cols = which(sapply(dt, class) == "logical")
    }
    for (this_col in these_cols) {
        ## this_val = as.data.frame(dt[, this_col, with = FALSE])[,1]
        this_val = dt[[this_col]]
        data.table::set(dt, j = this_col, value = na2false(this_val))
    }
    return(dt)
}
