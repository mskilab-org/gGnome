% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gGnome.R
\name{gWalk.subset}
\alias{gWalk.subset}
\alias{gWalk}
\title{sets metadata of gWalk object
(accessible through $dt accessor)}
\usage{
\S4method{\%&\%}{gWalk}(x, y)
}
\arguments{
\item{x}{a gWalk}
}
\value{
the subset of gWalks overlapping with the query
}
\description{
sets metadata of gWalk object
(accessible through $dt accessor)

sets metadata of gWalk object
(accessible through $dt accessor)

subset gWalk on overlaps
}
\author{
Marcin Imielinski

Julie Behr
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-gWalk-new}{\code{gWalk$new()}}
\item \href{#method-gWalk-set}{\code{gWalk$set()}}
\item \href{#method-gWalk-subset}{\code{gWalk$subset()}}
\item \href{#method-gWalk-dts}{\code{gWalk$dts()}}
\item \href{#method-gWalk-fix}{\code{gWalk$fix()}}
\item \href{#method-gWalk-rep}{\code{gWalk$rep()}}
\item \href{#method-gWalk-print}{\code{gWalk$print()}}
\item \href{#method-gWalk-json}{\code{gWalk$json()}}
\item \href{#method-gWalk-eval}{\code{gWalk$eval()}}
\item \href{#method-gWalk-disjoin}{\code{gWalk$disjoin()}}
\item \href{#method-gWalk-simplify}{\code{gWalk$simplify()}}
\item \href{#method-gWalk-mark}{\code{gWalk$mark()}}
\item \href{#method-gWalk-fitcn}{\code{gWalk$fitcn()}}
\item \href{#method-gWalk-gtrack}{\code{gWalk$gtrack()}}
\item \href{#method-gWalk-clone}{\code{gWalk$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-new"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$new(
  snode.id = NULL,
  sedge.id = NULL,
  grl = NULL,
  graph = NULL,
  meta = NULL,
  circular = NULL,
  disjoin = FALSE,
  drop = FALSE
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-set"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-set}{}}}
\subsection{Method \code{set()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$set(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-subset"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-subset}{}}}
\subsection{Method \code{subset()}}{
Allows for subsetting of the gWalk Object using bracket notation
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$subset(i)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-dts"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-dts}{}}}
\subsection{Method \code{dts()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$dts(ix = 1:self$length, makelists = TRUE)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-fix"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-fix}{}}}
\subsection{Method \code{fix()}}{
Returns a modified version of the gWalk with seqlevels "fixed", including
keeping only certain seqlevels, dropping certain seqlevels, and replacing seqlevels.

Warning: this may modify the walks including getting rid of nodes and edges (i.e. those outside
the retained seqlevels) and also change coordinates (ie move ranges that were previously on different
chromosomes to the same chromosome etc.).  Use with caution!

Default behavior is to replace 'chr', with ''.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$fix(pattern = NULL, replacement = NULL, drop = TRUE, seqlengths = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{pattern}}{character pattern to replace in seqlevels (used in a gsub, can have backreferences)}

\item{\code{replacement}}{character to replace pattern with (used in a gsub, can have backreferences)}

\item{\code{drop}}{character vector of seqlevels to drop or logical TRUE to drop all seqlevels that are unused  (TRUE)}

\item{\code{seqlengths}}{new seqlengths i.e. named integer vector of seqlevels to drop or embed graph into}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
current graph modified in place with additional nodes and edges, as specified by user
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-rep"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-rep}{}}}
\subsection{Method \code{rep()}}{
Creates "bubbles" in the underlying graph by replicating the gwalks.  Node replication
replicates edges going in and out of all replicated nodes.  If an edge connects a pair of replicated nodes
that edge will be replicated across all pairs of those replciated nodes.   Walk replication will create "longer bubbles"
with fewer edges getting replicated i.e. it will only replicate intra-walk edges within each walk replicate (but not between
separate walk replicates).

(note that this returns a copy of the gGraph that this node is linked to, i.e. it does not make 
this gWalk stale)

New graph keeps track of the parent node and edge ids in the original graph using node metadata parent.node.id
and edge metadata parent.edge.id i.e. the replicated nodes will be connected to the sources of the original nodes
and if replicated nodes connect to each other, then there will exist an edge connecting
all of their instances to each other.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$rep(times)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{times}}{scalar or vector of length self$length specifying how many times to replicate each of the nodes.}

\item{\code{nodes}}{= gNode object must point to a node in the graph, can also be an index of a node (but not a metadata expression), can also be a gWalk object}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
Returns a pointer to the new nodes
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-print"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-print}{}}}
\subsection{Method \code{print()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-json"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-json}{}}}
\subsection{Method \code{json()}}{
Dumps walk (and underlying graph)
to json file
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$json(
  filename = ".",
  save = TRUE,
  verbose = FALSE,
  annotations = NULL,
  nfields = NULL,
  efields = NULL,
  stack.gap = 1e+05,
  include.graph = TRUE,
  settings = list(y_axis = list(title = "copy number", visible = TRUE)),
  cid.field = NULL,
  no.y = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{filename}}{character path to save to}

\item{\code{save}}{whether to save or return list object representing json contents}

\item{\code{annotations}}{which graph annotations to dump to json}

\item{\code{nfields}}{which node fields to dump to json (NULL)}

\item{\code{efields}}{which edge fields to dump to json (NULL)}

\item{\code{stack.gap}}{this currently does not do anything, but is intended to control the spread of intervals on the y-axis}

\item{\code{include.graph}}{if set to TRUE, then the output will include the underlying graph for the gWalk (TRUE)}

\item{\code{settings}}{for details see gGraph json() module}

\item{\code{no.y}}{for details see gGraph json() module}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-eval"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-eval}{}}}
\subsection{Method \code{eval()}}{
Evaluates an expression on gWalk node or edge metadata
and returns a vector of length length(self) one per walk.

Ideally the expression should generate a scalar value for each walk.id, though will dedup and re-order
even if it does not.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$eval(x, node, edge)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{data.table style expression on node or edge metadata (will try each)}

\item{\code{node}}{data.table style expresion on node metadata}

\item{\code{edge}}{data.table style expresion on edge metadata}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-disjoin"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-disjoin}{}}}
\subsection{Method \code{disjoin()}}{
similar to gGraph disjoin but also updates according gWalk
and allows disjoining around a set of gRanges, essentially
adding "breaks" to the walks and allowing the ranges in those
walks to inherit metadata from overlapping GRanges
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$disjoin(
  gr = NULL,
  graph = NULL,
  by = NULL,
  collapse = TRUE,
  na.rm = TRUE,
  avg = FALSE,
  sep = ",",
  FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg)
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-simplify"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-simplify}{}}}
\subsection{Method \code{simplify()}}{
gWalk simplify will merge reference adjacent nodes in walk +/- "by field"
(i.e. which requires the nodes to additional match on by field)
and simplify the underlying graph.

Note that this only simplifies walks;
the associated graph may not be completely simplified 
since we retain all walk endpoints
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$simplify(
  by = NULL,
  na.rm = TRUE,
  avg = FALSE,
  sep = ",",
  FUN = default.agg.fun.generator(na.rm = na.rm, sep = sep, avg = avg)
)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-mark"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-mark}{}}}
\subsection{Method \code{mark()}}{
marks the nodes and edges in the graph corresponding to this walk
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$mark(...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-fitcn"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-fitcn}{}}}
\subsection{Method \code{fitcn()}}{
Given gWalks on a graph with node and edge
metadata "cn" attempts to find a non-negative integer
combination of those walks that generates those walks
minimizes the number of unique walks used (L0 norm).

Modifications to this objecti include (default behavior)
if min.alt = TRUE then solution will minimizes the number of
walks that contain an ALT edge.  If weight is provided
then the objective function be \eqn{\sum_i w_i*I[x_i>0]} where
w_i is the provided weight for walk i and I[x_i>0] is indicator
that is 1 if walk i has a nonzero copy number in the solution
and 0 otherwise. 

Complains if "cn" is not a metadata of nodes and edges
and if the graph is not "balanced" ie either
every node cn >= sum(edge cn) for both incoming
and outgoing edges.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$fitcn(
  trim = TRUE,
  weight = NULL,
  obs.mat = NULL,
  verbose = FALSE,
  min.alt = TRUE,
  edgeonly = FALSE,
  evolve = FALSE,
  n.sol = 2
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{weight}}{either a numeric metadata field of self or a  self$length integer weight}

\item{\code{verbose}}{logical flag whether to output CPLEX output}

\item{\code{min.alt}}{logical flag whether to make objective function}

\item{\code{edgeonly}}{logical}

\item{\code{evolve}}{logical}

\item{\code{n.sol}}{numeric}

\item{\code{logical}}{flag trim whether to "trim" the graph relative to the footprint of the walks}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-gtrack"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-gtrack}{}}}
\subsection{Method \code{gtrack()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$gtrack(name = NULL, stack.gap = 1e+05, ...)}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-gWalk-clone"></a>}}
\if{latex}{\out{\hypertarget{method-gWalk-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{gWalk$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
